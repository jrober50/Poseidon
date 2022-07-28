   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE FP_AndersonM_Module                                                          !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!   +101+    Fixed_Point_AndersonM                                               !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!

!*D*================================!
!                                   !
!           Dependencies            !
!                                   !
!===================================!
USE MPI

USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Variables_Mesh, &
            ONLY :  Num_R_Elements,             &
                    R_Inner,                    &
                    R_Outer

USE Variables_Derived, &
            ONLY :  LM_Length,                  &
                    Num_R_Nodes,                &
                    Var_Dim,                    &
                    Beta_Prob_Dim,              &
                    Prob_Dim

USE Variables_IO, &
            ONLY :  Frame_Update_Table,         &
                    Write_Flags,                &
                    Print_Flags


USE Poseidon_Parameters, &
            ONLY :  DEGREE,                     &
                    L_LIMIT,                    &
                    CUR_ITERATION,              &
                    MAX_ITERATIONS,             &
                    Convergence_Type,           &
                    Convergence_Type_Names,     &
                    Convergence_Criteria,       &
                    Convergence_Flag,           &
                    NUM_CFA_EQs,                &
                    CFA_Eq_Flags,               &
                    Verbose_Flag

USE Variables_FP,  &
            ONLY :  FP_Anderson_M

USE FP_Load_Vector_Module, &
            ONLY :  Calc_FP_Load_Vector,          &
                    Allocate_FP_Source_Variables,   &
                    Deallocate_FP_Source_Variables

USE IO_Print_Results, &
            ONLY :  Print_Results

USE IO_Convergence_Output,  &
            ONLY :  Output_Convergence_Reports

USE IO_Write_Final_Results, &
            ONLY :  Write_Final_Results

USE FP_System_Solvers_Module,   &
            ONLY :  Solve_FP_System,                &
                    Solve_FP_System_Beta

USE FP_Functions_Coeffs_Module, &
            ONLY :  Coeff_To_Vector,    &
                    Vector_To_Coeff

USE Timer_Routines_Module, &
            ONLY :  TimerStart,                     &
                    TimerStop

USE Timer_Variables_Module, &
            ONLY :  Timer_FP_Load_Vector,         &
                    Timer_FP_CFLF_Solve,            &
                    Timer_FP_Beta_Solve


IMPLICIT NONE





CONTAINS


!+101+###########################################################################!
!                                                                                !
!           Fixed_Point_AndersonM                                                !
!                                                                                !
!################################################################################!
SUBROUTINE Fixed_Point_AndersonM()

LOGICAL                                                 ::  CONVERGED
!LOGICAL                                                 ::  CONVERGED_Residual
INTEGER                                                 ::  i


INTEGER                                                 ::  M
INTEGER                                                 ::  LWORK
INTEGER                                                 ::  mk
INTEGER                                                 ::  INFO

REAL(KIND = idp), DIMENSION(1:4)                        :: timer


COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: Resid_Vector
COMPLEX(idp),DIMENSION(:,:), ALLOCATABLE                :: FVector
COMPLEX(idp),DIMENSION(:,:), ALLOCATABLE                :: GVector
COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: BVector
COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: UVector
COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: GVectorM
COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: FVectorM
COMPLEX(idp),DIMENSION(:,:), ALLOCATABLE                :: AMatrix
COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: Work
COMPLEX(idp),DIMENSION(:),   ALLOCATABLE                :: Alpha

LOGICAL                                                 :: PR = .FALSE.

M = FP_Anderson_M
LWORK = 2*M

ALLOCATE( Resid_Vector(1:Prob_Dim) )
ALLOCATE( FVector(1:Prob_Dim,1:M) )
ALLOCATE( GVectorM(1:Prob_Dim) )
ALLOCATE( FVectorM(1:Prob_Dim) )
ALLOCATE( BVector(1:Prob_Dim) )
ALLOCATE( GVector(1:Prob_Dim,1:M) )
ALLOCATE( Work(1:LWORK) )
ALLOCATE( Alpha(1:M) )
ALLOCATE( AMatrix(1:Prob_Dim,1:M))
ALLOCATE( UVector(1:Prob_Dim) )

CALL Allocate_FP_Source_Variables()

timer = 0.0_idp
CUR_ITERATION = 0
CONVERGED = .FALSE.

IF ( Verbose_Flag ) THEN
    PRINT*,"Begining Fixed Point Iterative Solve."
END IF



!
!   Begin Method
!

!IF ( (Write_Flags(5) == 1) .OR. (Write_Flags(5) == 3) ) THEN
!    PR = .TRUE.
!    WRITE(*,'(A)')"Initial Guess Values"
!    CALL Print_Results()
!    PRINT*," "
!END IF



!
!   Calculate Source Vector with u_0
!
CALL TimerStart(Timer_FP_Load_Vector)
CALL Calc_FP_Load_Vector()
CALL TimerStop(Timer_FP_Load_Vector)





DO WHILE ( .NOT. CONVERGED  .AND. Cur_Iteration < Max_Iterations)


    Timer(4) = MPI_WTime()

    Cur_Iteration = Cur_Iteration+1

    mk = MIN(Cur_Iteration, M)

    IF ( Cur_Iteration .NE. 1 ) THEN
        UVector = GVectorM
    ELSE
        CALL Coeff_To_Vector( UVector )
    END IF



    !
    !   Solve Systems
    !
    IF ( Verbose_Flag ) THEN
        PRINT*,"In Anderson FP loop, Before System Solves."
    END IF

    IF ( ANY( CFA_EQ_Flags(1:2) == 1) ) THEN
        CALL TimerStart(Timer_FP_CFLF_Solve)
        Call Solve_FP_System()
        CALL TimerStop(Timer_FP_CFLF_Solve)
    END IF
    
    IF ( ANY( CFA_EQ_Flags(3:5) == 1) ) THEN
        CALL TimerStart(Timer_FP_Beta_Solve)
        Call Solve_FP_System_Beta()
        CALL TimerStop(Timer_FP_Beta_Solve)
    END IF




    IF ( Verbose_Flag ) THEN
        PRINT*,"In Anderson FP loop, Before Update."
    END IF

    CALL Coeff_To_Vector( GVector(:,mk) )

    
    FVector(:,mk) = GVector(:,mk) - UVector(:)



    IF ( mk == 1 ) THEN

        GVectorM = GVector(:,mk)

    ELSE
!        PRINT*,"Here",mk,M

        BVector = -FVector(:,mk)

        AMatrix(:,1:mk-1) = FVector(:,1:mk-1) - SPREAD( FVector(:,mk), DIM=2, NCOPIES = mk-1)

!        PRINT*,"Before ZGELS"
        CALL ZGELS('N',Prob_Dim,mk-1,1,              &
                    AMatrix(:,1:mk-1), Prob_Dim,     &
                    BVector, Prob_Dim,                &
                    WORK, LWORK, INFO )

!        PRINT*,"After ZGELS"
        Alpha(1:mk-1) = BVector(1:mk-1)
        Alpha(mk)     = 1.0_idp - SUM( Alpha(1:mk-1) )


        GVectorM = 0.0_idp
        DO i = 1,mk

            GVectorM = GVectorM + Alpha(i)*GVector(:,i)

        END DO

    END IF

    FVectorM = GVectorM - UVector(:)
    Frame_Update_Table(CUR_Iteration,1) = MAXVAL( ABS( FVectorM ) )


    IF ( Verbose_Flag ) THEN
        PRINT*,"L_Inf(FVectorM) = ",MAXVAL(ABS(FVectorM))
    END IF

    IF ( ALL( ABS( FVectorM ) <= Convergence_Criteria ) ) THEN
        IF ( Verbose_Flag ) THEN
            PRINT*,"The Method has converged. The absolute update is less than the tolerance set. "
        END IF
        Converged = .TRUE.
        
    END IF
!    IF ( ALL( ABS( FVectorM ) <= Convergence_Criteria*ABS( GVectorM ) ) ) THEN
!        PRINT*,"The Method has converged. The relative update is less than the tolerance set. "
!        Converged = .TRUE.
!    END IF


    IF ( mk == M .AND. .NOT. Converged ) THEN
        GVector = CSHIFT( GVector, SHIFT = +1, DIM = 2)
        FVector = CSHIFT( FVector, SHIFT = +1, DIM = 2)
    END IF



    IF ( Cur_Iteration == Max_Iterations ) THEN
        PRINT*,"FP_Accelerated has reached the maximum number of allowed iterations. "
        Converged = .TRUE.
    END IF


!   Update Coefficient Vectors
    CALL Vector_To_Coeff( GVectorM )


    IF ( Verbose_Flag ) THEN
        PRINT*,"In Anderson FP loop, Before calculating Source Vector."
    END IF



!   Calculate Source Vector with new solution
    CALL TimerStart(Timer_FP_Load_Vector)
    CALL Calc_FP_Load_Vector()
    CALL TimerStop(Timer_FP_Load_Vector)



    IF ( PR ) THEN
        WRITE(*,'(A,I3.3,A)')'Iteration ',Cur_Iteration,' Results'
        CALL Print_Results()
        PRINT*," "
    END IF


    IF ( Verbose_Flag .EQV. .TRUE. ) THEN
        WRITE(*,'(A,1X,I3.3,/)') "End of Iteration",Cur_Iteration
    END IF
!    STOP


END DO ! Converged Loop


CALL Deallocate_FP_Source_Variables()




CALL Write_Final_Results()

CALL Output_Convergence_Reports()

END SUBROUTINE Fixed_Point_AndersonM








END MODULE FP_AndersonM_Module