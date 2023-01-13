   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE CFA_Method_Module                                                            !##!
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
USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Poseidon_Message_Routines_Module, &
            ONLY :  Run_Message,        &
                    Warning_Message

USE Poseidon_Parameters, &
            ONLY :  Degree,                     &
                    L_Limit,                    &
                    Cur_Iteration,              &
                    Max_Iterations,             &
                    Convergence_Criteria,       &
                    Verbose_Flag

USE Poseidon_IO_Parameters, &
            ONLY :  CFA_Var_Names


USE Variables_Mesh, &
            ONLY :  Num_R_Elements,             &
                    R_Inner,                    &
                    R_Outer

USE Variables_Derived, &
            ONLY :  Var_Dim,                    &
                    Prob_Dim

USE Variables_FP,  &
            ONLY :  FP_Anderson_M

USE CFA_Load_Vector_Module, &
            ONLY :  Calc_CFA_Load_Vector,          &
                    Allocate_CFA_Source_Variables,   &
                    Deallocate_CFA_Source_Variables

USE CFA_System_Solvers_Module, &
            ONLY :  Solve_CFA_Systems

USE IO_Print_Results, &
            ONLY :  Print_Results

USE IO_Convergence_Output,  &
            ONLY :  Output_Convergence_Reports

USE IO_Write_Final_Results, &
            ONLY :  Write_Final_Results

USE Functions_Coeffs_Module, &
            ONLY :  Coeff_To_Vector,    &
                    Vector_To_Coeff

USE Timer_Routines_Module, &
            ONLY :  TimerStart,                     &
                    TimerStop

USE Timer_Variables_Module, &
            ONLY :  Timer_CFA_Load_Vector


IMPLICIT NONE





CONTAINS


!+101+###########################################################################!
!                                                                                !
!           Fixed_Point_AndersonM                                                !
!                                                                                !
!################################################################################!
SUBROUTINE CFA_Method()

LOGICAL                                                 ::  CONVERGED
!LOGICAL                                                 ::  CONVERGED_Residual
INTEGER                                                 ::  i


INTEGER                                                 ::  M
INTEGER                                                 ::  LWORK
INTEGER                                                 ::  mk
INTEGER                                                 ::  INFO

CHARACTER(LEN = 300)                                    ::  Message

REAL(idp),DIMENSION(:),   ALLOCATABLE                :: Resid_Vector
REAL(idp),DIMENSION(:,:), ALLOCATABLE                :: FVector
REAL(idp),DIMENSION(:,:), ALLOCATABLE                :: GVector
REAL(idp),DIMENSION(:),   ALLOCATABLE                :: BVector
REAL(idp),DIMENSION(:),   ALLOCATABLE                :: UVector
REAL(idp),DIMENSION(:),   ALLOCATABLE                :: GVectorM
REAL(idp),DIMENSION(:),   ALLOCATABLE                :: FVectorM
REAL(idp),DIMENSION(:,:), ALLOCATABLE                :: AMatrix
REAL(idp),DIMENSION(:),   ALLOCATABLE                :: Work
REAL(idp),DIMENSION(:),   ALLOCATABLE                :: Alpha


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

CALL Allocate_CFA_Source_Variables()

CUR_ITERATION = 0
CONVERGED = .FALSE.


IF ( Verbose_Flag ) CALL Run_Message('Beginning Fixed Point Iterations.')


!
!   Begin Method
!




!
!   Calculate Source Vector with u_0
!
CALL TimerStart(Timer_CFA_Load_Vector)
CALL Calc_CFA_Load_Vector()
CALL TimerStop(Timer_CFA_Load_Vector)





DO WHILE ( .NOT. CONVERGED  .AND. Cur_Iteration < Max_Iterations)
    Cur_Iteration = Cur_Iteration+1

    IF ( Verbose_Flag ) THEN
        WRITE(Message,'(A,A,A,I2.2,A)')'Starting Iteration ',Cur_Iteration,'.'
        CALL Run_Message(TRIM(Message))
    END IF

    mk = MIN(Cur_Iteration, M)

    IF ( Cur_Iteration .NE. 1 ) THEN
        UVector = GVectorM
    ELSE
        CALL Coeff_To_Vector( UVector )
    END IF



    !
    !   Solve Systems
    !
    CALL Solve_CFA_Systems()

    CALL Coeff_To_Vector( GVector(:,mk) )
    
    FVector(:,mk) = GVector(:,mk) - UVector(:)



    IF ( mk == 1 ) THEN
        GVectorM = GVector(:,mk)
    ELSE


        BVector = -FVector(:,mk)
        AMatrix(:,1:mk-1) = FVector(:,1:mk-1) - SPREAD( FVector(:,mk), DIM=2, NCOPIES = mk-1)

        CALL ZGELS('N',Prob_Dim,mk-1,1,              &
                    AMatrix(:,1:mk-1), Prob_Dim,     &
                    BVector, Prob_Dim,                &
                    WORK, LWORK, INFO )

        IF ( INFO .NE. 0 ) THEN
            WRITE(Message,'(A,I1.1,A,I1.1)')'In CFA_Method_Module, ZGELS failed with INFO = ',INFO
            CALL WARNING_MESSAGE(Message)
        END IF

        Alpha(1:mk-1) = BVector(1:mk-1)
        Alpha(mk)     = 1.0_idp - SUM( Alpha(1:mk-1) )


        GVectorM = 0.0_idp
        DO i = 1,mk
            GVectorM = GVectorM + Alpha(i)*GVector(:,i)
        END DO

    END IF

    FVectorM = GVectorM - UVector(:)



    ! Check for Convergence
    CALL Convergence_Check( FVectorM, Cur_Iteration, Converged )





    IF ( mk == M .AND. .NOT. Converged ) THEN
        GVector = CSHIFT( GVector, SHIFT = +1, DIM = 2)
        FVector = CSHIFT( FVector, SHIFT = +1, DIM = 2)
    END IF




!   Update Coefficient Vectors
    CALL Vector_To_Coeff( GVectorM )



!   Calculate Source Vector with new solution
    CALL TimerStart(Timer_CFA_Load_Vector)
    CALL Calc_CFA_Load_Vector()
    CALL TimerStop(Timer_CFA_Load_Vector)



    IF ( Verbose_Flag ) THEN
        WRITE(Message,'(A,A,A,I2.2,A)')'Ending Iteration ',Cur_Iteration,'.'
        CALL Run_Message(TRIM(Message))
        WRITE(*,'()')
    END IF



END DO ! Converged Loop


CALL Deallocate_CFA_Source_Variables()




CALL Write_Final_Results()

CALL Output_Convergence_Reports()

END SUBROUTINE CFA_Method





 !+201+########################################################!
!                                                               !
!           Convergence_Check                                   !
!                                                               !
 !#############################################################!
SUBROUTINE Convergence_Check( Update, Iter, Flag )


REAL(idp),DIMENSION(Var_Dim), INTENT(IN)         :: Update
INTEGER,                         INTENT(IN)         :: Iter
LOGICAL,                         INTENT(INOUT)      :: Flag

CHARACTER(LEN = 300)                                ::  Message


IF ( Verbose_Flag ) THEN
    WRITE(Message,'(A,ES22.15,A)')'Maximum change in the coefficient vector : ',MAXVAL(ABS(Update)),'.'
    CALL Run_Message(TRIM(Message))
END IF



IF ( ALL( ABS( Update ) <= Convergence_Criteria ) ) THEN
    IF ( Verbose_Flag ) THEN
        CALL Run_Message("The Method has converged.")
        CALL Run_Message("The absolute update is less than the tolerance set. ")
    END IF
    Flag = .TRUE.
    
END IF


IF ( Iter == Max_Iterations ) THEN
    CALL Warning_Message('FP_Accelerated has reached the maximum number of allowed iterations.')
    Flag = .TRUE.
END IF




END SUBROUTINE Convergence_Check



END MODULE CFA_Method_Module
