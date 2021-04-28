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
!##!   +101+    Fixed_Point_Method                                                  !##!
!##!                                                                                !##!
!##!   +201+    Check_FP_Convergence                                                !##!
!##!   +202+    Calc_FP_Residual                                                    !##!
!##!                                                                                !##!
!##!   +301+    Solve_FP_System                                                     !##!
!##!   +302+                                                                        !##!
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
            ONLY :  Frame_Report_Flag,          &
                    Write_Report_Flag,          &
                    Write_Results_Flag,         &
                    Iter_Report_Num_Samples,    &
                    Iter_Report_File_ID,        &
                    Frame_Report_File_ID,       &
                    Iter_Time_Table,            &
                    Frame_Update_Table,         &
                    Frame_Residual_Table


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
                    Verbose_Flag



USE Variables_FP,  &
            ONLY :  Matrix_Format,              &
                    Linear_Solver,              &
                    FP_Source_Vector,           &
                    FP_Coeff_Vector,            &
                    FP_Update_Vector,           &
                    FP_Laplace_Vector,          &
                    FP_Residual_Vector,         &
                    FP_Laplace_Vector_Beta,     &
                    FP_Residual_Vector_Beta,    &
                    Laplace_Matrix_Full,        &
                    Laplace_Matrix_VAL,         &
                    Laplace_Matrix_ROW,         &
                    Laplace_Matrix_COL,         &
                    Laplace_Factored_VAL,       &
                    Laplace_Factored_ROW,       &
                    Laplace_Factored_COL,       &
                    Laplace_Matrix_Beta,        &
                    FP_Source_Vector_Beta,      &
                    FP_Coeff_Vector_Beta,       &
                    CFA_EQ_Flags,               &
                    CFA_Eq_Map,                 &
                    CFA_MAT_Map,                &
                    Laplace_NNZ,                &
                    Num_Matrices,               &
                    MCF_Flag,                   &
                    FP_Anderson_M

USE Functions_Mesh, &
            ONLY :  Create_Logarithmic_1D_Mesh,     &
                    Create_Uniform_1D_Mesh

USE FP_Source_Vector_Module, &
            ONLY :  Calc_FP_Source_Vector,          &
                    Allocate_FP_Source_Variables,   &
                    Deallocate_FP_Source_Variables

USE FP_Functions_BC,  &
            ONLY :  Dirichlet_BC,                   &
                    Dirichlet_BC_CCS,               &
                    Dirichlet_BC_CHOL,              &
                    Dirichlet_BC_Beta,              &
                    Neumann_BC,                     &
                    Neumann_BC_CCS

USE IO_FP_Linear_System, &
            ONLY :  Output_Laplace_Beta,            &
                    Output_Laplace

USE IO_Print_Results, &
            ONLY :  Print_Results

USE Poseidon_Cholesky_Module,   &
            ONLY :  Cholesky_Factorization,         &
                    CCS_Back_Substitution,          &
                    CCS_Forward_Substitution

USE Linear_Solvers_And_Preconditioners, &
            ONLY :  PRECOND_CONJ_GRAD_CCS,          &
                    JACOBI_CONDITIONING,            &
                    Jacobi_Conditioning_Beta

USE Poseidon_IO_Module, &
            ONLY :  Clock_In,                           &
                    OPEN_ITER_REPORT_FILE,              &
                    CLOSE_ITER_REPORT_FILE,             &
                    OUTPUT_FINAL_RESULTS

USE FP_Functions_Mapping, &
            ONLY :  FP_LM_Map,                      &
                    FP_Array_Map

USE FP_System_Solvers_Module,   &
            ONLY :  Solve_FP_System,                &
                    Solve_FP_System_Beta

USE FP_Resid_Convergence_Module,    &
            ONLY :  Check_FP_Convergence

IMPLICIT NONE





CONTAINS






!+101+###########################################################################!
!                                                                                !
!           Fixed_Point_Method                                                   !
!                                                                                !
!################################################################################!
SUBROUTINE Fixed_Point_AndersonM()

LOGICAL                                                 ::  CONVERGED
LOGICAL                                                 ::  CONVERGED_Residual
INTEGER                                                 ::  i, k, lm_loc, ui, lm
INTEGER                                                 ::  here, there

INTEGER                                                 ::  M
INTEGER                                                 ::  LWORK
INTEGER                                                 ::  mk
INTEGER                                                 ::  INFO

REAL(KIND = idp), DIMENSION(1:4)                        :: timer

COMPLEX(idp), DIMENSION(1:NUM_R_NODES)                  ::  WORK_VEC
COMPLEX(idp), DIMENSION(1:NUM_R_NODES,1:NUM_R_NODES)    ::  WORK_MAT


COMPLEX(idp),DIMENSION(:), ALLOCATABLE                  :: Resid_Vector
COMPLEX(idp),DIMENSION(:,:), ALLOCATABLE                :: FVector
COMPLEX(idp),DIMENSION(:,:), ALLOCATABLE                :: GVector
COMPLEX(idp),DIMENSION(:), ALLOCATABLE                  :: BVector
COMPLEX(idp),DIMENSION(:), ALLOCATABLE                  :: UVector
COMPLEX(idp),DIMENSION(:), ALLOCATABLE                  :: GVectorM
COMPLEX(idp),DIMENSION(:), ALLOCATABLE                  :: FVectorM
COMPLEX(idp),DIMENSION(:,:), ALLOCATABLE                :: AMatrix
COMPLEX(idp),DIMENSION(:), ALLOCATABLE                  :: Work
COMPLEX(idp),DIMENSION(:), ALLOCATABLE                  :: Alpha


M = FP_Anderson_M
LWORK = 2*M

ALLOCATE( Resid_Vector(1:Var_Dim) )
ALLOCATE( FVector(1:Var_Dim,1:M) )
ALLOCATE( GVectorM(1:Var_Dim) )
ALLOCATE( FVectorM(1:Var_Dim) )
ALLOCATE( BVector(1:Var_Dim) )
ALLOCATE( GVector(1:Var_Dim,1:M) )
ALLOCATE( Work(1:LWORK) )
ALLOCATE( Alpha(1:M) )
ALLOCATE( AMatrix(1:Var_Dim,1:M))
ALLOCATE( UVector(1:Var_Dim) )

CALL Allocate_FP_Source_Variables()

timer = 0.0_idp
CUR_ITERATION = 0
CONVERGED = .FALSE.



!PRINT*,"In AndersonM"
!
!   Begin Method
!

IF (Verbose_Flag .EQV. .TRUE. ) THEN
    PRINT*,"Initial Guess"
    CALL Print_Results()
    PRINT*," "
END IF



!
!   Calculate Source Vector with u_0
!
IF ( Verbose_Flag ) THEN
    PRINT*,"Before Anderson FP loop, Before calculating Source Vector."
END IF
timer(1) = MPI_WTime()
CALL Calc_FP_Source_Vector()
timer(2) = MPI_Wtime()
Call Clock_In(timer(2)-timer(1),3)


ui = 1

DO WHILE ( .NOT. CONVERGED  .AND. Cur_Iteration < Max_Iterations)



    Timer(4) = MPI_WTime()

    Cur_Iteration = Cur_Iteration+1

    mk = MIN(Cur_Iteration, M)

    IF ( Cur_Iteration .NE. 1 ) THEN
        UVector = GVectorM
    ELSE
        DO lm_loc = 1,LM_Length
!            PRINT*,FP_Coeff_Vector(:,lm_loc,ui)
            here = (lm_loc-1)*Num_R_Nodes + 1
            there = lm_loc*Num_R_Nodes
            UVector(here:there) = FP_Coeff_Vector(:,lm_loc,ui)
        END DO
    END IF


    !
    !   Solve Systems
    !
    IF ( Verbose_Flag ) THEN
        PRINT*,"In Anderson FP loop, Before System Solves."
    END IF
    timer(1) = MPI_Wtime()
    IF ( ANY( CFA_EQ_Flags(1:2) == 1) ) THEN
        Call Solve_FP_System()
    END IF
    timer(2) = MPI_WTime()
    IF ( ANY( CFA_EQ_Flags(3:5) == 1) ) THEN
        Call Solve_FP_System_Beta()
    END IF
    timer(3) = MPI_Wtime()

    CALL Clock_In(timer(2)-timer(1),4)
    CALL Clock_In(timer(3)-timer(2),5)





    IF ( Verbose_Flag ) THEN
        PRINT*,"In Anderson FP loop, Before Update."
    END IF




    DO lm_loc = 1,LM_Length
!        PRINT*,FP_Coeff_Vector(:,lm_loc,ui)
        here = (lm_loc-1)*Num_R_Nodes + 1
        there = lm_loc*Num_R_Nodes
        GVector(here:there,mk) = FP_Coeff_Vector(:,lm_loc,ui)
    END DO

    FVector(:,mk) = GVector(:,mk) - UVector(:)



    IF ( mk == 1 ) THEN

        GVectorM = GVector(:,mk)

    ELSE
!        PRINT*,"Here",mk,M

        BVector = -FVector(:,mk)

        AMatrix(:,1:mk-1) = FVector(:,1:mk-1) - SPREAD( FVector(:,mk), DIM=2, NCOPIES = mk-1)

!        PRINT*,"Before ZGELS"
        CALL ZGELS('N',Var_Dim,mk-1,1,              &
                    AMatrix(:,1:mk-1), Var_Dim,     &
                    BVector, Var_Dim,                &
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






    IF ( ALL( ABS( FVectorM ) <= Convergence_Criteria*ABS( GVectorM ) ) ) THEN
        PRINT*,"The Method has converged. The update is within the tolerance set. "
        Converged = .TRUE.
    END IF


    IF ( mk == M .AND. .NOT. Converged ) THEN
        GVector = CSHIFT( GVector, SHIFT = +1, DIM = 2)
        FVector = CSHIFT( FVector, SHIFT = +1, DIM = 2)
    END IF



    IF ( Cur_Iteration == Max_Iterations ) THEN
        PRINT*,"FP_Accelerated has reached the maximum number of allowed iterations. "
        Converged = .TRUE.
    END IF


    !
    !   Calculate Source Vector with u_k
    !




    DO lm_loc = 1,LM_Length
        here = (lm_loc-1)*Num_R_Nodes + 1
        there = lm_loc*Num_R_Nodes
        FP_Coeff_Vector(:,lm_loc,ui) = GVectorM(Here:There)
    END DO




    IF ( Verbose_Flag ) THEN
        PRINT*,"In Anderson FP loop, Before calculating Source Vector."
    END IF
!    PRINT*,"Before Calc_FP_Source_Vector in loop"
    timer(1) = MPI_Wtime()
    CALL Calc_FP_Source_Vector()





!    PRINT*,"Before Check_FP_Convergence"
    timer(2) = MPI_Wtime()
    Call Check_FP_Convergence(Converged_Residual)
    timer(3) = MPI_WTime()
    IF ( Converged_Residual ) THEN
        PRINT*,"The Method has converged. The residual is within the tolerance set. "
        Converged = .TRUE.
    END IF

    CALL Clock_In(timer(2)-timer(1),6)
    CALL Clock_In(timer(3)-timer(2),7)
    Call Clock_In(timer(3)-timer(4),8)

    IF (Verbose_Flag .EQV. .TRUE. ) THEN
        CALL Print_Results()
        PRINT*," "
    END IF



!    DO ui= 3,3
!    DO lm = 1,LM_Length
!        PRINT*,ui,lm
!        PRINT*,FP_Coeff_Vector(:,lm,ui)
!    END DO
!    END DO



    IF ( Verbose_Flag .EQV. .TRUE. ) THEN
        WRITE(*,'(A,1X,I3.3,/)') "End of Iteration",Cur_Iteration
    END IF
!    STOP




END DO ! Converged Loop






CALL Deallocate_FP_Source_Variables()




END SUBROUTINE Fixed_Point_AndersonM
















!+101+###########################################################################!
!                                                                                !
!           Fixed_Point_Method                                                   !
!                                                                                !
!################################################################################!
SUBROUTINE Fixed_Point_AndersonM_B()

LOGICAL                                                 ::  CONVERGED
LOGICAL                                                 ::  CONVERGED_Residual
INTEGER                                                 ::  i, k, lm_loc
INTEGER                                                 ::  here, there
INTEGER                                                 ::  ui, re, d, l


INTEGER                                                 ::  M
INTEGER                                                 ::  LWORK
INTEGER                                                 ::  mk
INTEGER                                                 ::  INFO

REAL(KIND = idp), DIMENSION(1:4)                        :: timer

COMPLEX(idp), DIMENSION(1:NUM_R_NODES)                  ::  WORK_VEC
COMPLEX(idp), DIMENSION(1:NUM_R_NODES,1:NUM_R_NODES)    ::  WORK_MAT


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



!PRINT*,"In AndersonM"
!
!   Begin Method
!

IF (Verbose_Flag .EQV. .TRUE. ) THEN
    PRINT*,"Initial Guess"
    CALL Print_Results()
    PRINT*," "
END IF



!
!   Calculate Source Vector with u_0
!
IF ( Verbose_Flag ) THEN
    PRINT*,"Before Anderson FP loop, Before calculating Source Vector."
END IF



timer(1) = MPI_WTime()
CALL Calc_FP_Source_Vector()
timer(2) = MPI_Wtime()
Call Clock_In(timer(2)-timer(1),3)



ui = 1

DO WHILE ( .NOT. CONVERGED  .AND. Cur_Iteration < Max_Iterations)


    Timer(4) = MPI_WTime()

    Cur_Iteration = Cur_Iteration+1

    mk = MIN(Cur_Iteration, M)

    IF ( Cur_Iteration .NE. 1 ) THEN
        UVector = GVectorM
    ELSE
        DO ui = 1,5
        DO lm_loc = 1,LM_Length
!            PRINT*,FP_Coeff_Vector(:,lm_loc,ui)
            here = (ui-1)*Var_Dim + (lm_loc-1)*Num_R_Nodes + 1
            there = (ui-1)*Var_Dim + lm_loc*Num_R_Nodes
            UVector(here:there) = FP_Coeff_Vector(:,lm_loc,ui)
        END DO
        END DO
    END IF



    !
    !   Solve Systems
    !
    IF ( Verbose_Flag ) THEN
        PRINT*,"In Anderson FP loop, Before System Solves."
    END IF
    timer(1) = MPI_Wtime()
    IF ( ANY( CFA_EQ_Flags(1:2) == 1) ) THEN
        Call Solve_FP_System()
    END IF
    timer(2) = MPI_WTime()
    IF ( ANY( CFA_EQ_Flags(3:5) == 1) ) THEN
        Call Solve_FP_System_Beta()
    END IF
    timer(3) = MPI_Wtime()

    CALL Clock_In(timer(2)-timer(1),4)
    CALL Clock_In(timer(3)-timer(2),5)




    IF ( Verbose_Flag ) THEN
        PRINT*,"In Anderson FP loop, Before Update."
    END IF



    DO ui = 1,5
    DO lm_loc = 1,LM_Length
!            PRINT*,FP_Coeff_Vector(:,lm_loc,ui)
        here = (ui-1)*Var_Dim + (lm_loc-1)*Num_R_Nodes + 1
        there = (ui-1)*Var_Dim + lm_loc*Num_R_Nodes
        GVector(here:there,mk) = FP_Coeff_Vector(:,lm_loc,ui)
    END DO
    END DO
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


    !
    !   Calculate Source Vector with u_k
    !


    DO ui = 1,5
    DO lm_loc = 1,LM_Length
!            PRINT*,FP_Coeff_Vector(:,lm_loc,ui)
        here = (ui-1)*Var_Dim + (lm_loc-1)*Num_R_Nodes + 1
        there = (ui-1)*Var_Dim + lm_loc*Num_R_Nodes
        FP_Coeff_Vector(:,lm_loc,ui) = GVectorM(here:there)
    END DO
    END DO




    IF ( Verbose_Flag ) THEN
        PRINT*,"In Anderson FP loop, Before calculating Source Vector."
    END IF
!    PRINT*,"Before Calc_FP_Source_Vector in loop"
    timer(1) = MPI_Wtime()
    CALL Calc_FP_Source_Vector()





!    PRINT*,"Before Check_FP_Convergence"
!    timer(2) = MPI_Wtime()
!    Call Check_FP_Convergence(Converged_Residual)
!    timer(3) = MPI_WTime()
!    IF ( Converged_Residual ) THEN
!        PRINT*,"The Method has converged. The residual is within the tolerance set. "
!        Converged = .TRUE.
!    END IF

    CALL Clock_In(timer(2)-timer(1),6)
    CALL Clock_In(timer(3)-timer(2),7)
    Call Clock_In(timer(3)-timer(4),8)


!    PRINT*,"Before PRint_Results"
    IF (Verbose_Flag .EQV. .TRUE. ) THEN
        CALL Print_Results()
        PRINT*," "
    END IF



!    DO ui= 3,3
!    DO lm = 1,LM_Length
!        PRINT*,ui,lm
!        PRINT*,FP_Coeff_Vector(:,lm,ui)
!    END DO
!    END DO



    IF ( Verbose_Flag .EQV. .TRUE. ) THEN
        WRITE(*,'(A,1X,I3.3,/)') "End of Iteration",Cur_Iteration
    END IF
!    STOP




END DO ! Converged Loop

CALL Deallocate_FP_Source_Variables()




END SUBROUTINE Fixed_Point_AndersonM_B





!!+101+###########################################################################!
!!                                                                                !
!!           Fixed_Point_Method                                                   !
!!                                                                                !
!!################################################################################!
!SUBROUTINE Fixed_Point_AndersonM_B()
!
!LOGICAL                                                 ::  CONVERGED
!LOGICAL                                                 ::  CONVERGED_Residual
!INTEGER                                                 ::  i, k, lm_loc
!INTEGER                                                 ::  here, there
!INTEGER                                                 ::  ui, re, d, l
!
!
!INTEGER                                                 ::  M
!INTEGER                                                 ::  LWORK
!INTEGER                                                 ::  mk
!INTEGER                                                 ::  INFO
!
!REAL(KIND = idp), DIMENSION(1:4)                        :: timer
!
!COMPLEX(idp), DIMENSION(1:NUM_R_NODES)                  ::  WORK_VEC
!COMPLEX(idp), DIMENSION(1:NUM_R_NODES,1:NUM_R_NODES)    ::  WORK_MAT
!
!
!COMPLEX(idp),DIMENSION(:), ALLOCATABLE                  :: Resid_Vector
!COMPLEX(idp),DIMENSION(:,:), ALLOCATABLE                :: FVector
!COMPLEX(idp),DIMENSION(:,:), ALLOCATABLE                :: GVector
!COMPLEX(idp),DIMENSION(:), ALLOCATABLE                  :: BVector
!COMPLEX(idp),DIMENSION(:), ALLOCATABLE                  :: UVector
!COMPLEX(idp),DIMENSION(:), ALLOCATABLE                  :: GVectorM
!COMPLEX(idp),DIMENSION(:), ALLOCATABLE                  :: FVectorM
!COMPLEX(idp),DIMENSION(:,:), ALLOCATABLE                :: AMatrix
!COMPLEX(idp),DIMENSION(:), ALLOCATABLE                  :: Work
!COMPLEX(idp),DIMENSION(:), ALLOCATABLE                  :: Alpha
!
!
!M = FP_Anderson_M
!LWORK = 2*M
!
!ALLOCATE( Resid_Vector(1:Prob_Dim) )
!ALLOCATE( FVector(1:Prob_Dim,1:M) )
!ALLOCATE( GVectorM(1:Prob_Dim) )
!ALLOCATE( FVectorM(1:Prob_Dim) )
!ALLOCATE( BVector(1:Prob_Dim) )
!ALLOCATE( GVector(1:Prob_Dim,1:M) )
!ALLOCATE( Work(1:LWORK) )
!ALLOCATE( Alpha(1:M) )
!ALLOCATE( AMatrix(1:Prob_Dim,1:M))
!ALLOCATE( UVector(1:Prob_Dim) )
!
!CALL Allocate_FP_Source_Variables()
!
!timer = 0.0_idp
!CUR_ITERATION = 0
!CONVERGED = .FALSE.
!
!
!
!!PRINT*,"In AndersonM"
!!
!!   Begin Method
!!
!
!IF (Verbose_Flag .EQV. .TRUE. ) THEN
!    PRINT*,"Initial Guess"
!    CALL Print_Results()
!    PRINT*," "
!END IF
!
!
!
!!
!!   Calculate Source Vector with u_0
!!
!IF ( Verbose_Flag ) THEN
!    PRINT*,"Before Anderson FP loop, Before calculating Source Vector."
!END IF
!timer(1) = MPI_WTime()
!CALL Calc_FP_Source_Vector()
!timer(2) = MPI_Wtime()
!Call Clock_In(timer(2)-timer(1),3)
!
!
!ui = 1
!
!DO WHILE ( .NOT. CONVERGED  .AND. Cur_Iteration < Max_Iterations)
!
!
!
!    Timer(4) = MPI_WTime()
!
!    Cur_Iteration = Cur_Iteration+1
!
!    mk = MIN(Cur_Iteration, M)
!
!    IF ( Cur_Iteration .NE. 1 ) THEN
!        UVector = GVectorM
!    ELSE
!
!        UVector = FP_Coeff_Vector
!
!!        DO ui = 1,5
!!        DO lm_loc = 1,LM_Length
!!!            PRINT*,FP_Coeff_Vector(:,lm_loc,ui)
!!            here = (ui-1)*Var_Dim + (lm_loc-1)*Num_R_Nodes + 1
!!            there = (ui-1)*Var_Dim + lm_loc*Num_R_Nodes
!!            print*,ui,lm_loc,here,there
!!            UVector(here:there) = FP_Coeff_Vector(:,lm_loc,ui)
!!        END DO
!!        END DO
!    END IF
!
!
!
!    !
!    !   Solve Systems
!    !
!    IF ( Verbose_Flag ) THEN
!        PRINT*,"In Anderson FP loop, Before System Solves."
!    END IF
!    timer(1) = MPI_Wtime()
!    IF ( ANY( CFA_EQ_Flags(1:2) == 1) ) THEN
!        Call Solve_FP_System()
!    END IF
!    timer(2) = MPI_WTime()
!    IF ( ANY( CFA_EQ_Flags(3:5) == 1) ) THEN
!        Call Solve_FP_System_Beta()
!    END IF
!    timer(3) = MPI_Wtime()
!
!    CALL Clock_In(timer(2)-timer(1),4)
!    CALL Clock_In(timer(3)-timer(2),5)
!
!
!
!
!
!    IF ( Verbose_Flag ) THEN
!        PRINT*,"In Anderson FP loop, Before Update."
!    END IF
!
!
!    GVector(:,mk) = FP_Coeff_Vector(:)
!!    DO ui = 1,5
!!    DO lm_loc = 1,LM_Length
!!!            PRINT*,FP_Coeff_Vector(:,lm_loc,ui)
!!        here = (ui-1)*Var_Dim + (lm_loc-1)*Num_R_Nodes + 1
!!        there = (ui-1)*Var_Dim + lm_loc*Num_R_Nodes
!!        GVector(here:there,mk) = FP_Coeff_Vector(:,lm_loc,ui)
!!    END DO
!!    END DO
!    FVector(:,mk) = GVector(:,mk) - UVector(:)
!
!
!
!    IF ( mk == 1 ) THEN
!
!        GVectorM = GVector(:,mk)
!
!    ELSE
!!        PRINT*,"Here",mk,M
!
!        BVector = -FVector(:,mk)
!
!        AMatrix(:,1:mk-1) = FVector(:,1:mk-1) - SPREAD( FVector(:,mk), DIM=2, NCOPIES = mk-1)
!
!!        PRINT*,"Before ZGELS"
!        CALL ZGELS('N',Prob_Dim,mk-1,1,              &
!                    AMatrix(:,1:mk-1), Prob_Dim,     &
!                    BVector, Prob_Dim,                &
!                    WORK, LWORK, INFO )
!
!!        PRINT*,"After ZGELS"
!        Alpha(1:mk-1) = BVector(1:mk-1)
!        Alpha(mk)     = 1.0_idp - SUM( Alpha(1:mk-1) )
!
!
!        GVectorM = 0.0_idp
!        DO i = 1,mk
!
!            GVectorM = GVectorM + Alpha(i)*GVector(:,i)
!
!        END DO
!
!    END IF
!
!    FVectorM = GVectorM - UVector(:)
!
!
!    IF ( ALL( ABS( FVectorM ) <= Convergence_Criteria ) ) THEN
!        PRINT*,"The Method has converged. The absolute update is less than the tolerance set. "
!        Converged = .TRUE.
!    END IF
!!    IF ( ALL( ABS( FVectorM ) <= Convergence_Criteria*ABS( GVectorM ) ) ) THEN
!!        PRINT*,"The Method has converged. The relative update is less than the tolerance set. "
!!        Converged = .TRUE.
!!    END IF
!
!
!    IF ( mk == M .AND. .NOT. Converged ) THEN
!        GVector = CSHIFT( GVector, SHIFT = +1, DIM = 2)
!        FVector = CSHIFT( FVector, SHIFT = +1, DIM = 2)
!    END IF
!
!
!
!    IF ( Cur_Iteration == Max_Iterations ) THEN
!        PRINT*,"FP_Accelerated has reached the maximum number of allowed iterations. "
!        Converged = .TRUE.
!    END IF
!
!
!    !
!    !   Calculate Source Vector with u_k
!    !
!
!
!    FP_Coeff_Vector(:) = GVectorM(:)
!
!!    DO ui = 1,5
!!    DO lm_loc = 1,LM_Length
!!!            PRINT*,FP_Coeff_Vector(:,lm_loc,ui)
!!        here = (ui-1)*Var_Dim + (lm_loc-1)*Num_R_Nodes + 1
!!        there = (ui-1)*Var_Dim + lm_loc*Num_R_Nodes
!!        FP_Coeff_Vector(:,lm_loc,ui) = GVectorM(here:there)
!!    END DO
!!    END DO
!
!
!
!    IF ( Verbose_Flag ) THEN
!        PRINT*,"In Anderson FP loop, Before calculating Source Vector."
!    END IF
!!    PRINT*,"Before Calc_FP_Source_Vector in loop"
!    timer(1) = MPI_Wtime()
!    CALL Calc_FP_Source_Vector()
!
!
!
!
!
!!    PRINT*,"Before Check_FP_Convergence"
!!    timer(2) = MPI_Wtime()
!!    Call Check_FP_Convergence(Converged_Residual)
!!    timer(3) = MPI_WTime()
!!    IF ( Converged_Residual ) THEN
!!        PRINT*,"The Method has converged. The residual is within the tolerance set. "
!!        Converged = .TRUE.
!!    END IF
!
!    CALL Clock_In(timer(2)-timer(1),6)
!    CALL Clock_In(timer(3)-timer(2),7)
!    Call Clock_In(timer(3)-timer(4),8)
!
!    IF (Verbose_Flag .EQV. .TRUE. ) THEN
!        CALL Print_Results()
!        PRINT*," "
!    END IF
!
!
!
!!    DO ui= 3,3
!!    DO lm = 1,LM_Length
!!        PRINT*,ui,lm
!!        PRINT*,FP_Coeff_Vector(:,lm,ui)
!!    END DO
!!    END DO
!
!
!
!    IF ( Verbose_Flag .EQV. .TRUE. ) THEN
!        WRITE(*,'(A,1X,I3.3,/)') "End of Iteration",Cur_Iteration
!    END IF
!!    STOP
!
!
!
!
!END DO ! Converged Loop
!
!
!
!
!
!
!CALL Deallocate_FP_Source_Variables()
!
!
!
!
!END SUBROUTINE Fixed_Point_AndersonM_B


END MODULE FP_AndersonM_Module
