   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE FP_Method_Module                                                             !##!
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
!##!                                                                                !##!
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

USE Poseidon_Numbers_Module, &
            ONLY :  pi

USE Units_Module, &
            ONLY :  C_Square,                   &
                    Centimeter,                 &
                    Meter,                      &
                    Second,                     &
                    GravPot_Units,              &
                    Shift_Units

USE Variables_Functions, &
            ONLY :  Potential_Solution,         &
                    Calc_3D_Values_At_Location

USE Variables_Mesh, &
            ONLY :  Num_R_Elements,             &
                    R_Inner,                    &
                    R_Outer

USE Variables_Derived, &
            ONLY :  LM_Length,                  &
                    Num_R_Nodes,                &
                    Beta_Prob_Dim

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


USE DRIVER_Parameters, &
            ONLY :  Driver_Test_Number

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

USE Variables_MPI, &
            ONLY :  myID_Poseidon


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

USE Functions_Matrix, &
            ONLY :  MVMULT_FULL,                &
                    MVMULT_FULL_SUB,            &
                    MVMULT_CCS

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
            ONLY :  FP_LM_Map

USE FP_System_Solvers_Module, &
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
SUBROUTINE Fixed_Point_Method()

LOGICAL                                                 ::  CONVERGED
INTEGER                                                 ::  i

REAL(KIND = idp), DIMENSION(1:4)                        :: timer


CALL Allocate_FP_Source_Variables()

timer = 0.0_idp
CUR_ITERATION = 1
CONVERGED = .FALSE.

IF (myID_Poseidon == 0 ) THEN
    IF ( Verbose_Flag .eqv. .TRUE. ) THEN
        CALL OUTPUT_GUESS(0, 0)
    END IF
END IF


!PRINT*,"RHS_Terms(:,:,3) = 0.0 in SubJacobian_Functions_1D. Message in FP_Method."
!
!FP_Source_Vector = 0.0_idp
!FP_Coeff_Vector = 0.0_idp
!FP_Laplace_Vector = 0.0_idp
!PRINT*,"FP_Source_Vector, FP_Coeff_Vector, FP_Laplace_Vector Zeroed in Poseidon_FP_Method.f90"

!
!   Begin Method
!


timer(1) = MPI_WTime()
CALL Calc_FP_Source_Vector()
timer(2) = MPI_WTime()
Call Clock_In(timer(2)-timer(1),3)


DO WHILE ( CONVERGED .EQV. .FALSE. )
    IF (myID_Poseidon == 0 ) THEN
!        CALL OPEN_ITER_REPORT_FILE(Cur_Iteration, myID_Poseidon)
    END IF


    timer(4) = MPI_Wtime()
    IF ( ANY( CFA_EQ_Flags(1:2) == 1) ) THEN
        Call Solve_FP_System()
    END IF
    timer(2) = MPI_WTime()
    IF ( ANY( CFA_EQ_Flags(3:5) == 1) ) THEN
        Call Solve_FP_System_Beta()
    END IF
    timer(3) = MPI_Wtime()

    CALL Clock_In(timer(2)-timer(4),4)
    CALL Clock_In(timer(3)-timer(2),5)

    timer(1) = MPI_Wtime()
    CALL Calc_FP_Source_Vector()
    timer(2) = MPI_Wtime()
    Call Check_FP_Convergence(Converged)
    timer(3) = MPI_WTime()

    CALL Clock_In(timer(2)-timer(1),6)
    CALL Clock_In(timer(3)-timer(2),7)
    
    timer(2) = MPI_Wtime()
    Call Clock_In(timer(2)-timer(4),8)


!    PRINT*,"Before Outputs"
    IF ( myID_Poseidon == 0 ) THEN
        IF ( Verbose_Flag .eqv. .TRUE. ) THEN
            CALL OUTPUT_ITERATION_REPORT(Cur_Iteration, myID_Poseidon)
!            PRINT*,"After OUTPUT_ITERATION_REPORT"
    !        CALL CLOSE_ITER_REPORT_FILE()
!            PRINT*,"After CLOSE_ITER_REPORT_FILE"
            !CALL OUTPUT_COEFFICIENT_VECTOR_FORTRAN()
        END IF
    END IF
    
    IF ( Verbose_Flag .EQV. .TRUE. ) THEN
        WRITE(*,'(A,1X,I3.3,/)') "End of Iteration",Cur_Iteration
    END IF

    Cur_Iteration = Cur_Iteration + 1
END DO ! Converged Loop
!PRINT*,"Outside Fixed-Point Iteration loop."


IF ( WRITE_RESULTS_FLAG == 1 ) THEN
!    CALL OUTPUT_FINAL_RESULTS()
END IF



CALL Deallocate_FP_Source_Variables()




END SUBROUTINE Fixed_Point_Method







!+101+###########################################################################!
!                                                                                !
!           Fixed_Point_Method                                                   !
!                                                                                !
!################################################################################!
SUBROUTINE Fixed_Point_Accelerated()

LOGICAL                                                 ::  CONVERGED
LOGICAL                                                 ::  CONVERGED_blank
INTEGER                                                 ::  i, k

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

ALLOCATE( Resid_Vector(1:Num_R_Nodes) )
ALLOCATE( FVector(1:Num_R_Nodes,1:M) )
ALLOCATE( GVectorM(1:Num_R_Nodes) )
ALLOCATE( FVectorM(1:Num_R_Nodes) )
ALLOCATE( BVector(1:Num_R_Nodes) )
ALLOCATE( GVector(1:Num_R_Nodes,1:M) )
ALLOCATE( Work(1:LWORK) )
ALLOCATE( Alpha(1:M) )
ALLOCATE( AMatrix(1:Num_R_Nodes,1:M))
ALLOCATE( UVector(1:Num_R_Nodes) )

CALL Allocate_FP_Source_Variables()

timer = 0.0_idp
CUR_ITERATION = 0
CONVERGED = .FALSE.




!
!   Begin Method
!


!
!   Calculate Source Vector with u_0
!
timer(1) = MPI_WTime()
CALL Calc_FP_Source_Vector()
timer(2) = MPI_Wtime()
Call Clock_In(timer(2)-timer(1),3)


DO WHILE ( .NOT. CONVERGED  .AND. Cur_Iteration < Max_Iterations)
    timer(4) = MPI_Wtime()
    Cur_Iteration = Cur_Iteration + 1

    mk = MIN(Cur_Iteration, M)

    UVector = FP_Coeff_Vector(:,0,CFA_EQ_Map(1))
    PRINT*,"UVector"
    PRINT*,UVector

    !
    !   Solve Systems
    !
    timer(1) = MPI_WTime()
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


    GVector(:,mk) = FP_Coeff_Vector(:,0,CFA_EQ_Map(1))
    FVector(:,mk) = GVector(:,mk) - UVector(:)
    PRINT*,"FVector"
    PRINT*,FVector

    IF ( mk == 1 ) THEN

        GVectorM = GVector(:,mk)

    ELSE

        BVector = -FVector(:,mk)

        AMatrix(:,1:mk-1) = FVector(:,1:mk-1) - SPREAD( FVector(:,mk), DIM=2, NCOPIES = mk-1)

        CALL ZGELS('N',Num_R_Nodes,mk-1,1,              &
                    AMatrix(:,1:mk-1), Num_R_Nodes,     &
                    BVector,Num_R_Nodes,                &
                    WORK, LWORK, INFO )


        Alpha(1:mk-1) = BVector(1:mk-1)
        Alpha(mk)     = 1.0_idp - SUM( Alpha(1:mk-1) )


        GVectorM = 0.0_idp
        DO i = 1,mk

            GVectorM = GVectorM + Alpha(i)*GVector(:,i)
    
        END DO

    END IF

    FVectorM = GVectorM - UVector(:)
    PRINT*,"GVectorM"
    PRINT*,GVectorM

    FP_Coeff_Vector(:,0,CFA_EQ_Map(1)) = GVectorM







    IF ( ALL( ABS( FVectorM ) <= Convergence_Criteria*ABS(FP_Coeff_Vector(:,0,CFA_EQ_Map(1))) ) ) THEN
        PRINT*,"The Method has converged."
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


    IF ( Verbose_Flag .EQV. .TRUE. ) THEN
        WRITE(*,'(A,1X,I3.3,/)') "End of Iteration",Cur_Iteration
    END IF
    

    !
    !   Calculate Source Vector with u_k
    !
    timer(1) = MPI_Wtime()
    CALL Calc_FP_Source_Vector()
    timer(2) = MPI_Wtime()
    Call Check_FP_Convergence(Converged_Blank)
    timer(3) = MPI_WTime()

    CALL Clock_In(timer(2)-timer(1),6)
    CALL Clock_In(timer(3)-timer(2),7)
    Call Clock_In(timer(3)-timer(4),8)


!    STOP
END DO ! Converged Loop



CALL Deallocate_FP_Source_Variables()




END SUBROUTINE Fixed_Point_Accelerated








!+501+##########################################################################!
!                                                                               !
!                   OUTPUT_ITERATION_REPORT                                     !
!                                                                               !
!###############################################################################!
SUBROUTINE OUTPUT_ITERATION_REPORT(Iter, Rank)

INTEGER, INTENT(IN)                 :: Iter, Rank

INTEGER, DIMENSION(0:1)                         ::  FILE_ID
INTEGER                                         ::  i, j
REAL(KIND = idp)                                ::  r, theta, phi, deltar
REAL(KIND = idp)                                ::  Analytic_Val, Solver_Val
REAL(KIND = idp)                                ::  Return_Psi, Return_AlphaPsi
REAL(KIND = idp)                                ::  Return_Beta1, Return_Beta2, Return_Beta3
REAL(KIND = idp)                                ::  PsiPot_Val, AlphaPsiPot_Val


REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Output_re,          &
                                                                Output_rc,          &
                                                                Output_dr

120 FORMAT (A61)
121 FORMAT (A1)
122 FORMAT (A41,I2.2)
123 FORMAT (A38,ES22.15)

109 FORMAT (A,I2.2,A,I2.2)
110 FORMAT (11X,A1,24X,A3,19X,A8,15X,A11,14X,A11,14X,A11)
111 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15)
112 FORMAT (A43,I2.2,A2,I2.2,A4)





FILE_ID = -1
IF ( FRAME_REPORT_FLAG == 1 ) THEN

    FILE_ID(0) = FRAME_REPORT_FILE_ID

END IF

IF (( WRITE_REPORT_FLAG == 2) .OR. (WRITE_REPORT_FLAG == 3) ) THEN
    FILE_ID(1) = ITER_REPORT_FILE_ID
END IF



! Write Title to File
DO i = 0,1

    IF ( FILE_ID(i) .NE. -1 ) THEN

        WRITE(FILE_ID(i),'(A)')"                                     Timing Results"
        WRITE(FILE_ID(i),'(A)')"            ============================================================="
        WRITE(FILE_ID(i),'(A)')" "
        WRITE(FILE_ID(i),123)"                    Initialize Time : ",ITER_TIME_TABLE(1)
        WRITE(FILE_ID(i),123)" Input/Communicate Source Data Time : ",ITER_TIME_TABLE(2)
        WRITE(FILE_ID(i),123)"     Input Boundary Conditions Time : ",ITER_TIME_TABLE(3)
        WRITE(FILE_ID(i),123)"        CFA_3D_Apply_BCs_Part1 Time : ",ITER_TIME_TABLE(4)
        WRITE(FILE_ID(i),120)"-------------------------------------------------------------"
        WRITE(FILE_ID(i),123)" ||     Calc_3D_Current_Values Time : ",ITER_TIME_TABLE(5)
        WRITE(FILE_ID(i),123)" ||    CREATE_3D_SubJcbn_Terms Time : ",ITER_TIME_TABLE(6)
        WRITE(FILE_ID(i),123)" ||       CREATE_3D_RHS_VECTOR Time : ",ITER_TIME_TABLE(7)
        WRITE(FILE_ID(i),123)"\  /     CREATE_3D_JCBN_MATRIX Time : ",ITER_TIME_TABLE(8)
        WRITE(FILE_ID(i),120)"-\/ ---------------------------------------------------------"
        WRITE(FILE_ID(i),123)"CREATE_3D_NONLAPLACIAN_STF_MAT Time : ",ITER_TIME_TABLE(9)
        WRITE(FILE_ID(i),123)"REDUCE_3D_NONLAPLACIAN_STF_MAT Time : ",ITER_TIME_TABLE(10)
        WRITE(FILE_ID(i),123)"FINISH_3D_NONLAPLACIAN_STF_MAT Time : ",ITER_TIME_TABLE(11)
        WRITE(FILE_ID(i),123)"          FINISH_3D_RHS_VECTOR Time : ",ITER_TIME_TABLE(12)
        WRITE(FILE_ID(i),123)"        CFA_3D_Apply_BCs_Part2 Time : ",ITER_TIME_TABLE(13)
        WRITE(FILE_ID(i),123)"                    CFA_Solver Time : ",ITER_TIME_TABLE(14)
        WRITE(FILE_ID(i),123)"        CFA_Coefficient_Update Time : ",ITER_TIME_TABLE(15)
        WRITE(FILE_ID(i),123)"   CFA_Coefficient_Share_PETSc Time : ",ITER_TIME_TABLE(16)
        WRITE(FILE_ID(i),123)"         CFA_Convergence_Check Time : ",ITER_TIME_TABLE(17)
        WRITE(FILE_ID(i),123)"               Total Iteration Time : ",ITER_TIME_TABLE(18)
        WRITE(FILE_ID(i),123)"             Poseidon_Dist_Sol Time : ",ITER_TIME_TABLE(19)
        WRITE(FILE_ID(i),120)"============================================================="
        WRITE(FILE_ID(i),121)" "
        WRITE(FILE_ID(i),121)" "
        WRITE(FILE_ID(i),121)" "
        WRITE(FILE_ID(i),121)" "

        WRITE(FILE_ID(i),'(A,I2.2,A)')"                        Iteration ",Iter," Results"
        WRITE(FILE_ID(i),'(A)')"            ============================================================="
        WRITE(FILE_ID(i),'(A)')" "
        WRITE(FILE_ID(i),110)"r","Analytic Potential","Psi Potential","AlphaPsi Potential","Beta Value1","Beta Value2","Beta Value3"

    END IF
END DO






IF ( Driver_Test_Number .NE. -1 ) THEN

    ALLOCATE( Output_re(0:ITER_REPORT_NUM_SAMPLES) )
    ALLOCATE( Output_rc(1:ITER_REPORT_NUM_SAMPLES) )
    ALLOCATE( Output_dr(1:ITER_REPORT_NUM_SAMPLES) )

    ! Write Results Table Header to Screen
    IF (( WRITE_REPORT_FLAG == 1) .OR. (WRITE_REPORT_FLAG == 3) ) THEN
        IF ( Rank == 0 ) THEN
            PRINT*,"+++++++++++++++++++ myID,",Rank," Iteration",Iter,"++++++++++++++++++++++"
            WRITE(*,110)"r","Analytic Potential","Psi Potential","AlphaPsi Potential","Beta Value1","Beta Value2","Beta Value3"
        END IF
    END IF


    
    IF ( R_OUTER/(R_Inner+1.0_idp) > 1E3 ) THEN
        CALL Create_Logarithmic_1D_Mesh( R_INNER, R_OUTER, ITER_REPORT_NUM_SAMPLES,     &
                                        output_re, output_rc, output_dr     )
    ELSE
        CALL Create_Uniform_1D_Mesh( R_INNER, R_OUTER, ITER_REPORT_NUM_SAMPLES,     &
                                     output_re, output_rc, output_dr     )
    END IF


    DO i = 0,ITER_REPORT_NUM_SAMPLES

        r = output_re(i)
        theta = pi/2.0_idp
        phi = pi/2.0_idp

        CALL Calc_3D_Values_At_Location( r, theta, phi,                              &
                                        Return_Psi, Return_AlphaPsi,                &
                                        Return_Beta1, Return_Beta2, Return_Beta3    )

        ! Determine the Newtonian Potential at the location r, theta, phi
        Analytic_Val = Potential_Solution(r,theta,phi)/GravPot_Units


        ! AlphaPsi_to_Pot   =   2*C_Square*(AlphaPsi - 1)
        ! Psi_to_Pot        =   2*C_Square*(1 - Psi)

        ! Calculate Conformal Factor value from Newtonian Potential
        PsiPot_Val = 2.0_idp*C_Square*(1.0_idp - Return_Psi)/GravPot_Units

        ! Calculate the product of the Conformal Factor and Lapse Function from Newtonian Potential
        AlphaPsiPot_Val = 2.0_idp*C_Square*(Return_AlphaPsi - 1.0_idp)/GravPot_Units


        ! Write Results to Screen
        IF (( WRITE_REPORT_FLAG == 1) .OR. (WRITE_REPORT_FLAG == 3) ) THEN
            IF ( Rank == 0 ) THEN
                WRITE(*,111) r/Centimeter,              &
                             Analytic_Val,              &
                             PsiPot_Val,                &
                             AlphaPsiPot_Val,           &
                             Return_Beta1/Shift_Units,  &
                             Return_Beta2,              &
                             Return_Beta3

                WRITE(*,111) r/Centimeter,              &
                             Analytic_Val,              &
                             PsiPot_Val,                &
                             AlphaPsiPot_Val,           &
                             Return_Beta1/Shift_Units,  &
                             Return_Beta2,              &
                             Return_Beta3
            END IF
        END IF

        ! Write Results to File
        DO j = 0,1
            IF ( FILE_ID(j) .NE. -1 ) THEN
                WRITE(FILE_ID(j),111) r/Centimeter,              &
                                      Analytic_Val,              &
                                      PsiPot_Val,                &
                                      AlphaPsiPot_Val,           &
                                      Return_Beta1/Shift_Units,  &
                                      Return_Beta2,              &
                                      Return_Beta3
            END IF
        END DO ! j Loop

    END DO  ! i Loop







ELSE

    ALLOCATE( Output_re(0:ITER_REPORT_NUM_SAMPLES) )
    ALLOCATE( Output_rc(1:ITER_REPORT_NUM_SAMPLES) )
    ALLOCATE( Output_dr(1:ITER_REPORT_NUM_SAMPLES) )

    ! Write Results Table Header to Screen
    IF (( WRITE_REPORT_FLAG == 1) .OR. (WRITE_REPORT_FLAG == 3) ) THEN
        IF ( Rank == 0 ) THEN
            PRINT*,"+++++++++++++++++++ myID,",Rank," Iteration",Iter,"++++++++++++++++++++++"
!            WRITE(*,110)"r","Psi Potential","AlphaPsi Potential","Beta Value1","Beta Value2","Beta Value3"
            WRITE(*,110)"r","Psi","AlphaPsi","Beta Value1","Beta Value2","Beta Value3"
        END IF
    END IF

    IF ( R_OUTER/(R_Inner+1.0_idp) > 1E3 ) THEN
        CALL Create_Logarithmic_1D_Mesh( R_INNER, R_OUTER, ITER_REPORT_NUM_SAMPLES,     &
                                        output_re, output_rc, output_dr     )
    ELSE
        CALL Create_Uniform_1D_Mesh( R_INNER, R_OUTER, ITER_REPORT_NUM_SAMPLES,     &
                                     output_re, output_rc, output_dr     )
    END IF

    DO i = 0,ITER_REPORT_NUM_SAMPLES

        r = output_re(i)
        theta = pi/2.0_idp
        phi = pi/2.0_idp


        CALL Calc_3D_Values_At_Location( r, theta, phi,                              &
                                        Return_Psi, Return_AlphaPsi,                &
                                        Return_Beta1, Return_Beta2, Return_Beta3    )

        ! AlphaPsi_to_Pot   =   2*C_Square*(AlphaPsi - 1)
        ! Psi_to_Pot        =   2*C_Square*(1 - Psi)

        ! Calculate Conformal Factor value from Newtonian Potential
        PsiPot_Val = 2.0_idp*C_Square*(1.0_idp - Return_Psi)/GravPot_Units

        ! Calculate the product of the Conformal Factor and Lapse Function from Newtonian Potential
        AlphaPsiPot_Val = 2.0_idp*C_Square*(Return_AlphaPsi - 1.0_idp)/GravPot_Units


        ! Write Results to Screen
        IF (( WRITE_REPORT_FLAG == 1) .OR. (WRITE_REPORT_FLAG == 3) ) THEN
            IF ( Rank == 0 ) THEN
!                WRITE(*,111) r/Centimeter,              &
!                             PsiPot_Val,                &
!                             AlphaPsiPot_Val,           &
!                             Return_Beta1/Shift_Units,  &
!                             Return_Beta2,              &
!                             Return_Beta3


                WRITE(*,111)  r/Centimeter,              &
                              Return_Psi,                &
                              Return_AlphaPsi,           &
                              Return_Beta1/Shift_Units,  &
                              Return_Beta2,              &
                              Return_Beta3
            END IF
        END IF

        ! Write Results to File
        DO j = 0,1
            IF ( FILE_ID(j) .NE. -1 ) THEN
                WRITE(FILE_ID(j),111) r/Centimeter,              &
                                      PsiPot_Val,                &
                                      AlphaPsiPot_Val,           &
                                      Return_Beta1/Shift_Units,  &
                                      Return_Beta2,              &
                                      Return_Beta3

!                WRITE(FILE_ID(j),111) r/Centimeter,              &
!                                      Return_Psi,                &
!                                      Return_AlphaPsi,           &
!                                      Return_Beta1/Shift_Units,  &
!                                      Return_Beta2,              &
!                                      Return_Beta3
            END IF
        END DO ! j Loop

    END DO  ! i Loop




END IF



!WRITE( FRAME_REPORT_FILE_ID, '(4/)')

END SUBROUTINE OUTPUT_ITERATION_REPORT


!+501+##########################################################################!
!                                                                               !
!                   OUTPUT_GUESS                                                !
!                                                                               !
!###############################################################################!
SUBROUTINE OUTPUT_GUESS(Iter, Rank)

INTEGER, INTENT(IN)                 :: Iter, Rank

INTEGER                                         ::  FILE_ID
INTEGER                                         ::  i
REAL(KIND = idp)                                ::  r, theta, phi, deltar
REAL(KIND = idp)                                ::  Analytic_Val, Solver_Val
REAL(KIND = idp)                                ::  Return_Psi, Return_AlphaPsi
REAL(KIND = idp)                                ::  Return_Beta1, Return_Beta2, Return_Beta3
REAL(KIND = idp)                                ::  PsiPot_Val, AlphaPsiPot_Val

REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Output_re,          &
                                                                Output_rc,          &
                                                                Output_dr


120 FORMAT (A61)
121 FORMAT (A1)
122 FORMAT (A41,I2.2)
123 FORMAT (A38,ES22.15)

109 FORMAT (A,I2.2,A,I2.2)
110 FORMAT (11X,A1,16X,A18,9X,A13,10X,A18,10X,A11,14X,A11,14X,A11)
111 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15)
112 FORMAT (A43,I2.2,A2,I2.2,A4)
113 FORMAT (11X,A1,14X,A13,10X,A18,14X,A11,14X,A11,14X,A11)


IF ( .FALSE. ) THEN
    ! Write Results Table Header to Screen
    IF (( WRITE_REPORT_FLAG == 1) .OR. (WRITE_REPORT_FLAG == 3) ) THEN
        IF ( Rank == 0 ) THEN
            PRINT*,"+++++++++++++++++++ myID,",Rank," Iteration",Iter,"++++++++++++++++++++++"
            WRITE(*,110)"r (cm)","Analytic Potential","Psi Potential","AlphaPsi Potential","Beta Value1","Beta Value2","Beta Value3"
        END IF
    END IF

    deltar = ( R_OUTER - R_INNER )/ REAL(ITER_REPORT_NUM_SAMPLES, KIND = idp)
    DO i = 0,ITER_REPORT_NUM_SAMPLES

        r = i*deltar + R_INNER
        theta = pi/2.0_idp
        phi = pi/2.0_idp


        CALL Calc_3D_Values_At_Location( r, theta, phi,                              &
                                        Return_Psi, Return_AlphaPsi,                &
                                        Return_Beta1, Return_Beta2, Return_Beta3    )

        ! Determine the Newtonian Potential at the location r, theta, phi
        Analytic_Val = Potential_Solution(r,theta,phi)/GravPot_Units


        ! AlphaPsi_to_Pot   =   2*C_Square*(AlphaPsi - 1)
        ! Psi_to_Pot        =   2*C_Square*(1 - Psi)

        ! Calculate Conformal Factor value from Newtonian Potential
        PsiPot_Val = 2.0_idp*C_Square*(1.0_idp - Return_Psi)/GravPot_Units

        ! Calculate the product of the Conformal Factor and Lapse Function from Newtonian Potential
        AlphaPsiPot_Val = 2.0_idp*C_Square*(Return_AlphaPsi - 1.0_idp)/GravPot_Units


        ! Write Results to Screen
        IF (( WRITE_REPORT_FLAG == 1) .OR. (WRITE_REPORT_FLAG == 3) ) THEN
            IF ( Rank == 0 ) THEN
                WRITE(*,111) r/Centimeter,              &
                              Analytic_Val,              &
                              PsiPot_Val,                &
                              AlphaPsiPot_Val,           &
                              Return_Beta1/Shift_Units,  &
                              Return_Beta2,              &
                              Return_Beta3
            END IF
        END IF


    END DO








ELSE

    ALLOCATE( Output_re(0:ITER_REPORT_NUM_SAMPLES) )
    ALLOCATE( Output_rc(1:ITER_REPORT_NUM_SAMPLES) )
    ALLOCATE( Output_dr(1:ITER_REPORT_NUM_SAMPLES) )

    ! Write Results Table Header to Screen
    IF (( WRITE_REPORT_FLAG == 1) .OR. (WRITE_REPORT_FLAG == 3) ) THEN
        IF ( Rank == 0 ) THEN
            PRINT*,"+++++++++++++++++++ myID,",Rank," Iteration",Iter,"++++++++++++++++++++++"
            WRITE(*,113)"r (cm)","Psi ","AlphaPsi ","Beta Value1","Beta Value2","Beta Value3"
        END IF
    END IF


    IF ( R_OUTER/(R_Inner+1.0_idp) > 1E3 ) THEN
        CALL Create_Logarithmic_1D_Mesh( R_INNER, R_OUTER, ITER_REPORT_NUM_SAMPLES,     &
                                        output_re, output_rc, output_dr     )
    ELSE
        CALL Create_Uniform_1D_Mesh( R_INNER, R_OUTER, ITER_REPORT_NUM_SAMPLES,     &
                                     output_re, output_rc, output_dr     )
    END IF


    DO i = 0,ITER_REPORT_NUM_SAMPLES

        r = output_re(i)
        theta = pi/2.0_idp
        phi = pi/2.0_idp


        CALL Calc_3D_Values_At_Location( r, theta, phi,                              &
                                        Return_Psi, Return_AlphaPsi,                &
                                        Return_Beta1, Return_Beta2, Return_Beta3    )



        ! Write Results to Screen
        IF (( WRITE_REPORT_FLAG == 1) .OR. (WRITE_REPORT_FLAG == 3) ) THEN
            IF ( Rank == 0 ) THEN
                WRITE(*,111) r/Centimeter,              &
                              Return_Psi,                &
                              Return_AlphaPsi,           &
                              Return_Beta1/Shift_Units,  &
                              Return_Beta2,              &
                              Return_Beta3
            END IF
        END IF


    END DO


END IF



END SUBROUTINE OUTPUT_GUESS










!+501+##########################################################################!
!                                                                               !
!                   Laplace_Test_Residual                                       !
!                                                                               !
!###############################################################################!
SUBROUTINE Laplace_Test_Residual()

REAL(KIND = idp)                                ::  deltar, r
REAL(KIND = idp)                                ::  theta, phi
REAL(KIND = idp)                                ::  C_One, C_Two

REAL(KIND = idp)                                ::  f
REAL(KIND = idp)                                ::  Return_Psi, Return_AlphaPsi
REAL(KIND = idp)                                ::  Return_Beta1, Return_Beta2, Return_Beta3

INTEGER                                         ::  Num_Samples,i


!C_One = -9.99999999999999E1_idp
!C_Two = -5.00000000000000E1_idp

C_One = -4.705882352941180E+01_idp
C_Two = -2.352941176470590E+01_idp

Num_Samples = 30


deltar = ( R_OUTER - R_INNER )/ REAL(Num_Samples, KIND = idp)

DO i = 0,Num_Samples

    r = i*deltar + R_INNER
    theta = pi/2.0_idp
    phi = pi/2.0_idp

    f = C_One*r + C_Two*(r)**(-2)

    CALL Calc_3D_Values_At_Location( r, theta, phi,                              &
                                  Return_Psi, Return_AlphaPsi,                &
                                  Return_Beta1, Return_Beta2, Return_Beta3    )

!    WRITE(*,*)r,Return_Beta1,f, abs(Return_Beta1-f), abs(Return_Beta1-f)/abs(f)

END DO  ! i Loop



END SUBROUTINE Laplace_Test_Residual










END MODULE FP_Method_Module
