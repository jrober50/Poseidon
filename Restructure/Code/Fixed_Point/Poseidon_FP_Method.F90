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
                    Num_R_Nodes

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
                    CONVERGENCE_CRITERIA,       &
                    CONVERGENCE_FLAG,           &
                    NUM_CFA_EQs

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
                    Laplace_Matrix_Full,        &
                    Laplace_Matrix_VAL,         &
                    Laplace_Matrix_ROW,         &
                    Laplace_Matrix_COL,         &
                    Laplace_Factored_VAL,       &
                    Laplace_Factored_ROW,       &
                    Laplace_Factored_COL,       &
                    CFA_Eq_Map,                 &
                    CFA_MAT_Map,                &
                    Laplace_NNZ,                &
                    Num_Matrices,               &
                    MCF_Flag

USE Functions_Matrix, &
            ONLY :  MVMULT_FULL,                &
                    MVMULT_CCS

USE FP_Source_Vector_Module, &
            ONLY :  Calc_FP_Source_Vector,          &
                    Allocate_FP_Source_Variables,   &
                    Deallocate_FP_Source_Variables

USE FP_Functions_BC,  &
            ONLY :  Dirichlet_BC,                   &
                    Dirichlet_BC_CCS,               &
                    Dirichlet_BC_CHOL,              &
                    Neumann_BC,                     &
                    Neumann_BC_CCS

USE Poseidon_Cholesky_Module,   &
            ONLY :  Cholesky_Factorization,         &
                    CCS_Back_Substitution,          &
                    CCS_Forward_Substitution

USE Linear_Solvers_And_Preconditioners, &
            ONLY :  PRECOND_CONJ_GRAD_CCS,           &
                    JACOBI_CONDITIONING

USE Poseidon_IO_Module, &
            ONLY :  Clock_In,                           &
                    OPEN_ITER_REPORT_FILE,              &
                    CLOSE_ITER_REPORT_FILE,             &
                    OUTPUT_FINAL_RESULTS

USE FP_Functions_Mapping, &
            ONLY :  FP_LM_Map


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

REAL(KIND = idp), DIMENSION(1:3)                        :: timer


CALL Allocate_FP_Source_Variables()

timer = 0.0_idp
CUR_ITERATION = 1
CONVERGED = .FALSE.

IF (myID_Poseidon == 0 ) THEN
    CALL OUTPUT_GUESS(0, 0)
END IF

PRINT*,"RHS_Terms(:,:,3) = 0.0 in SubJacobian_Functions_1D. Message in FP_Method."

!
!   Begin Method
!

!PRINT*,"Before First Calc_FP_Source_Vector"
CALL Calc_FP_Source_Vector()
!
!PRINT*,"Stop after first Calc_FP_Source_Vector()"
!STOP

DO WHILE ( CONVERGED .EQV. .FALSE. )
    IF (myID_Poseidon == 0 ) THEN
        CALL OPEN_ITER_REPORT_FILE(Cur_Iteration, myID_Poseidon)
    END IF


!    PRINT*,"Starting FP Iteration ",Cur_Iteration
    timer(1) = MPI_Wtime()
    
!    PRINT*,"Before Solve_FP_System"
    Call Solve_FP_System()
    
!    PRINT*,"Before Inner Calc_FP_Source_Vector"
    CALL Calc_FP_Source_Vector()

!    PRINT*,"Before Check_FP_Convergence"
    Call Check_FP_Convergence(Converged)

    IF ( myID_Poseidon == 0 ) THEN
        CALL OUTPUT_ITERATION_REPORT(Cur_Iteration, myID_Poseidon)
        CALL CLOSE_ITER_REPORT_FILE()
        !CALL OUTPUT_COEFFICIENT_VECTOR_FORTRAN()
    END IF



    Cur_Iteration = Cur_Iteration + 1
END DO ! Converged Loop




IF ( WRITE_RESULTS_FLAG == 1 ) THEN
    CALL OUTPUT_FINAL_RESULTS()
END IF



CALL Deallocate_FP_Source_Variables()

END SUBROUTINE Fixed_Point_Method







!+201+###########################################################################!
!                                                                                !
!           Check_FP_Convergence                                                 !
!                                                                                !
!################################################################################!
SUBROUTINE Check_FP_Convergence( CONVERGED )

LOGICAL, INTENT(INOUT)                              ::  CONVERGED

REAL(KIND = idp)                                    ::  Convergence_Stat


INTEGER                                             ::  ui, l


! Test Residual
CALL Calc_FP_Residual( Convergence_Stat )

Frame_Update_Table(CUR_ITERATION)   =    MAXVAL(ABS(FP_Update_Vector))
Frame_Residual_Table(CUR_ITERATION) =    Convergence_Stat



PRINT*,"Convergence Check, Iter ",CUR_ITERATION," - Residual : ",Convergence_Stat," Criteria : ",CONVERGENCE_CRITERIA
PRINT*,"MAX(ABS(FP_Update_Vector))",MAXVAL(ABS(FP_Update_Vector))
IF ( Convergence_Stat < CONVERGENCE_CRITERIA) THEN
    CONVERGED = .TRUE.
    CONVERGENCE_FLAG = 1
    PRINT*,"Solver has converged. Residual ",Convergence_Stat," < ",CONVERGENCE_CRITERIA

ELSEIF ( MAXVAL(ABS(FP_Update_Vector)) == 0.0_idp ) THEN
    CONVERGED = .TRUE.
    CONVERGENCE_FLAG = 1
    PRINT*,"Solver has converged. Fixed Point Solution change is zero. "


!
!   Has the Solver taken the max number of iterations allowed?
!
ELSEIF ( CUR_ITERATION > MAX_ITERATIONS-1 ) THEN
    CONVERGED = .TRUE.
    CONVERGENCE_FLAG = 2
    PRINT*,"Solver has not converged. Iterations has exceeded maximum, ", Max_Iterations

END IF


END SUBROUTINE Check_FP_Convergence









!+202+##########################################################################!
!                                                                               !
!           Calc_FP_Residual                                                    !
!                                                                               !
!###############################################################################!
SUBROUTINE Calc_FP_Residual( Convergence_Stat )

REAL(KIND = idp), INTENT(OUT)                       ::  Convergence_Stat


INTEGER                                             ::  Convergence_Type = 1

INTEGER                                             ::  ui, l, m, lm_loc, map_loc

COMPLEX(KIND = idp), DIMENSION(1:NUM_R_NODES)                               ::  WORK_VEC
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:)                            ::  WORK_MAT

! Multi Matrices and Coeff Vectors to form Laplace
IF ( Matrix_Format == 'Full' ) THEN

    ALLOCATE (WORK_MAT(1:NUM_R_NODES, 1:NUM_R_NODES))
    

    DO ui = 1,NUM_CFA_EQs
        DO l = 0,L_Limit
            DO m = -l,l
                lm_loc = FP_LM_Map(l,m)
                map_loc = CFA_MAT_Map(CFA_EQ_Map(ui))

                WORK_MAT = Laplace_Matrix_Full(:,:,l,map_loc)
                WORK_VEC = FP_Source_Vector(:,lm_loc,ui)
    
                CALL DIRICHLET_BC(WORK_MAT, WORK_VEC, l, m, ui)

                FP_Laplace_Vector(:,lm_loc,ui) = MVMULT_FULL( WORK_MAT,             &
                                                        FP_Coeff_Vector(:,lm_loc,CFA_EQ_Map(ui)), &
                                                        NUM_R_NODES, NUM_R_NODES                   )


                FP_Residual_Vector(:,lm_loc,ui) = FP_Laplace_Vector(:,lm_loc,ui) - Work_Vec
                



!                PRINT*,"Laplace", lm_loc, ui
!                PRINT*,FP_Laplace_Vector
!                PRINT*,"Matrix"
!                PRINT*,WORK_MAT
!                PRINT*,"Vector"
!                PRINT*,FP_Coeff_Vector(:,lm_loc,CFA_EQ_Map(ui))
!                PRINT*,"Source"
!                PRINT*,WORK_VEC

            END DO ! m
        END DO ! l
    END DO ! ui


ELSE IF ( Matrix_Format == 'CCS' ) THEN

    DO ui = 1,Num_Matrices
        DO l = 0,LM_LENGTH

            FP_Laplace_Vector(:,l,ui) = MVMULT_CCS( NUM_R_NODES,                                &
                                                    Laplace_NNZ,                                &
                                                    Laplace_Matrix_VAL(:,l,ui),                 &
                                                    Laplace_Matrix_COL(:,l),                    &
                                                    Laplace_Matrix_ROW(:,l),                    &
                                                    FP_Coeff_Vector(:,l,CFA_EQ_Map(ui))   )

        END DO ! l
    END DO ! ui

END IF




CALL Laplace_Test_Residual()







! Check Residual against Criteria

!
!   Has the Solver achieved convergence
!
SELECT CASE (Convergence_Type )
CASE(1)
    ! L_Inf Error
    Convergence_Stat = MAXVAL(ABS(FP_Residual_Vector))

CASE(2)
    ! L_One Error
    Convergence_Stat = SUM(ABS(FP_Residual_Vector))
CASE(3)
    ! L_Two Error


END SELECT


END SUBROUTINE Calc_FP_Residual








!+301+###########################################################################!
!                                                                                !
!           Call Solve_FP_System                                                 !
!                                                                                !
!################################################################################!
SUBROUTINE Solve_FP_System()

INTEGER                                                                     ::  INFO, LDAB
INTEGER, DIMENSION(NUM_R_NODES)                                             ::  IPIV


REAL(KIND = idp)                                                            ::  SCALE_FACTOR
COMPLEX(KIND = idp), DIMENSION(1:NUM_R_NODES)                               ::  WORK_VEC
COMPLEX(KIND = idp), DIMENSION(0:NUM_R_NODES)                               ::  WORK_VECB

COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:)                            ::  WORK_MAT
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)                              ::  WORK_ELEM_VAL

INTEGER                                                                     ::  NNZ

CHARACTER(LEN = 70)                                                         ::  FILE_NAME
INTEGER                                                                     ::  FILE_ID

INTEGER                                                                     ::  l, m, k, lm_loc, ui
INTEGER                                                                     ::  Guess_Flag
INTEGER                                                                     ::  Mat_Loc

IF (MATRIX_FORMAT =='FULL') THEN

    ALLOCATE (WORK_MAT(1:NUM_R_NODES, 1:NUM_R_NODES))

ELSE IF (MATRIX_FORMAT == 'CCS' ) THEN

    NNZ = NUM_R_ELEMENTS*(DEGREE + 1)*(DEGREE + 1) - NUM_R_ELEMENTS + 1
    ALLOCATE(WORK_ELEM_VAL(0:NNZ-1))

END IF



WORK_VECB = 0




IF (( LINEAR_SOLVER == "CHOL" ) .AND. (MCF_Flag == 0 ) ) THEN
    
    !
    !   This only needs to be done everytime the radial mesh is defined/redefined.
    !   This performs Cholesky factorization on the stiffness matrix and overwrites
    !   the stiffness matrix variables STF_ELEM_VAL, STF_ROW_IND, and STF_COL_PTR to
    !   represent the factorization matrix, L.  This matrix can then be reused to
    !   solve the linear system using forward and backward substitution.
    !

    CALL Cholesky_Factorization()
    MCF_Flag = 1


END IF






113 FORMAT (A,I2.2,A,I5.5,A)







IF (LINEAR_SOLVER =='Full') THEN
    !####################################!
    !
    !           Full Matrix Solver       !
    !
    !####################################!
    PRINT*,"In Full Solver", NUM_CFA_EQs
    DO ui = 1,NUM_CFA_EQs
        DO l = 0,L_LIMIT
            DO m = -l,l

                lm_loc = FP_LM_Map(l,m)
                mat_loc = CFA_Mat_Map(CFA_Eq_Map(ui))
                WORK_MAT = Laplace_Matrix_Full(:,:,l,mat_loc)
                WORK_VEC = FP_Source_Vector(:,lm_loc,ui)

!                PRINT*,"Before DIRICHLET_BC",ui,l,m
!                PRINT*,WORK_VEC

                CALL DIRICHLET_BC(WORK_MAT, WORK_VEC, l, m, ui)



                CALL NEUMANN_BC(l, WORK_VEC)




!                PRINT*,"Before JACOBI_CONDITIONING"
                CALL JACOBI_CONDITIONING(WORK_MAT, WORK_VEC)



                CALL ZGESV(NUM_R_NODES, 1, WORK_MAT, NUM_R_NODES, IPIV, WORK_VEC, NUM_R_NODES, INFO)
                IF (INFO > 0) THEN
                    print*,"DGESV has failed with INFO = ",INFO
                END IF




            !CALL PRECOND_CONJ_GRAD_FULL(WORK_MAT, WORK_VEC)



                lm_loc = FP_LM_Map(l,m)
                FP_Update_Vector(:,lm_loc,ui) = WORK_VEC(:)-FP_Coeff_Vector(:,lm_loc,CFA_EQ_Map(ui))
                FP_Coeff_Vector(:,lm_loc,CFA_EQ_Map(ui)) = WORK_VEC(:)

                
            END DO ! m Loop
        END DO ! l Loop
    END DO ! ui Loop


ELSE IF (LINEAR_SOLVER == 'CCS') THEN
            !#######################################################################!
            !                                                                       !
            !           CCS Preconditioned Conjugate Gradient Matrix Solver         !
            !                                                                       !
            !#######################################################################!
    PRINT*,"In CCS Solver"
    DO ui = 1,NUM_CFA_EQs
        DO l = 0,L_LIMIT

            DO m = -l,l

                Mat_Loc = CFA_MAT_Map(ui)
                WORK_ELEM_VAL(:) = Laplace_Matrix_VAL(:,l,Mat_Loc)
                WORK_VEC(:) = FP_Source_Vector(:,m,l)



                CALL DIRICHLET_BC_CCS(  NUM_R_NODES,                &
                                        Laplace_NNZ,                &
                                        l,                          &
                                        m,                          &
                                        WORK_ELEM_VAL,              &
                                        Laplace_Matrix_COL(:,l),    &
                                        Laplace_Matrix_ROW(:,l),    &
                                        WORK_VEC                    )



                CALL NEUMANN_BC_CCS(    NUM_R_NODES,                &
                                        Laplace_NNZ,                &
                                        l,                          &
                                        m,                          &
                                        WORK_ELEM_VAL,              &
                                        Laplace_Matrix_COL(:,l),    &
                                        Laplace_Matrix_ROW(:,l),    &
                                        WORK_VEC                    )




                GUESS_FLAG = 0  !!! 0 Means no guess, 1 means guess


                IF ( 1 == 1) THEN

                    WORK_VECB = WORK_VEC

                    CALL PRECOND_CONJ_GRAD_CCS(NUM_R_NODES,                 &
                                                Laplace_NNZ,                &
                                                WORK_ELEM_VAL,              &
                                                Laplace_Matrix_COL(:,l),    &
                                                Laplace_Matrix_ROW(:,l),    &
                                                WORK_VECB,                  &
                                                GUESS_FLAG,                 &
                                                WORK_VECB,                  &
                                                0                           )


                    GUESS_FLAG = 1

                END IF




                CALL PRECOND_CONJ_GRAD_CCS( NUM_R_NODES,                &
                                            Laplace_NNZ,                &
                                            WORK_ELEM_VAL,              &
                                            Laplace_Matrix_COL(:,l),     &
                                            Laplace_Matrix_ROW(:,l),     &
                                            WORK_VECB,                  &
                                            GUESS_FLAG,                 &
                                            WORK_VECB,                  &
                                            0                           )

                lm_loc = l*(l+1)+ m
                FP_Update_Vector(:,lm_loc,ui) = WORK_VEC(:)

            END DO ! m
        END DO ! l
    END DO ! ui


ELSE IF (LINEAR_SOLVER == "CHOL") THEN
            !#######################################################################!
            !                                                                       !
            !               CCS Cholesky Factorization Matrix Solver                !
            !                                                                       !
            !#######################################################################!
    PRINT*,"In CHOL Solver"
    DO ui = 1,NUM_CFA_EQs
        DO l = 0,L_LIMIT
            DO m = -l,l

                WORK_VEC = FP_Source_Vector(:,m,l)

                CALL DIRICHLET_BC_CHOL( NUM_R_NODES,                &
                                        Laplace_NNZ,                &
                                        l,                          &
                                        m,                          &
                                        Laplace_Factored_COL(:,l),  &
                                        Laplace_Factored_ROW(:,l),  &
                                        WORK_VEC                    )

                CALL NEUMANN_BC_CCS(    NUM_R_NODES,                &
                                        Laplace_NNZ,                &
                                        l,                          &
                                        m,                          &
                                        WORK_ELEM_VAL,              &
                                        Laplace_Factored_COL(:,l),  &
                                        Laplace_Factored_ROW(:,l),  &
                                        WORK_VEC                    )

                

                CALL CCS_Forward_Substitution(  NUM_R_NODES,                    &
                                                Laplace_NNZ,                    &
                                                Laplace_Factored_VAL(:,l,ui),   &
                                                Laplace_Factored_COL(:,l),      &
                                                Laplace_Factored_ROW(:,l),      &
                                                WORK_VEC                        )


                CALL CCS_Back_Substitution(     NUM_R_NODES,                    &
                                                Laplace_NNZ,                    &
                                                Laplace_Factored_VAL(:,l,ui),   &
                                                Laplace_Factored_COL(:,l),      &
                                                Laplace_Factored_ROW(:,l),      &
                                                WORK_VEC                        )


                lm_loc = l*(l+1)+ m
                FP_Update_Vector(:,lm_loc,ui) = WORK_VEC(:)

            END DO
        END DO
    END DO


END IF








IF (MATRIX_FORMAT =='FULL') THEN

    DEALLOCATE (WORK_MAT)

ELSE IF (MATRIX_FORMAT == 'CCS' ) THEN

    DEALLOCATE(WORK_ELEM_VAL)

END IF



END SUBROUTINE Solve_FP_System








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


120 FORMAT (A61)
121 FORMAT (A1)
122 FORMAT (A41,I2.2)
123 FORMAT (A38,ES22.15)

109 FORMAT (A,I2.2,A,I2.2)
110 FORMAT (11X,A1,16X,A18,9X,A13,10X,A18,10X,A11,14X,A11,14X,A11)
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

    ! Write Results Table Header to Screen
    IF (( WRITE_REPORT_FLAG == 1) .OR. (WRITE_REPORT_FLAG == 3) ) THEN
        IF ( Rank == 0 ) THEN
            PRINT*,"+++++++++++++++++++ myID,",Rank," Iteration",Iter,"++++++++++++++++++++++"
            WRITE(*,110)"r","Analytic Potential","Psi Potential","AlphaPsi Potential","Beta Value1","Beta Value2","Beta Value3"
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

    ! Write Results Table Header to Screen
    IF (( WRITE_REPORT_FLAG == 1) .OR. (WRITE_REPORT_FLAG == 3) ) THEN
        IF ( Rank == 0 ) THEN
            PRINT*,"+++++++++++++++++++ myID,",Rank," Iteration",Iter,"++++++++++++++++++++++"
            WRITE(*,110)"r","Psi Potential","AlphaPsi Potential","Beta Value1","Beta Value2","Beta Value3"
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
                                      PsiPot_Val,                &
                                      AlphaPsiPot_Val,           &
                                      Return_Beta1/Shift_Units,  &
                                      Return_Beta2,              &
                                      Return_Beta3
            END IF
        END DO ! j Loop

    END DO  ! i Loop




END IF



WRITE( FRAME_REPORT_FILE_ID, '(4/)')

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

120 FORMAT (A61)
121 FORMAT (A1)
122 FORMAT (A41,I2.2)
123 FORMAT (A38,ES22.15)

109 FORMAT (A,I2.2,A,I2.2)
110 FORMAT (11X,A1,16X,A18,9X,A13,10X,A18,10X,A11,14X,A11,14X,A11)
111 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15)
112 FORMAT (A43,I2.2,A2,I2.2,A4)




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


END SUBROUTINE OUTPUT_GUESS











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

    
    f = C_One*r*100 + C_Two*(r*100)**(-2)

    CALL Calc_3D_Values_At_Location( r, theta, phi,                              &
                                  Return_Psi, Return_AlphaPsi,                &
                                  Return_Beta1, Return_Beta2, Return_Beta3    )

!    WRITE(*,*)r,Return_Beta1,f, abs(Return_Beta1-f), abs(Return_Beta1-f)/abs(f)

END DO  ! i Loop



END SUBROUTINE Laplace_Test_Residual





END MODULE FP_Method_Module
