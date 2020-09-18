   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_FP_Method_Module                                                    !##!
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

USE Poseidon_Constants_Module, &
            ONLY :  idp, pi, fdp


USE Poseidon_Parameters, &
            ONLY :  CUR_ITERATION,              &
                    MAX_ITERATIONS,             &
                    CONVERGENCE_CRITERIA,       &
                    CONVERGENCE_FLAG

USE Poseidon_Variables_Module, &
            ONLY :  LM_LENGTH,                 &
                    NUM_R_NODES


USE Poseidon_FP_Variables_Module,  &
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
                    CFA_Eq_Map,                 &
                    Laplace_NNZ,                &
                    Num_Matrices

USE Poseidon_Matrix_Functions_Module, &
            ONLY :  MVMULT_FULL,                &
                    MVMULT_CCS

USE Poseidon_FP_Source_Vector_Module, &
            ONLY :  Calc_FP_Source_Vector,          &
                    Allocate_FP_Source_Variables,   &
                    Deallocate_FP_Source_Variables


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

!
!   Begin Method
!




CALL Calc_FP_Source_Vector()

DO WHILE ( CONVERGED .EQV. .FALSE. )

    timer(1) = MPI_Wtime()
    
    ! Call Solve_FP_System()
    
    ! Call Update_FP_Coefficients()

    CALL Calc_FP_Source_Vector()

    Call Check_FP_Convergence(Converged)

    Cur_Iteration = Cur_Iteration + 1
END DO ! Converged Loop



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



IF ( Convergence_Stat > CONVERGENCE_CRITERIA) THEN
    CONVERGED = .TRUE.
    CONVERGENCE_FLAG = 1

!
!   Has the Solver taken the max number of iterations allowed?
!
ELSEIF ( CUR_ITERATION > MAX_ITERATIONS ) THEN
    CONVERGED = .TRUE.
    CONVERGENCE_FLAG = 2


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

INTEGER                                             ::  ui, l


! Multi Matrices and Coeff Vectors to form Laplace
IF ( Matrix_Format == 'Full' ) THEN

    DO ui = 1,Num_Matrices
        DO l = 0,LM_LENGTH
    
            FP_Laplace_Vector(:,l,ui) = MVMULT_FULL( Laplace_Matrix_Full(:,:,l,ui),             &
                                                     FP_Coeff_Vector(:,l,CFA_EQ_Map(ui)), &
                                                     NUM_R_NODES, NUM_R_NODES                   )

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


 
! Add to FP_Source_Vector
FP_Residual_Vector = FP_Laplace_Vector+FP_Source_Vector



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
COMPLEX(KIND = idp), DIMENSION(0:NUM_R_NODES-1)                             ::  WORK_VEC
COMPLEX(KIND = idp), DIMENSION(0:NUM_R_NODES-1)                             ::  WORK_VECB

COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:)                            ::  WORK_MAT
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)                              ::  WORK_ELEM_VAL

INTEGER                                                                     ::  NNZ


INTEGER                                                                     ::  l, m, k
INTEGER                                                                     ::  Guess_Flag

IF (MATRIX_FORMAT =='FULL') THEN

    ALLOCATE (WORK_MAT(0:NUM_R_NODES-1, 0:NUM_R_NODES-1))

ELSE IF (MATRIX_FORMAT == 'CCS' ) THEN

    NNZ = NUM_R_ELEMENTS*(DEGREE + 1)*(DEGREE + 1) - NUM_R_ELEMENTS + 1
    ALLOCATE(WORK_ELEM_VAL(0:NNZ-1))

END IF



WORK_VECB = 0




IF (( LINEAR_SOLVER == "CHOL" ) .AND. (Matrix_Cholesky_Factorized_Flag .EQV. .FALSE.) ) THEN

    !
    !   This only needs to be done everytime the radial mesh is defined/redefined.
    !   This performs Cholesky factorization on the stiffness matrix and overwrites
    !   the stiffness matrix variables STF_ELEM_VAL, STF_ROW_IND, and STF_COL_PTR to
    !   represent the factorization matrix, L.  This matrix can then be reused to
    !   solve the linear system using forward and backward substitution.
    !

    CALL Cholesky_Factorization()
    Matrix_Cholesky_Factorized_Flag = .TRUE.


END IF






113 FORMAT (A,I2.2,A,I5.5,A)




DO l = 0,L_LIMIT

    DO m = -l,l


        IF (LINEAR_SOLVER =='FULL') THEN
            !####################################!
            !
            !           Full Matrix Solver       !
            !
            !####################################!


            WORK_MAT = STF_MAT(:,:,l)
            WORK_VEC = Source_Vector(:,m,l)




            CALL DIRICHLET_BC(WORK_MAT, WORK_VEC, l, m)

            CALL NEUMANN_BC(l, WORK_VEC)

            WRITE(FILE_NAME,113) "OUTPUT/Matrix_D",DEGREE,"_RE",NUM_R_ELEMENTS,".out"
            CALL MATRIX_TO_AIJ(REAL(WORK_MAT(:,:),KIND = idp), NUM_R_NODES, NUM_R_NODES, FILE_NAME)


            CALL JACOBI_CONDITIONING(WORK_MAT, WORK_VEC)



            WRITE(FILE_NAME,113) "OUTPUT/Matrix_D",DEGREE,"_RE",NUM_R_ELEMENTS,"_COND.out"
            CALL MATRIX_TO_AIJ(REAL(WORK_MAT(:,:),KIND = idp), NUM_R_NODES, NUM_R_NODES, FILE_NAME)


            CALL ZGESV(NUM_R_NODES, 1, WORK_MAT, NUM_R_NODES, IPIV, WORK_VEC, NUM_R_NODES, INFO)
            IF (INFO > 0) THEN
                print*,"DGESV has failed with INFO = ",INFO
            END IF




        !CALL PRECOND_CONJ_GRAD_FULL(WORK_MAT, WORK_VEC)









        ELSE IF (LINEAR_SOLVER == 'CCS') THEN
            !#######################################################################!
            !                                                                       !
            !           CCS Preconditioned Conjugate Gradient Matrix Solver         !
            !                                                                       !
            !#######################################################################!

            WORK_ELEM_VAL = STF_ELEM_VAL(:,l)
            WORK_VEC = Source_Vector(:,m,l)



            CALL DIRICHLET_BC_CCS(  NUM_R_NODES,   NNZ, l, m,                               &
                                    WORK_ELEM_VAL, STF_COL_PTR, STF_ROW_IND, WORK_VEC)



            CALL NEUMANN_BC_CCS(    NUM_R_NODES, NNZ, l, m,                                 &
                                    WORK_ELEM_VAL, STF_COL_PTR, STF_ROW_IND, WORK_VEC)





            GUESS_FLAG = 0  !!! 0 Means no guess, 1 means guess


            IF ( 1 == 1) THEN

                WORK_VECB = WORK_VEC

                CALL PRECOND_CONJ_GRAD_CCS(NUM_R_NODES, NNZ, WORK_ELEM_VAL, STF_COL_PTR, STF_ROW_IND, &
                                            WORK_VECB, GUESS_FLAG, WORK_VECB,0)


                GUESS_FLAG = 1

            END IF




            CALL PRECOND_CONJ_GRAD_CCS(NUM_R_NODES, NNZ, WORK_ELEM_VAL, STF_COL_PTR, STF_ROW_IND,       &
                                        WORK_VEC,GUESS_FLAG, WORK_VECB,0)





        ELSE IF (LINEAR_SOLVER == "CHOL") THEN
            !#######################################################################!
            !                                                                       !
            !               CCS Cholesky Factorization Matrix Solver                !
            !                                                                       !
            !#######################################################################!


            WORK_VEC = Source_Vector(:,m,l)

            CALL DIRICHLET_BC_CHOL(  NUM_R_NODES, STF_NNZ, l, m,                                &
                                    STF_COL_PTR, STF_ROW_IND, WORK_VEC)


            CALL NEUMANN_BC_CCS(    NUM_R_NODES, STF_NNZ, l, m,                                 &
                                    WORK_ELEM_VAL, STF_COL_PTR, STF_ROW_IND, WORK_VEC)


            

            CALL CCS_Forward_Substitution(NUM_R_NODES, STF_NNZ, STF_ELEM_VAL(:,l), STF_COL_PTR, STF_ROW_IND, WORK_VEC )



            CALL CCS_Back_Substitution(NUM_R_NODES, STF_NNZ, STF_ELEM_VAL(:,l), STF_COL_PTR, STF_ROW_IND, WORK_VEC )





        END IF



        Do k = 0,NUM_R_NODES - 1
            Coefficient_Vector(k,m,l) = WORK_VEC(k)
        END DO



    END DO
END DO





IF (MATRIX_FORMAT =='FULL') THEN

    DEALLOCATE (WORK_MAT)

ELSE IF (MATRIX_FORMAT == 'CCS' ) THEN

    DEALLOCATE(WORK_ELEM_VAL)

END IF



END SUBROUTINE Solve_FP_System








END MODULE Poseidon_FP_Method_Module
