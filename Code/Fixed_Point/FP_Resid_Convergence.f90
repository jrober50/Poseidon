   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE FP_Resid_Convergence_Module                                                  !##!
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

USE Parameters_Variable_Indices, &
            ONLY :  iU_CF,                        &
                    iU_LF,                        &
                    iU_S1,                        &
                    iU_S2,                        &
                    iU_S3,                        &
                    iVB_S


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
            ONLY :  Frame_Update_Table,         &
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
                    CFA_Eq_Flags,               &
                    Verbose_Flag

USE Variables_MPI, &
            ONLY :  myID_Poseidon


USE Variables_FP,  &
            ONLY :  FP_Update_Vector,           &
                    FP_Laplace_Vector,          &
                    FP_Residual_Vector,         &
                    FP_Laplace_Vector_Beta,     &
                    FP_Residual_Vector_Beta,    &
                    FP_Anderson_M

USE Variables_Vectors,  &
            ONLY :  cVA_Load_Vector,         &
                    cVB_Load_Vector,         &
                    cVA_Coeff_Vector,          &
                    cVB_Coeff_Vector

USE Variables_Matrices,  &
            ONLY :  Matrix_Format,              &
                    Linear_Solver,              &
                    Laplace_Matrix_Full,        &
                    Laplace_Matrix_VAL,         &
                    Laplace_Matrix_ROW,         &
                    Laplace_Matrix_COL,         &
                    Laplace_Factored_VAL,       &
                    Laplace_Factored_ROW,       &
                    Laplace_Factored_COL,       &
                    Laplace_Matrix_Beta,        &
                    Beta_MVL_Banded,            &
                    Beta_Diagonals,             &
                    Laplace_NNZ,                &
                    Factored_NNZ,               &
                    Beta_Bandwidth

USE Functions_Matrix, &
            ONLY :  MVMULT_FULL,                &
                    MVMULT_CCS

USE Functions_Mesh, &
            ONLY :  Create_Logarithmic_1D_Mesh,     &
                    Create_Uniform_1D_Mesh

USE FP_Load_Vector_Module, &
            ONLY :  Calc_FP_Load_Vector,          &
                    Allocate_FP_Source_Variables,   &
                    Deallocate_FP_Source_Variables

USE Matrix_Boundary_Condition_Routines,  &
            ONLY :  Dirichlet_BC,                   &
                    Dirichlet_BC_CCS,               &
                    Dirichlet_BC_CHOL,              &
                    Dirichlet_BC_Beta,              &
                    DIRICHLET_BC_Beta_Banded,       &
                    Neumann_BC,                     &
                    Neumann_BC_CCS

USE Matrix_Vector_Laplacian_Routines, &
            ONLY :  Jacobi_PC_MVL_Banded_Vector

USE Matrix_Cholesky_Factorization_Module,   &
            ONLY :  Cholesky_Factorization,         &
                    CCS_Back_Substitution,          &
                    CCS_Forward_Substitution

USE IO_FP_Linear_System, &
            ONLY :  Output_Laplace_Beta,            &
                    Output_Laplace

USE Linear_Solvers_And_Preconditioners, &
            ONLY :  JACOBI_CONDITIONING,            &
                    Jacobi_Conditioning_Beta

USE  Maps_Domain, &
            ONLY :  Map_To_lm




IMPLICIT NONE

CONTAINS


!+201+###########################################################################!
!                                                                                !
!           Check_FP_Convergence                                                 !
!                                                                                !
!################################################################################!
SUBROUTINE Check_FP_Convergence( CONVERGED )

LOGICAL, INTENT(INOUT)                              ::  CONVERGED


REAL(idp), DIMENSION(1:3,1:LM_LENGTH, 1:5)          ::  Resid_Data
INTEGER                                             ::  ui


! Test Residual

! CALL Calc_FP_Residual( Resid_Data )


DO ui = 1,5
    Frame_Update_Table(CUR_ITERATION,ui)      =  MAXVAL( ABS( FP_Update_Vector(:,:,ui) ) )
    Frame_Residual_Table(1,CUR_ITERATION,ui)  =  MAXVAL(Resid_Data(1,:,ui))
    Frame_Residual_Table(2,CUR_ITERATION,ui)  =  MAXVAL(Resid_Data(2,:,ui))
    Frame_Residual_Table(3,CUR_ITERATION,ui)  =  MAXVAL(Resid_Data(3,:,ui))
END DO




IF ( Verbose_Flag .EQV. .TRUE. ) THEN
    WRITE(*,'(/,A,I3.3)') "Convergence Check, Iteration ",CUR_ITERATION

    WRITE(*,'(A,A)')         "Type       :  ",Convergence_Type_Names(Convergence_Type)
    WRITE(*,'(A,ES22.15)')   "Residual   : ", MAXVAL(Resid_Data(Convergence_Type,:,:))
    WRITE(*,'(A,ES22.15)')   "Criteria   : ", Convergence_Criteria
    WRITE(*,'(A,ES22.15,/)') "Max Change : ", MAXVAL(Frame_Update_Table(CUR_ITERATION,:))
END IF

!
!   Has the Solver met the convergence criteria?
!
IF (  ALL( Resid_Data(Convergence_Type,:,:)  < CONVERGENCE_CRITERIA ) ) THEN
    CONVERGED = .TRUE.
    CONVERGENCE_FLAG = 1
ELSE
    CONVERGED = .FALSE.
END IF



END SUBROUTINE Check_FP_Convergence








!+202+##########################################################################!
!                                                                               !
!           Calc_FP_Residual                                                    !
!                                                                               !
!###############################################################################!
SUBROUTINE Calc_FP_Residual( Convergence_Stat )

REAL(KIND = idp), DIMENSION(1:3,1:LM_LENGTH,1:5), INTENT(OUT)       ::  Convergence_Stat


INTEGER                                             ::  ui, l, m, lm_loc

COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)      ::  WORK_VEC
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:)    ::  WORK_MAT

REAL(idp), DIMENSION(1:3,1:LM_LENGTH,1:5)                       ::  LOne_Norm
REAL(idp), DIMENSION(1:3,1:LM_LENGTH,1:5)                       ::  LTwo_Norm
REAL(idp), DIMENSION(1:3,1:LM_LENGTH,1:5)                       ::  LInf_Norm


LOne_Norm = 0.0_idp
LTwo_Norm = 0.0_idp
LInf_Norm = 0.0_idp


! Multi Matrices and Coeff Vectors to form Laplace
IF ( Matrix_Format == 'Full' ) THEN

    ALLOCATE( WORK_MAT(1:NUM_R_NODES, 1:NUM_R_NODES) )
    ALLOCATE( WORK_VEC(1:NUM_R_NODES) )
    
    

    DO ui = 1,2
        IF ( CFA_EQ_Flags(ui) == 1 ) THEN

            DO l = 0,L_Limit
                DO m = -l,l
                    lm_loc  = Map_To_lm(l,m)

                    WORK_MAT = Laplace_Matrix_Full(:,:,l)
                    WORK_VEC = cVA_Load_Vector(:,lm_loc,ui)
        
                    CALL DIRICHLET_BC(WORK_MAT, WORK_VEC, l, m, ui)

                    FP_Laplace_Vector(:,lm_loc,ui) = MVMULT_FULL( WORK_MAT,                             &
                                                            cVA_Coeff_Vector(:,lm_loc,ui),   &
                                                            NUM_R_NODES, NUM_R_NODES                    )


                    
                    FP_Residual_Vector(:,lm_loc,ui) = (FP_Laplace_Vector(:,lm_loc,ui) - Work_Vec)


                    LOne_Norm(1,lm_loc,ui) = LOne_Norm(1,lm_loc,ui) + SUM(ABS(FP_Laplace_Vector(:,lm_loc,ui) ) )
                    LOne_Norm(2,lm_loc,ui) = LOne_Norm(2,lm_loc,ui) + SUM(ABS(Work_Vec) )
                    LOne_Norm(3,lm_loc,ui) = LOne_Norm(3,lm_loc,ui) + SUM(ABS(FP_Residual_Vector(:,lm_loc,ui) ) )

                    LTwo_Norm(1,lm_loc,ui) = LTwo_Norm(1,lm_loc,ui)                                             &
                                 + SQRT(REAL(DOT_PRODUCT(FP_Laplace_Vector(:,lm_loc,ui),         &
                                                    FP_Laplace_Vector(:,lm_loc,ui) ), idp) )
                    LTwo_Norm(2,lm_loc,ui) = LTwo_Norm(2,lm_loc,ui) + SQRT(REAL(DOT_PRODUCT(Work_Vec,                &
                                                                   Work_Vec ), idp) )
                    LTwo_Norm(3,lm_loc,ui) = LTwo_Norm(3,lm_loc,ui)                                             &
                                 + SQRT(REAL(DOT_PRODUCT(FP_Residual_Vector(:,lm_loc,ui),        &
                                                    FP_Residual_Vector(:,lm_loc,ui) ), idp) )

                    LInf_Norm(1,lm_loc,ui) = MAXVAL( (/ LInf_Norm(1,lm_loc,ui), ABS(FP_Laplace_Vector(:,lm_loc,ui) ) /) )
                    LInf_Norm(2,lm_loc,ui) = MAXVAL( (/ LInf_Norm(2,lm_loc,ui), ABS(Work_Vec ) /) )
                    LInf_Norm(3,lm_loc,ui) = MAXVAL( (/ LInf_Norm(3,lm_loc,ui), ABS(FP_Residual_Vector(:,lm_loc,ui) ) /) )

    


                END DO ! m
            END DO ! l


        END IF
    END DO ! ui
    
    
    IF ( ANY(CFA_EQ_Flags(3:5) .EQ. 1) ) THEN
        DEALLOCATE( WORK_VEC )
        DEALLOCATE( WORK_MAT )
        ALLOCATE( WORK_VEC(1:Beta_Prob_Dim) )
        ALLOCATE( WORK_MAT(1:Beta_Prob_Dim,1:Beta_Prob_Dim) )

        WORK_MAT(:,:) = Laplace_Matrix_Beta(:,:)
        WORK_VEC(:) = cVB_Load_Vector(:,iVB_S)

        CALL DIRICHLET_BC_Beta(WORK_MAT, WORK_VEC)


        CALL Jacobi_Conditioning_Beta(Work_mat,work_Vec,Beta_Prob_Dim,Beta_Prob_Dim)


        CALL ZGEMV('N',                         &
                    Beta_Prob_Dim,              &
                    Beta_Prob_Dim,              &
                    1.0_idp,                    &
                    Work_Mat,                   &
                    Beta_Prob_Dim,              &
                    cVB_Coeff_Vector(:,iVB_S), &
                    1,                          &
                    0.0_idp,                    &
                    FP_Laplace_Vector_Beta,     &
                    1                       )


        FP_Residual_Vector_Beta = (FP_Laplace_Vector_Beta - Work_Vec)
    


        LOne_Norm(1,:,3:5) = SUM( ABS(FP_Laplace_Vector_Beta) )
        LTwo_Norm(1,:,3:5) = SQRT(REAL(DOT_PRODUCT(FP_Laplace_Vector_Beta,FP_Laplace_Vector_Beta) ) )
        LInf_Norm(1,:,3:5) = MAXVAL( ABS( FP_Laplace_Vector_Beta ) )

        LOne_Norm(2,:,3:5) = SUM(ABS(WORK_VEC))
        LTwo_Norm(2,:,3:5) = SQRT(REAL(DOT_PRODUCT(WORK_VEC,WORK_VEC) ) )
        LInf_Norm(2,:,3:5) = MAXVAL( ABS( WORK_VEC ) )

        LOne_Norm(3,:,3:5) = SUM(ABS(FP_Residual_Vector_Beta))
        LTwo_Norm(3,:,3:5) = SQRT(REAL(DOT_PRODUCT(FP_Residual_Vector_Beta,FP_Residual_Vector_Beta) ) )
        LInf_Norm(3,:,3:5) = MAXVAL( ABS( FP_Residual_Vector_Beta ) )


!        DO i = 1,Num_R_Nodes
!            PRINT*,FP_RESIDUAL_VECTOR_BETA(i)/WORK_VEC(i)
!        END DO

    END IF





!    PRINT*,"In Calc_FP_Residual, A"

ELSE IF ( Matrix_Format == 'CCS' ) THEN

    ALLOCATE( Work_Vec(Num_R_Nodes) )


    DO ui = 1,2

        IF ( CFA_EQ_Flags(ui) == 1 ) THEN
            DO l = 0,L_Limit
                DO m = -l,l


                    lm_loc = Map_To_lm(l,m)


                    WORK_VEC = -cVA_Load_Vector(:,lm_loc,ui)

                    CALL DIRICHLET_BC_CHOL( NUM_R_NODES,                &
                                            Factored_NNZ,               &
                                            l,                          &
                                            m,                          &
                                            Laplace_Factored_COL(:,l),  &
                                            Laplace_Factored_ROW(:,l),  &
                                            WORK_VEC,                   &
                                            ui                          )




                    CALL Calc_CCS_Laplacian_Vector( NUM_R_NODES,                            &
                                                    Factored_NNZ,                            &
                                                    Laplace_Factored_VAL(:,l),              &
                                                    Laplace_Factored_COL(:,l),              &
                                                    Laplace_Factored_ROW(:,l),              &
                                                    cVA_Coeff_Vector(:,lm_loc,ui),    &
                                                    FP_Laplace_Vector(:,lm_loc,ui) )


                    FP_Residual_Vector(:,lm_loc,ui) = (FP_Laplace_Vector(:,lm_loc,ui) - Work_Vec)


                    LOne_Norm(1,lm_loc,ui) = LOne_Norm(1,lm_loc,ui) + SUM(ABS(FP_Laplace_Vector(:,lm_loc,ui) ) )
                    LOne_Norm(2,lm_loc,ui) = LOne_Norm(2,lm_loc,ui) + SUM(ABS(Work_Vec) )
                    LOne_Norm(3,lm_loc,ui) = LOne_Norm(3,lm_loc,ui) + SUM(ABS(FP_Residual_Vector(:,lm_loc,ui) ) )


                    LTwo_Norm(1,lm_loc,ui) = LTwo_Norm(1,lm_loc,ui)                                             &
                                 + SQRT(REAL(DOT_PRODUCT(FP_Laplace_Vector(:,lm_loc,ui),         &
                                                    FP_Laplace_Vector(:,lm_loc,ui) ), idp) )
                    LTwo_Norm(2,lm_loc,ui) = LTwo_Norm(2,lm_loc,ui) + SQRT(REAL(DOT_PRODUCT(Work_Vec,                &
                                                                   Work_Vec ), idp) )
                    LTwo_Norm(3,lm_loc,ui) = LTwo_Norm(3,lm_loc,ui)                                             &
                                 + SQRT(REAL(DOT_PRODUCT(FP_Residual_Vector(:,lm_loc,ui),        &
                                                    FP_Residual_Vector(:,lm_loc,ui) ), idp) )

                    LInf_Norm(1,lm_loc,ui) = MAXVAL( (/ LInf_Norm(1,lm_loc,ui), ABS(FP_Laplace_Vector(:,lm_loc,ui) ) /) )
                    LInf_Norm(2,lm_loc,ui) = MAXVAL( (/ LInf_Norm(2,lm_loc,ui), ABS(Work_Vec ) /) )
                    LInf_Norm(3,lm_loc,ui) = MAXVAL( (/ LInf_Norm(3,lm_loc,ui), ABS(FP_Residual_Vector(:,lm_loc,ui) ) /) )


                END DO ! m
            END DO ! l
        END IF
    END DO ! ui







    DEALLOCATE(Work_Vec)
    ALLOCATE( Work_Vec(1:Beta_Prob_Dim ) )
    Work_Vec = cVB_Load_Vector(:,iVB_S)

    CALL DIRICHLET_BC_Beta_Banded(Beta_Prob_Dim, Work_Vec )

    CALL Jacobi_PC_MVL_Banded_Vector( Work_Vec )

    ! Shift Residual
    CALL ZGBMV('N',                         &
                Beta_Prob_Dim,              &
                Beta_Prob_Dim,              &
                Beta_Diagonals,             &
                2*Beta_Diagonals,           &
                1.0_idp,                    &
                Beta_MVL_Banded,            &
                3*Beta_Diagonals+1,         &
                cVB_Coeff_Vector(:,iVB_S), &
                1,                          &
                0.0_idp,                    &
                FP_Laplace_Vector_Beta,     &
                1                           )


   

!    DO l = 1,Beta_Prob_Dim
!        PRINT*,l,FP_Laplace_Vector_Beta(l),Work_Vec(l),FP_Laplace_Vector_Beta(l)-Work_Vec(l)
!    END DO
!     PRINT*,MAXLOC(abs(FP_Laplace_Vector_Beta(:))),MAXVAL( ABS( FP_Laplace_Vector_Beta ) )

    FP_Residual_Vector_Beta = FP_Laplace_Vector_Beta - Work_Vec
    

    LOne_Norm(1,:,3:5) = SUM( ABS(FP_Laplace_Vector_Beta) )
    LTwo_Norm(1,:,3:5) = SQRT(REAL(DOT_PRODUCT(FP_Laplace_Vector_Beta,FP_Laplace_Vector_Beta) ) )
    LInf_Norm(1,:,3:5) = MAXVAL( ABS( FP_Laplace_Vector_Beta ) )

    LOne_Norm(2,:,3:5) = SUM(ABS(WORK_VEC))
    LTwo_Norm(2,:,3:5) = SQRT(REAL(DOT_PRODUCT(WORK_VEC,WORK_VEC) ) )
    LInf_Norm(2,:,3:5) = MAXVAL( ABS( WORK_VEC ) )

    LOne_Norm(3,:,3:5) = SUM(ABS(FP_Residual_Vector_Beta))
    LTwo_Norm(3,:,3:5) = SQRT(REAL(DOT_PRODUCT(FP_Residual_Vector_Beta,FP_Residual_Vector_Beta) ) )
    LInf_Norm(3,:,3:5) = MAXVAL( ABS( FP_Residual_Vector_Beta ) )


END IF







! Convergence Check
!   Is the relative residual below a tolerance
!   norm/max(norm(laplace),norm(source)) <? Tolerance
!
! Here the relative residual is calculated.


Convergence_Stat = 0.0_idp
DO ui = 1,5

!    IF ( ui == 3) THEN
!        PRINT*, LOne_Norm(1,:,ui)
!        PRINT*, " "
!        PRINT*, LOne_norm(2,:,ui)
!        PRINT*," "
!        PRINT*, LOne_Norm(3,:,ui)
!        PRINT*," ------ "
!    END IF

    IF ( MAXVAL(LOne_Norm(1:2,:,ui)) .NE. 0.0_idp ) THEN
        Convergence_Stat(1,:, ui) = LOne_Norm(3,:,ui)/MAXVAL(LOne_Norm(1:2,:,ui))
    ELSE
        Convergence_Stat(1,:, ui) = LOne_Norm(3,:,ui)
    END IF


    IF ( MAXVAL(LTwo_Norm(1:2,:,ui)) .NE. 0.0_idp ) THEN
        Convergence_Stat(2,:, ui) = LTwo_Norm(3,:,ui)/MAXVAL(LTwo_Norm(1:2,:,ui))
    ELSE
        Convergence_Stat(2,:, ui) = LTwo_Norm(3,:,ui)
    END IF


    IF ( MAXVAL(LInf_Norm(1:2,:,ui)) .NE. 0.0_idp ) THEN
        Convergence_Stat(3,:, ui) = LInf_Norm(3,:,ui)/MAXVAL(LInf_Norm(1:2,:,ui))
    ELSE
        Convergence_Stat(3,:, ui) = LInf_Norm(3,:,ui)
    END IF

END DO


END SUBROUTINE Calc_FP_Residual





!+202+##########################################################################!
!                                                                               !
!           Calc_FP_Residual                                                    !
!                                                                               !
!###############################################################################!
SUBROUTINE Calc_CCS_Laplacian_Vector(N, NNZ, ELEM_VAL, COL_PTR, ROW_IND, VECT, Lap)


INTEGER, INTENT(IN)                                         :: N, NNZ
INTEGER, DIMENSION(0:N), INTENT(IN)                         :: COL_PTR
INTEGER, DIMENSION(0:NNZ - 1),INTENT(IN)                    :: ROW_IND

COMPLEX(KIND = idp), DIMENSION(0:NNZ - 1), INTENT(IN)       :: ELEM_VAL
COMPLEX(KIND = idp), DIMENSION(0:N-1), INTENT(IN)           :: VECT
COMPLEX(KIND = idp), DIMENSION(0:N-1), INTENT(OUT)          :: Lap

COMPLEX(KIND = idp), DIMENSION(0:N-1)                       :: Laplacian_A
COMPLEX(KIND = idp), DIMENSION(0:N-1)                       :: Laplacian_B

INTEGER                                                     :: i,j


! The Matrix that enters is the L matrix of Factorized matrix A = LL^t.

! To calculate the Laplacian, we need to calculate L*L^t*x = L * (L^t*x)


Laplacian_A = 0.0_idp
Laplacian_B = 0.0_idp


! L^t*x
DO i = 0,N-1
    DO j = COL_PTR(i),COL_PTR(i+1)-1

        Laplacian_A(i) = Laplacian_A(i) + ELEM_VAL(j)*Vect(ROW_IND(j))

    END DO ! j Loop
END DO ! i Loop



!   L*(L^t*x)
DO i = 0,N-1
    DO j = COL_PTR(i),COL_PTR(i+1)-1

        Laplacian_B(ROW_IND(j)) = Laplacian_B(ROW_IND(j)) + ELEM_VAL(j)*Laplacian_A(i)
    
    END DO ! j Loop
END DO ! i Loop



Lap = Laplacian_B

END SUBROUTINE Calc_CCS_Laplacian_Vector









END MODULE FP_Resid_Convergence_Module 