   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE FP_System_Solvers_Module                                                      !##!
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
                    Beta_Prob_Dim,              &
                    Prob_Dim

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
                    FP_Source_Vector_A,         &
                    FP_Source_Vector_B,         &
                    FP_Coeff_Vector_A,          &
                    FP_Coeff_Vector_B,          &
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
                    CFA_EQ_Flags,               &
                    CFA_Eq_Map,                 &
                    CFA_Var_Map,                &
                    CFA_MAT_Map,                &
                    Laplace_NNZ,                &
                    Factored_NNZ,               &
                    Num_Matrices,               &
                    MCF_Flag,                   &
                    FP_Anderson_M,              &
                    Beta_IPIV,                  &
                    Beta_Diagonals,             &
                    Beta_Bandwidth,             &
                    Beta_MVL_Banded,            &
                    Beta_Factorized_Flag

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
                    DIRICHLET_BC_Beta_Banded,       &
                    Neumann_BC,                     &
                    Neumann_BC_CCS

USE IO_FP_Linear_System, &
            ONLY :  Output_Laplace_Beta,            &
                    Output_Laplace

USE Poseidon_Cholesky_Module,   &
            ONLY :  CCS_Back_Substitution,          &
                    CCS_Forward_Substitution,       &
                    Cholesky_Factorization

USE Linear_Solvers_And_Preconditioners, &
            ONLY :  PRECOND_CONJ_GRAD_CCS,          &
                    JACOBI_CONDITIONING,            &
                    Jacobi_Conditioning_Beta

USE Poseidon_IO_Module, &
            ONLY :  Clock_In,                       &
                    OPEN_ITER_REPORT_FILE,          &
                    CLOSE_ITER_REPORT_FILE

USE FP_Functions_Mapping, &
            ONLY :  FP_Beta_Array_Map,              &
                    FP_Array_Map_TypeB,             &
                    FP_Array_Map

USE Functions_Domain_Maps, &
                ONLY :  Map_To_lm,                  &
                        Map_To_FEM_Node

USE FP_Factorize_Beta_Banded, &
            ONLY :  Factorize_Beta_Banded,          &
                    Jacobi_PC_MVL_Banded_Vector

IMPLICIT NONE





CONTAINS







!+301+###########################################################################!
!                                                                                !
!           Call Solve_FP_System                                                 !
!                                                                                !
!################################################################################!
SUBROUTINE Solve_FP_System()

INTEGER                                                                     ::  INFO
INTEGER, DIMENSION(NUM_R_NODES)                                             ::  IPIV


COMPLEX(KIND = idp), DIMENSION(1:NUM_R_NODES)                               ::  WORK_VEC
COMPLEX(KIND = idp), DIMENSION(0:NUM_R_NODES)                               ::  WORK_VECB

COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:)                            ::  WORK_MAT
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)                              ::  WORK_ELEM_VAL

INTEGER                                                                     ::  l, m, lm_loc, ui
REAL(idp)                                                                   ::  omega


REAL(KIND = idp), DIMENSION(1:4)                        :: timer

WORK_VECB = 0

IF ( Verbose_Flag ) THEN
    PRINT*,"In Anderson FP loop, In Solve_FP_System."
END IF

Omega = 2.0_idp/(1.0_idp + sin( pi/(Num_R_Nodes) ) )



IF (( LINEAR_SOLVER == "CHOL" ) .AND. (MCF_Flag == 0 ) ) THEN
    
    !
    !   This only needs to be done everytime the radial mesh is defined/redefined.
    !   This performs Cholesky factorization on the stiffness matrix and overwrites
    !   the stiffness matrix variables STF_ELEM_VAL, STF_ROW_IND, and STF_COL_PTR to
    !   represent the factorization matrix, L.  This matrix can then be reused to
    !   solve the linear system using forward and backward substitution.
    !

    timer(1) = MPI_Wtime()
    CALL Cholesky_Factorization()
    MCF_Flag = 1


    timer(2) = MPI_Wtime()
    CALL Clock_IN(timer(2)-timer(1), 11)
END IF




timer(1) = MPI_Wtime()
IF (LINEAR_SOLVER =='Full') THEN
    !####################################!
    !
    !           Full Matrix Solver       !
    !
    !####################################!
    ALLOCATE (WORK_MAT(1:NUM_R_NODES, 1:NUM_R_NODES))
!    PRINT*,"In Full Solver"
    DO ui = 1,2
    IF ( CFA_EQ_Flags(ui) == 1 ) THEN
        DO l = 0,L_LIMIT
        DO m = -l,l

            lm_loc   = Map_To_lm(l,m)
            WORK_MAT = Laplace_Matrix_Full(:,:,l)
            WORK_VEC = FP_Source_Vector_A(:,lm_loc,ui)


            CALL DIRICHLET_BC(WORK_MAT, WORK_VEC, l, m, ui)


            CALL NEUMANN_BC(l, WORK_VEC)

            
            CALL JACOBI_CONDITIONING(WORK_MAT, WORK_VEC)


            CALL ZGESV(NUM_R_NODES, 1, WORK_MAT, NUM_R_NODES, IPIV, WORK_VEC, NUM_R_NODES, INFO)
            IF (INFO > 0) THEN
                print*,"DGESV has failed with INFO = ",INFO
            END IF

        
            FP_Update_Vector(:,lm_loc,ui) = WORK_VEC(:)-FP_Coeff_Vector_A(:,lm_loc,CFA_EQ_Map(ui))
            FP_Coeff_Vector_A(:,lm_loc,ui) = WORK_VEC(:)



            
        END DO ! m Loop
        END DO ! l Loop
    END IF
    END DO ! ui Loop


ELSE IF (LINEAR_SOLVER == "CHOL") THEN
            !#######################################################################!
            !                                                                       !
            !               Cholesky Factorization Matrix Solver                !
            !                                                                       !
            !#######################################################################!


    ALLOCATE(WORK_ELEM_VAL(0:Factored_NNZ-1))
    
    DO ui = 1,2
    IF ( CFA_EQ_Flags(ui) == 1 ) THEN
        DO l = 0,L_LIMIT
        DO m = -l,l


            lm_loc = Map_To_lm(l,m)
            WORK_VEC = -FP_Source_Vector_A(:,lm_loc,ui)
            WORK_ELEM_VAL(:) = Laplace_Factored_VAL(:,l,CFA_Var_MAP(ui))

    

!                PRINT*,"Before Dirichelet_BC",ui
            CALL DIRICHLET_BC_CHOL( NUM_R_NODES,                &
                                    Factored_NNZ,               &
                                    l,                          &
                                    m,                          &
                                    Laplace_Factored_COL(:,l),  &
                                    Laplace_Factored_ROW(:,l),  &
                                    WORK_VEC,                   &
                                    ui                          )



!                PRINT*,"Before Neumann_BC",ui
            CALL NEUMANN_BC_CCS(    NUM_R_NODES,                &
                                    Factored_NNZ,               &
                                    l,                          &
                                    m,                          &
                                    WORK_ELEM_VAL,              &
                                    Laplace_Factored_COL(:,l),  &
                                    Laplace_Factored_ROW(:,l),  &
                                    WORK_VEC                    )


            CALL CCS_Forward_Substitution(  NUM_R_NODES,                    &
                                            Factored_NNZ,                   &
                                            WORK_ELEM_VAL,                  &
                                            Laplace_Factored_COL(:,l),      &
                                            Laplace_Factored_ROW(:,l),      &
                                            WORK_VEC                        )


            CALL CCS_Back_Substitution(     NUM_R_NODES,                    &
                                            Factored_NNZ,                   &
                                            WORK_ELEM_VAL,                  &
                                            Laplace_Factored_COL(:,l),      &
                                            Laplace_Factored_ROW(:,l),      &
                                            WORK_VEC                        )



            FP_Update_Vector(:,lm_loc,ui) = WORK_VEC(:)-FP_Coeff_Vector_A(:,lm_loc,ui)
            FP_Coeff_Vector_A( :,lm_loc,ui) = WORK_VEC(:)



        END DO  ! m Loop
        END DO  ! l Loop
    END IF
    END DO  ! ui Loop


END IF









IF (MATRIX_FORMAT =='FULL') THEN

    DEALLOCATE (WORK_MAT)

ELSE IF (MATRIX_FORMAT == 'CCS' ) THEN

    DEALLOCATE(WORK_ELEM_VAL)

END IF



timer(1) = MPI_Wtime()
CALL Clock_In(timer(1)-timer(2),16)

END SUBROUTINE Solve_FP_System
















!+301+###########################################################################!
!                                                                                !
!           Call Solve_FP_System                                                 !
!                                                                                !
!################################################################################!
SUBROUTINE Solve_FP_System_Beta()

INTEGER                                                         ::  INFO
INTEGER, DIMENSION(1:Beta_Prob_Dim)                             ::  IPIV


COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)                  ::  WORK_VEC
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:)                ::  WORK_MAT


INTEGER                                                         ::  ui, re, d, l
INTEGER                                                         ::  Here, There, Thither

REAL(idp), DIMENSION(1:4)                                       ::  timer


IF ( Verbose_Flag ) THEN
    PRINT*,"In Anderson FP loop, In Solve_FP_System_Beta."
END IF



!####################################!
!
!           Full Matrix Solver       !
!
!####################################!
IF (LINEAR_SOLVER =='Full') THEN



    ALLOCATE( WORK_VEC( 1:Beta_Prob_Dim ) )
    ALLOCATE( WORK_MAT( 1:Beta_Prob_Dim, 1:Beta_Prob_Dim ) )

    WORK_MAT(:,:) = Laplace_Matrix_Beta(:,:)
    WORK_VEC(:)   = FP_Source_Vector_B(:,iVB_S)

    

    CALL DIRICHLET_BC_Beta(WORK_MAT, WORK_VEC)


    CALL JACOBI_CONDITIONING_Beta(WORK_MAT, WORK_VEC, Beta_Prob_Dim, Beta_Prob_Dim)
    
    
!    CALL Calc_RCOND_Full( Work_Mat, RCOND )


    CALL ZGESV(Beta_Prob_Dim, 1, WORK_MAT, Beta_Prob_Dim, IPIV, WORK_VEC, Beta_Prob_Dim, INFO)
    IF (INFO .NE. 0) THEN
        print*,"ZGESV has failed with INFO = ",INFO
    END IF

    
    DO ui = 1,3
    DO re = 0,Num_R_Elements-1
    DO d = 0,Degree
    DO l = 1,LM_Length

        Here    = FP_Beta_Array_Map(re,d,ui,l)
        There   = Map_To_FEM_Node(re,d)
        Thither = FP_Array_Map_TypeB(ui,iVB_S,re,d,l)

        FP_Update_Vector(There,l,ui+2) = WORK_VEC(Here)-FP_Coeff_Vector_B(Thither,iVB_S)

    END DO  ! l
    END DO ! d
    END DO ! re
    END DO ! ui

    FP_Coeff_Vector_B(:,iVB_S) = WORK_VEC(:)


    DEALLOCATE( Work_Vec )
    DEALLOCATE( Work_Mat )




!####################################!
!
!         Sparse Matrix Solver       !
!
!####################################!
ELSE IF (LINEAR_SOLVER == "CHOL") THEN


    timer(1) = MPI_Wtime()
    IF ( .NOT. Beta_Factorized_Flag ) THEN
        
        CALL Factorize_Beta_Banded()

    END IF
    timer(2) = MPI_Wtime()
    
    CALL Clock_IN(timer(2)-timer(1), 12)


    ALLOCATE( WORK_VEC( 1:Beta_Prob_Dim ) )
!    ALLOCATE( WORK_MAT( 1:(3*Beta_Diagonals+1), 1:Beta_Prob_Dim ) )
    
!    Work_Mat = Beta_MVL_Banded
    Work_Vec = FP_Source_Vector_B(:,iVB_S)



    CALL DIRICHLET_BC_Beta_Banded(Beta_Prob_Dim, Work_Vec )



    CALL Jacobi_PC_MVL_Banded_Vector( Work_Vec )



!    PRINT*,"Before ZGBTRS "
    CALL ZGBTRS( 'N',                   &
                 Beta_Prob_Dim,         &
                 Beta_Diagonals,        &
                 Beta_Diagonals,        &
                 1,                     &
                 Beta_MVL_Banded,              &
                 3*Beta_Diagonals+1,    &
                 Beta_IPIV,             &
                 -Work_Vec,              &
                 Beta_Prob_Dim,         &
                 INFO                   )

    IF (INFO .NE. 0) THEN
        print*,"ZGBTRS has failed with INFO = ",INFO
    END IF



    DO ui = 1,3
    DO re = 0,Num_R_Elements-1
    DO d = 0,Degree
    DO l = 1,LM_Length

        Here    = FP_Beta_Array_Map(re,d,ui,l)
        There   = Map_To_FEM_Node(re,d)
        Thither = FP_Array_Map_TypeB(ui+2,iVB_S,re,d,l)

        FP_Update_Vector(There,l,ui+2) = WORK_VEC(Here)-FP_Coeff_Vector_B(Thither,iVB_S)

    END DO  ! l
    END DO ! d
    END DO ! re
    END DO ! ui

    FP_Coeff_Vector_B(:,iVB_S) = WORK_VEC(:)
    

    DEALLOCATE( Work_Vec )
!    DEALLOCATE( Work_Mat )


    timer(1) = MPI_Wtime()
    CALL Clock_In(timer(1)-timer(2),17)

END IF




!PRINT*,"End of Routine"

END SUBROUTINE Solve_FP_System_Beta










!+301+###########################################################################!
!                                                                                !
!           Calc_RCOND_Full                                                      !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_RCOND_Full( Work_Mat, RCOND )

REAL(idp), INTENT(OUT)                                                      ::  RCOND
COMPLEX(idp),   DIMENSION(1:Beta_Prob_Dim, 1:Beta_Prob_Dim), INTENT(IN)     ::  Work_Mat

COMPLEX(idp), DIMENSION(2*Beta_Prob_Dim)                                    ::  WORK
REAL(idp), DIMENSION(Beta_Prob_Dim)                                         ::  RWORK
REAL(idp)                                                                   ::  NORM
INTEGER                                                                     ::  INFO

INTEGER                                                                     ::  i


NORM = 0.0_idp
DO i = 1,Beta_Prob_Dim
    NORM = MAX( NORM, ABS(SUM(WORK_MAT(:,i) ) ) )
END DO


CALL ZGECON('1',    &
            Beta_Prob_Dim,      &
            Work_Mat,           &
            Beta_Prob_Dim,      &
            NORM,                &
            RCOND,              &
            WORK,               &
            RWORK,              &
            INFO    )
IF (INFO .NE. 0) THEN
    print*,"ZGECON has failed with INFO = ",INFO
ELSE
    IF ( Verbose_Flag ) THEN
        PRINT*,"Shift Matrix : RCOND = ",RCOND," Norm = ",Norm
    END IF
END IF



END SUBROUTINE Calc_RCOND_Full





END MODULE FP_System_Solvers_Module
