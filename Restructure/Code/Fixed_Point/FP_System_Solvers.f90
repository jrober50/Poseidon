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
                    Factored_NNZ,               &
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
                    SSOR_CONDITIONING,              &
                    Jacobi_Conditioning_Beta

USE Poseidon_IO_Module, &
            ONLY :  Clock_In,                           &
                    OPEN_ITER_REPORT_FILE,              &
                    CLOSE_ITER_REPORT_FILE,             &
                    OUTPUT_FINAL_RESULTS

USE FP_Functions_Mapping, &
            ONLY :  FP_LM_Map


IMPLICIT NONE





CONTAINS







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
REAL(idp)                                                                   ::  omega


REAL(KIND = idp), DIMENSION(1:4)                        :: timer

WORK_VECB = 0



timer(1) = MPI_Wtime()
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



timer(2) = MPI_Wtime()
PRINT*,"Cholesky Factorization Time",timer(2)-timer(1)

Omega = 2.0_idp/(1.0_idp + sin( pi/(Num_R_Nodes) ) )


113 FORMAT (A,I2.2,A,I5.5,A)







IF (LINEAR_SOLVER =='Full') THEN
    !####################################!
    !
    !           Full Matrix Solver       !
    !
    !####################################!
    ALLOCATE (WORK_MAT(1:NUM_R_NODES, 1:NUM_R_NODES))
!    PRINT*,"In Full Solver"
    DO ui = 1,2
!        PRINT*,"CFA_EQ_Flags(ui)",CFA_EQ_Flags(ui)
        IF ( CFA_EQ_Flags(ui) == 1 ) THEN
            DO l = 0,L_LIMIT
                DO m = -l,l

                    lm_loc = FP_LM_Map(l,m)
                    WORK_MAT = Laplace_Matrix_Full(:,:,l)
                    WORK_VEC = FP_Source_Vector(:,lm_loc,ui)
!                    PRINT*,"Work_Vec"
!                    PRINT*,WORK_VEC

                    CALL DIRICHLET_BC(WORK_MAT, WORK_VEC, l, m, ui)


                    CALL NEUMANN_BC(l, WORK_VEC)

                    
!                    CALL SSOR_CONDITIONING(WORK_MAT,WORK_VEC,Omega)
                    CALL JACOBI_CONDITIONING(WORK_MAT, WORK_VEC)

                    CALL ZGESV(NUM_R_NODES, 1, WORK_MAT, NUM_R_NODES, IPIV, WORK_VEC, NUM_R_NODES, INFO)
                    IF (INFO > 0) THEN
                        print*,"DGESV has failed with INFO = ",INFO
                    END IF

                    
                    lm_loc = FP_LM_Map(l,m)


!                    PRINT*,lm_Loc
!                    PRINT*,WORK_VEC
    
                    FP_Update_Vector(:,lm_loc,ui) = WORK_VEC(:)-FP_Coeff_Vector(:,lm_loc,CFA_EQ_Map(ui))
                    FP_Coeff_Vector(:,lm_loc,ui) = WORK_VEC(:)

    

                    
                END DO ! m Loop
            END DO ! l Loop
        END IF
    END DO ! ui Loop


ELSE IF (LINEAR_SOLVER == 'CCS') THEN
            !#######################################################################!
            !                                                                       !
            !           CCS Preconditioned Conjugate Gradient Matrix Solver         !
            !                                                                       !
            !#######################################################################!
    NNZ = NUM_R_ELEMENTS*(DEGREE + 1)*(DEGREE + 1) - NUM_R_ELEMENTS + 1
    ALLOCATE(WORK_ELEM_VAL(0:NNZ-1))

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

    PRINT*,"Factored_NNZ",Factored_NNZ
    ALLOCATE(WORK_ELEM_VAL(0:Factored_NNZ-1))
    
!    PRINT*,"In CHOL Solver"
    DO ui = 1,2

        IF ( CFA_EQ_Flags(ui) == 1 ) THEN
        DO l = 0,L_LIMIT
            DO m = -l,l


                lm_loc = FP_LM_Map(l,m)
                WORK_VEC = -FP_Source_Vector(:,lm_loc,ui)
                WORK_ELEM_VAL(:) = Laplace_Factored_VAL(:,l,ui)

               
!                PRINT*,"Before Dirichelet_BC"
                CALL DIRICHLET_BC_CHOL( NUM_R_NODES,                &
                                        Factored_NNZ,               &
                                        l,                          &
                                        m,                          &
                                        Laplace_Factored_COL(:,l),  &
                                        Laplace_Factored_ROW(:,l),  &
                                        WORK_VEC                    )

!                PRINT*,"Before Neumann_BC"
                CALL NEUMANN_BC_CCS(    NUM_R_NODES,                &
                                        Factored_NNZ,               &
                                        l,                          &
                                        m,                          &
                                        WORK_ELEM_VAL,              &
                                        Laplace_Factored_COL(:,l),  &
                                        Laplace_Factored_ROW(:,l),  &
                                        WORK_VEC                    )

                
               
!                PRINT*,"Before Forward Sub"
                CALL CCS_Forward_Substitution(  NUM_R_NODES,                    &
                                                Factored_NNZ,                   &
                                                WORK_ELEM_VAL,                  &
                                                Laplace_Factored_COL(:,l),      &
                                                Laplace_Factored_ROW(:,l),      &
                                                WORK_VEC                        )

!                PRINT*,"Before Backward Sub"
                CALL CCS_Back_Substitution(     NUM_R_NODES,                    &
                                                Factored_NNZ,                   &
                                                WORK_ELEM_VAL,                  &
                                                Laplace_Factored_COL(:,l),      &
                                                Laplace_Factored_ROW(:,l),      &
                                                WORK_VEC                        )

!                PRINT*,lm_Loc
!                PRINT*,WORK_VEC
                

                lm_loc = l*(l+1)+ m

                FP_Update_Vector(:,lm_loc,ui) = WORK_VEC(:)-FP_Coeff_Vector(:,lm_loc,CFA_EQ_Map(ui))
                FP_Coeff_Vector(:,lm_loc,ui) = WORK_VEC(:)


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



END SUBROUTINE Solve_FP_System
















!+301+###########################################################################!
!                                                                                !
!           Call Solve_FP_System                                                 !
!                                                                                !
!################################################################################!
SUBROUTINE Solve_FP_System_Beta()

INTEGER                                                                     ::  INFO
INTEGER, DIMENSION(1:Beta_Prob_Dim)                                         ::  IPIV


COMPLEX(KIND = idp), DIMENSION(1:Beta_Prob_Dim)                             ::  WORK_VEC
COMPLEX(KIND = idp), DIMENSION(1:Beta_Prob_Dim,1:Beta_Prob_Dim)             ::  WORK_MAT


INTEGER                                                                     ::  ui, re, d, l

INTEGER                                                                     ::  Here, There







IF (LINEAR_SOLVER =='Full') THEN
    !####################################!
    !
    !           Full Matrix Solver       !
    !
    !####################################!
!    PRINT*,"In Beta Solver"


    WORK_MAT(:,:) = Laplace_Matrix_Beta(:,:)
    WORK_VEC(:) = FP_Source_Vector_Beta(:)

!    PRINT*,"Before DIRICHLET_BC"



    CALL DIRICHLET_BC_Beta(WORK_MAT, WORK_VEC)


!    PRINT*,"After Dirichlet_BC"
!    Call Output_Laplace_Beta(Work_Mat, Beta_Prob_Dim, Beta_Prob_Dim, "C")
!    CALL NEUMANN_BC(l, WORK_VEC)

    


!    PRINT*,"Before JACOBI_CONDITIONING"
    CALL JACOBI_CONDITIONING_Beta(WORK_MAT, WORK_VEC, Beta_Prob_Dim, Beta_Prob_Dim)
    

!    Call Output_Laplace_Beta(Work_Mat, Beta_Prob_Dim, Beta_Prob_Dim)


!    PRINT*,Work_Vec

!    PRINT*,"Before ZGESV"
    CALL ZGESV(Beta_Prob_Dim, 1, WORK_MAT, Beta_Prob_Dim, IPIV, WORK_VEC, Beta_Prob_Dim, INFO)
    IF (INFO > 0) THEN
        print*,"DGESV has failed with INFO = ",INFO
    END IF



    FP_Coeff_Vector_Beta(:) = WORK_VEC(:)

!    PRINT*,"FP_Coeff_Vector_Beta"
!    PRINT*,FP_Coeff_Vector_Beta

    DO ui = 1,3
        DO re = 0,Num_R_Elements -1
            DO d = 0,Degree
                DO l = 1,LM_Length

                    Here = (re*Degree + d) * 3 * LM_Length  &
                         + (ui - 1) * LM_Length             &
                         + l

                    There = (re*Degree + d) + 1

!                    PRINT*,"Here",FP_Coeff_Vector_Beta(Here),Here, There, ui, re, d, l
!                    FP_Coeff_Vector(There,l-1,ui+2) = FP_Coeff_Vector_Beta(Here)
                    FP_Update_Vector(There,l,ui+2) = WORK_VEC(Here)-FP_Coeff_Vector(There,l,CFA_EQ_Map(ui))
                    FP_Coeff_Vector(There,l-1,ui+2) = Work_Vec(Here)


                END DO  ! l
            END DO ! d
        END DO ! re
        
    END DO ! ui



END IF





END SUBROUTINE Solve_FP_System_Beta








END MODULE FP_System_Solvers_Module
