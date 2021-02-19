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


!+201+###########################################################################!
!                                                                                !
!           Check_FP_Convergence                                                 !
!                                                                                !
!################################################################################!
SUBROUTINE Check_FP_Convergence( CONVERGED )

LOGICAL, INTENT(INOUT)                              ::  CONVERGED

REAL(idp), DIMENSION(1:3,0:LM_LENGTH-1, 1:5)        ::  Resid_Data
REAL(idp)                                           ::  Convergence_Stat

INTEGER                                             ::  ui, l


! Test Residual
!PRINT*,"Before Calc_FP_Residual"
CALL Calc_FP_Residual( Resid_Data )
!PRINT*,"After Calc_FP_Residual"



DO ui = 1,5
    Frame_Update_Table(CUR_ITERATION,ui)      =  MAXVAL( ABS( FP_Update_Vector(:,:,ui) ) )
    Frame_Residual_Table(1,CUR_ITERATION,ui)  =  MAXVAL(Resid_Data(1,:,ui))
    Frame_Residual_Table(2,CUR_ITERATION,ui)  =  MAXVAL(Resid_Data(2,:,ui))
    Frame_Residual_Table(3,CUR_ITERATION,ui)  =  MAXVAL(Resid_Data(3,:,ui))
END DO




IF ( Verbose_Flag .EQV. .TRUE. ) THEN
    WRITE(*,'(/,A,I3.3)') "Convergence Check, Iteration ",CUR_ITERATION

    WRITE(*,'(A,A)')       "Type     : ",Convergence_Type_Names(Convergence_Type)
    WRITE(*,'(A,ES22.15)') "Residual : ",Convergence_Stat, MAXLOC(Resid_Data(Convergence_Type,0,:) )
    WRITE(*,'(A,ES22.15)') "Criteria : ",Convergence_Criteria

    WRITE(*,'(A,ES22.15,/)') "Max Change : ",MAXVAL(ABS(FP_Update_Vector))
END IF

!
!   Has the Solver met the convergence criteria?
!
IF (  MAXVAL(Resid_Data(Convergence_Type,:,:) ) < CONVERGENCE_CRITERIA ) THEN
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

REAL(KIND = idp), DIMENSION(1:3,0:LM_LENGTH-1,1:5), INTENT(OUT)       ::  Convergence_Stat


INTEGER                                             ::  Convergence_Type = 3

INTEGER                                             ::  ui, l, m, lm_loc, map_loc, i

COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)      ::  WORK_VEC
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:)    ::  WORK_MAT

REAL(idp), DIMENSION(1:3,0:LM_LENGTH-1,1:5)                       ::  LOne_Norm
REAL(idp), DIMENSION(1:3,0:LM_LENGTH-1,1:5)                       ::  LTwo_Norm
REAL(idp), DIMENSION(1:3,0:LM_LENGTH-1,1:5)                       ::  LInf_Norm

! Multi Matrices and Coeff Vectors to form Laplace
IF ( Matrix_Format == 'Full' ) THEN

    ALLOCATE( WORK_MAT(1:NUM_R_NODES, 1:NUM_R_NODES) )
    ALLOCATE( WORK_VEC(1:NUM_R_NODES) )
    
    

    DO ui = 1,2
        IF ( CFA_EQ_Flags(ui) == 1 ) THEN

            LOne_Norm = 0.0_idp
            LTwo_Norm = 0.0_idp
            LInf_Norm = 0.0_idp

            DO l = 0,L_Limit
                DO m = -l,l
                    lm_loc = FP_LM_Map(l,m)
                    map_loc = CFA_MAT_Map(CFA_EQ_Map(ui))

                    WORK_MAT = Laplace_Matrix_Full(:,:,l)
                    WORK_VEC = FP_Source_Vector(:,lm_loc,ui)
        
                    CALL DIRICHLET_BC(WORK_MAT, WORK_VEC, l, m, ui)

                    FP_Laplace_Vector(:,lm_loc,ui) = MVMULT_FULL( WORK_MAT,                             &
                                                            FP_Coeff_Vector(:,lm_loc,CFA_EQ_Map(ui)),   &
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
        WORK_VEC(:) = FP_Source_Vector_Beta(:)

        CALL DIRICHLET_BC_Beta(WORK_MAT, WORK_VEC)

        FP_Laplace_Vector_Beta = MVMULT_FULL( WORK_MAT,                 &
                                                FP_Coeff_Vector_Beta(:),        &
                                                Beta_Prob_Dim, Beta_Prob_Dim    )


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


! Convergence Check
!   Is the relative residual below a tolerance
!   norm/max(norm(laplace),norm(source)) <? Tolerance
!
! Here the relative residual is calculated.


Convergence_Stat = 0.0_idp
DO ui = 1,5

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










END MODULE FP_Resid_Convergence_Module 
