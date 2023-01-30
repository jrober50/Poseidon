   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Allocation_XCFC_Linear_Systems                                               !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   Allocate_XCFC_Linear_Systems                                        !##!
!##!    +102+   Deallocate_XCFC_Linear_Systems                                      !##!
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
            ONLY : idp

USE Poseidon_Parameters, &
            ONLY :  DEGREE,                 &
                    L_LIMIT,                &
                    Max_Iterations

Use Variables_Derived, &
            ONLY :  Num_R_Nodes,            &
                    Num_R_Nodesp1,          &
                    LM_Length,              &
                    Var_Dim,                &
                    Prob_Dim,               &
                    iVB_Prob_Dim

USE Variables_FP,  &
            ONLY :  FP_Update_Vector,       &
                    FP_Residual_Vector,     &
                    FP_Laplace_Vector,      &
                    FP_Diagnostics_Flag,    &
                    FP_Iter_Matrix_Storage, &
                    FP_Iter_Load_Storage,   &
                    FP_Iteration_Log,       &
                    Resid_Norms,            &
                    Update_Norms
            

USE Variables_Vectors,  &
            ONLY :  cVA_Load_Vector,        &
                    cVB_Load_Vector,        &
                    cVA_Coeff_Vector,       &
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
                    zMB_Matrix_Banded,            &
                    zMB_Matrix_Diagonal,          &
                    iMB_Diagonals,              &
                    iMB_IPIV,                   &
                    Laplace_NNZ,                &
                    Factored_NNZ,               &
                    zMA_First_Col_Storage,       &
                    zMA_Last_Col_Storage,        &
                    zMB_First_Col_Storage,  &
                    zMB_Last_Col_Storage

USE Flags_Initialization_Module, &
            ONLY :  lPF_Init_Flags,         &
                    iPF_Init_Alloc_LinSys


IMPLICIT NONE


CONTAINS
 !+101+####################################################!
!                                                           !
!          Allocate_Poseidon_Variables                      !
!                                                           !
 !#########################################################!
SUBROUTINE Allocate_XCFC_Linear_Systems()


IF ( MATRIX_FORMAT == 'Full' ) THEN

    ALLOCATE( Laplace_Matrix_Full(1:NUM_R_NODES,1:NUM_R_NODES,0:L_LIMIT) )
    ALLOCATE( Laplace_Matrix_Beta(1:iVB_Prob_Dim,1:iVB_Prob_Dim) )


ELSEIF ( MATRIX_FORMAT == 'CCS' ) THEN

    ALLOCATE( Laplace_Matrix_VAL(0:Laplace_NNZ-1, 0:L_LIMIT) )
    ALLOCATE( Laplace_Matrix_ROW(0:Laplace_NNZ-1, 0:L_LIMIT) )
    ALLOCATE( Laplace_Matrix_COL(0:NUM_R_NODES, 0:L_LIMIT) )

    ALLOCATE( Laplace_Factored_VAL(0:Laplace_NNZ-1, 0:L_LIMIT) )
    ALLOCATE( Laplace_Factored_ROW(0:Laplace_NNZ-1, 0:L_LIMIT) )
    ALLOCATE( Laplace_Factored_COL(0:NUM_R_NODES, 0:L_LIMIT) )

    ALLOCATE( iMB_IPIV(1:iVB_Prob_Dim) )
    ALLOCATE( zMB_Matrix_Banded(1:(3*iMB_Diagonals+1), 1:iVB_Prob_Dim))
    ALLOCATE( zMB_Matrix_Diagonal(1:iVB_Prob_Dim) )

    ALLOCATE( zMA_First_Col_Storage(0:DEGREE,0:L_LIMIT)   )
    ALLOCATE( zMA_Last_Col_Storage(0:DEGREE,0:L_LIMIT)    )

    ALLOCATE( zMB_First_Col_Storage(1:LM_Length,0:DEGREE,1:6)   )
    ALLOCATE( zMB_Last_Col_Storage(1:LM_Length,0:DEGREE,1:6)    )

END IF


ALLOCATE( cVA_Load_Vector(1:NUM_R_NODES,1:LM_LENGTH,1:2)   )
ALLOCATE( cVB_Load_Vector(1:iVB_Prob_Dim,1:2) )

ALLOCATE( cVA_Coeff_Vector(1:NUM_R_NODES,1:LM_LENGTH,1:2) )
ALLOCATE( cVB_Coeff_Vector(1:iVB_Prob_Dim,1:2) )

IF ( FP_Diagnostics_Flag ) THEN
    ALLOCATE( FP_Update_Vector(1:NUM_R_NODES,1:LM_LENGTH,1:2) )
    ALLOCATE( FP_Laplace_Vector(1:NUM_R_NODES,1:LM_LENGTH,1:2) )
    ALLOCATE( FP_Residual_Vector(1:NUM_R_NODES,1:LM_LENGTH,1:2) )
    
    ALLOCATE(FP_Iter_Load_Storage(1:NUM_R_NODES,1:LM_LENGTH) )
    
    ALLOCATE( FP_Iteration_Log(1:2) )
    
    ALLOCATE( Resid_Norms(1:3,1:LM_LENGTH,1:Max_Iterations,1:8) )
    ALLOCATE( Update_Norms(1:3,1:Max_Iterations,1:8) )
    
    Resid_Norms = 0.0_idp
    Update_Norms = 0.0_idp
END IF

lPF_Init_Flags(iPF_Init_Alloc_LinSys) = .TRUE.

END SUBROUTINE Allocate_XCFC_Linear_Systems











 !+102+####################################################!
!                                                           !
!          Deallocate_XCFC_Linear_Systems                   !
!                                                           !
 !#########################################################!
SUBROUTINE Deallocate_XCFC_Linear_Systems()

IF ( MATRIX_FORMAT == 'Full' ) THEN
    DEALLOCATE( Laplace_Matrix_Full )
    DEALLOCATE( Laplace_Matrix_Beta )
ELSEIF ( MATRIX_FORMAT == 'CCS' ) THEN
    DEALLOCATE( Laplace_Matrix_VAL )
    DEALLOCATE( Laplace_Matrix_ROW )
    DEALLOCATE( Laplace_Matrix_COL )

    DEALLOCATE( Laplace_Factored_VAL )
    DEALLOCATE( Laplace_Factored_ROW )
    DEALLOCATE( Laplace_Factored_COL )

    DEALLOCATE( iMB_IPIV )
    DEALLOCATE( zMB_Matrix_Banded )
    DEALLOCATE( zMB_Matrix_Diagonal )

    DEALLOCATE( zMA_First_Col_Storage )
    DEALLOCATE( zMA_Last_Col_Storage )

    DEALLOCATE( zMB_First_Col_Storage )
    DEALLOCATE( zMB_Last_Col_Storage )
    
END IF

DEALLOCATE( cVA_Load_Vector )
DEALLOCATE( cVB_Load_Vector )

DEALLOCATE( cVA_Coeff_Vector )
DEALLOCATE( cVB_Coeff_Vector )

IF ( FP_Diagnostics_Flag ) THEN
    DEALLOCATE( FP_Update_Vector )
    DEALLOCATE( FP_Residual_Vector )
    DEALLOCATE( FP_Laplace_Vector )
    
    IF ( ALLOCATED(FP_Iter_Matrix_Storage)) THEN
        DEALLOCATE( FP_Iter_Matrix_Storage )
    END IF
    
    DEALLOCATE(FP_Iter_Load_Storage)
END IF

lPF_Init_Flags(iPF_Init_Alloc_LinSys) = .FALSE.

END SUBROUTINE Deallocate_XCFC_Linear_Systems






 !+103+####################################################!
!                                                           !
!          Reallocate_XCFC_Linear_Systems                   !
!                                                           !
 !#########################################################!
SUBROUTINE Reallocate_XCFC_Linear_Systems()

CALL Deallocate_XCFC_Linear_Systems()
CALL Allocate_XCFC_Linear_Systems()

END SUBROUTINE Reallocate_XCFC_Linear_Systems



END MODULE Allocation_XCFC_Linear_Systems

