   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_FP_Allocation_Module                                                !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains the functions and subroutines associated with the stiffness matrix !##!
!##!        involved in the linear system formed by the expansions. This includes   !##!
!##!        functions that build the matrix in various storage formats.             !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   Allocate_Poseidon_FP_Variables                                      !##!
!##!    +102+   Deallocate_Poseidon_FP_Variables                                    !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!



!*D*================================!
!                                   !
!           Dependencies            !
!                                   !
!===================================!
USE Poseidon_Constants_Module, &
            ONLY : idp

USE Poseidon_Parameters, &
            ONLY :  DOMAIN_DIM,             &
                    DEGREE,                 &
                    L_LIMIT,                &
                    NUM_CFA_Eqs

Use Poseidon_Variables_Module, &
            ONLY :  Num_R_Nodes,            &
                    Num_R_Nodesp1,          &
                    LM_Length


USE Poseidon_FP_Variables_Module, &
            ONLY :  Laplace_Matrix_Full,    &
                    Laplace_Matrix_VAL,     &
                    Laplace_Matrix_ROW,     &
                    Laplace_Matrix_COL,     &
                    Laplace_Factored_VAL,   &
                    Laplace_Factored_ROW,   &
                    Laplace_Factored_COL,   &
                    Laplace_NNZ,            &
                    FP_Source_Vector,       &
                    FP_Coeff_Vector,  &
                    FP_Update_Vector,       &
                    FP_Laplace_Vector,      &
                    FP_Residual_Vector,     &
                    Matrix_Format,          &
                    Num_Matrices




IMPLICIT NONE


CONTAINS



!+101+###########################################################################!
!                                                                                !
!                            Allocate_Poseidon_Variables                         !
!                                                                                !
!################################################################################!
SUBROUTINE Allocate_Poseidon_FP_Variables()

IF ( MATRIX_FORMAT == 'FULL' ) THEN

    ALLOCATE( Laplace_Matrix_Full(1:NUM_R_NODES,1:NUM_R_NODES,0:L_LIMIT,1:Num_Matrices) )

ELSEIF ( MATRIX_FORMAT == 'CCS' ) THEN

    ALLOCATE( Laplace_Matrix_VAL(1:Laplace_NNZ, 0:L_LIMIT,1:Num_Matrices) )
    ALLOCATE( Laplace_Matrix_ROW(1:Laplace_NNZ, 0:L_LIMIT) )
    ALLOCATE( Laplace_Matrix_COL(1:NUM_R_NODESp1, 0:L_LIMIT) )

    ALLOCATE( Laplace_Factored_VAL(1:Laplace_NNZ, 0:L_LIMIT,1:Num_Matrices) )
    ALLOCATE( Laplace_Factored_ROW(1:Laplace_NNZ, 0:L_LIMIT) )
    ALLOCATE( Laplace_Factored_COL(1:NUM_R_NODESp1, 0:L_LIMIT) )

END IF


ALLOCATE( FP_Source_Vector(1:NUM_R_NODESp1,0:LM_LENGTH,1:NUM_CFA_EQS)   )
ALLOCATE( FP_Coeff_Vector(1:NUM_R_NODESp1,0:LM_LENGTH,1:5)        )
ALLOCATE( FP_Update_Vector(1:NUM_R_NODESp1,0:LM_LENGTH,1:NUM_CFA_EQS)   )
ALLOCATE( FP_Laplace_Vector(1:NUM_R_NODESp1,0:LM_LENGTH,1:NUM_CFA_EQS)  )
ALLOCATE( FP_Residual_Vector(1:NUM_R_NODESp1,0:LM_LENGTH,1:NUM_CFA_EQS) )

END SUBROUTINE Allocate_Poseidon_FP_Variables











!+102+###########################################################################!
!                                                                                !
!                           Deallocate_Poseidon_Variables                        !
!                                                                                !
!################################################################################!
SUBROUTINE Deallocate_Poseidon_FP_Variables()

DEALLOCATE( Laplace_Matrix_Full )

DEALLOCATE( Laplace_Matrix_VAL )
DEALLOCATE( Laplace_Matrix_ROW )
DEALLOCATE( Laplace_Matrix_COL )

DEALLOCATE( FP_Source_Vector )
DEALLOCATE( FP_Coeff_Vector )
DEALLOCATE( FP_Update_Vector )

END SUBROUTINE Deallocate_Poseidon_FP_Variables







END MODULE Poseidon_FP_Allocation_Module

