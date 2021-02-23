   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Allocation_FP                                                                !##!
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
USE Poseidon_Kinds_Module, &
            ONLY : idp

USE Poseidon_Parameters, &
            ONLY :  DOMAIN_DIM,             &
                    DEGREE,                 &
                    L_LIMIT,                &
                    NUM_CFA_Eqs

Use Variables_Derived, &
            ONLY :  Num_R_Nodes,            &
                    Num_R_Nodesp1,          &
                    LM_Length,              &
                    Var_Dim,                &
                    Beta_Prob_Dim


USE Variables_FP, &
            ONLY :  Laplace_Matrix_Full,    &
                    Laplace_Matrix_Beta,    &
                    Laplace_Matrix_VAL,     &
                    Laplace_Matrix_ROW,     &
                    Laplace_Matrix_COL,     &
                    Laplace_Factored_VAL,   &
                    Laplace_Factored_ROW,   &
                    Laplace_Factored_COL,   &
                    Laplace_NNZ,            &
                    FP_Source_Vector,       &
                    FP_Source_Vector_Beta,  &
                    FP_Coeff_Vector,        &
                    FP_Coeff_Vector_Beta,   &
                    FP_Update_Vector,       &
                    FP_Laplace_Vector,      &
                    FP_Laplace_Vector_Beta, &
                    FP_Residual_Vector,     &
                    FP_Residual_Vector_Beta,&
                    Matrix_Format,          &
                    Num_Matrices,           &
                    First_Column_Storage,   &
                    Last_Column_Storage




IMPLICIT NONE


CONTAINS



!+101+###########################################################################!
!                                                                                !
!                            Allocate_Poseidon_Variables                         !
!                                                                                !
!################################################################################!
SUBROUTINE Allocate_FP()

IF ( MATRIX_FORMAT == 'Full' ) THEN

    ALLOCATE( Laplace_Matrix_Full(1:NUM_R_NODES,1:NUM_R_NODES,0:L_LIMIT) )
    ALLOCATE( Laplace_Matrix_Beta(1:Beta_Prob_Dim,1:Beta_Prob_Dim) )


ELSEIF ( MATRIX_FORMAT == 'CCS' ) THEN

    PRINT*,"Num_Matrices",Num_Matrices,Num_R_Nodes, L_LIMIT

    ALLOCATE( Laplace_Matrix_VAL(0:Laplace_NNZ-1, 0:L_LIMIT,1:Num_Matrices) )
    ALLOCATE( Laplace_Matrix_ROW(0:Laplace_NNZ-1, 0:L_LIMIT) )
    ALLOCATE( Laplace_Matrix_COL(0:NUM_R_NODES, 0:L_LIMIT) )

    ALLOCATE( Laplace_Factored_VAL(0:Laplace_NNZ-1, 0:L_LIMIT,1:Num_Matrices) )
    ALLOCATE( Laplace_Factored_ROW(0:Laplace_NNZ-1, 0:L_LIMIT) )
    ALLOCATE( Laplace_Factored_COL(0:NUM_R_NODES, 0:L_LIMIT) )

    ALLOCATE( First_Column_Storage(0:DEGREE,0:L_LIMIT,1:Num_Matrices)   )
    ALLOCATE( Last_Column_Storage(0:DEGREE,0:L_LIMIT,1:Num_Matrices)    )

END IF


ALLOCATE( FP_Source_Vector(1:NUM_R_NODES,0:LM_LENGTH-1,1:2)   )
ALLOCATE( FP_Source_Vector_Beta(1:Beta_Prob_Dim) )
ALLOCATE( FP_Coeff_Vector(1:NUM_R_NODES,0:LM_LENGTH-1,1:5)        )
ALLOCATE( FP_Coeff_Vector_Beta(1:Beta_Prob_Dim) )
ALLOCATE( FP_Update_Vector(1:NUM_R_NODES,0:LM_LENGTH-1,1:5)   )
ALLOCATE( FP_Laplace_Vector(1:NUM_R_NODES,0:LM_LENGTH-1,1:2)  )
ALLOCATE( FP_Laplace_Vector_Beta(1:Beta_Prob_Dim)  )
ALLOCATE( FP_Residual_Vector(1:NUM_R_NODES,0:LM_LENGTH-1,1:2) )
ALLOCATE( FP_Residual_Vector_Beta(1:Beta_Prob_Dim)  )





END SUBROUTINE Allocate_FP











!+102+###########################################################################!
!                                                                                !
!                           Deallocate_Poseidon_Variables                        !
!                                                                                !
!################################################################################!
SUBROUTINE Deallocate_FP()

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

    DEALLOCATE( First_Column_Storage )
    DEALLOCATE( Last_Column_Storage )
    
END IF

DEALLOCATE( FP_Source_Vector )
DEALLOCATE( FP_Source_Vector_Beta )
DEALLOCATE( FP_Coeff_Vector )
DEALLOCATE( FP_Coeff_Vector_Beta )
DEALLOCATE( FP_Update_Vector )
DEALLOCATE( FP_Laplace_Vector )
DEALLOCATE( FP_Laplace_Vector_Beta )
DEALLOCATE( FP_Residual_Vector )
DEALLOCATE( FP_Residual_Vector_Beta )





END SUBROUTINE Deallocate_FP







END MODULE Allocation_FP

