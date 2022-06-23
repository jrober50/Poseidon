   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Allocation_CFA_Linear_Systems                                                !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   Allocate_CFA_Linear_Systems                                         !##!
!##!    +102+   Deallocate_CFA_Linear_Systems                                       !##!
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
            ONLY :  Degree,                 &
                    L_Limit

Use Variables_Derived, &
            ONLY :  Num_R_Nodes,            &
                    Num_R_Nodesp1,          &
                    LM_Length,              &
                    Var_Dim,                &
                    Prob_Dim,               &
                    Beta_Prob_Dim


USE Variables_Vectors, &
            ONLY :  cVA_Source_Vector,      &
                    cVB_Source_Vector,      &
                    cVA_Coeff_Vector,      &
                    cVB_Coeff_Vector
                    

USE Variables_Matrices, &
            ONLY :  Num_Matrices,           &
                    Matrix_Format,          &
                    First_Column_Storage,   &
                    Last_Column_Storage,    &
                    First_Column_Beta_Storage,   &
                    Last_Column_Beta_Storage,   &
                    Laplace_Matrix_Full,    &
                    Laplace_Matrix_Beta,    &
                    Laplace_Matrix_VAL,     &
                    Laplace_Matrix_ROW,     &
                    Laplace_Matrix_COL,     &
                    Laplace_Factored_VAL,   &
                    Laplace_Factored_ROW,   &
                    Laplace_Factored_COL,   &
                    Laplace_NNZ,            &
                    Beta_IPIV,              &
                    Beta_Diagonals,         &
                    Beta_Bandwidth,         &
                    Beta_MVL_Banded,        &
                    Beta_MVL_Diagonal

USE Variables_FP, &
            ONLY :  FP_Update_Vector,       &
                    FP_Laplace_Vector,      &
                    FP_Laplace_Vector_Beta, &
                    FP_Residual_Vector,     &
                    FP_Residual_Vector_Beta
                    




IMPLICIT NONE


CONTAINS



!+101+###########################################################################!
!                                                                                !
!                            Allocate_Poseidon_Variables                         !
!                                                                                !
!################################################################################!
SUBROUTINE Allocate_CFA_Linear_Systems

IF ( MATRIX_FORMAT == 'Full' ) THEN

    ALLOCATE( Laplace_Matrix_Full(1:NUM_R_NODES,1:NUM_R_NODES,0:L_LIMIT) )
    ALLOCATE( Laplace_Matrix_Beta(1:Beta_Prob_Dim,1:Beta_Prob_Dim) )


ELSEIF ( MATRIX_FORMAT == 'CCS' ) THEN

    ALLOCATE( Laplace_Matrix_VAL(0:Laplace_NNZ-1, 0:L_LIMIT,1:Num_Matrices) )
    ALLOCATE( Laplace_Matrix_ROW(0:Laplace_NNZ-1, 0:L_LIMIT) )
    ALLOCATE( Laplace_Matrix_COL(0:NUM_R_NODES, 0:L_LIMIT) )

    ALLOCATE( Laplace_Factored_VAL(0:Laplace_NNZ-1, 0:L_LIMIT,1:Num_Matrices) )
    ALLOCATE( Laplace_Factored_ROW(0:Laplace_NNZ-1, 0:L_LIMIT) )
    ALLOCATE( Laplace_Factored_COL(0:NUM_R_NODES, 0:L_LIMIT) )

    ALLOCATE( Beta_IPIV(1:Beta_Prob_Dim) )
    ALLOCATE( Beta_MVL_Banded(1:(3*Beta_Diagonals+1), 1:Beta_Prob_Dim))
    ALLOCATE( Beta_MVL_Diagonal(1:Beta_Prob_Dim) )

    ALLOCATE( First_Column_Storage(0:DEGREE,0:L_LIMIT,1:Num_Matrices)   )
    ALLOCATE( Last_Column_Storage(0:DEGREE,0:L_LIMIT,1:Num_Matrices)    )

    ALLOCATE( First_Column_Beta_Storage(1:LM_Length,0:DEGREE,1:3)   )
    ALLOCATE( Last_Column_Beta_Storage(1:LM_Length,0:DEGREE,1:3)    )

END IF


ALLOCATE( cVA_Source_Vector(1:NUM_R_NODES,1:LM_LENGTH,1:2)   )
ALLOCATE( cVB_Source_Vector(1:Beta_Prob_Dim,1:2) )

ALLOCATE( cVA_Coeff_Vector(1:NUM_R_NODES,1:LM_LENGTH,1:5)         )
ALLOCATE( cVB_Coeff_Vector(1:Beta_Prob_Dim,1:2) )

ALLOCATE( FP_Update_Vector(1:NUM_R_NODES,1:LM_LENGTH,1:5)  )

ALLOCATE( FP_Laplace_Vector(1:NUM_R_NODES,1:LM_LENGTH,1:2)  )
ALLOCATE( FP_Laplace_Vector_Beta(1:Beta_Prob_Dim)  )

ALLOCATE( FP_Residual_Vector(1:NUM_R_NODES,1:LM_LENGTH,1:5)  )
ALLOCATE( FP_Residual_Vector_Beta(1:Beta_Prob_Dim)  )

!ALLOCATE( FP_Coeff_Vector(1:Prob_Dim) )
!ALLOCATE( FP_Update_Vector(1:Prob_Dim) )


END SUBROUTINE Allocate_CFA_Linear_Systems











!+102+###########################################################################!
!                                                                                !
!                           Deallocate_Poseidon_Variables                        !
!                                                                                !
!################################################################################!
SUBROUTINE Deallocate_CFA_Linear_Systems()

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

    DEALLOCATE( Beta_IPIV )
    DEALLOCATE( Beta_MVL_Banded )
    DEALLOCATE( Beta_MVL_Diagonal )

    DEALLOCATE( First_Column_Storage )
    DEALLOCATE( Last_Column_Storage )

    DEALLOCATE( First_Column_Beta_Storage )
    DEALLOCATE( Last_Column_Beta_Storage )
    
END IF

DEALLOCATE( cVA_Source_Vector )
DEALLOCATE( cVB_Source_Vector )
DEALLOCATE( cVA_Coeff_Vector )
DEALLOCATE( cVB_Coeff_Vector )
DEALLOCATE( FP_Update_Vector )
DEALLOCATE( FP_Laplace_Vector )
DEALLOCATE( FP_Laplace_Vector_Beta )
DEALLOCATE( FP_Residual_Vector )
DEALLOCATE( FP_Residual_Vector_Beta )





END SUBROUTINE Deallocate_CFA_Linear_Systems







END MODULE Allocation_CFA_Linear_Systems

