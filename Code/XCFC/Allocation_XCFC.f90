   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Allocation_XCFC                                                                !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
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
            ONLY :  DEGREE,                 &
                    L_LIMIT

Use Variables_Derived, &
            ONLY :  Num_R_Nodes,            &
                    Num_R_Nodesp1,          &
                    LM_Length,              &
                    Var_Dim,                &
                    Prob_Dim,               &
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
                    Beta_IPIV,              &
                    Beta_Diagonals,         &
                    Beta_Bandwidth,         &
                    Beta_MVL_Banded,        &
                    Beta_MVL_Diagonal,      &
                    FP_Source_Vector_A,     &
                    FP_Source_Vector_B,     &
                    FP_Coeff_Vector_A,      &
                    FP_Coeff_Vector_B,      &
                    FP_Update_Vector,       &
                    FP_Laplace_Vector,      &
                    FP_Laplace_Vector_Beta, &
                    FP_Laplace_Vector_X,    &
                    FP_Residual_Vector,     &
                    FP_Residual_Vector_Beta,&
                    FP_Residual_Vector_X,   &
                    Matrix_Format,          &
                    Num_Matrices,           &
                    First_Column_Storage,   &
                    Last_Column_Storage,    &
                    First_Column_Beta_Storage,   &
                    Last_Column_Beta_Storage




IMPLICIT NONE


CONTAINS



!+101+###########################################################################!
!                                                                                !
!                            Allocate_Poseidon_Variables                         !
!                                                                                !
!################################################################################!
SUBROUTINE Allocate_XCFC()


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

    ALLOCATE( First_Column_Beta_Storage(1:LM_Length,0:DEGREE,1:6)   )
    ALLOCATE( Last_Column_Beta_Storage(1:LM_Length,0:DEGREE,1:6)    )

END IF


ALLOCATE( FP_Source_Vector_A(1:NUM_R_NODES,1:LM_LENGTH,1:2)   )
ALLOCATE( FP_Source_Vector_B(1:Beta_Prob_Dim,1:2) )

ALLOCATE( FP_Coeff_Vector_A(1:NUM_R_NODES,1:LM_LENGTH,1:2) )
ALLOCATE( FP_Coeff_Vector_B(1:Beta_Prob_Dim,1:2) )

ALLOCATE( FP_Update_Vector(1:NUM_R_NODES,1:LM_LENGTH,1:8)  )

ALLOCATE( FP_Laplace_Vector(1:NUM_R_NODES,1:LM_LENGTH,1:2)  )
ALLOCATE( FP_Laplace_Vector_Beta(1:Beta_Prob_Dim)  )
ALLOCATE( FP_Laplace_Vector_X(1:Beta_Prob_Dim)  )

ALLOCATE( FP_Residual_Vector(1:NUM_R_NODES,1:LM_LENGTH,1:5)  )
ALLOCATE( FP_Residual_Vector_Beta(1:Beta_Prob_Dim)  )
ALLOCATE( FP_Residual_Vector_X(1:Beta_Prob_Dim)  )

!ALLOCATE( FP_Coeff_Vector(1:Prob_Dim) )
!ALLOCATE( FP_Update_Vector(1:Prob_Dim) )


END SUBROUTINE Allocate_XCFC











!+102+###########################################################################!
!                                                                                !
!                           Deallocate_Poseidon_Variables                        !
!                                                                                !
!################################################################################!
SUBROUTINE Deallocate_XCFC()

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

DEALLOCATE( FP_Source_Vector_A )
DEALLOCATE( FP_Source_Vector_B )

DEALLOCATE( FP_Coeff_Vector_A )
DEALLOCATE( FP_Coeff_Vector_B )

DEALLOCATE( FP_Update_Vector )
DEALLOCATE( FP_Laplace_Vector )
DEALLOCATE( FP_Laplace_Vector_Beta )
DEALLOCATE( FP_Laplace_Vector_X )

DEALLOCATE( FP_Residual_Vector )
DEALLOCATE( FP_Residual_Vector_Beta )
DEALLOCATE( FP_Residual_Vector_X )




END SUBROUTINE Deallocate_XCFC







END MODULE Allocation_XCFC

