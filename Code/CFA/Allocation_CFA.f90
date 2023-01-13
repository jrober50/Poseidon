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
                    iVB_Prob_Dim


USE Variables_Vectors, &
            ONLY :  dVA_Load_Vector,      &
                    dVB_Load_Vector,      &
                    dVA_Coeff_Vector,      &
                    dVB_Coeff_Vector
                    

USE Variables_Matrices, &
            ONLY :  Matrix_Format,          &
                    dMA_First_Col_Storage,   &
                    dMA_Last_Col_Storage,    &
                    dMB_First_Col_Storage,   &
                    dMB_Last_Col_Storage,   &
                    Laplace_Matrix_Full,    &
                    Laplace_Matrix_Beta,    &
                    Laplace_Matrix_VAL,     &
                    Laplace_Matrix_ROW,     &
                    Laplace_Matrix_COL,     &
                    Laplace_Factored_VAL,   &
                    Laplace_Factored_ROW,   &
                    Laplace_Factored_COL,   &
                    Laplace_NNZ,            &
                    iMB_IPIV,               &
                    iMB_Diagonals,          &
                    dMB_Matrix_Banded,        &
                    dMB_Matrix_Diagonal

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
    ALLOCATE( Laplace_Matrix_Beta(1:iVB_Prob_Dim,1:iVB_Prob_Dim) )


ELSEIF ( MATRIX_FORMAT == 'CCS' ) THEN

    ALLOCATE( Laplace_Matrix_VAL(0:Laplace_NNZ-1, 0:L_LIMIT) )
    ALLOCATE( Laplace_Matrix_ROW(0:Laplace_NNZ-1, 0:L_LIMIT) )
    ALLOCATE( Laplace_Matrix_COL(0:NUM_R_NODES, 0:L_LIMIT) )

    ALLOCATE( Laplace_Factored_VAL(0:Laplace_NNZ-1, 0:L_LIMIT) )
    ALLOCATE( Laplace_Factored_ROW(0:Laplace_NNZ-1, 0:L_LIMIT) )
    ALLOCATE( Laplace_Factored_COL(0:NUM_R_NODES, 0:L_LIMIT) )

    ALLOCATE( iMB_IPIV(1:iVB_Prob_Dim) )
    ALLOCATE( dMB_Matrix_Banded(1:(3*iMB_Diagonals+1), 1:iVB_Prob_Dim))
    ALLOCATE( dMB_Matrix_Diagonal(1:iVB_Prob_Dim) )

    ALLOCATE( dMA_First_Col_Storage(0:DEGREE,0:L_LIMIT)   )
    ALLOCATE( dMA_Last_Col_Storage(0:DEGREE,0:L_LIMIT)    )

    ALLOCATE( dMB_First_Col_Storage(1:LM_Length,0:DEGREE,1:3)   )
    ALLOCATE( dMB_Last_Col_Storage(1:LM_Length,0:DEGREE,1:3)    )

END IF


ALLOCATE( dVA_Load_Vector(1:NUM_R_NODES,1:LM_LENGTH,1:2)   )
ALLOCATE( dVB_Load_Vector(1:iVB_Prob_Dim,1:2) )

ALLOCATE( dVA_Coeff_Vector(1:NUM_R_NODES,1:LM_LENGTH,1:5)         )
ALLOCATE( dVB_Coeff_Vector(1:iVB_Prob_Dim,1:2) )

ALLOCATE( FP_Update_Vector(1:NUM_R_NODES,1:LM_LENGTH,1:5)  )

ALLOCATE( FP_Laplace_Vector(1:NUM_R_NODES,1:LM_LENGTH,1:2)  )
ALLOCATE( FP_Laplace_Vector_Beta(1:iVB_Prob_Dim)  )

ALLOCATE( FP_Residual_Vector(1:NUM_R_NODES,1:LM_LENGTH,1:5)  )
ALLOCATE( FP_Residual_Vector_Beta(1:iVB_Prob_Dim)  )

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

    DEALLOCATE( iMB_IPIV )
    DEALLOCATE( dMB_Matrix_Banded )
    DEALLOCATE( dMB_Matrix_Diagonal )

    DEALLOCATE( dMA_First_Col_Storage )
    DEALLOCATE( dMA_Last_Col_Storage )

    DEALLOCATE( dMB_First_Col_Storage )
    DEALLOCATE( dMB_Last_Col_Storage )
    
END IF

DEALLOCATE( dVA_Load_Vector )
DEALLOCATE( dVB_Load_Vector )
DEALLOCATE( dVA_Coeff_Vector )
DEALLOCATE( dVB_Coeff_Vector )
DEALLOCATE( FP_Update_Vector )
DEALLOCATE( FP_Laplace_Vector )
DEALLOCATE( FP_Laplace_Vector_Beta )
DEALLOCATE( FP_Residual_Vector )
DEALLOCATE( FP_Residual_Vector_Beta )





END SUBROUTINE Deallocate_CFA_Linear_Systems







END MODULE Allocation_CFA_Linear_Systems

