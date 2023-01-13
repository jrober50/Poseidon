   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Allocation_Poisson_Linear_System                                      !##!
!##!                                                                         !##!
!##!_________________________________________________________________________!##!
!##!                                                                         !##!
!##!                                                                         !##!
!##!=========================================================================!##!
!##!                                                                         !##!
!##!    Contains:                                                            !##!
!##!                                                                         !##!
!##!                                                                         !##!
!###############################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !##########################################################################!


!*D*================================!
!                                   !
!           Dependencies            !
!                                   !
!===================================!
USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Poseidon_Parameters, &
            ONLY :  Degree,             &
                    L_Limit

USE Variables_Vectors, &
            ONLY :  dVA_Coeff_Vector,   &
                    dVA_Load_Vector

USE Variables_Mesh, &
            ONLY :  Num_R_Elements,     &
                    Num_T_Elements,     &
                    Num_P_Elements

USE Variables_FP,  &
            ONLY :  FP_Update_Vector

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
                    Laplace_NNZ,                &
                    dMA_First_Col_Storage,       &
                    dMA_Last_Col_Storage

USE Variables_Quadrature, &
            ONLY :  Local_Quad_DOF

USE Variables_Derived, &
            ONLY :  LM_Length,          &
                    Num_R_Nodes

USE Flags_Initialization_Module, &
            ONLY :  lPF_Init_Flags,         &
                    iPF_Init_Alloc_LinSys

IMPLICIT NONE


CONTAINS



 !+101+############################################################!
!                                                                   !
!          Allocate_Poisson_Linear_System            			    !
!                                                                   !
 !#################################################################!
SUBROUTINE Allocate_Poisson_Linear_System



IF ( MATRIX_FORMAT == 'Full' ) THEN

    ALLOCATE( Laplace_Matrix_Full(1:NUM_R_NODES,1:NUM_R_NODES,0:L_LIMIT) )


ELSEIF ( MATRIX_FORMAT == 'CCS' ) THEN

    ALLOCATE( Laplace_Matrix_VAL(0:Laplace_NNZ-1, 0:L_LIMIT) )
    ALLOCATE( Laplace_Matrix_ROW(0:Laplace_NNZ-1, 0:L_LIMIT) )
    ALLOCATE( Laplace_Matrix_COL(0:NUM_R_NODES, 0:L_LIMIT) )

    ALLOCATE( Laplace_Factored_VAL(0:Laplace_NNZ-1, 0:L_LIMIT) )
    ALLOCATE( Laplace_Factored_ROW(0:Laplace_NNZ-1, 0:L_LIMIT) )
    ALLOCATE( Laplace_Factored_COL(0:NUM_R_NODES, 0:L_LIMIT) )

    ALLOCATE( dMA_First_Col_Storage(0:DEGREE,0:L_LIMIT)   )
    ALLOCATE( dMA_Last_Col_Storage(0:DEGREE,0:L_LIMIT)    )


END IF


ALLOCATE( dVA_Load_Vector(1:NUM_R_NODES,1:LM_LENGTH,1)   )
ALLOCATE( dVA_Coeff_Vector(1:NUM_R_NODES,1:LM_LENGTH,1) )


ALLOCATE( FP_Update_Vector(1:NUM_R_NODES,1:LM_LENGTH,1:2) )

lPF_Init_Flags(iPF_Init_Alloc_LinSys) = .TRUE.




END SUBROUTINE Allocate_Poisson_Linear_System





!+101+##########################################################################!
!                                                                               !
!                                                                     !
!                                                                               !
!###############################################################################!
SUBROUTINE Deallocate_Poisson_Linear_System

IF ( MATRIX_FORMAT == 'Full' ) THEN
    DEALLOCATE( Laplace_Matrix_Full )
ELSEIF ( MATRIX_FORMAT == 'CCS' ) THEN
    DEALLOCATE( Laplace_Matrix_VAL )
    DEALLOCATE( Laplace_Matrix_ROW )
    DEALLOCATE( Laplace_Matrix_COL )

    DEALLOCATE( Laplace_Factored_VAL )
    DEALLOCATE( Laplace_Factored_ROW )
    DEALLOCATE( Laplace_Factored_COL )

    DEALLOCATE( dMA_First_Col_Storage )
    DEALLOCATE( dMA_Last_Col_Storage )

END IF

DEALLOCATE( dVA_Load_Vector )
DEALLOCATE( dVA_Coeff_Vector )

DEALLOCATE( FP_Update_Vector )

lPF_Init_Flags(iPF_Init_Alloc_LinSys) = .FALSE.

END SUBROUTINE Deallocate_Poisson_Linear_System








END MODULE Allocation_Poisson_Linear_System
