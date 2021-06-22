   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Allocation_NR                                                              !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
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
USE Poseidon_Kinds_Module, &
            ONLY : idp



USE Variables_NR, &
            ONLY :  NR_Update_Vector,              &
                    NR_Coeff_Vector,         &
                    Block_RHS_Vector,           &
                    Block_STF_Mat,              &
                    Block_Elem_STF_MatVec

USE Variables_Derived, &
            ONLY :  Prob_Dim,                   &
                    Block_Prob_Dim,             &
                    SubShell_Prob_Dim,          &
                    Elem_Prob_Dim_Sqr,          &
                    Num_Off_Diagonals

USE Variables_MPI,  &
            ONLY :  Num_R_Elems_Per_Block



IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!                            Allocate_Mesh                                      !
!                                                                               !
!###############################################################################!
SUBROUTINE Allocate_NR()


ALLOCATE( NR_Update_Vector(0:PROB_DIM-1 ) )
ALLOCATE( NR_Coeff_Vector(0:PROB_DIM-1 ) )


ALLOCATE( Block_RHS_Vector( 1:Block_PROB_DIM ) )
ALLOCATE( Block_STF_Mat( 0:2*NUM_OFF_DIAGONALS, 0:SUBSHELL_PROB_DIM-1) )

ALLOCATE( Block_Elem_STF_MatVec(0:ELEM_PROB_DIM_SQR-1 ,0:NUM_R_ELEMS_PER_BLOCK-1) )


END SUBROUTINE Allocate_NR











!+102+##########################################################################!
!                                                                               !
!                           Deallocate_Mesh                                     !
!                                                                               !
!###############################################################################!
SUBROUTINE Deallocate_NR()


DEALLOCATE( NR_Update_Vector )
DEALLOCATE( NR_Coeff_Vector )

DEALLOCATE( Block_RHS_Vector )
DEALLOCATE( Block_STF_Mat )

DEALLOCATE( Block_Elem_STF_MatVec )

END SUBROUTINE Deallocate_NR





END MODULE Allocation_NR


