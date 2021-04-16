   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Allocation_Tables                                                            !##!
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

USE Poseidon_Parameters, &
            ONLY :  Degree,                 &
                    L_Limit

USE Variables_Derived, &
            ONLY :  LM_Length

USE Variables_MPI, &
            ONLY :  Num_R_Elems_Per_Block,  &
                    Num_T_Elems_Per_Block,  &
                    Num_P_Elems_Per_Block

USE Variables_Quadrature, &
            ONLY :  Num_R_Quad_Points,      &
                    Num_T_Quad_Points,      &
                    Num_P_Quad_Points,      &
                    Num_TP_Quad_Points

USE Variables_Tables, &
            ONLY :  Ylm_Table_Block,        &
                    Ylm_Values,             &
                    Ylm_dt_Values,          &
                    Ylm_dp_Values,          &
                    Ylm_CC_Values,          &
                    Ylm_CC_dt_Values,       &
                    Ylm_CC_dp_Values,       &
                    Lagrange_Poly_Table,    &
                    LPT_LPT,                &
                    M_Values




IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!                            Allocate_Mesh                                      !
!                                                                               !
!###############################################################################!
SUBROUTINE Allocate_Tables()



ALLOCATE( M_VALUES(0:L_LIMIT) )



ALLOCATE( Ylm_Values(       1:LM_Length,                                            &
                            1:NUM_TP_QUAD_POINTS,                                   &
                            0:NUM_T_ELEMS_PER_BLOCK-1, 0:NUM_P_ELEMS_PER_BLOCK-1)   )

ALLOCATE( Ylm_dt_Values(    1:LM_Length,                                            &
                            1:NUM_TP_QUAD_POINTS,                                   &
                            0:NUM_T_ELEMS_PER_BLOCK-1, 0:NUM_P_ELEMS_PER_BLOCK-1)   )

ALLOCATE( Ylm_dp_Values(    1:LM_Length,                                            &
                            1:NUM_TP_QUAD_POINTS,                                   &
                            0:NUM_T_ELEMS_PER_BLOCK-1, 0:NUM_P_ELEMS_PER_BLOCK-1)   )

ALLOCATE( Ylm_CC_Values(    1:NUM_TP_QUAD_POINTS,                                   &
                            1:LM_Length,                                            &
                            0:NUM_T_ELEMS_PER_BLOCK-1, 0:NUM_P_ELEMS_PER_BLOCK-1)   )

ALLOCATE( Ylm_CC_dt_Values( 1:NUM_TP_QUAD_POINTS,                                   &
                            1:LM_Length,                                            &
                            0:NUM_T_ELEMS_PER_BLOCK-1, 0:NUM_P_ELEMS_PER_BLOCK-1)   )

ALLOCATE( Ylm_CC_dp_Values( 1:NUM_TP_QUAD_POINTS,                                   &
                            1:LM_Length,                                            &
                            0:NUM_T_ELEMS_PER_BLOCK-1, 0:NUM_P_ELEMS_PER_BLOCK-1)   )


ALLOCATE( Lagrange_Poly_Table(0:DEGREE, 1:NUM_R_QUAD_POINTS, 0:2)   )
ALLOCATE( LPT_LPT( 1:NUM_R_QUAD_POINTS,0:DEGREE,0:DEGREE,0:1,0:2)       )


END SUBROUTINE Allocate_Tables











!+102+##########################################################################!
!                                                                               !
!                           Deallocate_Mesh                                     !
!                                                                               !
!###############################################################################!
SUBROUTINE Deallocate_Tables()


DEALLOCATE( Ylm_Values )
DEALLOCATE( Ylm_dt_Values )
DEALLOCATE( Ylm_dp_Values )

DEALLOCATE( Ylm_CC_Values )
DEALLOCATE( Ylm_CC_dt_Values )
DEALLOCATE( Ylm_CC_dp_Values )

DEALLOCATE( M_Values )

DEALLOCATE( Lagrange_Poly_Table )
DEALLOCATE( LPT_LPT )

END SUBROUTINE Deallocate_Tables





END MODULE Allocation_Tables


