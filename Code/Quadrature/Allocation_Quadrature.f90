   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Allocation_Quadrature                                                        !##!
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

USE Poseidon_Message_Routines_Module, &
            ONLY :  Init_Message

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag

USE Variables_Quadrature, &
            ONLY :  Num_R_Quad_Points,      &
                    Num_T_Quad_Points,      &
                    Num_P_Quad_Points,      &
                    Num_TP_Quad_Points,     &
                    Int_R_Locations,        &
                    Int_R_Weights,          &
                    Int_T_Locations,        &
                    Int_T_Weights,          &
                    Int_P_Locations,        &
                    Int_P_Weights,          &
                    Int_TP_Weights,         &
                    Local_Node_Locations

USE Flags_Initialization_Module, &
            ONLY :  lPF_Init_Quad_Flags,    &
                    iPF_Init_Quad_Alloc



IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!                            Allocate_Mesh                                      !
!                                                                               !
!###############################################################################!
SUBROUTINE Allocate_Quadrature()

IF ( Verbose_Flag ) CALL Init_Message('Allocating Quadrature Variables.')

ALLOCATE(INT_R_LOCATIONS(1:NUM_R_QUAD_POINTS), INT_R_WEIGHTS(1:NUM_R_QUAD_POINTS))
ALLOCATE(INT_T_LOCATIONS(1:NUM_T_QUAD_POINTS), INT_T_WEIGHTS(1:NUM_T_QUAD_POINTS))
ALLOCATE(INT_P_LOCATIONS(1:NUM_P_QUAD_POINTS), INT_P_WEIGHTS(1:NUM_P_QUAD_POINTS))

ALLOCATE( INT_TP_WEIGHTS(1:NUM_TP_QUAD_POINTS) )


lPF_Init_Quad_Flags(iPF_Init_Quad_Alloc) = .TRUE.

END SUBROUTINE Allocate_Quadrature











!+102+##########################################################################!
!                                                                               !
!                           Deallocate_Mesh                                     !
!                                                                               !
!###############################################################################!
SUBROUTINE Deallocate_Quadrature()


DEALLOCATE(INT_R_LOCATIONS, INT_R_WEIGHTS)
DEALLOCATE(INT_T_LOCATIONS, INT_T_WEIGHTS)
DEALLOCATE(INT_P_LOCATIONS, INT_P_WEIGHTS)

DEALLOCATE(INT_TP_WEIGHTS )
!DEALLOCATE( LOCAL_NODE_LOCATIONS )

lPF_Init_Quad_Flags(iPF_Init_Quad_Alloc) = .FALSE.

END SUBROUTINE Deallocate_Quadrature





END MODULE Allocation_Quadrature


