   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Allocation_Mesh                                                              !##!
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

USE Variables_Mesh, &
            ONLY :  Num_R_Elements,         &
                    Num_T_Elements,         &
                    Num_P_Elements,         &
                    rlocs,                  &
                    tlocs,                  &
                    plocs,                  &
                    drlocs,                 &
                    dtlocs,                 &
                    dplocs


USE Flags_Initialization_Module, &
            ONLY :  lPF_Init_Mesh_Flags,    &
                    iPF_Init_Mesh_Alloc



IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!                            Allocate_Mesh                                      !
!                                                                               !
!###############################################################################!
SUBROUTINE Allocate_Mesh()

IF ( Verbose_Flag ) CALL Init_Message('Allocating Mesh Variables.')

PRINT*,Num_R_Elements
PRINT*,Num_T_Elements
PRint*,Num_P_Elements

ALLOCATE(rlocs(0:NUM_R_ELEMENTS))
ALLOCATE(tlocs(0:NUM_T_ELEMENTS))
ALLOCATE(plocs(0:NUM_P_ELEMENTS))

ALLOCATE(drlocs(0:NUM_R_ELEMENTS-1))
ALLOCATE(dtlocs(0:NUM_T_ELEMENTS-1))
ALLOCATE(dplocs(0:NUM_P_ELEMENTS-1))

lPF_Init_Mesh_Flags(iPF_Init_Mesh_Alloc) = .TRUE.


END SUBROUTINE Allocate_Mesh





!+102+##########################################################################!
!                                                                               !
!                           Deallocate_Mesh                                     !
!                                                                               !
!###############################################################################!
SUBROUTINE Deallocate_Mesh()

IF ( Verbose_Flag ) CALL Init_Message('Deallocating Mesh Variables.')

DEALLOCATE(rlocs)
DEALLOCATE(tlocs)
DEALLOCATE(plocs)


DEALLOCATE(drlocs)
DEALLOCATE(dtlocs)
DEALLOCATE(dplocs)

lPF_Init_Mesh_Flags(iPF_Init_Mesh_Alloc) = .FALSE.

END SUBROUTINE Deallocate_Mesh




!+102+##########################################################################!
!                                                                               !
!                           Reallocate_Mesh                                     !
!                                                                               !
!###############################################################################!
SUBROUTINE Reallocate_Mesh()

CALL Deallocate_Mesh()
CALL Allocate_Mesh()

END SUBROUTINE Reallocate_Mesh


END MODULE Allocation_Mesh

