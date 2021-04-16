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






IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!                            Allocate_Mesh                                      !
!                                                                               !
!###############################################################################!
SUBROUTINE Allocate_Mesh()



ALLOCATE(rlocs(0:NUM_R_ELEMENTS))
ALLOCATE(tlocs(0:NUM_T_ELEMENTS))
ALLOCATE(plocs(0:NUM_P_ELEMENTS))

ALLOCATE(drlocs(0:NUM_R_ELEMENTS-1))
ALLOCATE(dtlocs(0:NUM_T_ELEMENTS-1))
ALLOCATE(dplocs(0:NUM_P_ELEMENTS-1))




END SUBROUTINE Allocate_Mesh











!+102+##########################################################################!
!                                                                               !
!                           Deallocate_Mesh                                     !
!                                                                               !
!###############################################################################!
SUBROUTINE Deallocate_Mesh()


DEALLOCATE(rlocs)
DEALLOCATE(tlocs)
DEALLOCATE(plocs)


DEALLOCATE(drlocs)
DEALLOCATE(dtlocs)
DEALLOCATE(dplocs)

END SUBROUTINE Deallocate_Mesh





END MODULE Allocation_Mesh

