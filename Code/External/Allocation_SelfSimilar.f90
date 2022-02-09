   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Allocation_SelfSimilar                                                       !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   Allocate_Yahil_Variables                                            !##!
!##!    +102+   Deallocate_Yahil_Variables                                          !##!
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
            ONLY :  idp

USE Variables_Mesh, &
            ONLY :  Num_R_Elements

USE Variables_Yahil, &
            ONLY :  SELFSIM_R_VALS,     &
                    SELFSIM_POT_VALS,   &
                    SELFSIM_SHIFT_VALs, &
                    SelfSim_Allocated




IMPLICIT NONE


CONTAINS



!+101+###########################################################################!
!                                                                                !
!                            Allocate_Poseidon_Variables                         !
!                                                                                !
!################################################################################!
SUBROUTINE Allocate_SelfSim(Num_Entries)

INTEGER, INTENT(IN)                             ::  Num_Entries


IF ( SelfSim_Allocated ) THEN
    WRITE(*,'(A)')'Warning attempting to reallocate Selfsim Variables.'
ELSE
    ALLOCATE( SELFSIM_R_VALS(0:NUM_ENTRIES) )
    ALLOCATE( SELFSIM_POT_VALS(0:NUM_ENTRIES) )
    ALLOCATE( SELFSIM_SHIFT_VALS(0:Num_R_Elements+1) )
    SelfSim_Allocated = .TRUE.
END IF

END SUBROUTINE Allocate_SelfSim











!+102+###########################################################################!
!                                                                                !
!                           Deallocate_Poseidon_Variables                        !
!                                                                                !
!################################################################################!
SUBROUTINE Deallocate_SelfSim()

IF ( SelfSim_Allocated ) THEN

    DEALLOCATE( SELFSIM_R_VALS )
    DEALLOCATE( SELFSIM_POT_VALS )
    DEALLOCATE( SELFSIM_SHIFT_VALS )
    SelfSim_Allocated = .FALSE.

END IF

END SUBROUTINE Deallocate_SelfSim







END MODULE Allocation_SelfSimilar


