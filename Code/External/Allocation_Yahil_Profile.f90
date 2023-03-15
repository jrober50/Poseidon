   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Allocation_Yahil_Profile                                                     !##!
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
USE Variables_Mesh, &
            ONLY :  Num_R_Elements

USE Variables_External, &
            ONLY :  SelfSim_R_Vals,     &
                    Selfsim_Pot_Vals,   &
                    Selfsim_Shift_Vals, &
                    SelfSim_Allocated
                    
USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag




IMPLICIT NONE


CONTAINS



!+101+###########################################################################!
!                                                                                !
!                            Allocate_Poseidon_Variables                         !
!                                                                                !
!################################################################################!
SUBROUTINE Allocate_Yahil_Profile(Num_Entries)

INTEGER, INTENT(IN)                             ::  Num_Entries


IF ( SelfSim_Allocated ) THEN
    IF ( Verbose_Flag ) THEN
        WRITE(*,'(A)')'Warning attempting to reallocate Selfsim Variables.'
    END IF
ELSE
    ALLOCATE( SELFSIM_R_VALS(0:NUM_ENTRIES) )
    ALLOCATE( SELFSIM_POT_VALS(0:NUM_ENTRIES) )
    ALLOCATE( SELFSIM_SHIFT_VALS(0:Num_R_Elements+1) )
    SelfSim_Allocated = .TRUE.
END IF

END SUBROUTINE Allocate_Yahil_Profile











!+102+###########################################################################!
!                                                                                !
!                           Deallocate_Poseidon_Variables                        !
!                                                                                !
!################################################################################!
SUBROUTINE Deallocate_Yahil_Profile()

IF ( SelfSim_Allocated ) THEN

    DEALLOCATE( SELFSIM_R_VALS )
    DEALLOCATE( SELFSIM_POT_VALS )
    DEALLOCATE( SELFSIM_SHIFT_VALS )
    SelfSim_Allocated = .FALSE.

END IF

END SUBROUTINE Deallocate_Yahil_Profile







END MODULE Allocation_Yahil_Profile


