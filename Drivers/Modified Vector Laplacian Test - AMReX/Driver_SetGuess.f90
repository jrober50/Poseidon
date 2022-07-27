   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Driver_SetGuess_Module                                                !##!
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

USE Poseidon_Units_Module, &
            ONLY :  C_Square

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag

USE Variables_Functions, &
            ONLY :  Potential_Solution

USE Poseidon_Interface_Initial_Guess, &
            ONLY :  Poseidon_Initialize_Flat_Guess


IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!     Driver_SetGuess                                                			!
!                                                                               !
!###############################################################################!
SUBROUTINE Driver_SetGuess()




IF ( Verbose_Flag ) THEN
    WRITE(*,'(A)')"In Driver, Setting Guess."
END IF

CALL Poseidon_Initialize_Flat_Guess()

END SUBROUTINE Driver_SetGuess



END MODULE Driver_SetGuess_Module
