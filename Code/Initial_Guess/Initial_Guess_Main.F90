   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Poseidon_Initial_Guess_Module                                         !##!
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

USE IG_Input_Native_Module, &
            ONLY :  IG_Input_Native,               &
                    IG_Input_Native_Caller

USE IG_Flat_Guess_Module, &
            ONLY :  IG_Init_Flat_Guess

IMPLICIT NONE




INTERFACE Poseidon_Input_Initial_Guess
    MODULE PROCEDURE IG_Input_Native
    MODULE PROCEDURE IG_Input_Native_Caller
END INTERFACE Poseidon_Input_Initial_Guess


INTERFACE Poseidon_Initialize_Flat_Guess
    MODULE PROCEDURE IG_Init_Flat_Guess
END INTERFACE Poseidon_Initialize_Flat_Guess


CONTAINS



END MODULE Poseidon_Initial_Guess_Module
