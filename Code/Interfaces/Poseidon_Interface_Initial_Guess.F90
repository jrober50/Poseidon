   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Poseidon_Interface_Initial_Guess                                      !##!
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
            ONLY :  IG_Input_XCFC_Native,           &
                    IG_Input_XCFC_Native_Caller,    &
                    IG_Input_Poisson_Native,        &
                    IG_Input_Poisson_Native_Caller

USE IG_Flat_Guess_Module, &
            ONLY :  IG_Init_Flat_Guess

IMPLICIT NONE




INTERFACE Poseidon_Input_Initial_Guess
    MODULE PROCEDURE IG_Input_XCFC_Native
    MODULE PROCEDURE IG_Input_XCFC_Native_Caller
    MODULE PROCEDURE IG_Input_Poisson_Native
    MODULE PROCEDURE IG_Input_Poisson_Native_Caller
END INTERFACE Poseidon_Input_Initial_Guess


INTERFACE Poseidon_Initialize_Flat_Guess
    MODULE PROCEDURE IG_Init_Flat_Guess
END INTERFACE Poseidon_Initialize_Flat_Guess


CONTAINS



END MODULE Poseidon_Interface_Initial_Guess
