   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Poseidon_Interface_Source_Input                                       !##!
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

USE Source_Input_AMReX_Module, &
            ONLY :  Poseidon_Input_Sources_AMReX,               &
                    Poseidon_Input_Sources_AMReX_Caller,        &
                    Poseidon_Input_Sources_Part1_AMReX,         &
                    Poseidon_Input_Sources_Part1_AMReX_Caller

USE Source_Input_Native_Module, &
            ONLY :  Poseidon_Input_Sources_Native,               &
                    Poseidon_Input_Sources_Native_Caller,        &
                    Poseidon_Input_Sources_Part1_Native,         &
                    Poseidon_Input_Sources_Part1_Native_Caller,  &
                    Poseidon_Input_Sources_Part2_Native,         &
                    Poseidon_Input_Sources_Part2_Native_Caller


IMPLICIT NONE




INTERFACE Poseidon_Input_Sources
    MODULE PROCEDURE Poseidon_Input_Sources_Native
    MODULE PROCEDURE Poseidon_Input_Sources_Native_Caller
    MODULE PROCEDURE Poseidon_Input_Sources_AMReX
    MODULE PROCEDURE Poseidon_Input_Sources_AMReX_Caller
END INTERFACE Poseidon_Input_Sources


INTERFACE Poseidon_Input_Sources_Part1
    MODULE PROCEDURE Poseidon_Input_Sources_Part1_Native
    MODULE PROCEDURE Poseidon_Input_Sources_Part1_Native_Caller
    MODULE PROCEDURE Poseidon_Input_Sources_Part1_AMReX
    MODULE PROCEDURE Poseidon_Input_Sources_Part1_AMReX_Caller
END INTERFACE Poseidon_Input_Sources_Part1


INTERFACE Poseidon_Input_Sources_Part2
    MODULE PROCEDURE Poseidon_Input_Sources_Part2_Native
    MODULE PROCEDURE Poseidon_Input_Sources_Part2_Native_Caller
    MODULE PROCEDURE Poseidon_Input_Sources_AMReX
    MODULE PROCEDURE Poseidon_Input_Sources_AMReX_Caller
END INTERFACE Poseidon_Input_Sources_Part2




CONTAINS

END MODULE Poseidon_Interface_Source_Input
