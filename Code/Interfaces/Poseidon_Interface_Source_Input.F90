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

USE Source_Input_XCFC_AMReX_Module, &
            ONLY :  Poseidon_Input_Sources_XCFC_AMReX,              &
                    Poseidon_Input_Sources_XCFC_AMReX_Caller,       &
                    Poseidon_Input_Sources_Part1_AMReX,             &
                    Poseidon_Input_Sources_Part1_AMReX_Caller

USE Source_Input_XCFC_Native_Module, &
            ONLY :  Poseidon_Input_Sources_XCFC_Native,             &
                    Poseidon_Input_Sources_XCFC_Native_Caller,      &
                    Poseidon_Input_Sources_Part1_Native,            &
                    Poseidon_Input_Sources_Part1_Native_Caller

USE Source_Input_Newtonian_Native_Module, &
            ONLY :  Poseidon_Input_Sources_Newtonian_Native,        &
                    Poseidon_Input_Sources_Newtonian_Native_Caller

USE Source_Input_Newtonian_AMReX_Module, &
            ONLY :  Poseidon_Input_Sources_Newtonian_AMReX,        &
                    Poseidon_Input_Sources_Newtonian_AMReX_Caller


IMPLICIT NONE




INTERFACE Poseidon_Input_Sources
    MODULE PROCEDURE Poseidon_Input_Sources_XCFC_Native
    MODULE PROCEDURE Poseidon_Input_Sources_XCFC_Native_Caller
    MODULE PROCEDURE Poseidon_Input_Sources_XCFC_AMReX
    MODULE PROCEDURE Poseidon_Input_Sources_XCFC_AMReX_Caller
    MODULE PROCEDURE Poseidon_Input_Sources_Newtonian_Native
    MODULE PROCEDURE Poseidon_Input_Sources_Newtonian_Native_Caller
!    MODULE PROCEDURE Poseidon_Input_Sources_Newtonian_AMReX
!    MODULE PROCEDURE Poseidon_Input_Sources_Newtonian_AMReX_Caller
END INTERFACE Poseidon_Input_Sources


INTERFACE Poseidon_Input_Sources_Part1
    MODULE PROCEDURE Poseidon_Input_Sources_Part1_Native
    MODULE PROCEDURE Poseidon_Input_Sources_Part1_Native_Caller
    MODULE PROCEDURE Poseidon_Input_Sources_Part1_AMReX
    MODULE PROCEDURE Poseidon_Input_Sources_Part1_AMReX_Caller
END INTERFACE Poseidon_Input_Sources_Part1


INTERFACE Poseidon_Input_Sources_Part2
    MODULE PROCEDURE Poseidon_Input_Sources_XCFC_Native
    MODULE PROCEDURE Poseidon_Input_Sources_XCFC_Native_Caller
    MODULE PROCEDURE Poseidon_Input_Sources_XCFC_AMReX
    MODULE PROCEDURE Poseidon_Input_Sources_XCFC_AMReX_Caller
END INTERFACE Poseidon_Input_Sources_Part2




CONTAINS

END MODULE Poseidon_Interface_Source_Input
