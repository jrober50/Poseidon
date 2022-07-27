   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Poseidon_Interface_Return_Routines                                    !##!
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


USE Parameters_Variable_Indices, &
            ONLY :  iU_CF,                  &
                    iU_LF,                  &
                    iU_S1,                  &
                    iU_S2,                  &
                    iU_S3,                  &
                    iVB_S


USE Poseidon_Return_Routines_CF, &
            ONLY :  Poseidon_Return_CF_Native,          &
                    Poseidon_Return_CF_Native_Caller,   &
                    Poseidon_Return_CF_AMReX,           &
                    Poseidon_Return_CF_AMReX_Caller

USE Poseidon_Return_Routines_LF, &
            ONLY :  Poseidon_Return_LF_Native,          &
                    Poseidon_Return_LF_Native_Caller,   &
                    Poseidon_Return_LF_AMReX,           &
                    Poseidon_Return_LF_AMReX_Caller

USE Poseidon_Return_Routines_SV, &
            ONLY :  Poseidon_Return_SV_Native,          &
                    Poseidon_Return_SV_Native_Caller,   &
                    Poseidon_Return_SV_AMReX,           &
                    Poseidon_Return_SV_AMReX_Caller

USE Poseidon_Return_Routines_XV, &
            ONLY :  Poseidon_Return_XV_Native,          &
                    Poseidon_Return_XV_Native_Caller,   &
                    Poseidon_Return_XV_AMReX,           &
                    Poseidon_Return_XV_AMReX_Caller

USE Poseidon_Return_Routines_Kij, &
            ONLY :  Poseidon_Return_Kij_Native,          &
                    Poseidon_Return_Kij_Native_Caller,   &
                    Poseidon_Return_Kij_AMReX,           &
                    Poseidon_Return_Kij_AMReX_Caller


USE Poseidon_Return_Routines_All, &
            ONLY :  Poseidon_Return_All_Native,          &
                    Poseidon_Return_All_Native_Caller,   &
                    Poseidon_Return_All_AMReX,           &
                    Poseidon_Return_All_AMReX_Caller


USE Poseidon_Return_Routines_Potential, &
            ONLY :  Poseidon_Return_Potential_Native,          &
                    Poseidon_Return_Potential_Native_Caller,   &
                    Poseidon_Return_Potential_AMReX,           &
                    Poseidon_Return_Potential_AMReX_Caller


IMPLICIT NONE


INTERFACE Poseidon_Return_Conformal_Factor
    MODULE PROCEDURE Poseidon_Return_CF_Native
    MODULE PROCEDURE Poseidon_Return_CF_Native_Caller
    MODULE PROCEDURE Poseidon_Return_CF_AMReX
    MODULE PROCEDURE Poseidon_Return_CF_AMReX_Caller
END INTERFACE Poseidon_Return_Conformal_Factor

INTERFACE Poseidon_Return_Lapse_Function
    MODULE PROCEDURE Poseidon_Return_LF_Native
    MODULE PROCEDURE Poseidon_Return_LF_Native_Caller
    MODULE PROCEDURE Poseidon_Return_LF_AMReX
    MODULE PROCEDURE Poseidon_Return_LF_AMReX_Caller
END INTERFACE Poseidon_Return_Lapse_Function

INTERFACE Poseidon_Return_Shift_Vector
    MODULE PROCEDURE Poseidon_Return_SV_Native
    MODULE PROCEDURE Poseidon_Return_SV_Native_Caller
    MODULE PROCEDURE Poseidon_Return_SV_AMReX
    MODULE PROCEDURE Poseidon_Return_SV_AMReX_Caller
END INTERFACE Poseidon_Return_Shift_Vector

INTERFACE Poseidon_Return_X_Vector
    MODULE PROCEDURE Poseidon_Return_XV_Native
    MODULE PROCEDURE Poseidon_Return_XV_Native_Caller
    MODULE PROCEDURE Poseidon_Return_XV_AMReX
    MODULE PROCEDURE Poseidon_Return_XV_AMReX_Caller
END INTERFACE Poseidon_Return_X_Vector

INTERFACE Poseidon_Return_Extrinsic_Curvature
    MODULE PROCEDURE Poseidon_Return_Kij_Native
    MODULE PROCEDURE Poseidon_Return_Kij_Native_Caller
    MODULE PROCEDURE Poseidon_Return_Kij_AMReX
    MODULE PROCEDURE Poseidon_Return_Kij_AMReX_Caller
END INTERFACE Poseidon_Return_Extrinsic_Curvature

INTERFACE Poseidon_Return_All
    MODULE PROCEDURE Poseidon_Return_All_Native
    MODULE PROCEDURE Poseidon_Return_All_Native_Caller
    MODULE PROCEDURE Poseidon_Return_All_AMReX
    MODULE PROCEDURE Poseidon_Return_All_AMReX_Caller
END INTERFACE Poseidon_Return_All





INTERFACE Poseidon_Return_Newtonian_Potential
    MODULE PROCEDURE Poseidon_Return_Potential_Native
    MODULE PROCEDURE Poseidon_Return_Potential_Native_Caller
    MODULE PROCEDURE Poseidon_Return_Potential_AMReX
    MODULE PROCEDURE Poseidon_Return_Potential_AMReX_Caller
END INTERFACE Poseidon_Return_Newtonian_Potential


CONTAINS





END MODULE Poseidon_Interface_Return_Routines 
