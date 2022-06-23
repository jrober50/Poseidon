   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Poseidon_Return_Routines_Module                                       !##!
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
USE Return_Functions_FP, &
            ONLY :  Calc_Var_At_Location_Type_A,        &
                    Calc_Var_At_Location_Type_B,        &
                    Calc_Values_Here_Type_A,            &
                    Calc_Values_Here_Type_B
                    

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




INTERFACE Calc_Var_At_Location
    MODULE PROCEDURE Calc_Var_At_Location_Type_A
    MODULE PROCEDURE Calc_Var_At_Location_Type_B
END INTERFACE Calc_Var_At_Location


INTERFACE Calc_Values_Here_XCFC
    MODULE PROCEDURE Calc_Values_Here_Type_A
    MODULE PROCEDURE Calc_Values_Here_Type_B
END INTERFACE Calc_Values_Here_XCFC



CONTAINS





END MODULE Poseidon_Return_Routines_Module
