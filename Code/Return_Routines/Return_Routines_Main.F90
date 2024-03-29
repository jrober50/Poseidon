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
                    


IMPLICIT NONE



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
