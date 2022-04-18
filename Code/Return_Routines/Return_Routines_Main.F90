   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Return_Routines_Main                                                  !##!
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



!+101+##########################################################################!
!                                                                               !
!                                                     				!
!                                                                               !
!###############################################################################!




END MODULE Return_Routines_Main
