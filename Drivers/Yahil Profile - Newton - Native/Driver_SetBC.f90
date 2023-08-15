   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Driver_SetBC_Module                                                   !##!
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

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag

USE Poseidon_Message_Routines_Module, &
            ONLY :  Driver_Init_Message

USE Poseidon_Units_Module, &
            ONLY :  C_Square

USE Variables_Functions, &
            ONLY :  Potential_Solution

USE Poseidon_Interface_Boundary_Conditions, &
            ONLY :  Poseidon_Set_Uniform_Boundary_Conditions

USE Variables_Mesh, &
            ONLY :  R_Outer


IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!   Driver_SetBC                                                                !
!                                                                               !
!###############################################################################!
SUBROUTINE Driver_SetBC( )


CHARACTER(LEN=1)                                ::  Inner_BC_Type
CHARACTER(LEN=1)                                ::  Outer_BC_Type
REAL(idp)                                       ::  Inner_BC_Value
REAL(idp)                                       ::  Outer_BC_Value

IF ( Verbose_Flag ) CALL Driver_Init_Message('Calculating boundary conditions.')



Inner_BC_Type = "N"
Outer_BC_Type = "D"

Inner_BC_Value = 0.0_idp
Outer_BC_Value = Potential_Solution(R_Outer, 0.0_idp, 0.0_idp)



CALL Poseidon_Set_Uniform_Boundary_Conditions("I", Inner_BC_Type, Inner_BC_Value)
CALL Poseidon_Set_Uniform_Boundary_Conditions("O", Outer_BC_Type, Outer_BC_Value)




END SUBROUTINE Driver_SetBC




END MODULE Driver_SetBC_Module

