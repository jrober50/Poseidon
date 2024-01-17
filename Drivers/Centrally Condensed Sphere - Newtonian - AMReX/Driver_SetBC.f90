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
            ONLY :  Grav_Constant_G,    &
                    Speed_of_Light,     &
                    C_Square,           &
                    GR_Source_Scalar,   &
                    Centimeter,         &
                    Second,             &
                    Millisecond,         &
                    Erg,                &
                    Gram

USE Variables_Mesh, &
            ONLY :  R_Outer

USE Poseidon_Interface_Boundary_Conditions, &
            ONLY :  Poseidon_Set_Uniform_Boundary_Conditions

USE External_CCS_Solution_Module, &
            ONLY :  CCS_Potential


IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!   Driver_SetBC                                                                !
!                                                                               !
!###############################################################################!
SUBROUTINE Driver_SetBC( )

REAL(idp)                                               ::  Psi_BC
REAL(idp)                                               ::  AlphaPsi_BC
REAL(idp)                                               ::  Shift_Vector_BC

CHARACTER(LEN=1)                        ::  INNER_BC_TYPES
CHARACTER(LEN=1)                        ::  OUTER_BC_TYPES
REAL(idp)                               ::  INNER_BC_VALUES
REAL(idp)                              ::  OUTER_BC_VALUES

REAL(idp)                                               ::  Potential


IF ( Verbose_Flag ) CALL Driver_Init_Message('Calculating boundary conditions.')



Potential = CCS_Potential(R_Outer)



Psi_BC          = Potential




INNER_BC_TYPES = "N"
OUTER_BC_TYPES = "D"


INNER_BC_VALUES = 0.0_idp
OUTER_BC_VALUES = Psi_BC 

IF ( Verbose_Flag ) CALL Driver_Init_Message('Setting boundary conditions.')



CALL Poseidon_Set_Uniform_Boundary_Conditions("I", INNER_BC_TYPES, INNER_BC_VALUES)
CALL Poseidon_Set_Uniform_Boundary_Conditions("O", OUTER_BC_TYPES, OUTER_BC_VALUES)



END SUBROUTINE Driver_SetBC



END MODULE Driver_SetBC_Module

