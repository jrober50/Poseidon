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
            ONLY :  C_Square,       &
                    Centimeter

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
SUBROUTINE Driver_SetBC( NE, x_e )

INTEGER,   DIMENSION(3),       INTENT(IN)       ::  NE
REAL(idp), DIMENSION(0:NE(1)), INTENT(IN)       ::  x_E

REAL(idp)                                       ::  Psi_BC
REAL(idp)                                       ::  AlphaPsi_BC
REAL(idp)                                       ::  Shift_Vector_BC
CHARACTER(LEN=1), DIMENSION(1:5)                ::  INNER_BC_TYPES, OUTER_BC_TYPES
REAL(idp), DIMENSION(1:5)                       ::  INNER_BC_VALUES, OUTER_BC_VALUES

IF ( Verbose_Flag ) CALL Driver_Init_Message('Calculating boundary conditions.')




Psi_BC = 1.0_idp    &
       - 0.5_idp*Potential_Solution(R_Outer*Centimeter, 0.0_idp, 0.0_idp)/C_Square

AlphaPsi_BC = 1.0_idp    &
            + 0.5_idp*Potential_Solution(R_Outer*Centimeter, 0.0_idp, 0.0_idp)/C_Square



INNER_BC_TYPES = (/"N", "N","N","N","N"/)
OUTER_BC_TYPES = (/"D", "D","D","D","D"/)


INNER_BC_VALUES = (/0.0_idp, 0.0_idp, 0.0_idp, 0.0_idp, 0.0_idp /)
OUTER_BC_VALUES = (/Psi_BC,  AlphaPsi_BC, 0.0_idp, 0.0_idp, 0.0_idp /)
!    OUTER_BC_VALUES = (/Psi_BC,  AlphaPsi_BC, Shift_Vector_BC, 0.0_idp, 0.0_idp /)



IF ( Verbose_Flag ) CALL Driver_Init_Message('Setting boundary conditions.')

CALL Poseidon_Set_Uniform_Boundary_Conditions("I", INNER_BC_TYPES, INNER_BC_VALUES)
CALL Poseidon_Set_Uniform_Boundary_Conditions("O", OUTER_BC_TYPES, OUTER_BC_VALUES)



END SUBROUTINE Driver_SetBC




END MODULE Driver_SetBC_Module

