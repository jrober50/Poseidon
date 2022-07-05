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


CHARACTER(LEN=1)                                ::  INNER_BC_TYPE
CHARACTER(LEN=1)                                ::  OUTER_BC_TYPE
REAL(idp)                                       ::  INNER_BC_VALUE
REAL(idp)                                       ::  OUTER_BC_VALUE

IF ( Verbose_Flag ) CALL Driver_Init_Message('Calculating boundary conditions.')


INNER_BC_TYPE = "N"
OUTER_BC_TYPE = "D"

INNER_BC_VALUE = 0.0_idp
OUTER_BC_VALUE = Potential_Solution(R_Outer, 0.0_idp, 0.0_idp)


CALL Poseidon_Set_Uniform_Boundary_Conditions("I", INNER_BC_TYPE, INNER_BC_VALUE)
CALL Poseidon_Set_Uniform_Boundary_Conditions("O", OUTER_BC_TYPE, OUTER_BC_VALUE)




END SUBROUTINE Driver_SetBC




END MODULE Driver_SetBC_Module

