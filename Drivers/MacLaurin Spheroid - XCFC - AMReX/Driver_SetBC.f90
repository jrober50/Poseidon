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

USE Poseidon_Numbers_Module, &
            ONLY :  pi
            
USE Poseidon_Units_Module, &
            ONLY :  GravPot_Units

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag

USE Poseidon_Message_Routines_Module, &
            ONLY :  Driver_Init_Message

USE Poseidon_Interface_Boundary_Conditions, &
            ONLY :  Poseidon_Set_Uniform_Boundary_Conditions

USE Poseidon_Units_Module, &
            ONLY :  C_Square,           &
                    Centimeter,         &
                    Gram,               &
                    Grav_Constant_G

USE Variables_Mesh, &
            ONLY :  R_Outer

USE Variables_External, &
            ONLY :  MLS_SemiMinor,    &
                    MLS_SemiMajor,    &
                    MLS_Rho

USE External_MLS_Solution_Module, &
            ONLY :  MacLaurin_Potential_Sub

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

CHARACTER(LEN=1), DIMENSION(1:5)                        ::  INNER_BC_TYPES
CHARACTER(LEN=1), DIMENSION(1:5)                        ::  OUTER_BC_TYPES
REAL(idp), DIMENSION(1:5)                               ::  INNER_BC_VALUES
REAL(idp), DIMENSION(1:5)                               ::  OUTER_BC_VALUES


REAL(idp)                                               ::  Potential


IF ( Verbose_Flag ) CALL Driver_Init_Message('Calculating boundary conditions.')


IF ( MLS_SemiMajor == MLS_SemiMinor ) THEN
    
    Potential = -((4.0_idp/3.0_idp)*pi*MLS_SemiMajor**3*MLS_Rho)*Grav_Constant_G/R_Outer        &
              * gram / centimeter
ELSE
    CALL MacLaurin_Potential_Sub(R_Outer, 0.5_idp*pi, 0.0_idp, Potential)
END IF

!print*,"Centimeter",Centimeter
!print*,"R_Outer",R_Outer
!PRINT*,"Potential",Potential,Potential/GravPot_Units
!print*,"C_Square",C_Square

Psi_BC      = 1.0_idp - 0.5_idp*Potential/C_Square
AlphaPsi_BC = 1.0_idp + 0.5_idp*Potential/C_Square
Shift_Vector_BC = 0.0_idp

!PRINT*,"Psi_BC",Psi_BC

!STOP
INNER_BC_TYPES = (/"N", "N","N","N","N"/)
OUTER_BC_TYPES = (/"D", "D","D","D","D"/)


INNER_BC_VALUES = (/0.0_idp, 0.0_idp, 0.0_idp, 0.0_idp, 0.0_idp /)
OUTER_BC_VALUES = (/Psi_BC,  AlphaPsi_BC, Shift_Vector_BC, 0.0_idp, 0.0_idp /)

IF ( Verbose_Flag ) CALL Driver_Init_Message('Setting boundary conditions.')


CALL Poseidon_Set_Uniform_Boundary_Conditions("I", INNER_BC_TYPES, INNER_BC_VALUES)
CALL Poseidon_Set_Uniform_Boundary_Conditions("O", OUTER_BC_TYPES, OUTER_BC_VALUES)



END SUBROUTINE Driver_SetBC








END MODULE Driver_SetBC_Module


