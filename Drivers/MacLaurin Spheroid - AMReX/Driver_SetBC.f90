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

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag

USE Poseidon_Main_Module, &
            ONLY :  Poseidon_CFA_Set_Uniform_Boundary_Conditions

USE Poseidon_Units_Module, &
            ONLY :  C_Square,       &
                    Set_Units,      &
                    Shift_Units,    &
                    Centimeter,     &
                    Second,         &
                    Gram,            &
                    Grav_Constant_G

USE Variables_External, &
            ONLY :  MacLaurin_SemiMinor,    &
                    MacLaurin_SemiMajor,    &
                    MacLaurin_Ecc,          &
                    MacLaurin_SphereType,   &
                    MacLaurin_Rho

USE MacLaurin_Module, &
            ONLY :  MacLaurin_Potential_Sub_B

IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!   Driver_SetBC                                                                !
!                                                                               !
!###############################################################################!
SUBROUTINE Driver_SetBC( Domain_Edge )

REAL(idp), DIMENSION(2)                                 ::  Domain_Edge

REAL(idp)                                               ::  Psi_BC
REAL(idp)                                               ::  AlphaPsi_BC
REAL(idp)                                               ::  Shift_Vector_BC

CHARACTER(LEN=1), DIMENSION(1:5)                        ::  INNER_BC_TYPES
CHARACTER(LEN=1), DIMENSION(1:5)                        ::  OUTER_BC_TYPES
REAL(idp), DIMENSION(1:5)                               ::  INNER_BC_VALUES
REAL(idp), DIMENSION(1:5)                               ::  OUTER_BC_VALUES


REAL(idp)                                               ::  Potential


IF ( Verbose_Flag ) THEN
    WRITE(*,'(A)')"In Driver, Setting Boundary Conditions."
END IF


IF ( MacLaurin_SemiMajor == MacLaurin_SemiMinor ) THEN
    
    Potential = -((4.0_idp/3.0_idp)*pi*MacLaurin_SemiMajor**3*MacLaurin_Rho)*Grav_Constant_G/Domain_Edge(2)        &
              * gram / centimeter
ELSE
    CALL MacLaurin_Potential_Sub_B(Domain_Edge(2), 0.5_idp*pi, 0.0_idp,Potential)
END IF

Psi_BC      = 1.0_idp - 0.5_idp*Potential/C_Square
AlphaPsi_BC = 1.0_idp + 0.5_idp*Potential/C_Square
Shift_Vector_BC = 0.0_idp



INNER_BC_TYPES = (/"N", "N","N","N","N"/)
OUTER_BC_TYPES = (/"D", "D","D","D","D"/)


INNER_BC_VALUES = (/0.0_idp, 0.0_idp, 0.0_idp, 0.0_idp, 0.0_idp /)
OUTER_BC_VALUES = (/Psi_BC,  AlphaPsi_BC, Shift_Vector_BC, 0.0_idp, 0.0_idp /)



CALL Poseidon_CFA_Set_Uniform_Boundary_Conditions("I", INNER_BC_TYPES, INNER_BC_VALUES)
CALL Poseidon_CFA_Set_Uniform_Boundary_Conditions("O", OUTER_BC_TYPES, OUTER_BC_VALUES)



END SUBROUTINE Driver_SetBC








END MODULE Driver_SetBC_Module


