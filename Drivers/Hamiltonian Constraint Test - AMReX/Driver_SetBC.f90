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

USE Poseidon_Numbers_Module, &
            ONLY :  pi

USE Poseidon_Units_Module, &
            ONLY :  C_Square

USE Variables_Mesh, &
            ONLY :  R_Outer

USE Variables_Functions, &
            ONLY :  Potential_Solution

USE Variables_External, &
            ONLY :  HCT_Alpha,      &
                    HCT_Star_Radius

USE Poseidon_Interface_BC_Input, &
            ONLY :  Poseidon_Set_Uniform_Boundary_Conditions


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

REAL(idp)                                               ::  Beta
REAL(idp)                                               ::  C
REAL(idp)                                               ::  rho_o
REAL(idp)                                               ::  uaR
REAL(idp)                                               ::  fofalpha

IF ( Verbose_Flag ) CALL Driver_Init_Message('Calculating boundary conditions.')


fofalpha            =  HCT_Alpha**5/(1.0_idp+HCT_Alpha*HCT_Alpha)**3

rho_o               =  (3.0_idp/(2.0_idp*pi*HCT_Star_Radius*HCT_Star_Radius) )*fofalpha*fofalpha
uaR                 =  sqrt(HCT_Alpha/((1.0_idp+HCT_Alpha*HCT_Alpha)*HCT_Star_Radius))
C                   =  1.0_idp/sqrt(sqrt( (2.0_idp/3.0_idp)*pi*rho_o  ) )
Beta                =  (C*uaR-1.0_idp)*HCT_Star_Radius

Psi_BC          = 1.0250_idp
!Psi_BC          = HCT_Solution( R_Outer, HCT_Alpha, Beta, C, HCT_Star_Radius )
AlphaPsi_BC     = 1.0_idp
Shift_Vector_BC = 0.0_idp



INNER_BC_TYPES = (/"N", "N","N","N","N"/)
OUTER_BC_TYPES = (/"D", "D","D","D","D"/)


INNER_BC_VALUES = (/0.0_idp, 0.0_idp, 0.0_idp, 0.0_idp, 0.0_idp /)
OUTER_BC_VALUES = (/Psi_BC,  AlphaPsi_BC, Shift_Vector_BC, 0.0_idp, 0.0_idp /)


IF ( Verbose_Flag ) CALL Driver_Init_Message('Setting boundary conditions.')


CALL Poseidon_Set_Uniform_Boundary_Conditions("I", INNER_BC_TYPES, INNER_BC_VALUES)
CALL Poseidon_Set_Uniform_Boundary_Conditions("O", OUTER_BC_TYPES, OUTER_BC_VALUES)





END SUBROUTINE Driver_SetBC







!############################################################!
!#                                                          #!
!#                   HCT_Solution Function                  #!
!#                                                          #!
!############################################################!
REAL FUNCTION HCT_Solution( r, Alpha, Beta, C, Star_Radius )

REAL(idp), INTENT(IN)                      ::   r
REAL(idp), INTENT(IN)                      ::   Alpha
REAL(idp), INTENT(IN)                      ::   Beta
REAL(idp), INTENT(IN)                      ::   C
REAL(idp), INTENT(IN)                      ::   Star_Radius


IF ( r .LE. Star_Radius ) THEN

    HCT_Solution = C*SQRT( (Alpha*Star_Radius)/(r*r + (Alpha*Star_Radius)**2 ) )

ELSE

    HCT_Solution = Beta/r + 1.0_idp

END IF



END FUNCTION HCT_Solution

END MODULE Driver_SetBC_Module

