   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE External_HCT_Solution_Module                                    	     !##!
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
            ONLY : idp

USE Poseidon_Numbers_Module, &
            ONLY :  pi

USE Poseidon_Units_Module, &
            ONLY :  Centimeter

USE Variables_External, &
            ONLY :  HCT_Alpha,          &
                    HCT_Star_Radius,    &
                    HCT_Rhoo,           &
                    HCT_Beta,           &
                    HCT_C

IMPLICIT NONE


CONTAINS



 !+101+####################################################!
!                                                           !
!          Set_HCT_Test_Params                         	    !
!                                                           !
 !#########################################################!
SUBROUTINE Set_HCT_Test_Params( Alpha_In,           &
                                Star_Radius_In      )

REAL(idp),  INTENT(IN)          ::  Alpha_In
REAL(idp),  INTENT(IN)          ::  Star_Radius_In

HCT_Alpha       = Alpha_In
HCT_Star_Radius = Star_Radius_In*Centimeter

HCT_Rhoo = (3.0_idp*FofAlpha(HCT_Alpha)**2) / (2.0_idp*pi*HCT_Star_Radius*HCT_Star_Radius)
HCT_C    = (2.0_idp/3.0_idp*pi*HCT_Rhoo)**(-1.0/4.0)
HCT_Beta = (HCT_C*UsubAlpha(HCT_Star_Radius) - 1.0_idp) * HCT_Star_Radius

END SUBROUTINE Set_HCT_Test_Params







 !+201+####################################################!
!                                                           !
!          UsubAlpha                                        !
!                                                           !
 !#########################################################!
FUNCTION HCT_Solution(r)


REAL(idp), INTENT(IN)       ::  r

REAL(idp)                   ::  r_wUnits
REAL(idp)                   ::  HCT_Solution

r_wUnits = r/Centimeter

IF ( r .LE. HCT_Star_Radius ) THEN
    HCT_Solution = HCT_C * UsubAlpha(r)
ELSE
    HCT_Solution = HCT_Beta/r + 1.0_idp
END IF


END FUNCTION HCT_Solution





 !+301+####################################################!
!                                                           !
!          UsubAlpha                                        !
!                                                           !
 !#########################################################!
FUNCTION UsubAlpha(r)


REAL(idp), INTENT(IN)       ::  r

REAL(idp)                   ::  UsubAlpha


UsubAlpha = sqrt(HCT_Alpha*HCT_Star_Radius              &
                /(r*r + HCT_Alpha*HCT_Alpha*HCT_Star_Radius*HCT_Star_Radius) )

END FUNCTION UsubAlpha



 !+302+####################################################!
!                                                           !
!          UsubAlpha                                        !
!                                                           !
 !#########################################################!
FUNCTION FofAlpha(Alpha)


REAL(idp), INTENT(IN)       ::  Alpha

REAL(idp)                   ::  FofAlpha

FofAlpha = Alpha**5/(1 + Alpha**2)**3


END FUNCTION FofAlpha


END MODULE External_HCT_Solution_Module
