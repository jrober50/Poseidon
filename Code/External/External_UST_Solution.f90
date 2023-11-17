   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE External_UST_Solution_Module                                             !##!
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

USE Poseidon_Units_Module, &
            ONLY :  Grav_Constant_G

USE Poseidon_Numbers_Module, &
            ONLY :  pi

USE Poseidon_Units_Module, &
            ONLY :  Gram,               &
                    Centimeter,         &
                    Second,             &
                    GravPot_Units

USE Variables_External, &
            ONLY :  UST_Rhoo,           &
                    UST_StarRadius

IMPLICIT NONE


CONTAINS



 !+101+####################################################!
!                                                           !
!          Set_HCT_Test_Params                                 !
!                                                           !
 !#########################################################!
SUBROUTINE Set_UST_Test_Params( Density_In,           &
                                Star_Radius_In      )

REAL(idp),  INTENT(IN)          ::  Density_In
REAL(idp),  INTENT(IN)          ::  Star_Radius_In

UST_Rhoo        = Density_In/(Gram/Centimeter**3)
UST_StarRadius = Star_Radius_In/Centimeter

!UST_Rhoo        = Density_In
!UST_StarRadius = Star_Radius_In

END SUBROUTINE Set_UST_Test_Params







 !+201+####################################################!
!                                                           !
!        UST_Solution                                       !
!                                                           !
 !#########################################################!
FUNCTION UST_Solution(r)


REAL(idp), INTENT(IN)       ::  r

REAL(idp)                   ::  r_wUnits
REAL(idp)                   ::  R_woUnits
REAL(idp)                   ::  UST_Solution

r_wUnits = r/Centimeter
R_woUnits = UST_StarRadius*Centimeter
!r_wUnits = r
IF ( r .LE. R_woUnits ) THEN
    UST_Solution = -(2.0_idp/3.0_idp)* pi                   &
                 * UST_Rhoo *  Grav_Constant_G              &
                 * (3.0_idp*R_woUnits*R_woUnits - r*r)
ELSE
    UST_Solution = -(4.0_idp/3.0_idp)* pi                   &
                 * UST_Rhoo *  Grav_Constant_G              &
                 * R_woUnits*R_woUnits*R_woUnits            &
                 / r
END IF

UST_Solution = UST_Solution!/(centimeter**3/(second*second))


END FUNCTION UST_Solution




END MODULE External_UST_Solution_Module

