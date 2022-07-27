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
            ONLY :  Gram, Centimeter, Second

USE Variables_External, &
            ONLY :  UST_Rhoo,          &
                    UST_Star_Radius

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
UST_Star_Radius = Star_Radius_In/Centimeter

!UST_Rhoo        = Density_In
!UST_Star_Radius = Star_Radius_In

END SUBROUTINE Set_UST_Test_Params







 !+201+####################################################!
!                                                           !
!        UST_Solution                                       !
!                                                           !
 !#########################################################!
FUNCTION UST_Solution(r)


REAL(idp), INTENT(IN)       ::  r

REAL(idp)                   ::  r_wUnits
REAL(idp)                   ::  UST_Solution

r_wUnits = r/Centimeter
!r_wUnits = r

IF ( r_wUnits .LE. UST_Star_Radius ) THEN
    UST_Solution = -(2.0_idp/3.0_idp)* pi                                       &
                 * UST_Rhoo *  Grav_Constant_G                                  &
                 * (3.0_idp*UST_Star_Radius*UST_Star_Radius - r*r)
ELSE
    UST_Solution = -(4.0_idp/3.0_idp)* pi                               &
                 * UST_Rhoo *  Grav_Constant_G                          &
                 * UST_Star_Radius*UST_Star_Radius*UST_Star_Radius      &
                 / r
END IF

UST_Solution = UST_Solution!/(centimeter**3/(second*second))


END FUNCTION UST_Solution




END MODULE External_UST_Solution_Module

