   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE External_CCS_Solution_Module                                          !##!
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
            ONLY :  Grav_Constant_G,        &
                    Gram,                  &
                    Centimeter,             &
                    C_Square,               &
                    Density_Units

USE Variables_External, &
            ONLY :  CCS_SurfaceRadius,      &
                    CCS_CoreDensity,        &
                    CCS_CoreRadius
            

IMPLICIT NONE

PUBLIC  :: CCS_Density, CCS_Potential

CONTAINS



 !+201+####################################################!
!                                                           !
!        CCS_Density                                        !
!                                                           !
 !#########################################################!
FUNCTION CCS_Density(r)


REAL(idp), INTENT(IN)       ::  r

REAL(idp)                   ::  CCS_Density

REAL(idp)                   ::  r_wUnits
REAL(idp)                   ::  roverrc



r_wUnits = r/Centimeter

roverrc = r/CCS_CoreRadius

IF ( r .LE. CCS_SurfaceRadius ) THEN

    CCS_Density = CCS_CoreDensity/(1.0_idp + roverrc*roverrc)

ELSE

    CCS_Density = 0.0_idp
    
END IF


END FUNCTION CCS_Density





 !+201+####################################################!
!                                                           !
!        CCS_Potential                                      !
!                                                           !
 !#########################################################!
FUNCTION CCS_Potential(r)


REAL(idp), INTENT(IN)       ::  r

REAL(idp)                   ::  CCS_Potential

REAL(idp)                   ::  r_wUnits
REAL(idp)                   ::  fourpiGrho
REAL(idp)                   ::  roverrc
REAL(idp)                   ::  bigroverrc

r_wUnits = r/Centimeter

fourpiGrho = 4.0_idp*pi*Grav_Constant_G*CCS_CoreDensity

roverrc = r/CCS_CoreRadius
bigroverrc = CCS_SurfaceRadius/CCS_CoreRadius

IF ( r .LE. CCS_SurfaceRadius ) THEN

    CCS_Potential = -fourpiGrho*CCS_CoreRadius*CCS_CoreRadius                      &
                  * ( 1.0_idp - atan(roverrc)/roverrc                                    &
                    - 0.5_idp*log((1.0_idp+roverrc*roverrc)/(1.0_idp+bigroverrc*bigroverrc) ) )

ELSE

    CCS_Potential = -fourpiGrho*CCS_CoreRadius*CCS_CoreRadius*CCS_CoreRadius/r  &
                  * (bigroverrc - atan(bigroverrc))
                  
!    print*,bigroverrc - atan(bigroverrc),fourpiGrho,CCS_CoreRadius*CCS_CoreRadius*CCS_CoreRadius

END IF


END FUNCTION CCS_Potential





END MODULE External_CCS_Solution_Module
