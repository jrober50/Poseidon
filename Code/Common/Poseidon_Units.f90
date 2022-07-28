   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_Units_Module                                                        !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains the functions, and subroutines used to define and convert to and   !##!
!##!    from different unit systems to that which Poseidon uses.                    !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!

                            !*D*================================!
                            !                                   !
                            !           Dependencies            !
                            !                                   !
                            !===================================!
USE Poseidon_Kinds_Module, &
            ONLY :  idp






IMPLICIT NONE

REAL(idp)               ::  Speed_of_Light_MKS  = 2.99792458e8_idp
REAL(idp)               ::  Grav_Constant_MKS   = 6.673e-11_idp

REAL(idp)               ::  Grav_Constant_G
REAL(idp)               ::  Speed_of_Light
REAL(idp)               ::  C_Square
REAL(idp)               ::  GR_Source_Scalar
REAL(idp)               ::  Meter, Kilometer, Centimeter
REAL(idp)               ::  Gram, Kilogram, SolarMass
REAL(idp)               ::  Second, Millisecond
REAL(idp)               ::  Joule, Erg, Newton
REAL(idp)               ::  GravPot_Units, Shift_Units
REAL(idp)               ::  E_Units, S_Units, Si_Units

REAL(idp)               ::  Solar_Mass



CONTAINS



 !+101+####################################################!
!                                                           !
!          Set_Units                                        !
!                                                           !
 !#########################################################!
SUBROUTINE Set_Units( Units_Flag )


CHARACTER(LEN = 1)                          :: Units_Flag

IF ( Units_Flag == "C" ) THEN
!   CGS (centimeter-gram-second ) system
!
!   Mass    ->  g
!   Length  ->  cm
!   time    ->  s
!
!   c = 2.99792458e10 cm / s
!   G = 6.67428e-8    cm^3 / ( g s^2 )

!    IF ( Verbose_Flag ) THEN
!        PRINT*,"Setting CGS Units"
!    END IF
    ! Length
    Centimeter  = 1.0_idp
    Meter       = 1.0E+2_idp * Centimeter
    Kilometer   = 1.0E+5_idp * Centimeter

    ! Mass
    Gram        = 1.0_idp
    Kilogram    = 1.0E+3_idp * Gram
    SolarMass   = 1.98892E30_idp * Kilogram

    ! Time
    Second      = 1.0_idp
    Millisecond  = 1.0E-3_idp * Second


    Grav_Constant_G = Grav_Constant_MKS * (Meter*Meter*Meter)/(Kilogram * Second * Second)
    Speed_of_Light  = Speed_of_Light_MKS * (Meter / Second )
    
    


ELSE IF ( Units_Flag == "S" ) THEN
    !   SI Units / MKS
    !
    !   Mass    ->  kg
    !   Length  ->  m
    !   time    ->  s
    !
    !   c = 2.99792458e8    m / s
    !   G = 6.67408e-11     m^2 / ( kg s^2 )

    ! Length
    Meter       = 1.0_idp
    Centimeter  = 1.0E-2_idp * Meter
    Kilometer   = 1.0E+3_idp * Meter

    ! Mass
    Kilogram    = 1.0_idp
    Gram        = 1.0E-3_idp * Kilogram
    SolarMass   = 1.98892E30_idp * Kilogram

    ! Time
    Second      = 1.0_idp
    Millisecond  = 1.0E-3_idp * Second

    Grav_Constant_G = Grav_Constant_MKS * (Meter*Meter*Meter)/(Kilogram * Second * Second)
    Speed_of_Light  = Speed_of_Light_MKS * (Meter / Second )




ELSE IF ( Units_Flag == "G" ) THEN
    !   Geometrized Units -> Poseidon Units
    !
    !   Mass    ->  g
    !   Length  ->  m
    !   time    ->  s
    !
    !   c = 1
    !   G = 1

!    IF ( Verbose_Flag ) THEN
!        PRINT*,"Setting Geometrized Units"
!    END IF

    Grav_Constant_G = 1.0_idp
    Speed_of_Light = 1.0_idp

    ! Length
    Meter       = 1.0_idp
    Centimeter  = 1.0E-2_idp * Meter
    Kilometer   = 1.0E+3_idp * Meter

    ! Time
    Second      = (Speed_of_Light_MKS/Speed_of_Light) * Meter
    Millisecond  = 1.0E-3_idp * Second

    ! Mass
    Kilogram    = (Grav_Constant_MKS/Grav_Constant_G) * Meter**3/Second**2
    Gram        = 1.0E-3_idp * Kilogram
    SolarMass   = 1.98892E30_idp * Kilogram

ELSE IF ( Units_Flag == "U" ) THEN
    ! Unitless
!    IF ( Verbose_Flag ) THEN
!        PRINT*,"Setting Unitless Units"
!    END IF

    Grav_Constant_G = 1.0_idp
    Speed_of_Light = 1.0_idp

    ! Length
    Meter       = 1.0_idp
    Centimeter  = 1.0_idp
    Kilometer   = 1.0_idp

    ! Time
    Second      = 1.0_idp
    Millisecond  = 1.0_idp

    ! Mass
    Kilogram    = 1.0_idp
    Gram        = 1.0_idp
    SolarMass   = 1.0_idp



END IF



C_Square = Speed_of_Light*Speed_of_Light
GR_Source_Scalar = Grav_Constant_G/(C_SQUARE*C_SQUARE)

Joule           = Kilogram * (Meter/Second)**2
Erg             = Gram * ( Centimeter / Second )**2
Newton          = Joule / Meter
GravPot_Units   = (Meter*Meter)/(Second*Second)
Shift_Units     = Centimeter/Second

E_Units         = Erg/Centimeter**3
S_Units         = Erg/Centimeter**3
Si_Units        = Gram/(Second*Centimeter**2)

Solar_Mass      = 1.988435E+30 * Kilogram

END SUBROUTINE Set_Units




END MODULE Poseidon_Units_Module
