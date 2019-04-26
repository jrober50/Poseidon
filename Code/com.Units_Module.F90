   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Units_Module                                                                 !##!
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
USE Poseidon_Constants_Module, &
            ONLY :  idp, pi,            &
                    Grav_Constant_G,    &
                    Speed_of_Light,     &
                    C_Square,           &
                    GR_Source_Scalar,   &
                    Meter,              &
                    Gram,               &
                    Kilogram,           &
                    Second






IMPLICIT NONE





                        !*F&S*==========================================!
                        !                                               !
                        !           Functions & Subroutines             !
                        !                                               !
                        !===============================================!
CONTAINS


SUBROUTINE Set_Units( Units_Flag )


CHARACTER(LEN = 1)                          :: Units_Flag


Second = 1.0_idp


IF ( Units_Flag == "C" ) THEN
!   CGS (centimeter-gram-second ) system
!
!   Mass    ->  g
!   Length  ->  cm
!   time    ->  s
!
!   c = 2.99792458e10 cm / s
!   G = 6.67428e-8    cm^3 / ( g s^2 )


    Meter = 100.0_idp     ! centimeters

    Gram = 1.0_idp        ! gram
    Kilogram = 1000_idp   ! grams


    !Potential_Units = "cm/s^2"

    Grav_Constant_G = 6.67408E-11 * (Meter*Meter*Meter)/(Kilogram * Second * Second)

    Speed_of_Light = 2.99792458E8 * (Meter / Second )
  

ELSE IF ( Units_Flag == "S" ) THEN
    !   SI Units / MKS
    !
    !   Mass    ->  kg
    !   Length  ->  m
    !   time    ->  s
    !
    !   c = 2.99792458e8    m / s
    !   G = 6.67408e-11     m^2 / ( kg s^2 )


    Meter = 1.0_idp       ! meter

    Gram = .0010_idp      ! kilograms
    Kilogram = 1.0_idp    ! kilogram


    !Potential_Units = "m/s^2"

    Grav_Constant_G = 6.67408E-11 * (Meter*Meter*Meter)/(Kilogram * Second * Second)

    Speed_of_Light = 2.99792458E8 * (Meter / Second )




ELSE IF ( Units_Flag == "G" ) THEN
    !   Geometrized Units -> Poseidon Units
    !
    !   Mass    ->  g
    !   Length  ->  m
    !   time    ->  s
    !
    !   c = 1
    !   G = 1

    Meter = 1.0_idp       ! meter

    Gram = 1.0_idp        ! gram
    Kilogram = 1000.0_idp ! grams


    Grav_Constant_G = 1.0_idp
    Speed_of_Light = 1.0_idp


    PRINT*,"Grav_Constant_G",Grav_Constant_G
    PRINT*,"Speed_of_Light",Speed_of_Light
    PRINT*," "
    PRINT*," "

END IF



C_Square = Speed_of_Light*Speed_of_Light
GR_Source_Scalar = Grav_Constant_G/(C_SQUARE*C_SQUARE)



END SUBROUTINE Set_Units




END MODULE Units_Module
