   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE External_MLS_Solution_Module                                                 !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
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

USE Poseidon_Units_Module, &
            ONLY :  Grav_Constant_G,        &
                    Gram,                  &
                    Centimeter,             &
                    C_Square

USE Poseidon_Numbers_Module, &
            ONLY :  pi

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag

USE Variables_MPI, &
            ONLY :  myID_Poseidon,          &
                    MasterID_Poseidon

USE Variables_Functions, &
            ONLY :  Potential_Solution

USE Variables_External, &
            ONLY :  MLS_SemiMinor,          &
                    MLS_SemiMajor,          &
                    MLS_Ecc,                &
                    MLS_SphereType,         &
                    MLS_Rho,                &
                    iMLS_Oblate,            &
                    iMLS_Prolate

USE Maps_X_Space, &
            ONLY :  Map_To_X_Space,         &
                    Map_From_X_Space

USE Maps_Quadrature,   &
            ONLY :  Quad_Map

USE Flags_IO_Module, &
            ONLY :  lPF_IO_Flags,           &
                    iPF_IO_Print_Setup

USE External_MLS_Profile_Module, &
            ONLY :  MacLaurin_Root_Finder,  &
                    Calc_MLS_ABCs


IMPLICIT NONE





CONTAINS


!+201+##########################################################
!
!   MacLaurin_Potential - Calulates the potential at a given location for the MacLaurin Spheriod
!
!################################################################
FUNCTION MacLaurin_Potential(r, theta, phi)


REAL(KIND = idp), INTENT(IN)                                ::  r, theta, phi


REAL(KIND = idp)                                            ::  MacLaurin_Potential

REAL(KIND = idp)                                            ::  X, Y, Z
REAL(KIND = idp)                                            ::  XX, YY, ZZ

REAL(idp)                                                   ::  A, B, C
REAL(idp)                                                   ::  AA, BB, CC

REAL(KIND = idp)                                            ::  VAL, &
                                                                eccsqr, eccmod, &
                                                                Parta, Partb, Partc, &
                                                                lambda, POTENTIAL


CALL Calc_MLS_ABCs( A, B, C, AA, BB, CC )


eccsqr = MLS_Ecc*MLS_Ecc



!   Calculate Cartesian Coordinates
X = r * cos(phi) * sin(theta)
Y = r * sin(phi) * sin(theta)
Z = r * cos(theta)

XX = X*X
YY = Y*Y
ZZ = Z*Z


VAL = XX/AA + YY/BB + ZZ/CC
!
!   VAL > 1 -> point is outside.  VAL <= 1 -> point is inside
!
IF (VAL .LE. 1.0_idp) THEN

    !
    !   Equations for potential inside homegenous spheroid
    !
    !       REF:    Ellipsoidal Figures of Equilibrium - Chandrasekhar - Chapter 3
    !
    IF ( MLS_SphereType == iMLS_Oblate ) THEN

        eccmod = sqrt(1- eccsqr)/(eccsqr*MLS_Ecc)*ASIN(MLS_Ecc)

        Parta = eccmod - (1 - eccsqr)/eccsqr

        Partb = 2/eccsqr - 2*eccmod


        POTENTIAL = pi*Grav_Constant_G*MLS_Rho * ((2*AA - XX - YY)*Parta + (CC - ZZ)*Partb)


    ELSE IF (MLS_SphereType == iMLS_Prolate) THEN

        eccmod = (1 - eccsqr)/(eccsqr*MLS_Ecc) * LOG((1+MLS_Ecc)/(1-MLS_Ecc))

        Parta = eccmod - 2*(1 - eccsqr)/eccsqr

        Partb = 1/eccsqr - 0.5_idp*eccmod


        POTENTIAL = pi*(2*CC - YY - ZZ)*Partb + (AA - XX)*Parta


    END IF


ELSE IF (VAL > 1.0_idp) THEN


    !
    !   Equations for potential outside homegenous spheroid
    !
    !       REF:    Ellipsoidal Figures of Equilibrium - Chandrasekhar - Chapter 3 -
    !
    lambda = MacLaurin_Root_Finder(r, theta, phi, AA, CC)

    Parta = ATAN(sqrt((AA - CC)/(CC + lambda)))

    Partb = 1/sqrt(AA - CC)

    Partc = Partb*Partb*Partb

    IF ( MLS_SphereType == iMLS_Oblate ) THEN




        POTENTIAL = -pi * MLS_Rho * Grav_Constant_G * AA * C                                &
                  * ( 2 * Partb * Parta                                                 &
                     - (XX + YY) * Partc                                                &
                        * (Parta - ( SQRT(AA-CC)*SQRT(CC + lambda) )/(AA + lambda) )    &
                     + 2*(ZZ) * Partc * ( Parta - sqrt(AA - CC)/sqrt(CC + lambda) ) )


    ELSE IF ( MLS_SphereType == iMLS_Prolate ) THEN


        POTENTIAL = A * CC                                                              &
                  * ( - Parta                                                           &
                      - Partb                                                           &
                      + ( XX )                                                          &
                        * ( (1.5*A - 0.5*C + lambda)*Partb*Partb + Parta/(A - C) )      &
                      + (YY + ZZ)                                                       &
                        * ( (1.5*C - 0.5*A + lambda)*Partc*Partc + Parta/(A - C) ) )



    END IF


END IF

MacLaurin_Potential  = POTENTIAL

END FUNCTION MacLaurin_Potential





!+201+##########################################################
!
!   MacLaurin_Potential - Calulates the potential at a given location for the MacLaurin Spheriod
!
!################################################################
SUBROUTINE MacLaurin_Potential_Sub(r, theta, phi, pot)


REAL(KIND = idp), INTENT(IN)                                ::  r, theta, phi
REAL(KIND = idp), INTENT(OUT)                               ::  Pot


REAL(KIND = idp)                                            ::  X, Y, Z
REAL(KIND = idp)                                            ::  XX, YY, ZZ

REAL(idp)                                                   ::  A, B, C
REAL(idp)                                                   ::  AA, BB, CC

REAL(KIND = idp)                                            ::  VAL, &
                                                                rsqr,                &
                                                                eccsqr, eccmod, &
                                                                Parta, Partb, Partc, &
                                                                lambda, POTENTIAL

CALL Calc_MLS_ABCs( A, B, C, AA, BB, CC )


IF ( MLS_SemiMajor == MLS_SemiMinor ) THEN

    IF ( r >= MLS_SemiMajor ) THEN

!        POTENTIAL = -GM/r

    ELSE

    END IF

ELSE

    eccsqr = MLS_Ecc*MLS_Ecc


    !   Calculate Cartesian Coordinates
    X = r * cos(phi) * sin(theta)
    Y = r * sin(phi) * sin(theta)
    Z = r * cos(theta)

    XX = X*X
    YY = Y*Y
    ZZ = Z*Z

    
    VAL = XX/AA + YY/BB + ZZ/CC

    !
    !   VAL > 1 -> point is outside.
    !   VAL < 1 -> point is inside.
    !   VAL = 1 -> point on surface.
    !
!    PRINT*,"VAL",Val
    IF (VAL .LE. 1.0_idp) THEN

        !
        !   Equations for potential inside homegenous spheroid
        !
        !       REF:    Ellipsoidal Figures of Equilibrium - Chandrasekhar - Chapter 3
        !
        IF ( MLS_SphereType == iMLS_Oblate ) THEN

            eccmod = sqrt(1- eccsqr)/(eccsqr*MLS_Ecc)*ASIN(MLS_Ecc)

            Parta = eccmod - (1 - eccsqr)/eccsqr

            Partb = 2/eccsqr - 2*eccmod


            POTENTIAL = -pi*Grav_Constant_G*MLS_Rho * ((2*AA - XX - YY)*Parta + (CC - ZZ)*Partb)
!            PRINT*,"HERERERERERE"
!            PRINT*,"A & C: ",AA,CC
!            print*,"X, Y, & Z: ",XX,YY,ZZ
!            PRINT*,"Parts: ",parta,partb

        ELSE IF (MLS_SphereType == iMLS_Prolate) THEN

            eccmod = (1 - eccsqr)/(eccsqr*MLS_Ecc) * LOG((1+MLS_Ecc)/(1-MLS_Ecc))

            Parta = eccmod - 2*(1 - eccsqr)/eccsqr

            Partb = 1/eccsqr - 0.5_idp*eccmod


            POTENTIAL = pi*(2*CC - YY - ZZ)*Partb + (AA - XX)*Parta


        END IF


    ELSE IF (VAL > 1.0_idp) THEN


        !
        !   Equations for potential outside homegenous spheroid
        !
        !       REF:    Ellipsoidal Figures of Equilibrium - Chandrasekhar - Chapter 3 -
        !
        lambda = MacLaurin_Root_Finder(r, theta, phi, AA, CC)

        rsqr = XX + YY

        Parta = (sqrt(1-eccsqr)/(eccsqr*MLS_Ecc))*asin(MLS_Ecc) - (1-eccsqr)/eccsqr;

        Partb = 2/eccsqr - 2*((1-eccsqr)/(MLS_Ecc*eccsqr))*asin(MLS_Ecc);

        Partc = pi/sqrt(AA-CC) - 2/sqrt(AA-CC)*atan(sqrt((CC+lambda)/(AA-CC)));


        IF ( MLS_SphereType == iMLS_Oblate ) THEN

            POTENTIAL = -pi*Grav_Constant_G*MLS_Rho                                  &
                        * AA * C                                                    &
                        * ( (1 + rsqr/(2*(CC-AA)) - ZZ/(CC-AA))*Partc                    &
                            - rsqr*sqrt(CC+lambda)/((CC-AA)*(AA+lambda))              &
                            - ZZ*( 2/((AA+lambda)*sqrt(CC+lambda))                  &
                                   - 2*sqrt(CC+lambda)/((CC-AA)*(AA+lambda)) ) );


        ELSE IF ( MLS_SphereType == iMLS_Prolate ) THEN


            POTENTIAL = A * CC                                                              &
                      * ( - Parta                                                           &
                          - Partb                                                           &
                          + ( XX )                                                          &
                            * ( (1.5*A - 0.5*C + lambda)*Partb*Partb + Parta/(A - C) )      &
                          + (YY + ZZ)                                                       &
                            * ( (1.5*C - 0.5*A + lambda)*Partc*Partc + Parta/(A - C) ) )



        END IF


    END IF

END IF

Pot  = POTENTIAL

END Subroutine MacLaurin_Potential_Sub




!+201+##########################################################
!
!   MacLaurin_Potential - Calulates the potential at a given location for the MacLaurin Spheriod
!
!################################################################
SUBROUTINE MacLaurin_Potential_Sub_B(r, theta, phi, pot)


REAL(KIND = idp), INTENT(IN)                                ::  r, theta, phi
REAL(KIND = idp), INTENT(OUT)                               ::  Pot


REAL(KIND = idp)                                            ::  X, Y, Z
REAL(KIND = idp)                                            ::  XX, YY, ZZ

REAL(idp)                                                   ::  A, B, C
REAL(idp)                                                   ::  AA, BB, CC

REAL(KIND = idp)                                            ::  VAL, &
                                                                rsqr,                &
                                                                eccsqr, eccmod, &
                                                                Parta, Partb, Partc, &
                                                                lambda, POTENTIAL

CALL Calc_MLS_ABCs( A, B, C, AA, BB, CC )

PRINT*,"loc: ",r,theta,phi
print*,"ABCs:",A,B,C,AA,BB,CC
PRINT*,"Major/Minor",MLS_SemiMajor,MLS_SemiMinor


IF ( MLS_SemiMajor == MLS_SemiMinor ) THEN

    IF ( r >= MLS_SemiMajor ) THEN

!        POTENTIAL = -GM/r

    ELSE

    END IF

ELSE
    print*,"MLS_Ecc",MLS_Ecc
    eccsqr = MLS_Ecc*MLS_Ecc

    !   Calculate Cartesian Coordinates
    X = r * cos(phi) * sin(theta)
    Y = r * sin(phi) * sin(theta)
    Z = r * cos(theta)

    XX = X*X
    YY = Y*Y
    ZZ = Z*Z


    VAL = XX/AA + YY/BB + ZZ/CC
    !
    !   VAL > 1 -> point is outside.
    !   VAL < 1 -> point is inside.
    !   VAL = 1 -> point on surface.
    !

    PRINT*,"Val",Val
    IF (VAL .LE. 1.0_idp) THEN

        !
        !   Equations for potential inside homegenous spheroid
        !
        !       REF:    Ellipsoidal Figures of Equilibrium - Chandrasekhar - Chapter 3
        !
        IF ( MLS_SphereType == iMLS_Oblate) THEN

            eccmod = sqrt(1- eccsqr)/(eccsqr*MLS_Ecc)*ASIN(MLS_Ecc)

            Parta = eccmod - (1 - eccsqr)/eccsqr

            Partb = 2/eccsqr - 2*eccmod


            POTENTIAL = -pi*Grav_Constant_G*MLS_Rho * ((2*AA - XX - YY)*Parta + (CC - ZZ)*Partb)

        ELSE IF (MLS_SphereType == iMLS_Prolate) THEN

            eccmod = (1 - eccsqr)/(eccsqr*MLS_Ecc) * LOG((1+MLS_Ecc)/(1-MLS_Ecc))

            Parta = eccmod - 2*(1 - eccsqr)/eccsqr

            Partb = 1/eccsqr - 0.5_idp*eccmod


            POTENTIAL = pi*Grav_Constant_G*MLS_Rho*(2*CC - YY - ZZ)*Partb + (AA - XX)*Parta


        END IF


    ELSE IF (VAL > 1.0_idp) THEN

        !
        !   Equations for potential outside homegenous spheroid
        !
        !       REF:    Ellipsoidal Figures of Equilibrium - Chandrasekhar - Chapter 3 -
        !
        lambda = MacLaurin_Root_Finder(r, theta, phi, AA, CC)

        rsqr = XX + YY

        Parta = (sqrt(1-eccsqr)/(eccsqr*MLS_Ecc))*asin(MLS_Ecc) - (1-eccsqr)/eccsqr;

        Partb = 2/eccsqr - 2*((1-eccsqr)/(MLS_Ecc*eccsqr))*asin(MLS_Ecc);

        Partc = pi/sqrt(AA-CC) - 2/sqrt(AA-CC)*atan(sqrt((CC+lambda)/(AA-CC)));


        IF ( MLS_SphereType == iMLS_Oblate ) THEN




            POTENTIAL = -pi*Grav_Constant_G*MLS_Rho                                  &
                        * AA * C                                                    &
                        * ( (1 + rsqr/(2*(CC-AA)) - ZZ/(CC-AA))*Partc                    &
                            - rsqr*sqrt(CC+lambda)/((CC-AA)*(AA+lambda))              &
                            - ZZ*( 2/((AA+lambda)*sqrt(CC+lambda))                  &
                                   - 2*sqrt(CC+lambda)/((CC-AA)*(AA+lambda)) ) );


        ELSE IF ( MLS_SphereType == iMLS_Prolate) THEN


            POTENTIAL = -pi*Grav_Constant_G*MLS_Rho                                   &
                      * A * CC                                                              &
                      * ( - Parta                                                           &
                          - Partb                                                           &
                          + ( XX )                                                          &
                            * ( (1.5*A - 0.5*C + lambda)*Partb*Partb + Parta/(A - C) )      &
                          + (YY + ZZ)                                                       &
                            * ( (1.5*C - 0.5*A + lambda)*Partc*Partc + Parta/(A - C) ) )



        END IF


    END IF

END IF

Pot  = POTENTIAL

END Subroutine MacLaurin_Potential_Sub_B










END MODULE External_MLS_Solution_Module

