   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE MacLaurin_Module                                                             !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   Initialize_MacLaurin_Sources                                        !##!
!##!                                                                                !##!
!##!    +201+   Test_Solution_MacLaurin                                             !##!
!##!                                                                                !##!
!##!    +301+   MacLaurin_Radius                                                    !##!
!##!    +302+   MacLaurin_Root_Finder                                               !##!
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

USE Units_Module, &
                        ONLY :  Grav_Constant_G,        &
                                Gram,                  &
                                Centimeter,             &
                                C_Square

USE Poseidon_Numbers_Module, &
                        ONLY :  pi

USE Poseidon_Parameters, &
                        ONLY :  Verbose_Flag

USE Functions_Mapping, &
                        ONLY :  Map_To_X_Space,         &
                                Map_From_X_Space


USE Variables_Functions, &
                        ONLY :  Potential_Solution

IMPLICIT NONE

INTEGER                                                     :: Spheroid_Type_Flag
CHARACTER(LEN=7)                                            :: Spheroid_Name
REAL(idp)                                                   :: SemiMajor_Axis
REAL(idp)                                                   :: SemiMinor_Axis
REAL(idp)                                                   :: Density

REAL(idp)                                                   ::  A, B, C
REAL(idp)                                                   ::  AA, BB, CC


CONTAINS



!+101+###########################################################################!
!                                                                                !
!                  Initialize_MacLaurin_Sources                                  !
!                                                                                !
!################################################################################!
SUBROUTINE Initialize_MacLaurin_Sources( Rho_O, Spheroid_Type,              &
                                         SemiMajor, SemiMinor,              &
                                         Num_Quad,                          &
                                         R_Quad, T_Quad, P_Quad,            &
                                         Left_Limit, Right_Limit,           &
                                         Num_Elem,                          &
                                         R_Locs, T_Locs, P_Locs,            &
                                         Output_E, Output_S, Output_Si      )


REAL(idp),  INTENT(IN)                                                      :: Rho_O
CHARACTER(LEN=7), INTENT(IN)                                                :: Spheroid_Type

REAL(idp),  INTENT(IN)                                                      :: SemiMajor
REAL(idp),  INTENT(IN)                                                      :: SemiMinor

INTEGER,    INTENT(IN), DIMENSION( 1:3 )                                    :: Num_Quad
REAL(idp),  INTENT(IN), DIMENSION( 1:Num_Quad(1) )                          :: R_Quad
REAL(idp),  INTENT(IN), DIMENSION( 1:Num_Quad(2) )                          :: T_Quad
REAL(idp),  INTENT(IN), DIMENSION( 1:Num_Quad(3) )                          :: P_Quad
REAL(idp),  INTENT(IN)                                                      :: Left_Limit
REAL(idp),  INTENT(IN)                                                      :: Right_Limit

INTEGER,    INTENT(IN), DIMENSION( 1:3 )                                    :: Num_Elem
REAL(idp),  INTENT(IN), DIMENSION( 1:Num_Elem(1)+1 )                        :: R_Locs
REAL(idp),  INTENT(IN), DIMENSION( 1:Num_Elem(2)+1 )                        :: T_Locs
REAL(idp),  INTENT(IN), DIMENSION( 1:Num_Elem(3)+1 )                        :: P_Locs

REAL(idp),  INTENT(OUT),DIMENSION( 1:Num_Quad(1)*Num_Quad(2)*Num_Quad(3),   &
                                   0:Num_Elem(1)-1,                         &
                                   0:Num_Elem(2)-1,                         &
                                   0:Num_Elem(3)-1  )                       ::  Output_E,   &
                                                                                Output_S

REAL(idp), INTENT(OUT), DIMENSION( 1:Num_Quad(1)*Num_Quad(2)*Num_Quad(3),   &
                                   0:Num_Elem(1)-1,                         &
                                   0:Num_Elem(2)-1,                         &
                                   0:Num_Elem(3)-1,                         &
                                   1:3                  )                   ::  Output_Si



INTEGER                                                     ::  re, te, pe
INTEGER                                                     ::  rd, td, pd
INTEGER                                                     ::  l, m, d
INTEGER                                                     ::  Tot_Quad

REAL(idp),  DIMENSION( 1:Num_Quad(1) )                      ::  RQ_Loc
REAL(idp),  DIMENSION( 1:Num_Quad(2) )                      ::  TQ_Loc
REAL(idp),  DIMENSION( 1:Num_Quad(3) )                      ::  PQ_Loc

REAL(idp),  DIMENSION( 1:Num_Quad(1) )                      ::  rsqr

REAL(idp),  DIMENSION( 1:Num_Quad(2) )                      ::  sinsqr_t
REAL(idp),  DIMENSION( 1:Num_Quad(2) )                      ::  cossqr_t

REAL(idp),  DIMENSION( 1:Num_Quad(3) )                      ::  sinsqr_p
REAL(idp),  DIMENSION( 1:Num_Quad(3) )                      ::  cossqr_p


REAL(idp),  DIMENSION( 1:Num_Quad(1) )                      ::  Cur_R_Locs
REAL(idp),  DIMENSION( 1:Num_Quad(2) )                      ::  Cur_T_Locs
REAL(idp),  DIMENSION( 1:Num_Quad(3) )                      ::  Cur_P_Locs

INTEGER                                                     ::  Here

REAL(idp)                                                   ::  Energy
REAL(idp)                                                   ::  Pressure
REAL(idp)                                                   ::  Spec_Ent

REAL(idp)                                                   ::  Value

Tot_Quad = Num_Quad(1) * Num_Quad(2) * Num_Quad(3)
SemiMajor_Axis = SemiMajor
SemiMinor_Axis = SemiMinor


RQ_Loc = Map_To_X_Space(Left_Limit, Right_Limit, R_Quad )
TQ_Loc = Map_To_X_Space(Left_Limit, Right_Limit, T_Quad )
PQ_Loc = Map_To_X_Space(Left_Limit, Right_Limit, P_Quad )


IF ( Spheroid_Type == 'Prolate') THEN
    Spheroid_Type_Flag  = 2
    Spheroid_Name       = 'Prolate'

    A = SemiMajor_Axis
    B = SemiMinor_Axis
    C = B
    
ELSE
    Spheroid_Type_Flag  = 1
    Spheroid_Name       = 'Oblate '

    A = SemiMajor_Axis
    B = A
    C = SemiMinor_Axis

END IF


AA = A*A
BB = B*B
CC = C*C


IF ( Verbose_Flag ) THEN
    PRINT*,"Initalizing Maclaurin Spheroid with the following parameters"
    WRITE(*,'(A,A)')      "Spheroid Type  : ", Spheroid_Name
    WRITE(*,'(A,ES12.5)') "Semimajor Axis : ", SemiMajor_Axis
    WRITE(*,'(A,ES12.5)') "Semiminor Axis : ", SemiMinor_Axis
    WRITE(*,'(A,ES12.5)') "Density        : ", Rho_O
END IF


Density  = Rho_O
Pressure = 0.0_idp
Energy   = 0.0_idp

!Pressure = kappa * Density**Gamma
!Energy   = Pressure / (Gamma - 1.0_dip )
Spec_Ent = C_Square + (Energy + Pressure)/Density



DO re = 1,Num_Elem(1)
DO te = 1,Num_Elem(2)
DO pe = 1,Num_Elem(3)

    Cur_R_Locs(:) = Map_From_X_Space(R_locs(re),R_locs(re+1),RQ_Loc)
    Cur_T_Locs(:) = Map_From_X_Space(T_locs(te),T_locs(te+1),TQ_Loc)
    Cur_P_Locs(:) = Map_From_X_Space(P_locs(pe),P_locs(pe+1),PQ_Loc)

    rsqr(:)     = Cur_R_Locs(:) * Cur_R_Locs(:)
    sinsqr_t(:) = SIN( Cur_T_Locs(:) ) * SIN( CUR_T_LOCS(:) )
    cossqr_t(:) = COS( Cur_T_Locs(:) ) * COS( CUR_T_LOCS(:) )
    sinsqr_p(:) = SIN( Cur_P_Locs(:) ) * SIN( CUR_P_LOCS(:) )
    cossqr_p(:) = COS( Cur_P_Locs(:) ) * COS( CUR_P_LOCS(:) )


    DO rd = 1,Num_Quad(1)
    DO td = 1,Num_Quad(2)
    DO pd = 1,Num_Quad(3)

        Here = (rd-1) * Num_Quad(3) * Num_Quad(2)       &
             + (td-1) * Num_Quad(3)                     &
             + pd

        Value = ( rsqr(rd) * cossqr_p(pd) * sinsqr_t(td) ) / AA      &
              + ( rsqr(rd) * sinsqr_p(pd) * sinsqr_t(td) ) / BB      &
              + ( rsqr(rd) * cossqr_t(td) ) / CC


        
        IF ( Value .LE. 1.0_idp ) THEN
            Output_E(Here,re-1,te-1,pe-1) = Density * Spec_Ent - Pressure
        ELSE
            Output_E(Here,re-1,te-1,pe-1) = 0.0_idp
        END IF

    END DO ! pd
    END DO ! td
    END DO ! rd
END DO ! pe
END DO ! te
END DO ! re


Potential_Solution => MacLaurin_Potential


END SUBROUTINE Initialize_MacLaurin_Sources








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
REAL(KIND = idp)                                            ::  VAL, &
                                                                rsqr, psqr, &
                                                                ecc, eccsqr, eccmod, &
                                                                Parta, Partb, Partc, &
                                                                lambda, POTENTIAL,h



!   Calculate Eccentricity !
ecc =   sqrt(1 - CC/AA)
eccsqr = ecc*ecc

PRINT*,"Eccs",ecc,eccsqr
PRINT*,"AA,BB,CC",AA,BB,CC

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
    IF ( Spheroid_Type_Flag == 1 ) THEN

        eccmod = sqrt(1- eccsqr)/(eccsqr*ecc)*ASIN(ecc)

        Parta = eccmod - (1 - eccsqr)/eccsqr

        Partb = 2/eccsqr - 2*eccmod


        POTENTIAL = pi*Grav_Constant_G*Density * ((2*AA - XX - YY)*Parta + (CC - ZZ)*Partb)


    ELSE IF (Spheroid_Type_Flag == 2) THEN

        eccmod = (1 - eccsqr)/(eccsqr*ecc) * LOG((1+ecc)/(1-ecc))

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
    lambda = MacLaurin_Root_Finder(r, theta, phi)

    Parta = ATAN(sqrt((AA - CC)/(CC + lambda)))

    Partb = 1/sqrt(AA - CC)

    Partc = Partb*Partb*Partb



    PRINT*,"PArts",Parta,partb, Lambda

    IF ( Spheroid_Type_Flag == 1 ) THEN




        POTENTIAL = -pi * Density * Grav_Constant_G * AA * C                                &
                  * ( 2 * Partb * Parta                                                 &
                     - (XX + YY) * Partc                                                &
                        * (Parta - ( SQRT(AA-CC)*SQRT(CC + lambda) )/(AA + lambda) )    &
                     + 2*(ZZ) * Partc * ( Parta - sqrt(AA - CC)/sqrt(CC + lambda) ) )


    ELSE IF ( Spheroid_Type_Flag == 2 ) THEN


        POTENTIAL = A * CC                                                              &
                  * ( - Parta                                                           &
                      - Partb                                                           &
                      + ( XX )                                                          &
                        * ( (1.5*A - 0.5*C + lambda)*Partb*Partb + Parta/(A - C) )      &
                      + (YY + ZZ)                                                       &
                        * ( (1.5*C - 0.5*A + lambda)*Partc*Partc + Parta/(A - C) ) )



    END IF


END IF

PRINT*,"At the end"
MacLaurin_Potential  = POTENTIAL

Print*,"at the end 2"
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
REAL(KIND = idp)                                            ::  VAL, &
                                                                rsqr, psqr, &
                                                                ecc, eccsqr, eccmod, &
                                                                Parta, Partb, Partc, &
                                                                lambda, POTENTIAL,h



!   Calculate Eccentricity !
ecc =   sqrt(1 - CC/AA)
eccsqr = ecc*ecc


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
    IF ( Spheroid_Type_Flag == 1 ) THEN

        eccmod = sqrt(1- eccsqr)/(eccsqr*ecc)*ASIN(ecc)

        Parta = eccmod - (1 - eccsqr)/eccsqr

        Partb = 2/eccsqr - 2*eccmod


        POTENTIAL = -pi*Grav_Constant_G*Density * ((2*AA - XX - YY)*Parta + (CC - ZZ)*Partb)


    ELSE IF (Spheroid_Type_Flag == 2) THEN

        eccmod = (1 - eccsqr)/(eccsqr*ecc) * LOG((1+ecc)/(1-ecc))

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
    lambda = MacLaurin_Root_Finder(r, theta, phi)

    rsqr = XX + YY

    Parta = (sqrt(1-eccsqr)/(eccsqr*ecc))*asin(ecc) - (1-eccsqr)/eccsqr;

    Partb = 2/eccsqr - 2*((1-eccsqr)/(ecc*eccsqr))*asin(ecc);

    Partc = pi/sqrt(AA-CC) - 2/sqrt(AA-CC)*atan(sqrt((CC+lambda)/(AA-CC)));


    IF ( Spheroid_Type_Flag == 1 ) THEN




        POTENTIAL = -pi*Grav_Constant_G*Density                                  &
                    * AA * C                                                    &
                    * ( (1 + rsqr/(2*(CC-AA)) - ZZ/(CC-AA))*Partc                    &
                        - rsqr*sqrt(CC+lambda)/((CC-AA)*(AA+lambda))              &
                        - ZZ*( 2/((AA+lambda)*sqrt(CC+lambda))                  &
                               - 2*sqrt(CC+lambda)/((CC-AA)*(AA+lambda)) ) );


    ELSE IF ( Spheroid_Type_Flag == 2 ) THEN


        POTENTIAL = A * CC                                                              &
                  * ( - Parta                                                           &
                      - Partb                                                           &
                      + ( XX )                                                          &
                        * ( (1.5*A - 0.5*C + lambda)*Partb*Partb + Parta/(A - C) )      &
                      + (YY + ZZ)                                                       &
                        * ( (1.5*C - 0.5*A + lambda)*Partc*Partc + Parta/(A - C) ) )



    END IF


END IF

Pot  = POTENTIAL

END Subroutine MacLaurin_Potential_Sub




!+301+###########################################################################!
!                                                                                !
!                               MacLaurin_Radius                                 !
!                                                                                !
!################################################################################!
PURE FUNCTION MacLaurin_Radius( theta, phi )


REAL(idp),  INTENT(IN)          :: theta, phi


REAL(idp)                       ::  RADIUS_HERE

REAL(idp)                       ::  MacLaurin_Radius



IF (Spheroid_Type_Flag == 1 ) THEN

    RADIUS_HERE = A * C                                     &
                / SQRT( CC * SIN(theta) * SIN(theta)        &
                      + AA * COS(theta) * COS(theta)     )

ELSE IF (Spheroid_Type_Flag == 2) THEN

    RADIUS_HERE = A * B                                                     &
                * SQRT( BB                                                  &
                        * COS(theta) * COS(theta) * SIN(phi) * SIN(phi)     &
                      + AA                                                  &
                        * ( SIN(phi) * SIN(phi) * SIN(theta) * SIN(theta)   &
                            + COS(theta) * COS(theta)   )   )

END IF


MacLaurin_Radius = RADIUS_HERE

END FUNCTION MacLaurin_Radius




!+302+#################################################################
!
!   MacLaurin_Root_Finder
!
!#######################################################################
PURE FUNCTION MacLaurin_Root_Finder( r, theta, phi )


REAL(KIND = idp),   INTENT(IN)                      ::  r, theta, phi

REAL(KIND = idp)                                    ::  MacLaurin_Root_Finder

REAL(KIND = idp)                                    ::  rsqr, xsqr, ysqr, zsqr

REAL(KIND = idp)                                    ::  LenA, LenB

!! r^2  !!
rsqr = r * r


!! x^2 = r^2 * cos(theta)^2 * sin(phi)^2 !!
xsqr = rsqr * cos(phi) * cos(phi)* sin(theta) * sin(theta)

!! y^2 = r^2 * sin(theta)^2 * sin(phi)^2 !!
ysqr = rsqr * sin(phi) * sin(phi)* sin(theta) * sin(theta)

!! z^2 = r^2 * cos(phi)^2  !!
zsqr = rsqr * cos(theta) * cos(theta)



LenA = AA + CC - xsqr - ysqr - zsqr


IF ( Spheroid_Type_Flag == 1 ) THEN

    LenB = AA*(CC - zsqr) - CC*(xsqr + ysqr)

ELSE IF ( Spheroid_Type_Flag == 2 ) THEN

    LenB = CC*(AA - xsqr) - AA*(ysqr + zsqr)

END IF


MacLaurin_Root_Finder = 0.5 * ( - LenA + sqrt(LenA*LenA - 4*LenB))




END FUNCTION MacLaurin_Root_Finder











END MODULE MacLaurin_Module
