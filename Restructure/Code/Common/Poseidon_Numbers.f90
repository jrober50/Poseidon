MODULE Poseidon_Numbers_Module

USE Poseidon_Kinds_Module, &
            ONLY : idp

IMPLICIT NONE

REAL(KIND = idp), PUBLIC, PARAMETER     :: pi = 4.0_idp*ATAN(1.0_idp)
REAL(KIND = idp), PUBLIC, PARAMETER     :: twopi = 2.0_idp*pi
REAL(KIND = idp), PUBLIC, PARAMETER     :: eps = EPSILON(1.0_idp)


REAL(KIND = idp), PUBLIC, PARAMETER     ::  OneThird = 1.0_idp/3.0_idp
REAL(KIND = idp), PUBLIC, PARAMETER     ::  TwoThirds = 2.0_idp/3.0_idp
REAL(KIND = idp), PUBLIC, PARAMETER     ::  FourThirds = 4.0_idp/3.0_idp
REAL(KIND = idp), PUBLIC, PARAMETER     ::  OneThirtySecond = 1.0_idp/32.0_idp



END MODULE  Poseidon_Numbers_Module
