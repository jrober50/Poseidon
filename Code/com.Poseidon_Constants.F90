MODULE Poseidon_Constants_Module



IMPLICIT NONE


INTEGER, PARAMETER              :: fdp = selected_int_kind(10)
INTEGER, PARAMETER              :: idp = KIND(1.D0)
REAL(KIND = idp), PARAMETER     :: pi = 4.0_idp*ATAN(1.0_idp)
REAL(KIND = idp), PARAMETER     :: twopi = 2.0_idp*pi
REAL(KIND = idp), PARAMETER     :: eps = EPSILON(1.0_idp)


REAL(KIND = idp)                ::  OneThird = 1.0_idp/3.0_idp
REAL(KIND = idp)                ::  TwoThirds = 2.0_idp/3.0_idp
REAL(KIND = idp)                ::  FourThirds = 4.0_idp/3.0_idp
REAL(KIND = idp)                ::  OneThirtySecond = 1.0_idp/32.0_idp



END MODULE  Poseidon_Constants_Module
