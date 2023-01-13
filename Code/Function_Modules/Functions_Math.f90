   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Functions_Math                                                               !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   Lagrange_Poly                                                       !##!
!##!    +102+   Lagrange_Poly_Deriv                                                 !##!
!##!    +103+   Lagrange_Second_Deriv                                               !##!
!##!                                                                                !##!
!##!    +201+   Legendre_Poly                                                       !##!
!##!    +202+   Norm_Factor                                                         !##!
!##!    +203+   Factorial                                                           !##!
!##!                                                                                !##!
!##!    +301+   Spherical_Harmonic                                                  !##!
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

USE Poseidon_Numbers_Module, &
            ONLY :  pi




IMPLICIT NONE

CONTAINS



 !+101+################################################################!
!                                                                       !
!   Lagrange_Poly - Calculates the value of the Lagrange Polynomial     !
!                                                                       !
!                                                                       !
 !#####################################################################!
PURE FUNCTION Lagrange_Poly(x, Ord, xlocs)

REAL(idp),  INTENT(IN)                          :: x
INTEGER,    INTENT(IN)                          :: Ord
REAL(idp),  INTENT(IN), DIMENSION(0:Ord)        :: xlocs


INTEGER                                         :: i,j
REAL(idp), DIMENSION(0:Ord)                     :: tmp
REAL(idp), DIMENSION(0:Ord)                     :: Lagrange_Poly



tmp(:) = 1.0_idp

DO j = 0,Ord
    DO i = 0,Ord
        IF (i .NE. j) THEN
            tmp(i) = tmp(i) * (x - xlocs(j))/(xlocs(i) - xlocs(j))
        END IF
    END DO
END DO



Lagrange_Poly = tmp


END FUNCTION Lagrange_Poly





 !+102+################################################################!
!                                                                       !
!   Lagrange_Poly_Deriv - Calculates the value of the first derivative  !
!                    of Lagrange Polynomials.                           !
!                                                                       !
!                                                                       !
 !#####################################################################!
PURE FUNCTION Lagrange_Poly_Deriv(x, Ord, xlocs)

INTEGER, INTENT(IN)                                 :: Ord
REAL(idp), INTENT(IN)                        :: x
REAL(idp), INTENT(IN), DIMENSION(0:Ord)      :: xlocs



INTEGER                                             :: i,j,k
REAL(idp), DIMENSION(0:Ord)                  :: tmp,tmp2
REAL(idp), DIMENSION(0:Ord)                  :: Lagrange_Poly_Deriv


tmp(:) = 1.0_idp
tmp2(:) = 0.0_idp

DO j = 0,Ord
    DO i = 0,Ord
        IF (i .NE. j) THEN

            tmp(i) = 1/(xlocs(i) - xlocs(j))

            DO k = 0,Ord

                IF (k .NE. i .AND. k .NE. j) THEN
                    tmp(i) = tmp(i) * (x - xlocs(k))/(xlocs(i) - xlocs(k))
                END IF

            END DO
            tmp2(i) = tmp2(i) + tmp(i)
        END IF
    END DO
END DO



Lagrange_Poly_Deriv = tmp2


END FUNCTION Lagrange_Poly_Deriv





 !+103+################################################################!
!                                                                       !
!   Lagrange_Second_Deriv - Calculates the value of the first derivative  !
!                    of Lagrange Polynomials.                           !
!                                                                       !
!                                                                       !
 !#####################################################################!
PURE FUNCTION Lagrange_Second_Deriv(x, Ord, xlocs)

INTEGER, INTENT(IN)                                 :: Ord
REAL(idp), INTENT(IN)                        :: x
REAL(idp), INTENT(IN), DIMENSION(0:Ord)      :: xlocs



INTEGER                                             :: i, j, k, l
REAL(idp), DIMENSION(0:Ord)                  :: tmp,tmp2
REAL(idp), DIMENSION(0:Ord)                  :: Lagrange_Second_Deriv





IF (Ord == 1) THEN

    LAGRANGE_SECOND_DERIV = 0.0_idp


ELSE



    tmp(:) = 1.0_idp
    tmp2(:) = 0.0_idp

    DO i = 0,Ord

        DO k = 0,Ord

            IF (k .NE. i) THEN


                DO l = 0,Ord

                    IF (l .NE. i .AND. l .NE. k) THEN

                        tmp(i) = (1/(xlocs(i) - xlocs(k)))*(1/(xlocs(i) - xlocs(l)))


                        DO j = 0,Ord

                            IF (j .NE. i .AND. j .NE. k .AND. j .NE. l) THEN

                                tmp(i) = tmp(i) * (x - xlocs(j))/(xlocs(i) - xlocs(j))

                            END IF

                        END DO  ! j Loop

                        tmp2(i) = tmp2(i) + tmp(i)

                    END IF

                END DO ! l Loop

            END IF

        END DO  ! k Loop

    END DO  ! i Loop



    Lagrange_Second_Deriv = tmp2


END IF

END FUNCTION Lagrange_Second_Deriv





 !+201+################################################################!
!                                                                       !
!   Legendre_Poly - Calculates the value of the Legendre Polynomial     !
!                   P^m_l(cos(theta))                                   !
!                                                                       !
 !#####################################################################!
PURE FUNCTION Legendre_Poly(l,m,num_points,theta)


!  Returns array that conatins the values P^m_n(cos(theta)) for n = 0:l
!  If m > n then returns 0 as the poly doesn't exist there.

INTEGER, INTENT(IN)                                 :: l,m, num_points
REAL(idp),  INTENT(IN), DIMENSION(1:num_points)     :: theta
REAL(idp), DIMENSION(1:num_points)                  :: Legendre_Poly

INTEGER                                             :: i, n
REAL(idp)                                           :: factor, normfactor
REAL(idp), DIMENSION(1:num_points)                  :: costheta, sqrfactor
REAL(idp), DIMENSION(0:l,1:num_points)              :: Plm


n = abs(m)
costheta = cos(theta)
sqrfactor = sqrt(1.0_idp - costheta*costheta)



Plm(:,:) = 0.0_idp




IF (n <= l) THEN

    Plm(n,:) = 1.0_idp
    factor = 1.0_idp

    DO i = 1,n
        Plm(n,:) = -Plm(n,:)*factor*sqrfactor(:)
        factor = factor + 2.0_idp
    END DO
END IF



IF (n + 1 <= l) THEN
    Plm(n+1,:) = costheta(:)*(2.0_idp*n + 1.0_idp) * Plm(n,:)
END IF



DO i = n+2,l

    Plm(i,:) = ( ( 2*i - 1)*costheta(:)*Plm(i-1,:) + (-i - n +1)*Plm(i-2,:) )/(i - n)

END DO





If (m < 0) THEN
    normfactor = ((-1)**n)*(1.0_idp*Factorial(l-n))/(1.0_idp*Factorial(l+n))

    Legendre_Poly(:) = normfactor*Plm(l,:)
ELSE
    Legendre_Poly(:) = Plm(l,:)
END IF




END FUNCTION Legendre_Poly



 !+201+################################################################!
!                                                                       !
!   Legendre_Poly - Calculates the value of the Legendre Polynomial     !
!                   P^m_l(cos(theta))                                   !
!                                                                       !
 !#####################################################################!
FUNCTION Legendre_Poly_Array(l,m,num_points,theta)


!  Returns array that conatins the values P^m_n(cos(theta)) for n = 0:l
!  If m > n then returns 0 as the poly doesn't exist there.

INTEGER, INTENT(IN)                                     :: l,m, num_points
REAL(idp),  INTENT(IN), DIMENSION(1:num_points)  :: theta
REAL(idp), DIMENSION(0:l,1:num_points)           :: Legendre_Poly_Array

INTEGER                                         :: i, n
REAL(idp)                                :: factor, normfactor
REAL(idp), DIMENSION(1:num_points)       :: costheta, sqrfactor
REAL(idp), DIMENSION(0:l,1:num_points)   :: Plm


n = abs(m)
costheta = DCOS(theta)
sqrfactor = sqrt(1.0_idp - costheta*costheta)



Plm(:,:) = 0.0_idp

IF (n <= l) THEN

    Plm(n,:) = 1.0_idp
    factor = 1.0_idp

    DO i = 1,n
        Plm(n,:) = -Plm(n,:)*factor*sqrfactor(:)
        factor = factor + 2.0_idp
    END DO
END IF

IF (n + 1 <= l) THEN
    Plm(n+1,:) = costheta(:)*(2.0_idp*n + 1.0_idp) * Plm(n,:)
END IF


DO i = n+2,l

    Plm(i,:) = ( ( 2*i - 1)*costheta(:)*Plm(i-1,:) + (-i - n +1)*Plm(i-2,:) )/(i - n)

END DO



DO i = 0,l
    If (m < 0) THEN
        normfactor = ((-1)**n)*(1.0_idp*Factorial(l-n))/(1.0_idp*Factorial(l+n))

        Legendre_Poly_Array(l,:) = normfactor*Plm(l,:)
    ELSE
        Legendre_Poly_Array(l,:) = Plm(l,:)
    END IF
END DO ! i Loop


END FUNCTION Legendre_Poly_Array






 !+202+####################################################!
!                                                           !
!       Norm_Factor - Calculate normalization factor for    !
!                       spherical harmonics                 !
!                                                           !
 !#########################################################!
PURE ELEMENTAL FUNCTION Norm_Factor(l,m)

INTEGER, INTENT(IN)         :: l,m
REAL(idp)            :: Norm_Factor

REAL(idp)            ::  Real_L
REAL(idp)            ::  Real_M

Real_L = REAL( l, Kind = idp)
Real_M = REAL( m, Kind = idp)


!Norm_Factor = 1.0_idp
Norm_Factor = sqrt( ( 2.0_idp * Real_L + 1.0_idp)               &
                    * REAL( Factorial(l-m), Kind = idp )        &
                    /( 4.0_idp                                  &
                        * pi                                    &
                        * REAL(Factorial(l+m), Kind = idp) )    &
                    )


END FUNCTION Norm_Factor





 !+203+####################################################!
!                                                           !
!       Sqrt_Factor - Calculate normalization factor for    !
!                       spherical harmonic derivatives      !
!                                                           !
 !#########################################################!
PURE ELEMENTAL FUNCTION Sqrt_Factor(l,m)

INTEGER, INTENT(IN)         :: l,m
REAL(idp)            :: Sqrt_Factor

REAL(idp)            ::  Real_L
REAL(idp)            ::  Real_M

Real_L = REAL( l, Kind = idp)
Real_M = REAL( m, Kind = idp)

!Norm_Factor = 1.0_idp
Sqrt_Factor = sqrt( ( 2.0_idp * Real_L + 1.0_idp ) / ( 2.0_idp * Real_L + 1.0_idp )     &
                     * ( Real_L - Real_M ) * (Real_L + Real_M ) )

END FUNCTION Sqrt_Factor




 !+204+################################################################!
!                                                                       !
!       Factorial - Calulcates N!, for use in Legendre_Poly             !
!                                                                       !
 !#####################################################################!
PURE ELEMENTAL REAL(idp) FUNCTION Factorial(n)

INTEGER, INTENT(IN) :: n


INTEGER             :: i
REAL(idp)  :: fact


fact = 1
DO i = 2,n
    fact = fact*i
END DO




Factorial = fact


END FUNCTION Factorial




 !+301+################################################################!
!                                                                       !
!   Spherical_Harmonic - Calculates the value of the spherical harmonic,!
!                       Y^M_L(Theta, Phi). Uses Legendre_Poly           !
!                                                                       !
!           Output - 2L Value array - (Real, Imaginary)                 !
!                                                                       !
 !#####################################################################!
PURE FUNCTION Spherical_Harmonic(l,m,theta,phi)




INTEGER, INTENT(IN)                         :: l,m
REAL(idp), INTENT(IN)                :: theta, phi

REAL(idp), DIMENSION(0:0)            :: Plm
COMPLEX(idp)                         :: Spherical_Harmonic


Plm = Legendre_Poly(l,m,1,[theta])
Spherical_Harmonic = Norm_Factor(l,m)*Plm(0)*EXP(CMPLX(0,m*phi, KIND = idp))


END FUNCTION Spherical_Harmonic



!+301+################################################################!
!                                                                       !
!   Spherical_Harmonic - Calculates the value of the spherical harmonic,!
!                       Y^M_L(Theta, Phi). Uses Legendre_Poly           !
!                                                                       !
!           Output - 2L Value array - (Real, Imaginary)                 !
!                                                                       !
!#####################################################################!
PURE FUNCTION Spherical_Harmonic_dt(l,m,theta,phi)


COMPLEX(idp)                        ::  Spherical_Harmonic_dt

INTEGER, INTENT(IN)                 ::  l,m
REAL(idp), INTENT(IN)               ::  theta, phi

REAL(idp), DIMENSION(0:0)           ::  Plml
REAL(idp), DIMENSION(0:0)           ::  Plmlm1
COMPLEX(idp)                        ::  Spherical_Harmonicl
COMPLEX(idp)                        ::  Spherical_Harmoniclm1

REAL(idp)                           ::  Real_L
REAL(idp)                           ::  Cot_Val
REAL(idp)                           ::  Csc_Val
COMPLEX(idp)                        ::  Sqrt_Term

Plml   = Legendre_Poly(l,m,1,[theta])
Spherical_Harmonicl   = Norm_Factor(l,m)*Plml(0)*EXP(CMPLX(0,m*phi, KIND = idp))


IF ( l .NE. 0 ) THEN
    Plmlm1 = Legendre_Poly(l-1,m,1,[theta])
    Spherical_Harmoniclm1 = Norm_Factor(l-1,m)*Plmlm1(0)*EXP(CMPLX(0,m*phi, KIND = idp))
ELSE
    Plmlm1 = 0.0_idp
    Spherical_Harmoniclm1 = 0.0_idp
END IF


Real_L = REAL(l, idp)

Cot_Val = 1.0_idp/DTAN(theta)
Csc_Val = 1.0_idp/DSIN(theta)

Sqrt_Term = sqrt( CMPLX( (2.0_idp * Real_L + 1.0_idp)/(2.0_idp*Real_L - 1.0_idp),0.0_idp,idp ) )    &
          * sqrt( CMPLX( (l-m)*(l+m),0.0_idp, idp ) )

Spherical_Harmonic_dt = Real_L * Cot_Val * Spherical_Harmonicl                     &
                      - Sqrt_Term * Csc_Val * Spherical_Harmoniclm1


END FUNCTION Spherical_Harmonic_dt



!+301+################################################################!
!                                                                       !
!   Spherical_Harmonic - Calculates the value of the spherical harmonic,!
!                       Y^M_L(Theta, Phi). Uses Legendre_Poly           !
!                                                                       !
!           Output - 2L Value array - (Real, Imaginary)                 !
!                                                                       !
!#####################################################################!
PURE FUNCTION Spherical_Harmonic_dp(l,m,theta,phi)


COMPLEX(idp)                        ::  Spherical_Harmonic_dp

INTEGER, INTENT(IN)                 ::  l,m
REAL(idp), INTENT(IN)               ::  theta, phi

REAL(idp), DIMENSION(0:0)           ::  Plm
COMPLEX(idp)                        ::  Spherical_Harmonic

Plm = Legendre_Poly(l,m,1,[theta])
Spherical_Harmonic = Norm_Factor(l,m)*Plm(0)*EXP(CMPLX(0,m*phi, KIND = idp))


Spherical_Harmonic_dp = CMPLX(0,m,idp) * Spherical_Harmonic


END FUNCTION Spherical_Harmonic_dp








 !+301+################################################################!
!                                                                       !
!   Spherical_Harmonic - Calculates the value of the spherical harmonic,!
!                       Y^M_L(Theta, Phi). Uses Legendre_Poly           !
!                                                                       !
!           Output - 2L Value array - (Real, Imaginary)                 !
!                                                                       !
 !#####################################################################!
PURE FUNCTION Real_Spherical_Harmonic(l,m,theta,phi)




INTEGER, INTENT(IN)                         :: l,m
REAL(idp), INTENT(IN)                       :: theta, phi

REAL(idp), DIMENSION(0:0)                   ::  Plm
REAL(idp)                                   ::  Nlm
REAL(idp)                                   ::  Am
REAL(idp)                                   :: Real_Spherical_Harmonic


Nlm = Norm_Factor(l,m)
Plm = Legendre_Poly(l,abs(m),1,[theta])

IF ( m < 0 ) THEN
    Am = sqrt(2.0_idp)*DSIN(abs(m)*phi)
ELSE IF ( m == 0 ) THEN
    Am = 1.0_idp
ELSE
    Am = sqrt(2.0_idp)*DCOS(m*phi)
END IF

Real_Spherical_Harmonic = Nlm*Plm(0)*Am


END FUNCTION Real_Spherical_Harmonic



 !+301+################################################################!
!                                                                       !
!   Spherical_Harmonic - Calculates the value of the spherical harmonic,!
!                       Y^M_L(Theta, Phi). Uses Legendre_Poly           !
!                                                                       !
!           Output - 2L Value array - (Real, Imaginary)                 !
!                                                                       !
 !#####################################################################!
PURE FUNCTION Real_Spherical_Harmonic_dt(l,m,theta,phi)




INTEGER, INTENT(IN)                         :: l,m
REAL(idp), INTENT(IN)                       :: theta, phi

REAL(idp), DIMENSION(0:0)                   ::  Plm
REAL(idp), DIMENSION(0:0)                   ::  Plm1
REAL(idp)                                   ::  Nlm
REAL(idp)                                   ::  Am
REAL(idp)                                   ::  Real_Spherical_Harmonic_dt

REAL(idp)                                   ::  Cotan
REAL(idp)                                   ::  Cosec

Nlm = Norm_Factor(l,m)
Plm = Legendre_Poly(l,abs(m),1,[theta])
IF ( abs(m) > l-1 ) THEN
    Plm1 = 0.0_idp
ELSE
    Plm1 = Legendre_Poly(l-1,abs(m),1,[theta])
END IF

Cotan = DCOS(theta)/DSIN(Theta)
Cosec = 1.0_idp/DSIN(theta)


IF ( m < 0 ) THEN
    Am = sqrt(2.0_idp)*DSIN(abs(m)*phi)
ELSE IF ( m == 0 ) THEN
    Am = 1.0_idp
ELSE
    Am = sqrt(2.0_idp)*DCOS(m*phi)
END IF


Real_Spherical_Harmonic_dt = Nlm                   &
                           * ( REAL(l,kind=idp)*Cotan*Plm(0) - REAL(l+m,Kind=idp)*Cosec*Plm1(0) )  &
                           * Am


END FUNCTION Real_Spherical_Harmonic_dt





 !+301+################################################################!
!                                                                       !
!   Spherical_Harmonic - Calculates the value of the spherical harmonic,!
!                       Y^M_L(Theta, Phi). Uses Legendre_Poly           !
!                                                                       !
!           Output - 2L Value array - (Real, Imaginary)                 !
!                                                                       !
 !#####################################################################!
PURE FUNCTION Real_Spherical_Harmonic_dp(l,m,theta,phi)




INTEGER, INTENT(IN)                         :: l,m
REAL(idp), INTENT(IN)                       :: theta, phi

REAL(idp), DIMENSION(0:0)                   ::  Plm
REAL(idp)                                   ::  Nlm
REAL(idp)                                   ::  Am
REAL(idp)                                   ::  Real_Spherical_Harmonic_dp


Nlm = Norm_Factor(l,m)
Plm = Legendre_Poly(l,abs(m),1,[theta])

IF ( m < 0 ) THEN
    Am = -sqrt(2.0_idp)*abs(m)*DCOS(abs(m)*phi)
ELSE IF ( m == 0 ) THEN
    Am = 0.0_idp
ELSE
    Am = sqrt(2.0_idp)*m*DSIN(m*phi)
END IF


Real_Spherical_Harmonic_dp = Nlm*Plm(0)*Am


END FUNCTION Real_Spherical_Harmonic_dp


END MODULE Functions_Math
