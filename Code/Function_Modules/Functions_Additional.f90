   !#################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Functions_Additional                                                         !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains functions and subroutine to calculate special functions such as    !##!
!##!    the Lagrange and Legendre Polynomials and Spherical Harmonics, initialize   !##!
!##!    quadratures, perform other useful operations                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +105+   Spherical_Harmonic                                                  !##!
!##!                                                                                !##!
!##!    +201+   Norm_Factor                                                         !##!
!##!    +202+   Factorial                                                           !##!
!##!    +203+   Factorial_Division                                                  !##!
!##!                                                                                !##!
!##!    +401+   MVMULT_CCS                                                          !##!
!##!    +402+   MVMULT_FULL                                                         !##!
!##!    +403+   SVVAD                                                               !##!
!##!                                                                                !##!
!##!    +601+   SphericalHarmonic_OrthogonalityTest                                 !##!
!##!                                                                                !##!
!##!    +801+   Generate_Defined_Coarse_Mesh                                        !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !#################################################################################!

!*D*================================!
!                                   !
!           Dependencies            !
!                                   !
!===================================!
USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Poseidon_Numbers_Module, &
            ONLY :  pi, eps

USE Functions_Math, &
            ONLY :  Legendre_Poly,      &
                    Norm_Factor,        &
                    Spherical_Harmonic

USE Functions_Mapping,  &
            ONLY :  Map_From_X_Space

USE Functions_Quadrature, &
            ONLY : Initialize_LGL_Quadrature


IMPLICIT NONE





CONTAINS






 !+202+################################################################!
!                                                                       !
!       Factorial - Calulcates N!, for use in Legendre_Poly             !
!                                                                       !
 !#####################################################################!
PURE ELEMENTAL REAL(KIND = idp) FUNCTION Factorial(n)

INTEGER, INTENT(IN) :: n


INTEGER             :: i
REAL(KIND = idp)  :: fact


fact = 1
DO i = 2,n
    fact = fact*i
END DO




Factorial = fact


END FUNCTION Factorial







 !+203+################################################################!
!                                                                       !
!       Factorial_Division - Calulcates M!/N!, for use in Legendre_Poly !
!                                                                       !
 !#####################################################################!
PURE ELEMENTAL REAL(KIND = idp) FUNCTION Factorial_Division(M,N)

INTEGER, INTENT(IN) :: N,M


INTEGER             :: i
REAL(KIND = idp)  :: fact



IF (M .GE. N) THEN

    fact = N
    DO i = N+1,M
        fact = fact*i
    END DO

    IF (fact .EQ. 0) THEN
        fact = 1
    END IF

END IF


IF (M < N) THEN

    fact = M

    IF (fact .EQ. 0) THEN
        fact = 1
    END IF

    DO i = M+1, N
        fact = fact*i
    END DO




    fact = 1/fact

END IF

Factorial_Division = fact


END FUNCTION Factorial_Division











!+601+##########################################################################!
!                                                                               !
!       SphericalHarmonic_OrthogonalityTest - Shows Orthogonality of Spherical  !
!                                       Harmonic Functions                      !
!                                                                               !
 !##############################################################################!
SUBROUTINE SphericalHarmonic_OrthogonalityTest(L_LIMIT, L_SPECIFIC, M_SPECIFIC)


!REAL(KIND = idp)                                :: BC_Integralb
!COMPLEX(KIND = idp)                                :: SphereHarmonic_OrthogTest



INTEGER, INTENT(IN)                             :: L_LIMIT, L_SPECIFIC, M_SPECIFIC


INTEGER                                         ::  l, m
INTEGER                                         ::  te, td, pe, pd


INTEGER                                         ::  T_Degree, P_Degree

INTEGER                                         ::  T_SUB_ELEMENTS, P_SUB_ELEMENTS


REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)     ::  P_xlocs, P_locs, P_weights, &
                                                    T_xlocs, T_locs, T_weights, &
                                                    Sub_tlocs, Sub_plocs

REAL(KIND = idp)                                ::  deltasubt, deltasubp, deltat, deltap


COMPLEX(KIND = idp)                             :: TMP_VALUE


T_SUB_ELEMENTS = 32
P_SUB_ELEMENTS = 32

T_Degree = 5
P_Degree = 5




ALLOCATE(P_xlocs(0:P_Degree), P_locs(0:P_Degree), P_weights(0:P_Degree))
ALLOCATE(T_xlocs(0:T_Degree), T_locs(0:T_Degree), T_weights(0:T_Degree))
ALLOCATE(Sub_tlocs(0:T_SUB_ELEMENTS),Sub_plocs(0:P_SUB_ELEMENTS))

CALL Initialize_LGL_Quadrature(P_Degree, P_xlocs, P_weights)
CALL Initialize_LGL_Quadrature(T_Degree, T_xlocs, T_weights)

deltasubt = (pi)/REAL(T_SUB_ELEMENTS)
deltasubp = (2*pi)/REAL(P_SUB_ELEMENTS)

DO te = 0, T_SUB_ELEMENTS
    Sub_tlocs(te) = 0.0_idp + te*deltasubt
END DO

DO pe = 0,P_SUB_ELEMENTS
    Sub_plocs(pe) = 0.0_idp + pe*deltasubp
END DO




PRINT*,"L_SPECIFIC, l, M_SPECIFIC, m, Integral Value"
DO l = 0,L_LIMIT
    DO m = -l,l

        TMP_VALUE = 0.0_idp

        DO te = 0, T_SUB_ELEMENTS - 1

            T_locs = Map_From_X_Space(Sub_tlocs(te), Sub_tlocs(te+1), T_xlocs)
            deltat = Sub_tlocs(te+1)-Sub_tlocs(te)

            DO pe = 0, P_SUB_ELEMENTS - 1

                P_locs = Map_From_X_Space(Sub_plocs(pe), Sub_plocs(pe+1), P_xlocs)
                deltap = Sub_plocs(pe+1) - Sub_plocs(pe)

                DO td = 0, T_Degree

                    DO pd = 0, P_Degree

                        TMP_VALUE = TMP_VALUE + T_weights(td)*P_weights(pd)                                             &
                                                * Spherical_Harmonic(L_SPECIFIC, M_SPECIFIC, T_locs(td), P_locs(pd))    & ! Ylm
                                                * (-1.0_idp**M)*Spherical_Harmonic(L, -M, T_locs(td), P_locs(pd))  & ! Ylm* = (-1)^m Yl-m
                                                * sin(T_locs(td))                                                       &
                                                * deltat/2.0_idp*deltap/2.0_idp


                    END DO
                END DO
            END DO
        END DO

        PRINT*, L_SPECIFIC, l, M_SPECIFIC, m, TMP_VALUE


    END DO
END DO



END SUBROUTINE SphericalHarmonic_OrthogonalityTest















END MODULE Functions_Additional
