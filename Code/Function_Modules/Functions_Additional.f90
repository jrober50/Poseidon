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
!##!                                                                                !##!
!##!    +202+   Factorial                                                           !##!
!##!    +203+   Factorial_Division                                                  !##!
!##!                                                                                !##!
!##!    +401+   MVMULT_CCS                                                          !##!
!##!    +402+   MVMULT_FULL                                                         !##!
!##!    +403+   SVVAD                                                               !##!
!##!                                                                                !##!
!##!                                                                                !##!
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
            ONLY :  Legendre_Poly

USE Maps_X_Space,  &
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
INTEGER             :: fact
REAL(idp)           :: tmp


IF (M .GE. N) THEN

    fact = N
    DO i = N+1,M
        fact = fact*i
    END DO

    IF (fact == 0) THEN
        fact = 1
    END IF

    tmp = REAL(fact, KIND = idp)
END IF


IF (M < N) THEN

    fact = M

    IF (fact == 0) THEN
        fact = 1
    END IF

    DO i = M+1, N
        fact = fact*i
    END DO

    
    tmp = 1.0_idp/REAL(fact, KIND = idp)
END IF

Factorial_Division = tmp


END FUNCTION Factorial_Division




















END MODULE Functions_Additional
