   !#################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_Additional_Functions_Module                                         !##!
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
!##!    +101+   Lagrange_Poly                                                       !##!
!##!    +102+   Lagrange_First_Deriv                                                !##!
!##!    +103+   Lagrange_Second_Deriv                                               !##!
!##!    +104+   Legendre_Poly                                                       !##!
!##!    +105+   Spherical_Harmonic                                                  !##!
!##!                                                                                !##!
!##!    +201+   Norm_Factor                                                         !##!
!##!    +202+   Factorial                                                           !##!
!##!    +203+   Factorial_Division                                                  !##!
!##!    +204+   Power                                                               !##!
!##!                                                                                !##!
!##!    +301+   Map_To_X_Space                                                      !##!
!##!    +302+   Map_From_X_Space                                                    !##!
!##!                                                                                !##!
!##!    +401+   MVMULT_CCS                                                          !##!
!##!    +402+   MVMULT_FULL                                                         !##!
!##!    +403+   SVVAD                                                               !##!
!##!                                                                                !##!
!##!    +501+   Initialize_LGL_Quadrature                                           !##!
!##!    +502+   Initialize_LGL_Quadrature_Locations                                 !##!
!##!    +503+   Initialize_LG_Quadrature                                            !##!
!##!    +504+   Initialize_LG_Quadrature_Locations                                  !##!
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
USE Poseidon_Constants_Module, &
            ONLY :  idp, pi, eps


USE Poseidon_Parameters, &
            ONLY :  DOMAIN_DIM,             &
                    DEGREE,                 &
                    L_LIMIT,                &
                    NUM_CFA_VARS,           &
                    STF_MAPPING_FLAG

USE Poseidon_Variables_Module,    &
            ONLY :  NUM_R_ELEMENTS,         &
                    VAR_DIM,                &
                    NUM_OFF_DIAGONALS,      &
                    ULM_LENGTH,             &
                    LM_LENGTH,              &
                    Coefficient_Vector,     &
                    M_VALUES,               &
                    rlocs,                  &
                    tlocs,                  &
                    plocs,                  &
                    Matrix_Location

IMPLICIT NONE





CONTAINS
 !+101+################################################################!
!                                                                       !
!   Lagrange_Poly - Calculates the value of the Lagrange Polynomial     !
!                                                                       !
!                                                                       !
 !#####################################################################!
PURE FUNCTION Lagrange_Poly(x, Ord, xlocs)


INTEGER, INTENT(IN)                                 :: Ord
REAL(KIND = idp), INTENT(IN)                        :: x
REAL(KIND = idp), INTENT(IN), DIMENSION(0:Ord)      :: xlocs


INTEGER                                             :: i,j
REAL(KIND = idp), DIMENSION(0:Ord)                  :: tmp
REAL(KIND = idp), DIMENSION(0:Ord)                  :: Lagrange_Poly



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
REAL(KIND = idp), INTENT(IN)                        :: x
REAL(KIND = idp), INTENT(IN), DIMENSION(0:Ord)      :: xlocs



INTEGER                                             :: i,j,k
REAL(KIND = idp), DIMENSION(0:Ord)                  :: tmp,tmp2
REAL(KIND = idp), DIMENSION(0:Ord)                  :: Lagrange_Poly_Deriv


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
REAL(KIND = idp), INTENT(IN)                        :: x
REAL(KIND = idp), INTENT(IN), DIMENSION(0:Ord)      :: xlocs



INTEGER                                             :: i, j, k, l
REAL(KIND = idp), DIMENSION(0:Ord)                  :: tmp,tmp2
REAL(KIND = idp), DIMENSION(0:Ord)                  :: Lagrange_Second_Deriv





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

                        !print*,"tmp A ",i, k, l, tmp(i)

                        DO j = 0,Ord

                            IF (j .NE. i .AND. j .NE. k .AND. j .NE. l) THEN

                                tmp(i) = tmp(i) * (x - xlocs(j))/(xlocs(i) - xlocs(j))
                                !print*,"tmp B ",i, k, l, j, tmp(i)

                            END IF

                        END DO  ! j Loop

                        tmp2(i) = tmp2(i) + tmp(i)

                        !print*,"tmp A",j, i, tmp2(i)

                    END IF

                END DO ! l Loop

            END IF

        END DO  ! k Loop

    END DO  ! i Loop



    Lagrange_Second_Deriv = tmp2


END IF

END FUNCTION Lagrange_Second_Deriv





 !+104+################################################################!
!                                                                       !
!   Legendre_Poly - Calculates the value of the Legendre Polynomial     !
!                   P^m_l(cos(theta))                                   !
!                                                                       !
 !#####################################################################!
PURE FUNCTION Legendre_Poly(l,m,num_points,theta)


!  Returns array that conatins the values P^m_n(cos(theta)) for n = 0:l
!  If m > n then returns 0 as the poly doesn't exist there.

INTEGER, INTENT(IN)                                     :: l,m, num_points
REAL(KIND = idp),  INTENT(IN), DIMENSION(1:num_points)  :: theta
REAL(KIND = idp), DIMENSION(1:num_points)               :: Legendre_Poly, Legendre_Polyc

INTEGER                                         :: i,j,n
REAL(KIND = idp)                                :: factor, normfactor
REAL(KIND = idp), DIMENSION(1:num_points)       :: costheta, sqrfactor
REAL(KIND = idp), DIMENSION(0:l,1:num_points)   :: Plm


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










 !+105+################################################################!
!                                                                       !
!   Spherical_Harmonic - Calculates the value of the spherical harmonic,!
!                       Y^M_L(Theta, Phi). Uses Legendre_Poly           !
!                                                                       !
!           Output - 2L Value array - (Real, Imaginary)                 !
!                                                                       !
 !#####################################################################!
PURE FUNCTION Spherical_Harmonic(l,m,theta,phi)




INTEGER, INTENT(IN)                         :: l,m
REAL(KIND = idp), INTENT(IN)                :: theta, phi


INTEGER                                     :: abs_m
REAL(KIND = idp), DIMENSION(0:0)            :: Plm
COMPLEX(KIND = idp)                         :: Spherical_Harmonic


Plm = Legendre_Poly(l,m,1,[theta])


Spherical_Harmonic = Norm_Factor(l,m)*Plm(0)*EXP(CMPLX(0,m*phi, KIND = idp))


END FUNCTION Spherical_Harmonic








 !+201+####################################################!
!                                                           !
!       Norm_Factor - Calculate normalization factor for    !
!                       spherical harmonics                 !
!                                                           !
 !#########################################################!
PURE ELEMENTAL FUNCTION Norm_Factor(l,m)

INTEGER, INTENT(IN)         :: l,m
REAL(KIND = idp)            :: Norm_Factor


!Norm_Factor = 1.0_idp
Norm_Factor = sqrt( ( ( 2.0_idp * l + 1.0_idp) * Factorial(l - m))/( 4.0_idp * pi * Factorial(l + m)))
!Norm_Factor = sqrt(  (( 2.0_idp * l + 1.0_idp)/( 4.0_idp * pi)) * Factorial_Division(l - m, l + m) )





END FUNCTION Norm_Factor






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










 !+204+################################################################!
!                                                                       !
!       POWER - Simple x**y function for testing code, TEST_LGL         !
!                                                                       !
 !#####################################################################!
PURE ELEMENTAL FUNCTION POWER(x,y)

INTEGER, INTENT(IN)             :: y
REAL(KIND = idp), INTENT(IN)    :: x

INTEGER                         :: i
REAL(KIND = idp)                :: tmp
REAL(KIND = idp)                :: POWER

tmp = 1.0_idp


IF (y < 0) THEN

    DO i = 1,abs(y)

        tmp = tmp/x

    END DO

ELSE

    DO i = 1,y
        tmp = tmp*x
    END DO

END IF

POWER = tmp

END FUNCTION POWER











!+301+##########################################################!
!                                                               !
!      Map_To_X_Space - maps r value between ra, and rb to x    !
!                   space such that x in [-1,1].                !
!                                                               !
!###############################################################!
PURE ELEMENTAL FUNCTION Map_To_X_Space(ra, rb, r)

REAL(KIND = idp)                            ::  Map_To_X_Space
REAL(KIND = idp), intent(in)                ::  ra, rb
REAL(KIND = idp), intent(in)                ::  r

Map_To_X_Space = (2.0_idp*(r - ra))/(rb - ra) - 1.0_idp


END FUNCTION Map_To_X_Space







!+302+##########################################################!
!                                                               !
!      Map_From_X_Space - maps x value between -1, and 1 to r   !
!                   space such that r in [ra,rb].               !
!                                                               !
!###############################################################!
PURE ELEMENTAL FUNCTION Map_From_X_Space(ra, rb, x)

REAL(KIND = idp)                            ::  Map_From_X_Space
REAL(KIND = idp), intent(in)                ::  ra, rb
REAL(KIND = idp), intent(in)                ::  x

Map_From_X_Space = ((rb - ra) * (x+1.0_idp))/2.0_idp  + ra


END FUNCTION Map_From_X_Space










!+401+##########################################################!
!                                                               !
!     MVMULT_CCS: Multiply a vector by a matrix in CCS format   !
!                                                               !
!###############################################################!
FUNCTION MVMULT_CCS(N, NNZ, ELEM_VAL, COL_PTR, ROW_IND, VECT)


INTEGER, INTENT(IN)                                         :: N, NNZ
INTEGER, DIMENSION(0:N), INTENT(IN)                         :: COL_PTR
INTEGER, DIMENSION(0:NNZ - 1),INTENT(IN)                    :: ROW_IND

REAL(KIND = idp), DIMENSION(0:NNZ - 1), INTENT(IN)          :: ELEM_VAL
REAL(KIND = idp), DIMENSION(0:N-1), INTENT(IN)              :: VECT


REAL(KIND = idp), DIMENSION(0:N-1)                          :: MVMULT_CCS

INTEGER                                                     :: i, j

MVMULT_CCS = 0.0_idp

DO i = 0,N-1


    DO j = COL_PTR(i),COL_PTR(i+1)-1

        MVMULT_CCS(ROW_IND(j)) =  MVMULT_CCS(ROW_IND(j)) + ELEM_VAL(j)*VECT(i)

    END DO


END DO



END FUNCTION MVMULT_CCS




!+402+##########################################################!
!                                                               !
!     MVMULT_FULL: Multiply a vector by a matrix.               !
!                                                               !
!###############################################################!
FUNCTION MVMULT_FULL(A, V, N, M)


INTEGER, INTENT(IN)                                         :: N, M
REAL(KIND = idp), INTENT(IN), DIMENSION(1:M)                :: V
REAL(KIND = idp), INTENT(IN), DIMENSION(1:N,1:M)            :: A


REAL(KIND = idp), DIMENSION(1:N)                            :: MVMULT_FULL


REAL(KIND = idp), DIMENSION(1:N)                            :: Sol

INTEGER                                                     :: i,j, io, jo

INTEGER                                                     :: Blocksize

Sol = 0.0_idp


Blocksize = MIN( 100, MAX(M,N) )

!$OMP PARALLEL                              &
!$OMP SHARED( N, M, Sol, A, V)              &
!$OMP PRIVATE( i, j, io, jo )


!$OMP DO SIMD SCHEDULE(Dynamic)
Do i = 1,N, Blocksize

    Do j = 1,M, Blocksize

        DO io = i,MIN(i+Blocksize-1, N)

            DO jo = j,MIN(j+Blocksize-1,M)

                Sol(io) = Sol(io) + A(io,jo)*V(jo)

            END DO

        END DO

    END DO

END DO
!$OMP END DO SIMD

!$OMP END PARALLEL



MVMULT_FULL = Sol

END FUNCTION MVMULT_FULL










!+403+##############################################!
!                                                   !
!     SVVAD - Scaled Vector-Vector Addition         !
!               A = alpha*A + B                     !
!                                                   !
!###################################################!
FUNCTION SVVADD(N, alpha, A, B )

INTEGER, INTENT(IN)                                 :: N
REAL(KIND = idp), INTENT(IN)                        :: alpha

REAL(KIND = idp), DIMENSION(0:N-1), INTENT(IN)      :: B
REAL(KIND = idp), DIMENSION(0:N-1), INTENT(IN)      :: A

REAL(KIND = idp), DIMENSION(0:N-1)                  :: SVVADD

INTEGER                                             :: i

DO i = 0,N-1

    SVVADD(i) = alpha*A(i) + B(i)

END DO


END FUNCTION SVVADD













 !+501+####################################################################!
!                                                                           !
!   Initialize_LGL_Quadrature - Calculate the Legendre-Gauss-Lobatto        !
!               Quadrature node locations and weights.                      !
!                                                                           !
 !#########################################################################!
SUBROUTINE Initialize_LGL_Quadrature(Ord, xloc, weights)


INTEGER, INTENT(IN)                                     ::  Ord
REAL(KIND = idp), INTENT(INOUT), DIMENSION(0:Ord)       ::  xloc, weights


INTEGER                                                 ::  i, j, Op1
REAL(KIND = idp), DIMENSION(0:Ord)                      ::  xtmp, xlast
REAL(KIND = idp), DIMENSION(0:Ord,0:Ord)                ::  VanderM



Op1 = Ord +1

!! Set Initial Guess at x locations !!
DO i = 0,Ord
    xtmp(i) = cos(i*pi/Ord)
END DO


!! Initialize VanderM and xlast !!
VanderM(:,:) = 0.0_idp
xlast(:) = 42

!! Loop until converged on machine precision !!
DO WHILE (MAXVAL(ABS(xtmp(:) - xlast(:))) > eps)

    !! Set xlast to newest versions of xtmp !!
    xlast(:) = xtmp(:)


    VanderM(:,0) = 1
    VanderM(:,1) = xtmp(:)

    !! Calculate Elements of Vandermonde Matrix for Legendre Polynomials !!
    DO i = 1,Ord-1
        VanderM(:,i+1) = ((2*(i+1)-1)*xtmp(:)*VanderM(:,i) - (i)*VanderM(:,i-1))/(i+1)
    END DO

    !! Update xtmp using derivative of Legendre Polynomials !!
    xtmp(:) = xlast(:) - ( xtmp(:)*VanderM(:,Ord) - VanderM(:,Ord-1) )/(Op1 * VanderM(:,Ord))


END DO


DO i = 0,Ord
    xloc(i) = xtmp(Ord-i)
    weights(i) = 2/(Ord*Op1 * VanderM(Ord-i,Ord)*VanderM(Ord-i,Ord))
END DO





END SUBROUTINE Initialize_LGL_Quadrature





 !+502+########################################################################!
!                                                                               !
!   Initialize_LGL_Quadrature_Locations - Calculates the Legendre-Gauss-Lobatto !
!               Quadrature node locations. Does not calculate the weights       !
!               needed for integration, use "Initialize_LGL_Quadrature", +501+, !
!               for locations and weights.                                      !
!                                                                               !
 !#############################################################################!
PURE FUNCTION Initialize_LGL_Quadrature_Locations(Ord)


INTEGER, INTENT(IN)                                     ::  Ord
REAL(KIND = idp),  DIMENSION(0:Ord)                     ::  Initialize_LGL_Quadrature_Locations


INTEGER                                                 ::  i, j, Op1
REAL(KIND = idp), DIMENSION(0:Ord)                      ::  xtmp, xlast
REAL(KIND = idp), DIMENSION(0:Ord,0:Ord)                ::  VanderM



Op1 = Ord +1

!! Set Initial Guess at x locations !!
DO i = 0,Ord
    xtmp(i) = cos(i*pi/Ord)
END DO


!! Initialize VanderM and xlast !!
VanderM(:,:) = 0.0_idp
xlast(:) = 42

!! Loop until converged on machine precision !!
DO WHILE (MAXVAL(ABS(xtmp(:) - xlast(:))) > eps)

    !! Set xlast to newest versions of xtmp !!
    xlast(:) = xtmp(:)


    VanderM(:,0) = 1
    VanderM(:,1) = xtmp(:)

    !! Calculate Elements of Vandermonde Matrix for Legendre Polynomials !!
    DO i = 1,Ord-1
        VanderM(:,i+1) = ((2*(i+1)-1)*xtmp(:)*VanderM(:,i) - (i)*VanderM(:,i-1))/(i+1)
    END DO

    !! Update xtmp using derivative of Legendre Polynomials !!
    xtmp(:) = xlast(:) - ( xtmp(:)*VanderM(:,Ord) - VanderM(:,Ord-1) )/(Op1 * VanderM(:,Ord))


END DO


DO i = 0,Ord
    Initialize_LGL_Quadrature_Locations(i) = xtmp(Ord-i)
END DO





END FUNCTION Initialize_LGL_Quadrature_Locations




!+503+##################################################################!
!                                                                       !
!   Initialize_LG_Quadrature - Calculate the Legendre-Gauss Quadrature  !
!                              node locations and weights.              !
!                                                                       !
!#######################################################################!
SUBROUTINE Initialize_LG_Quadrature(Ord, xloc, weights)

INTEGER, INTENT(IN)                                     ::  Ord
REAL(KIND = idp), INTENT(INOUT), DIMENSION(1:Ord)       ::  xloc, weights

INTEGER                                                 :: i, j, m
REAL(KIND = idp)                                        :: p1, p2, p3, pp, z, z1


m = (Ord + 1)/2

DO i = 1,m

    z = cos(pi * (i-0.25_idp)/(Ord + 0.5_idp))
    z1 = 42.0_idp

    DO WHILE ( ABS(z - z1) .GT. eps)
        p1 = 1.0_idp
        p2 = 0.0_idp

        DO j = 1, Ord

            p3 = p2
            p2 = p1
            p1 = ((2.0_idp*j - 1.0_idp)*z*p2 - (j - 1.0_idp)*p3)/j

        END DO


        pp = Ord*(z*p1 - p2)/(z*z-1.0_idp)
        z1 = z
        z = z1 - p1/pp


    END DO

    xloc(i) = -z
    xloc(Ord-i+1) = +z

    weights(i) = 2.0_idp/((1.0_idp - z*z)*pp*pp)
    weights(Ord-i+1) = weights(i)

END DO

END SUBROUTINE Initialize_LG_Quadrature





 !+504+########################################################################!
!                                                                               !
!   Initialize_LG_Quadrature_Locations - Calculates the Legendre-Gauss          !
!               Quadrature node locations. Does not calculate the weights       !
!               needed for integration, use "Initialize_LG_Quadrature", +503+,  !
!               for locations and weights.                                      !
!                                                                               !
 !#############################################################################!
PURE FUNCTION Initialize_LG_Quadrature_Locations(Ord)

INTEGER, INTENT(IN)                                     ::  Ord


REAL(KIND = idp), DIMENSION(1:Ord)                      ::  Initialize_LG_Quadrature_Locations

INTEGER                                                 :: i, j, m
REAL(KIND = idp)                                        :: p1, p2, p3, pp, z, z1


m = (Ord + 1)/2

DO i = 1,m

    z = cos(pi * (i-0.25_idp)/(Ord + 0.5_idp))
    z1 = 42.0_idp

    DO WHILE ( ABS(z - z1) .GT. eps)
        p1 = 1.0_idp
        p2 = 0.0_idp

        DO j = 1, Ord

            p3 = p2
            p2 = p1
            p1 = ((2.0_idp*j - 1.0_idp)*z*p2 - (j - 1.0_idp)*p3)/j

        END DO


        pp = Ord*(z*p1 - p2)/(z*z-1.0_idp)
        z1 = z
        z = z1 - p1/pp


    END DO

    Initialize_LG_Quadrature_Locations(i) = -z
    Initialize_LG_Quadrature_Locations(Ord-i+1) = +z


END DO

!PRINT*,"!!!!!!!!!!!!!!!"
!PRINT*,xloc
!PRINT*," "
!PRINT*,weights
!PRINT*,"!!!!!!!!!!!!!!!"


END FUNCTION Initialize_LG_Quadrature_Locations









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
                                                * POWER(-1.0_idp, M)*Spherical_Harmonic(L, -M, T_locs(td), P_locs(pd))  & ! Ylm* = (-1)^m Yl-m
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







!+701+##################################################################################!
!                                                                                       !
!       GENERATE_DEFINED_MESH - Generate the values for the mesh sent in using          !
!                                predefined values for the width of each element.       !
!                                                                                       !
!---------------------------------------------------------------------------------------!
!                                                                                       !
!   Input:  Mesh_Start - Single Real number defining inner boundary location.           !
!                                                                                       !
!           Number_of_Elements - Single Integer value defining the number of elements   !
!                                   in the mesh.                                        !
!                                                                                       !
!           Element_Width_Vector - Real valued vector of length(1:Number_of_Elements)   !
!                                       containing Real numbers defining the width of   !
!                                        each element.                                  !
!                                                                                       !
!---------------------------------------------------------------------------------------!
!                                                                                       !
!   Output: Mesh - Real Vector,length(0:Number_of_Elements) that on output contains     !
!                       values describing the element edge locations.                   !
!                                                                                       !
!#######################################################################################!
SUBROUTINE Generate_Defined_Mesh(Number_of_Elements, Mesh_Start, Element_Width_Vector, Mesh)


INTEGER, INTENT(IN)                                                 ::  Number_of_Elements
REAL(KIND = idp), INTENT(IN)                                        ::  Mesh_Start
REAL(KIND = idp), DIMENSION(1:Number_of_Elements), INTENT(IN)       ::  Element_Width_Vector

REAL(KIND = idp), DIMENSION(0:Number_of_Elements), INTENT(OUT)      ::  Mesh

INTEGER                                                             ::  i


mesh(0) = Mesh_Start
DO i = 1,Number_of_Elements


    mesh(i) = mesh(i-1) + Element_Width_Vector(i)


END DO


END SUBROUTINE Generate_Defined_Mesh






!+802+##################################################################################!
!                                                                                       !
!       GENERATE_DEFINED_COARSE_MESH                                                    !
!                                                                                       !
!#######################################################################################!
SUBROUTINE Generate_Defined_Coarse_Mesh(Input_Number_of_Elements, Output_Number_of_Elements,    &
                                        Coarsen_Factor, Mesh_Start, Element_Width_Vector, Mesh)


INTEGER, INTENT(IN)                                                     ::  Input_Number_of_Elements
INTEGER, INTENT(IN)                                                     ::  Output_Number_of_Elements
INTEGER, INTENT(IN)                                                     ::  Coarsen_Factor
REAL(KIND = idp), INTENT(IN)                                            ::  Mesh_Start
REAL(KIND = idp), DIMENSION(1:Input_Number_of_Elements), INTENT(IN)     ::  Element_Width_Vector

REAL(KIND = idp), DIMENSION(0:Output_Number_of_Elements), INTENT(OUT)   ::  Mesh

INTEGER                                                                 ::  i

mesh(0) = Mesh_Start
DO i = 1,Output_Number_of_Elements


    mesh(i) = mesh(i-1) + SUM(Element_Width_Vector((i-1)*Coarsen_Factor+1:i*Coarsen_Factor))
         

END DO


END SUBROUTINE Generate_Defined_Coarse_Mesh













!+701+###########################################################################!
!                                                                                !
!                  Calc_3D_Values_At_Location          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_3D_Values_At_Location( r, theta, phi, Return_Psi, Return_AlphaPsi,  &
                                        Return_Beta1, Return_Beta2, Return_Beta3    )


REAL(KIND = idp), INTENT(IN)                                ::  r, theta, phi
REAL(KIND = idp), INTENT(INOUT)                             ::  Return_Psi,         &
                                                                Return_AlphaPsi,    &
                                                                Return_Beta1,       &
                                                                Return_Beta2,       &
                                                                Return_Beta3



COMPLEX(KIND = idp), DIMENSION(1:5)                         ::  Tmp_U_Value


REAL(KIND = idp)                                                ::  r_tmp
REAL(KIND = idp), DIMENSION(0:DEGREE)                           ::  LagP

INTEGER                                                         ::  re, l, m, d, u


INTEGER                                                         :: Current_Location



COMPLEX(KIND = idp)                                             ::  TMP_VALUE_A
REAL(KIND = idp), DIMENSION(0:DEGREE)                           ::  xlocP, weightP


Tmp_U_Value = 0.0_idp

IF ( r == rlocs(0) ) THEN

    DO u = 1,NUM_CFA_VARS
        DO l = 0,L_Limit
            DO m = -M_VALUES(l),M_VALUES(l)
        


                TMP_VALUE_A = Spherical_Harmonic(l,m,theta,phi)
                Current_Location =  Matrix_Location( u, l, m, 0, 0 )
                Tmp_U_Value(u) = Tmp_U_Value(u) + Coefficient_Vector(Current_Location) * TMP_VALUE_A
                
                
            END DO ! m Loop
        END DO  ! l Loop
    END DO  ! u Loop



ELSE

    CALL Initialize_LGL_Quadrature(DEGREE,xlocP,weightP)

    DO re = 0,NUM_R_ELEMENTS-1

        IF ( r > rlocs(re) .AND. r <= rlocs(re+1) ) THEN

            r_tmp = Map_To_X_Space(rlocs(re),rlocs(re+1),r)
            LagP = Lagrange_Poly(r_tmp,DEGREE,xlocP)

            DO u = 1,NUM_CFA_VARS
                DO l = 0,L_Limit
                    DO m = -M_VALUES(l),M_VALUES(l)
                        DO d = 0,DEGREE


                            TMP_VALUE_A = Spherical_Harmonic(l,m,theta,phi) * LagP(d)
                            Current_Location = Matrix_Location( u, l, m, re, d )
                            Tmp_U_Value(u) = Tmp_U_Value(u) + Coefficient_Vector( Current_Location ) * TMP_VALUE_A

                        END DO  !   d Loop
                    END DO  !   m Loop
                END DO  !   l Loop
            END DO ! u Loop

            EXIT
        END IF

    END DO

    IF ( r > rlocs(NUM_R_ELEMENTS) ) THEN

        DO u = 1,NUM_CFA_VARS
            DO l = 0,L_Limit
                DO m = -M_VALUES(l),M_VALUES(l)


                    Current_Location = Matrix_Location( u, l, m, NUM_R_ELEMENTS-1, DEGREE )
                    TMP_VALUE_A = Spherical_Harmonic(l,m,theta,phi)
                    Tmp_U_Value(u) = Tmp_U_Value(u) + Coefficient_Vector( Current_Location ) * TMP_VALUE_A


                END DO  !   m Loop
            END DO  !   l Loop
        END DO ! u Loop
    END IF

END IF



Return_Psi      = REAL(Tmp_U_Value(1), KIND = idp)
Return_AlphaPsi = REAL(Tmp_U_Value(2), KIND = idp)
Return_Beta1    = REAL(Tmp_U_Value(3), KIND = idp)
Return_Beta2    = REAL(Tmp_U_Value(4), KIND = idp)
Return_Beta3    = REAL(Tmp_U_Value(5), KIND = idp)





END SUBROUTINE Calc_3D_Values_At_Location











END MODULE Poseidon_Additional_Functions_Module
