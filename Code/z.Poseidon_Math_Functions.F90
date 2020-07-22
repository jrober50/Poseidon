   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_Math_Functions_Module                                               !##!
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
USE Poseidon_Constants_Module, &
            ONLY :  idp, pi


USE Poseidon_Parameters, &
            ONLY :  DOMAIN_DIM,             &
                    DEGREE,                 &
                    L_LIMIT,                &
                    NUM_BLOCK_THETA_ROWS,   &
                    NUM_T_ELEMS_PER_BLOCK,  &
                    NUM_P_ELEMS_PER_BLOCK,  &
                    NUM_R_QUAD_POINTS,      &
                    NUM_T_QUAD_POINTS,      &
                    NUM_P_QUAD_POINTS

USE Poseidon_Variables_Module, &
            ONLY :  NUM_R_ELEMENTS,         &
                    NUM_T_ELEMENTS,         &
                    NUM_P_ELEMENTS,         &
                    NUM_TP_QUAD_POINTS,     &
                    rlocs, tlocs, plocs,    &
                    NUM_R_NODES,            &
                    INT_R_LOCATIONS,        &
                    INT_T_LOCATIONS,        &
                    INT_P_LOCATIONS,        &
                    VAR_DIM,                &
                    PROB_DIM,               &
                    ULM_LENGTH,             &
                    Ylm_Table_Block,        &
                    Ylm_Values,             &
                    Ylm_dt_Values,          &
                    Ylm_dp_Values,          &
                    Ylm_CC_Values,          &
                    Ylm_CC_dt_Values,       &
                    Ylm_CC_dp_Values,       &
                    LM_LENGTH,              &
                    myID_Shell,             &
                    Lagrange_Poly_Table,    &
                    LPT_LPT,                &
                    Coefficient_Vector,     &
                    LM_Location,            &
                    Matrix_Location,        &
                    M_VALUES,               &
                    R_Inner,                &
                    R_Outer

USE Poseidon_Quadrature_Module, &
            ONLY :  Initialize_LGL_Quadrature_Locations

USE Poseidon_Mapping_Functions_Module, &
            ONLY :  Map_From_X_Space


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





 !+201+################################################################!
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






 !+202+####################################################!
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


END FUNCTION Norm_Factor




 !+203+################################################################!
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
REAL(KIND = idp), INTENT(IN)                :: theta, phi


INTEGER                                     :: abs_m
REAL(KIND = idp), DIMENSION(0:0)            :: Plm
COMPLEX(KIND = idp)                         :: Spherical_Harmonic


Plm = Legendre_Poly(l,m,1,[theta])


Spherical_Harmonic = Norm_Factor(l,m)*Plm(0)*EXP(CMPLX(0,m*phi, KIND = idp))


END FUNCTION Spherical_Harmonic




END MODULE Poseidon_Math_Functions_Module
