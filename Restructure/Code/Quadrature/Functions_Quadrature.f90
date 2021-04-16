   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Functions_Quadrature                                                         !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   Initialize_LGL_Quadrature                                           !##!
!##!    +102+   Initialize_LGL_Quadrature_Locations                                 !##!
!##!                                                                                !##!
!##!    +201+   Initialize_LG_Quadrature                                            !##!
!##!    +202+   Initialize_LG_Quadrature_Locations                                  !##!
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
            ONLY :  pi, eps


IMPLICIT NONE


CONTAINS


 !+101+####################################################################!
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





 !+102+########################################################################!
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









!+201+##################################################################!
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










!+201+##################################################################!
!                                                                       !
!   Initialize_Trapezoid_Quadrature - Calculate the Trapezoid Rule      !
!                              quadrature node locations and weights.   !
!                                                                       !
!#######################################################################!
SUBROUTINE Initialize_Trapezoid_Quadrature(Ord, xloc, weights)

INTEGER, INTENT(IN)                                     ::  Ord
REAL(KIND = idp), INTENT(INOUT), DIMENSION(1:Ord)       ::  xloc, weights

INTEGER                                                 ::  j

DO j = 0,Ord-1
    xloc(j+1) = -pi + (2.0_idp*pi/Ord) * j
END DO

weights(:) = 2.0_idp*pi/Ord



END SUBROUTINE Initialize_Trapezoid_Quadrature


END MODULE Functions_Quadrature
