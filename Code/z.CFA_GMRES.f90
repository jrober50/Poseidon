   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE CFA_GMRES_Module                                                             !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!

USE MPI

USE OMP_LIB


USE Poseidon_Constants_Module, &
                ONLY :  idp

USE Poseidon_Variables_Module, &
                ONLY :  Prob_Dim,               &
                        Coefficient_Vector

USE GMRES_Functions_Module, &
                ONLY :  CDOT_PRODUCT,   &
                        CNORM2,         &
                        ROTMAT

USE GMRES_MultJV_Module, &
                ONLY : GMRES_MultJV





IMPLICIT NONE


INTEGER, PARAMETER                              :: Iter_Max     = 10
INTEGER, PARAMETER                              :: Iter_Restart = 5
REAL(KIND = idp), PARAMETER                     :: Tolerance    = 1.0E-15


CONTAINS




  !+101+############################################################################!
 !                                                                                   !
 !                       CFA_GMRES                                                   !
 !                                                                                   !
  !#################################################################################!
 SUBROUTINE CFA_GMRES( Eq_Flag_Array )

 INTEGER, DIMENSION(1:5), INTENT(IN)                ::  Eq_Flag_Array



INTEGER                                             :: Converge_Flag

REAL(KIND = idp)                                    :: Error
REAL(KIND = idp)                                    :: bnorm2


INTEGER                                             :: Iter_Cur
INTEGER                                             :: i, k

INTEGER, DIMENSION(:), ALLOCATABLE                  :: IPIV
INTEGER                                             :: INFO

COMPLEX(KIND = idp), DIMENSION(:), ALLOCATABLE      :: Residual
COMPLEX(KIND = idp), DIMENSION(:), ALLOCATABLE      :: RHS_Vector
COMPLEX(KIND = idp)                                 :: tmp


COMPLEX(KIND = idp), DIMENSION(:,:), ALLOCATABLE    :: V
COMPLEX(KIND = idp), DIMENSION(:,:), ALLOCATABLE    :: H, H_tmp
COMPLEX(KIND = idp), DIMENSION(:), ALLOCATABLE      :: cs
COMPLEX(KIND = idp), DIMENSION(:), ALLOCATABLE      :: sn
COMPLEX(KIND = idp), DIMENSION(:), ALLOCATABLE      :: cs_tmp
COMPLEX(KIND = idp), DIMENSION(:), ALLOCATABLE      :: y
COMPLEX(KIND = idp), DIMENSION(:), ALLOCATABLE      :: x, xb
COMPLEX(KIND = idp), DIMENSION(:), ALLOCATABLE      :: e1
COMPLEX(KIND = idp), DIMENSION(:), ALLOCATABLE      :: s
COMPLEX(KIND = idp), DIMENSION(:), ALLOCATABLE      :: s_tmp

COMPLEX(KIND = idp), DIMENSION(:), ALLOCATABLE      :: Alt

COMPLEX(KIND = idp), DIMENSION(:), ALLOCATABLE      :: w
                       


ALLOCATE( V(1:Prob_Dim,1:Iter_Restart+1) )
ALLOCATE( H(1:Iter_Restart+1,1:Iter_Restart+1))
ALLOCATE( H_tmp(1:Iter_Restart+1,1:Iter_Restart+1))
ALLOCATE( IPIV(1:Iter_Restart+1))
ALLOCATE( cs(1:Iter_Restart) )
ALLOCATE( sn(1:Iter_Restart) )
ALLOCATE( cs_tmp(1:2) )
ALLOCATE( y(1:Iter_Restart) )
ALLOCATE( e1(1:Prob_Dim) )
ALLOCATE( s(1:Iter_Restart+1) )
ALLOCATE( s_tmp(1:Iter_Restart+1) )
ALLOCATE( w(1:Prob_Dim) )
ALLOCATE( x(1:Prob_Dim) )   ! Coefficient Vector
ALLOCATE( xb(1:Prob_Dim) )  ! Coefficient Vector

ALLOCATE( Alt(1:Prob_Dim) )


ALLOCATE( Residual(1:Prob_Dim))
ALLOCATE( RHS_Vector(1:Prob_Dim) )

e1 = 0.0_idp
e1(1) = 1.0_idp
cs = 0.0_idp
sn = 0.0_idp


!Residual = Equation(x, Eq_Flag)
bnorm2 = cnorm2( Residual)
IF (bnorm2 == 0.0_idp) THEN
    bnorm2 = 1.0_idp
END IF
Error  = cnorm2(Residual)/bnorm2


DO Iter_Cur = 1,Iter_Max

    ! Calculate Residual !
    ! r = b-Ax
    CALL GMRES_MultJV( Residual, V(:,i), Coefficient_Vector, Prob_Dim)
!    Residual = Equation(x, Eq_Flag) - Residual


    V(:,1) = Residual/cnorm2(Residual)
    s(1:Prob_Dim) = cnorm2(Residual)*e1(1:Prob_Dim)


    DO i = 1,Iter_Restart       ! Construct orthonomral basis via Gram-Schmidt

        ! Calculate w
        ! w = A*v
        CALL GMRES_MultJV( w, V(:,i), Coefficient_Vector, Prob_Dim)



        DO k = 1,i
            H(k,i) = DOT_PRODUCT(V(:,k),w(:))
            w = w - H(k,i)*V(:,k)
        END DO ! k Loop


        H(i+1,i) = cnorm2( w )
        V(:,i+1) = w/H(i+1,i)


        CALL Givens_Rotation( i, Iter_Restart, cs, sn, H )

        ! form i-th rotation matrix
        ! cs and sn
        cs_tmp = Rotmat( H(i,i), H(i+1,i) )
        cs(i) = cs_tmp(1)
        sn(i) = cs_tmp(2)
        tmp = cs(i)*s(i)
        s(i+1) = -DCONJG(sn(i))*s(i)
        s(i) = tmp



        H(i,i) = cs(i)*H(i,i) + sn(i)*H(i+1,i)
        H(i+1,i) = 0.0_idp
        Error = abs( s(i+1) )/ bnorm2


        IF ( Error <= Tolerance) THEN

!           y = H(1:i,1:i) / s(1:i)  replaced by the following LAPACK call
            H_tmp(1:i,1:i) = H(1:i,1:i)
            s_tmp(1:i) = s(1:i)

            CALL ZGESV( i, 1, H_tmp, Iter_Restart+1, IPIV, s_tmp, Iter_Restart+1, INFO )
            y(1:i) = s_tmp(1:i)


!           x = x + V(:,1:i)*y(1:i)
            CALL ZGEMV( 'N', Prob_Dim, i, 1.0_idp, V(:,1:i), Prob_Dim, y, 1, 1.0_idp, Coefficient_Vector, 1 )

            EXIT
        END IF

    END DO ! i Loop

    IF ( Error <= Tolerance ) THEN
        EXIT
    END IF


!   y = H(1:Iter_Restart,1:Iter_Restart) / s(1:Iter_Restart)  replaced by the following LAPACK call
    H_tmp(1:Iter_Restart,1:Iter_Restart) = H(1:Iter_Restart,1:Iter_Restart)
    s_tmp(1:Iter_Restart) = s(1:Iter_Restart)
    CALL ZGESV( Iter_Restart, 1, H_tmp, Iter_Restart+1, IPIV, s_tmp, Iter_Restart+1, INFO )
    y = s_tmp


!    PRINT*,"Before ZGEMV 2"
!    x = x + V(:,1:Iter_Restart)*y
    CALL ZGEMV( 'N', Prob_Dim, Iter_Restart, 1.0_idp,   &
                    V(:,1:Iter_Restart), Prob_Dim, y,   &
                    1, 1.0_idp, Coefficient_Vector, 1 )


!    PRINT*,"Before MultJV 2"
    ! Calculate Residual !
    ! r = b-Ax
    CALL GMRES_MultJV( Residual, V(:,Iter_Restart), Coefficient_Vector, Prob_Dim)
!    Residual = Residual - Equation(x, Eq_Flag)


    s(i+1) = cnorm2(Residual)
    Error = s(i+1)/bnorm2
    IF ( Error <= Tolerance ) THEN
            EXIT
    END IF

END DO ! Iter_Cur Loop



END SUBROUTINE CFA_GMRES





  !+201+############################################################################!
 !                                                                                   !
 !                       Givens_Rotation                                             !
 !                                                                                   !
  !#################################################################################!
SUBROUTINE Givens_Rotation(Iter, Size, cs, sn, H )

INTEGER, INTENT(IN)                                                 :: Iter
INTEGER, INTENT(IN)                                                 :: Size
COMPLEX( KIND = idp ), DIMENSION(1:Size), INTENT(IN)                :: cs
COMPLEX( KIND = idp ), DIMENSION(1:Size), INTENT(IN)                :: sn
COMPLEX( KIND = idp ), DIMENSION(1:Size+1,1:Size+1), INTENT(INOUT)  :: H


INTEGER                                                             :: k
COMPLEX( KIND = idp )                                               :: tmp

DO k = 1,Iter-1        ! Apply Givens rotation

    tmp = cs(k)*H(k,Iter) + sn(k)*H(k+1,Iter)
    H(k+1,Iter) = -DCONJG(sn(k))*H(k,Iter) + cs(k)*H(k+1,Iter)
    H(k,Iter) = tmp

END DO ! k Loop



END SUBROUTINE Givens_Rotation




END MODULE CFA_GMRES_Module
