   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE GMRES_Functions_Module                                                       !##!
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

USE Poseidon_Constants_Module, &
                ONLY :  idp


IMPLICIT NONE

CONTAINS



!+101+############################################################################!
!                                                                                   !
!                       CDOT_PRODUCT                                                !
!                                                                                   !
 !#################################################################################!
FUNCTION CDOT_PRODUCT(x, y )

COMPLEX(KIND = idp)                                                     :: CDOT_PRODUCT
COMPLEX(KIND = idp), DIMENSION(:),INTENT(IN)                            :: x
COMPLEX(KIND = idp), DIMENSION(:),INTENT(IN)                            :: y


INTEGER                                                     :: Length_x, Length_y,i

Length_x = size(x)
Length_y = size(y)

If (Length_x .NE. Length_y ) THEN
    PRINT*,"Error in CDOT_PRODUCT, Vectors are different lengths. "
    PRINT*,"Length Vector x, ",Length_x
    PRINT*,"Length Vector y, ",Length_y
ELSE
    CDOT_PRODUCT = 0.0_idp
    DO i = 1,Length_x
        CDOT_PRODUCT = x(i)*DCONJG(y(i))
    END DO

END IF

END FUNCTION CDOT_PRODUCT









 !+102+############################################################################!
!                                                                                   !
!                       CNORM2                                                      !
!                                                                                   !
 !#################################################################################!
FUNCTION CNORM2( Vec )

COMPLEX(Kind = idp), DIMENSION(:), INTENT(IN)       :: Vec
REAL(KIND = idp )                                   :: CNORM2

CNORM2 = NORM2( (/NORM2(REAL(vec, KIND = idp)), NORM2(AIMAG(vec))/) )


END FUNCTION CNORM2






 !+103+############################################################################!
!                                                                                   !
!                       ROTMAT                                                      !
!                                                                                   !
 !#################################################################################!
FUNCTION ROTMAT( a, b )

COMPLEX(KIND = idp), DIMENSION(1:2)         :: Rotmat
COMPLEX(KIND = idp), INTENT(IN)             :: a, b

Rotmat(1) = a/sqrt(zabs(a)*zabs(a) + b*b)
Rotmat(2) = b/sqrt(zabs(a)*zabs(a) + b*b)

END FUNCTION Rotmat









END MODULE GMRES_Functions_Module

