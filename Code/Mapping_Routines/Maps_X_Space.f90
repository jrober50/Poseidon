   !#################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Maps_X_Space                                                                 !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   Map_To_X_Space                                                      !##!
!##!    +102+   Map_From_X_Space                                                    !##!
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
            ONLY : idp

USE Variables_Quadrature, &
            ONLY :  xLeftLimit,     &
                    xRightLimit

IMPLICIT NONE


CONTAINS




!+101+##########################################################!
!                                                               !
!      Map_To_X_Space - maps r value between ra, and rb to x    !
!                   space such that x in [-1,1].                !
!                                                               !
!###############################################################!
PURE ELEMENTAL FUNCTION Map_To_X_Space( rLeftLimit,     &
                                        rRightLimit,    &
                                        r               ) RESULT( x )

REAL(idp), INTENT(IN)                   ::  rLeftLimit
REAL(idp), INTENT(IN)                   ::  rRightLimit
REAL(idp), INTENT(IN)                   ::  r

REAL(idp)                               ::  x

x = Map_Between_Spaces( r,                  &
                        rLeftLimit,         &
                        rRightLimit,        &
                        xLeftLimit,         &
                        xRightLimit         )


END FUNCTION Map_To_X_Space







!+102+##########################################################!
!                                                               !
!      Map_From_X_Space - maps x value between -1, and 1 to r   !
!                   space such that r in [ra,rb].               !
!                                                               !
!###############################################################!
PURE ELEMENTAL FUNCTION Map_From_X_Space( rLeftLimit,   &
                                          rRightLimit,  &
                                          x             ) RESULT( r )

REAL(idp), INTENT(IN)                   ::  rLeftLimit
REAL(idp), INTENT(IN)                   ::  rRightLimit
REAL(idp), INTENT(IN)                   ::  x

REAL(idp)                               ::  r


r = Map_Between_Spaces( x,                  &
                        xLeftLimit,         &
                        xRightLimit,        &
                        rLeftLimit,         &
                        rRightLimit         )


END FUNCTION Map_From_X_Space







!+201+##########################################################!
!                                                               !
!      Map_To_General_X_Space - maps r value between sa, and sb !
!                               to x-space such that x is in the!
!                               range [da,db].                  !
!                                                               !
!###############################################################!
PURE ELEMENTAL FUNCTION Map_Between_Spaces( rs,         &
                                            sLeft,      &
                                            sRight,     &
                                            dLeft,      &
                                            dRight      ) RESULT( rd )

REAL(idp), INTENT(IN)                   ::  rs
REAL(idp), INTENT(IN)                   ::  sLeft
REAL(idp), INTENT(IN)                   ::  sRight
REAL(idp), INTENT(IN)                   ::  dLeft
REAL(idp), INTENT(IN)                   ::  dRight

REAL(idp)                               ::  rd

rd = ((dRight-dLeft)/(sRight-sLeft))*(rs-sLeft) + dLeft

END FUNCTION Map_Between_Spaces









END MODULE Maps_X_Space

