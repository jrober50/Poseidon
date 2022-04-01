   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Driver_SetGuess_Module                                                !##!
!##!                                                                         !##!
!##!_________________________________________________________________________!##!
!##!                                                                         !##!
!##!                                                                         !##!
!##!=========================================================================!##!
!##!                                                                         !##!
!##!    Contains:                                                            !##!
!##!                                                                         !##!
!##!                                                                         !##!
!###############################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !##########################################################################!


!*D*================================!
!                                   !
!           Dependencies            !
!                                   !
!===================================!


USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Poseidon_Units_Module, &
            ONLY :  C_Square

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag

USE Variables_Functions, &
            ONLY :  Potential_Solution

USE Initial_Guess_Module, &
            ONLY :  Poseidon_Input_Guess


IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!     Driver_SetGuess                                                			!
!                                                                               !
!###############################################################################!
SUBROUTINE Driver_SetGuess( NE, NQ,                     &
                            dx_c, x_e,                  &
                            R_Quad, T_Quad, P_Quad,     &
                            LeftLimit, RightLimit,      &
                            Guess_Type                  )

INTEGER, INTENT(IN), DIMENSION(3)                       ::  NE
INTEGER, INTENT(IN), DIMENSION(3)                       ::  NQ

REAL(idp), INTENT(IN), DIMENSION(1:NE(1))               ::  dx_c
REAL(idp), INTENT(IN), DIMENSION(0:NE(1))               ::  x_e

REAL(idp), INTENT(IN), DIMENSION(1:NQ(1))               ::  R_Quad
REAL(idp), INTENT(IN), DIMENSION(1:NQ(2))               ::  T_Quad
REAL(idp), INTENT(IN), DIMENSION(1:NQ(3))               ::  P_Quad

REAL(idp), INTENT(IN)                                   ::  LeftLimit
REAL(idp), INTENT(IN)                                   ::  RightLimit

INTEGER,   INTENT(IN)                                   ::  Guess_Type


INTEGER                                                 ::  Num_DOF

REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              ::  Psi_Guess
REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              ::  AlphaPsi_Guess
REAL(idp), DIMENSION(:,:,:,:,:), ALLOCATABLE            ::  Beta_Guess




IF ( Verbose_Flag ) THEN
    WRITE(*,'(A)')"In Driver, Setting Guess."
END IF

Num_DOF = NQ(1)*NQ(2)*NQ(3)

ALLOCATE( Psi_Guess(1:Num_DOF,0:NE(1)-1,0:NE(2)-1, 0:NE(3)-1 )  )
ALLOCATE( AlphaPsi_Guess(1:Num_DOF,0:NE(1)-1,0:NE(2)-1, 0:NE(3)-1 )  )
ALLOCATE( Beta_Guess(1:Num_DOF,0:NE(1)-1,0:NE(2)-1, 0:NE(3)-1, 3 )  )

Psi_Guess = 1.0_idp
AlphaPsi_Guess = 1.0_idp
Beta_Guess = 0.0_idp



CALL Poseidon_Input_Guess(  Psi_Guess,          &
                            AlphaPsi_Guess,     &
                            Beta_Guess,         &
                            NE(1),              &
                            NE(2),              &
                            NE(3),              &
                            NQ(1),              &
                            NQ(2),              &
                            NQ(3),              &
                            R_Quad,             &
                            T_Quad,             &
                            P_Quad,             &
                            LeftLimit,          &
                            RightLimit          )



DEALLOCATE( Psi_Guess )
DEALLOCATE( AlphaPsi_Guess )
DEALLOCATE( Beta_Guess )

END SUBROUTINE Driver_SetGuess



END MODULE Driver_SetGuess_Module
