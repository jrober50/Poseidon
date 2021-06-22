   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Initial_Guess_Module                                                      !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
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

USE Poseidon_Parameters, &
            ONLY :  Method_Flag

USE FP_Initial_Guess_Module, &
            ONLY : FP_Input_Guess

USE NR_Initial_Guess_Module, &
            ONLY : NR_Input_Guess

IMPLICIT NONE

CONTAINS



!+101+###########################################################################!
!                                                                                !
!               Poseidon_Input_Guess                                             !
!                                                                                !
!################################################################################!
SUBROUTINE Poseidon_Input_Guess(Psi_Guess,                                  &
                                AlphaPsi_Guess,                             &
                                Beta_Guess,                                 &
                                Input_RE, Input_TE, Input_PE,               &
                                Input_RQ, Input_TQ, Input_PQ,               &
                                Input_R_Quad, Input_T_Quad, Input_P_Quad,   &
                                Left_Limit, Right_Limit                     )


REAL(idp), DIMENSION(1:Input_RQ*Input_TQ*Input_PQ,1:Input_RE,1:Input_TE,1:Input_PE),        &
                                                                        INTENT(IN)  :: Psi_Guess

REAL(idp), DIMENSION(1:Input_RQ*Input_TQ*Input_PQ,1:Input_RE,1:Input_TE,1:Input_PE),        &
INTENT(IN)  :: AlphaPsi_Guess

REAL(idp), DIMENSION(1:Input_RQ*Input_TQ*Input_PQ,1:Input_RE,1:Input_TE,1:Input_PE,1:3),        &
INTENT(IN)  :: Beta_Guess

INTEGER, INTENT(IN)                                                                 :: Input_RE
INTEGER, INTENT(IN)                                                                 :: Input_TE
INTEGER, INTENT(IN)                                                                 :: Input_PE

INTEGER, INTENT(IN)                                                                 :: Input_RQ
INTEGER, INTENT(IN)                                                                 :: Input_TQ
INTEGER, INTENT(IN)                                                                 :: Input_PQ

REAL(idp), DIMENSION(1:Input_RQ), INTENT(IN)                                        :: Input_R_Quad
REAL(idp), DIMENSION(1:Input_TQ), INTENT(IN)                                        :: Input_T_Quad
REAL(idp), DIMENSION(1:Input_PQ), INTENT(IN)                                        :: Input_P_Quad

REAL(idp), INTENT(IN)                                                               :: Left_Limit
REAL(idp), INTENT(IN)                                                               :: Right_Limit




IF ( Method_Flag == 1 ) THEN

    CALL NR_Input_Guess( Psi_Guess,                                  &
                         AlphaPsi_Guess,                             &
                         Beta_Guess,                                 &
                         Input_RE, Input_TE, Input_PE,               &
                         Input_RQ, Input_TQ, Input_PQ,               &
                         Input_R_Quad, Input_T_Quad, Input_P_Quad,   &
                         Left_Limit, Right_Limit                     )


ELSE IF ( Method_Flag == 2 ) THEN

    CALL FP_Input_Guess( Psi_Guess,                                  &
                         AlphaPsi_Guess,                             &
                         Beta_Guess,                                 &
                         Input_RE, Input_TE, Input_PE,               &
                         Input_RQ, Input_TQ, Input_PQ,               &
                         Input_R_Quad, Input_T_Quad, Input_P_Quad,   &
                         Left_Limit, Right_Limit                     )

END IF





END SUBROUTINE Poseidon_Input_Guess



















END MODULE Initial_Guess_Module
