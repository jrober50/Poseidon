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

USE Variables_Mesh, &
        ONLY :  NUM_R_ELEMENTS,             &
                NUM_T_ELEMENTS,             &
                NUM_P_ELEMENTS

USE Variables_Quadrature, &
        ONLY :  NUM_R_QUAD_POINTS,          &
                NUM_T_QUAD_POINTS,          &
                NUM_P_QUAD_POINTS,         &
                INT_R_LOCATIONS,            &
                INT_T_LOCATIONS,            &
                INT_P_LOCATIONS


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


ELSE IF ( Method_Flag >= 2 ) THEN

    CALL FP_Input_Guess( Psi_Guess,                                  &
                         AlphaPsi_Guess,                             &
                         Beta_Guess,                                 &
                         Input_RE, Input_TE, Input_PE,               &
                         Input_RQ, Input_TQ, Input_PQ,               &
                         Input_R_Quad, Input_T_Quad, Input_P_Quad,   &
                         Left_Limit, Right_Limit                     )

END IF





END SUBROUTINE Poseidon_Input_Guess










!+101+###########################################################################!
!                                                                                !
!               Poseidon_Input_Guess                                             !
!                                                                                !
!################################################################################!
SUBROUTINE Poseidon_Init_FlatGuess()


REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              ::  Psi_Guess
REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              ::  AlphaPsi_Guess
REAL(idp), DIMENSION(:,:,:,:,:), ALLOCATABLE            ::  Beta_Guess


REAL(idp)                                                               :: Left_Limit
REAL(idp)                                                               :: Right_Limit
INTEGER                                                                 :: Num_DOF



Num_DOF =NUM_R_QUAD_POINTS*NUM_T_QUAD_POINTS*NUM_P_QUAD_POINTS

ALLOCATE( Psi_Guess(1:Num_DOF, 0:Num_R_Elements-1, 0:Num_T_Elements-1, 0:Num_P_Elements-1 ) )
ALLOCATE( AlphaPsi_Guess(1:Num_DOF, 0:Num_R_Elements-1, 0:Num_T_Elements-1, 0:Num_P_Elements-1 ) )
ALLOCATE( Beta_Guess(1:Num_DOF, 0:Num_R_Elements-1, 0:Num_T_Elements-1, 0:Num_P_Elements-1,1:3 ) )


Left_Limit  = -0.50_idp
Right_Limit = +0.50_idp


Psi_Guess = 1.0_idp
AlphaPsi_Guess = 1.0_idp
Beta_Guess = 0.0_idp


IF ( Method_Flag == 1 ) THEN

    CALL NR_Input_Guess( Psi_Guess,                                  &
                         AlphaPsi_Guess,                             &
                         Beta_Guess,                                 &
                         Num_R_Elements, Num_R_Elements, Num_R_Elements,               &
                         NUM_R_QUAD_POINTS, NUM_T_QUAD_POINTS, NUM_P_QUAD_POINTS,               &
                         INT_R_LOCATIONS, INT_T_LOCATIONS, INT_P_LOCATIONS,   &
                         Left_Limit, Right_Limit                     )


ELSE IF ( Method_Flag >= 2 ) THEN

    CALL FP_Input_Guess( Psi_Guess,                                  &
                         AlphaPsi_Guess,                             &
                         Beta_Guess,                                 &
                         Num_R_Elements, Num_R_Elements, Num_R_Elements,               &
                         NUM_R_QUAD_POINTS, NUM_T_QUAD_POINTS, NUM_P_QUAD_POINTS,               &
                         INT_R_LOCATIONS, INT_T_LOCATIONS, INT_P_LOCATIONS,   &
                         Left_Limit, Right_Limit                     )

END IF




END SUBROUTINE Poseidon_Init_FlatGuess













END MODULE Initial_Guess_Module
