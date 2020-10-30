   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Variables_Derived                                                            !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains parameters used to define the running of Poseidon.                 !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!

USE Poseidon_Kinds_Module, &
                ONLY : idp


IMPLICIT NONE

INTEGER, PUBLIC                             :: NUM_R_NODES
INTEGER, PUBLIC                             :: NUM_R_NODESp1
INTEGER, PUBLIC                             :: BLOCK_NUM_R_NODES
INTEGER, PUBLIC                             :: SUBSHELL_NUM_R_NODES

INTEGER, PUBLIC                             :: VAR_DIM
INTEGER, PUBLIC                             :: ELEM_VAR_DIM
INTEGER, PUBLIC                             :: BLOCK_VAR_DIM
INTEGER, PUBLIC                             :: SUBSHELL_VAR_DIM

INTEGER, PUBLIC                             :: PROB_DIM
INTEGER, PUBLIC                             :: ELEM_PROB_DIM
INTEGER, PUBLIC                             :: ELEM_PROB_DIM_SQR
INTEGER, PUBLIC                             :: BLOCK_PROB_DIM
INTEGER, PUBLIC                             :: SUBSHELL_PROB_DIM

INTEGER, PUBLIC                             :: LM_LENGTH
INTEGER, PUBLIC                             :: ULM_LENGTH

INTEGER, PUBLIC                             :: NUM_OFF_DIAGONALS


END MODULE Variables_Derived

