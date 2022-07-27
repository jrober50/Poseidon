   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Flags_Run_Check_Module                                                       !##!
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


! Allocation Flags
INTEGER,    PUBLIC, PARAMETER       ::  iPF_RC_Num_Flags    = 3

INTEGER,    PUBLIC, PARAMETER       ::  iPF_RC_Ready        = 1
INTEGER,    PUBLIC, PARAMETER       ::  iPF_RC_Init         = 2
INTEGER,    PUBLIC, PARAMETER       ::  iPF_RC_Source       = 3
INTEGER,    PUBLIC, PARAMETER       ::  iPF_RC_BC           = 4
INTEGER,    PUBLIC, PARAMETER       ::  iPF_RC_Guess        = 5


LOGICAL,    PUBLIC, DIMENSION(1:iPF_RC_Num_Flags)    ::  lPF_RC_Flags






CONTAINS




END MODULE Flags_Run_Check_Module




