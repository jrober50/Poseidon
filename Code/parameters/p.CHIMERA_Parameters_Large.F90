   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE CHIMERA_Parameters                                                           !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains the top level subroutines needed to inialize, run, and close       !##!
!##!        Poseidon.                                                               !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!


USE Poseidon_Constants_Module, &
            ONLY :  idp, pi, fdp


IMPLICIT NONE



INTEGER, PARAMETER                     ::   CHIMERA_R_ELEMS = 256
INTEGER, PARAMETER                     ::   CHIMERA_T_ELEMS = 32
INTEGER, PARAMETER                     ::   CHIMERA_P_ELEMS = 16

INTEGER, PARAMETER                     ::   CHIMERA_R_INPUT_NODES = 1
INTEGER, PARAMETER                     ::   CHIMERA_T_INPUT_NODES = 1
INTEGER, PARAMETER                     ::   CHIMERA_P_INPUT_NODES = 1

INTEGER, PARAMETER                     ::   CHIMERA_LEFT_LIMIT = -0.5_idp
INTEGER, PARAMETER                     ::   CHIMERA_RIGHT_LIMIT = 0.5_idp

INTEGER, PARAMETER                     ::   CHIMERA_DIMENSION = 3

INTEGER, PARAMETER                     ::   CHIMERA_PROCS = 512
INTEGER, PARAMETER                     ::   CHIMERA_y_PROCS = 32
INTEGER, PARAMETER                     ::   CHIMERA_z_PROCS = 16






END MODULE CHIMERA_Parameters
