   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Variables_MPI                                                             !##!
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

INTEGER, PUBLIC                                             ::  POSEIDON_COMM_WORLD
INTEGER, PUBLIC                                             ::  POSEIDON_COMM_SHELL
INTEGER, PUBLIC                                             ::  POSEIDON_COMM_DIST
INTEGER, PUBLIC                                             ::  POSEIDON_COMM_PETSC

INTEGER, PUBLIC                                             ::  nPROCS_POSEIDON
INTEGER, PUBLIC                                             ::  nPROCS_PETSC
INTEGER, PUBLIC                                             ::  nPROCS_SHELL

INTEGER, PUBLIC                                             ::  ierr,                   &
                                                                NUM_BLOCKS_THETA,       &
                                                                NUM_BLOCKS_PHI


INTEGER, PUBLIC                                             ::  myShell = -1

INTEGER, PUBLIC                                             ::  myID_Poseidon   =  0
INTEGER, PUBLIC                                             ::  myID_Shell      = -1
INTEGER, PUBLIC                                             ::  myID_SubShell   = -1
INTEGER, PUBLIC                                             ::  myID_Dist       = -1
INTEGER, PUBLIC                                             ::  myID_PETSc      = -1

INTEGER, PUBLIC                                             ::  Local_Length = -1

INTEGER, PUBLIC, PARAMETER                                  ::  MasterID_Poseidon = 0

INTEGER, PUBLIC                      ::  NUM_SHELLS                 = 1
INTEGER, PUBLIC                      ::  NUM_SUBSHELLS_PER_SHELL    = 1
INTEGER, PUBLIC                      ::  NUM_BLOCKS                 = 1
INTEGER, PUBLIC                      ::  NUM_SUBSHELLS              = 1

INTEGER, PUBLIC                      ::  NUM_BLOCKS_PER_SHELL       = 1
INTEGER, PUBLIC                      ::  NUM_BLOCK_THETA_ROWS       = 1
INTEGER, PUBLIC                      ::  NUM_BLOCK_PHI_COLUMNS      = 1

INTEGER, PUBLIC                      ::  NUM_R_ELEMS_PER_SHELL     ! NUM_SHELS*RAD_ELEMS_PER_SHELL = CHIMERA_R_ELEM
INTEGER, PUBLIC                      ::  NUM_R_ELEMS_PER_SUBSHELL   ! NUM_R_ELEMS_PER_SUBSHELL = NUM_R_ELEMS_PER_BLOCK/NUM_SUBSHELLS_PER_SHELL

INTEGER, PUBLIC                      ::  NUM_R_ELEMS_PER_BLOCK
INTEGER, PUBLIC                      ::  NUM_T_ELEMS_PER_BLOCK     ! NUM_BLOCK_THETA_ROW*NUM_T_ELEMS_PER_BLOCK = CHIMERA_T_ELEM
INTEGER, PUBLIC                      ::  NUM_P_ELEMS_PER_BLOCK     ! NUM_BLOCK_PHI_COLUMNS*NUM_P_ELEMS_PER_BLOCK = CHIMERA_P_ELEM


INTEGER                              ::  Ratio_BNDLperBLCK,      &
                                         Ratio_T_BNDLperBLCK,    &
                                         Ratio_P_BNDLperBLCK



END MODULE Variables_MPI


