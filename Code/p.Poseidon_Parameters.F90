   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_Parameters                                                          !##!
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


USE Poseidon_Constants_Module, &
            ONLY :  idp, pi, fdp


IMPLICIT NONE




INTEGER                     ::  DOMAIN_DIM
INTEGER                     ::  DEGREE
INTEGER                     ::  L_LIMIT

INTEGER                     ::  nPROCS_POSEIDON

INTEGER, PARAMETER          ::  NUM_CFA_VARS = 5

INTEGER, PARAMETER          ::  DATA_DIST_MODE = 4
INTEGER, PARAMETER          ::  SOL_DIST_SCHEME = 2
INTEGER, PARAMETER          ::  STF_MAPPING_FLAG = 2

INTEGER                     ::  NUM_SHELLS
INTEGER                     ::  NUM_SUBSHELLS_PER_SHELL


INTEGER                     ::  NUM_BLOCKS_PER_SHELL
INTEGER                     ::  NUM_BLOCK_THETA_ROWS
INTEGER                     ::  NUM_BLOCK_PHI_COLUMNS


INTEGER                     ::  NUM_R_ELEMS_PER_SHELL     ! NUM_SHELS*RAD_ELEMS_PER_SHELL = CHIMERA_R_ELEM
INTEGER                     ::  NUM_R_ELEMS_PER_SUBSHELL  ! NUM_R_ELEMS_PER_SUBSHELL = NUM_R_ELEMS_PER_BLOCK/NUM_SUBSHELLS_PER_SHELL


INTEGER                     ::  NUM_R_ELEMS_PER_BLOCK
INTEGER                     ::  NUM_T_ELEMS_PER_BLOCK     ! NUM_BLOCK_THETA_ROW*NUM_T_ELEMS_PER_BLOCK = CHIMERA_T_ELEM
INTEGER                     ::  NUM_P_ELEMS_PER_BLOCK     ! NUM_BLOCK_PHI_COLUMNS*NUM_P_ELEMS_PER_BLOCK = CHIMERA_P_ELEM


INTEGER                     ::  NUM_R_QUAD_POINTS
INTEGER                     ::  NUM_T_QUAD_POINTS
INTEGER                     ::  NUM_P_QUAD_POINTS

INTEGER                     ::  NUM_QUAD_DOF

INTEGER                     ::  R_COARSEN_FACTOR    = 1
INTEGER                     ::  T_COARSEN_FACTOR    = 1
INTEGER                     ::  P_COARSEN_FACTOR    = 1

INTEGER                     ::  NUM_BLOCKS
INTEGER                     ::  NUM_SUBSHELLS


INTEGER                     ::  POSEIDON_FRAME           = 0
INTEGER                     ::  CUR_ITERATION
INTEGER                     ::  MAX_ITERATIONS      = 5
INTEGER                     ::  CONVERGENCE_FLAG    = 0

REAL(KIND =idp)             ::  CONVERGENCE_CRITERIA

INTEGER                     ::  INITIAL_GUESS_FLAG

INTEGER                     ::  OUTPUT_SETUP_TABLE_FLAG     = 1

INTEGER                     ::  WRITE_TIMETABLE_FLAG        = 0
INTEGER                     ::  WRITE_REPORT_FLAG           = 0
INTEGER                     ::  ITER_REPORT_NUM_SAMPLES     = 20

INTEGER                     ::  WRITE_RESULTS_FLAG          = 1
INTEGER                     ::  WRITE_RESULTS_R_SAMPS       = 1000
INTEGER                     ::  WRITE_RESULTS_T_SAMPS       = 1
INTEGER                     ::  WRITE_RESULTS_P_SAMPS       = 1

INTEGER                     ::  WRITE_SOURCES_FLAG          = 1


INTEGER                     ::  RUN_REPORT_FILE_ID
INTEGER                     ::  ITER_REPORT_FILE_ID
INTEGER                     ::  FRAME_REPORT_FILE_ID

INTEGER                     ::  OUTPUT_MATRIX_FLAG          = 0
INTEGER                     ::  OUTPUT_RHS_VECTOR_FLAG      = 0
INTEGER                     ::  OUTPUT_UPDATE_VECTOR_FLAG   = 0

INTEGER                     ::  NEW_PETSC_SOLVER_FLAG
INTEGER                     ::  SOLVER_TYPE_FLAG            = 1 ! 1 = Regular N-R (Default), 2 = Jacobian-Free GMRES

LOGICAL                     ::  POSEIDON_INITIALIZED_FLAG   = .FALSE.

INTEGER                     ::  Ratio_T_BNDLperBLCK
INTEGER                     ::  Ratio_P_BNDLperBLCK
INTEGER                     ::  Ratio_BNDLperBLCK


END MODULE Poseidon_Parameters
