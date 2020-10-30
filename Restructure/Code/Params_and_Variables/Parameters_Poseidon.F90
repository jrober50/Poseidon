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



USE Poseidon_Kinds_Module, &
            ONLY : idp, fdp

USE Poseidon_Numbers_Module, &
            ONLY : pi

IMPLICIT NONE




INTEGER, PUBLIC                     ::  DOMAIN_DIM
INTEGER, PUBLIC                     ::  DEGREE
INTEGER, PUBLIC                     ::  L_LIMIT

INTEGER, PUBLIC                     ::  NUM_CFA_VARS = 5
INTEGER, PUBLIC                     ::  NUM_CFA_EQs

INTEGER, PUBLIC, PARAMETER          ::  DATA_DIST_MODE = 4
INTEGER, PUBLIC, PARAMETER          ::  SOL_DIST_SCHEME = 2
INTEGER, PUBLIC, PARAMETER          ::  STF_MAPPING_FLAG = 2



INTEGER, PUBLIC                     ::  POSEIDON_FRAME      = 0
INTEGER, PUBLIC                     ::  CUR_ITERATION
INTEGER, PUBLIC                     ::  MAX_ITERATIONS      = 5
INTEGER, PUBLIC                     ::  CONVERGENCE_FLAG    = 0

REAL(KIND =idp), PUBLIC             ::  CONVERGENCE_CRITERIA = 1.0E-8_idp

INTEGER, PUBLIC                     ::  INITIAL_GUESS_FLAG


INTEGER, PUBLIC                     ::  NEW_PETSC_SOLVER_FLAG
INTEGER, PUBLIC                     ::  Solver_Type            = 1 ! 1 = Regular N-R (Default), 2 = Jacobian-Free GMRES
CHARACTER(LEN = 32), PUBLIC         ::  Solver_Name



LOGICAL, PUBLIC                     ::  POSEIDON_INITIALIZED_FLAG   = .FALSE.

INTEGER, PUBLIC                     ::  Ratio_T_BNDLperBLCK
INTEGER, PUBLIC                     ::  Ratio_P_BNDLperBLCK
INTEGER, PUBLIC                     ::  Ratio_BNDLperBLCK


Logical, PUBLIC                     ::  Verbose_Flag

END MODULE Poseidon_Parameters
