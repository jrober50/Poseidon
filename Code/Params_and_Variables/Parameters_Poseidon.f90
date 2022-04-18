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
INTEGER, DIMENSION(1:5)             ::  CFA_EQ_Flags


INTEGER, PUBLIC, PARAMETER          ::  DATA_DIST_MODE = 4
INTEGER, PUBLIC, PARAMETER          ::  SOL_DIST_SCHEME = 2
INTEGER, PUBLIC, PARAMETER          ::  STF_MAPPING_FLAG = 2



INTEGER, PUBLIC                     ::  POSEIDON_FRAME          = 1
INTEGER, PUBLIC                     ::  CUR_ITERATION           = 1
INTEGER, PUBLIC                     ::  MAX_ITERATIONS          = 20
INTEGER, PUBLIC, PARAMETER          ::  MAX_ITERATIONS_Default  = 20

INTEGER, PUBLIC                     ::  CONVERGENCE_FLAG     = 0
INTEGER, PUBLIC                     ::  Convergence_Type     = 2
REAL(KIND =idp), PUBLIC             ::  CONVERGENCE_CRITERIA = 1.0E-8_idp
REAL(KIND =idp), PUBLIC, PARAMETER  ::  CONVERGENCE_CRITERIA_Default = 1.0E-8_idp

INTEGER, PUBLIC                     ::  INITIAL_GUESS_FLAG


INTEGER, PUBLIC                     ::  NEW_PETSC_SOLVER_FLAG
INTEGER, PUBLIC                     ::  Method_Flag          = 1 ! 1 = Regular N-R (Default), 2 = Jacobian-Free GMRES



LOGICAL, PUBLIC                     ::  POSEIDON_INITIALIZED_FLAG   = .FALSE.

INTEGER, PUBLIC                     ::  Ratio_T_BNDLperBLCK
INTEGER, PUBLIC                     ::  Ratio_P_BNDLperBLCK
INTEGER, PUBLIC                     ::  Ratio_BNDLperBLCK


LOGICAL, PUBLIC                     ::  Verbose_Flag = .FALSE.

LOGICAL, PUBLIC                     ::  Poisson_Mode = .FALSE.

LOGICAL, PUBLIC                     ::  Poseidon_Remesh_Flag = .TRUE.

CHARACTER(17), DIMENSION(3), PUBLIC, PARAMETER        ::  &
Convergence_Type_Names = ['L_One Error Norm ',     &
                          'L_Two Error Norm ',     &
                          'L_Inf Error Norm '      ]

END MODULE Poseidon_Parameters
