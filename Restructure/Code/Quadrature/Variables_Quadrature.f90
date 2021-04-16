   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Variables_Quadrature                                                             !##!
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

INTEGER, PUBLIC                      ::  NUM_R_QUAD_POINTS
INTEGER, PUBLIC                      ::  NUM_T_QUAD_POINTS
INTEGER, PUBLIC                      ::  NUM_P_QUAD_POINTS
INTEGER, PUBLIC                      ::  NUM_TP_QUAD_POINTS

INTEGER, PUBLIC                      ::  NUM_QUAD_DOF

REAL(KIND = idp), PUBLIC, ALLOCATABLE, DIMENSION(:)         ::  INT_R_LOCATIONS, INT_R_WEIGHTS
REAL(KIND = idp), PUBLIC, ALLOCATABLE, DIMENSION(:)         ::  INT_T_LOCATIONS, INT_T_WEIGHTS
REAL(KIND = idp), PUBLIC, ALLOCATABLE, DIMENSION(:)         ::  INT_P_LOCATIONS, INT_P_WEIGHTS

REAL(KIND = idp), PUBLIC, ALLOCATABLE, DIMENSION(:)         ::  INT_TP_WEIGHTS

REAL(KIND = idp), PUBLIC, ALLOCATABLE, DIMENSION(:)         ::  LOCAL_NODE_LOCATIONS



END MODULE Variables_Quadrature

