   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Variables_Mesh                                                               !##!
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

INTEGER, PUBLIC                     :: NUM_R_ELEMENTS   !! # of Radial Elements
INTEGER, PUBLIC                     :: NUM_T_ELEMENTS   !! # of Theta Elements
INTEGER, PUBLIC                     :: NUM_P_ELEMENTS   !! # of Phi Elements

INTEGER, PUBLIC                     :: NUM_LOC_R_ELEMENTS   !! # of Radial Elements
INTEGER, PUBLIC                     :: NUM_LOC_T_ELEMENTS   !! # of Theta Elements
INTEGER, PUBLIC                     :: NUM_LOC_P_ELEMENTS   !! # of Phi Elements,

REAL(KIND = idp)                                            ::  R_INNER
REAL(KIND = idp)                                            ::  R_OUTER

REAL(KIND = idp), PUBLIC , ALLOCATABLE, DIMENSION(:)        ::  rlocs
REAL(KIND = idp), PUBLIC , ALLOCATABLE, DIMENSION(:)        ::  tlocs
REAL(KIND = idp), PUBLIC , ALLOCATABLE, DIMENSION(:)        ::  plocs

REAL(KIND = idp), PUBLIC , ALLOCATABLE, DIMENSION(:)        ::  drlocs
REAL(KIND = idp), PUBLIC , ALLOCATABLE, DIMENSION(:)        ::  dtlocs
REAL(KIND = idp), PUBLIC , ALLOCATABLE, DIMENSION(:)        ::  dplocs


INTEGER, PUBLIC                      ::  R_COARSEN_FACTOR    = 1
INTEGER, PUBLIC                      ::  T_COARSEN_FACTOR    = 1
INTEGER, PUBLIC                      ::  P_COARSEN_FACTOR    = 1


LOGICAL, PUBLIC, DIMENSION(3)       ::  locs_Set
LOGICAL, PUBLIC, DIMENSION(3)       ::  dlocs_Set


LOGICAL                                             ::  RADIAL_MESH_SET_FLAG = .FALSE.
LOGICAL                                             ::  THETA_MESH_SET_FLAG = .FALSE.
LOGICAL                                             ::  PHI_MESH_SET_FLAG = .FALSE.

END MODULE Variables_Mesh
