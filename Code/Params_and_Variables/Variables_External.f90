   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Variables_External                                                           !##!
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


REAL(KIND=idp), PUBLIC                                  ::  SELFSIM_T
REAL(KIND=idp), PUBLIC                                  ::  SELFSIM_START_T
REAL(KIND=idp), PUBLIC                                  ::  SELFSIM_END_T
REAL(KIND=idp), PUBLIC                                  ::  SELFSIM_NUM_FRAMES
REAL(KIND=idp), PUBLIC                                  ::  SELFSIM_KAPPA
REAL(KIND=idp), PUBLIC                                  ::  SELFSIM_GAMMA
REAL(KIND=idp), PUBLIC                                  ::  SELFSIM_ECC

INTEGER, PUBLIC                                         :: NUM_ENTRIES
REAL(KIND=idp), PUBLIC, ALLOCATABLE, DIMENSION(:)       :: SELFSIM_R_VALS
REAL(KIND=idp), PUBLIC, ALLOCATABLE, DIMENSION(:)       :: SELFSIM_POT_VALS
REAL(KIND=idp), PUBLIC, ALLOCATABLE, DIMENSION(:)       :: SELFSIM_SHIFT_VALS

INTEGER, PUBLIC                                         :: SELFSIM_V_SWITCH = 0

LOGICAL, PUBLIC                                         ::  SelfSim_Allocated = .FALSE.

INTEGER                                                 ::  OUTPUT_PRIMATIVES_FLAG = 1

REAL(KIND=idp), PUBLIC                                  ::  Central_E



INTEGER                                                 ::  MLS_RE

REAL(idp)                                               ::  MLS_SemiMinor
REAL(idp)                                               ::  MLS_SemiMajor
REAL(idp)                                               ::  MLS_Ecc
REAL(idp)                                               ::  MLS_Rho

INTEGER                                                 ::  MLS_SphereType
CHARACTER(LEN=7),   DIMENSION(2),   PARAMETER           ::  MLS_SphereName = &
                                                            ['Oblate ',  &
                                                             'Prolate'   ]
                                                             
INTEGER,                            PARAMETER           ::  iMLS_Oblate  = 1
INTEGER,                            PARAMETER           ::  iMLS_Prolate = 2




REAL(idp)                                               ::  MVL_Boundary_Value
REAL(idp)                                               ::  MVL_C1
REAL(idp)                                               ::  MVL_C2


REAL(idp)                                               ::  HCT_Alpha
REAL(idp)                                               ::  HCT_Star_Radius
REAL(idp)                                               ::  HCT_Rhoo
REAL(idp)                                               ::  HCT_Beta
REAL(idp)                                               ::  HCT_C


REAL(idp)                                               ::  UST_Rhoo
REAL(idp)                                               ::  UST_StarRadius


REAL(idp),  PUBLIC,     ALLOCATABLE,    DIMENSION(:)    :: CFLD_Update
REAL(idp),  PUBLIC,     ALLOCATABLE,    DIMENSION(:)    :: CFLD_Residual
INTEGER,    PUBLIC                                      :: CFLD_MaxIters  = 1
INTEGER,    PUBLIC                                      :: CFLD_Iters     = 0
REAL(idp),  PUBLIC                                      :: CFLD_Tolerance = 0.0_idp


REAL(idp)                                               ::  CCS_SurfaceRadius
REAL(idp)                                               ::  CCS_CoreDensity
REAL(idp)                                               ::  CCS_CoreRadius


END MODULE Variables_External



