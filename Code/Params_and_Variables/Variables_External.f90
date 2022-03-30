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

INTEGER                                                 ::  OUTPUT_PRIMATIVES_FLAG = 0

REAL(KIND=idp), PUBLIC                                  ::  Central_E



INTEGER                                                 ::  MacLaurin_RE

REAL(idp)                                               ::  MacLaurin_SemiMinor
REAL(idp)                                               ::  MacLaurin_SemiMajor
REAL(idp)                                               ::  MacLaurin_Ecc
REAL(idp)                                               ::  MacLaurin_Rho

CHARACTER(LEN=1)                                        ::  MacLaurin_SphereType


END MODULE Variables_External



