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


!
!   CHIMERA Set-up Variables
!
INTEGER                     ::  CHIMERA_DIMENSION   =   3   ! Default = 3

INTEGER                     ::  CHIMERA_R_ELEMS
INTEGER                     ::  CHIMERA_C_ELEMS
INTEGER                     ::  CHIMERA_T_ELEMS
INTEGER                     ::  CHIMERA_P_ELEMS

INTEGER                     ::  CHIMERA_R_INPUT_NODES
INTEGER                     ::  CHIMERA_T_INPUT_NODES
INTEGER                     ::  CHIMERA_P_INPUT_NODES

! x-space limits
REAL(KIND=idp)              ::  CHIMERA_LEFT_LIMIT
REAL(KIND=idp)              ::  CHIMERA_RIGHT_LIMIT

! r-space limits
REAL(KIND=idp)              ::  CHIMERA_INNER_RADIUS
REAL(KIND=idp)              ::  CHIMERA_CORE_RADIUS
REAL(KIND=idp)              ::  CHIMERA_OUTER_RADIUS


! Mesh Type
INTEGER                     ::  CHIMERA_Mesh_Type

! MPI variables
INTEGER                     ::  CHIMERA_PROCS
INTEGER                     ::  CHIMERA_y_PROCS
INTEGER                     ::  CHIMERA_z_PROCS
INTEGER                     ::  myID
INTEGER                     ::  myID_theta
INTEGER                     ::  myID_phi
INTEGER                     ::  nPROCS

INTEGER                     ::  ij_ray_dim
INTEGER                     ::  ik_ray_dim

INTEGER                     ::  CHIMERA_SOLVER_TYPE = 2     ! Default = 2 -> CFA

INTEGER                     ::  MPI_COMM_XY
INTEGER                     ::  MPI_COMM_XZ
INTEGER                     ::  MPI_COMM_GRID
INTEGER                     ::  ngrid = 1

REAL(KIND=idp), ALLOCATABLE, DIMENSION(:)      ::  CHIMERA_R_LOCS
REAL(KIND=idp), ALLOCATABLE, DIMENSION(:)      ::  CHIMERA_T_LOCS
REAL(KIND=idp), ALLOCATABLE, DIMENSION(:)      ::  CHIMERA_P_LOCS

REAL(KIND=idp), ALLOCATABLE, DIMENSION(:)      ::  CHIMERA_Delta_R
REAL(KIND=idp), ALLOCATABLE, DIMENSION(:)      ::  CHIMERA_Delta_T
REAL(KIND=idp), ALLOCATABLE, DIMENSION(:)      ::  CHIMERA_Delta_P

REAL(KIND=idp), ALLOCATABLE, DIMENSION(:)      ::  Enclosed_Mass


REAL(KIND=idp), ALLOCATABLE, DIMENSION(:)      ::  CHIMERA_T_LOCS_LOCAL
REAL(KIND=idp), ALLOCATABLE, DIMENSION(:)      ::  CHIMERA_P_LOCS_LOCAL

REAL(KIND=idp), ALLOCATABLE, DIMENSION(:)      ::  CHIMERA_Potential

REAL(KIND=idp), ALLOCATABLE, DIMENSION(:,:,:)  ::  CHIMERA_E
REAL(KIND=idp), ALLOCATABLE, DIMENSION(:,:,:)  ::  CHIMERA_S
REAL(KIND=idp), ALLOCATABLE, DIMENSION(:,:,:,:)::  CHIMERA_Si

INTEGER                    ::  Ratio_T_BNDLperBLCK
INTEGER                    ::  Ratio_P_BNDLperBLCK
INTEGER                    ::  Ratio_BNDLperBLCK

!===================================================================!
!                                                                   !
!   Analytic Solution Interfaces                                    !
!                                                                   !
!===================================================================!
ABSTRACT INTERFACE
    FUNCTION Solution_Function_Pointer(r, theta, phi)
        REAL(KIND = KIND(1.D0))                         ::  Solution_Function_Pointer
        REAL(KIND = KIND(1.D0)), INTENT(IN)             ::  r, theta, phi
    END FUNCTION Solution_Function_Pointer
END INTERFACE


ABSTRACT INTERFACE
    FUNCTION Shift_Function_Pointer(r, r_locs, num_r_elems)
        REAL(KIND = KIND(1.D0))                         ::  Shift_Function_Pointer
        REAL(KIND = KIND(1.D0)), INTENT(IN)             ::  r
        INTEGER, INTENT(IN)                             ::  num_r_elems
        REAL(KIND = KIND(1.D0)), DIMENSION(0:num_r_elems), INTENT(IN) :: r_locs
    END FUNCTION Shift_Function_Pointer
END INTERFACE


PROCEDURE(Solution_Function_Pointer), POINTER           ::   Analytic_Solution => NULL()
PROCEDURE(Shift_Function_Pointer), POINTER              ::   Shift_Solution => NULL()



!
!   Spherically Symmetric Variables
!
INTEGER                                                 :: POWER_A
REAL(KIND = idp)                                        :: RHO_O


!
!   SelfSimilar Variables
!
REAL(KIND=idp)              ::  SELFSIM_T
REAL(KIND=idp)              ::  SELFSIM_KAPPA
REAL(KIND=idp)              ::  SELFSIM_GAMMA

INTEGER                                     :: NUM_ENTRIES
REAL(KIND=idp), ALLOCATABLE, DIMENSION(:)   :: SELFSIM_R_VALS
REAL(KIND=idp), ALLOCATABLE, DIMENSION(:)   :: SELFSIM_POT_VALS
REAL(KIND=idp), ALLOCATABLE, DIMENSION(:)   :: SELFSIM_SHIFT_VALS






END MODULE CHIMERA_Parameters
