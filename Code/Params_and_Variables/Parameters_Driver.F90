   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE DRIVER_Parameters                                                           !##!
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



USE Poseidon_Kinds_Module, &
            ONLY : idp, fdp

USE Poseidon_Numbers_Module, &
            ONLY : pi
IMPLICIT NONE


!
!   DRIVER Set-up Variables
!
INTEGER                     ::  DRIVER_DIMENSION   =   3   ! Default = 3
INTEGER                     ::  DRIVER_TEST_NUMBER =   -1

INTEGER                     ::  DRIVER_R_ELEMS
INTEGER                     ::  DRIVER_C_ELEMS
INTEGER                     ::  DRIVER_T_ELEMS
INTEGER                     ::  DRIVER_P_ELEMS

INTEGER                     ::  DRIVER_R_INPUT_NODES
INTEGER                     ::  DRIVER_T_INPUT_NODES
INTEGER                     ::  DRIVER_P_INPUT_NODES

! x-space limits
REAL(KIND=idp)              ::  DRIVER_LEFT_LIMIT
REAL(KIND=idp)              ::  DRIVER_RIGHT_LIMIT

! r-space limits
REAL(KIND=idp)              ::  DRIVER_INNER_RADIUS
REAL(KIND=idp)              ::  DRIVER_CORE_RADIUS
REAL(KIND=idp)              ::  DRIVER_OUTER_RADIUS

INTEGER                     ::  SOLVER_MODE             ! 1 = NR, 2 = FP
INTEGER, DIMENSION(1:5)     ::  CFA_EQs_Flag_Vector

! Mesh Type
INTEGER                     ::  DRIVER_Mesh_Type
REAL(KIND=idp)              ::  DRIVER_Zoom

! MPI variables
INTEGER                     ::  DRIVER_PROCS
INTEGER                     ::  DRIVER_y_PROCS
INTEGER                     ::  DRIVER_z_PROCS
INTEGER                     ::  myID
INTEGER                     ::  myID_theta
INTEGER                     ::  myID_phi
INTEGER                     ::  nPROCS

INTEGER                     ::  ij_ray_dim
INTEGER                     ::  ik_ray_dim

INTEGER                     ::  DRIVER_SOLVER_TYPE = 2     ! Default = 2 -> CFA

INTEGER                     ::  MPI_COMM_XY
INTEGER                     ::  MPI_COMM_XZ
INTEGER                     ::  MPI_COMM_GRID
INTEGER                     ::  ngrid = 1


INTEGER                     ::  DRIVER_FIRST_GUESS_FLAG
INTEGER                     ::  DRIVER_SUBSEQUENT_GUESS_FLAG

INTEGER                     ::  OUTPUT_PRIMATIVES_FLAG

INTEGER                     ::  DRIVER_FRAME = 1
INTEGER                     ::  DRIVER_START_FRAME
INTEGER                     ::  DRIVER_END_FRAME
INTEGER                     ::  DRIVER_TOTAL_FRAMES

REAL(KIND=idp), ALLOCATABLE, DIMENSION(:)           ::  DRIVER_R_LOCS
REAL(KIND=idp), ALLOCATABLE, DIMENSION(:)           ::  DRIVER_T_LOCS
REAL(KIND=idp), ALLOCATABLE, DIMENSION(:)           ::  DRIVER_P_LOCS

REAL(KIND=idp), ALLOCATABLE, DIMENSION(:)           ::  DRIVER_Delta_R
REAL(KIND=idp), ALLOCATABLE, DIMENSION(:)           ::  DRIVER_Delta_T
REAL(KIND=idp), ALLOCATABLE, DIMENSION(:)           ::  DRIVER_Delta_P

REAL(KIND=idp), ALLOCATABLE, DIMENSION(:)           ::  Enclosed_Mass


REAL(KIND=idp), ALLOCATABLE, DIMENSION(:)           ::  DRIVER_T_LOCS_LOCAL
REAL(KIND=idp), ALLOCATABLE, DIMENSION(:)           ::  DRIVER_P_LOCS_LOCAL

REAL(KIND=idp), ALLOCATABLE, DIMENSION(:)           ::  DRIVER_POTENTIAL
REAL(KIND=idp), ALLOCATABLE, DIMENSION(:)           ::  DRIVER_SHIFT_VAL
REAL(KIND=idp), ALLOCATABLE, DIMENSION(:,:,:)       ::  DRIVER_LAPSE_VAL

REAL(KIND=idp), ALLOCATABLE, DIMENSION(:,:,:,:)     ::  DRIVER_E
REAL(KIND=idp), ALLOCATABLE, DIMENSION(:,:,:,:)     ::  DRIVER_S
REAL(KIND=idp), ALLOCATABLE, DIMENSION(:,:,:,:,:)   ::  DRIVER_Si



!===================================================================!
!                                                                   !
!   Analytic Solution Interfaces                                    !
!                                                                   !
!===================================================================!
LOGICAL                                                 ::   Potential_Sol_Flag
LOGICAL                                                 ::   Shift_Sol_Flag

!
!   Spherically Symmetric Variables
!
INTEGER                                                 :: POWER_A
REAL(KIND = idp)                                        :: RHO_O


INTEGER                                                 ::  CHIMERA_START_FRAME
INTEGER                                                 ::  CHIMERA_END_FRAME








END MODULE DRIVER_Parameters
