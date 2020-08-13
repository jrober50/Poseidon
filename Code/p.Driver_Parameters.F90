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


USE Poseidon_Constants_Module, &
            ONLY :  idp, pi, fdp


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

INTEGER                     ::  SOURCE_OUTPUT_FLAG
INTEGER                     ::  RESULTS_OUTPUT_FLAG
INTEGER                     ::  RUN_REPORT_FLAG
INTEGER                     ::  FRAME_REPORT_FLAG

INTEGER                     ::  DRIVER_FIRST_GUESS_FLAG
INTEGER                     ::  DRIVER_SUBSEQUENT_GUESS_FLAG


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


PROCEDURE(Solution_Function_Pointer), POINTER           ::   Potential_Solution => NULL()
PROCEDURE(Shift_Function_Pointer), POINTER              ::   Shift_Solution => NULL()
LOGICAL                                                 ::   Potential_Sol_Flag
LOGICAL                                                 ::   Shift_Sol_Flag

!
!   Spherically Symmetric Variables
!
INTEGER                                                 :: POWER_A
REAL(KIND = idp)                                        :: RHO_O


INTEGER                                                 ::  CHIMERA_START_FRAME
INTEGER                                                 ::  CHIMERA_END_FRAME


!
!   SelfSimilar Variables
!
REAL(KIND=idp)                                          ::  SELFSIM_T
REAL(KIND=idp)                                          ::  SELFSIM_START_T
REAL(KIND=idp)                                          ::  SELFSIM_END_T
REAL(KIND=idp)                                          ::  SELFSIM_NUM_FRAMES
REAL(KIND=idp)                                          ::  SELFSIM_KAPPA
REAL(KIND=idp)                                          ::  SELFSIM_GAMMA
REAL(KIND=idp)                                          ::  SELFSIM_ECC

INTEGER                                                 :: NUM_ENTRIES
REAL(KIND=idp), ALLOCATABLE, DIMENSION(:)               :: SELFSIM_R_VALS
REAL(KIND=idp), ALLOCATABLE, DIMENSION(:)               :: SELFSIM_POT_VALS
REAL(KIND=idp), ALLOCATABLE, DIMENSION(:)               :: SELFSIM_SHIFT_VALS

INTEGER                                                 :: SELFSIM_V_SWITCH






INTEGER, DIMENSION(:), ALLOCATABLE                      ::  Iteration_History
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  Driver_Iter_Time_Table
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  Driver_Frame_Time_Table
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE             ::  Driver_Run_Time_Table






END MODULE DRIVER_Parameters
