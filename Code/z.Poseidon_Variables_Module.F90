   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_Variables_Module                                                    !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains the global variables that store the stiffness matrix, source       !##!
!##!           vector, and coefficient vector.  Also location of global varaibles   !##!
!##!           used to define the parameters of the solve.                          !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!


USE Poseidon_Constants_Module, &
            ONLY : idp, pi



IMPLICIT NONE

!###############################################################################################!
!                                                                                               !
!                                           Solver Parameters                                   !
!                                                                                               !
!###############################################################################################!

!!! Specify the Number of Dimensions in the Domain



!!! Specify the Number of Elements in Each Dimension !!!
INTEGER                                     :: NUM_R_ELEMENTS   !! # of Radial Elements
INTEGER                                     :: NUM_T_ELEMENTS   !! # of Theta Elements
INTEGER                                     :: NUM_P_ELEMENTS   !! # of Phi Elements



CHARACTER(LEN = 3)                          :: PHYSICS_TYPE = 'CFA'     ! NEW - Newton , CFA - CFA





!###############################################################################################!
!                                                                                               !
!                                           Problem Parameters                                  !
!                                                                                               !
!###############################################################################################!

!!! Specify Radial Domain !!!
REAL(KIND = idp)                            :: R_INNER
REAL(KIND = idp)                            :: R_OUTER





!###############################################################################################!
!                                                                                               !
!                                           Derived Parameters                                  !
!                                                                                               !
!###############################################################################################!

INTEGER                                     :: NUM_R_NODES
INTEGER                                     :: BLOCK_NUM_R_NODES
INTEGER                                     :: SUBSHELL_NUM_R_NODES

INTEGER                                     :: VAR_DIM
INTEGER                                     :: ELEM_VAR_DIM
INTEGER                                     :: BLOCK_VAR_DIM
INTEGER                                     :: SUBSHELL_VAR_DIM

INTEGER                                     :: PROB_DIM
INTEGER                                     :: ELEM_PROB_DIM
INTEGER                                     :: ELEM_PROB_DIM_SQR
INTEGER                                     :: BLOCK_PROB_DIM
INTEGER                                     :: SUBSHELL_PROB_DIM


INTEGER                                     ::  NUM_TP_QUAD_POINTS



INTEGER                                     :: LM_LENGTH
INTEGER                                     :: ULM_LENGTH

INTEGER, DIMENSION(:), ALLOCATABLE          :: M_VALUES



INTEGER                                                     ::  Ratio_BNDLperBLCK,      &
                                                                Ratio_T_BNDLperBLCK,    &
                                                                Ratio_P_BNDLperBLCK


ABSTRACT INTERFACE
    PURE FUNCTION Matrix_Location_Pointer(ui, l, m, re, d)
        INTEGER                                             ::  Matrix_Location_Pointer
        INTEGER, INTENT(IN)                                 ::  ui, l, m, re, d
    END FUNCTION Matrix_Location_Pointer
END INTERFACE

PROCEDURE(Matrix_Location_Pointer), POINTER                 ::   Matrix_Location => NULL()


ABSTRACT INTERFACE
    PURE FUNCTION LM_Location_Pointer(l, m)
        INTEGER                                             ::  LM_Location_Pointer
        INTEGER, INTENT(IN)                                 ::  l, m
    END FUNCTION LM_Location_Pointer
END INTERFACE

PROCEDURE(LM_Location_Pointer), POINTER                     ::   LM_Location => NULL()




!###############################################################################################!
!                                                                                               !
!                                           Global Variables                                    !
!                                                                                               !
!###############################################################################################!


!===================================================================!
!                                                                   !
!   Initialization Flags                                            !
!                                                                   !
!===================================================================!
LOGICAL                                             ::  Stiffness_Matrix_Initialized_Flag = .FALSE.
LOGICAL                                             ::  Test_Space_Allocated_Flag = .FALSE.
LOGICAL                                             ::  FirstCall_Flag = .TRUE.
LOGICAL                                             ::  Matrix_Cholesky_Factorized_Flag = .FALSE.


!===================================================================!
!                                                                   !
!   Integration Variables                                           !
!                                                                   !
!===================================================================!
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)         ::  INT_R_LOCATIONS, INT_R_WEIGHTS
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)         ::  INT_T_LOCATIONS, INT_T_WEIGHTS
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)         ::  INT_P_LOCATIONS, INT_P_WEIGHTS

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)         ::  INT_TP_WEIGHTS

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)         ::  LOCAL_NODE_LOCATIONS




!===================================================================!
!                                                                   !
!   Mesh Variables                                                  !
!                                                                   !
!===================================================================!
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)         ::  rlocs
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)         ::  tlocs
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)         ::  plocs


LOGICAL                                             ::  RADIAL_MESH_SET_FLAG = .FALSE.
LOGICAL                                             ::  THETA_MESH_SET_FLAG = .FALSE.
LOGICAL                                             ::  PHI_MESH_SET_FLAG = .FALSE.







!===================================================================!
!                                                                   !
!   Ylm Table                                                       !
!                                                                   !
!===================================================================!
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:,:,:,:)        :: Ylm_Table_Block


COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:,:)          :: Ylm_Values
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:,:)          :: Ylm_dt_Values
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:,:)          :: Ylm_dp_Values

COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:,:)          :: Ylm_CC_Values
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:,:)          :: Ylm_CC_DT_Values
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:,:)          :: Ylm_CC_DP_Values





!===================================================================!
!                                                                   !
!   Lagrange_Poly Table                                             !
!                                                                   !
!===================================================================!
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)                 :: Lagrange_Poly_Table
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:,:,:)             :: LPT_LPT



!===================================================================!
!                                                                   !
!   Coefficient Variable                                            !
!                                                                   !
!===================================================================!
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)                  :: Coefficient_Vector
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)                  :: Update_Vector







!===================================================================!
!                                                                   !
!   RHS Vector Variable                                          !
!                                                                   !
!===================================================================!
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)  :: RHS_Vector
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)  :: Block_RHS_Vector







!===================================================================!
!                                                                   !
!   Stiffness Matrix Variables                                      !
!                                                                   !
!===================================================================!
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:)    ::  Block_STF_Mat
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:)    ::  BLOCK_ELEM_STF_MATVEC

INTEGER                                             ::  NUM_OFF_DIAGONALS



!!! Newtonian Variables
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)     ::  STF_MAT_Integrals


!!! STF_MAT in full matrix form !!!
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)     ::  STF_MAT



!!! STF_MAT in CCS form !!!
INTEGER                                             ::  STF_NNZ
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)       ::  STF_ELEM_VAL
INTEGER, ALLOCATABLE, DIMENSION(:)                  ::  STF_COL_PTR, STF_ROW_IND


!===================================================================!
!                                                                   !
!   Source Function Variables                                       !
!                                                                   !
!===================================================================!
ABSTRACT INTERFACE
    FUNCTION Source_Function_Pointer(r, theta, phi)
        REAL(KIND = KIND(1.D0))                         ::  Source_Function_Pointer
        REAL(KIND = KIND(1.D0)), INTENT(IN)             ::  r, theta, phi
    END FUNCTION Source_Function_Pointer
END INTERFACE

PROCEDURE(Source_Function_Pointer), POINTER             ::  Source_Function => NULL()








!===================================================================!
!                                                                   !
!   Source Vectors                                                  !
!                                                                   !
!===================================================================!
REAL(KIND = idp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE   ::  Source_Term_Coefficients

REAL(KIND = idp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE   ::  Block_Source_E
REAL(KIND = idp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE   ::  Block_Source_S
REAL(KIND = idp), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE ::  Block_Source_Si












!===================================================================!
!                                                                   !
!   Boundary Condtion Variables                                     !
!                                                                   !
!===================================================================!


!!!   Specify Boundary Conditions  !!!
!!                                  !!
!!  DRCH, D  -> Dirichlet           !!
!!  NEUM, N  -> Neumann !CCS Works! !!
!!                                  !!


!!! INNER !!!
CHARACTER(LEN = 4)                              :: INNER_BC_TYPE
LOGICAL                                         :: INNER_BC_SET_FLAG = .FALSE.
LOGICAL                                         :: INNER_UNIFORM_DIR_BC_FLAG

COMPLEX(KIND = idp)                             :: INNER_DIR_BC_INPUT
COMPLEX(KIND = idp)                             :: INNER_NEU_BC_INPUT


!!! OUTER !!!
CHARACTER(LEN = 4)                              :: OUTER_BC_TYPE
LOGICAL                                         :: OUTER_BC_SET_FLAG = .FALSE.
LOGICAL                                         :: OUTER_UNIFORM_DIR_BC_FLAG

COMPLEX(KIND = idp)                             :: OUTER_DIR_BC_INPUT
COMPLEX(KIND = idp)                             :: OUTER_NEU_BC_INPUT





CHARACTER(LEN = 1), DIMENSION(1:5)              :: INNER_CFA_BC_TYPE
REAL(KIND = idp), DIMENSION(1:5)                :: INNER_CFA_BC_VALUES


CHARACTER(LEN = 1), DIMENSION(1:5)              :: OUTER_CFA_BC_TYPE
REAL(KIND = idp), DIMENSION(1:5)                :: OUTER_CFA_BC_VALUES







!===================================================================!
!                                                                   !
!   MPI & Parallelism Variables                                     !
!                                                                   !
!===================================================================!

INTEGER                                                     ::  POSEIDON_COMM_WORLD
INTEGER                                                     ::  POSEIDON_COMM_SHELL
INTEGER                                                     ::  POSEIDON_COMM_DIST
INTEGER                                                     ::  POSEIDON_COMM_PETSC

INTEGER                                                     ::  nPROCS_PETSC
INTEGER                                                     ::  nPROCS_SHELL


INTEGER                                                     ::  POSEIDON_DATA_DIST_MODE

INTEGER                                                     ::  ierr,                   &
                                                                NUM_BLOCKS_THETA,       &
                                                                NUM_BLOCKS_PHI


INTEGER                                                     ::  myShell = -1

INTEGER                                                     ::  myID_Poseidon
INTEGER                                                     ::  myID_Shell      = -1
INTEGER                                                     ::  myID_SubShell   = -1
INTEGER                                                     ::  myID_Dist       = -1
INTEGER                                                     ::  myID_PETSc      = -1


INTEGER                                                     ::  Local_Length = -1


INTEGER                                                     ::  Local_Theta_Start,      &
                                                                Local_Theta_Range,      &
                                                                Local_Phi_Start,        &
                                                                Local_Phi_Range






!###############################################################################################!
!                                                                                               !
!                                     Output/Timing Variables                                   !
!                                                                                               !
!###############################################################################################!
INTEGER                                                     ::  Num_Timer_Calls=25
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Iter_Time_Table
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Frame_Time_Table
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Run_Time_Table

INTEGER                                                     ::  Total_Run_Iters=1

REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Frame_Convergence_Table

INTEGER, DIMENSION(:), ALLOCATABLE                          ::  Iteration_Histogram



END MODULE Poseidon_Variables_Module
