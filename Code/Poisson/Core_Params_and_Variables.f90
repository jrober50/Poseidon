   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Global_Variables_And_Parameters                                              !##!
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


USE constants, &
            ONLY : idp, pi



IMPLICIT NONE

!###############################################################################!
!                                                                               !
!                           Solver Parameters                                   !
!                                                                               !
!###############################################################################!

INTEGER                                     :: NUM_R_ELEMENTS
INTEGER                                     :: NUM_T_ELEMENTS
INTEGER                                     :: NUM_P_ELEMENTS


INTEGER                                     :: DEGREE
INTEGER                                     :: L_LIMIT
INTEGER                                     :: LM_Length



CHARACTER(LEN = 4)                         :: MATRIX_FORMAT = 'CCS'



REAL(idp)                            :: R_INNER
REAL(idp)                            :: R_OUTER



!!!            Specify Source Function         !!!
!!                                              !!
!!  #                    Source                 !!
!!                                              !!
!!  1                   rho_0 r^a               !!
!!  2       rho_0 r^a : rho_0 = 0 for r > R_m   !!
!!  3               MacLaurin Spheroid          !!
!!  4              User Defined Source?         !!
!!                                              !!

INTEGER                                         :: Source_Function_Flag = 1





!!!             Specify Additional Problem Variables           !!!
!!                                                              !!
!!  POWER_A     - power in source functions 1 & 2               !!
!!  RHO_O       - Density constant in Source functions 1 & 2    !!
!!  POINT_MASS  - mass within inner boundary                    !!
!!


INTEGER                                         :: POWER_A = 0
REAL(idp)                                :: RHO_O = 1.0_idp
REAL(idp)                                :: STAR_SURFACE = 0.50_idp





CHARACTER(LEN = 7), PARAMETER                   :: SPHEROID_TYPE = 'OBLATE'    !  Choices : OBLATE, PROLATE

REAL(idp)                                :: R_MIN
REAL(idp)                                :: R_MAX










!###############################################################################################!
!                                                                                               !
!                                           Derived Parameters                                  !
!                                                                                               !
!###############################################################################################!

INTEGER                                     :: NUM_R_NODES  != DEGREE*NUM_R_ELEMENTS+1








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
LOGICAL                                             ::  Test_Run_Flag = .FALSE.



LOGICAL                                             ::  Matrix_Cholesky_Factorized_Flag = .FALSE.


!===================================================================!
!                                                                   !
!   Solver Variables                                                !
!                                                                   !
!===================================================================!
CHARACTER(LEN = 4)                                  ::  LINEAR_SOLVER = "CHOL"


!===================================================================!
!                                                                   !
!   Mesh Variables                                                  !
!                                                                   !
!===================================================================!
REAL(idp), ALLOCATABLE, DIMENSION(:)                :: rlocs
REAL(idp), ALLOCATABLE, DIMENSION(:)                :: tlocs
REAL(idp), ALLOCATABLE, DIMENSION(:)                :: plocs

INTEGER                                             :: Discontinuity_Flag

LOGICAL                                             ::  RADIAL_MESH_SET_FLAG = .FALSE.
LOGICAL                                             ::  THETA_MESH_SET_FLAG = .FALSE.
LOGICAL                                             ::  PHI_MESH_SET_FLAG = .FALSE.






!===================================================================!
!                                                                   !
!   Coefficient Variable                                            !
!                                                                   !
!===================================================================!
COMPLEX(idp), ALLOCATABLE, DIMENSION(:,:,:)     :: Coefficient_Vector









!===================================================================!
!                                                                   !
!   Source Vector Variable                                          !
!                                                                   !
!===================================================================!
COMPLEX(idp), ALLOCATABLE, DIMENSION(:,:)     :: Source_Vector










!===================================================================!
!                                                                   !
!   Stiffness Matrix Variables                                      !
!                                                                   !
!===================================================================!
ABSTRACT INTERFACE
    SUBROUTINE Subroutine_No_Input()
    END SUBROUTINE Subroutine_No_Input
END INTERFACE

PROCEDURE(Subroutine_No_Input), POINTER             :: St => NULL()


REAL(idp), ALLOCATABLE, DIMENSION(:,:,:)            :: STF_MAT_Integrals


!!! STF_MAT in full matrix form !!!
REAL(idp), ALLOCATABLE, DIMENSION(:,:,:)            :: STF_MAT



!!! STF_MAT in CCS form !!!
INTEGER                                             :: STF_NNZ
REAL(idp), ALLOCATABLE, DIMENSION(:,:)              :: STF_ELEM_VAL
INTEGER, ALLOCATABLE, DIMENSION(:)                  :: STF_COL_PTR, STF_ROW_IND











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


REAL(idp), DIMENSION(:,:,:,:)    , ALLOCATABLE          ::  Test_Source_Input
REAL(idp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE          ::  Source_Term_Coefficients
REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              ::  Source_Terms

INTEGER,DIMENSION(1:3)                                  ::  Source_Degrees
INTEGER                                                 ::  Num_Quad_DOF








!===================================================================!
!                                                                   !
!   Analytic_Solution                                               !
!                                                                   !
!===================================================================!
ABSTRACT INTERFACE
    FUNCTION Solution_Function_Pointer(r, theta, phi)
        REAL(KIND = KIND(1.D0))                         ::  Solution_Function_Pointer
        REAL(KIND = KIND(1.D0)), INTENT(IN)             ::  r, theta, phi
    END FUNCTION Solution_Function_Pointer
END INTERFACE

PROCEDURE(Solution_Function_Pointer), POINTER          ::   Analytic_Solution => NULL()









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

COMPLEX(idp)                                    :: INNER_DIR_BC_INPUT
COMPLEX(idp)                                    :: INNER_NEU_BC_INPUT


!!! OUTER !!!
CHARACTER(LEN = 4)                              :: OUTER_BC_TYPE
LOGICAL                                         :: OUTER_BC_SET_FLAG = .FALSE.
LOGICAL                                         :: OUTER_UNIFORM_DIR_BC_FLAG

COMPLEX(idp)                                    :: OUTER_DIR_BC_INPUT
COMPLEX(idp)                                    :: OUTER_NEU_BC_INPUT





!!! Used to store origonal matrix to allow for Dirichlet BC updates !!!
!!!
!!! Currently allocated in CHOLESKY_FACTORIZATION subroutine in Linear_Solvers_And_Preconditioners.f90
REAL(idp), ALLOCATABLE, DIMENSION(:,:)   ::  First_Column_Storage
REAL(idp), ALLOCATABLE, DIMENSION(:,:)   ::  Last_Column_Storage





!===================================================================!
!                                                                   !
!   Unit Variables                                                  !
!                                                                   !
!===================================================================!

REAL(idp)                                   ::  Meter,              &
                                                Kilometer,          &
                                                Gram



REAL(idp)                                   ::  Grav_Constant_G,    &
                                                Speed_of_Light






END MODULE Global_Variables_And_Parameters
