   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Variables_Flags                                                              !##!
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

!  Initialization Flags
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Num_Flags          = 6

INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Units_Set          = 1
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Matrices           = 2
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_IO_Params          = 3
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Quadrature         = 4
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_M_T_And_Gen_Vars   = 5
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Method_Vars        = 6
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Caller_Vars        = 7

LOGICAL,    PUBLIC, DIMENSION(1:iPF_Init_Num_Flags)     ::  lPF_Init_Flags

LOGICAL, PUBLIC         ::  Stiffness_Matrix_Initialized_Flag = .FALSE.
LOGICAL, PUBLIC         ::  Test_Space_Allocated_Flag = .FALSE.
LOGICAL, PUBLIC         ::  FirstCall_Flag = .TRUE.
LOGICAL, PUBLIC         ::  Matrix_Cholesky_Factorized_Flag = .FALSE.



!  Initialization Flags - Matrices
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Matrices_Num_Flags         = 4

INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Matrices_Type_A            = 1
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Matrices_Type_A_Cholesky   = 2
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Matrices_Type_B_Matrix     = 3
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Matrices_Type_B_LU         = 4

LOGICAL,    PUBLIC, DIMENSION(1:iPF_Init_Matrices_Num_Flags)    ::  lPF_Init_Matrices_Flags


!  Initialization Flags - Quadrature
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Quad_Num_Flags = 2

INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Quad_Vars      = 1
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Quad_Init      = 2

LOGICAL,    PUBLIC, DIMENSION(1:iPF_Init_Quad_Num_Flags)    ::  lPF_Init_Quad_Flags


!  Initialization Flags - Maps, Tables and General Variables (MTGV)
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_MTGV_Num_Flags = 3

INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_MTGV_Eq_Maps    = 1
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_MTGV_Derived    = 2
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_MTGV_Tables     = 3

LOGICAL,    PUBLIC, DIMENSION(1:iPF_Init_MTGV_Num_Flags)     ::  lPF_Init_MTGV_Flags




! Allocation Flags
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Alloc_Num_Flags         = 3

INTEGER,    PUBLIC, PARAMETER       ::  iPF_Alloc_Variables         = 1
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Alloc_XCFC_Variables    = 2
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Alloc_XCFC_Source_Vars  = 3

LOGICAL,    PUBLIC, DIMENSION(1:iPF_Alloc_Num_Flags)    ::  lPF_Alloc_Flags

LOGICAL, PUBLIC         ::  Poseidon_Variables_Allocated_Flag = .FALSE.





!  Initial_Guess_Flags
LOGICAL, PUBLIC         ::  Initial_Guess_Flag = .FALSE.
LOGICAL, PUBLIC         ::  Flat_Guess_Flag = .TRUE.



!  Boundary Condition Flags
INTEGER,    PUBLIC, PARAMETER       ::  iPF_BC_Num_Flags         = 3

INTEGER,    PUBLIC, PARAMETER       ::  iPF_BC_Set          = 1
INTEGER,    PUBLIC, PARAMETER       ::  iPF_BC_Inner_Set    = 2
INTEGER,    PUBLIC, PARAMETER       ::  iPF_BC_Outer_Set    = 3

LOGICAL,    PUBLIC, DIMENSION(1:iPF_BC_Num_Flags)    ::  lPF_BC_Flags



LOGICAL, PUBLIC         ::  Inner_BC_Set_Flag = .FALSE.
LOGICAL, PUBLIC         ::  Outer_BC_Set_Flag = .FALSE.


!  Debug Flags
!LOGICAL, PUBLIC         ::  Debug_Flag



CONTAINS





END MODULE Variables_Flags


