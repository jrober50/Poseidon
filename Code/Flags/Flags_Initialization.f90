   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Flags_Initialization_Module                                                  !##!
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
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Num_Flags          = 7

INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Units_Set          = 1
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Matrices           = 2
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_IO_Params          = 3
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Quadrature         = 4
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_M_T_And_Gen_Vars   = 5
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Method_Vars        = 6
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Caller_Vars        = 7

LOGICAL,    PUBLIC, DIMENSION(1:iPF_Init_Num_Flags)     ::  lPF_Init_Flags





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




CONTAINS








END MODULE Flags_Initialization_Module


