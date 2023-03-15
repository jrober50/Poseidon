   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Timer_Variables_Module                                               !##!
!##!                                                                         !##!
!##!_________________________________________________________________________!##!
!##!                                                                         !##!
!##!                                                                         !##!
!##!=========================================================================!##!
!##!                                                                         !##!
!##!    Contains:                                                            !##!
!##!                                                                         !##!
!##!                                                                         !##!
!###############################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !##########################################################################!


!*D*================================!
!                                   !
!           Dependencies            !
!                                   !
!===================================!
USE Poseidon_Kinds_Module, &
            ONLY :  idp

IMPLICIT NONE

REAL(idp), PUBLIC   ::  Timer_Total
REAL(idp), PUBLIC   ::  Timer_Poseidon

REAL(idp), PUBLIC   ::  Timer_Core_Init_Test_Problem
REAL(idp), PUBLIC   ::  Timer_Core_PrintResults
REAL(idp), PUBLIC   ::  Timer_Core_Utilities

REAL(idp), PUBLIC   ::  Timer_Initialization_Total
REAL(idp), PUBLIC   ::  Timer_Initialization_Core
REAL(idp), PUBLIC   ::  Timer_Initialization_XCFC
REAL(idp), PUBLIC   ::  Timer_Initialization_Newtonian

REAL(idp), PUBLIC   ::  Timer_Poisson_SourceInput
REAL(idp), PUBLIC   ::  Timer_Poisson_SourceInput_PartA
REAL(idp), PUBLIC   ::  Timer_Poisson_SourceInput_PartB

REAL(idp), PUBLIC   ::  Timer_Poisson_SourceVector
REAL(idp), PUBLIC   ::  Timer_Poisson_SourceVector_SubParts
REAL(idp), PUBLIC   ::  Timer_Poisson_SourceVector_Main
REAL(idp), PUBLIC   ::  Timer_Poisson_LinearSolve


REAL(idp), PUBLIC   ::  Timer_GR_SourceInput
REAL(idp), PUBLIC   ::  Timer_GR_SourceInput_PartA
REAL(idp), PUBLIC   ::  Timer_GR_SourceInput_PartB

REAL(idp), PUBLIC   ::  Timer_Initialization


REAL(idp), PUBLIC   ::  Timer_Poisson_Initialization
REAL(idp), PUBLIC   ::  Timer_Poisson_Matrix_Init

REAL(idp), PUBLIC   ::  Timer_Matrix_Init
REAL(idp), PUBLIC   ::  Timer_Matrix_Radial_Terms
REAL(idp), PUBLIC   ::  Timer_Matrix_Angular_Terms
REAL(idp), PUBLIC   ::  Timer_Matrix_Laplace_Init
REAL(idp), PUBLIC   ::  Timer_Matrix_MVL_Init

REAL(idp), PUBLIC   ::  Timer_Matrix_Factorization
REAL(idp), PUBLIC   ::  Timer_Matrix_Cholesky

REAL(idp), PUBLIC   ::  Timer_Banded_Factorization

REAL(idp), PUBLIC   ::  Timer_X
REAL(idp), PUBLIC   ::  Timer_X_SourceVector
REAL(idp), PUBLIC   ::  Timer_X_LinearSolve

REAL(idp), PUBLIC   ::  Timer_Shift
REAL(idp), PUBLIC   ::  Timer_Shift_SourceVector
REAL(idp), PUBLIC   ::  Timer_Shift_LinearSolve

REAL(idp), PUBLIC   ::  Timer_Lapse
REAL(idp), PUBLIC   ::  Timer_Lapse_SourceVector
REAL(idp), PUBLIC   ::  Timer_Lapse_LinearSolve

REAL(idp), PUBLIC   ::  Timer_ConFactor
REAL(idp), PUBLIC   ::  Timer_ConFactor_SourceVector
REAL(idp), PUBLIC   ::  Timer_ConFactor_LinearSolve


REAL(idp), PUBLIC   ::  Timer_Remesh
REAL(idp), PUBLIC   ::  Timer_Remesh_MakeCopies
REAL(idp), PUBLIC   ::  Timer_Remesh_FillTotal
REAL(idp), PUBLIC   ::  Timer_Remesh_MakeLambdaArray
REAL(idp), PUBLIC   ::  Timer_Remesh_FillTypeA
REAL(idp), PUBLIC   ::  Timer_Remesh_FillX
REAL(idp), PUBLIC   ::  Timer_Remesh_FillS
REAL(idp), PUBLIC   ::  Timer_Remesh_DestroyCopies

REAL(idp), PUBLIC   ::  Timer_CFA_Load_Vector


REAL(idp), PUBLIC   ::  Timer_Driver_SetSource
REAL(idp), PUBLIC   ::  Timer_Driver_SetSource_InitTest
REAL(idp), PUBLIC   ::  Timer_Driver_SetSource_Scale
REAL(idp), PUBLIC   ::  Timer_Driver_SetSource_SetSource
REAL(idp), PUBLIC   ::  Timer_Driver_SetBC
REAL(idp), PUBLIC   ::  Timer_Driver_SetGuess
REAL(idp), PUBLIC   ::  Timer_Driver_Run
REAL(idp), PUBLIC   ::  Timer_Driver_Extra

REAL(idp), PUBLIC   ::  Timer_Plm_Init
REAL(idp), PUBLIC   ::  Timer_Ylm_Init


CONTAINS


END MODULE Timer_Variables_Module
