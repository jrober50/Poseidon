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

REAL(idp), PUBLIC   ::  Timer_XCFC_Initialization
REAL(idp), PUBLIC   ::  Timer_XCFC_Matrix_Init

REAL(idp), PUBLIC   ::  Timer_Poisson_Initialization
REAL(idp), PUBLIC   ::  Timer_Poisson_Matrix_Init

REAL(idp), PUBLIC   ::  Timer_XCFC_Matrix_Factorization
REAL(idp), PUBLIC   ::  Timer_XCFC_Matrix_Cholesky
REAL(idp), PUBLIC   ::  Timer_XCFC_Banded_Factorization

REAL(idp), PUBLIC   ::  Timer_XCFC_X
REAL(idp), PUBLIC   ::  Timer_XCFC_X_SourceVector
REAL(idp), PUBLIC   ::  Timer_XCFC_X_LinearSolve

REAL(idp), PUBLIC   ::  Timer_XCFC_Shift
REAL(idp), PUBLIC   ::  Timer_XCFC_Shift_SourceVector
REAL(idp), PUBLIC   ::  Timer_XCFC_Shift_LinearSolve

REAL(idp), PUBLIC   ::  Timer_XCFC_Lapse
REAL(idp), PUBLIC   ::  Timer_XCFC_Lapse_SourceVector
REAL(idp), PUBLIC   ::  Timer_XCFC_Lapse_LinearSolve

REAL(idp), PUBLIC   ::  Timer_XCFC_ConFactor
REAL(idp), PUBLIC   ::  Timer_XCFC_ConFactor_SourceVector
REAL(idp), PUBLIC   ::  Timer_XCFC_ConFactor_LinearSolve


REAL(idp), PUBLIC   ::  Timer_FP_Initialization
REAL(idp), PUBLIC   ::  Timer_FP_Matrix_Init

REAL(idp), PUBLIC   ::  Timer_FP_Load_Vector
REAL(idp), PUBLIC   ::  Timer_FP_CFLF_Solve
REAL(idp), PUBLIC   ::  Timer_FP_Beta_Solve

REAL(idp), PUBLIC   ::  Timer_FP_Matrix_Cholesky
REAL(idp), PUBLIC   ::  Timer_FP_Banded_Factorization

REAL(idp), PUBLIC   ::  Timer_Driver_SetSource
REAL(idp), PUBLIC   ::  Timer_Driver_SetSource_InitTest
REAL(idp), PUBLIC   ::  Timer_Driver_SetSource_Scale
REAL(idp), PUBLIC   ::  Timer_Driver_SetSource_SetSource
REAL(idp), PUBLIC   ::  Timer_Driver_SetBC
REAL(idp), PUBLIC   ::  Timer_Driver_SetGuess
REAL(idp), PUBLIC   ::  Timer_Driver_Run
REAL(idp), PUBLIC   ::  Timer_Driver_Extra


CONTAINS


END MODULE Timer_Variables_Module
