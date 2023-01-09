   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Memory_Variables_Module                                               !##!
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



#ifdef POSEIDON_MEMORY_FLAG

INTEGER, PUBLIC   ::  Memory_HWM

INTEGER, PUBLIC   ::  Memory_Start
INTEGER, PUBLIC   ::  Memory_After_MPI_Init
INTEGER, PUBLIC   ::  Memory_After_AMReX_Init
INTEGER, PUBLIC   ::  Memory_After_AMReX_Finalize
INTEGER, PUBLIC   ::  Memory_End

INTEGER, PUBLIC   ::  Memory_Loop_Start
INTEGER, PUBLIC   ::  Memory_Loop_Before_Init
INTEGER, PUBLIC   ::  Memory_Loop_After_Init
INTEGER, PUBLIC   ::  Memory_Loop_Before_Run
INTEGER, PUBLIC   ::  Memory_Loop_After_Run
INTEGER, PUBLIC   ::  Memory_Loop_Before_Close
INTEGER, PUBLIC   ::  Memory_Loop_End


INTEGER, PUBLIC   ::  Memory_FineMask_Before_LVB
INTEGER, PUBLIC   ::  Memory_FineMask_After_LVB

INTEGER, PUBLIC   ::  Memory_Method_Start
INTEGER, PUBLIC   ::  Memory_Method_X_Between
INTEGER, PUBLIC   ::  Memory_Method_Before_X_Load
INTEGER, PUBLIC   ::  Memory_Method_Before_Bnd_Factorize
INTEGER, PUBLIC   ::  Memory_Method_After_Bnd_Factorize



INTEGER, PUBLIC   ::  Memory_Method_Before_CF

INTEGER, PUBLIC   ::  Memory_Method_Before_CF_LoadVector1
INTEGER, PUBLIC   ::  Memory_Method_Before_CF_LoadVector2
INTEGER, PUBLIC   ::  Memory_Method_Before_CF_LoadVector3
INTEGER, PUBLIC   ::  Memory_Method_Before_CF_LoadVector4
INTEGER, PUBLIC   ::  Memory_Method_After_CF_LoadVector1


INTEGER, PUBLIC   ::  Memory_Method_Before_CF_Solve1

INTEGER, PUBLIC   ::  Memory_Method_Before_Chol_Fact
INTEGER, PUBLIC   ::  Memory_Method_After_Chol_Fact

INTEGER, PUBLIC   ::  Memory_Method_Before_CF_Solve2
INTEGER, PUBLIC   ::  Memory_Method_Before_CF_Solve3
INTEGER, PUBLIC   ::  Memory_Method_After_CF_Solve1
INTEGER, PUBLIC   ::  Memory_Method_After_CF_Solve2
INTEGER, PUBLIC   ::  Memory_Method_After_CF_Solve3
INTEGER, PUBLIC   ::  Memory_Method_After_CF_FixedPoint
INTEGER, PUBLIC   ::  Memory_Method_After_CF_DeallocWork

INTEGER, PUBLIC   ::  Memory_Method_Before_LF

INTEGER, PUBLIC   ::  Memory_Method_Before_LF_LoadVector1
INTEGER, PUBLIC   ::  Memory_Method_Before_LF_LoadVector2
INTEGER, PUBLIC   ::  Memory_Method_Before_LF_LoadVector3
INTEGER, PUBLIC   ::  Memory_Method_Before_LF_LoadVector4
INTEGER, PUBLIC   ::  Memory_Method_After_LF_LoadVector1
INTEGER, PUBLIC   ::  Memory_Method_Before_LF_Solve1
INTEGER, PUBLIC   ::  Memory_Method_Before_LF_Solve2
INTEGER, PUBLIC   ::  Memory_Method_Before_LF_Solve3
INTEGER, PUBLIC   ::  Memory_Method_After_LF_Solve1
INTEGER, PUBLIC   ::  Memory_Method_After_LF_Solve2
INTEGER, PUBLIC   ::  Memory_Method_After_LF_Solve3
INTEGER, PUBLIC   ::  Memory_Method_After_LF_FixedPoint
INTEGER, PUBLIC   ::  Memory_Method_After_LF_DeallocWork

INTEGER, PUBLIC   ::  Memory_Method_Before_SV
INTEGER, PUBLIC   ::  Memory_Method_End


INTEGER, PUBLIC   ::  Memory_Before_Init
INTEGER, PUBLIC   ::  Memory_After_Init

INTEGER, PUBLIC   ::  Memory_Before_Close
INTEGER, PUBLIC   ::  Memory_After_Close

INTEGER, PUBLIC   ::  Memory_Loop_After_Source
INTEGER, PUBLIC   ::  Memory_Loop_After_SetBC


#endif

CONTAINS


END MODULE Memory_Variables_Module

