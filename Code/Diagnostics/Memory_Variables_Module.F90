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

INTEGER, PUBLIC   ::  Memory_Max

INTEGER, PUBLIC   ::  Memory_Start
INTEGER, PUBLIC   ::  Memory_End

INTEGER, PUBLIC   ::  Memory_Loop_Start
INTEGER, PUBLIC   ::  Memory_Loop_Before_Init
INTEGER, PUBLIC   ::  Memory_Loop_Before_Run
INTEGER, PUBLIC   ::  Memory_Loop_After_Run
INTEGER, PUBLIC   ::  Memory_Loop_Before_Close
INTEGER, PUBLIC   ::  Memory_Loop_End


INTEGER, PUBLIC   ::  Memory_FineMask_Before_LVB
INTEGER, PUBLIC   ::  Memory_FineMask_After_LVB

#endif

CONTAINS


END MODULE Memory_Variables_Module

