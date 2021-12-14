   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Initilaization_Poseidon_Core                                          !##!
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


IMPLICIT NONE


CONTAINS



 !+101+########################################################!
!                                                               !
!          Initialize_Poseidon_Core                             !
!                                                               !
 !#############################################################!
SUBROUTINE Initialize_Poseidon_Core()


ALLOCATE( ITER_TIME_TABLE(1:NUM_TIMER_CALLS) )
ALLOCATE( FRAME_TIME_TABLE(1:NUM_TIMER_CALLS) )
ALLOCATE( RUN_TIME_TABLE(1:NUM_TIMER_CALLS) )

ITER_TIME_TABLE = 0.0_idp
FRAME_TIME_TABLE = 0.0_idp
RUN_TIME_TABLE = 0.0_idp

ALLOCATE( Frame_Update_Table(1:MAX_ITERATIONS,1:5) )
ALLOCATE( Frame_Residual_Table(1:3, 1:MAX_ITERATIONS,1:5) )
ALLOCATE( Iteration_Histogram(1:MAX_ITERATIONS) )

Frame_Update_Table = 0.0_idp
Frame_Residual_Table = 0.0_idp
Iteration_Histogram = 0

END SUBROUTINE Initialize_Poseidon_Core




END MODULE Initilaization_Poseidon_Core
