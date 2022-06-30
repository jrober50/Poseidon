   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Poisson_Main_Module                                                   !##!
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

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag

USE Poseidon_Message_Routines_Module, &
            ONLY :  Run_Message

USE Poisson_Load_Vector, &
            ONLY :  Calculate_Poisson_Load_Vector

USE XCFC_System_Solvers_TypeA_Module, &
            ONLY :  XCFC_Solve_System_TypeA

USE Timer_Routines_Module, &
            ONLY :  TimerStart,                     &
                    TimerStop

USE Timer_Variables_Module, &
            ONLY :  Timer_Poisson_SourceVector,        &
                    Timer_Poisson_LinearSolve

IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!                                                     				!
!                                                                               !
!###############################################################################!
SUBROUTINE Poisson_Solve()

IF ( Verbose_Flag ) CALL Run_Message("Begining Newtonian Gravity Poisson Solve.")


!!! Generate Src Vector !!!
CALL TimerStart( Timer_Poisson_SourceVector )
CALL Calculate_Poisson_Load_Vector()
CALL TimerStop( Timer_Poisson_SourceVector )



!!! Calculate Solution Coefficients !!!
CALL TimerStart( Timer_Poisson_LinearSolve )
CALL XCFC_Solve_System_TypeA(1)
CALL TimerStop( Timer_Poisson_LinearSolve )




END SUBROUTINE Poisson_Solve






END MODULE Poisson_Main_Module
