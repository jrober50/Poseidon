   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Newtonian_Main_Module                                                   !##!
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

USE Parameters_Variable_Indices, &
            ONLY :  iU_NP

USE Poseidon_Message_Routines_Module, &
            ONLY :  Run_Message

USE Load_Vector_Newtonian_Module, &
            ONLY :  Create_Load_Vector_Newtonian

USE Linear_System_Solvers_TypeA_Module, &
            ONLY :  Solve_Linear_System_TypeA

USE Timer_Routines_Module, &
            ONLY :  TimerStart,                     &
                    TimerStop

USE Timer_Variables_Module, &
            ONLY :  Timer_Newtonian_LoadVector,        &
                    Timer_Newtonian_LinearSolve

IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!                                                     				!
!                                                                               !
!###############################################################################!
SUBROUTINE Newtonian_Method()

IF ( Verbose_Flag ) CALL Run_Message("Begining Newtonian Gravity Solve.")


!!! Generate Src Vector !!!
CALL TimerStart( Timer_Newtonian_LoadVector )
CALL Create_Load_Vector_Newtonian()
CALL TimerStop( Timer_Newtonian_LoadVector )




!!! Calculate Solution Coefficients !!!
CALL TimerStart( Timer_Newtonian_LinearSolve )
CALL Solve_Linear_System_TypeA(iU_NP)
CALL TimerStop( Timer_Newtonian_LinearSolve )


END SUBROUTINE Newtonian_Method






END MODULE Newtonian_Main_Module
