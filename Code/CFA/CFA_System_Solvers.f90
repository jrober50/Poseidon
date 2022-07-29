   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE CFA_System_Solvers_Module                                       	     !##!
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

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag,               &
                    Eq_Flags

USE Poseidon_Message_Routines_Module, &
            ONLY :  Run_Message

USE Parameters_Variable_Indices, &
            ONLY :  iVB_X,                      &
                    iVB_S,                      &
                    iU_CF,                      &
                    iU_LF,                      &
                    iU_S1,                      &
                    iU_S2,                      &
                    iU_S3,                      &
                    iU_X1,                      &
                    iU_X2,                      &
                    iU_X3

USE Linear_System_Solvers_TypeA_Module, &
            ONLY :  Solve_Linear_System_TypeA

USE Linear_System_Solvers_TypeB_Module, &
            ONLY :  Solve_Linear_System_TypeB


USE Timer_Routines_Module, &
            ONLY :  TimerStart,                     &
                    TimerStop

USE Timer_Variables_Module, &
ONLY :  Timer_Lapse_LinearSolve,        &
        Timer_ConFactor_LinearSolve,    &
        Timer_Shift_LinearSolve

IMPLICIT NONE


CONTAINS



 !+101+####################################################!
!                                                           !
!          Solve_CFA_Systems                           	    !
!                                                           !
 !#########################################################!
SUBROUTINE Solve_CFA_Systems()

INTEGER, DIMENSION(3)                                   ::  iU
INTEGER                                                 ::  iVB



IF ( Eq_Flags(iU_CF) == 1 ) THEN
    CALL TimerStart( Timer_ConFactor_LinearSolve )
    CALL Solve_Linear_System_TypeA(iU_CF)
    CALL TimerStop( Timer_ConFactor_LinearSolve )
END IF




IF ( Eq_Flags(iU_LF) == 1 ) THEN
    CALL TimerStart( Timer_Lapse_LinearSolve )
    CALL Solve_Linear_System_TypeA(iU_LF)
    CALL TimerStop( Timer_Lapse_LinearSolve )
END IF




IF ( ANY(Eq_Flags(iU_S1:iU_S3) == 1) ) THEN
    
    iU = [iU_S1, iU_S2, iU_S3]
    iVB = iVB_S
    CALL TimerStart( Timer_Shift_LinearSolve )
    CALL Solve_Linear_System_TypeB( iU, iVB )
    CALL TimerStop(  Timer_Shift_LinearSolve )

END IF


END SUBROUTINE Solve_CFA_Systems








END MODULE CFA_System_Solvers_Module
