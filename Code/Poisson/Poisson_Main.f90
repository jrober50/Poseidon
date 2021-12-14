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

USE Poisson_Source_Vector, &
            ONLY :  Calculate_Poisson_Source_Vector

USE Poisson_Linear_Solve_Module, &
            ONLY : Poisson_Linear_Solve


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

IF ( Verbose_Flag ) THEN
    PRINT*,"Begining Newtonian Gravity Poisson Solve."
END IF


!!! Generate Src Vector !!!
CALL TimerStart( Timer_Poisson_SourceVector )
CALL Calculate_Poisson_Source_Vector()
CALL TimerStop( Timer_Poisson_SourceVector )



!!! Calculate Solution Coefficients !!!
CALL TimerStart( Timer_Poisson_LinearSolve )
CALL Poisson_Linear_Solve()
CALL TimerStop( Timer_Poisson_LinearSolve )




END SUBROUTINE Poisson_Solve






END MODULE Poisson_Main_Module
