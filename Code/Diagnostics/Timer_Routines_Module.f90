   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Timer_Routines_Module                                                !##!
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

USE Timer_Variables_Module

USE Flags_IO_Module, &
            ONLY :  lPF_IO_Flags,           &
                    iPF_IO_Write_Timetable

USE Timer_IO_Module, &
            ONLY :  Output_Time_Report

USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: &
  I8 => INT64



IMPLICIT NONE






CONTAINS


!+301+##########################################################################!
!                                                                               !
!           Init_Timers                                                          !
!                                                                               !
!###############################################################################!
SUBROUTINE Init_Timers

Timer_Total                         = 0.0_idp
Timer_Poseidon                      = 0.0_idp

Timer_Core_Init_Test_Problem        = 0.0_idp
Timer_Core_PrintResults             = 0.0_idp
Timer_Core_Utilities                = 0.0_idp


Timer_Initialization_Total           = 0.0_idp
Timer_Initialization_Core           = 0.0_idp
Timer_Initialization_XCFC           = 0.0_idp

Timer_Newtonian_Matrix_Init           = 0.0_idp
Timer_Newtonian_SourceInput           = 0.0_idp
Timer_Newtonian_SourceInput_PartA     = 0.0_idp
Timer_Newtonian_SourceInput_PartB     = 0.0_idp

Timer_Newtonian_LoadVector          = 0.0_idp
Timer_Newtonian_LoadVector_SubParts = 0.0_idp
Timer_Newtonian_LoadVector_Main     = 0.0_idp
Timer_Newtonian_LinearSolve           = 0.0_idp




Timer_GR_SourceInput                = 0.0_idp
Timer_GR_SourceInput_PartA          = 0.0_idp
Timer_GR_SourceInput_PartB          = 0.0_idp

Timer_Initialization           = 0.0_idp
Timer_Matrix_Init              = 0.0_idp
Timer_Matrix_Cholesky          = 0.0_idp
Timer_Banded_Factorization     = 0.0_idp

Timer_X                        = 0.0_idp
Timer_X_SourceVector           = 0.0_idp
Timer_X_LinearSolve            = 0.0_idp

Timer_Shift                    = 0.0_idp
Timer_Shift_SourceVector       = 0.0_idp
Timer_Shift_LinearSolve        = 0.0_idp

Timer_Lapse                    = 0.0_idp
Timer_Lapse_SourceVector       = 0.0_idp
Timer_Lapse_LinearSolve        = 0.0_idp

Timer_ConFactor                = 0.0_idp
Timer_ConFactor_SourceVector   = 0.0_idp
Timer_ConFactor_LinearSolve    = 0.0_idp

Timer_Remesh                    = 0.0_idp
Timer_Remesh_MakeCopies         = 0.0_idp
Timer_Remesh_FillTotal          = 0.0_idp
Timer_Remesh_MakeLambdaArray    = 0.0_idp
Timer_Remesh_FillTypeA          = 0.0_idp
Timer_Remesh_FillX              = 0.0_idp
Timer_Remesh_FillS              = 0.0_idp
Timer_Remesh_DestroyCopies      = 0.0_idp


Timer_Driver_SetSource              = 0.0_idp
Timer_Driver_SetSource_InitTest     = 0.0_idp
Timer_Driver_SetSource_SetSource    = 0.0_idp
Timer_Driver_SetSource_Scale        = 0.0_idp
Timer_Driver_SetBC                  = 0.0_idp
Timer_Driver_SetGuess               = 0.0_idp
Timer_Driver_Run                    = 0.0_idp
Timer_Driver_Extra                  = 0.0_idp

CALL TimerStart( Timer_Total )


END SUBROUTINE Init_Timers





!+302+##########################################################################!
!                                                                               !
!           Print_Timers                                                        !
!                                                                               !
!###############################################################################!
SUBROUTINE Finalize_Timers()


CALL TimerStop( Timer_Total )


Timer_Poseidon = Timer_Total - Timer_Core_Init_Test_Problem


! Matrix factorizations occur at two different times.
! Here the seperate factorization time are added to create a total.
Timer_Matrix_Factorization = Timer_Matrix_Cholesky        &
                                + Timer_Banded_Factorization

Timer_Initialization_Total = Timer_Initialization_Core      &
                           + Timer_Initialization_XCFC


! The first call to the X system and CF system contain factorization calls.
! Therefore we subtract those times out to get the time for just the source create
! and linear solve(s)
Timer_X            = Timer_X                      &
                        - Timer_Banded_Factorization

Timer_X_LinearSolve = Timer_X_LinearSolve         &
                         - Timer_Banded_Factorization

Timer_ConFactor    = Timer_ConFactor              &
                        - Timer_Matrix_Cholesky




CALL Output_Time_Report


END SUBROUTINE Finalize_Timers




!+101+##########################################################################!
!                                                                               !
!      TimersStart                                               				!
!                                                                               !
!###############################################################################!
SUBROUTINE TimerStart( Timer )

REAL(idp), INTENT(INOUT)            :: Timer

Timer = Timer - TimerWtime()

END SUBROUTINE TimerStart




!+102+##########################################################################!
!                                                                               !
!      TimersStop                                                               !
!                                                                               !
!###############################################################################!
SUBROUTINE TimerStop( Timer )

REAL(idp), INTENT(INOUT)            :: Timer

Timer = Timer + TimerWtime()

END SUBROUTINE TimerStop




!+201+##########################################################################!
!                                                                               !
!        TimersWtime                                                            !
!                                                                               !
!###############################################################################!
REAL(idp) FUNCTION TimerWtime()

INTEGER(I8)                ::  Clock_Read
INTEGER(I8)                ::  Clock_Rate
INTEGER(I8)                ::  Clock_Max

CALL SYSTEM_CLOCK( Clock_Read, Clock_Rate, Clock_Max)
TimerWtime = REAL( Clock_Read, idp )/REAL( Clock_Rate, idp )

RETURN

END FUNCTION TimerWtime










END MODULE Timer_Routines_Module
