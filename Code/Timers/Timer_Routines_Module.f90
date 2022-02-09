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

USE Variables_IO, &
            ONLY :  Report_Flags,       &
                    iRF_Time

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
Timer_Core_Initialization           = 0.0_idp


Timer_Poisson_Matrix_Init           = 0.0_idp
Timer_Poisson_SourceInput           = 0.0_idp
Timer_Poisson_SourceInput_PartA     = 0.0_idp
Timer_Poisson_SourceInput_PartB     = 0.0_idp

Timer_Poisson_SourceVector          = 0.0_idp
Timer_Poisson_SourceVector_SubParts = 0.0_idp
Timer_Poisson_SourceVector_Main     = 0.0_idp
Timer_Poisson_LinearSolve           = 0.0_idp




Timer_GR_SourceInput                = 0.0_idp
Timer_GR_SourceInput_PartA          = 0.0_idp
Timer_GR_SourceInput_PartB          = 0.0_idp

Timer_XCFC_Initialization           = 0.0_idp
Timer_XCFC_Matrix_Init              = 0.0_idp
Timer_XCFC_Matrix_Cholesky          = 0.0_idp
Timer_XCFC_Type_B_Factorization     = 0.0_idp

Timer_XCFC_X                        = 0.0_idp
Timer_XCFC_X_SourceVector           = 0.0_idp
Timer_XCFC_X_LinearSolve            = 0.0_idp

Timer_XCFC_Shift                    = 0.0_idp
Timer_XCFC_Shift_SourceVector       = 0.0_idp
Timer_XCFC_Shift_LinearSolve        = 0.0_idp

Timer_XCFC_Lapse                    = 0.0_idp
Timer_XCFC_Lapse_SourceVector       = 0.0_idp
Timer_XCFC_Lapse_LinearSolve        = 0.0_idp

Timer_XCFC_ConFactor                = 0.0_idp
Timer_XCFC_ConFactor_SourceVector   = 0.0_idp
Timer_XCFC_ConFactor_LinearSolve    = 0.0_idp


Timer_FP_Initialization             = 0.0_idp

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



IF ( Report_Flags(iRF_Time) > 0 ) THEN

    CALL Output_Time_Report

END IF

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
