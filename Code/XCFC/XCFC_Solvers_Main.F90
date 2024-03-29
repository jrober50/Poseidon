   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE XCFC_Solvers_Main_Module                                                     !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   XCFC_X_Solve                                                        !##!
!##!    +201+   XCFC_ConFactor_Solve                                                !##!
!##!    +301+   XCFC_Lapse_Solve                                                    !##!
!##!    +401+   XCFC_Shift_Solve                                                    !##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!


!*D*================================!
!                                   !
!           Dependencies            !
!                                   !
!===================================!
USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag

USE Poseidon_Message_Routines_Module, &
            ONLY :  Run_Message
            
USE Poseidon_Bailout_Module, &
            ONLY :  Poseidon_Bailout

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

USE Variables_Mesh, &
            ONLY :  Num_R_Elements,             &
                    Num_T_Elements,             &
                    Num_P_Elements

USE XCFC_Fixed_Point_Module, &
            ONLY :  XCFC_Fixed_Point

USE Load_Vector_XCFC_TypeB_Module , &
            ONLY :  Create_Load_Vector_XCFC_TypeB

USE Linear_System_Solvers_TypeB_Module, &
            ONLY :  Solve_Linear_System_TypeB

USE Timer_Routines_Module, &
            ONLY :  TimerStart,                     &
                    TimerStop

USE Timer_Variables_Module, &
            ONLY :  Timer_X,                   &
                    Timer_X_SourceVector,      &
                    Timer_X_LinearSolve,       &
                    Timer_Shift,               &
                    Timer_Shift_SourceVector,  &
                    Timer_Shift_LinearSolve,   &
                    Timer_Lapse,               &
                    Timer_ConFactor

USE IO_Print_Results, &
            ONLY :  Print_Single_Var_Results,       &
                    Print_Vector_Var_Results

#ifdef POSEIDON_MEMORY_FLAG
USE Poseidon_Memory_Routines, &
            ONLY :  Poseidon_Mark_Memory

USE Memory_Variables_Module, &
            ONLY :  Memory_Method_Before_X_Load,    &
                    Memory_Method_X_Between
                    
                    
#endif


IMPLICIT NONE

CONTAINS

!+201+##########################################################################!
!                                                                               !
!                       Solve_X_System                                             !
!                                                                               !
!###############################################################################!
SUBROUTINE XCFC_X_Solve()

INTEGER, DIMENSION(3)                                   ::  iU
INTEGER                                                 ::  iVB

INTEGER, DIMENSION(3)                                   ::  iEU
INTEGER, DIMENSION(3)                                   ::  iEL



IF ( Verbose_Flag ) THEN
    WRITE(*,'()')
    CALL Run_Message('Beginning X System.')
END IF
CALL TimerStart( Timer_X )


iU = [iU_X1, iU_X2, iU_X3]
iVB = iVB_X
iEL = [0, 0, 0]
iEU = [Num_R_Elements-1,Num_T_Elements-1,Num_P_Elements-1]

!PRINT*,"Before Calc_Source"

#ifdef POSEIDON_MEMORY_FLAG
CALL Poseidon_Mark_Memory(Memory_Method_Before_X_Load)
PRINT*,"Before X Load Vector                 : ",Memory_Method_Before_X_Load
#endif



CALL TimerStart( Timer_X_SourceVector )
CALL Create_Load_Vector_XCFC_TypeB( iU, iVB, iEU, iEL )
CALL TimerStop(  Timer_X_SourceVector )


#ifdef POSEIDON_MEMORY_FLAG
CALL Poseidon_Mark_Memory(Memory_Method_X_Between)
PRINT*,"Between X Vector and Solve           : ",Memory_Method_X_Between
#endif

!PRINT*,"Before Solve"
CALL TimerStart( Timer_X_LinearSolve )
CALL Solve_Linear_System_TypeB( iU, iVB )
CALL TimerStop(  Timer_X_LinearSolve )




CALL TimerStop( Timer_X )


!CALL Print_Single_Var_Results( iU_X1, iVB_X )
!CALL Print_Vector_Var_Results( iVB_X )
!IF ( .TRUE. ) CALL Poseidon_Bailout("at the end of XCFC_x_Solve")


END SUBROUTINE XCFC_X_Solve








!+201+##########################################################################!
!                                                                               !
!                       XCFC_ConFactor_Solve                                  !
!                                                                               !
!###############################################################################!
SUBROUTINE XCFC_ConFactor_Solve()


IF ( Verbose_Flag ) THEN
    WRITE(*,'()')
    CALL Run_Message('Beginning Conformal Factor System.')
END IF

CALL TimerStart( Timer_ConFactor )
CALL XCFC_Fixed_Point(iU_CF)
CALL TimerStop( Timer_ConFactor )


!CALL Print_Single_Var_Results( iU_CF )
!IF ( .TRUE. ) CALL Poseidon_Bailout("Stopping at the end of XCFC_ConFactor_Solve")


END SUBROUTINE XCFC_ConFactor_Solve







!+301+##########################################################################!
!                                                                               !
!                       XCFC_ConFactor_Solve                                  !
!                                                                               !
!###############################################################################!
SUBROUTINE XCFC_Lapse_Solve()

IF ( Verbose_Flag ) THEN
    WRITE(*,'()')
    CALL Run_Message('Beginning Lapse Function System.')
END IF


CALL TimerStart( Timer_Lapse )

CALL XCFC_Fixed_Point(iU_LF)

CALL TimerStop( Timer_Lapse )


!CALL Print_Single_Var_Results( iU_LF )
!IF ( .TRUE. ) CALL Poseidon_Bailout("Stopping at the end of XCFC_Lapse_Solve")

END SUBROUTINE XCFC_Lapse_Solve




!+401+##########################################################################!
!                                                                               !
!                       Solve_Shift_System                                      !
!                                                                               !
!###############################################################################!
SUBROUTINE XCFC_Shift_Solve()

INTEGER, DIMENSION(3)                               ::  iU
INTEGER                                             ::  iVB

INTEGER, DIMENSION(3)                               ::  iEU
INTEGER, DIMENSION(3)                               ::  iEL


IF ( Verbose_Flag ) THEN
    WRITE(*,'()')
    CALL Run_Message('Beginning Shift System.')
END IF



CALL TimerStart( Timer_Shift )


iU = [iU_S1, iU_S2, iU_S3]
iVB = iVB_S
iEL = [0, 0, 0]
iEU = [Num_R_Elements-1,Num_T_Elements-1,Num_P_Elements-1]

CALL TimerStart( Timer_Shift_SourceVector )
CALL Create_Load_Vector_XCFC_TypeB( iU, iVB, iEU, iEL )
CALL TimerStop(  Timer_Shift_SourceVector )

CALL TimerStart( Timer_Shift_LinearSolve )
CALL Solve_Linear_System_TypeB( iU, iVB )
CALL TimerStop(  Timer_Shift_LinearSolve )


CALL TimerStop( Timer_Shift )

!CALL Print_Single_Var_Results( iU_S1, iVB_S )
!IF ( .TRUE. ) STOP "at the end of XCFC_Shift_Solve"

END SUBROUTINE XCFC_Shift_Solve




END MODULE XCFC_Solvers_Main_Module

