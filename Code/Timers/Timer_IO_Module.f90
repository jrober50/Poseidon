   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Timer_IO_Module                                                      !##!
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
            ONLY :  Poisson_Mode

USE Timer_Variables_Module


IMPLICIT NONE


CONTAINS



 !+101+########################################################!
!                                                               !
!          Output_Time_Report                         		    !
!                                                               !
 !#############################################################!
SUBROUTINE Output_Time_Report()


IF ( Poisson_Mode ) THEN
    CALL Output_Poisson_Time_Report()
ELSE
    CALL Output_XCFC_Time_Report()
END IF


END SUBROUTINE Output_Time_Report








 !+201+########################################################!
!                                                               !
!          Output_Poisson_Time_Report                           !
!                                                               !
 !#############################################################!
SUBROUTINE Output_Poisson_Time_Report()



CALL Output_Header()


WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  '- Source Input                   :', Timer_Poisson_SourceInput, ' s'
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  '   - Source Input Part A         :', Timer_Poisson_SourceInput_PartA, ' s'
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  '   - Source Input Part B         :', Timer_Poisson_SourceInput_PartB, ' s'
WRITE(*,*)



!
! Initiliazation
!   + Matrix Construction
!
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  '- Poseidon Initialization        :', Timer_Core_Initialization, ' s     '
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  '   - Matrix Construction         :', Timer_Poisson_Matrix_Init,' s     '
WRITE(*,*)




!
! Source Vector Construction
!   +SubParts
!   +Main
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  '- Source Vector Construction     :', Timer_Poisson_SourceVector, ' s'
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  '   - Source Vector Subparts      :', Timer_Poisson_SourceVector_SubParts,' s'
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  '   - Source Vector Main          :', Timer_Poisson_SourceVector_Main, ' s'
WRITE(*,*)






! Linear Solver
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  '- Linear Solve                   :', Timer_Poisson_LinearSolve, ' s'
WRITE(*,*)





WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  '- Print Results                  :', Timer_Core_PrintResults, ' s'
WRITE(*,*)
WRITE(*,*)
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  'Missing Time                     :'      &
        ,Timer_Poseidon                     &
        - Timer_Core_Initialization         &
        - Timer_Poisson_SourceVector        &
        - Timer_Poisson_LinearSolve         &
        - Timer_Core_PrintResults, ' s'



END SUBROUTINE Output_Poisson_Time_Report







 !+202+########################################################!
!                                                               !
!          Output_XCFC_Time_Report                              !
!                                                               !
 !#############################################################!
SUBROUTINE Output_XCFC_Time_Report()

CALL Output_Header()


WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  '- Source Input                   :', Timer_GR_SourceInput, ' s'
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  '   - Source Input Part A         :', Timer_GR_SourceInput_PartA, ' s'
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  '   - Source Input Part B         :', Timer_GR_SourceInput_PartB, ' s'
WRITE(*,*)



!
! Initiliazation
!   + Matrix Construction
!
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  '- Poseidon Initialization        :', Timer_Core_Initialization, ' s     '
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  '   - Matrix Construction         :', Timer_XCFC_Matrix_Init,' s     '
WRITE(*,*)


WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  '- Cholesky Factorization         :', Timer_XCFC_Matrix_Cholesky,' s     '
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  '- Type B Matrix Factorization    :', Timer_XCFC_Type_B_Factorization,' s     '
WRITE(*,*)




!
! XCFC Variable X
!   +SourceVector
!   +LinearSolve
WRITE(*,'(7X,A,5X,ES12.6E2,A)')&
  '- XCFC, X                        :', Timer_XCFC_X,' s'
WRITE(*,'(10X,A,5X,ES12.6E2,A)') &
  '- Source Vector Construction  :', Timer_XCFC_X_SourceVector,' s'
WRITE(*,'(10X,A,5X,ES12.6E2,A)') &
  '- Linear Solve                :', Timer_XCFC_X_LinearSolve-Timer_XCFC_Type_B_Factorization,' s'
WRITE(*,*)



!
! XCFC Variable Conformal Factor
!   +SourceVector
!   +LinearSolve
WRITE(*,'(7X,A,5X,ES12.6E2,A)')&
  '- XCFC, Conformal Factor         :', Timer_XCFC_ConFactor,' s'
WRITE(*,'(10X,A,5X,ES12.6E2,A)') &
  '- Source Vector Construction  :', Timer_XCFC_ConFactor_SourceVector,' s'
WRITE(*,'(10X,A,5X,ES12.6E2,A)') &
  '- Linear Solve                :', Timer_XCFC_ConFactor_LinearSolve, ' s'
WRITE(*,*)



!
! XCFC Variable Lapse Function
!   +SourceVector
!   +LinearSolve
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  '- XCFC, Lapse Function           :', Timer_XCFC_Lapse,' s'
WRITE(*,'(10X,A,5X,ES12.6E2,A)') &
  '- Source Vector Construction  :', Timer_XCFC_Lapse_SourceVector,' s'
WRITE(*,'(10X,A,5X,ES12.6E2,A)') &
  '- Linear Solve                :', Timer_XCFC_Lapse_LinearSolve, ' s'
WRITE(*,*)



!
! XCFC Variable Shift
!   +SourceVector
!   +LinearSolve
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  '- XCFC Variable, Shift Vector    :',Timer_XCFC_Shift
WRITE(*,'(10X,A,5X,ES12.6E2,A)') &
  '- Source Vector Construction  :', Timer_XCFC_Shift_SourceVector, ' s'
WRITE(*,'(10X,A,5X,ES12.6E2,A)') &
  '- Linear Solve                :', Timer_XCFC_Shift_LinearSolve, ' s'
WRITE(*,*)


CALL Output_Footer()






END SUBROUTINE Output_XCFC_Time_Report




 !+202+########################################################!
!                                                               !
!          Output_Header                                        !
!                                                               !
 !#############################################################!
SUBROUTINE Output_Header()

REAL(idp)           :: Timer_Accountable
!REAL(idp)           :: Timer_Driver_Accountable

Timer_Accountable = Timer_GR_SourceInput            &
                    + Timer_Core_Initialization     &
                    + Timer_XCFC_Matrix_Cholesky    &
                    + Timer_XCFC_X                  &
                    + Timer_XCFC_ConFactor          &
                    + Timer_XCFC_Lapse              &
                    + Timer_XCFC_Shift              &
                    + Timer_Core_Utilities          &
                    + Timer_Core_PrintResults

WRITE(*,*)
WRITE(*,'(10X,A)')'    Timers    '
WRITE(*,'(10X,A)')'--------------'
WRITE(*,*)
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  'Program Total Time               :', Timer_Total, ' s'
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  'Initialize Test Problem          :', Timer_Core_Init_Test_Problem, ' s'
WRITE(*,*)

WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  'Poseidon Total Run Time          :', Timer_Poseidon, ' s'
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  'Print Results Time               :', Timer_Core_PrintResults, ' s'
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  'Utilities Time                   :', Timer_Core_Utilities, ' s'
WRITE(*,*)


!WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
!  'Driver - Set Source              :', Timer_Driver_SetSource, ' s'
!WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
!  '  Driver - Set Source - InitTest :', Timer_Driver_SetSource_InitTest, ' s'
!WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
!  '  Driver - Set Source - Scale    :', Timer_Driver_SetSource_Scale, ' s'
!WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
!  '  Driver - Set Source - Input    :', Timer_Driver_SetSource_SetSource, ' s'
!WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
!  'Driver - Set BC                  :', Timer_Driver_SetBC, ' s'
!WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
!  'Driver - Set Guess               :', Timer_Driver_SetGuess, ' s'
!WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
!  'Driver - Run                     :', Timer_Driver_Run, ' s'
!WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
!  'Driver - Extra                   :', Timer_Driver_Extra, ' s'
!
!Timer_Driver_Accountable = Timer_Driver_SetSource    &
!                        + Timer_Driver_SetBC        &
!                        + Timer_Driver_SetGuess     &
!                        + Timer_Driver_Run          &
!                        + Timer_Driver_Extra
!
!WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
!'Driver - Total                   :', Timer_Driver_Accountable, ' s'
!WRITE(*,*)



END SUBROUTINE Output_Header




 !+202+########################################################!
!                                                               !
!          Output_Header                                        !
!                                                               !
 !#############################################################!
SUBROUTINE Output_Footer()

REAL(idp)           :: Timer_Accountable

Timer_Accountable = Timer_GR_SourceInput            &
                    + Timer_Core_Initialization     &
                    + Timer_XCFC_Matrix_Cholesky    &
                    + Timer_XCFC_X                  &
                    + Timer_XCFC_ConFactor          &
                    + Timer_XCFC_Lapse              &
                    + Timer_XCFC_Shift              &
                    + Timer_Core_Utilities          &
                    + Timer_Core_PrintResults

WRITE(*,*)
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  'Accountable Time                 :', Timer_Accountable, ' s'
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  'Unaccounted Time                 :', Timer_Poseidon - Timer_Accountable, ' s'
WRITE(*,*)





END SUBROUTINE Output_Footer









END MODULE Timer_IO_Module
