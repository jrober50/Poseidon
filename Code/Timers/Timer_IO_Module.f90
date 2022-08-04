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

USE Poseidon_IO_Parameters, &
            ONLY :  Method_Names,           &
                    Poseidon_Reports_Dir

USE Variables_IO, &
            ONLY :  File_Suffix

USE Poseidon_File_Routines_Module, &
            ONLY :  Open_New_File

USE Timer_Variables_Module

USE Flags_IO_Module, &
            ONLY :  lPF_IO_Flags,               &
                    iPF_IO_Print_TimeTable,     &
                    iPF_IO_Write_TimeTable

USE Flags_Core_Module, &
            ONLY :  iPF_Core_Flags,         &
                    iPF_Core_Method_Mode,   &
                    iPF_Core_Method_Newtonian

IMPLICIT NONE


CONTAINS



 !+101+########################################################!
!                                                               !
!          Output_Time_Report                         		    !
!                                                               !
 !#############################################################!
SUBROUTINE Output_Time_Report()


IF ( iPF_Core_Flags(iPF_Core_Method_Mode) == iPF_Core_Method_Newtonian ) THEN
    CALL Output_Poisson_Time_Report()
ELSE
    CALL Output_XCFC_Time_Report()
    CALL Write_XCFC_Time_Report()
END IF


END SUBROUTINE Output_Time_Report








 !+201+########################################################!
!                                                               !
!          Output_Poisson_Time_Report                           !
!                                                               !
 !#############################################################!
SUBROUTINE Output_Poisson_Time_Report()


IF ( lPF_IO_Flags(iPF_IO_Print_TimeTable) )  THEN

    CALL Output_Header()

    101 FORMAT (7X,A,5X,ES12.6E2,A)



    WRITE(*,101) '- Source Input                    :', Timer_Poisson_SourceInput, ' s'
    WRITE(*,101) '   - Source Input Part A          :', Timer_Poisson_SourceInput_PartA, ' s'
    WRITE(*,101) '   - Source Input Part B          :', Timer_Poisson_SourceInput_PartB, ' s'
    WRITE(*,*)




    WRITE(*,101) '- Poseidon Initialization, Total   :', Timer_Initialization_Total, ' s     '
    WRITE(*,101) '   - Poseidon Initialization, Core :', Timer_Initialization_Core, ' s     '
    WRITE(*,*)



    WRITE(*,101) '- Matrix Construction              :', Timer_Matrix_Init,' s     '
    WRITE(*,*)



    !
    ! Source Vector Construction
    !   +SubParts
    !   +Main
    WRITE(*,101) '- Source Vector Construction      :', Timer_Poisson_SourceVector, ' s'
    WRITE(*,101) '   - Source Vector Subparts       :', Timer_Poisson_SourceVector_SubParts,' s'
    WRITE(*,101) '   - Source Vector Main           :', Timer_Poisson_SourceVector_Main, ' s'
    WRITE(*,*)






    ! Linear Solver
    WRITE(*,101) '- Linear Solve                    :', Timer_Poisson_LinearSolve, ' s'
    WRITE(*,*)





    WRITE(*,101) '- Print Results                   :', Timer_Core_PrintResults, ' s'
    WRITE(*,*)
    WRITE(*,*)
END IF


END SUBROUTINE Output_Poisson_Time_Report







 !+202+########################################################!
!                                                               !
!          Output_XCFC_Time_Report                              !
!                                                               !
 !#############################################################!
SUBROUTINE Output_XCFC_Time_Report()

IF ( lPF_IO_Flags(iPF_IO_Print_TimeTable) )  THEN

    CALL Output_Header()


    101 FORMAT (7X,A,5X,ES12.6E2,A)


    WRITE(*,101) '- Source Input                     :', Timer_GR_SourceInput, ' s'
    !WRITE(*,101) '   - Source Input Part A           :', Timer_GR_SourceInput_PartA, ' s'
    !WRITE(*,101) '   - Source Input Part B           :', Timer_GR_SourceInput_PartB, ' s'
    WRITE(*,*)


    WRITE(*,101) '- Poseidon Initialization, Total   :', Timer_Initialization_Total, ' s     '
    WRITE(*,101) '   - Poseidon Initialization, Core :', Timer_Initialization_Core, ' s     '
    WRITE(*,101) '   - Poseidon Initialization, XCFC :', Timer_Initialization_XCFC, ' s     '
    WRITE(*,*)



    WRITE(*,101) '- Matrix Construction              :', Timer_Matrix_Init,' s     '
    WRITE(*,*)

    WRITE(*,101) '- Matrix Factorization, Total      :', Timer_Matrix_Factorization, ' s     '
    WRITE(*,101) '   - Cholesky Factorization        :', Timer_Matrix_Cholesky,' s     '
    WRITE(*,101) '   - Banded Matrix Factorization   :', Timer_Banded_Factorization,' s     '
    WRITE(*,*)



    WRITE(*,101) '- XCFC, X                          :', Timer_X,' s'
    WRITE(*,101) '   - Source Vector Construction    :', Timer_X_SourceVector,' s'
    WRITE(*,101) '   - Linear Solve                  :', Timer_X_LinearSolve,' s'
    WRITE(*,*)


    WRITE(*,101) '- XCFC, Conformal Factor           :', Timer_ConFactor,' s'
    WRITE(*,101) '   - Source Vector Construction    :', Timer_ConFactor_SourceVector,' s'
    WRITE(*,101) '   - Linear Solve                  :', Timer_ConFactor_LinearSolve, ' s'
    WRITE(*,*)



    WRITE(*,101) '- XCFC, Lapse Function             :', Timer_Lapse,' s'
    WRITE(*,101) '   - Source Vector Construction    :', Timer_Lapse_SourceVector,' s'
    WRITE(*,101) '   - Linear Solve                  :', Timer_Lapse_LinearSolve, ' s'
    WRITE(*,*)


    WRITE(*,101) '- XCFC Variable, Shift Vector      :', Timer_Shift
    WRITE(*,101) '   - Source Vector Construction    :', Timer_Shift_SourceVector, ' s'
    WRITE(*,101) '   - Linear Solve                  :', Timer_Shift_LinearSolve, ' s'
    WRITE(*,*)


    CALL Output_Footer()

END IF




END SUBROUTINE Output_XCFC_Time_Report




 !+202+########################################################!
!                                                               !
!          Output_Header                                        !
!                                                               !
 !#############################################################!
SUBROUTINE Output_Header()



101 FORMAT (7X,A,5X,ES12.6E2,A)



WRITE(*,*)
WRITE(*,'(10X,A)')'    Timers    '
WRITE(*,'(10X,A)')'--------------'
WRITE(*,*)
WRITE(*,101) '- Program Total Time               :', Timer_Total, ' s'
WRITE(*,101) '   - Initialize Test Problem       :', Timer_Core_Init_Test_Problem, ' s'
WRITE(*,101) '   - Poseidon Total Run Time       :', Timer_Poseidon, ' s'
WRITE(*,*)





END SUBROUTINE Output_Header




 !+202+########################################################!
!                                                               !
!          Output_Header                                        !
!                                                               !
 !#############################################################!
SUBROUTINE Output_Footer()

REAL(idp)           :: Timer_Accountable
REAL(idp)           :: Timer_One_Time_Costs
REAL(idp)           :: Timer_Multi_Time_Costs

Timer_Accountable = Timer_GR_SourceInput                &
                    + Timer_Initialization_Total        &
                    + Timer_Matrix_Factorization   &
                    + Timer_X                      &
                    + Timer_ConFactor              &
                    + Timer_Lapse                  &
                    + Timer_Shift                  &
                    + Timer_Core_Utilities              &
                    + Timer_Core_PrintResults


Timer_One_Time_Costs = Timer_Initialization_Total       &
                     + Timer_Matrix_Factorization  &
                     + Timer_Core_Utilities             &
                     + Timer_Core_PrintResults

Timer_Multi_Time_Costs = Timer_X                      &
                        + Timer_ConFactor              &
                        + Timer_Lapse                  &
                        + Timer_Shift

101 FORMAT (7X,A,5X,ES12.6E2,A)


WRITE(*,101) '- Print Results Time               :', Timer_Core_PrintResults, ' s'
WRITE(*,101) '- Utilities Time                   :', Timer_Core_Utilities, ' s'
WRITE(*,*)


WRITE(*,101) '- Accountable Time                 :', Timer_Accountable, ' s'
WRITE(*,101) '- Unaccounted Time                 :', Timer_Poseidon - Timer_Accountable, ' s'
WRITE(*,*)

WRITE(*,101) '- Single Call Time                 :', Timer_One_Time_Costs,' s'
WRITE(*,101) '- Multiple Call Time               :', Timer_Multi_Time_Costs,' s'
WRITE(*,*)





END SUBROUTINE Output_Footer















 !+202+########################################################!
!                                                               !
!          Write_XCFC_Time_Report                              !
!                                                               !
 !#############################################################!
SUBROUTINE Write_XCFC_Time_Report()

CHARACTER(LEN = 100)                                ::  Report_Name
INTEGER                                             ::  Suggested_Number = 400
INTEGER                                             ::  File_ID


IF ( lPF_IO_Flags(iPF_IO_Write_TimeTable) )  THEN

    WRITE(Report_Name,'(A,A,A,A)') Poseidon_Reports_Dir,"Timer_Report_",trim(File_Suffix),".out"
    CALL Open_New_File( Report_Name, File_ID, Suggested_Number)

    WRITE(File_ID,'(ES12.6E2)') Timer_Poseidon
    WRITE(File_ID,'(ES12.6E2)') Timer_Initialization_Total
    WRITE(File_ID,'(ES12.6E2)') Timer_Matrix_Init
    WRITE(File_ID,'(ES12.6E2)') Timer_Matrix_Cholesky
    WRITE(File_ID,'(ES12.6E2)') Timer_Banded_Factorization
    WRITE(File_ID,'(ES12.6E2)') Timer_X
    WRITE(File_ID,'(ES12.6E2)') Timer_X_SourceVector
    WRITE(File_ID,'(ES12.6E2)') Timer_X_LinearSolve
    WRITE(File_ID,'(ES12.6E2)') Timer_ConFactor
    WRITE(File_ID,'(ES12.6E2)') Timer_ConFactor_SourceVector
    WRITE(File_ID,'(ES12.6E2)') Timer_ConFactor_LinearSolve
    WRITE(File_ID,'(ES12.6E2)') Timer_Lapse
    WRITE(File_ID,'(ES12.6E2)') Timer_Lapse_SourceVector
    WRITE(File_ID,'(ES12.6E2)') Timer_Lapse_LinearSolve
    WRITE(File_ID,'(ES12.6E2)') Timer_Shift
    WRITE(File_ID,'(ES12.6E2)') Timer_Shift_SourceVector
    WRITE(File_ID,'(ES12.6E2)') Timer_Shift_LinearSolve

    CLOSE(File_ID)

END IF

END SUBROUTINE Write_XCFC_Time_Report





END MODULE Timer_IO_Module
