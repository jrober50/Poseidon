   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE FP_IO_Module                                                                 !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
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

USE Poseidon_File_Routines_Module, &
        ONLY :  Open_New_File


USE Poseidon_IO_Parameters, &
        ONLY :  Poseidon_Reports_Dir

USE Variables_IO, &
        ONLY :  Run_Time_Table,         &
                File_Suffix,            &
                Report_Flags,           &
                iRF_Time

USE Timer_IO_Module, &
        ONLY :  Output_Time_Report

USE Variables_FP, &
        ONLY :  FPTT_Names

!*F&S*==========================================!
!                                               !
!           Functions & Subroutines             !
!                                               !
!===============================================!
CONTAINS


!+101+###########################################################################!
!                                                                                !
!           Output_FP_TimeTable                                                  !
!                                                                                !
!################################################################################!
SUBROUTINE Output_FP_Timetable()

INTEGER                                             ::  File_ID

CHARACTER(LEN = 100)                                ::  File_Name



IF ( Report_Flags(iRF_Time) == 1 ) THEN
    CALL Output_Time_Report()
!    WRITE(*,'(A)')'------------- Fixed Point Timetable -----------'
!    WRITE(*,'(A,A,ES12.4E3,A)')FPTT_Names(1), ' : ',Run_Time_Table(1), ' s'
!    WRITE(*,'(A,A,ES12.4E3,A)')FPTT_Names(3), ' : ',Run_Time_Table(3), ' s'
!    WRITE(*,'(A,A,ES12.4E3,A)')FPTT_Names(11),' : ',Run_Time_Table(11),' s'
!    WRITE(*,'(A,A,ES12.4E3,A)')FPTT_Names(4), ' : ',Run_Time_Table(4), ' s'
!    WRITE(*,'(A,A,ES12.4E3,A)')FPTT_Names(12),' : ',Run_Time_Table(12),' s'
!    WRITE(*,'(A,A,ES12.4E3,A)')FPTT_Names(5), ' : ',Run_Time_Table(5), ' s'
!    WRITE(*,'(A,A,ES12.4E3,A)')FPTT_Names(8), ' : ',Run_Time_Table(8), ' s'
!    WRITE(*,'(/,/)')
END IF



IF ( Report_Flags(5) > 1 ) THEN
    WRITE(File_Name, '(A,A,A,A)') Poseidon_Reports_Dir,"FP_Timetable_",trim(File_Suffix),".out"

    CALL OPEN_NEW_FILE( trim(File_Name), File_ID, 50)

    WRITE( File_ID,*) Run_Time_Table

    CLOSE( UNIT = File_ID)

END IF

END SUBROUTINE Output_FP_Timetable








END MODULE FP_IO_Module
