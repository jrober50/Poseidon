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
!##!    +101+   OPEN_ITER_REPORT_FILE                                               !##!
!##!    +102+   OUTPUT_ITER_TIMETABLE                                               !##!
!##!    +103+   CLOSE_ITER_REPORT_FILE                                              !##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!    +301+   OPEN_FRAME_REPORT_FILE                                              !##!
!##!    +303+   CLOSE_FRAME_REPORT_FILE                                             !##!
!##!                                                                                !##!
!##!    +401+   OUTPUT_FINAL_RESULTS                                                !##!
!##!    +402+   OUTPUT_PETSC_REPORT                                                 !##!
!##!    +403+   OUTPUT_STF_ELEM_BLOCK_MATRIX                                        !##!
!##!                                                                                !##!
!##!    +501+   CLOCK_IN                                                            !##!
!##!    +502+   OPEN_NEW_FILE                                                       !##!
!##!                                                                                !##!
!##!    +601+   OUTPUT_POSEIDON_SOURCES                                             !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!

!*D*================================!
!                                   !
!           Dependencies            !
!                                   !
!===================================!

USE Poseidon_IO_Module, &
        ONLY :  Open_New_File


USE Poseidon_IO_Parameters, &
        ONLY :  Poseidon_Reports_Dir

USE Variables_IO, &
        ONLY :  Run_Time_Table,         &
                File_Suffix

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


WRITE(File_Name, '(A,A,A,A)') Poseidon_Reports_Dir,"FP_Timetable_",trim(File_Suffix),".out"

CALL OPEN_NEW_FILE( trim(File_Name), File_ID, 50)

WRITE( File_ID,*) Run_Time_Table

CLOSE( UNIT = File_ID)



END SUBROUTINE Output_FP_Timetable








END MODULE FP_IO_Module
