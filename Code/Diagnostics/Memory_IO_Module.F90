   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Memory_IO_Module                                                      !##!
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

USE Memory_Variables_Module

IMPLICIT NONE


CONTAINS



 !+201+########################################################!
!                                                               !
!          Output_Poisson_Time_Report                           !
!                                                               !
 !#############################################################!
SUBROUTINE Output_Poseidon_Memory_Loop_Report( Loop )

INTEGER,        INTENT(IN)                  ::  Loop

CHARACTER(LEN = 100)                        ::  Report_Name
INTEGER                                     ::  Suggested_Number = 475
INTEGER                                     ::  File_ID


#ifdef POSEIDON_MEMORY_FLAG

WRITE(Report_Name,'(A,A,A,A,I3.3,A)') Poseidon_Reports_Dir,            &
                                    "Memory_Report_",               &
                                    trim(File_Suffix),              &
                                    "_Loop_",                       &
                                    Loop,                  &
                                    ".out"

                            

CALL Open_New_File( Report_Name, File_ID, Suggested_Number)

WRITE(File_ID,'(I16)') Memory_Loop_Start
WRITE(File_ID,'(I16)') Memory_Loop_Before_Init
WRITE(File_ID,'(I16)') Memory_Loop_After_Init
WRITE(File_ID,'(I16)') Memory_Loop_Before_Run

WRITE(File_ID,'(I16)') Memory_Method_Start
WRITE(File_ID,'(I16)') Memory_Method_Before_CF

WRITE(File_ID,'(I16)') Memory_Method_Before_CF_LoadVector
WRITE(File_ID,'(I16)') Memory_Method_After_CF_LoadVector
WRITE(File_ID,'(I16)') Memory_Method_After_CF_FixedPoint
WRITE(File_ID,'(I16)') Memory_Method_After_CF_DeallocWork

WRITE(File_ID,'(I16)') Memory_Method_Before_LF

WRITE(File_ID,'(I16)') Memory_Method_Before_LF_LoadVector
WRITE(File_ID,'(I16)') Memory_Method_After_LF_LoadVector
WRITE(File_ID,'(I16)') Memory_Method_After_LF_FixedPoint
WRITE(File_ID,'(I16)') Memory_Method_After_LF_DeallocWork


WRITE(File_ID,'(I16)') Memory_Method_Before_SV
WRITE(File_ID,'(I16)') Memory_Method_End

WRITE(File_ID,'(I16)') Memory_Loop_After_Run
WRITE(File_ID,'(I16)') Memory_Loop_Before_Close
WRITE(File_ID,'(I16)') Memory_Loop_End

CLOSE( File_ID)


#endif

END SUBROUTINE Output_Poseidon_Memory_Loop_Report



 !+201+########################################################!
!                                                               !
!          Output_Poisson_Time_Report                           !
!                                                               !
 !#############################################################!
SUBROUTINE Output_Poseidon_Memory_Total_Report( )


CHARACTER(LEN = 100)                        ::  Report_Name
INTEGER                                     ::  Suggested_Number = 475
INTEGER                                     ::  File_ID


#ifdef POSEIDON_MEMORY_FLAG

WRITE(Report_Name,'(A,A,A,A,A,A)') Poseidon_Reports_Dir,            &
                                    "Memory_Report_",               &
                                    trim(File_Suffix),              &
                                    "_Total.out"
CALL Open_New_File( Report_Name, File_ID, Suggested_Number)

WRITE(File_ID,'(I16)') Memory_Start
WRITE(File_ID,'(I16)') Memory_After_MPI_Init
WRITE(File_ID,'(I16)') Memory_After_AMReX_Init
WRITE(File_ID,'(I16)') Memory_After_AMReX_Finalize
WRITE(File_ID,'(I16)') Memory_End


CLOSE( File_ID)

#endif

END SUBROUTINE Output_Poseidon_Memory_Total_Report




END MODULE Memory_IO_Module


