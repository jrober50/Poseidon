   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE IO_Fixed_Point_Diagnostics                                                   !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!


!*D*================================!
!                                   !
!           Dependencies            !
!                                   !
!===================================!


USE Poseidon_Kinds_Module, &
            ONLY : idp

USE Poseidon_Parameters, &
            ONLY :  Cur_Iteration,          &
                    Max_Iterations,         &
                    Convergence_Criteria

USE Variables_IO, &
            ONLY :  Frame_Update_Table,     &
                    Frame_Residual_Table,   &
                    File_Suffix

USE Poseidon_IO_Parameters, &
            ONLY :  Poseidon_Reports_Dir

USE Poseidon_File_Routines_Module, &
            ONLY :  Open_New_File

USE Flags_IO_Module, &
            ONLY :  lPF_IO_Flags,               &
                    iPF_IO_Write_FP_Diagnostics
                    
USE Variables_FP, &
            ONLY :  Resid_Norms,                &
                    Update_Norms,               &
                    FP_Iteration_Log

IMPLICIT NONE




!*F&S*==========================================!
!                                               !
!           Functions & Subroutines             !
!                                               !
!===============================================!
CONTAINS



!+101+##########################################################################!
!                                                                               !
!                    Output_Convergence_Reports                                 !
!                                                                               !
!###############################################################################!
SUBROUTINE Output_FP_Diagnostics()

CHARACTER(LEN = 100)                                        ::  Iter_Name
CHARACTER(LEN = 200), DIMENSION(:), ALLOCATABLE             ::  Resid_Name
CHARACTER(LEN = 200), DIMENSION(:), ALLOCATABLE             ::  Update_Name

INTEGER                                                     ::  Iter_ID
INTEGER, DIMENSION(:), ALLOCATABLE                          ::  Resid_IDs
INTEGER, DIMENSION(:), ALLOCATABLE                          ::  Update_IDs
INTEGER                                                     ::  Num_Files

INTEGER                                                     ::  i, iU, iter

116 FORMAT (A,A,A,A)



IF ( lPF_IO_Flags(iPF_IO_Write_FP_Diagnostics) ) THEN

    Num_Files = 2
    ALLOCATE( Resid_Name(1:Num_Files) )
    ALLOCATE( Update_Name(1:Num_Files) )
    ALLOCATE( Resid_IDs(1:Num_Files) )
    ALLOCATE( Update_IDs(1:Num_Files) )

    WRITE(Iter_Name,116) Poseidon_Reports_Dir,"Iterations_",TRIM(File_Suffix),".out"

    WRITE(Resid_Name(1),116) Poseidon_Reports_Dir,"Residual_ConFactor_",TRIM(File_Suffix),".out"
    WRITE(Resid_Name(2),116) Poseidon_Reports_Dir,"Residual_Lapse_",TRIM(File_Suffix),".out"
!    WRITE(Resid_Name(3),116) Poseidon_Reports_Dir,"Residual_Beta1_",TRIM(File_Suffix),".out"
!    WRITE(Resid_Name(4),116) Poseidon_Reports_Dir,"Residual_Beta2_",TRIM(File_Suffix),".out"
!    WRITE(Resid_Name(5),116) Poseidon_Reports_Dir,"Residual_Beta3_",TRIM(File_Suffix),".out"

    WRITE(Update_Name(1),116) Poseidon_Reports_Dir,"Update_ConFactor_",TRIM(File_Suffix),".out"
    WRITE(Update_Name(2),116) Poseidon_Reports_Dir,"Update_Lapse_",TRIM(File_Suffix),".out"
!    WRITE(Update_Name(3),116) Poseidon_Reports_Dir,"Update_Beta1_",TRIM(File_Suffix),".out"
!    WRITE(Update_Name(4),116) Poseidon_Reports_Dir,"Update_Beta2_",TRIM(File_Suffix),".out"
!    WRITE(Update_Name(5),116) Poseidon_Reports_Dir,"Update_Beta3_",TRIM(File_Suffix),".out"


    Iter_ID = 920
    Resid_IDs = [(900 + i, i=1,Num_Files)]
    Update_IDs = [(910 + i, i=1,Num_Files)]

    CALL OPEN_NEW_FILE( Iter_Name, Iter_ID, Iter_ID )
    DO i = 1,Num_Files
        CALL OPEN_NEW_FILE( Resid_Name(i), Resid_IDs(i), Resid_IDs(i) )
        CALL OPEN_NEW_FILE( Update_Name(i), Update_IDs(i), Update_IDs(i) )
    END DO

    WRITE(Iter_ID,*)FP_Iteration_Log, Max_Iterations, Convergence_Criteria



    DO iU = 1,2
        DO iter = 1,Max_Iterations
            WRITE(Resid_IDs(iU),*) Resid_Norms(3,1,iter,iU)
        END DO
    END DO



    DO iU = 1,2
        DO iter = 1,Max_Iterations
            WRITE(Update_IDs(iU),*) Update_Norms(3,iter,iU)
        END DO
    END DO


    CLOSE( Unit = Iter_ID )
    DO i = 1,Num_Files
        CLOSE( Unit = Resid_IDs(i) )
        CLOSE( Unit = Update_IDs(i) )
    END DO

END IF

END SUBROUTINE  Output_FP_Diagnostics









END MODULE IO_Fixed_Point_Diagnostics
