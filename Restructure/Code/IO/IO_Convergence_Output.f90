   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE IO_Convergence_Output                                                        !##!
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

USE Units_Module, &
                    ONLY :  Centimeter,             &
                            Shift_Units

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

USE Poseidon_IO_Module, &
                    ONLY :  Open_New_File


IMPLICIT NONE




!*F&S*==========================================!
!                                               !
!           Functions & Subroutines             !
!                                               !
!===============================================!
CONTAINS



!+101+##########################################################################!
!                                                                               !
!                    Output_Convergence_Data                                    !
!                                                                               !
!###############################################################################!
SUBROUTINE Output_Convergence_Data()

CHARACTER(LEN = 100)                                        ::  Iter_Name
CHARACTER(LEN = 100), DIMENSION(:), ALLOCATABLE             ::  Resid_Name
CHARACTER(LEN = 100), DIMENSION(:), ALLOCATABLE             ::  Update_Name

INTEGER                                                     ::  Iter_ID
INTEGER, DIMENSION(:), ALLOCATABLE                          ::  Resid_IDs
INTEGER, DIMENSION(:), ALLOCATABLE                          ::  Update_IDs
INTEGER                                                     ::  Num_Files

INTEGER                                                     ::  i, ui, iter

116 FORMAT (A,A,A,A)

Num_Files = 5
ALLOCATE( Resid_Name(1:Num_Files) )
ALLOCATE( Update_Name(1:Num_Files) )
ALLOCATE( Resid_IDs(1:Num_Files) )
ALLOCATE( Update_IDs(1:Num_Files) )

WRITE(Iter_Name,116) Poseidon_Reports_Dir,"Iterations_",TRIM(File_Suffix),".out"

WRITE(Resid_Name(1),116) Poseidon_Reports_Dir,"Residual_ConFactor_",TRIM(File_Suffix),".out"
WRITE(Resid_Name(2),116) Poseidon_Reports_Dir,"Residual_Lapse_",TRIM(File_Suffix),".out"
WRITE(Resid_Name(3),116) Poseidon_Reports_Dir,"Residual_Beta1_",TRIM(File_Suffix),".out"
WRITE(Resid_Name(4),116) Poseidon_Reports_Dir,"Residual_Beta2_",TRIM(File_Suffix),".out"
WRITE(Resid_Name(5),116) Poseidon_Reports_Dir,"Residual_Beta3_",TRIM(File_Suffix),".out"

WRITE(Update_Name(1),116) Poseidon_Reports_Dir,"Update_ConFactor_",TRIM(File_Suffix),".out"
WRITE(Update_Name(2),116) Poseidon_Reports_Dir,"Update_Lapse_",TRIM(File_Suffix),".out"
WRITE(Update_Name(3),116) Poseidon_Reports_Dir,"Update_Beta1_",TRIM(File_Suffix),".out"
WRITE(Update_Name(4),116) Poseidon_Reports_Dir,"Update_Beta2_",TRIM(File_Suffix),".out"
WRITE(Update_Name(5),116) Poseidon_Reports_Dir,"Update_Beta3_",TRIM(File_Suffix),".out"

Iter_ID = 920
Resid_IDs = [(900 + i, i=1,Num_Files)]
Update_IDs = [(910 + i, i=1,Num_Files)]

CALL OPEN_NEW_FILE( Iter_Name, Iter_ID, Iter_ID )
DO i = 1,Num_Files
    CALL OPEN_NEW_FILE( Resid_Name(i), Resid_IDs(i), Resid_IDs(i) )
    CALL OPEN_NEW_FILE( Update_Name(i), Update_IDs(i), Update_IDs(i) )
END DO

WRITE(Iter_ID,*)Cur_Iteration, Max_Iterations, Convergence_Criteria

DO ui = 1,5
    DO iter = 1,Cur_Iteration
        WRITE(Resid_IDs(ui),*) Frame_Residual_Table(:,iter,ui)
    END DO
END DO

DO ui = 1,5
    DO iter = 1,Cur_Iteration
        WRITE(Update_IDs(ui),*) Frame_Update_Table(iter,ui)
    END DO
END DO


CLOSE( Unit = Iter_ID )
DO i = 1,Num_Files
    CLOSE( Unit = Resid_IDs(i) )
    CLOSE( Unit = Update_IDs(i) )
END DO



END SUBROUTINE Output_Convergence_Data









END MODULE IO_Convergence_Output
