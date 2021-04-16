   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Variables_IO                                                                 !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains parameters used to define the running of Poseidon.                 !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!

USE Poseidon_Kinds_Module, &
                ONLY : idp


IMPLICIT NONE



INTEGER                     ::  OUTPUT_SETUP_TABLE_FLAG     = 1


INTEGER                     ::  WRITE_TIMETABLE_FLAG        = 0     ! 0=Off,1 = Screen, 2 = File, 3 = Both
INTEGER                     ::  WRITE_REPORT_FLAG           = 1
INTEGER                     ::  WRITE_RESULTS_FLAG          = 1


INTEGER                     ::  ITER_REPORT_NUM_SAMPLES     = 20

INTEGER                     ::  WRITE_RESULTS_R_SAMPS       = 256
INTEGER                     ::  WRITE_RESULTS_T_SAMPS       = 256
INTEGER                     ::  WRITE_RESULTS_P_SAMPS       = 1

INTEGER                     ::  WRITE_SOURCES_FLAG          = 0


INTEGER                     ::  SOURCE_OUTPUT_FLAG
INTEGER                     ::  RESULTS_OUTPUT_FLAG
INTEGER                     ::  RUN_REPORT_FLAG
INTEGER                     ::  FRAME_REPORT_FLAG


! Report File IDs
INTEGER                     ::  RUN_REPORT_FILE_ID          = -1
INTEGER                     ::  ITER_REPORT_FILE_ID         = -1
INTEGER                     ::  FRAME_REPORT_FILE_ID        = -1


! Output Linear System to File Flags,       0 = Off, 1 = On
INTEGER                     ::  OUTPUT_MATRIX_FLAG          = 0
INTEGER                     ::  OUTPUT_RHS_VECTOR_FLAG      = 0
INTEGER                     ::  OUTPUT_UPDATE_VECTOR_FLAG   = 0



INTEGER                                                     ::  Num_Timer_Calls=25
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Iter_Time_Table
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Frame_Time_Table
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Run_Time_Table

INTEGER                                                     ::  Total_Run_Iters=1

REAL(KIND = idp), DIMENSION(:,:), ALLOCATABLE               ::  Frame_Update_Table
REAL(KIND = idp), DIMENSION(:,:,:), ALLOCATABLE             ::  Frame_Residual_Table

INTEGER, DIMENSION(:), ALLOCATABLE                          ::  Iteration_Histogram


CHARACTER(LEN=40)                                           ::  File_Suffix

END MODULE Variables_IO
