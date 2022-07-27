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




INTEGER                     ::  ITER_REPORT_NUM_SAMPLES     = 20

INTEGER                     ::  WRITE_RESULTS_R_SAMPS       = 256
INTEGER                     ::  WRITE_RESULTS_T_SAMPS       = 256
INTEGER                     ::  WRITE_RESULTS_P_SAMPS       = 1






INTEGER, PARAMETER                  :: N_RF         = 6

INTEGER, PARAMETER                  :: iRF_Run      = 1
INTEGER, PARAMETER                  :: iRF_Frame    = 2
INTEGER, PARAMETER                  :: iRF_Iter     = 3
INTEGER, PARAMETER                  :: iRF_Setup    = 4
INTEGER, PARAMETER                  :: iRF_Time     = 5
INTEGER, PARAMETER                  :: iRF_Converge = 6

INTEGER, DIMENSION(N_RF)            :: Report_Flags = 0
INTEGER, DIMENSION(N_RF)            :: Report_IDs   = -1
CHARACTER(LEN=14), DIMENSION(N_RF)  :: Report_Names = [ 'Run           ',   &
                                                        'Frame         ',   &
                                                        'Iteration     ',   &
                                                        'Setup         ',   &
                                                        'Timetable     ',   &
                                                        'Convergence   '    ]




INTEGER, PARAMETER                  :: N_PF         = 2

INTEGER, PARAMETER                  :: iPF_Cond     = 1
INTEGER, PARAMETER                  :: iPF_Results  = 1

INTEGER, DIMENSION(N_PF)            :: Print_Flags = 2
CHARACTER(LEN=16), DIMENSION(N_PF)  :: Print_Names = [  'Condition Number', &
                                                        'Results         ']



INTEGER, PARAMETER                  :: N_WF         = 8

INTEGER, PARAMETER                  :: iWF_Matrix   = 1
INTEGER, PARAMETER                  :: iWF_RHS      = 2
INTEGER, PARAMETER                  :: iWF_Update   = 3
INTEGER, PARAMETER                  :: iWF_Source   = 4
INTEGER, PARAMETER                  :: iWF_Results  = 5
INTEGER, PARAMETER                  :: iWF_Coeffs   = 6
INTEGER, PARAMETER                  :: iWF_Mesh     = 7
INTEGER, PARAMETER                  :: iWF_Cond     = 8

INTEGER, DIMENSION(N_WF)            :: Write_Flags = 0
CHARACTER(LEN=16), DIMENSION(N_WF)  :: Write_Names = [  'Matrix(ces)     ', &
                                                        'RHS Vector(s)   ', &
                                                        'Update Vector(s)', &
                                                        'Source(s)       ', &
                                                        'Results         ', &
                                                        'Coefficients    ', &
                                                        'Mesh            ', &
                                                        'Condition Number']



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
