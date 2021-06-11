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
INTEGER, DIMENSION(N_RF)            :: Report_Flags = 0
INTEGER, DIMENSION(N_RF)            :: Report_IDs   = -1
CHARACTER(LEN=14), DIMENSION(N_RF)  :: Report_Names = [ 'Run           ',   &
                                                        'Frame         ',   &
                                                        'Iteration     ',   &
                                                        'Setup         ',   &
                                                        'Timetable     ',   &
                                                        'Convergence   '    ]




INTEGER, PARAMETER          :: N_PF                         = 5
INTEGER, DIMENSION(N_PF)    :: Print_Flags




INTEGER, PARAMETER                  :: N_WF        = 7
INTEGER, DIMENSION(N_WF)            :: Write_Flags = 0
CHARACTER(LEN=16), DIMENSION(N_WF)  :: Write_Names = [  'Matrix(ces)     ', &
                                                        'RHS Vector(s)   ', &
                                                        'Update Vector(s)', &
                                                        'Source(s)       ', &
                                                        'Results         ', &
                                                        'Coefficients    ', &
                                                        'Mesh            '  ]



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
