   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE IO_Frame_Report                                                              !##!
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


USE Poseidon_Kinds_Module, &
                    ONLY : idp

USE Poseidon_Numbers_Module, &
                    ONLY : pi

USE Units_Module, &
                    ONLY :  C_Square,       &
                            Gram,           &
                            Centimeter,     &
                            Kilometer,      &
                            Erg,            &
                            Second,         &
                            GravPot_Units,  &
                            Shift_Units


USE Driver_Parameters,  &
                    ONLY  : DRIVER_TEST_NUMBER,                             &
                            DRIVER_FRAME,                                   &
                            DRIVER_INNER_RADIUS,                            &
                            DRIVER_OUTER_RADIUS,                            &
                            DRIVER_TOTAL_FRAMES,                            &
                            myID,                                           &
                            Driver_R_Input_Nodes

USE Poseidon_Parameters, &
                    ONLY :  DEGREE,                 &
                            L_LIMIT,                &
                            Domain_Dim,             &
                            CUR_ITERATION,          &
                            Poseidon_Frame,         &
                            Convergence_Flag,       &
                            Max_Iterations


USE Variables_Mesh, &
                    ONLY :  NUM_R_ELEMENTS,         &
                            rlocs,                  &
                            R_Inner,                &
                            R_Outer

USE Variables_IO, &
                    ONLY :  Iter_Report_File_ID,    &
                            Frame_Report_File_ID,   &
                            Run_Report_File_ID,     &
                            Write_Report_Flag,      &
                            Write_Results_Flag,     &
                            Write_Timetable_Flag,   &
                            Write_Sources_Flag,     &
                            Run_Report_Flag,        &
                            Results_Output_Flag,    &
                            Iter_Time_Table,        &
                            Frame_Time_Table,       &
                            Run_Time_Table,         &
                            Write_Results_R_Samps,  &
                            Write_Results_T_Samps,  &
                            Write_Results_P_Samps,  &
                            Total_Run_Iters,        &
                            Iter_Report_Num_Samples,&
                            Iter_Time_Table,        &
                            Iteration_Histogram


USE Variables_Yahil, &
                    ONLY :  SelfSim_T

USE Variables_Functions, &
                    ONLY :  Potential_Solution,     &
                            Shift_Solution

USE Functions_Results, &
                    ONLY :  Calc_1D_CFA_Values

USE Functions_Quadrature, &
                    ONLY :  Initialize_LG_Quadrature_Locations


USE Functions_Mesh, &
                    ONLY :  Create_Logarithmic_1D_Mesh

USE Poseidon_IO_Parameters, &
                    ONLY :  Poseidon_Reports_Dir,                           &
                            Poseidon_IterReports_Dir,                       &
                            Poseidon_Objects_Dir,                           &
                            Poseidon_LinSys_Dir,                            &
                            Poseidon_Results_Dir,                           &
                            Poseidon_Sources_Dir

USE Functions_Info,   &
                    ONLY  : PQ_ITERATIONS_MAX,  &
                            PQ_ITERATIONS_HIST, &
                            PQ_TIMETABLE_FRAME, &
                            PQ_TIMETABLE_RUN


IMPLICIT NONE


CHARACTER(LEN = 20), PARAMETER    :: Filename_Format_A = "(A,A)"
CHARACTER(LEN = 20), PARAMETER    :: Filename_Format_B = "(A,A,I2.2,A,I2.2,A)"


!*F&S*==========================================!
!                                               !
!           Functions & Subroutines             !
!                                               !
!===============================================!
CONTAINS

