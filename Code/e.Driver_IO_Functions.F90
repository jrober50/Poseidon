   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Driver_IO_Functions_Module                                                   !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +201+   OPEN_RUN_REPORT_FILE                                                !##!
!##!    +202+   OUTPUT_RUN_REPORT                                                   !##!
!##!    +203+   CLOSE_RUN_REPORT_FILE                                               !##!
!##!                                                                                !##!
!##!    +301+   OPEN_FRAME_REPORT_FILE                                              !##!
!##!    +302+   OUTPUT_FRAME_REPORT                                                 !##!
!##!    +303+   CLOSE_FRAME_REPORT_FILE                                             !##!
!##!                                                                                !##!
!##!    +502+   OPEN_NEW_FILE                                                       !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!
USE Poseidon_Constants_Module, &
        ONLY :  idp, pi, fdp

USE Units_Module, &
        ONLY :  C_Square,       &
                Gram,           &
                Centimeter,     &
                Meter,          &
                Kilometer,      &
                Second,         &
                GravPot_Units,  &
                Shift_Units

USE Driver_Parameters,  &
        ONLY  : Iteration_History,          &
                SOURCE_OUTPUT_FLAG,                             &
                RESULTS_OUTPUT_FLAG,                            &
                RUN_REPORT_FLAG,                                &
                FRAME_REPORT_FLAG,                              &
                DRIVER_TEST_NUMBER,                             &
                DRIVER_FRAME,                                   &
                DRIVER_INNER_RADIUS,                            &
                DRIVER_OUTER_RADIUS,                            &
                DRIVER_TOTAL_FRAMES,                            &
                Potential_Solution,                              &
                myID,                                           &
                Iteration_History

USE Poseidon_IO_Parameters, &
        ONLY :  Poseidon_Sources_Dir

USE Poseidon_Parameters, &
        ONLY :  Poseidon_Frame

USE Poseidon_Info_Module,   &
        ONLY  : PQ_ITERATIONS_MAX,  &
                PQ_ITERATIONS_HIST, &
                PQ_TIMETABLE_FRAME, &
                PQ_TIMETABLE_RUN

USE Driver_Additional_Functions_Module, &
        ONLY  : Calc_3D_Values_At_Location


USE Poseidon_Parameters, &
                        ONLY :  DOMAIN_DIM,                 &
                                DEGREE,                     &
                                L_LIMIT,                    &
                                NUM_SHELLS,                 &
                                NUM_SUBSHELLS,              &
                                nPROCS_POSEIDON,            &
                                STF_MAPPING_FLAG,           &
                                NUM_R_ELEMS_PER_BLOCK,      &
                                NUM_R_ELEMS_PER_SUBSHELL,   &
                                CUR_ITERATION,              &
                                MAX_ITERATIONS,             &
                                CONVERGENCE_CRITERIA,       &
                                CONVERGENCE_FLAG,           &
                                WRITE_REPORT_FLAG,          &
                                OUTPUT_MATRIX_FLAG,         &
                                OUTPUT_RHS_VECTOR_FLAG,     &
                                ITER_REPORT_NUM_SAMPLES,    &
                                ITER_REPORT_FILE_ID,        &
                                FRAME_REPORT_FILE_ID,       &
                                RUN_REPORT_FILE_ID

USE Poseidon_Variables_Module, &
                        ONLY :  NUM_R_ELEMENTS,             &
                                NUM_T_ELEMENTS,             &
                                NUM_P_ELEMENTS,             &
                                rlocs, tlocs, plocs,        &
                                NUM_R_NODES,                &
                                Coefficient_Vector,         &
                                Block_RHS_Vector,           &
                                myID_Poseidon,              &
                                myID_PETSC,                 &
                                myID_Shell,                 &
                                myShell,                    &
                                BLOCK_ELEM_STF_MATVEC,      &
                                POSEIDON_COMM_WORLD,        &
                                POSEIDON_COMM_PETSC,        &
                                POSEIDON_COMM_SHELL,        &
                                PROB_DIM,                   &
                                Block_Prob_Dim,             &
                                SUBSHELL_PROB_DIM,          &
                                Local_Length,               &
                                ULM_LENGTH,                 &
                                LM_LENGTH,                  &
                                Update_Vector,              &
                                NUM_OFF_DIAGONALS,          &
                                R_INNER, R_OUTER,           &
                                Total_Run_Iters,            &
                                ITER_TIME_TABLE,            &
                                FRAME_CONVERGENCE_TABLE,    &
                                Iteration_Histogram,        &
                                Matrix_Location

USE Mesh_Module, &
                        ONLY  : Create_Logarithmic_1D_Mesh

USE Poseidon_IO_Parameters, &
                    ONLY :  Poseidon_Reports_Dir,                           &
                            Poseidon_IterReports_Dir,                       &
                            Poseidon_Objects_Dir,                           &
                            Poseidon_LinSys_Dir,                            &
                            Poseidon_Results_Dir,                           &
                            Poseidon_Sources_Dir

IMPLICIT NONE

CONTAINS

 !+201+############################################################################!
!                                                                                   !
!                     OPEN_RUN_REPORT_FILE                                          !
!                                                                                   !
 !#################################################################################!
SUBROUTINE OPEN_RUN_REPORT_FILE()

CHARACTER(LEN = 39)                                     ::  FILE_NAME
INTEGER                                                 ::  istat
LOGICAL                                                 ::  FLAG, OK


IF ( RUN_REPORT_FLAG == 1 ) THEN
    WRITE(FILE_NAME,'(A,A)') Poseidon_Reports_Dir,"Run_Report.out"
    CALL OPEN_NEW_FILE( FILE_NAME, RUN_REPORT_FILE_ID )
END IF


END SUBROUTINE OPEN_RUN_REPORT_FILE


!+202+##########################################################################!
!                                                                               !
!                   OUTPUT_RUN_REPORT                                           !
!                                                                               !
!###############################################################################!
SUBROUTINE OUTPUT_RUN_REPORT()


INTEGER                                         ::  FILE_ID
INTEGER                                         ::  i
REAL(KIND = idp)                                ::  r, theta, phi, deltar
REAL(KIND = idp), DIMENSION(0:ITER_REPORT_NUM_SAMPLES)  ::  x_e
REAL(KIND = idp), DIMENSION(1:ITER_REPORT_NUM_SAMPLES)  ::  x_c, dx_c

REAL(KIND = idp)                                ::  Analytic_Val, Solver_Val
REAL(KIND = idp)                                ::  Return_Psi, Return_AlphaPsi
REAL(KIND = idp)                                ::  Return_Beta1, Return_Beta2, Return_Beta3
REAL(KIND = idp)                                ::  PsiPot_Val, AlphaPsiPot_Val


INTEGER                                         ::  Max_Iters
INTEGER, DIMENSION(:), ALLOCATABLE              ::  Hist_Iters

REAL(KIND = idp), DIMENSION(1:25)               ::  Time_Table

120 FORMAT (A61)
121 FORMAT (A1)
122 FORMAT (A41,I2.2)
123 FORMAT (A38,ES22.15)

109 FORMAT (A,I2.2,A,I2.2)
110 FORMAT (11X,A1,18X,A13,10X,A18,10X,A11,14X,A11,14X,A11)
111 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15)
112 FORMAT (A43,I2.2,A2,I2.2,A4)

FILE_ID = RUN_REPORT_FILE_ID


CALL PQ_ITERATIONS_MAX( Max_Iters )
ALLOCATE( Hist_Iters(1:Max_Iters) )
CALL PQ_ITERATIONS_HIST( Hist_Iters, Max_Iters )

CALL PQ_TIMETABLE_RUN( Time_Table )

! Write Timetable to File
IF ( RUN_REPORT_FLAG == 1 ) THEN

    WRITE(FILE_ID,'(A)')"                                Average Timing Results"
    WRITE(FILE_ID,'(A)')"            ============================================================="
    WRITE(FILE_ID,'(A)')" "
    WRITE(FILE_ID,123)"                    Initialize Time : ",TIME_TABLE(1)
    WRITE(FILE_ID,123)" Input/Communicate Source Data Time : ",TIME_TABLE(2)
    WRITE(FILE_ID,123)"     Input Boundary Conditions Time : ",TIME_TABLE(3)
    WRITE(FILE_ID,123)"        CFA_3D_Apply_BCs_Part1 Time : ",TIME_TABLE(4)
    WRITE(FILE_ID,120)"-------------------------------------------------------------"
    WRITE(FILE_ID,123)" ||     Calc_3D_Current_Values Time : ",TIME_TABLE(5)
    WRITE(FILE_ID,123)" ||    CREATE_3D_SubJcbn_Terms Time : ",TIME_TABLE(6)
    WRITE(FILE_ID,123)" ||       CREATE_3D_RHS_VECTOR Time : ",TIME_TABLE(7)
    WRITE(FILE_ID,123)"\  /     CREATE_3D_JCBN_MATRIX Time : ",TIME_TABLE(8)
    WRITE(FILE_ID,120)"-\/ ---------------------------------------------------------"
    WRITE(FILE_ID,123)"CREATE_3D_NONLAPLACIAN_STF_MAT Time : ",TIME_TABLE(9)
    WRITE(FILE_ID,123)"REDUCE_3D_NONLAPLACIAN_STF_MAT Time : ",TIME_TABLE(10)
    WRITE(FILE_ID,123)"FINISH_3D_NONLAPLACIAN_STF_MAT Time : ",TIME_TABLE(11)
    WRITE(FILE_ID,123)"          FINISH_3D_RHS_VECTOR Time : ",TIME_TABLE(12)
    WRITE(FILE_ID,123)"        CFA_3D_Apply_BCs_Part2 Time : ",TIME_TABLE(13)
    WRITE(FILE_ID,123)"                    CFA_Solver Time : ",TIME_TABLE(14)
    WRITE(FILE_ID,123)"        CFA_Coefficient_Update Time : ",TIME_TABLE(15)
    WRITE(FILE_ID,123)"   CFA_Coefficient_Share_PETSc Time : ",TIME_TABLE(16)
    WRITE(FILE_ID,123)"         CFA_Convergence_Check Time : ",TIME_TABLE(17)
    WRITE(FILE_ID,123)"               Total Iteration Time : ",TIME_TABLE(18)
    WRITE(FILE_ID,123)"             Poseidon_Dist_Sol Time : ",TIME_TABLE(19)
    WRITE(FILE_ID,120)"============================================================="
    WRITE(FILE_ID,121)" "
    WRITE(FILE_ID,121)" "
    WRITE(FILE_ID,121)" "
    WRITE(FILE_ID,121)" "



    

    WRITE(FILE_ID,'(A)')"                                 Number of Iterations"
    WRITE(FILE_ID,'(A)')"            ============================================================="
    WRITE(FILE_ID,'(A)')""
    WRITE(FILE_ID,'(A)')"Iterations   Times Occurred"
    DO i = 1,Max_Iters
        WRITE(FILE_ID,'(A,I2.2,A,I2.2)')"    ",i,"   |         ",Hist_Iters(i)
    END DO
    WRITE(FILE_ID,'(A)')" "
    WRITE(FILE_ID,'(A)')" "


END IF


! Write Timetable to Screen
IF (( WRITE_REPORT_FLAG == 1) .OR. (WRITE_REPORT_FLAG == 3) ) THEN
    WRITE(*,'(A)')" "
    WRITE(*,'(A)')" "
    WRITE(*,'(A)')"            ============================================================="
    WRITE(*,'(A)')"                                Average Timing Results"
    WRITE(*,'(A)')"            ============================================================="
    WRITE(*,'(A)')" "
    WRITE(*,123)"                    Initialize Time : ",TIME_TABLE(1)
    WRITE(*,123)" Input/Communicate Source Data Time : ",TIME_TABLE(2)
    WRITE(*,123)"     Input Boundary Conditions Time : ",TIME_TABLE(3)
    WRITE(*,123)"        CFA_3D_Apply_BCs_Part1 Time : ",TIME_TABLE(4)
    WRITE(*,120)"-------------------------------------------------------------"
    WRITE(*,123)" ||     Calc_3D_Current_Values Time : ",TIME_TABLE(5)
    WRITE(*,123)" ||    CREATE_3D_SubJcbn_Terms Time : ",TIME_TABLE(6)
    WRITE(*,123)" ||       CREATE_3D_RHS_VECTOR Time : ",TIME_TABLE(7)
    WRITE(*,123)"\  /     CREATE_3D_JCBN_MATRIX Time : ",TIME_TABLE(8)
    WRITE(*,120)"-\/ ---------------------------------------------------------"
    WRITE(*,123)"CREATE_3D_NONLAPLACIAN_STF_MAT Time : ",TIME_TABLE(9)
    WRITE(*,123)"REDUCE_3D_NONLAPLACIAN_STF_MAT Time : ",TIME_TABLE(10)
    WRITE(*,123)"FINISH_3D_NONLAPLACIAN_STF_MAT Time : ",TIME_TABLE(11)
    WRITE(*,123)"          FINISH_3D_RHS_VECTOR Time : ",TIME_TABLE(12)
    WRITE(*,123)"        CFA_3D_Apply_BCs_Part2 Time : ",TIME_TABLE(13)
    WRITE(*,123)"                    CFA_Solver Time : ",TIME_TABLE(14)
    WRITE(*,123)"        CFA_Coefficient_Update Time : ",TIME_TABLE(15)
    WRITE(*,123)"   CFA_Coefficient_Share_PETSc Time : ",TIME_TABLE(16)
    WRITE(*,123)"         CFA_Convergence_Check Time : ",TIME_TABLE(17)
    WRITE(*,123)"               Total Iteration Time : ",TIME_TABLE(18)
    WRITE(*,123)"             Poseidon_Dist_Sol Time : ",TIME_TABLE(19)
    WRITE(*,120)"============================================================="
    WRITE(*,121)" "
    WRITE(*,121)" "
    WRITE(*,121)" "
    WRITE(*,121)" "

END IF




! Write Results Table Header to Screen
IF (( WRITE_REPORT_FLAG == 1) .OR. (WRITE_REPORT_FLAG == 3) ) THEN
    PRINT*,"++++++++++++++++++++++++++ Sample Run Results ++++++++++++++++++++++++++"
    WRITE(*,110)"r","Psi Potential","AlphaPsi Potential","Beta Value1","Beta Value2","Beta Value3"
END IF


Call Create_Logarithmic_1D_Mesh( DRIVER_INNER_RADIUS, DRIVER_OUTER_RADIUS, ITER_REPORT_NUM_SAMPLES,     &
                                 x_e, x_c, dx_c )


!deltar = ( DRIVER_OUTER_RADIUS - DRIVER_INNER_RADIUS )/ REAL(ITER_REPORT_NUM_SAMPLES, KIND = idp)
DO i = 0,ITER_REPORT_NUM_SAMPLES

!    r = i*deltar + DRIVER_INNER_RADIUS
    r = x_e(i)*Centimeter
    theta = pi/8.0_idp
    phi = pi/2.0_idp
    

    CALL Calc_3D_Values_At_Location( r, theta, phi,                              &
                                    Return_Psi, Return_AlphaPsi,                &
                                    Return_Beta1, Return_Beta2, Return_Beta3    )


    ! Calculate Conformal Factor value from Newtonian Potential
    PsiPot_Val = 2.0_idp*C_Square*(1.0_idp - Return_Psi)/GravPot_Units

    ! Calculate the product of the Conformal Factor and Lapse Function from Newtonian Potential
    AlphaPsiPot_Val = 2.0_idp*C_Square*(Return_AlphaPsi - 1.0_idp)/GravPot_Units


    ! Write Results to Screen
    IF (( WRITE_REPORT_FLAG == 1) .OR. (WRITE_REPORT_FLAG == 3) ) THEN
        WRITE(*,111) r/Centimeter,              &
                     PsiPot_Val,                &
                     AlphaPsiPot_Val,           &
                     Return_Beta1/Shift_Units,  &
                     Return_Beta2,              &
                     Return_Beta3
    END IF

    ! Write Results to File
    IF ( RUN_REPORT_FLAG == 1 ) THEN
        WRITE(FILE_ID,111) r/Centimeter,              &
                           PsiPot_Val,                &
                           AlphaPsiPot_Val,           &
                           Return_Beta1/Shift_Units,  &
                           Return_Beta2,              &
                           Return_Beta3
    END IF

END DO


DEALLOCATE( Hist_Iters )

END SUBROUTINE OUTPUT_RUN_REPORT





 !+203+############################################################################!
!                                                                                   !
!                    CLOSE_RUN_REPORT_FILE                                          !
!                                                                                   !
 !#################################################################################!
SUBROUTINE CLOSE_RUN_REPORT_FILE()

IF ( RUN_REPORT_FLAG == 1 ) THEN
    CLOSE( UNIT = RUN_REPORT_FILE_ID )
END IF

END SUBROUTINE CLOSE_RUN_REPORT_FILE

















 !+301+############################################################################!
!                                                                                   !
!                     OPEN_FRAME_REPORT_FILE                                          !
!                                                                                   !
 !#################################################################################!
SUBROUTINE OPEN_FRAME_REPORT_FILE(Frame)

INTEGER, INTENT(IN)                                     ::  Frame

CHARACTER(LEN = 70)                                     ::  FILE_NAME
INTEGER                                                 ::  istat
LOGICAL                                                 ::  FLAG, OK

109 FORMAT (A,I2.2)



IF ( myID == 0 ) THEN
    IF ( FRAME_REPORT_FLAG == 1 ) THEN

        WRITE(FILE_NAME,'(A,A)') Poseidon_IterReports_Dir,"Frame_Report_SCRATCH.out"
        CALL OPEN_NEW_FILE( FILE_NAME, FRAME_REPORT_FILE_ID )

    END IF ! WRITE_REPORT_FLAGS
END IF ! myID == 0


END SUBROUTINE OPEN_FRAME_REPORT_FILE




!+302+##########################################################################!
!                                                                               !
!                   OUTPUT_ITERATION_REPORT                                     !
!                                                                               !
!###############################################################################!
SUBROUTINE OUTPUT_FRAME_REPORT(FRAME)

INTEGER, INTENT(IN)                 :: FRAME

INTEGER                                                 ::  FILE_ID, FILE_IDb
INTEGER                                                 ::  i
REAL(KIND = idp)                                        ::  r, theta, phi, deltar
REAL(KIND = idp), DIMENSION(0:ITER_REPORT_NUM_SAMPLES)  ::  x_e
REAL(KIND = idp), DIMENSION(1:ITER_REPORT_NUM_SAMPLES)  ::  x_c, dx_c


REAL(KIND = idp)                                        ::  Analytic_Val, Solver_Val
REAL(KIND = idp)                                        ::  Return_Psi, Return_AlphaPsi
REAL(KIND = idp)                                        ::  Return_Beta1, Return_Beta2, Return_Beta3
REAL(KIND = idp)                                        ::  PsiPot_Val, AlphaPsiPot_Val

CHARACTER(LEN = 70)                                     ::  FILE_NAME

CHARACTER(LEN = 300)                                    ::  Line
INTEGER                                                 ::  io_stat
INTEGER                                                 ::  ITER_REPORT_NUM_SAMPLES = 20

REAL(KIND = idp), DIMENSION(1:25)                       ::  Time_Table


120 FORMAT (12X,A)
121 FORMAT (A1)
122 FORMAT (A,I2.2)
123 FORMAT (12X,A38,ES22.15)

109 FORMAT (A,I2.2,A,I2.2)
110 FORMAT (11X,A1,16X,A18,9X,A13,10X,A18,10X,A11,14X,A11,14X,A11)
111 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15)
112 FORMAT (A43,I2.2,A2,I2.2,A4)

132 FORMAT (A,A,I2.2,A)

142 FORMAT (16X,A,10X,A)
143 FORMAT (19X,I2.2,15X,ES22.15)



CALL PQ_TIMETABLE_FRAME( Time_Table )

WRITE(FILE_NAME,132)Poseidon_IterReports_Dir,"Frame_Report_",Frame,".out"
CALL OPEN_NEW_FILE( FILE_NAME, FILE_ID )


IF ( myID == 0 ) THEN
    ! Write Title to File
    IF ( FRAME_REPORT_FLAG == 1 ) THEN

        WRITE(FILE_ID,109)"                                Report for Frame ",Frame
        WRITE(FILE_ID,'(A)')"-----------------------------------------------------------------------------------------"
        WRITE(FILE_ID,'(2/)')


        IF ( Cur_Iteration == 2 ) THEN
            WRITE(FILE_ID,'(A,I2.2,A)')"The Newton-Raphson solver exited after ",Cur_Iteration-1," Iteration."
        ELSE
            WRITE(FILE_ID,'(A,I2.2,A)')"The Newton-Raphson solver exited after ",Cur_Iteration-1," Iterations."
        END IF


        IF ( CONVERGENCE_FLAG == 1 ) THEN

            WRITE(FILE_ID,'(A)')"The solution converged within the tolerance set."

        ELSE IF ( CONVERGENCE_FLAG == 2 ) THEN

            WRITE(FILE_ID,'(A)')"The solution did not converge within the maximum number of iterations allowed."
            WRITE(FILE_ID,'(A,I2.2)')"The maximum number of iterations allowed is ",Max_Iterations

        END IF
        WRITE(FILE_ID,'(2/)')




        WRITE(FILE_ID,'(A)')"                              Average Timing Results"
        WRITE(FILE_ID,'(A)')"            ============================================================="
        WRITE(FILE_ID,'(A)')" "
        WRITE(FILE_ID,123)"                    Initialize Time : ",TIME_TABLE(1)
        WRITE(FILE_ID,123)" Input/Communicate Source Data Time : ",TIME_TABLE(2)
        WRITE(FILE_ID,123)"     Input Boundary Conditions Time : ",TIME_TABLE(3)
        WRITE(FILE_ID,123)"        CFA_3D_Apply_BCs_Part1 Time : ",TIME_TABLE(4)
        WRITE(FILE_ID,120)"-------------------------------------------------------------"
        WRITE(FILE_ID,123)" ||     Calc_3D_Current_Values Time : ",TIME_TABLE(5)
        WRITE(FILE_ID,123)" ||    CREATE_3D_SubJcbn_Terms Time : ",TIME_TABLE(6)
        WRITE(FILE_ID,123)" ||       CREATE_3D_RHS_VECTOR Time : ",TIME_TABLE(7)
        WRITE(FILE_ID,123)"\  /     CREATE_3D_JCBN_MATRIX Time : ",TIME_TABLE(8)
        WRITE(FILE_ID,120)"-\/ ---------------------------------------------------------"
        WRITE(FILE_ID,123)"CREATE_3D_NONLAPLACIAN_STF_MAT Time : ",TIME_TABLE(9)
        WRITE(FILE_ID,123)"REDUCE_3D_NONLAPLACIAN_STF_MAT Time : ",TIME_TABLE(10)
        WRITE(FILE_ID,123)"FINISH_3D_NONLAPLACIAN_STF_MAT Time : ",TIME_TABLE(11)
        WRITE(FILE_ID,123)"          FINISH_3D_RHS_VECTOR Time : ",TIME_TABLE(12)
        WRITE(FILE_ID,123)"        CFA_3D_Apply_BCs_Part2 Time : ",TIME_TABLE(13)
        WRITE(FILE_ID,123)"                    CFA_Solver Time : ",TIME_TABLE(14)
        WRITE(FILE_ID,123)"        CFA_Coefficient_Update Time : ",TIME_TABLE(15)
        WRITE(FILE_ID,123)"   CFA_Coefficient_Share_PETSc Time : ",TIME_TABLE(16)
        WRITE(FILE_ID,123)"         CFA_Convergence_Check Time : ",TIME_TABLE(17)
        WRITE(FILE_ID,123)"               Total Iteration Time : ",TIME_TABLE(18)
        WRITE(FILE_ID,123)"             Poseidon_Dist_Sol Time : ",TIME_TABLE(19)
        WRITE(FILE_ID,120)"============================================================="
        WRITE(FILE_ID,121)" "
        WRITE(FILE_ID,121)" "
        WRITE(FILE_ID,121)" "
        WRITE(FILE_ID,121)" "



        WRITE(FILE_ID,'(A)')"                                  Frame's Convergence"
        WRITE(FILE_ID,'(A)')"            ============================================================="
        WRITE(FILE_ID,'(A)')" "
        WRITE(FILE_ID,142)"Iteration","MAXVAL(ABS(Update_Vector))"
        DO i = 1,CUR_ITERATION-1
            WRITE(FILE_ID,143)i,FRAME_CONVERGENCE_TABLE(i)
        END DO
        WRITE(FILE_ID,'(2/)')




        WRITE(FILE_ID,'(A)')"                                   Frame's Results"
        WRITE(FILE_ID,'(A)')"            ============================================================="
        WRITE(FILE_ID,'(A)')" "
        WRITE(FILE_ID,110)"r","Analytic Potential","Psi Potential","AlphaPsi Potential","Beta Value1","Beta Value2","Beta Value3"



        Call Create_Logarithmic_1D_Mesh( DRIVER_INNER_RADIUS, DRIVER_OUTER_RADIUS, ITER_REPORT_NUM_SAMPLES,     &
                                         x_e, x_c, dx_c )
!        deltar = ( DRIVER_OUTER_RADIUS - DRIVER_INNER_RADIUS )/ REAL(ITER_REPORT_NUM_SAMPLES, KIND = idp)

        DO i = 0,ITER_REPORT_NUM_SAMPLES

 !          r = i*deltar + DRIVER_OUTER_RADIUS
            r = x_e(i)*Centimeter
            theta = pi/2.0_idp
            phi = pi/2.0_idp


            CALL Calc_3D_Values_At_Location( r, theta, phi,                              &
                                            Return_Psi, Return_AlphaPsi,                &
                                            Return_Beta1, Return_Beta2, Return_Beta3    )

            ! Determine the Newtonian Potential at the location r, theta, phi
            Analytic_Val = Potential_Solution(r,theta,phi)/GravPot_Units


            ! AlphaPsi_to_Pot   =   2*C_Square*(AlphaPsi - 1)
            ! Psi_to_Pot        =   2*C_Square*(1 - Psi)

            ! Calculate Conformal Factor value from Newtonian Potential
            PsiPot_Val = 2.0_idp*C_Square*(Return_Psi - 1.0_idp)/GravPot_Units

            ! Calculate the product of the Conformal Factor and Lapse Function from Newtonian Potential
            AlphaPsiPot_Val = 2.0_idp*C_Square*(1.0_idp - Return_AlphaPsi)/GravPot_Units

            ! Write Results to File
            WRITE(FILE_ID,111)  r/Centimeter,               &
                                Analytic_Val,               &
                                PsiPot_Val,                 &
                                AlphaPsiPot_Val,            &
                                Return_Beta1/Shift_Units,   &
                                Return_Beta2,               &
                                Return_Beta3

        END DO
        WRITE( FILE_ID, '(4/)')




        ! Copy from Scratch File into Frame Report !
        REWIND( UNIT=FRAME_REPORT_FILE_ID )
        DO
            READ(FRAME_REPORT_FILE_ID,'(A)', iostat=io_stat) Line
            IF (io_stat /= 0 ) THEN
                EXIT
            END IF
            WRITE(FILE_ID, '(A)') Line
        END DO


    END IF ! FRAME_REPORT_FLAG == 1

    CLOSE( UNIT = FILE_ID)
    CLOSE( UNIT = FRAME_REPORT_FILE_ID, STATUS='delete' )

END IF ! myID == 0

END SUBROUTINE OUTPUT_FRAME_REPORT





 !+303+############################################################################!
!                                                                                   !
!                    CLOSE_RUN_REPORT_FILE                                          !
!                                                                                   !
 !#################################################################################!
SUBROUTINE CLOSE_FRAME_REPORT_FILE()


CLOSE( UNIT = FRAME_REPORT_FILE_ID, STATUS='delete' )


END SUBROUTINE CLOSE_FRAME_REPORT_FILE






 !+501+############################################################################!
!                                                                                   !
!                     OPEN_NEW_FILE                                                 !
!                                                                                   !
 !#################################################################################!
SUBROUTINE OPEN_NEW_FILE(File_Name, File_Number)



CHARACTER(LEN = *), INTENT(IN)                          ::  File_Name
INTEGER,            INTENT(INOUT)                       ::  File_Number

INTEGER                                                 ::  Temp_Number
INTEGER                                                 ::  istat
LOGICAL                                                 ::  FLAG, OP, EX
LOGICAL                                                 ::  UNIT_FLAG, NAME_FLAG


UNIT_FLAG = .FALSE.
NAME_FLAG = .FALSE.


!  Assigned an unused number, and assign it to new file
FLAG = .TRUE.
Temp_Number = 421
DO WHILE (FLAG)
    INQUIRE( UNIT = Temp_Number, OPENED = OP )

    IF ( OP ) THEN
        Temp_Number = Temp_Number + 1
    ELSE
        File_Number = Temp_Number
        FLAG = .FALSE.
        UNIT_FLAG = .TRUE.
    END IF
END DO




! Check if file already exists !
!INQUIRE( FILE = File_Name, EXIST = EX )
!IF ( EX ) THEN
!    PRINT*,"File ",File_Name," is already opened"
!
!ELSE
!    PRINT*,File_Name," is not already opened."
!    NAME_FLAG = .TRUE.
!END IF



! Open New File
IF ( UNIT_FLAG  ) THEN

    OPEN( Unit = File_Number, File = File_Name, IOSTAT = istat )
    IF ( istat .NE. 0 ) THEN

        PRINT*,"WARNING: Could not open Run Report File at ", File_Name

    END IF
END IF


END SUBROUTINE OPEN_NEW_FILE













 !+701+############################################################################!
!                                                                                   !
!                     OUTPUT_ITERATION_HISTORY                                      !
!                                                                                   !
 !#################################################################################!
SUBROUTINE OUTPUT_ITERATION_HISTORY()

INTEGER                                         ::  i
INTEGER                                         ::  FILE_ID
CHARACTER(LEN = 65)                             ::  FILE_NAME





WRITE(FILE_NAME,'(A,A)')Poseidon_IterReports_Dir,"Iteration_Histogram.out"
CALL OPEN_NEW_FILE( FILE_NAME, FILE_ID )

DO i = 1,DRIVER_TOTAL_FRAMES

    WRITE(FILE_ID,'(I3.3,A,I3.3)')i,"   ",Iteration_History(i)

END DO

CLOSE( FILE_ID )

END SUBROUTINE OUTPUT_ITERATION_HISTORY







SUBROUTINE OUTPUT_PRIMATIVES( Density, Velocity, Num_Entries )

REAL(KIND = idp), DIMENSION(1:Num_Entries), INTENT(IN)  :: Density
REAL(KIND = idp), DIMENSION(1:Num_Entries), INTENT(IN)  :: Velocity
INTEGER,                                    INTENT(IN)  :: Num_Entries

CHARACTER(LEN = 100), DIMENSION(:), ALLOCATABLE             ::  Filenames
INTEGER, DIMENSION(:), ALLOCATABLE                          ::  File_IDs
INTEGER                                                     ::  Num_Files = 2
INTEGER                                                     ::  i

116 FORMAT (A,A,I5.5,A)

ALLOCATE( Filenames(1:Num_Files) )
ALLOCATE( File_IDs(1:Num_Files) )

WRITE(Filenames(1),116) Poseidon_Sources_Dir,"Sources_Density_",Poseidon_Frame,".out"
WRITE(Filenames(2),116) Poseidon_Sources_Dir,"Sources_Velocity_",Poseidon_Frame,".out"

! Open Files
File_IDs = [(161 + i, i =1,Num_Files)]
DO i = 1,Num_Files
    CALL OPEN_NEW_FILE( Filenames(i), File_IDs(i) )
END DO

! Write to Files
DO i = 1,Num_Entries
    WRITE(File_IDs(1),*)Density(i) /(Gram/Centimeter**3)
    WRITE(FILE_IDs(2),*)Velocity(i)/(Centimeter/Second)
END DO

! Close Files
DO i = 1,Num_Files
    CLOSE( Unit = File_IDs(i))
END DO



END SUBROUTINE OUTPUT_PRIMATIVES





END MODULE Driver_IO_Functions_Module
