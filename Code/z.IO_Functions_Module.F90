   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE IO_Functions_Module                                                          !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   OUTPUT_ITER_TIMETABLE                                                    !##!
!##!    +102+   CLOCK_IN                                                            !##!
!##!                                                                                !##!
!##!    +201+   OPEN_RUN_REPORT_FILE                                                !##!
!##!    +202+   CLOSE_RUN_REPORT_FILE                                               !##!
!##!                                                                                !##!
!##!    +301+   OPEN_ITER_REPORT_FILE                                               !##!
!##!    +302+   CLOSE_ITER_REPORT_FILE                                              !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!

!*D*================================!
!                                   !
!           Dependencies            !
!                                   !
!===================================!



USE Poseidon_Constants_Module, &
                    ONLY : idp, pi

USE Poseidon_Parameters, &
                    ONLY :  DOMAIN_DIM,                                     &
                            DEGREE,                                         &
                            L_LIMIT,                                        &
                            NUM_R_QUAD_POINTS,                              &
                            NUM_T_QUAD_POINTS,                              &
                            NUM_P_QUAD_POINTS,                              &
                            NUM_R_ELEMS_PER_BLOCK,                          &
                            CUR_ITERATION,                                  &
                            WRITE_TIMETABLE_FLAG,                           &
                            WRITE_RESULTS_FLAG,                             &
                            RUN_REPORT_FILE_ID,                             &
                            ITER_REPORT_FILE_ID,                            &
                            WRITE_REPORT_FLAG

USE CHIMERA_Parameters,  &
                    ONLY :  Analytic_Solution,                              &
                            Shift_Solution

USE Global_Variables_And_Parameters, &
                    ONLY :  NUM_R_NODES,                                    &
                            NUM_R_ELEMENTS,                                 &
                            RHS_Vector,                                     &
                            Coefficient_Vector,                             &
                            rlocs, tlocs, plocs,                            &
                            R_INNER, R_OUTER,                               &
                            INT_R_LOCATIONS,                                &
                            INT_T_LOCATIONS,                                &
                            INT_P_LOCATIONS,                                &
                            VAR_DIM,                                        &
                            NUM_OFF_DIAGONALS,                              &
                            Num_Timer_Calls,                                &
                            ITER_TIME_TABLE,                                &
                            RUN_TIME_TABLE
 

USE Additional_Functions_Module, &
                    ONLY :  Map_From_x_Space,                               &
                            Initialize_LGL_Quadrature_Locations,            &
                            Lagrange_Poly,                                  &
                            Lagrange_Poly_Deriv,                            &
                            Lagrange_Second_Deriv,                          &
                            CFA_Matrix_Map





USE Mesh_Module, &
                    ONLY :  Create_Uniform_1D_Mesh,                         &
                            Create_Logarithmic_1D_Mesh,                     &
                            Create_Split_1D_Mesh


IMPLICIT NONE


!*F&S*==========================================!
!                                               !
!           Functions & Subroutines             !
!                                               !
!===============================================!
CONTAINS







 !+101+############################################################################!
!                                                                                   !
!                           OUTPUT_ITER_TIMETABLE                                        !
!                                                                                   !
 !#################################################################################!
SUBROUTINE OUTPUT_ITER_TIMETABLE(Ident)


INTEGER,INTENT(IN)              ::  Ident

CHARACTER(LEN = 23)             ::  FILENAME
INTEGER                         ::  FILE_ID

110 FORMAT (A61)
111 FORMAT (A1)
112 FORMAT (A41,I2.2)
113 FORMAT (A38,ES22.15)

IF (( WRITE_TIMETABLE_FLAG == 1 ) .OR. ( WRITE_TIMETABLE_FLAG == 3 ) ) THEN

    WRITE(*,110)"============================================================="
    WRITE(*,111)" "
    WRITE(*,112)"                 Time Table for Process ",Ident
    WRITE(*,111)" "
    WRITE(*,110)"============================================================="
    WRITE(*,113)"                    Initialize Time : ",ITER_TIME_TABLE(1)
    WRITE(*,113)" Input/Communicate Source Data Time : ",ITER_TIME_TABLE(2)
    WRITE(*,113)"     Input Boundary Conditions Time : ",ITER_TIME_TABLE(3)
    WRITE(*,113)"        CFA_3D_Apply_BCs_Part1 Time : ",ITER_TIME_TABLE(4)
    WRITE(*,110)"-------------------------------------------------------------"
    WRITE(*,113)" ||     Calc_3D_Current_Values Time : ",ITER_TIME_TABLE(5)
    WRITE(*,113)" ||    CREATE_3D_SubJcbn_Terms Time : ",ITER_TIME_TABLE(6)
    WRITE(*,113)" ||       CREATE_3D_RHS_VECTOR Time : ",ITER_TIME_TABLE(7)
    WRITE(*,113)"\  /     CREATE_3D_JCBN_MATRIX Time : ",ITER_TIME_TABLE(8)
    WRITE(*,110)"-\/ ---------------------------------------------------------"
    WRITE(*,113)"CREATE_3D_NONLAPLACIAN_STF_MAT Time : ",ITER_TIME_TABLE(9)
    WRITE(*,113)"REDUCE_3D_NONLAPLACIAN_STF_MAT Time : ",ITER_TIME_TABLE(10)
    WRITE(*,113)"FINISH_3D_NONLAPLACIAN_STF_MAT Time : ",ITER_TIME_TABLE(11)
    WRITE(*,113)"          FINISH_3D_RHS_VECTOR Time : ",ITER_TIME_TABLE(12)
    WRITE(*,113)"        CFA_3D_Apply_BCs_Part2 Time : ",ITER_TIME_TABLE(13)
    WRITE(*,113)"                    CFA_Solver Time : ",ITER_TIME_TABLE(14)
    WRITE(*,113)"        CFA_Coefficient_Update Time : ",ITER_TIME_TABLE(15)
    WRITE(*,113)"   CFA_Coefficient_Share_PETSc Time : ",ITER_TIME_TABLE(16)
    WRITE(*,113)"         CFA_Convergence_Check Time : ",ITER_TIME_TABLE(17)
    WRITE(*,113)"               Total Iteration Time : ",ITER_TIME_TABLE(18)
    WRITE(*,113)"             Poseidon_Dist_Sol Time : ",ITER_TIME_TABLE(19)
    WRITE(*,110)"============================================================="


END IF


IF (( WRITE_TIMETABLE_FLAG == 2 ) .OR. ( WRITE_TIMETABLE_FLAG == 3 ) ) THEN

    FILE_ID = 42
    WRITE(FILENAME,'(A,I2.2,A)')'OUTPUT/Timetable_',Ident,'.out'

    OPEN(unit = FILE_ID,file = FILENAME)

    WRITE(42,110)"============================================================="
    WRITE(42,111)" "
    WRITE(42,112)"                 Time Table for Process ",Ident
    WRITE(42,111)" "
    WRITE(42,110)"============================================================="
    WRITE(42,113)"                    Initialize Time : ",ITER_TIME_TABLE(1)
    WRITE(42,113)" Input/Communicate Source Data Time : ",ITER_TIME_TABLE(2)
    WRITE(42,113)"     Input Boundary Conditions Time : ",ITER_TIME_TABLE(3)
    WRITE(42,113)"        CFA_3D_Apply_BCs_Part1 Time : ",ITER_TIME_TABLE(4)
    WRITE(42,110)"-------------------------------------------------------------"
    WRITE(42,113)" ||     Calc_3D_Current_Values Time : ",ITER_TIME_TABLE(5)
    WRITE(42,113)" ||    CREATE_3D_SubJcbn_Terms Time : ",ITER_TIME_TABLE(6)
    WRITE(42,113)" ||       CREATE_3D_RHS_VECTOR Time : ",ITER_TIME_TABLE(7)
    WRITE(42,113)"\  /     CREATE_3D_JCBN_MATRIX Time : ",ITER_TIME_TABLE(8)
    WRITE(42,110)"-\/ ---------------------------------------------------------"
    WRITE(42,113)"CREATE_3D_NONLAPLACIAN_STF_MAT Time : ",ITER_TIME_TABLE(9)
    WRITE(42,113)"REDUCE_3D_NONLAPLACIAN_STF_MAT Time : ",ITER_TIME_TABLE(10)
    WRITE(42,113)"FINISH_3D_NONLAPLACIAN_STF_MAT Time : ",ITER_TIME_TABLE(11)
    WRITE(42,113)"          FINISH_3D_RHS_VECTOR Time : ",ITER_TIME_TABLE(12)
    WRITE(42,113)"        CFA_3D_Apply_BCs_Part2 Time : ",ITER_TIME_TABLE(13)
    WRITE(42,113)"                    CFA_Solver Time : ",ITER_TIME_TABLE(14)
    WRITE(42,113)"        CFA_Coefficient_Update Time : ",ITER_TIME_TABLE(15)
    WRITE(42,113)"   CFA_Coefficient_Share_PETSc Time : ",ITER_TIME_TABLE(16)
    WRITE(42,113)"         CFA_Convergence_Check Time : ",ITER_TIME_TABLE(17)
    WRITE(42,113)"               Total Iteration Time : ",ITER_TIME_TABLE(18)
    WRITE(42,113)"             Poseidon_Dist_Sol Time : ",ITER_TIME_TABLE(19)
    WRITE(42,110)"============================================================="




END IF



END SUBROUTINE OUTPUT_ITER_TIMETABLE







 !+102+############################################################################!
!                                                                                   !
!                           CLOCK_IN                                                !
!                                                                                   !
 !#################################################################################!
SUBROUTINE CLOCK_IN( Time, Ident )


REAL(KIND = idp), INTENT(IN)                 ::   Time
INTEGER, INTENT(IN)                          ::   Ident


! Add Time to Iteration Time Table
ITER_TIME_TABLE(Ident) = Time


! Add Times to the Run Time Table
IF (( Ident <= 3 ) .OR. (Ident == 19) ) THEN
    ! The 3 first events and last event only happen once !
    RUN_TIME_TABLE(Ident) = Time
ELSE
    ! Add This Iterations Time to Running Average for the Run
    RUN_TIME_TABLE(Ident) = (RUN_TIME_TABLE(ident)*(Cur_Iteration-1) + Time)/Cur_Iteration
END IF


END SUBROUTINE CLOCK_IN







 !+201+############################################################################!
!                                                                                   !
!                     OPEN_RUN_REPORT_FILE                                          !
!                                                                                   !
 !#################################################################################!
SUBROUTINE OPEN_RUN_REPORT_FILE()

CHARACTER(LEN = 22)                                     ::  FILE_NAME
INTEGER                                                 ::  istat
LOGICAL                                                 ::  FLAG, OK


FLAG = .TRUE.
RUN_REPORT_FILE_ID = 421
DO WHILE (FLAG)
    INQUIRE( UNIT = RUN_REPORT_FILE_ID, OPENED = OK)
    IF ( OK ) THEN
        RUN_REPORT_FILE_ID = RUN_REPORT_FILE_ID + 1
    ELSE
        FLAG = .FALSE.
    END IF
END DO


WRITE(FILE_NAME,'(A)')"OUTPUT/Run_Report.out"
OPEN( Unit = RUN_REPORT_FILE_ID, File = FILE_NAME, IOSTAT = istat )
IF ( istat .NE. 0 ) THEN

    PRINT*,"WARNING: Could not open Run Report File at 'OUTPUT/Run_Report.out'."

END IF


END SUBROUTINE OPEN_RUN_REPORT_FILE





 !+202+############################################################################!
!                                                                                   !
!                    CLOSE_RUN_REPORT_FILE                                          !
!                                                                                   !
 !#################################################################################!
SUBROUTINE CLOSE_RUN_REPORT_FILE()


CLOSE( UNIT = RUN_REPORT_FILE_ID )


END SUBROUTINE CLOSE_RUN_REPORT_FILE














 !+301+############################################################################!
!                                                                                   !
!                     OPEN_RUN_REPORT_FILE                                          !
!                                                                                   !
 !#################################################################################!
SUBROUTINE OPEN_ITER_REPORT_FILE(Iteration, Rank)

INTEGER, INTENT(IN)                                     ::  Iteration, Rank

CHARACTER(LEN = 53)                                     ::  FILE_NAME
INTEGER                                                 ::  istat
LOGICAL                                                 ::  FLAG, OK

109 FORMAT (A,I2.2,A,I2.2)
112 FORMAT (A43,I2.2,A2,I2.2,A4)


IF (( WRITE_REPORT_FLAG == 2) .OR. (WRITE_REPORT_FLAG == 3) ) THEN

    FLAG = .TRUE.
    ITER_REPORT_FILE_ID = 842
    DO WHILE (FLAG)
        INQUIRE( UNIT = ITER_REPORT_FILE_ID, OPENED = OK)
        IF ( OK ) THEN
            ITER_REPORT_FILE_ID = ITER_REPORT_FILE_ID + 1
        ELSE
            FLAG = .FALSE.
        END IF
    END DO



    WRITE(FILE_NAME,112)"OUTPUT/Iteration_Reports/Iteration_Report_P",Rank,"_I",Iteration,".out"
    OPEN( Unit = ITER_REPORT_FILE_ID, File = FILE_NAME, IOSTAT = istat )
    IF ( istat .NE. 0 ) THEN

        PRINT*,"WARNING: Could not open Run Report File at",FILE_NAME

    END IF




    WRITE(ITER_REPORT_FILE_ID,109)"                      Report for Iteration ",Iteration," by Process ",Rank
    WRITE(ITER_REPORT_FILE_ID,'(A)')"-----------------------------------------------------------------------------------------"
    WRITE(ITER_REPORT_FILE_ID,'(A)')" "

END IF


END SUBROUTINE OPEN_ITER_REPORT_FILE




 !+302+############################################################################!
!                                                                                   !
!                    CLOSE_ITER_REPORT_FILE                                          !
!                                                                                   !
 !#################################################################################!
SUBROUTINE CLOSE_ITER_REPORT_FILE()


CLOSE( UNIT = ITER_REPORT_FILE_ID )


END SUBROUTINE CLOSE_ITER_REPORT_FILE








END MODULE IO_Functions_Module
