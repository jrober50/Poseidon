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
!##!    +101+   OUTPUT_ITER_TIMETABLE                                               !##!
!##!    +102+   CLOCK_IN                                                            !##!
!##!                                                                                !##!
!##!    +201+   OPEN_RUN_REPORT_FILE                                                !##!
!##!    +202+   CLOSE_RUN_REPORT_FILE                                               !##!
!##!                                                                                !##!
!##!    +301+   OPEN_ITER_REPORT_FILE                                               !##!
!##!    +302+   CLOSE_ITER_REPORT_FILE                                              !##!
!##!                                                                                !##!
!##!    +401+   OUTPUT_FINAL_RESULTS                                                !##!
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
                    ONLY : idp, pi, C_Square

USE Poseidon_Parameters, &
                    ONLY :  DOMAIN_DIM,                                     &
                            DEGREE,                                         &
                            L_LIMIT,                                        &
                            NUM_R_QUAD_POINTS,                              &
                            NUM_T_QUAD_POINTS,                              &
                            NUM_P_QUAD_POINTS,                              &
                            NUM_R_ELEMS_PER_BLOCK,                          &
                            CUR_ITERATION,                                  &
                            MAX_ITERATIONS,                                 &
                            CONVERGENCE_FLAG,                               &
                            WRITE_TIMETABLE_FLAG,                           &
                            WRITE_RESULTS_FLAG,                             &
                            ITER_REPORT_NUM_SAMPLES,                        &
                            RUN_REPORT_FILE_ID,                             &
                            ITER_REPORT_FILE_ID,                            &
                            FRAME_REPORT_FILE_ID,                           &
                            WRITE_REPORT_FLAG

USE DRIVER_Parameters,  &
                    ONLY :  Analytic_Solution,                              &
                            Shift_Solution,                                 &
                            myID,                                           &
                            SOURCE_OUTPUT_FLAG,                             &
                            RESULTS_OUTPUT_FLAG,                            &
                            RUN_REPORT_FLAG,                                &
                            FRAME_REPORT_FLAG,                              &
                            DRIVER_TEST_NUMBER,                             &
                            DRIVER_FRAME

USE Poseidon_Variables_Module, &
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
                            Total_Run_Iters,                                &
                            Num_Timer_Calls,                                &
                            ITER_TIME_TABLE,                                &
                            FRAME_TIME_TABLE,                               &
                            FRAME_CONVERGENCE_TABLE,                        &
                            RUN_TIME_TABLE
 

USE Additional_Functions_Module, &
                    ONLY :  Map_From_x_Space,                               &
                            Initialize_LGL_Quadrature_Locations,            &
                            Lagrange_Poly,                                  &
                            Lagrange_Poly_Deriv,                            &
                            Lagrange_Second_Deriv,                          &
                            CFA_Matrix_Map

USE Additional_Functions_Module,     &
                    ONLY : Calc_3D_Values_At_Location


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


! Add Times to the Frame Time Table
IF (( Ident <= 3 ) .OR. (Ident == 19) ) THEN
    ! The 3 first events and last event only happen once !
    FRAME_TIME_TABLE(Ident) = Time
ELSE
    ! Add This Iterations Time to Running Average for the Run
    FRAME_TIME_TABLE(Ident) = (FRAME_TIME_TABLE(ident)*(Cur_Iteration-1) + Time)/Cur_Iteration
END IF



! Add Times to the Frame Time Table
IF (( Ident <= 3 ) .OR. (Ident == 19) ) THEN
    ! The 3 first events and last event only happen once per Frame !
    RUN_TIME_TABLE(Ident) = (RUN_TIME_TABLE(ident)*(DRIVER_FRAME-1) + Time)/DRIVER_FRAME
ELSE
    ! Add This Iterations Time to Running Average for the Run
    RUN_TIME_TABLE(Ident) = (RUN_TIME_TABLE(ident)*(Total_Run_Iters-1) + Time)/Total_Run_Iters
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


IF ( RUN_REPORT_FLAG == 1 ) THEN
    WRITE(FILE_NAME,'(A)')"OUTPUT/Run_Report.out"
    CALL OPEN_NEW_FILE( FILE_NAME, RUN_REPORT_FILE_ID )
END IF


END SUBROUTINE OPEN_RUN_REPORT_FILE





 !+202+############################################################################!
!                                                                                   !
!                    CLOSE_RUN_REPORT_FILE                                          !
!                                                                                   !
 !#################################################################################!
SUBROUTINE CLOSE_RUN_REPORT_FILE()

IF ( RUN_REPORT_FLAG == 1 ) THEN
    CLOSE( UNIT = RUN_REPORT_FILE_ID )
END IF

END SUBROUTINE CLOSE_RUN_REPORT_FILE






 !+201+############################################################################!
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



    WRITE(FILE_NAME,112)"OUTPUT/Iteration_Reports/Iteration_Report_P",Rank,"_I",Iteration,".out"
    CALL OPEN_NEW_FILE( FILE_NAME, ITER_REPORT_FILE_ID )



    WRITE(ITER_REPORT_FILE_ID,109)"                      Report for Iteration ",Iteration," by Process ",Rank
    WRITE(ITER_REPORT_FILE_ID,'(A)')"-----------------------------------------------------------------------------------------"
    WRITE(ITER_REPORT_FILE_ID,'(A)')" "

END IF


IF ( FRAME_REPORT_FLAG == 1 ) THEN

    WRITE(FRAME_REPORT_FILE_ID,'(A)')"-----------------------------------------------------------------------------------------"
    WRITE(FRAME_REPORT_FILE_ID,109)"                      Report for Iteration ",Iteration," by Process ",Rank
    WRITE(FRAME_REPORT_FILE_ID,'(A)')"-----------------------------------------------------------------------------------------"
    WRITE(FRAME_REPORT_FILE_ID,'(3/)')

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




 !+401+############################################################################!
!                                                                                   !
!                           OUTPUT_FINAL_RESULTS                                    !
!                                                                                   !
 !#################################################################################!
SUBROUTINE OUTPUT_FINAL_RESULTS()


INTEGER                                                     ::  NUM_SAMPLES


CHARACTER(LEN = 50), DIMENSION(:), ALLOCATABLE              ::  Filenames
INTEGER, DIMENSION(:), ALLOCATABLE                          ::  File_IDs
INTEGER                                                     ::  Num_Files

CHARACTER(LEN = 50)                                         ::  filenamea,          &
                                                                filenameb,          &
                                                                filenamec,          &
                                                                filenamed,          &
                                                                filenamee

INTEGER                                                     ::  file_ida,           &
                                                                file_idb,           &
                                                                file_idc,           &
                                                                file_idd

REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Output_re,          &
                                                                Output_rc,          &
                                                                Output_dr

REAL(KIND = idp)                                            ::  Return_Psi,         &
                                                                Return_AlphaPsi,    &
                                                                Return_Beta1,       &
                                                                Return_Beta2,       &
                                                                Return_Beta3


REAL(KIND = idp)                                            ::  deltar, r,          &
                                                                theta, phi,         &
                                                                Analytic_Val,       &
                                                                Solver_Val,         &
                                                                Solver_Valb,        &
                                                                Error_Val

REAL(KIND = idp)                                            ::  csqr
INTEGER                                                     ::  i, j, k

INTEGER                                                     ::  NUM_THETA_RAYS,     &
                                                                NUM_PHI_RAYS,       &
                                                                NUM_RADIAL_SAMPLES

REAL(KIND = idp)                                            ::  DELTA_THETA,        &
                                                                THETA_VAL

REAL(KIND = idp)                                            ::  DELTA_PHI,        &
                                                                PHI_VAL

REAL(KIND = idp), DIMENSION(:,:,:), ALLOCATABLE             ::  Lapse_Holder,       &
                                                                ConForm_Holder
REAL(KIND = idp), DIMENSION(:,:,:,:), ALLOCATABLE           ::  Shift_Holder
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  R_Holder,           &
                                                                T_Holder,           &
                                                                P_Holder


110 FORMAT (11X,A1,16X,A18,9X,A13,10X,A18,10X,A10)                      !!! Output Header

111 FORMAT (11X,A1,24X,A3,19X,A8,15X,A11,14X,A11,14X,A11)             !!! Output Header for Results file
112 FORMAT (11X,A1,16X,A18,9x,A14)                                      !!! Output Header for Analytic Solution file

113 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15)               !!! Output
114 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15)    !!! Output for Results file
115 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15)                                     !!! Output for Analytic Solution file


IF ( RESULTS_OUTPUT_FLAG == 1 ) THEN
    NUM_SAMPLES = 1000

    ALLOCATE( Output_re(0:NUM_SAMPLES) )
    ALLOCATE( Output_rc(1:NUM_SAMPLES) )
    ALLOCATE( Output_dr(1:NUM_SAMPLES) )

    IF ( DRIVER_TEST_NUMBER == 1 ) THEN
        ! Open Results File
        file_ida = 42
        filenamea = "OUTPUT/Results.out"
        CALL OPEN_NEW_FILE( filenamea, file_ida)
        WRITE(file_ida,111)"r","Psi","AlphaPsi","Beta1 Value","Beta2 Value","Beta3 Value"

        ! Open Solution File
        file_idb = 43
        filenameb = "OUTPUT/Solution.out"
        CALL OPEN_NEW_FILE( filenameb, file_idb)
        WRITE(file_idb,112)"r","Analytic Potential","Beta1 Solution"


        ! Set Output locations
        theta = pi/2.0_idp
        phi = pi/2.0_idp

        CALL Create_Logarithmic_1D_Mesh( R_INNER, R_OUTER, NUM_SAMPLES,        &
                                        output_re, output_rc, output_dr                 )



        DO i = 1,NUM_SAMPLES

            CALL Calc_3D_Values_At_Location( output_rc(i), theta, phi,                   &
                                             Return_Psi, Return_AlphaPsi,                &
                                             Return_Beta1, Return_Beta2, Return_Beta3    )

            Analytic_Val = Analytic_Solution(output_rc(i),theta,phi)
            Solver_Val = 2.0_idp*csqr*(1.0_idp - Return_Psi)
            Solver_Valb = 2.0_idp*csqr*(Return_AlphaPsi - 1.0_idp)
            Error_Val = ABS((Analytic_Val - Solver_Val)/Analytic_Val)


           WRITE(42,114) output_rc(i), Return_Psi, Return_AlphaPsi, Return_Beta1,Return_Beta2,Return_Beta3
           WRITE(43,115) output_rc(i),Analytic_Val, Shift_Solution(output_rc(i),rlocs,NUM_R_ELEMENTS)



        END DO


        ! Close Files
        CLOSE( Unit = file_ida)
        CLOSE( Unit = file_idb)


    ELSE IF ( DRIVER_TEST_NUMBER == 5 ) THEN

        ! Open Results File
        file_ida = 42
        WRITE(filenamea,'(A31,I5.5,A4)')"OUTPUT/CHIMERA_RESULTS/Results_",DRIVER_FRAME,".out"
        CALL OPEN_NEW_FILE( filenamea, file_ida)
        WRITE(file_ida,111)"r","Psi","AlphaPsi","Beta1 Value","Beta2 Value","Beta3 Value"

        ! Open Solution File
        file_idb = 43
        WRITE(filenameb,'(A32,I5.5,A4)')"OUTPUT/CHIMERA_RESULTS/Solution_",DRIVER_FRAME,".out"
        CALL OPEN_NEW_FILE( filenameb, file_idb)
        WRITE(file_idb,112)"r","Analytic Potential","Beta1 Solution"


        ! Set Output locations
        theta = pi/2.0_idp
        phi = pi/2.0_idp

        CALL Create_Logarithmic_1D_Mesh( R_INNER, R_OUTER, NUM_SAMPLES,        &
                                        output_re, output_rc, output_dr                 )



        DO i = 1,NUM_SAMPLES

            CALL Calc_3D_Values_At_Location( output_rc(i), theta, phi,                   &
                                             Return_Psi, Return_AlphaPsi,                &
                                             Return_Beta1, Return_Beta2, Return_Beta3    )

            Analytic_Val = Analytic_Solution(output_rc(i),theta,phi)
            Solver_Val = 2.0_idp*csqr*(1.0_idp - Return_Psi)
            Solver_Valb = 2.0_idp*csqr*(Return_AlphaPsi - 1.0_idp)
            Error_Val = ABS((Analytic_Val - Solver_Val)/Analytic_Val)


           WRITE(42,114) output_rc(i), Return_Psi, Return_AlphaPsi, Return_Beta1,Return_Beta2,Return_Beta3
           WRITE(43,115) output_rc(i),Analytic_Val, Shift_Solution(output_rc(i),rlocs,NUM_R_ELEMENTS)



        END DO


        ! Close Files
        CLOSE( Unit = file_ida)
        CLOSE( Unit = file_idb)







    ELSE IF ( DRIVER_TEST_NUMBER == 2 ) THEN

        Num_Files = 9

        ALLOCATE( Filenames(1:Num_Files) )
        ALLOCATE( File_IDs(1:Num_Files) )


        WRITE(Filenames(1),'(A37,I5.5,A4)')"OUTPUT/CHIMERA_RESULTS/Results_Lapse_",DRIVER_FRAME,".out"
        WRITE(Filenames(2),'(A41,I5.5,A4)')"OUTPUT/CHIMERA_RESULTS/Results_ConFactor_",DRIVER_FRAME,".out"
        WRITE(Filenames(3),'(A37,I5.5,A4)')"OUTPUT/CHIMERA_RESULTS/Results_Beta1_",DRIVER_FRAME,".out"
        WRITE(Filenames(4),'(A37,I5.5,A4)')"OUTPUT/CHIMERA_RESULTS/Results_Beta2_",DRIVER_FRAME,".out"
        WRITE(Filenames(5),'(A37,I5.5,A4)')"OUTPUT/CHIMERA_RESULTS/Results_Beta3_",DRIVER_FRAME,".out"
        WRITE(Filenames(6),'(A45)')"OUTPUT/CHIMERA_RESULTS/Results_Dimensions.out"
        WRITE(Filenames(7),'(A32,I5.5,A4)')"OUTPUT/CHIMERA_RESULTS/R_VALUES_",DRIVER_FRAME,".out"
        WRITE(Filenames(8),'(A32,I5.5,A4)')"OUTPUT/CHIMERA_RESULTS/T_VALUES_",DRIVER_FRAME,".out"
        WRITE(Filenames(9),'(A32,I5.5,A4)')"OUTPUT/CHIMERA_RESULTS/P_VALUES_",DRIVER_FRAME,".out"



        File_IDs = [(141 + i, i=1,Num_Files)]
        DO i = 1,Num_Files
            CALL OPEN_NEW_FILE( Filenames(i), File_IDs(i) )
        END DO





        ! Create Output Spacing
        NUM_PHI_RAYS = 1
        IF ( NUM_PHI_RAYS == 1 ) THEN
            DELTA_PHI = pi/2.0_idp
!            DELTA_PHI = 0.0_idp
        ELSE
            DELTA_PHI = 2.0_idp*pi/(NUM_PHI_RAYS-1)
        END IF

        NUM_THETA_RAYS = 1
        DELTA_THETA = pi/(NUM_THETA_RAYS+1)

        NUM_RADIAL_SAMPLES = 1000
        CALL Create_Logarithmic_1D_Mesh( R_INNER, R_OUTER, NUM_RADIAL_SAMPLES,  &
                                         output_re, output_rc, output_dr        )

        ALLOCATE( Lapse_Holder(1:NUM_PHI_RAYS, 1:NUM_THETA_RAYS, 1:NUM_RADIAL_SAMPLES) )
        ALLOCATE( ConForm_Holder(1:NUM_PHI_RAYS, 1:NUM_THETA_RAYS, 1:NUM_RADIAL_SAMPLES) )
        ALLOCATE( Shift_Holder(1:3,1:NUM_PHI_RAYS, 1:NUM_THETA_RAYS, 1:NUM_RADIAL_SAMPLES) )
        ALLOCATE( R_Holder(1:NUM_RADIAL_SAMPLES) )
        ALLOCATE( T_Holder(1:NUM_THETA_RAYS) )
        ALLOCATE( P_Holder(1:NUM_PHI_RAYS) )



        ! Calculate Output
        DO k = 1,NUM_PHI_RAYS
            DO j = 1,NUM_THETA_RAYS
                DO i = 1,NUM_RADIAL_SAMPLES

                    PHI_VAL = k*DELTA_PHI
                    THETA_VAL = (j-1)*DELTA_THETA

                    CALL Calc_3D_Values_At_Location( output_rc(i), THETA_VAL, PHI_VAL,                   &
                                                     Return_Psi, Return_AlphaPsi,                &
                                                     Return_Beta1, Return_Beta2, Return_Beta3    )

                    Lapse_Holder(k,j,i) = Return_AlphaPsi/Return_Psi
                    ConForm_Holder(k,j,i) = Return_Psi
                    Shift_Holder(1:3,k,j,i) = (/ Return_Beta1, Return_Beta2, Return_Beta3 /)

                    R_Holder(i) = output_rc(i)
                    T_Holder(j) = THETA_VAL
                    P_Holder(k) = PHI_VAL


                END DO ! i Loop
            END DO ! j Loop
        END DO ! k Loop


        ! Write Output Location Files
        WRITE(File_IDs(7),*)R_Holder
        WRITE(File_IDs(8),*)T_Holder
        WRITE(File_IDs(9),*)P_Holder


        ! Write Output Value Files
        DO k = 1,NUM_PHI_RAYS
            DO j = 1,NUM_THETA_RAYS


                WRITE(File_IDs(1),*)Lapse_Holder(k,j,:)
                WRITE(File_IDs(2),*)ConForm_Holder(k,j,:)
                WRITE(File_IDs(3),*)Shift_Holder(1,k,j,:)
                WRITE(File_IDs(4),*)Shift_Holder(2,k,j,:)
                WRITE(File_IDs(5),*)Shift_Holder(3,k,j,:)

            END DO ! j Loop
        END DO ! k Loop



        ! Close Files
        DO i = 1,Num_Files
            CLOSE( Unit = File_IDs(i))
        END DO


    END IF


END IF

END SUBROUTINE OUTPUT_FINAL_RESULTS






 !+301+############################################################################!
!                                                                                   !
!                     OPEN_FRAME_REPORT_FILE                                          !
!                                                                                   !
 !#################################################################################!
SUBROUTINE OPEN_FRAME_REPORT_FILE(Frame)

INTEGER, INTENT(IN)                                     ::  Frame

CHARACTER(LEN = 53)                                     ::  FILE_NAME
INTEGER                                                 ::  istat
LOGICAL                                                 ::  FLAG, OK

109 FORMAT (A,I2.2)



IF ( myID == 0 ) THEN
    IF ( FRAME_REPORT_FLAG == 1 ) THEN

        WRITE(FILE_NAME,'(A)')"OUTPUT/Iteration_Reports/Frame_Report_SCRATCH.out"
        CALL OPEN_NEW_FILE( FILE_NAME, FRAME_REPORT_FILE_ID )

    END IF ! WRITE_REPORT_FLAGS
END IF ! myID == 0


END SUBROUTINE OPEN_FRAME_REPORT_FILE




!+501+##########################################################################!
!                                                                               !
!                   OUTPUT_ITERATION_REPORT                                     !
!                                                                               !
!###############################################################################!
SUBROUTINE OUTPUT_FRAME_REPORT(FRAME)

INTEGER, INTENT(IN)                 :: FRAME

INTEGER                                         ::  FILE_ID, FILE_IDb
INTEGER                                         ::  i
REAL(KIND = idp)                                ::  r, theta, phi, deltar
REAL(KIND = idp)                                ::  Analytic_Val, Solver_Val
REAL(KIND = idp)                                ::  Return_Psi, Return_AlphaPsi
REAL(KIND = idp)                                ::  Return_Beta1, Return_Beta2, Return_Beta3
REAL(KIND = idp)                                ::  PsiPot_Val, AlphaPsiPot_Val

CHARACTER(LEN = 53)                             ::  FILE_NAME

CHARACTER(LEN = 300)                            ::  Line
INTEGEr                                         ::  io_stat


120 FORMAT (12X,A61)
121 FORMAT (A1)
122 FORMAT (A41,I2.2)
123 FORMAT (12X,A38,ES22.15)

109 FORMAT (A,I2.2,A,I2.2)
110 FORMAT (11X,A1,16X,A18,9X,A13,10X,A18,10X,A11,14X,A11,14X,A11)
111 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15)
112 FORMAT (A43,I2.2,A2,I2.2,A4)

132 FORMAT (A38,I2.2,A4)

142 FORMAT (16X,A,10X,A)
143 FORMAT (19X,I2.2,15X,ES22.15)


WRITE(FILE_NAME,132)"OUTPUT/Iteration_Reports/Frame_Report_",Frame,".out"
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
        WRITE(FILE_ID,123)"                    Initialize Time : ",FRAME_TIME_TABLE(1)
        WRITE(FILE_ID,123)" Input/Communicate Source Data Time : ",FRAME_TIME_TABLE(2)
        WRITE(FILE_ID,123)"     Input Boundary Conditions Time : ",FRAME_TIME_TABLE(3)
        WRITE(FILE_ID,123)"        CFA_3D_Apply_BCs_Part1 Time : ",FRAME_TIME_TABLE(4)
        WRITE(FILE_ID,120)"-------------------------------------------------------------"
        WRITE(FILE_ID,123)" ||     Calc_3D_Current_Values Time : ",FRAME_TIME_TABLE(5)
        WRITE(FILE_ID,123)" ||    CREATE_3D_SubJcbn_Terms Time : ",FRAME_TIME_TABLE(6)
        WRITE(FILE_ID,123)" ||       CREATE_3D_RHS_VECTOR Time : ",FRAME_TIME_TABLE(7)
        WRITE(FILE_ID,123)"\  /     CREATE_3D_JCBN_MATRIX Time : ",FRAME_TIME_TABLE(8)
        WRITE(FILE_ID,120)"-\/ ---------------------------------------------------------"
        WRITE(FILE_ID,123)"CREATE_3D_NONLAPLACIAN_STF_MAT Time : ",FRAME_TIME_TABLE(9)
        WRITE(FILE_ID,123)"REDUCE_3D_NONLAPLACIAN_STF_MAT Time : ",FRAME_TIME_TABLE(10)
        WRITE(FILE_ID,123)"FINISH_3D_NONLAPLACIAN_STF_MAT Time : ",FRAME_TIME_TABLE(11)
        WRITE(FILE_ID,123)"          FINISH_3D_RHS_VECTOR Time : ",FRAME_TIME_TABLE(12)
        WRITE(FILE_ID,123)"        CFA_3D_Apply_BCs_Part2 Time : ",FRAME_TIME_TABLE(13)
        WRITE(FILE_ID,123)"                    CFA_Solver Time : ",FRAME_TIME_TABLE(14)
        WRITE(FILE_ID,123)"        CFA_Coefficient_Update Time : ",FRAME_TIME_TABLE(15)
        WRITE(FILE_ID,123)"   CFA_Coefficient_Share_PETSc Time : ",FRAME_TIME_TABLE(16)
        WRITE(FILE_ID,123)"         CFA_Convergence_Check Time : ",FRAME_TIME_TABLE(17)
        WRITE(FILE_ID,123)"               Total Iteration Time : ",FRAME_TIME_TABLE(18)
        WRITE(FILE_ID,123)"             Poseidon_Dist_Sol Time : ",FRAME_TIME_TABLE(19)
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




        deltar = ( R_OUTER - R_INNER )/ REAL(ITER_REPORT_NUM_SAMPLES, KIND = idp)

        DO i = 0,ITER_REPORT_NUM_SAMPLES

            r = i*deltar + R_INNER
            theta = pi/2.0_idp
            phi = pi/2.0_idp


            CALL Calc_3D_Values_At_Location( r, theta, phi,                              &
                                            Return_Psi, Return_AlphaPsi,                &
                                            Return_Beta1, Return_Beta2, Return_Beta3    )

            ! Determine the Newtonian Potential at the location r, theta, phi
            Analytic_Val = Analytic_Solution(r,theta,phi)


            ! AlphaPsi_to_Pot   =   2*C_Square*(AlphaPsi - 1)
            ! Psi_to_Pot        =   2*C_Square*(1 - Psi)

            ! Calculate Conformal Factor value from Newtonian Potential
            PsiPot_Val = 2.0_idp*C_Square*(1.0_idp - Return_Psi)

            ! Calculate the product of the Conformal Factor and Lapse Function from Newtonian Potential
            AlphaPsiPot_Val = 2.0_idp*C_Square*(Return_AlphaPsi - 1.0_idp)

            ! Write Results to File
            WRITE(FILE_ID,111) r,Analytic_Val,PsiPot_Val,AlphaPsiPot_Val,Return_Beta1,Return_Beta2,Return_Beta3

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

 !+202+############################################################################!
!                                                                                   !
!                    CLOSE_RUN_REPORT_FILE                                          !
!                                                                                   !
 !#################################################################################!
SUBROUTINE CLOSE_FRAME_REPORT_FILE()


CLOSE( UNIT = FRAME_REPORT_FILE_ID, STATUS='delete' )


END SUBROUTINE CLOSE_FRAME_REPORT_FILE




!+501+##########################################################################!
!                                                                               !
!                   OUTPUT_ITERATION_REPORT                                     !
!                                                                               !
!###############################################################################!
SUBROUTINE OUTPUT_RUN_REPORT()


INTEGER                                         ::  FILE_ID
INTEGER                                         ::  i
REAL(KIND = idp)                                ::  r, theta, phi, deltar
REAL(KIND = idp)                                ::  Analytic_Val, Solver_Val
REAL(KIND = idp)                                ::  Return_Psi, Return_AlphaPsi
REAL(KIND = idp)                                ::  Return_Beta1, Return_Beta2, Return_Beta3
REAL(KIND = idp)                                ::  PsiPot_Val, AlphaPsiPot_Val

120 FORMAT (A61)
121 FORMAT (A1)
122 FORMAT (A41,I2.2)
123 FORMAT (A38,ES22.15)

109 FORMAT (A,I2.2,A,I2.2)
110 FORMAT (11X,A1,18X,A13,10X,A18,10X,A11,14X,A11,14X,A11)
111 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15)
112 FORMAT (A43,I2.2,A2,I2.2,A4)

FILE_ID = RUN_REPORT_FILE_ID



! Write Timetable to File
IF ( RUN_REPORT_FLAG == 1 ) THEN

    WRITE(FILE_ID,'(A)')"                                Average Timing Results"
    WRITE(FILE_ID,'(A)')"            ============================================================="
    WRITE(FILE_ID,'(A)')" "
    WRITE(FILE_ID,123)"                    Initialize Time : ",RUN_TIME_TABLE(1)
    WRITE(FILE_ID,123)" Input/Communicate Source Data Time : ",RUN_TIME_TABLE(2)
    WRITE(FILE_ID,123)"     Input Boundary Conditions Time : ",RUN_TIME_TABLE(3)
    WRITE(FILE_ID,123)"        CFA_3D_Apply_BCs_Part1 Time : ",RUN_TIME_TABLE(4)
    WRITE(FILE_ID,120)"-------------------------------------------------------------"
    WRITE(FILE_ID,123)" ||     Calc_3D_Current_Values Time : ",RUN_TIME_TABLE(5)
    WRITE(FILE_ID,123)" ||    CREATE_3D_SubJcbn_Terms Time : ",RUN_TIME_TABLE(6)
    WRITE(FILE_ID,123)" ||       CREATE_3D_RHS_VECTOR Time : ",RUN_TIME_TABLE(7)
    WRITE(FILE_ID,123)"\  /     CREATE_3D_JCBN_MATRIX Time : ",RUN_TIME_TABLE(8)
    WRITE(FILE_ID,120)"-\/ ---------------------------------------------------------"
    WRITE(FILE_ID,123)"CREATE_3D_NONLAPLACIAN_STF_MAT Time : ",RUN_TIME_TABLE(9)
    WRITE(FILE_ID,123)"REDUCE_3D_NONLAPLACIAN_STF_MAT Time : ",RUN_TIME_TABLE(10)
    WRITE(FILE_ID,123)"FINISH_3D_NONLAPLACIAN_STF_MAT Time : ",RUN_TIME_TABLE(11)
    WRITE(FILE_ID,123)"          FINISH_3D_RHS_VECTOR Time : ",RUN_TIME_TABLE(12)
    WRITE(FILE_ID,123)"        CFA_3D_Apply_BCs_Part2 Time : ",RUN_TIME_TABLE(13)
    WRITE(FILE_ID,123)"                    CFA_Solver Time : ",RUN_TIME_TABLE(14)
    WRITE(FILE_ID,123)"        CFA_Coefficient_Update Time : ",RUN_TIME_TABLE(15)
    WRITE(FILE_ID,123)"   CFA_Coefficient_Share_PETSc Time : ",RUN_TIME_TABLE(16)
    WRITE(FILE_ID,123)"         CFA_Convergence_Check Time : ",RUN_TIME_TABLE(17)
    WRITE(FILE_ID,123)"               Total Iteration Time : ",RUN_TIME_TABLE(18)
    WRITE(FILE_ID,123)"             Poseidon_Dist_Sol Time : ",RUN_TIME_TABLE(19)
    WRITE(FILE_ID,120)"============================================================="
    WRITE(FILE_ID,121)" "
    WRITE(FILE_ID,121)" "
    WRITE(FILE_ID,121)" "
    WRITE(FILE_ID,121)" "


    WRITE(FILE_ID,'(A)')"                                 Convergence Results"
    WRITE(FILE_ID,'(A)')"            ============================================================="
    WRITE(FILE_ID,'(A)')""
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
    WRITE(FILE_ID,'(A)')" "
    WRITE(FILE_ID,'(A)')" "


    WRITE(FILE_ID,'(A)')"                                 Sample of Final Results"
    WRITE(FILE_ID,'(A)')"            ============================================================="
    WRITE(FILE_ID,'(A)')" "
    WRITE(FILE_ID,110)"r","Psi Potential","AlphaPsi Potential","Beta Value1","Beta Value2","Beta Value3"
END IF


! Write Timetable to Screen
IF (( WRITE_REPORT_FLAG == 1) .OR. (WRITE_REPORT_FLAG == 3) ) THEN
    WRITE(*,'(A)')" "
    WRITE(*,'(A)')" "
    WRITE(*,'(A)')"            ============================================================="
    WRITE(*,'(A)')"                                Average Timing Results"
    WRITE(*,'(A)')"            ============================================================="
    WRITE(*,'(A)')" "
    WRITE(*,123)"                    Initialize Time : ",RUN_TIME_TABLE(1)
    WRITE(*,123)" Input/Communicate Source Data Time : ",RUN_TIME_TABLE(2)
    WRITE(*,123)"     Input Boundary Conditions Time : ",RUN_TIME_TABLE(3)
    WRITE(*,123)"        CFA_3D_Apply_BCs_Part1 Time : ",RUN_TIME_TABLE(4)
    WRITE(*,120)"-------------------------------------------------------------"
    WRITE(*,123)" ||     Calc_3D_Current_Values Time : ",RUN_TIME_TABLE(5)
    WRITE(*,123)" ||    CREATE_3D_SubJcbn_Terms Time : ",RUN_TIME_TABLE(6)
    WRITE(*,123)" ||       CREATE_3D_RHS_VECTOR Time : ",RUN_TIME_TABLE(7)
    WRITE(*,123)"\  /     CREATE_3D_JCBN_MATRIX Time : ",RUN_TIME_TABLE(8)
    WRITE(*,120)"-\/ ---------------------------------------------------------"
    WRITE(*,123)"CREATE_3D_NONLAPLACIAN_STF_MAT Time : ",RUN_TIME_TABLE(9)
    WRITE(*,123)"REDUCE_3D_NONLAPLACIAN_STF_MAT Time : ",RUN_TIME_TABLE(10)
    WRITE(*,123)"FINISH_3D_NONLAPLACIAN_STF_MAT Time : ",RUN_TIME_TABLE(11)
    WRITE(*,123)"          FINISH_3D_RHS_VECTOR Time : ",RUN_TIME_TABLE(12)
    WRITE(*,123)"        CFA_3D_Apply_BCs_Part2 Time : ",RUN_TIME_TABLE(13)
    WRITE(*,123)"                    CFA_Solver Time : ",RUN_TIME_TABLE(14)
    WRITE(*,123)"        CFA_Coefficient_Update Time : ",RUN_TIME_TABLE(15)
    WRITE(*,123)"   CFA_Coefficient_Share_PETSc Time : ",RUN_TIME_TABLE(16)
    WRITE(*,123)"         CFA_Convergence_Check Time : ",RUN_TIME_TABLE(17)
    WRITE(*,123)"               Total Iteration Time : ",RUN_TIME_TABLE(18)
    WRITE(*,123)"             Poseidon_Dist_Sol Time : ",RUN_TIME_TABLE(19)
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




deltar = ( R_OUTER - R_INNER )/ REAL(ITER_REPORT_NUM_SAMPLES, KIND = idp)
DO i = 0,ITER_REPORT_NUM_SAMPLES

    r = i*deltar + R_INNER
    theta = pi/2.0_idp
    phi = pi/2.0_idp


    CALL Calc_3D_Values_At_Location( r, theta, phi,                              &
                                    Return_Psi, Return_AlphaPsi,                &
                                    Return_Beta1, Return_Beta2, Return_Beta3    )



    ! Calculate Conformal Factor value from Newtonian Potential
    PsiPot_Val = 2.0_idp*C_Square*(1.0_idp - Return_Psi)

    ! Calculate the product of the Conformal Factor and Lapse Function from Newtonian Potential
    AlphaPsiPot_Val = 2.0_idp*C_Square*(Return_AlphaPsi - 1.0_idp)


    ! Write Results to Screen
    IF (( WRITE_REPORT_FLAG == 1) .OR. (WRITE_REPORT_FLAG == 3) ) THEN
        WRITE(*,111) r,PsiPot_Val,AlphaPsiPot_Val,Return_Beta1,Return_Beta2,Return_Beta3
    END IF

    ! Write Results to File
    IF ( RUN_REPORT_FLAG == 1 ) THEN
        WRITE(FILE_ID,111) r,PsiPot_Val,AlphaPsiPot_Val,Return_Beta1,Return_Beta2,Return_Beta3
    END IF

END DO



END SUBROUTINE OUTPUT_RUN_REPORT



!+101+###########################################################################!
!                                                                                !
!                  OUTPUT_PETSC_REPORT                                           !
!                                                                                !
!################################################################################!
SUBROUTINE OUTPUT_PETSC_REPORT( matset_time, solve_time, Iter_Count,   &
                                rtol, abstol, dtol, maxits             )

REAL(KIND = idp),   INTENT(IN)                          ::  matset_time, solve_time
INTEGER,            INTENT(IN)                          ::  Iter_Count, maxits
REAL(KIND = idp),   INTENT(IN)                          ::  rtol, abstol, dtol

INTEGER, DIMENSION(0:1)                                 ::  FILE_ID
INTEGER                                                 ::  i

FILE_ID = -1
IF ( FRAME_REPORT_FLAG == 1 ) THEN

    FILE_ID(0) = FRAME_REPORT_FILE_ID

END IF

IF (( WRITE_REPORT_FLAG == 2) .OR. (WRITE_REPORT_FLAG == 3) ) THEN
    FILE_ID(1) = ITER_REPORT_FILE_ID
END IF


DO i = 0,1
    IF ( FILE_ID(i) .NE. -1 ) THEN

        WRITE(FILE_ID(i),'(A)')"                                     PETSc Results"
        WRITE(FILE_ID(i),'(A)')"            ============================================================="
        WRITE(FILE_ID(i),'(A)')" "
        WRITE(FILE_ID(i),'(A,ES22.15)')"Set Matrix Values Time ",matset_time
        WRITE(FILE_ID(i),'(A)')"Solver Type : KSPPREONLY, PCILU"
        WRITE(FILE_ID(i),'(A,ES22.15)')"Solve Time ",solve_time
        WRITE(FILE_ID(i),'(A,I2.2)')"Solve Iterations :",Iter_Count
        WRITE(FILE_ID(i),'(A,ES22.15,ES22.15,ES22.15,I2.2)')"Solve Tolerances :",rtol, abstol, dtol, maxits
        WRITE(FILE_ID(i),'(A)')" "
        WRITE(FILE_ID(i),'(A)')" "
        WRITE(FILE_ID(i),'(A)')" "

    END IF
END DO


END SUBROUTINE OUTPUT_PETSC_REPORT







END MODULE IO_Functions_Module
