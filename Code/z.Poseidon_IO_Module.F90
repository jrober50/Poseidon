   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_IO_Module                                                           !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   OPEN_ITER_REPORT_FILE                                               !##!
!##!    +102+   OUTPUT_ITER_TIMETABLE                                               !##!
!##!    +103+   CLOSE_ITER_REPORT_FILE                                              !##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!    +301+   OPEN_FRAME_REPORT_FILE                                              !##!
!##!    +302+   OUTPUT_FRAME_REPORT                                                 !##!
!##!    +303+   CLOSE_FRAME_REPORT_FILE                                             !##!
!##!                                                                                !##!
!##!    +401+   OUTPUT_FINAL_RESULTS                                                !##!
!##!    +402+   OUTPUT_PETSC_REPORT                                                 !##!
!##!    +403+   OUTPUT_STF_ELEM_BLOCK_MATRIX                                        !##!
!##!                                                                                !##!
!##!    +501+   CLOCK_IN                                                            !##!
!##!    +502+   OPEN_NEW_FILE                                                       !##!
!##!                                                                                !##!
!##!    +601+   OUTPUT_POSEIDON_SOURCES                                             !##!
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

USE Units_Module, &
                    ONLY :  C_Square,       &
                            Gram,           &
                            Centimeter,     &
                            Kilometer,      &
                            Erg,            &
                            Second,         &
                            GravPot_Units,  &
                            Shift_Units

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
                            WRITE_RESULTS_R_SAMPS,                          &
                            WRITE_RESULTS_T_SAMPS,                          &
                            WRITE_RESULTS_P_SAMPS,                          &
                            ITER_REPORT_NUM_SAMPLES,                        &
                            RUN_REPORT_FILE_ID,                             &
                            ITER_REPORT_FILE_ID,                            &
                            FRAME_REPORT_FILE_ID,                           &
                            WRITE_REPORT_FLAG,                              &
                            WRITE_SOURCES_FLAG,                             &
                            OUTPUT_RHS_VECTOR_FLAG,                         &
                            Poseidon_Frame

USE DRIVER_Parameters,  &
                    ONLY :  Potential_Solution,                              &
                            Shift_Solution,                                 &
                            myID,                                           &
                            SOURCE_OUTPUT_FLAG,                             &
                            RESULTS_OUTPUT_FLAG,                            &
                            RUN_REPORT_FLAG,                                &
                            FRAME_REPORT_FLAG,                              &
                            DRIVER_TEST_NUMBER,                             &
                            Driver_Frame,                                   &
                            SELFSIM_T

USE Poseidon_Variables_Module, &
                    ONLY :  NUM_R_NODES,                                    &
                            NUM_R_ELEMENTS,                                 &
                            RHS_Vector,                                     &
                            Coefficient_Vector,                             &
                            Update_Vector,                                  &
                            rlocs, tlocs, plocs,                            &
                            R_INNER, R_OUTER,                               &
                            INT_R_LOCATIONS,                                &
                            INT_T_LOCATIONS,                                &
                            INT_P_LOCATIONS,                                &
                            VAR_DIM,                                        &
                            NUM_OFF_DIAGONALS,                              &
                            SUBSHELL_PROB_DIM,                              &
                            Block_PROB_DIM,                                 &
                            BLOCK_RHS_VECTOR,                               &
                            BLOCK_ELEM_STF_MATVEC,                          &
                            ELEM_PROB_DIM_SQR,                              &
                            Total_Run_Iters,                                &
                            Num_Timer_Calls,                                &
                            ITER_TIME_TABLE,                                &
                            FRAME_TIME_TABLE,                               &
                            FRAME_CONVERGENCE_TABLE,                        &
                            RUN_TIME_TABLE,                                 &
                            BLOCK_STF_MAT
 


USE Poseidon_Math_Functions_Module, &
                    ONLY :  Lagrange_Poly

USE Poseidon_Calculate_Results_Module, &
                    ONLY :  Calc_3D_Values_At_Location

USE Poseidon_Quadrature_Module, &
                    ONLY :  Initialize_LG_Quadrature_Locations


USE Poseidon_Mesh_Module, &
                    ONLY :  Create_Uniform_1D_Mesh,                         &
                            Create_Logarithmic_1D_Mesh,                     &
                            Create_Split_1D_Mesh

USE Poseidon_IO_Parameters, &
                    ONLY :  Poseidon_Reports_Dir,                           &
                            Poseidon_IterReports_Dir,                       &
                            Poseidon_Objects_Dir,                           &
                            Poseidon_LinSys_Dir,                            &
                            Poseidon_Results_Dir,                           &
                            Poseidon_Sources_Dir




IMPLICIT NONE


CHARACTER(LEN = 20), PARAMETER    :: Filename_Format_A = "(A,A)"
CHARACTER(LEN = 20), PARAMETER    :: Filename_Format_B = "(A,A,I2.2,A,I2.2,A)"


!*F&S*==========================================!
!                                               !
!           Functions & Subroutines             !
!                                               !
!===============================================!
CONTAINS



 !+101+############################################################################!
!                                                                                   !
!                     OPEN_ITER_REPORT_FILE                                         !
!                                                                                   !
 !#################################################################################!
SUBROUTINE OPEN_ITER_REPORT_FILE(Iteration, Rank)

INTEGER, INTENT(IN)                                     ::  Iteration, Rank

CHARACTER(LEN = 70)                                     ::  FILE_NAME
INTEGER                                                 ::  istat
LOGICAL                                                 ::  FLAG, OK

109 FORMAT (A,I2.2,A,I2.2)
112 FORMAT (A,I2.2,A,I2.2,A)
113 FORMAT (A,A,I2.2,A,I2.2,A)

IF (( WRITE_REPORT_FLAG == 2) .OR. (WRITE_REPORT_FLAG == 3) ) THEN

    WRITE(FILE_NAME,113) Poseidon_IterReports_Dir,"Iteration_Report_P",Rank,"_I",Iteration,".out"
    CALL OPEN_NEW_FILE( FILE_NAME, ITER_REPORT_FILE_ID, 100 )



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





 !+102+############################################################################!
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

     
    WRITE(FILENAME,'(A,I2.2,A)')'OUTPUT/Timetable_',Ident,'.out'

    OPEN(unit = FILE_ID,file = FILENAME)
    CALL OPEN_NEW_FILE( FILENAME, File_ID, 50 )

    WRITE(FILE_ID,110)"============================================================="
    WRITE(FILE_ID,111)" "
    WRITE(FILE_ID,112)"                 Time Table for Process ",Ident
    WRITE(FILE_ID,111)" "
    WRITE(FILE_ID,110)"============================================================="
    WRITE(FILE_ID,113)"                    Initialize Time : ",ITER_TIME_TABLE(1)
    WRITE(FILE_ID,113)" Input/Communicate Source Data Time : ",ITER_TIME_TABLE(2)
    WRITE(FILE_ID,113)"     Input Boundary Conditions Time : ",ITER_TIME_TABLE(3)
    WRITE(FILE_ID,113)"        CFA_3D_Apply_BCs_Part1 Time : ",ITER_TIME_TABLE(4)
    WRITE(FILE_ID,110)"-------------------------------------------------------------"
    WRITE(FILE_ID,113)" ||     Calc_3D_Current_Values Time : ",ITER_TIME_TABLE(5)
    WRITE(FILE_ID,113)" ||    CREATE_3D_SubJcbn_Terms Time : ",ITER_TIME_TABLE(6)
    WRITE(FILE_ID,113)" ||       CREATE_3D_RHS_VECTOR Time : ",ITER_TIME_TABLE(7)
    WRITE(FILE_ID,113)"\  /     CREATE_3D_JCBN_MATRIX Time : ",ITER_TIME_TABLE(8)
    WRITE(FILE_ID,110)"-\/ ---------------------------------------------------------"
    WRITE(FILE_ID,113)"CREATE_3D_NONLAPLACIAN_STF_MAT Time : ",ITER_TIME_TABLE(9)
    WRITE(FILE_ID,113)"REDUCE_3D_NONLAPLACIAN_STF_MAT Time : ",ITER_TIME_TABLE(10)
    WRITE(FILE_ID,113)"FINISH_3D_NONLAPLACIAN_STF_MAT Time : ",ITER_TIME_TABLE(11)
    WRITE(FILE_ID,113)"          FINISH_3D_RHS_VECTOR Time : ",ITER_TIME_TABLE(12)
    WRITE(FILE_ID,113)"        CFA_3D_Apply_BCs_Part2 Time : ",ITER_TIME_TABLE(13)
    WRITE(FILE_ID,113)"                    CFA_Solver Time : ",ITER_TIME_TABLE(14)
    WRITE(FILE_ID,113)"        CFA_Coefficient_Update Time : ",ITER_TIME_TABLE(15)
    WRITE(FILE_ID,113)"   CFA_Coefficient_Share_PETSc Time : ",ITER_TIME_TABLE(16)
    WRITE(FILE_ID,113)"         CFA_Convergence_Check Time : ",ITER_TIME_TABLE(17)
    WRITE(FILE_ID,113)"               Total Iteration Time : ",ITER_TIME_TABLE(18)
    WRITE(FILE_ID,113)"             Poseidon_Dist_Sol Time : ",ITER_TIME_TABLE(19)
    WRITE(FILE_ID,110)"============================================================="




END IF



END SUBROUTINE OUTPUT_ITER_TIMETABLE




 !+103+############################################################################!
!                                                                                   !
!                    CLOSE_ITER_REPORT_FILE                                          !
!                                                                                   !
 !#################################################################################!
SUBROUTINE CLOSE_ITER_REPORT_FILE()


CLOSE( UNIT = ITER_REPORT_FILE_ID )


END SUBROUTINE CLOSE_ITER_REPORT_FILE





! !+103+############################################################################!
!!                                                                                   !
!!                    OUTPUT_ITERATION_RESULTS                                          !
!!                                                                                   !
! !#################################################################################!
!SUBROUTINE OUTPUT_ITERATION_RESULTS()
!
!
!
!deltar = ( R_OUTER - R_INNER )/ REAL(ITER_REPORT_NUM_SAMPLES, KIND = idp)
!DO i = 0,ITER_REPORT_NUM_SAMPLES
!
!    r = i*deltar + R_INNER
!    theta = pi/2.0_idp
!    phi = pi/2.0_idp
!
!
!    CALL Calc_3D_Values_At_Location( r, theta, phi,                              &
!                                    Return_Psi, Return_AlphaPsi,                &
!                                    Return_Beta1, Return_Beta2, Return_Beta3    )
!
!    ! Determine the Newtonian Potential at the location r, theta, phi
!    Analytic_Val = Potential_Solution(r,theta,phi)
!
!
!    ! AlphaPsi_to_Pot   =   2*C_Square*(AlphaPsi - 1)
!    ! Psi_to_Pot        =   2*C_Square*(1 - Psi)
!
!    ! Calculate Conformal Factor value from Newtonian Potential
!    PsiPot_Val = 2.0_idp*C_Square*(Return_Psi - 1.0_idp)
!
!    ! Calculate the product of the Conformal Factor and Lapse Function from Newtonian Potential
!    AlphaPsiPot_Val = 2.0_idp*C_Square*(1.0_idp - Return_AlphaPsi)
!
!
!    ! Write Results to Screen
!    IF (( WRITE_REPORT_FLAG == 1) .OR. (WRITE_REPORT_FLAG == 3) ) THEN
!        IF ( Rank == 0 ) THEN
!            WRITE(*,111) r,Analytic_Val,PsiPot_Val,AlphaPsiPot_Val,Return_Beta1,Return_Beta2,Return_Beta3
!        END IF
!    END IF
!
!    ! Write Results to File
!    DO j = 0,1
!        IF ( FILE_ID(j) .NE. -1 ) THEN
!            WRITE(FILE_ID(j),111) r,Analytic_Val,PsiPot_Val,AlphaPsiPot_Val,Return_Beta1,Return_Beta2,Return_Beta3
!        END IF
!    END DO ! j Loop
!
!END DO  ! i Loop
!
!END SUBROUTINE OUTPUT_ITERATION_RESULTS()







 !+401+############################################################################!
!                                                                                   !
!                           OUTPUT_FINAL_RESULTS                                    !
!                                                                                   !
 !#################################################################################!
SUBROUTINE OUTPUT_FINAL_RESULTS()


INTEGER                                                     ::  NUM_SAMPLES


CHARACTER(LEN = 100), DIMENSION(:), ALLOCATABLE             ::  Filenames
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

INTEGER                                                     ::  i, j, k, ORD, Here

INTEGER                                                     ::  NUM_THETA_RAYS,     &
                                                                NUM_PHI_RAYS,       &
                                                                NUM_RADIAL_SAMPLES

REAL(KIND = idp)                                            ::  deltar_overtwo

REAL(KIND = idp)                                            ::  DELTA_THETA,        &
                                                                THETA_VAL

REAL(KIND = idp)                                            ::  DELTA_PHI,        &
                                                                PHI_VAL

REAL(KIND = idp), DIMENSION(:,:,:), ALLOCATABLE             ::  Lapse_Holder,       &
                                                                ConForm_Holder
REAL(KIND = idp), DIMENSION(:,:,:,:), ALLOCATABLE           ::  Shift_Holder
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  R_Holder,           &
                                                                T_Holder,           &
                                                                P_Holder,           &
                                                                xlocs, cur_r_locs


110 FORMAT (11X,A1,16X,A18,9X,A13,10X,A18,10X,A10)                              !!! Output Header

111 FORMAT (11X,A1,24X,A3,19X,A8,15X,A11,14X,A11,14X,A11)                       !!! Output Header for Results file
112 FORMAT (11X,A1,16X,A18,9x,A14)                                              !!! Output Header for Analytic Solution file

113 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15)               !!! Output
114 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15)    !!! Output for Results file
115 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15)                                     !!! Output for Analytic Solution file

116 FORMAT (A,A,I5.5,A4)
117 FORMAT (A,A)


IF ( WRITE_RESULTS_FLAG == 1 ) THEN
    NUM_SAMPLES = 1000



    IF ( DRIVER_TEST_NUMBER == 1 ) THEN

        ALLOCATE( Output_re(0:NUM_SAMPLES) )
        ALLOCATE( Output_rc(1:NUM_SAMPLES) )
        ALLOCATE( Output_dr(1:NUM_SAMPLES) )
        ! Open Results File
        WRITE(filenamea,117) Poseidon_Results_Dir,"Results.out"
        CALL OPEN_NEW_FILE( filenamea, file_ida, 200)
        WRITE(file_ida,111)"r (cm)","Psi","AlphaPsi","Beta1 Value","Beta2 Value","Beta3 Value"

        ! Open Solution File
        WRITE(filenameb,117) Poseidon_Results_Dir,"Solution.out"
        CALL OPEN_NEW_FILE( filenameb, file_idb,200)
        WRITE(file_idb,112)"r (cm)","Analytic Potential","Beta1 Solution"


        ! Set Output locations
        theta = pi/2.0_idp
        phi = pi/2.0_idp

        CALL Create_Logarithmic_1D_Mesh( R_INNER, R_OUTER, NUM_SAMPLES,        &
                                        output_re, output_rc, output_dr                 )



        DO i = 1,NUM_SAMPLES

            CALL Calc_3D_Values_At_Location( output_rc(i), theta, phi,                   &
                                             Return_Psi, Return_AlphaPsi,                &
                                             Return_Beta1, Return_Beta2, Return_Beta3    )

            Analytic_Val = Potential_Solution(output_rc(i),theta,phi)/GravPot_Units
            Solver_Val = 2.0_idp*C_Square*(1.0_idp - Return_Psi)/GravPot_Units
            Solver_Valb = 2.0_idp*C_Square*(Return_AlphaPsi - 1.0_idp)/GravPot_Units
            Error_Val = ABS((Analytic_Val - Solver_Val)/Analytic_Val)


           WRITE(42,114) output_rc(i)*Centimeter, Return_Psi, Return_AlphaPsi, Return_Beta1/Shift_Units,Return_Beta2,Return_Beta3
           WRITE(43,115) output_rc(i)*Centimeter, Analytic_Val, Shift_Solution(output_rc(i),rlocs,NUM_R_ELEMENTS)/Shift_Units



        END DO


        ! Close Files
        CLOSE( Unit = file_ida)
        CLOSE( Unit = file_idb)


    ELSE IF ( DRIVER_TEST_NUMBER == 5 ) THEN

        ALLOCATE( Output_re(0:NUM_SAMPLES) )
        ALLOCATE( Output_rc(1:NUM_SAMPLES) )
        ALLOCATE( Output_dr(1:NUM_SAMPLES) )

        ! Open Results File
        WRITE(filenamea,'(A,A,I5.5,A4)')Poseidon_Results_Dir,"Results_",Poseidon_Frame,".out"
        CALL OPEN_NEW_FILE( filenamea, file_ida, 200)
        WRITE(file_ida,111)"r","Psi","AlphaPsi","Beta1 Value","Beta2 Value","Beta3 Value"

        ! Open Solution File
        WRITE(filenameb,'(A,A,I5.5,A4)')Poseidon_Results_Dir,"Solution_",Poseidon_Frame,".out"
        CALL OPEN_NEW_FILE( filenameb, file_idb, 200)
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

            Analytic_Val = Potential_Solution(output_rc(i),theta,phi)/GravPot_Units
            Solver_Val = 2.0_idp*C_Square*(1.0_idp - Return_Psi)/GravPot_Units
            Solver_Valb = 2.0_idp*C_Square*(Return_AlphaPsi - 1.0_idp)/GravPot_Units
            Error_Val = ABS((Analytic_Val - Solver_Val)/Analytic_Val)


           WRITE(42,114) output_rc(i)/Centimeter, Return_Psi, Return_AlphaPsi, Return_Beta1/Shift_Units,Return_Beta2,Return_Beta3
           WRITE(43,115) output_rc(i)/Centimeter,Analytic_Val, Shift_Solution(output_rc(i),rlocs,NUM_R_ELEMENTS)/Shift_Units



        END DO


        ! Close Files
        CLOSE( Unit = file_ida)
        CLOSE( Unit = file_idb)







    ELSE
    

        Num_Files = 9

        ALLOCATE( Filenames(1:Num_Files) )
        ALLOCATE( File_IDs(1:Num_Files) )

        WRITE(Filenames(1),116) Poseidon_Results_Dir,"Results_Lapse_",Poseidon_Frame,".out"
        WRITE(Filenames(2),116) Poseidon_Results_Dir,"Results_ConFactor_",Poseidon_Frame,".out"
        WRITE(Filenames(3),116) Poseidon_Results_Dir,"Results_Beta1_",Poseidon_Frame,".out"
        WRITE(Filenames(4),116) Poseidon_Results_Dir,"Results_Beta2_",Poseidon_Frame,".out"
        WRITE(Filenames(5),116) Poseidon_Results_Dir,"Results_Beta3_",Poseidon_Frame,".out"
        WRITE(Filenames(6),'(A,A)') Poseidon_Results_Dir,"Results_Dimensions.out"
        WRITE(Filenames(7),116) Poseidon_Results_Dir,"Results_Radial_Locs_",Poseidon_Frame,".out"
        WRITE(Filenames(8),116) Poseidon_Results_Dir,"Results_Theta_Locs_",Poseidon_Frame,".out"
        WRITE(Filenames(9),116) Poseidon_Results_Dir,"Results_Phi_Locs_",Poseidon_Frame,".out"


        File_IDs = [(141 + i, i=1,Num_Files)]
        DO i = 1,Num_Files
            CALL OPEN_NEW_FILE( Filenames(i), File_IDs(i), 200 )
        END DO




        IF ( .FALSE. ) THEN
            NUM_RADIAL_SAMPLES = WRITE_RESULTS_R_SAMPS

            ! Create Radial Spacing !
            ALLOCATE( Output_rc(1:NUM_RADIAL_SAMPLES) )
            ALLOCATE( Output_dr(1:NUM_RADIAL_SAMPLES) )
            ALLOCATE( Output_re(0:NUM_RADIAL_SAMPLES) )
            CALL Create_Logarithmic_1D_Mesh( R_INNER, R_OUTER, NUM_RADIAL_SAMPLES,  &
                                             output_re, output_rc, output_dr        )
        ELSE
            ORD = 2
            NUM_RADIAL_SAMPLES = NUM_R_ELEMENTS*ORD
            ALLOCATE( Output_rc(1:NUM_RADIAL_SAMPLES) )
            ALLOCATE(xlocs(1:ORD))
            ALLOCATE(cur_r_locs(1:ORD))
            xlocs =  Initialize_LG_Quadrature_Locations(ORD)
            DO k = 1,NUM_R_ELEMENTS
                deltar_overtwo = (rlocs(k) - rlocs(k-1))/2.0_idp
                cur_r_locs(:) = deltar_overtwo*(xlocs(:) + 1.0_idp) + rlocs(k-1)
                DO j = 1,ORD
                    here = (k-1)*ord + j
                    Output_rc(here) = cur_r_locs(j)
                END DO
            END DO

        END IF


        ! Create Output Spacing
        ! Pull Number of Samples From Parameters !

        NUM_THETA_RAYS = WRITE_RESULTS_T_SAMPS
        NUM_PHI_RAYS = WRITE_RESULTS_P_SAMPS

        !  Create Phi Spacing !
        IF ( NUM_PHI_RAYS == 1 ) THEN
            DELTA_PHI = pi/2.0_idp
!            DELTA_PHI = 0.0_idp
        ELSE
            DELTA_PHI = 2.0_idp*pi/(NUM_PHI_RAYS-1)
        END IF

        ! Create Theta Spacing
        IF ( NUM_THETA_RAYS == 1 ) THEN
            DELTA_THETA = pi
!            DELTA_PHI = 0.0_idp
        ELSE
            DELTA_THETA = pi/(NUM_THETA_RAYS-1)
        END IF




        ! Allocate Data Holders !
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
        WRITE(File_IDs(7),*)R_Holder/Centimeter
        WRITE(File_IDs(8),*)T_Holder
        WRITE(File_IDs(9),*)P_Holder


        ! Write Output Value Files
        DO k = 1,NUM_PHI_RAYS
            DO j = 1,NUM_THETA_RAYS


                WRITE(File_IDs(1),*)Lapse_Holder(k,j,:)
                WRITE(File_IDs(2),*)ConForm_Holder(k,j,:)
                WRITE(File_IDs(3),*)Shift_Holder(1,k,j,:)/Shift_Units
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




!+402+###########################################################################!
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
        WRITE(FILE_ID(i),'(A,ES22.15,ES22.15,ES22.15,I6.5)')"Solve Tolerances :",rtol, abstol, dtol, maxits
        WRITE(FILE_ID(i),'(A)')" "
        WRITE(FILE_ID(i),'(A)')" "
        WRITE(FILE_ID(i),'(A)')" "

    END IF
END DO


END SUBROUTINE OUTPUT_PETSC_REPORT








 !+501+############################################################################!
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
    RUN_TIME_TABLE(Ident) = (RUN_TIME_TABLE(ident)*(Poseidon_Frame-1) + Time)/Poseidon_Frame
ELSE
    ! Add This Iterations Time to Running Average for the Run
    RUN_TIME_TABLE(Ident) = (RUN_TIME_TABLE(ident)*(Total_Run_Iters-1) + Time)/Total_Run_Iters
END IF


END SUBROUTINE CLOCK_IN



 !+501+############################################################################!
!                                                                                   !
!                     OPEN_NEW_FILE                                                 !
!                                                                                   !
 !#################################################################################!
SUBROUTINE OPEN_NEW_FILE(File_Name, File_Number, Suggested_Number)



CHARACTER(LEN = *), INTENT(IN)                          ::  File_Name
INTEGER,            INTENT(INOUT)                       ::  File_Number
INTEGER, OPTIONAL,  INTENT(IN)                          ::  Suggested_Number

INTEGER                                                 ::  Temp_Number
INTEGER                                                 ::  istat = 0
LOGICAL                                                 ::  FLAG, OP, EX
LOGICAL                                                 ::  UNIT_FLAG, NAME_FLAG


UNIT_FLAG = .FALSE.
NAME_FLAG = .FALSE.


!  Assigned an unused number, and assign it to new file
FLAG = .TRUE.
IF ( Present(Suggested_Number) ) THEN
    Temp_Number = Suggested_Number
ELSE
    Temp_Number = 1000
END IF
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

        PRINT*,"WARNING: Could not open file at ", File_Name, istat

    END IF
END IF


END SUBROUTINE OPEN_NEW_FILE






 !+501+############################################################################!
!                                                                                   !
!                     OPEN_EXISTING_FILE                                            !
!                                                                                   !
 !#################################################################################!
SUBROUTINE OPEN_EXISTING_FILE(File_Name, File_Number, istat)



CHARACTER(LEN = *), INTENT(IN)                          ::  File_Name
INTEGER,            INTENT(INOUT)                       ::  File_Number
INTEGER,            INTENT(INOUT)                       ::  istat


INTEGER                                                 ::  Temp_Number

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

    OPEN(unit=File_Number, file=File_Name, status='old', &
         iostat=istat, action='read', position='rewind')
    IF ( istat .NE. 0 ) THEN

        PRINT*,"WARNING: Could not open file at ", File_Name

    END IF
END IF


END SUBROUTINE OPEN_EXISTING_FILE










 !+601+############################################################################!
!                                                                                   !
!                     OUTPUT_POSIEDON_SOURCES                                       !
!                                                                                   !
 !#################################################################################!
SUBROUTINE OUTPUT_POSEIDON_SOURCES_1D( Local_E, Local_S, Local_Si,                         &
                                    Local_RE_Dim, Local_TE_Dim, Local_PE_Dim,           &
                                    Local_RQ_Dim, Local_TQ_Dim, Local_PQ_Dim,           &
                                    Input_R_Quad, Input_T_Quad, Input_P_Quad,           &
                                    Left_Limit, Right_Limit                             )


REAL(KIND = idp), INTENT(IN), DIMENSION(    1:Local_RQ_Dim*Local_TQ_Dim*Local_PQ_Dim,       &
                                            0:Local_RE_Dim-1,                               &
                                            0:Local_TE_Dim-1,                               &
                                            0:Local_PE_Dim-1  )             ::  Local_E,    &
                                                                                Local_S

REAL(KIND = idp), INTENT(IN), DIMENSION(    1:Local_RQ_Dim*Local_TQ_Dim*Local_PQ_Dim,       &
                                            0:Local_RE_Dim-1,                               &
                                            0:Local_TE_Dim-1,                               &
                                            0:Local_PE_Dim-1,                               &
                                            1:DOMAIN_DIM                )   :: Local_Si





INTEGER, INTENT(IN)                                                     ::  Local_RE_Dim,   &
                                                                            Local_TE_Dim,   &
                                                                            Local_PE_Dim,   &
                                                                            Local_RQ_Dim,   &
                                                                            Local_TQ_Dim,   &
                                                                            Local_PQ_Dim


REAL(KIND = idp), DIMENSION(1:Local_RQ_Dim), INTENT(IN)                 ::  Input_R_Quad
REAL(KIND = idp), DIMENSION(1:Local_TQ_Dim), INTENT(IN)                 ::  Input_T_Quad
REAL(KIND = idp), DIMENSION(1:Local_PQ_Dim), INTENT(IN)                 ::  Input_P_Quad

REAL(KIND = idp), DIMENSION(1:Local_RQ_Dim)                             ::  CUR_R_LOCS

REAL(KIND = idp), INTENT(IN)                                            ::  Left_Limit,     &
                                                                            Right_Limit

INTEGER                                                                 ::  i, re, te, pe, rd, td, pd

CHARACTER(LEN = 100), DIMENSION(:), ALLOCATABLE                         ::  Filenames
INTEGER, DIMENSION(:), ALLOCATABLE                                      ::  File_IDs
INTEGER                                                                 ::  Num_Files

REAL(KIND = idp)                                                        ::  Delta_X, Dr_Over_Dx

REAL(KIND = idp)                                                        ::  E_units, S_units, Si_units

E_Units = Erg/Centimeter**3
S_Units = Erg/Centimeter**3
Si_Units = Gram/(Second*Centimeter**2)



116 FORMAT (A,A,I5.5,A)



IF ( WRITE_SOURCES_FLAG == 1 ) THEN
    Num_Files = 6

    ALLOCATE( Filenames(1:Num_Files) )
    ALLOCATE( File_IDs(1:Num_Files) )

    WRITE(Filenames(1),116) Poseidon_Sources_Dir,"Sources_E_",Poseidon_Frame,".out"
    WRITE(Filenames(2),116) Poseidon_Sources_Dir,"Sources_S_",Poseidon_Frame,".out"
    WRITE(Filenames(3),116) Poseidon_Sources_Dir,"Sources_S1_",Poseidon_Frame,".out"
    WRITE(Filenames(4),'(A,A)') Poseidon_Sources_Dir,"Sources_Dimensions.out"
    WRITE(Filenames(5),116) Poseidon_Sources_Dir,"Sources_Radial_Locs_",Poseidon_Frame,".out"
    
    WRITE(Filenames(6),116) Poseidon_Sources_Dir,"Sources_Time_",Poseidon_Frame,".out"

    !WRITE(Filenames(6),116) Poseidon_Sources_Dir,"Sources_Theta_Locs_",Poseidon_Frame,".out"
    !WRITE(Filenames(7),116) Poseidon_Sources_Dir,"Sources_Phi_Locs_",Poseidon_Frame,".out"



    File_IDs = [(161 + i, i=1,Num_Files)]
    DO i = 1,Num_Files
        CALL OPEN_NEW_FILE( Filenames(i), File_IDs(i) )
    END DO

    

    WRITE(File_IDs(4),* )Local_RE_Dim, Local_TE_Dim, Local_PE_DIM
    WRITE(File_IDs(4),* )Local_RQ_Dim, Local_TQ_Dim, Local_PQ_DIM
    WRITE(File_IDs(6),* )SELFSIM_T


    Delta_X = Right_Limit - Left_Limit

    DO re = 0,NUM_R_ELEMENTS-1

        Dr_Over_Dx = (rlocs(re+1) - rlocs(re))/Delta_X
        CUR_R_LOCS(:) = Dr_Over_Dx * (INPUT_R_QUAD(:)-Left_Limit) + rlocs(re)


        DO rd = 1,Local_RQ_Dim
!            PRINT*,Local_E(rd,re,0,0),Local_S(rd,re,0,0)
            WRITE(File_IDs(5),*) CUR_R_LOCS(rd)/Centimeter
            WRITE(File_IDs(1),*) Local_E(rd,re,0,0)/E_Units
            WRITE(File_IDs(2),*) Local_S(rd,re,0,0)/S_Units
            WRITE(File_IDs(3),*) Local_Si(rd,re,0,0,1)/Si_Units
        END DO
    END DO







    ! Close Files
    DO i = 1,Num_Files
        CLOSE( Unit = File_IDs(i))
    END DO

END IF

END SUBROUTINE OUTPUT_POSEIDON_SOURCES_1D















!+401+###########################################################################!
!                                                                                !
!              Calc_Shift_1D                                                  !
!                                                                                !
!################################################################################!
SUBROUTINE Write_Shift_1D( Shift_Vector, NUM_R_ELEM, R_LOCS, SelfSim_T )

REAL(KIND = idp), DIMENSION( 0:NUM_R_ELEM ),INTENT( IN )          ::  Shift_Vector

INTEGER,               INTENT( IN )                               ::  NUM_R_ELEM

REAL(KIND = idp), DIMENSION( 0:NUM_R_ELEM ),  INTENT(IN)          ::  r_locs
REAL(KIND = idp)                                                  ::  SelfSim_T


INTEGER                                                           ::  re

INTEGER                                                           ::  File_id

CHARACTER( LEN = 43 )                                             ::  FileName

CHARACTER( LEN = 20 )                                             ::  FileDir
CHARACTER( LEN = 13 )                                             ::  FilePre
CHARACTER( LEN = 4 )                                              ::  FileExt


FileDir = "Shift_Vector_Output/"
FilePre = "Shift_Vector_"
FileExt = ".out"

WRITE( FileName, '(A,A,F6.4,A)' ) FileDir,FilePre,SelfSim_T,FileExt


PRINT*,Filename

File_id = 42

OPEN(Unit = File_id, file = FileName)

DO re = 0,NUM_R_ELEM

   WRITE(File_id,'(ES22.15,1X,ES22.15)') r_locs(re), Shift_Vector(re)

END DO


CLOSE(Unit = File_id)

END SUBROUTINE Write_Shift_1D








SUBROUTINE OPEN_FILE_INQUISITION()


INTEGER                 :: i
INTEGER                 :: i_min
INTEGER                 :: i_max
LOGICAL                 :: OP
CHARACTER(len=150)        :: name_of_file

i_min = 1
i_max = 1000

DO i = i_min,i_max
    OP = .FALSE.
    INQUIRE( UNIT = i, OPENED = OP, name=name_of_file )
    IF ( OP .EQV. .TRUE. ) THEN
        PRINT*,"File ",TRIM(name_of_file)," with Unit ",i," is still open."
    END IF
END DO


END SUBROUTINE OPEN_FILE_INQUISITION





END MODULE Poseidon_IO_Module
