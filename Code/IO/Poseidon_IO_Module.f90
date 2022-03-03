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


USE Poseidon_Kinds_Module, &
            ONLY : idp

USE Poseidon_Numbers_Module, &
            ONLY : pi

USE Poseidon_Units_Module, &
            ONLY :  C_Square,       &
                    Gram,           &
                    Centimeter,     &
                    Kilometer,      &
                    Erg,            &
                    Second,         &
                    GravPot_Units,  &
                    E_Units,        &
                    S_Units,        &
                    Si_Units,       &
                    Shift_Units

USE Parameters_Variable_Indices, &
            ONLY :  iU_CF,      &
                    iU_LF,      &
                    iU_S1,      &
                    iU_S2,      &
                    iU_S3,      &
                    iU_X1,      &
                    iU_X2,      &
                    iU_X3,      &
                    iVB_S,      &
                    iVB_X

USE Poseidon_Parameters, &
            ONLY :  DEGREE,                 &
                    L_LIMIT,                &
                    Num_CFA_Vars,           &
                    Domain_Dim,             &
                    CUR_ITERATION,          &
                    Poseidon_Frame,         &
                    Convergence_Type,       &
                    Convergence_Flag,       &
                    Max_Iterations

USE Variables_MPI, &
            ONLY :  myID_Poseidon

USE Variables_Derived, &
            ONLY :  Num_R_Nodes


USE Variables_Mesh, &
            ONLY :  NUM_R_ELEMENTS,         &
                    rlocs,                  &
                    tlocs,                  &
                    R_Inner,                &
                    R_Outer

USE Variables_IO, &
            ONLY :  Write_Flags,            &
                    Report_Flags,           &
                    Report_IDs,             &
                    Iter_Time_Table,        &
                    Frame_Time_Table,       &
                    Run_Time_Table,         &
                    Write_Results_R_Samps,  &
                    Write_Results_T_Samps,  &
                    Write_Results_P_Samps,  &
                    Total_Run_Iters,        &
                    Iter_Report_Num_Samples,&
                    Iter_Time_Table,        &
                    Frame_Residual_Table,   &
                    Frame_Update_Table,     &
                    Iteration_Histogram,    &
                    File_Suffix,            &
                    iWF_Source,             &
                    iWF_Results,            &
                    iRF_Run,                &
                    iRF_Frame,              &
                    iRF_Time,               &
                    iRF_Iter


USE Variables_Yahil, &
            ONLY :  SelfSim_T

USE Variables_Functions, &
            ONLY :  Potential_Solution,         &
                    Shift_Solution,             &
                    Calc_3D_Values_at_Location, &
                    Calc_1D_CFA_Values,         &
                    Calc_Var_At_Loc_A,          &
                    Calc_Var_At_Loc_B

USE Functions_Quadrature, &
            ONLY :  Initialize_LG_Quadrature_Locations,     &
                    Initialize_LGL_Quadrature_Locations


USE Functions_Mesh, &
            ONLY :  Create_Logarithmic_1D_Mesh,             &
                    Create_Uniform_1D_Mesh

USE Poseidon_IO_Parameters, &
            ONLY :  Poseidon_Reports_Dir,                           &
                    Poseidon_IterReports_Dir,                       &
                    Poseidon_Objects_Dir,                           &
                    Poseidon_Mesh_Dir,                              &
                    Poseidon_LinSys_Dir,                            &
                    Poseidon_Results_Dir,                           &
                    Poseidon_Sources_Dir,                           &
                    CFA_ShortVars

USE Functions_Info,   &
            ONLY  : PQ_ITERATIONS_MAX

USE IO_File_Routines_Module, &
            ONLY :  Open_New_File,                  &
                    Open_Existing_File

USE Variables_FP, &
            ONLY :  CFA_Eq_Flags



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

109 FORMAT (A,I2.2,A,I2.2)
113 FORMAT (A,A,I2.2,A,I2.2,A)

IF (( Report_Flags(iRF_Iter) == 2) .OR. (Report_Flags(iRF_Iter) == 3) ) THEN

    WRITE(FILE_NAME,113) Poseidon_IterReports_Dir,"Iteration_Report_P",Rank,"_I",Iteration,".out"
    CALL OPEN_NEW_FILE( FILE_NAME, Report_IDs(iRF_Iter), 100 )



    WRITE(Report_IDs(iRF_Iter),109)"                      Report for Iteration ",Iteration," by Process ",Rank
    WRITE(Report_IDs(iRF_Iter),'(A)')"-----------------------------------------------------------------------------------------"
    WRITE(Report_IDs(iRF_Iter),'(A)')" "

END IF


IF ( Report_Flags(iRF_Frame) == 1 ) THEN

    WRITE(Report_IDs(iRF_Frame),'(A)')"-----------------------------------------------------------------------------------------"
    WRITE(Report_IDs(iRF_Frame),109)"                      Report for Iteration ",Iteration," by Process ",Rank
    WRITE(Report_IDs(iRF_Frame),'(A)')"-----------------------------------------------------------------------------------------"
    WRITE(Report_IDs(iRF_Frame),'(3/)')

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


IF (( Report_Flags(iRF_Time) == 1 ) .OR. ( Report_Flags(iRF_Time) == 3 ) ) THEN

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


IF (( Report_Flags(iRF_Time) == 2 ) .OR. ( Report_Flags(iRF_Time) == 3 ) ) THEN

     
    WRITE(FILENAME,'(A,I2.2,A)')'OUTPUT/Timetable_',Ident,'.out'

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


CLOSE( UNIT = Report_IDs(iRF_Iter) )

END SUBROUTINE CLOSE_ITER_REPORT_FILE













 !+401+############################################################################!
!                                                                                   !
!                           OUTPUT_FINAL_RESULTS                                    !
!                                                                                   !
 !#################################################################################!
SUBROUTINE OUTPUT_FINAL_RESULTS()


CHARACTER(LEN = 100), DIMENSION(:), ALLOCATABLE             ::  Filenames
INTEGER, DIMENSION(:), ALLOCATABLE                          ::  File_IDs
INTEGER                                                     ::  Num_Files


REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Output_re,          &
                                                                Output_rc,          &
                                                                Output_dr

REAL(KIND = idp)                                            ::  Return_Psi,         &
                                                                Return_AlphaPsi,    &
                                                                Return_Beta1,       &
                                                                Return_Beta2,       &
                                                                Return_Beta3

INTEGER                                                     ::  i, j, k, Here

INTEGER                                                     ::  NUM_THETA_RAYS,     &
                                                                NUM_PHI_RAYS,       &
                                                                NUM_RADIAL_SAMPLES


REAL(KIND = idp)                                            ::  DELTA_THETA,        &
                                                                THETA_VAL

REAL(KIND = idp)                                            ::  DELTA_PHI,        &
                                                                PHI_VAL


REAL(KIND = idp), DIMENSION(:,:,:,:), ALLOCATABLE           ::  Var_Holder
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  R_Holder,           &
                                                                T_Holder,           &
                                                                P_Holder

REAL(idp), DIMENSION(6)                          ::  Units
REAL(idp)                                                   ::  CF

116 FORMAT (A,A,A,A,A,A)



IF ( Write_Flags(iWF_Results) > 1 ) THEN


    Num_Files = 4 + SUM(CFA_Eq_Flags)+1

    ALLOCATE( Filenames(1:Num_Files) )
    ALLOCATE( File_IDs(1:Num_Files) )




    WRITE(Filenames(1),116) Poseidon_Results_Dir,"Results_Dimensions_",TRIM(File_Suffix),".out"
    WRITE(Filenames(2),116) Poseidon_Results_Dir,"Results_Radial_Locs_",TRIM(File_Suffix),".out"
    WRITE(Filenames(3),116) Poseidon_Results_Dir,"Results_Theta_Locs_",TRIM(File_Suffix),".out"
    WRITE(Filenames(4),116) Poseidon_Results_Dir,"Results_Phi_Locs_",TRIM(File_Suffix),".out"
    Here = 5
    DO i = 1,5
        IF ( CFA_Eq_Flags(i) == 1 ) THEN
            Here = 5 + i - 1
            WRITE(Filenames(Here),116)                &
                    Poseidon_Results_Dir,"Results_",TRIM(CFA_ShortVars(i)),"_",TRIM(File_Suffix),".out"
        END IF
    END DO
    WRITE(Filenames(Num_Files),116)                &
    Poseidon_Results_Dir,"Results_",TRIM(CFA_ShortVars(6)),"_",TRIM(File_Suffix),".out"






!        DO i = 0,PROB_DIM-1
!            PRINT*,Coefficient_Vector(i)
!        END DO

    DO i = 1,Num_Files
        CALL OPEN_NEW_FILE( Filenames(i), File_IDs(i), 200 )
    END DO


    NUM_RADIAL_SAMPLES = WRITE_RESULTS_R_SAMPS
    NUM_THETA_RAYS = WRITE_RESULTS_T_SAMPS
    NUM_PHI_RAYS = WRITE_RESULTS_P_SAMPS
    

    ! Create Radial Spacing !
    ALLOCATE( Output_rc(1:NUM_RADIAL_SAMPLES) )
    ALLOCATE( Output_dr(1:NUM_RADIAL_SAMPLES) )
    ALLOCATE( Output_re(0:NUM_RADIAL_SAMPLES) )

    IF ( R_OUTER/(R_Inner+1.0_idp) > 1E3 ) THEN
        CALL Create_Logarithmic_1D_Mesh( R_INNER, R_OUTER, NUM_RADIAL_SAMPLES,     &
                                         output_re, output_rc, output_dr     )
    ELSE
        CALL Create_Uniform_1D_Mesh( R_INNER, R_OUTER, NUM_RADIAL_SAMPLES,     &
                                     output_re, output_rc, output_dr     )
    END IF

    ! Create Output Spacing
    ! Pull Number of Samples From Parameters !

    !  Create Phi Spacing !
    IF ( NUM_PHI_RAYS == 1 ) THEN
        DELTA_PHI = pi/2.0_idp
!        DELTA_PHI = 0.0_idp
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
    ALLOCATE( Var_Holder(1:NUM_PHI_RAYS, 1:NUM_THETA_RAYS, 1:NUM_RADIAL_SAMPLES,1:6) )
    ALLOCATE( R_Holder(1:NUM_RADIAL_SAMPLES) )
    ALLOCATE( T_Holder(1:NUM_THETA_RAYS) )
    ALLOCATE( P_Holder(1:NUM_PHI_RAYS) )


    ! Calculate Output
    DO k = 1,NUM_PHI_RAYS
    DO j = 1,NUM_THETA_RAYS
    DO i = 1,NUM_RADIAL_SAMPLES



        PHI_VAL = k*DELTA_PHI
        THETA_VAL = (j-1)*DELTA_THETA

!        CALL Calc_3D_Values_At_Location( output_re(i), THETA_VAL, PHI_VAL,           &
!                                         Return_Psi, Return_AlphaPsi,                &
!                                         Return_Beta1, Return_Beta2, Return_Beta3    )
!!
!!
!        Var_Holder(k,j,i,1) = Return_Psi
!        IF ( Return_Psi == 0.0_idp ) THEN
!            Var_Holder(k,j,i,2) = 0.0_idp
!        ELSE
!            Var_Holder(k,j,i,2) = Return_AlphaPsi/Return_Psi
!        END IF
!        Var_Holder(k,j,i,3:5) = (/ Return_Beta1, Return_Beta2, Return_Beta3 /)


        CF = Calc_Var_At_Loc_A(Output_re(i),theta_Val,phi_val, iU_CF)

        Var_Holder(k,j,i,1) = CF
        Var_Holder(k,j,i,2) = Calc_Var_At_Loc_A(Output_re(i),theta_Val,phi_val, iU_LF)/CF
        Var_Holder(k,j,i,3) = Calc_Var_At_Loc_B(Output_re(i),theta_Val,phi_val, iU_S1, iVB_S)
        Var_Holder(k,j,i,4) = Calc_Var_At_Loc_B(Output_re(i),theta_Val,phi_val, iU_S2, iVB_S)
        Var_Holder(k,j,i,5) = Calc_Var_At_Loc_B(Output_re(i),theta_Val,phi_val, iU_S3, iVB_S)
        Var_Holder(k,j,i,6) = Calc_Var_At_Loc_B(Output_re(i),theta_Val,phi_val, iU_X1, iVB_X)

        R_Holder(i) = output_re(i)
        T_Holder(j) = THETA_VAL
        P_Holder(k) = PHI_VAL



    END DO ! i Loop
    END DO ! j Loop
    END DO ! k Loop



    ! Write Output Location Files
    WRITE(File_IDs(1),*)Num_Radial_Samples, Num_Theta_Rays, Num_Phi_Rays
    WRITE(File_IDs(2),*)R_Holder/Centimeter
    WRITE(File_IDs(3),*)T_Holder
    WRITE(File_IDs(4),*)P_Holder


    
    
    Units = 1.0_idp
    Units(3) = Shift_Units
    ! Write Output Value Files

    Here = 5
    DO i = 1,5
        IF ( CFA_Eq_Flags(i) == 1 ) THEN
            DO k = 1,NUM_PHI_RAYS
            DO j = 1,NUM_THETA_RAYS
    
                WRITE(File_IDs(Here),*)Var_Holder(k,j,:,i)/Units(i)

            END DO ! j Loop
            END DO ! k Loop
            Here = Here + 1
        END IF
    END DO


    DO k = 1,NUM_PHI_RAYS
    DO j = 1,NUM_THETA_RAYS

        WRITE(File_IDs(Num_Files),*)Var_Holder(k,j,:,6)/Units(3)

    END DO ! j Loop
    END DO ! k Loop


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
IF ( Report_Flags(iRF_Frame) == 1 ) THEN

    FILE_ID(0) = Report_IDs(iRF_Frame)

END IF

IF (( Report_Flags(iRF_Iter) == 2) .OR. (Report_Flags(iRF_Iter) == 3) ) THEN
    FILE_ID(1) = Report_IDs(iRF_Iter)
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

INTEGER                                                                 ::  i, re, rd

CHARACTER(LEN = 100), DIMENSION(:), ALLOCATABLE                         ::  Filenames
INTEGER, DIMENSION(:), ALLOCATABLE                                      ::  File_IDs
INTEGER                                                                 ::  Num_Files

REAL(KIND = idp)                                                        ::  Delta_X, Dr_Over_Dx





116 FORMAT (A,A,A,A)

IF ( Write_Flags(iWF_Source) == 1 ) THEN
    Num_Files = 6

    ALLOCATE( Filenames(1:Num_Files) )
    ALLOCATE( File_IDs(1:Num_Files) )

    WRITE(Filenames(1),116) Poseidon_Sources_Dir,"Sources_E_",trim(File_Suffix),".out"
    WRITE(Filenames(2),116) Poseidon_Sources_Dir,"Sources_S_",trim(File_Suffix),".out"
    WRITE(Filenames(3),116) Poseidon_Sources_Dir,"Sources_S1_",trim(File_Suffix),".out"
    WRITE(Filenames(4),116) Poseidon_Sources_Dir,"Sources_Dimensions_",trim(File_Suffix),".out"
    WRITE(Filenames(5),116) Poseidon_Sources_Dir,"Sources_Radial_Locs_",trim(File_Suffix),".out"
    
    WRITE(Filenames(6),116) Poseidon_Sources_Dir,"Sources_Time_",trim(File_Suffix),".out"

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





 !+601+############################################################################!
!                                                                                   !
!                     OUTPUT_POSIEDON_SOURCES                                       !
!                                                                                   !
 !#################################################################################!
SUBROUTINE OUTPUT_POSEIDON_SOURCES_3D( Local_E, Local_S, Local_Si,                         &
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
REAL(KIND = idp), DIMENSION(1:Local_TQ_Dim)                             ::  CUR_T_LOCS


REAL(KIND = idp), INTENT(IN)                                            ::  Left_Limit,     &
                                                                            Right_Limit

INTEGER                                                                 ::  i, re, te, pe,  &
                                                                            rd, td, pd, tpd

CHARACTER(LEN = 100), DIMENSION(:), ALLOCATABLE                         ::  Filenames
INTEGER, DIMENSION(:), ALLOCATABLE                                      ::  File_IDs
INTEGER                                                                 ::  Num_Files

REAL(KIND = idp)                                                        ::  Delta_X, Dr_Over_Dx
REAL(KIND = idp)                                                        ::  Dt_Over_DX


INTEGER                                                                 ::  iSF_E     = 1
INTEGER                                                                 ::  iSF_S     = 2
INTEGER                                                                 ::  iSF_S1    = 3
INTEGER                                                                 ::  iSF_X1    = 4
INTEGER                                                                 ::  iSF_Dim   = 5
INTEGER                                                                 ::  iSF_RLocs = 6
INTEGER                                                                 ::  iSF_TLocs = 7
INTEGER                                                                 ::  iSF_PLocs = 8
INTEGER                                                                 ::  iSF_Time  = 9



116 FORMAT (A,A,A,A)


IF ( Write_Flags(iWF_Source) == 2 ) THEN
    Num_Files = 9

    ALLOCATE( Filenames(1:Num_Files) )
    ALLOCATE( File_IDs(1:Num_Files) )




    WRITE(Filenames(iSF_E),116) Poseidon_Sources_Dir,"Sources_E_",trim(File_Suffix),".out"
    WRITE(Filenames(iSF_S),116) Poseidon_Sources_Dir,"Sources_S_",trim(File_Suffix),".out"
    WRITE(Filenames(iSF_S1),116) Poseidon_Sources_Dir,"Sources_S1_",trim(File_Suffix),".out"
    WRITE(Filenames(iSF_X1),116) Poseidon_Sources_Dir,"Sources_S1_",trim(File_Suffix),".out"
    WRITE(Filenames(iSF_Dim),116) Poseidon_Sources_Dir,"Sources_Dimensions_",trim(File_Suffix),".out"
    WRITE(Filenames(iSF_Rlocs),116) Poseidon_Sources_Dir,"Sources_Radial_Locs_",trim(File_Suffix),".out"
    WRITE(Filenames(iSF_Tlocs),116) Poseidon_Sources_Dir,"Sources_Theta_Locs_",trim(File_Suffix),".out"
    WRITE(Filenames(iSF_Plocs),116) Poseidon_Sources_Dir,"Sources_Phi_Locs_",trim(File_Suffix),".out"

    WRITE(Filenames(iSF_Time),116) Poseidon_Sources_Dir,"Sources_Time_",trim(File_Suffix),".out"

    !WRITE(Filenames(6),116) Poseidon_Sources_Dir,"Sources_Theta_Locs_",Poseidon_Frame,".out"
    !WRITE(Filenames(7),116) Poseidon_Sources_Dir,"Sources_Phi_Locs_",Poseidon_Frame,".out"



    File_IDs = [(161 + i, i=1,Num_Files)]
    DO i = 1,Num_Files
        CALL OPEN_NEW_FILE( Filenames(i), File_IDs(i) )
    END DO

    

    WRITE(File_IDs(4),* )Local_RE_Dim, Local_TE_Dim, Local_PE_DIM
    WRITE(File_IDs(4),* )Local_RQ_Dim, Local_TQ_Dim, Local_PQ_DIM
    WRITE(File_IDs(8),* )SELFSIM_T


    Delta_X = Right_Limit - Left_Limit
    
!    DO re = 0,Local_RE_Dim-1
!    DO te = 0,Local_TE_Dim-1
!    DO pe = 0,Local_PE_Dim-1
!
!        DO rd = 1,Local_RQ_Dim
!        DO td = 1,Local_TQ_Dim
!        DO pd = 1,Local_PQ_Dim
!
!            tpd = (rd-1)*Local_TQ_Dim*Local_PQ_Dim      &
!                + (td-1)*Local_PQ_Dim                   &
!                + pd
!
!
!            WRITE(File_IDs(1),*) Local_E(tpd,re,te,pe)/E_Units
!            WRITE(File_IDs(2),*) Local_S(tpd,re,te,pe)/S_Units
!            WRITE(File_IDs(3),*) Local_Si(tpd,re,te,pe,1)/Si_Units
!
!        END DO  ! pd
!        END DO  ! td
!        END DO  ! rd
!
!
!    END DO  ! PE
!    END DO  ! TE
!    END DO  ! RE

    DO pe = 0,Local_PE_Dim-1
    DO pd = 1,Local_PQ_Dim
    DO te = 0,Local_TE_Dim-1
    DO td = 1,Local_TQ_Dim
    DO re = 0,Local_RE_Dim-1

        DO rd = 1,Local_RQ_Dim

            tpd = (rd-1)*Local_TQ_Dim*Local_PQ_Dim      &
                + (td-1)*Local_PQ_Dim                   &
                + pd

            
            WRITE(File_IDs(1),*) Local_E(tpd,re,te,pe)/E_Units
            WRITE(File_IDs(2),*) Local_S(tpd,re,te,pe)/S_Units
            WRITE(File_IDs(3),*) Local_Si(tpd,re,te,pe,1)/Si_Units

        END DO  ! pd
        END DO  ! td
        END DO  ! rd


    END DO  ! PE
    END DO  ! TE
    END DO  ! RE

        


    DO re = 0,Local_RE_Dim-1
        Dr_Over_Dx = (rlocs(re+1) - rlocs(re))/Delta_X
        CUR_R_LOCS(:) = Dr_Over_Dx * (INPUT_R_QUAD(:)-Left_Limit) + rlocs(re)
        DO rd = 1,Local_RQ_Dim
            WRITE(File_IDs(5),*) CUR_R_LOCS(rd)/Centimeter
        END DO
    END DO

    DO te = 0,Local_TE_Dim-1
        Dt_Over_Dx = (tlocs(te+1) - tlocs(te))/Delta_X
        CUR_T_LOCS(:) = Dt_Over_Dx * (Input_T_Quad(:)-Left_Limit) + tlocs(te)
        DO td = 1,Local_TQ_Dim
            WRITE(File_IDs(6),*) CUR_T_LOCS(td)
        END DO
    END DO

    DO pe = 0,Local_PE_Dim-1
    DO pd = 1,Local_PQ_Dim
        WRITE(File_IDs(7),*) Input_P_Quad(pd)
    END DO
    END DO
    



    ! Close Files
    DO i = 1,Num_Files
        CLOSE( Unit = File_IDs(i))
    END DO

END IF




END SUBROUTINE OUTPUT_POSEIDON_SOURCES_3D









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







 !+201+############################################################################!
!                                                                                   !
!                     OPEN_FILE_INQUISITION                                         !
!                                                                                   !
 !#################################################################################!
SUBROUTINE OPEN_FILE_INQUISITION()


INTEGER                 :: i
INTEGER                 :: i_min
INTEGER                 :: i_max
LOGICAL                 :: OP
CHARACTER(len=150)        :: name_of_file

i_min = 1
i_max = 3000

DO i = i_min,i_max
    OP = .FALSE.
    INQUIRE( UNIT = i, OPENED = OP, name=name_of_file )
    IF ( OP .EQV. .TRUE. ) THEN
        PRINT*,"File ",TRIM(name_of_file)," with Unit ",i," is still open."
    END IF
END DO


END SUBROUTINE OPEN_FILE_INQUISITION



 !+201+############################################################################!
!                                                                                   !
!                     OPEN_RUN_REPORT_FILE                                          !
!                                                                                   !
 !#################################################################################!
SUBROUTINE OPEN_RUN_REPORT_FILE()

CHARACTER(LEN = 39)                                     ::  FILE_NAME

IF ( Report_Flags(iRF_Run) == 1 ) THEN
    WRITE(FILE_NAME,'(A,A)') Poseidon_Reports_Dir,"Run_Report.out"
    CALL OPEN_NEW_FILE( FILE_NAME, Report_IDs(iRF_Run), 10 )
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
REAL(KIND = idp)                                ::  r, theta, phi
REAL(KIND = idp), DIMENSION(0:ITER_REPORT_NUM_SAMPLES)  ::  x_e
REAL(KIND = idp), DIMENSION(1:ITER_REPORT_NUM_SAMPLES)  ::  x_c, dx_c

REAL(KIND = idp)                                ::  Return_Psi, Return_AlphaPsi
REAL(KIND = idp)                                ::  Return_Beta1, Return_Beta2, Return_Beta3
REAL(KIND = idp)                                ::  PsiPot_Val, AlphaPsiPot_Val


INTEGER                                         ::  Max_Iters


120 FORMAT (A61)
121 FORMAT (A1)
123 FORMAT (A38,ES22.15)

110 FORMAT (11X,A1,18X,A13,10X,A18,10X,A11,14X,A11,14X,A11)
111 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15)

FILE_ID = Report_IDs(iRF_Run)



CALL PQ_ITERATIONS_MAX( Max_Iters )


! Write Timetable to File
IF ( Report_Flags(iRF_Run) == 1 ) THEN

    WRITE(FILE_ID,'(A)')"                                Average Timing Results"
    WRITE(FILE_ID,'(A)')"            ============================================================="
    WRITE(FILE_ID,'(A)')" "
    WRITE(FILE_ID,123)"                    Initialize Time : ",Run_Time_Table(1)
    WRITE(FILE_ID,123)" Input/Communicate Source Data Time : ",Run_Time_Table(2)
    WRITE(FILE_ID,123)"     Input Boundary Conditions Time : ",Run_Time_Table(3)
    WRITE(FILE_ID,123)"        CFA_3D_Apply_BCs_Part1 Time : ",Run_Time_Table(4)
    WRITE(FILE_ID,120)"-------------------------------------------------------------"
    WRITE(FILE_ID,123)" ||     Calc_3D_Current_Values Time : ",Run_Time_Table(5)
    WRITE(FILE_ID,123)" ||    CREATE_3D_SubJcbn_Terms Time : ",Run_Time_Table(6)
    WRITE(FILE_ID,123)" ||       CREATE_3D_RHS_VECTOR Time : ",Run_Time_Table(7)
    WRITE(FILE_ID,123)"\  /     CREATE_3D_JCBN_MATRIX Time : ",Run_Time_Table(8)
    WRITE(FILE_ID,120)"-\/ ---------------------------------------------------------"
    WRITE(FILE_ID,123)"CREATE_3D_NONLAPLACIAN_STF_MAT Time : ",Run_Time_Table(9)
    WRITE(FILE_ID,123)"REDUCE_3D_NONLAPLACIAN_STF_MAT Time : ",Run_Time_Table(10)
    WRITE(FILE_ID,123)"FINISH_3D_NONLAPLACIAN_STF_MAT Time : ",Run_Time_Table(11)
    WRITE(FILE_ID,123)"          FINISH_3D_RHS_VECTOR Time : ",Run_Time_Table(12)
    WRITE(FILE_ID,123)"        CFA_3D_Apply_BCs_Part2 Time : ",Run_Time_Table(13)
    WRITE(FILE_ID,123)"                    CFA_Solver Time : ",Run_Time_Table(14)
    WRITE(FILE_ID,123)"        CFA_Coefficient_Update Time : ",Run_Time_Table(15)
    WRITE(FILE_ID,123)"   CFA_Coefficient_Share_PETSc Time : ",Run_Time_Table(16)
    WRITE(FILE_ID,123)"         CFA_Convergence_Check Time : ",Run_Time_Table(17)
    WRITE(FILE_ID,123)"               Total Iteration Time : ",Run_Time_Table(18)
    WRITE(FILE_ID,123)"             Poseidon_Dist_Sol Time : ",Run_Time_Table(19)
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
        WRITE(FILE_ID,'(A,I2.2,A,I2.2)')"    ",i,"   |         ",Iteration_Histogram(i)
    END DO
    WRITE(FILE_ID,'(A)')" "
    WRITE(FILE_ID,'(A)')" "


END IF



! Write Timetable to Screen
IF (( Report_Flags(iRF_Iter) == 1) .OR. (Report_Flags(iRF_Iter) == 3) ) THEN
    WRITE(*,'(A)')" "
    WRITE(*,'(A)')" "
    WRITE(*,'(A)')"            ============================================================="
    WRITE(*,'(A)')"                                Average Timing Results"
    WRITE(*,'(A)')"            ============================================================="
    WRITE(*,'(A)')" "
    WRITE(*,123)"                    Initialize Time : ",Run_Time_Table(1)
    WRITE(*,123)" Input/Communicate Source Data Time : ",Run_Time_Table(2)
    WRITE(*,123)"     Input Boundary Conditions Time : ",Run_Time_Table(3)
    WRITE(*,123)"        CFA_3D_Apply_BCs_Part1 Time : ",Run_Time_Table(4)
    WRITE(*,120)"-------------------------------------------------------------"
    WRITE(*,123)" ||     Calc_3D_Current_Values Time : ",Run_Time_Table(5)
    WRITE(*,123)" ||    CREATE_3D_SubJcbn_Terms Time : ",Run_Time_Table(6)
    WRITE(*,123)" ||       CREATE_3D_RHS_VECTOR Time : ",Run_Time_Table(7)
    WRITE(*,123)"\  /     CREATE_3D_JCBN_MATRIX Time : ",Run_Time_Table(8)
    WRITE(*,120)"-\/ ---------------------------------------------------------"
    WRITE(*,123)"CREATE_3D_NONLAPLACIAN_STF_MAT Time : ",Run_Time_Table(9)
    WRITE(*,123)"REDUCE_3D_NONLAPLACIAN_STF_MAT Time : ",Run_Time_Table(10)
    WRITE(*,123)"FINISH_3D_NONLAPLACIAN_STF_MAT Time : ",Run_Time_Table(11)
    WRITE(*,123)"          FINISH_3D_RHS_VECTOR Time : ",Run_Time_Table(12)
    WRITE(*,123)"        CFA_3D_Apply_BCs_Part2 Time : ",Run_Time_Table(13)
    WRITE(*,123)"                    CFA_Solver Time : ",Run_Time_Table(14)
    WRITE(*,123)"        CFA_Coefficient_Update Time : ",Run_Time_Table(15)
    WRITE(*,123)"   CFA_Coefficient_Share_PETSc Time : ",Run_Time_Table(16)
    WRITE(*,123)"         CFA_Convergence_Check Time : ",Run_Time_Table(17)
    WRITE(*,123)"               Total Iteration Time : ",Run_Time_Table(18)
    WRITE(*,123)"             Poseidon_Dist_Sol Time : ",Run_Time_Table(19)
    WRITE(*,120)"============================================================="
    WRITE(*,121)" "
    WRITE(*,121)" "
    WRITE(*,121)" "
    WRITE(*,121)" "

END IF


Report_Flags(iRF_Iter) = 1
PRINT*,"Report_Flags(iRF_Iter)",Report_Flags(iRF_Iter)
! Write Results Table Header to Screen
IF (( Report_Flags(iRF_Iter) == 1) .OR. (Report_Flags(iRF_Iter) == 3) ) THEN
    WRITE(*,'(A)')"++++++++++++++++++++++++++ Sample Run Results ++++++++++++++++++++++++++"
    WRITE(*,110)"r","Psi Potential","AlphaPsi Potential","Beta Value1","Beta Value2","Beta Value3"

    IF ( Report_Flags(iRF_Run) == 1 ) THEN
        WRITE(FILE_ID,'(A)')"++++++++++++++++++++++++++ Sample Run Results ++++++++++++++++++++++++++"
        WRITE(FILE_ID,110)"r","Psi Potential","AlphaPsi Potential","Beta Value1","Beta Value2","Beta Value3"
    END IF

    Call Create_Logarithmic_1D_Mesh( R_Inner, R_Outer, ITER_REPORT_NUM_SAMPLES,     &
                                     x_e, x_c, dx_c )


    !deltar = ( R_Outer - R_Inner )/ REAL(ITER_REPORT_NUM_SAMPLES, KIND = idp)
    DO i = 0,ITER_REPORT_NUM_SAMPLES

    !    r = i*deltar + R_Inner
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
        WRITE(*,111) r/Centimeter,              &
                     PsiPot_Val,                &
                     AlphaPsiPot_Val,           &
                     Return_Beta1/Shift_Units,  &
                     Return_Beta1,              &
                     Return_Beta3

        ! Write Results to File
        IF ( Report_Flags(iRF_Run) == 1 ) THEN
            WRITE(FILE_ID,111) r/Centimeter,              &
                               PsiPot_Val,                &
                               AlphaPsiPot_Val,           &
                               Return_Beta1/Shift_Units,  &
                               Return_Beta2,              &
                               Return_Beta3
        END IF

    END DO

END IF





END SUBROUTINE OUTPUT_RUN_REPORT





 !+203+############################################################################!
!                                                                                   !
!                    CLOSE_RUN_REPORT_FILE                                          !
!                                                                                   !
 !#################################################################################!
SUBROUTINE CLOSE_RUN_REPORT_FILE()

IF ( Report_Flags(iRF_Run) == 1 ) THEN
    CLOSE( UNIT = Report_IDs(iRF_Run) )
END IF

END SUBROUTINE CLOSE_RUN_REPORT_FILE

















 !+301+############################################################################!
!                                                                                   !
!                     OPEN_FRAME_REPORT_FILE                                          !
!                                                                                   !
 !#################################################################################!
SUBROUTINE OPEN_FRAME_REPORT_FILE(Frame)

INTEGER, INTENT(IN)                                     :: frame

CHARACTER(LEN = 70)                                     ::  FILE_NAME



IF ( myID_Poseidon == 0 ) THEN
    IF ( Report_Flags(iRF_Frame) == 1 ) THEN

        WRITE(FILE_NAME,'(A,A)') Poseidon_IterReports_Dir,"Frame_Report_SCRATCH.out"
        CALL OPEN_NEW_FILE( FILE_NAME, Report_IDs(iRF_Frame) )

    END IF ! Report_Flags(3)S
END IF ! myID_Poseidon == 0


END SUBROUTINE OPEN_FRAME_REPORT_FILE




!+302+##########################################################################!
!                                                                               !
!                   OUTPUT_ITERATION_REPORT                                     !
!                                                                               !
!###############################################################################!
SUBROUTINE OUTPUT_FRAME_REPORT()

!INTEGER, INTENT(IN)                 :: FRAME

INTEGER                                                 ::  FILE_ID
INTEGER                                                 ::  i
REAL(KIND = idp)                                        ::  r, theta, phi
REAL(KIND = idp), DIMENSION(0:ITER_REPORT_NUM_SAMPLES)  ::  x_e
REAL(KIND = idp), DIMENSION(1:ITER_REPORT_NUM_SAMPLES)  ::  x_c, dx_c


REAL(KIND = idp)                                        ::  Analytic_Val
REAL(KIND = idp)                                        ::  Return_Psi, Return_AlphaPsi
REAL(KIND = idp)                                        ::  Return_Beta1, Return_Beta2, Return_Beta3
REAL(KIND = idp)                                        ::  PsiPot_Val, AlphaPsiPot_Val

CHARACTER(LEN = 300)                                     ::  FILE_NAME

CHARACTER(LEN = 300)                                    ::  Line
INTEGER                                                 ::  io_stat
INTEGER                                                 ::  ITER_REPORT_NUM_SAMPLES = 20


120 FORMAT (12X,A)
123 FORMAT (12X,A38,ES22.15)

109 FORMAT (A,I2.2,A,I2.2)
110 FORMAT (11X,A1,16X,A18,9X,A13,10X,A18,10X,A11,14X,A11,14X,A11)
111 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15)

132 FORMAT (A,A,A,A)

142 FORMAT (16X,A,10X,A,15X,A)
143 FORMAT (19X,I2.2,15X,ES22.15,10X,ES22.15)




WRITE(FILE_NAME,132)Poseidon_IterReports_Dir,"Frame_Report_",trim(File_Suffix),".out"
CALL OPEN_NEW_FILE( trim(FILE_NAME), FILE_ID )

IF ( myID_Poseidon == 0 ) THEN
    ! Write Title to File
    IF ( Report_Flags(iRF_FRame) == 1 ) THEN

        WRITE(FILE_ID,109)"                                           Frame Report"
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
        WRITE(FILE_ID,123)"                    Initialize Time : ",Frame_Time_Table(1)
        WRITE(FILE_ID,123)" Input/Communicate Source Data Time : ",Frame_Time_Table(2)
        WRITE(FILE_ID,123)"     Input Boundary Conditions Time : ",Frame_Time_Table(3)
        WRITE(FILE_ID,123)"        CFA_3D_Apply_BCs_Part1 Time : ",Frame_Time_Table(4)
        WRITE(FILE_ID,120)"-------------------------------------------------------------"
        WRITE(FILE_ID,123)" ||     Calc_3D_Current_Values Time : ",Frame_Time_Table(5)
        WRITE(FILE_ID,123)" ||    CREATE_3D_SubJcbn_Terms Time : ",Frame_Time_Table(6)
        WRITE(FILE_ID,123)" ||       CREATE_3D_RHS_VECTOR Time : ",Frame_Time_Table(7)
        WRITE(FILE_ID,123)"\  /     CREATE_3D_JCBN_MATRIX Time : ",Frame_Time_Table(8)
        WRITE(FILE_ID,120)"-\/ ---------------------------------------------------------"
        WRITE(FILE_ID,123)"CREATE_3D_NONLAPLACIAN_STF_MAT Time : ",Frame_Time_Table(9)
        WRITE(FILE_ID,123)"REDUCE_3D_NONLAPLACIAN_STF_MAT Time : ",Frame_Time_Table(10)
        WRITE(FILE_ID,123)"FINISH_3D_NONLAPLACIAN_STF_MAT Time : ",Frame_Time_Table(11)
        WRITE(FILE_ID,123)"          FINISH_3D_RHS_VECTOR Time : ",Frame_Time_Table(12)
        WRITE(FILE_ID,123)"        CFA_3D_Apply_BCs_Part2 Time : ",Frame_Time_Table(13)
        WRITE(FILE_ID,123)"                    CFA_Solver Time : ",Frame_Time_Table(14)
        WRITE(FILE_ID,123)"        CFA_Coefficient_Update Time : ",Frame_Time_Table(15)
        WRITE(FILE_ID,123)"   CFA_Coefficient_Share_PETSc Time : ",Frame_Time_Table(16)
        WRITE(FILE_ID,123)"         CFA_Convergence_Check Time : ",Frame_Time_Table(17)
        WRITE(FILE_ID,123)"               Total Iteration Time : ",Frame_Time_Table(18)
        WRITE(FILE_ID,123)"             Poseidon_Dist_Sol Time : ",Frame_Time_Table(19)
        WRITE(FILE_ID,120)"============================================================="
        WRITE(FILE_ID,'(A)')" "
        WRITE(FILE_ID,'(A)')" "
        WRITE(FILE_ID,'(A)')" "
        WRITE(FILE_ID,'(A)')" "



        WRITE(FILE_ID,'(A)')"                                  Frame's Convergence"
        WRITE(FILE_ID,'(A)')"            ============================================================="
        WRITE(FILE_ID,'(A)')" "
        WRITE(FILE_ID,142)"Iteration","MAXVAL(ABS(Update_Vector))","Residual"
        DO i = 1,CUR_ITERATION-1
            WRITE(FILE_ID,143)i,MAXVAL(Frame_Update_Table(i,:)),MAXVAL(Frame_Residual_Table(i,Convergence_Type,:))
        END DO
        WRITE(FILE_ID,'(2/)')




        WRITE(FILE_ID,'(A)')"                                   Frame's Results"
        WRITE(FILE_ID,'(A)')"            ============================================================="
        WRITE(FILE_ID,'(A)')" "
        WRITE(FILE_ID,110)"r","Analytic Potential","Psi Potential","AlphaPsi Potential","Beta Value1","Beta Value2","Beta Value3"



        Call Create_Logarithmic_1D_Mesh( R_Inner, R_Outer, ITER_REPORT_NUM_SAMPLES,     &
                                         x_e, x_c, dx_c )
!        deltar = ( R_Outer - R_Inner )/ REAL(ITER_REPORT_NUM_SAMPLES, KIND = idp)

        DO i = 0,ITER_REPORT_NUM_SAMPLES

 !          r = i*deltar + R_Outer
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
        REWIND( UNIT=Report_IDs(iRF_Frame) )
        DO
            READ(Report_IDs(iRF_Frame),'(A)', iostat=io_stat) Line
            IF (io_stat /= 0 ) THEN
                EXIT
            END IF
            WRITE(FILE_ID, '(A)') Line
        END DO


    END IF ! Report_Flags(2) == 1

    CLOSE( UNIT = FILE_ID)
    CLOSE( UNIT = Report_IDs(iRF_Frame), STATUS='delete' )

END IF ! myID_Poseidon == 0

END SUBROUTINE OUTPUT_FRAME_REPORT





 !+303+############################################################################!
!                                                                                   !
!                    CLOSE_RUN_REPORT_FILE                                          !
!                                                                                   !
 !#################################################################################!
SUBROUTINE CLOSE_FRAME_REPORT_FILE()


CLOSE( UNIT = Report_IDs(iRF_Frame), STATUS='delete' )


END SUBROUTINE CLOSE_FRAME_REPORT_FILE






 




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

DO i = 1,Poseidon_Frame

    WRITE(FILE_ID,'(I3.3,A,I3.3)')i,"   ",Iteration_Histogram(i)

END DO

CLOSE( FILE_ID )

END SUBROUTINE OUTPUT_ITERATION_HISTORY






 !+701+############################################################################!
!                                                                                   !
!                     OUTPUT_PRIMATIVES                                      !
!                                                                                   !
 !#################################################################################!
SUBROUTINE OUTPUT_PRIMATIVES( Density, Velocity, RadLocs, Num_Entries )

REAL(KIND = idp), DIMENSION(1:Num_Entries), INTENT(IN)  :: Density
REAL(KIND = idp), DIMENSION(1:Num_Entries), INTENT(IN)  :: Velocity
REAL(KIND = idp), DIMENSION(1:Num_Entries), INTENT(IN)  :: RadLocs
INTEGER,                                    INTENT(IN)  :: Num_Entries

CHARACTER(LEN = 100), DIMENSION(:), ALLOCATABLE             ::  Filenames
INTEGER, DIMENSION(:), ALLOCATABLE                          ::  File_IDs
INTEGER                                                     ::  Num_Files = 3
INTEGER                                                     ::  i

116 FORMAT (A,A,A,A)

ALLOCATE( Filenames(1:Num_Files) )
ALLOCATE( File_IDs(1:Num_Files) )

WRITE(Filenames(1),116) Poseidon_Sources_Dir,"Primatives_Density_",trim(File_Suffix),".out"
WRITE(Filenames(2),116) Poseidon_Sources_Dir,"Primatives_Velocity_",trim(File_Suffix),".out"
WRITE(Filenames(3),116) Poseidon_Sources_Dir,"Primatives_rlocs_",trim(File_Suffix),".out"

! Open Files
File_IDs = [(161 + i, i =1,Num_Files)]
DO i = 1,Num_Files
    CALL OPEN_NEW_FILE( Filenames(i), File_IDs(i) )
END DO

! Write to Files
DO i = 1,Num_Entries
    WRITE(File_IDs(1),*)Density(i) /(Gram/Centimeter**3)
    WRITE(FILE_IDs(2),*)Velocity(i)/(Centimeter/Second)
    WRITE(FILE_IDs(3),*)RadLocs(i)/Centimeter
END DO

! Close Files
DO i = 1,Num_Files
    CLOSE( Unit = File_IDs(i))
END DO



END SUBROUTINE OUTPUT_PRIMATIVES






 !+701+############################################################################!
!                                                                                   !
!                     OUTPUT_YAHIL_PRIMATIVES                                      !
!                                                                                   !
 !#################################################################################!
SUBROUTINE OUTPUT_YAHIL_PRIMATIVES( Density, Velocity, Num_Entries )

REAL(KIND = idp), DIMENSION(1:Num_Entries), INTENT(IN)  :: Density
REAL(KIND = idp), DIMENSION(1:Num_Entries), INTENT(IN)  :: Velocity
INTEGER,                                    INTENT(IN)  :: Num_Entries

CHARACTER(LEN = 100), DIMENSION(:), ALLOCATABLE             ::  Filenames
INTEGER, DIMENSION(:), ALLOCATABLE                          ::  File_IDs
INTEGER                                                     ::  Num_Files = 2
INTEGER                                                     ::  i

116 FORMAT (A,A,A,A)

ALLOCATE( Filenames(1:Num_Files) )
ALLOCATE( File_IDs(1:Num_Files) )

WRITE(Filenames(1),116) Poseidon_Sources_Dir,"Sources_DX_",trim(File_Suffix),".out"
WRITE(Filenames(2),116) Poseidon_Sources_Dir,"Sources_VX_",trim(File_Suffix),".out"

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



END SUBROUTINE OUTPUT_YAHIL_PRIMATIVES












!+202+###########################################################################!
!                                                                                !
!                   OUTPUT_JACOBIAN_MATRIX                                       !
!                                                                                !
!################################################################################!
SUBROUTINE Output_Mesh( Mesh, row_in, flag )


INTEGER,                        INTENT(IN)                  ::  row_in
REAL(idp), DIMENSION(1:Row_In), INTENT(IN)               ::  Mesh


CHARACTER(LEN = 1 ), INTENT(IN), OPTIONAL                   ::  flag

CHARACTER(LEN = 300)                                        ::  FILE_NAME
CHARACTER(LEN = 300)                                        ::  FILE_NAMEb
CHARACTER(LEN = 40)                                         ::  fmt


INTEGER                                                     ::  rows

INTEGER                                                     ::  FILE_ID
INTEGEr                                                     ::  i

CHARACTER(LEN = 29)        :: Poseidon_Mesh_Dir      = "Poseidon_Output/Objects/Mesh/"

100 FORMAT (A,A,A,A,A,A)
101 FORMAT (A,A,A,A)

fmt = '(ES24.16E3,SP,ES24.16E3,"i")'


rows = size( Mesh, 1 )

IF ( present(flag) ) THEN
    WRITE(FILE_NAMEb,100) Poseidon_Mesh_Dir,"Mesh_Radial_Dim_",trim(File_Suffix),"_",flag,".out"
ELSE
    WRITE(FILE_NAMEb,101) Poseidon_Mesh_Dir,"Mesh_Radial_Dim_",trim(File_Suffix),".out"
END IF


CALL OPEN_NEW_FILE( trim(FILE_NAMEb), FILE_ID, 300 )
WRITE(FILE_ID,*) rows
CLOSE(FILE_ID)



IF ( present(flag) ) THEN
    WRITE(FILE_NAME,100) Poseidon_Mesh_Dir,"Mesh_Radial_Loc_",trim(File_Suffix),"_",flag,".out"
ELSE
    WRITE(FILE_NAME,100) Poseidon_Mesh_Dir,"Mesh_Radial_Loc_",trim(File_Suffix),".out"
END IF

CALL OPEN_NEW_FILE( trim(FILE_NAME), FILE_ID, 300 )
DO i = 1,rows
    WRITE(FILE_ID,fmt) Mesh(i)
END DO
CLOSE( FILE_ID )



END SUBROUTINE Output_Mesh




!+202+###########################################################################!
!                                                                                !
!                   Output_Nodal_Mesh                                       !
!                                                                                !
!################################################################################!
SUBROUTINE Output_Nodal_Mesh( Mesh, row_in, flag )


INTEGER,                        INTENT(IN)                  ::  row_in
REAL(idp), DIMENSION(1:Row_In), INTENT(IN)               ::  Mesh


CHARACTER(LEN = 1 ), INTENT(IN), OPTIONAL                   ::  flag

CHARACTER(LEN = 300)                                        ::  FILE_NAME
CHARACTER(LEN = 300)                                        ::  FILE_NAMEb
CHARACTER(LEN = 40)                                         ::  fmt

REAL(KIND = idp)                                            ::  DROT
REAL(KIND = idp), DIMENSION(0:Degree)                       ::  Local_Locations
REAL(KIND = idp), DIMENSION(0:Degree)                       ::  Cur_R_Locs

INTEGER                                                     ::  rows

INTEGER                                                     ::  FILE_ID
INTEGER                                                     ::  re, d

CHARACTER(LEN = 29)        :: Poseidon_Mesh_Dir      = "Poseidon_Output/Objects/Mesh/"

100 FORMAT (A,A,A,A,A,A)
101 FORMAT (A,A,A,A)

fmt = '(ES24.16E3,SP,ES24.16E3,"i")'


rows = size( Mesh, 1 )

IF ( present(flag) ) THEN
    WRITE(FILE_NAMEb,100) Poseidon_Mesh_Dir,"Mesh_Nodal_Dim_",trim(File_Suffix),"_",flag,".out"
ELSE
    WRITE(FILE_NAMEb,101) Poseidon_Mesh_Dir,"Mesh_Nodal_Dim_",trim(File_Suffix),".out"
END IF


CALL OPEN_NEW_FILE( trim(FILE_NAMEb), FILE_ID, 300 )
WRITE(FILE_ID,*) Num_R_Nodes
CLOSE(FILE_ID)




Local_Locations = Initialize_LGL_Quadrature_Locations(Degree)


IF ( present(flag) ) THEN
    WRITE(FILE_NAME,100) Poseidon_Mesh_Dir,"Mesh_Nodal_Loc_",trim(File_Suffix),"_",flag,".out"
ELSE
    WRITE(FILE_NAME,100) Poseidon_Mesh_Dir,"Mesh_Nodal_Loc_",trim(File_Suffix),".out"
END IF

CALL OPEN_NEW_FILE( trim(FILE_NAME), FILE_ID, 300 )
WRITE(File_ID,fmt) Mesh(1)
DO re = 1,NUM_R_ELEMENTS

    DROT = 0.50_idp*(Mesh(re+1) - Mesh(re))
    CUR_R_LOCS(:) = DROT * (Local_Locations(:)+1.0_idp) + Mesh(re)

    DO d = 1,Degree
        WRITE(FILE_ID,fmt) Cur_R_Locs(d)
    END DO
END DO
CLOSE( FILE_ID )



END SUBROUTINE Output_Nodal_Mesh



END MODULE Poseidon_IO_Module
