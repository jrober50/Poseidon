   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE IO_Write_Final_Results                                                !##!
!##!                                                                         !##!
!##!_________________________________________________________________________!##!
!##!                                                                         !##!
!##!                                                                         !##!
!##!=========================================================================!##!
!##!                                                                         !##!
!##!    Contains:                                                            !##!
!##!                                                                         !##!
!##!                                                                         !##!
!###############################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !##########################################################################!


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
                    Max_Iterations,         &
                    CFA_Eq_Flags

USE Variables_MPI, &
            ONLY :  myID_Poseidon

USE Variables_Derived, &
            ONLY :  Num_R_Nodes


USE Variables_Mesh, &
            ONLY :  NUM_R_ELEMENTS,         &
                    Num_T_Elements,         &
                    Num_P_Elements,         &
                    rlocs,                  &
                    tlocs,                  &
                    plocs,                  &
                    R_Inner,                &
                    R_Outer

USE Variables_IO, &
            ONLY :  Write_Flags,            &
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


USE Variables_External, &
            ONLY :  SelfSim_T

USE Variables_Functions, &
            ONLY :  Potential_Solution,         &
                    Shift_Solution,             &
                    Calc_3D_Values_at_Location, &
                    Calc_1D_CFA_Values

USE Poseidon_Return_Routines_Module, &
            ONLY :  Calc_Var_at_Location

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

USE Poseidon_File_Routines_Module, &
            ONLY :  Open_New_File,                  &
                    Open_Existing_File


USE Flags_IO_Module, &
            ONLY :  lPF_IO_Flags,               &
                    iPF_IO_Write_Results


IMPLICIT NONE



CONTAINS



 !+401+####################################################################!
!                                                                           !
!          Write_Final_Results                                              !
!                                                                           !
 !#########################################################################!
SUBROUTINE Write_Final_Results( Output_Locations_Flag, CFA_Eq_Overide )


INTEGER,    INTENT(IN), OPTIONAL                            ::  Output_Locations_Flag
INTEGER,    INTENT(IN), OPTIONAL,   DIMENSION(1:5)          ::  CFA_Eq_Overide


CHARACTER(LEN = 100), DIMENSION(:), ALLOCATABLE             ::  Filenames
INTEGER, DIMENSION(:), ALLOCATABLE                          ::  File_IDs
INTEGER                                                     ::  Num_Files


INTEGER                                                     ::  i, Here

INTEGER, DIMENSION(1:5)                                     ::  CFA_Eq_Flag_Used
INTEGER                                                     ::  Output_Locations_Mode

116 FORMAT (A,A,A,A,A,A)

IF ( PRESENT(CFA_Eq_Overide) ) THEN
    CFA_Eq_Flag_Used = CFA_Eq_Overide
ELSE
    CFA_Eq_Flag_Used = CFA_Eq_Flags
END IF

IF ( PRESENT(Output_Locations_Flag) ) THEN
    Output_Locations_Mode = Output_Locations_Flag
ELSE
    Output_Locations_Mode = 2  ! Default is radially nodal, angularly uniform.
END IF



IF ( lPF_IO_Flags(iPF_IO_Write_Results) ) THEN


    Num_Files = 4 + SUM(CFA_Eq_Flag_Used)+1

    ALLOCATE( Filenames(1:Num_Files) )
    ALLOCATE( File_IDs(1:Num_Files) )




    WRITE(Filenames(1),116) Poseidon_Results_Dir,"Results_Dimensions_",TRIM(File_Suffix),".out"
    WRITE(Filenames(2),116) Poseidon_Results_Dir,"Results_Radial_Locs_",TRIM(File_Suffix),".out"
    WRITE(Filenames(3),116) Poseidon_Results_Dir,"Results_Theta_Locs_",TRIM(File_Suffix),".out"
    WRITE(Filenames(4),116) Poseidon_Results_Dir,"Results_Phi_Locs_",TRIM(File_Suffix),".out"
    Here = 5
    DO i = 1,5
        IF ( CFA_Eq_Flag_Used(i) == 1 ) THEN
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


    IF ( .TRUE. ) THEN
        CALL Write_Final_Results_Samples( Num_Files, File_IDs, CFA_Eq_Flag_Used, Output_Locations_Mode )
    ELSE
!        iEL = [0, 0, 0]
!        iEU = [Num_R_Elements-1,Num_T_Elements-1,Num_P_Elements-1]
    END IF


    DO i = 1,Num_Files
        CLOSE( File_IDs(i))
    END DO

END IF

END SUBROUTINE Write_Final_Results







 !+401+####################################################################!
!                                                                           !
!          Write_Final_Results_Samples                                      !
!                                                                           !
 !#########################################################################!
SUBROUTINE Write_Final_Results_Samples( Num_Files, File_IDs, CFA_Eq_Flags,OL_Flag )

INTEGER,                            INTENT(IN)              ::  Num_Files
INTEGER, DIMENSION(1:Num_Files),    INTENT(IN)              ::  File_IDs
INTEGER, DIMENSION(1:5),            INTENT(IN)              ::  CFA_Eq_Flags
INTEGER,    INTENT(IN), OPTIONAL                            ::  OL_Flag


REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Output_re,          &
                                                                Output_rc,          &
                                                                Output_dr

INTEGER                                                     ::  i, j, k, Here

INTEGER                                                     ::  NUM_THETA_RAYS,     &
                                                                NUM_PHI_RAYS,       &
                                                                NUM_RADIAL_SAMPLES


REAL(KIND = idp)                                            ::  DELTA_THETA
REAL(KIND = idp)                                            ::  DELTA_PHI


REAL(KIND = idp), DIMENSION(:,:,:,:), ALLOCATABLE           ::  Var_Holder
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  R_Holder,           &
                                                                T_Holder,           &
                                                                P_Holder

REAL(idp)                                                   ::  DROT

REAL(idp), DIMENSION(0:Degree)                              ::  Node_X_Locs

REAL(idp), DIMENSION(6)                                     ::  Units
REAL(idp)                                                   ::  CF


IF ( OL_Flag == 1 ) THEN


    NUM_RADIAL_SAMPLES = WRITE_RESULTS_R_SAMPS
    NUM_THETA_RAYS = WRITE_RESULTS_T_SAMPS
    NUM_PHI_RAYS = WRITE_RESULTS_P_SAMPS

ELSE


    NUM_RADIAL_SAMPLES = Num_R_Elements*Degree + 1
    NUM_THETA_RAYS = Num_T_Elements
    NUM_PHI_RAYS = Num_P_Elements

END IF

! Allocate Data Holders !
ALLOCATE( Var_Holder(1:NUM_PHI_RAYS, 1:NUM_THETA_RAYS, 1:NUM_RADIAL_SAMPLES,1:6) )
ALLOCATE( R_Holder(1:NUM_RADIAL_SAMPLES) )
ALLOCATE( T_Holder(1:NUM_THETA_RAYS) )
ALLOCATE( P_Holder(1:NUM_PHI_RAYS) )




! Create Output Spacing
! Pull Number of Samples From Parameters !
IF ( OL_Flag == 1 ) THEN

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


    DO i = 1,NUM_RADIAL_SAMPLES
        R_Holder(i) = output_re(i)
    END DO ! i Loop

    DO j = 1,NUM_THETA_RAYS
        T_Holder(j) = (j-1)*DELTA_THETA
    END DO ! j

    DO k = 1,NUM_PHI_RAYS
        P_Holder(k) = k*DELTA_PHI
    END DO ! k



ELSE
    Node_X_Locs = Initialize_LGL_Quadrature_Locations(DEGREE)

    R_Holder(1) = R_Inner
    DO i = 0,Num_R_Elements-1
    DO j = 1,Degree
        k = i*Degree + j + 1
        DROT = 0.5_idp * (rlocs(i+1) - rlocs(i))

        R_Holder(k) = DROT * (Node_X_Locs(j)+1.0_idp) + rlocs(i)

    END DO
    END DO

    DO i = 0,Num_T_Elements-1
        T_Holder = 0.5_idp * (tlocs(i+1) + tlocs(i))
    END DO


    DO i = 0,Num_P_Elements-1
        P_Holder = 0.5_idp * (plocs(i+1) + plocs(i))
    END DO



END IF



! Calculate Output
DO k = 1,NUM_PHI_RAYS
DO j = 1,NUM_THETA_RAYS
DO i = 1,NUM_RADIAL_SAMPLES

    CF = Calc_Var_at_Location(R_Holder(i),T_Holder(j),P_Holder(k), iU_CF)
    
    Var_Holder(k,j,i,1) = CF
    Var_Holder(k,j,i,2) = Calc_Var_at_Location(R_Holder(i),T_Holder(j),P_Holder(k), iU_LF)/CF
    Var_Holder(k,j,i,3) = Calc_Var_at_Location(R_Holder(i),T_Holder(j),P_Holder(k), iU_S1, iVB_S)
    Var_Holder(k,j,i,4) = Calc_Var_at_Location(R_Holder(i),T_Holder(j),P_Holder(k), iU_S2, iVB_S)
    Var_Holder(k,j,i,5) = Calc_Var_at_Location(R_Holder(i),T_Holder(j),P_Holder(k), iU_S3, iVB_S)
    Var_Holder(k,j,i,6) = Calc_Var_at_Location(R_Holder(i),T_Holder(j),P_Holder(k), iU_X1, iVB_X)

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


END SUBROUTINE Write_Final_Results_Samples


























END MODULE IO_Write_Final_Results
