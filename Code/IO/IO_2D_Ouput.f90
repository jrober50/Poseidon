   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE IO_2D_Output                                                                 !##!
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
                    ONLY :  Iter_Time_Table,        &
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
                            Write_Flags


USE Variables_Yahil, &
                    ONLY :  SelfSim_T

USE Variables_Functions, &
                    ONLY :  Potential_Solution,         &
                            Shift_Solution,             &
                            Calc_3D_Values_at_Location, &
                            Calc_1D_CFA_Values

USE Functions_Quadrature, &
                    ONLY :  Initialize_LG_Quadrature_Locations


USE Functions_Mesh, &
                    ONLY :  Create_Logarithmic_1D_Mesh

USE Poseidon_IO_Parameters, &
                    ONLY :  Poseidon_Reports_Dir,                           &
                            Poseidon_IterReports_Dir,                       &
                            Poseidon_Objects_Dir,                           &
                            Poseidon_Mesh_Dir,                              &
                            Poseidon_LinSys_Dir,                            &
                            Poseidon_Results_Dir,                           &
                            Poseidon_Sources_Dir

USE Functions_Info,   &
                    ONLY  : PQ_ITERATIONS_MAX


IMPLICIT NONE


CHARACTER(LEN = 20), PARAMETER    :: Filename_Format_A = "(A,A)"
CHARACTER(LEN = 20), PARAMETER    :: Filename_Format_B = "(A,A,I2.2,A,I2.2,A)"


!*F&S*==========================================!
!                                               !
!           Functions & Subroutines             !
!                                               !
!===============================================!
CONTAINS






 !+401+############################################################################!
!                                                                                   !
!                           Output_2D_Results                                       !
!                                                                                   !
 !#################################################################################!
SUBROUTINE Output_2D_Results()


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

116 FORMAT (A,A,A,A)
117 FORMAT (A,A)


IF ( Write_Flags(5) == 1 ) THEN
    NUM_SAMPLES = 1000


    Num_Files = 9

    ALLOCATE( Filenames(1:Num_Files) )
    ALLOCATE( File_IDs(1:Num_Files) )


    WRITE(Filenames(1),116) Poseidon_Results_Dir,"Results_Lapse_",TRIM(File_Suffix),".out"
    WRITE(Filenames(2),116) Poseidon_Results_Dir,"Results_ConFactor_",TRIM(File_Suffix),".out"
    WRITE(Filenames(3),116) Poseidon_Results_Dir,"Results_Beta1_",TRIM(File_Suffix),".out"
    WRITE(Filenames(4),116) Poseidon_Results_Dir,"Results_Beta2_",TRIM(File_Suffix),".out"
    WRITE(Filenames(5),116) Poseidon_Results_Dir,"Results_Beta3_",TRIM(File_Suffix),".out"
    WRITE(Filenames(6),116) Poseidon_Results_Dir,"Results_Dimensions_",TRIM(File_Suffix),".out"
    WRITE(Filenames(7),116) Poseidon_Results_Dir,"Results_Radial_Locs_",TRIM(File_Suffix),".out"
    WRITE(Filenames(8),116) Poseidon_Results_Dir,"Results_Theta_Locs_",TRIM(File_Suffix),".out"
    WRITE(Filenames(9),116) Poseidon_Results_Dir,"Results_Phi_Locs_",TRIM(File_Suffix),".out"


    File_IDs = [(141 + i, i=1,Num_Files)]
    DO i = 1,Num_Files
        CALL OPEN_NEW_FILE( Filenames(i), File_IDs(i), 200 )
    END DO





    NUM_RADIAL_SAMPLES = WRITE_RESULTS_R_SAMPS

    ! Create Radial Spacing !
    ALLOCATE( Output_rc(1:NUM_RADIAL_SAMPLES) )
    ALLOCATE( Output_dr(1:NUM_RADIAL_SAMPLES) )
    ALLOCATE( Output_re(0:NUM_RADIAL_SAMPLES) )
    CALL Create_Logarithmic_1D_Mesh( R_INNER, R_OUTER, NUM_RADIAL_SAMPLES,  &
                                     output_re, output_rc, output_dr        )


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

                IF ( Return_Psi == 0.0_idp ) THEN
                    Lapse_Holder(k,j,i) = 0.0_idp
                ELSE
                    Lapse_Holder(k,j,i) = Return_AlphaPsi/Return_Psi
                END IF
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






END SUBROUTINE Output_2D_Results






END MODULE Poseidon_IO_Module

