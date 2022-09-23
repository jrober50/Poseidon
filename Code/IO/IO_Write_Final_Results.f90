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
            ONLY :  Degree,                 &
                    L_Limit,                &
                    Eq_Flags

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

USE Poseidon_File_Routines_Module, &
            ONLY :  Open_New_File,                  &
                    Open_Existing_File

USE Flags_IO_Module, &
            ONLY :  lPF_IO_Flags,               &
                    iPF_IO_Write_Results

USE Flags_Core_Module, &
            ONLY :  iPF_Core_Flags,         &
                    iPF_Core_Method_Mode,   &
                    iPF_Core_Method_Newtonian


IMPLICIT NONE



CONTAINS



 !+401+####################################################################!
!                                                                           !
!          Write_Final_Results                                              !
!                                                                           !
 !#########################################################################!
SUBROUTINE Write_Final_Results( Output_Locations_Flag,      &
                                U_Flag_Option,              &
                                U_Overide,                  &
                                Kij_Flag_Option,            &
                                X_Flag_Option               )


INTEGER,    INTENT(IN), OPTIONAL                            ::  Output_Locations_Flag
INTEGER,    INTENT(IN), OPTIONAL                            ::  U_Flag_Option
INTEGER,    INTENT(IN), OPTIONAL,   DIMENSION(1:5)          ::  U_Overide
INTEGER,    INTENT(IN), OPTIONAL                            ::  Kij_Flag_Option
INTEGER,    INTENT(IN), OPTIONAL                            ::  X_Flag_Option

INTEGER, DIMENSION(4), ALLOCATABLE                          ::  mFile_IDs
INTEGER                                                     ::  mNum_Files

INTEGER, DIMENSION(:), ALLOCATABLE                          ::  uFile_IDs
INTEGER                                                     ::  uNum_Files

INTEGER, DIMENSION(:), ALLOCATABLE                          ::  xFile_IDs
INTEGER                                                     ::  xNum_Files

INTEGER, DIMENSION(:), ALLOCATABLE                          ::  kFile_IDs
INTEGER                                                     ::  kNum_Files


INTEGER                                                     ::  i, Here

INTEGER, DIMENSION(1:5)                                     ::  U_Flag_Used
INTEGER                                                     ::  Output_Locations_Mode

INTEGER                                                     ::  R_Dim
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  R_Holder

INTEGER                                                     ::  T_Dim
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  T_Holder

INTEGER                                                     ::  P_Dim
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  P_Holder

116 FORMAT (A,A,A,A,A,A)


IF ( lPF_IO_Flags(iPF_IO_Write_Results) ) THEN



    IF ( PRESENT(U_Flag_Option) ) THEN

        IF ( PRESENT(CFA_Eq_Overide) ) THEN
            U_Flag_Used = CFA_Eq_Overide
        ELSE
            IF ( iPF_Core_Flags(iPF_Core_Method_Mode) == iPF_Core_Method_Newtonian ) THEN
                U_Flag_Used = [1,0,0,0,0]
            ELSE
                U_Flag_Used = Eq_Flags
            END IF
        END IF

        uNum_Files  = SUM(CFA_Eq_Flag_Used)     ! Main metric variables

    ELSE
        uNum_Files  = 0
    END IF



    IF ( PRESENT(X_Flag_Option) ) THEN
        xNum_Files  = 1                ! Temporary
    ELSE
        xNum_Files  = 0
    END IF



    IF ( PRESENT(Kij_Flag_Option) ) THEN
        kNum_Files  = 1                ! Temporary
    ELSE
        kNum_Files  = 0
    END IF


    IF ( ANY((uNum_Files,xNum_Files,kNum_Files)==1 ) ) THEN
        mNum_Files  = 4
    ELSE
        mNum_Files  = 0
    END IF





    IF ( mNum_Files .NE. 0 ) THEN



        IF ( PRESENT(Output_Locations_Flag) ) THEN
            Output_Locations_Mode = Output_Locations_Flag
        ELSE
            Output_Locations_Mode = 2  ! Default is radially nodal, angularly uniform.
        END IF




        ALLOCATE( uFile_IDs(1:uNum_Files) )
        ALLOCATE( xFile_IDs(1:xNum_Files) )
        ALLOCATE( kFile_IDs(1:kNum_Files) )
        ALLOCATE( mFile_IDs(1:mNum_Files) )


        CALL Create_Final_Results_Files(uNum_Files, uFile_IDs,      &
                                        xNum_Files, xFile_IDs,      &
                                        kNum_Files, kFile_IDs,      &
                                        mNum_Files, mFile_IDs,      &
                                        U_Flag                      )




        !
        !   Write Results
        !
        !   Base Metric Variables
        IF ( .TRUE. ) THEN
            CALL Write_Final_Results_Samples( Num_Files, File_IDs, CFA_Eq_Flag_Used, Output_Locations_Mode )
        END IF



        !
        !   Close Files
        !
        CALL Close_Final_Results_Files( uNum_Files, uFile_IDs,      &
                                        xNum_Files, xFile_IDs,      &
                                        kNum_Files, kFile_IDs,      &
                                        mNum_Files, mFile_IDs       )

    END IF
END IF


END SUBROUTINE Write_Final_Results










 !+401+####################################################################!
!                                                                           !
!          Create_Final_Results_Files                                       !
!                                                                           !
 !#########################################################################!
SUBROUTINE Create_Final_Results_Files(  uNum_Files, uFile_IDs,      &
                                        xNum_Files, xFile_IDs,      &
                                        kNum_Files, kFile_IDs,      &
                                        mNum_Files, mFile_IDs,      &
                                        U_Flag                      )


INTEGER,                DIMENSION(1:uNum_Files),    INTENT(INOUT)   ::  uFile_IDs
INTEGER,                                            INTENT(IN)      ::  uNum_Files

INTEGER,                DIMENSION(1:uNum_Files),    INTENT(INOUT)   ::  xFile_IDs
INTEGER,                                            INTENT(IN)      ::  xNum_Files

INTEGER,                DIMENSION(1:uNum_Files),    INTENT(INOUT)   ::  kFile_IDs
INTEGER,                                            INTENT(IN)      ::  kNum_Files

INTEGER,                DIMENSION(1:uNum_Files),    INTENT(INOUT)   ::  mFile_IDs
INTEGER,                                            INTENT(IN)      ::  mNum_Files

INTEGER,                DIMENSION(1:5),             INTENT(IN)      ::  U_Flag


INTEGER                                                             ::  i

CHARACTER(LEN = 500),   DIMENSION(uNum_Files)                       ::  uFilenames
CHARACTER(LEN = 500),   DIMENSION(xNum_Files)                       ::  xFilenames
CHARACTER(LEN = 500),   DIMENSION(kNum_Files)                       ::  kFilenames
CHARACTER(LEN = 500),   DIMENSION(mNum_Files)                       ::  mFilenames

116 FORMAT (A,A,A,A,A,A)

!
!   Create Filenames
!


!   Dimension and Location Files
WRITE(mFilenames(1),116) Poseidon_Results_Dir,"Results_Dimensions_",TRIM(File_Suffix),".out"
WRITE(mFilenames(2),116) Poseidon_Results_Dir,"Results_Radial_Locs_",TRIM(File_Suffix),".out"
WRITE(mFilenames(3),116) Poseidon_Results_Dir,"Results_Theta_Locs_",TRIM(File_Suffix),".out"
WRITE(mFilenames(4),116) Poseidon_Results_Dir,"Results_Phi_Locs_",TRIM(File_Suffix),".out"



!   Base Metric Variables
DO i = 1,5
    IF ( U_Flag(i) == 1 ) THEN
        WRITE(uFilenames(Here),116)                &
                Poseidon_Results_Dir,"Results_",TRIM(CFA_ShortVars(i)),"_",TRIM(File_Suffix),".out"
    END IF
END DO



!   X Vector
IF ( xNum_Files .GE. 1 ) THEN

    DO i = 1,xNum_Files
        WRITE(xFilenames(i),116)                &
                Poseidon_Results_Dir,"Results_",TRIM(CFA_ShortVars(5+i)),"_",TRIM(File_Suffix),".out"
    END DO

    DO i = 1,xNum_Files
        CALL OPEN_NEW_FILE( xFilenames(i), xFile_IDs(i), 210 )
    END DO

END IF



!   Kij
IF ( kNum_Files .GE. 1 ) THEN

    DO i = 1,kNum_Files
        WRITE(kFilenames(i),116)                &
                Poseidon_Results_Dir,"Results_",TRIM(Kij_ShortVars(i)),"_",TRIM(File_Suffix),".out"
    END DO

    DO i = 1,kNum_Files
        CALL OPEN_NEW_FILE( kFilenames(i), kFile_IDs(i), 220 )
    END DO


END IF



END SUBROUTINE Create_Final_Results_Files






   
 !+401+####################################################################!
!                                                                           !
!          Close_Final_Results_Files                                        !
!                                                                           !
 !#########################################################################!
SUBROUTINE Close_Final_Results_Files(   uNum_Files, uFile_IDs,      &
                                        xNum_Files, xFile_IDs,      &
                                        kNum_Files, kFile_IDs,      &
                                        mNum_Files, mFile_IDs       )


INTEGER,                DIMENSION(1:uNum_Files),    INTENT(INOUT)   ::  uFile_IDs
INTEGER,                                            INTENT(IN)      ::  uNum_Files

INTEGER,                DIMENSION(1:uNum_Files),    INTENT(INOUT)   ::  xFile_IDs
INTEGER,                                            INTENT(IN)      ::  xNum_Files

INTEGER,                DIMENSION(1:uNum_Files),    INTENT(INOUT)   ::  kFile_IDs
INTEGER,                                            INTENT(IN)      ::  kNum_Files

INTEGER,                DIMENSION(1:uNum_Files),    INTENT(INOUT)   ::  mFile_IDs
INTEGER,                                            INTENT(IN)      ::  mNum_Files

INTEGER                                                             ::  i

!   Dimension and Location Files
DO i = 1,mNum_Files
    CLOSE ( mFiles_IDs(i) )
END DO


! Base Metric Files Files
DO i = 1,uNum_Files
    CLOSE( uFile_IDs(i))
END DO

!   X Vector
IF ( xNum_Files .GE. 1 ) THEN
    DO i = 1,xNum_Files
        CLOSE( xFile_IDs(i))
    END DO
END IF

!   Kij
IF ( kNum_Files .GE. 1 ) THEN
    DO i = 1,kNum_Files
        CLOSE( kFile_IDs(i))
    END DO
END IF


END SUBROUTINE Close_Final_Results_Files












 !+401+####################################################################!
!                                                                           !
!          Write_Final_Results_Samples                                      !
!                                                                           !
 !#########################################################################!
SUBROUTINE Write_Final_Results_Samples( Num_Files, File_IDs, CFA_Eq_Flags, OL_Flag )

INTEGER,                            INTENT(IN)              ::  Num_Files
INTEGER, DIMENSION(1:Num_Files),    INTENT(IN)              ::  File_IDs
INTEGER, DIMENSION(1:5),            INTENT(IN)              ::  CFA_Eq_Flags
INTEGER,    INTENT(IN), OPTIONAL                            ::  OL_Flag


REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Output_re,          &
                                                                Output_rc,          &
                                                                Output_dr

INTEGER                                                     ::  i, j, k, Here



REAL(KIND = idp), DIMENSION(:,:,:,:), ALLOCATABLE           ::  Var_Holder
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  R_Holder,           &
                                                                T_Holder,           &
                                                                P_Holder

REAL(idp)                                                   ::  DROT

REAL(idp), DIMENSION(0:Degree)                              ::  Node_X_Locs

REAL(idp), DIMENSION(6)                                     ::  Units
REAL(idp)                                                   ::  CF




! Allocate Data Holders !
ALLOCATE( Var_Holder(1:NUM_PHI_RAYS, 1:NUM_THETA_RAYS, 1:NUM_RADIAL_SAMPLES,1:6) )

ALLOCATE( R_Holder(1:NUM_RADIAL_SAMPLES) )
ALLOCATE( T_Holder(1:NUM_THETA_RAYS) )
ALLOCATE( P_Holder(1:NUM_PHI_RAYS) )





! Calculate Output
DO k = 1,NUM_PHI_RAYS
DO j = 1,NUM_THETA_RAYS
DO i = 1,NUM_RADIAL_SAMPLES


    IF ( ANY(CFA_Eq_Flags(1:2) == 1) ) THEN
        CF = Calc_Var_at_Location(R_Holder(i),T_Holder(j),P_Holder(k), iU_CF)
        
        Var_Holder(k,j,i,1) = CF
    END IF

    IF ( CFA_Eq_Flags(2) == 1 ) THEN
        Var_Holder(k,j,i,2) = Calc_Var_at_Location(R_Holder(i),T_Holder(j),P_Holder(k), iU_LF)/CF
    END IF

    IF ( CFA_Eq_Flags(3) == 1 ) THEN
        Var_Holder(k,j,i,3) = Calc_Var_at_Location(R_Holder(i),T_Holder(j),P_Holder(k), iU_S1, iVB_S)
    END IF


    IF ( CFA_Eq_Flags(4) == 1 ) THEN
        Var_Holder(k,j,i,4) = Calc_Var_at_Location(R_Holder(i),T_Holder(j),P_Holder(k), iU_S2, iVB_S)
    END IF

    IF ( CFA_Eq_Flags(5) == 1 ) THEN
        Var_Holder(k,j,i,5) = Calc_Var_at_Location(R_Holder(i),T_Holder(j),P_Holder(k), iU_S3, iVB_S)
    END IF

    IF ( .FALSE. ) THEN
        Var_Holder(k,j,i,6) = Calc_Var_at_Location(R_Holder(i),T_Holder(j),P_Holder(k), iU_X1, iVB_X)
    END IF

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


























 !+401+####################################################################!
!                                                                           !
!          Write_Final_Results_Samples                                      !
!                                                                           !
 !#########################################################################!
SUBROUTINE Write_Final_Results_All( Num_Files, File_IDs, CFA_Eq_Flags, OL_Flag )

INTEGER,                            INTENT(IN)              ::  Num_Files
INTEGER, DIMENSION(1:Num_Files),    INTENT(IN)              ::  File_IDs
INTEGER, DIMENSION(1:5),            INTENT(IN)              ::  CFA_Eq_Flags
INTEGER,    INTENT(IN), OPTIONAL                            ::  OL_Flag


REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Output_re,          &
                                                                Output_rc,          &
                                                                Output_dr

INTEGER                                                     ::  i, j, k, Here

INTEGER                                                     ::  R_Dim
INTEGER                                                     ::  T_Dim
INTEGER                                                     ::  P_Dim

REAL(KIND = idp), DIMENSION(:,:,:,:),   ALLOCATABLE         ::  uVar_Holder
REAL(KIND = idp), DIMENSION(:,:,:,:),   ALLOCATABLE         ::  xVar_Holder
REAL(KIND = idp), DIMENSION(:,:,:,:),   ALLOCATABLE         ::  kVar_Holder

REAL(KIND = idp), DIMENSION(:),         ALLOCATABLE         ::  R_Holder
REAL(KIND = idp), DIMENSION(:),         ALLOCATABLE         ::  T_Holder
REAL(KIND = idp), DIMENSION(:),         ALLOCATABLE         ::  P_Holder

REAL(idp), DIMENSION(6)                                     ::  Units
REAL(idp)                                                   ::  CF






CALL Create_Final_Results_Locations(OL_Flag,                        &
                                    R_Dim, T_Dim, P_Dim,            &
                                    R_Holder, T_Holder, P_Holder    )


ALLOCATE( uVar_Holder(1:P_Dim,1:T_Dim,1:R_Dim,5) )
ALLOCATE( xVar_Holder(1:P_Dim,1:T_Dim,1:R_Dim,xNum_Files) )
ALLOCATE( kVar_Holder(1:P_Dim,1:T_Dim,1:R_Dim,kNum_Files) )

! Calculate Output
DO k = 1,P_Dim
DO j = 1,T_Dim
DO i = 1,R_Dim


    IF ( ANY(CFA_Eq_Flags(1:2) == 1) ) THEN
        CF = Calc_Var_at_Location(R_Holder(i),T_Holder(j),P_Holder(k), iU_CF)
        uVar_Holder(k,j,i,1) = CF
    END IF

    IF ( CFA_Eq_Flags(2) == 1 ) THEN
        uVar_Holder(k,j,i,2) = Calc_Var_at_Location(R_Holder(i),T_Holder(j),P_Holder(k), iU_LF)/CF
    END IF

    IF ( CFA_Eq_Flags(3) == 1 ) THEN
        uVar_Holder(k,j,i,3) = Calc_Var_at_Location(R_Holder(i),T_Holder(j),P_Holder(k), iU_S1, iVB_S)
    END IF


    IF ( CFA_Eq_Flags(4) == 1 ) THEN
        uVar_Holder(k,j,i,4) = Calc_Var_at_Location(R_Holder(i),T_Holder(j),P_Holder(k), iU_S2, iVB_S)
    END IF

    IF ( CFA_Eq_Flags(5) == 1 ) THEN
        uVar_Holder(k,j,i,5) = Calc_Var_at_Location(R_Holder(i),T_Holder(j),P_Holder(k), iU_S3, iVB_S)
    END IF


END DO ! i Loop
END DO ! j Loop
END DO ! k Loop


IF (xNum_Files .NE. 0 ) THEN

    IF ( CFA_Eq_Flags(3) == 1 ) THEN
        xVar_Holder(k,j,i,1) = Calc_Var_at_Location(R_Holder(i),T_Holder(j),P_Holder(k), iU_S1, iVB_S)
    END IF

END IF



IF (kNum_Files .NE. 0 ) THEN

    

END IF





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




END SUBROUTINE Write_Final_Results_All




 !+401+####################################################################!
!                                                                           !
!          Create_Final_Results_Locations                                   !
!                                                                           !
 !#########################################################################!
SUBROUTINE Create_Final_Results_Locations(  OL_Flag,                        &
                                            R_Dim, T_Dim, P_Dim,            &
                                            R_Holder, T_Holder, P_Holder    )



INTEGER,                                            INTENT(IN)      ::  OL_Flag
INTEGER,                                            INTENT(INOUT)   ::  R_Dim
INTEGER,                                            INTENT(INOUT)   ::  T_Dim
INTEGER,                                            INTENT(INOUT)   ::  P_Dim

INTEGER,                DIMENSION(:), ALLOCATABLE,  INTENT(INOUT)   ::  R_Holder
INTEGER,                DIMENSION(:), ALLOCATABLE,  INTENT(INOUT)   ::  T_Holder
INTEGER,                DIMENSION(:), ALLOCATABLE,  INTENT(INOUT)   ::  P_Holder

REAL(idp)                                                   ::  DROT
REAL(idp)                                                   ::  dt
REAL(idp)                                                   ::  dp


REAL(idp), DIMENSION(0:Degree)                              ::  Node_X_Locs


IF ( OL_Flag == 1 ) THEN

    R_Dim = WRITE_RESULTS_R_SAMPS
    T_Dim = WRITE_RESULTS_T_SAMPS
    P_Dim = WRITE_RESULTS_P_SAMPS

ELSE

    R_Dim = Num_R_Elements*Degree + 1
    T_Dim = Num_T_Elements
    P_Dim = Num_P_Elements

END IF




IF ( OL_Flag == 1 ) THEN

    ! Create Radial Spacing !
    ALLOCATE( Output_rc(1:R_Dim) )
    ALLOCATE( Output_dr(1:R_Dim) )
    ALLOCATE( Output_re(0:R_Dim) )


    IF ( R_OUTER/(R_Inner+1.0_idp) > 1E3 ) THEN
        CALL Create_Logarithmic_1D_Mesh( R_INNER, R_OUTER, R_Dim,           &
                                         output_re, output_rc, output_dr    )
    ELSE
        CALL Create_Uniform_1D_Mesh( R_INNER, R_OUTER, R_Dim,           &
                                     output_re, output_rc, output_dr    )
    END IF

    !  Create Phi Spacing !
    IF ( P_Dim == 1 ) THEN
        dp = pi/2.0_idp
    ELSE
        dp = 2.0_idp*pi/(P_Dim-1)
    END IF

    ! Create Theta Spacing
    IF ( T_Dim == 1 ) THEN
        dt = pi
    ELSE
        dt = pi/(T_Dim-1)
    END IF




    DO i = 1,R_Dim
        R_Holder(i) = output_re(i)
    END DO ! i Loop

    DO j = 1,T_Dim
        T_Holder(j) = (j-1)*dt
    END DO ! j

    DO k = 1,P_Dim
        P_Holder(k) = k*dp
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
        T_Holder(i+1) = 0.5_idp * (tlocs(i+1) + tlocs(i))
    END DO


    DO i = 0,Num_P_Elements-1
        P_Holder(i+1) = 0.5_idp * (plocs(i+1) + plocs(i))
    END DO


END IF

END SUBROUTINE Create_Final_Results_Locations

















END MODULE IO_Write_Final_Results
