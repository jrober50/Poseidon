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

USE Poseidon_Message_Routines_Module, &
            ONLY :  Warning_Message

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
                    
USE Poseidon_Interface_Return_Routines, &
            ONLY :  Poseidon_Return_All

USE Variables_MPI, &
            ONLY :  myID_Poseidon

USE Variables_Derived, &
            ONLY :  Num_R_Nodes


USE Variables_Mesh, &
            ONLY :  Num_R_Elements,         &
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
                    
USE Variables_Interface, &
            ONLY :  Caller_NQ,                      &
                    Caller_Quad_DOF,                &
                    Caller_xL,                      &
                    Caller_RQ_xlocs,                &
                    Caller_TQ_xlocs,                &
                    Caller_PQ_xlocs
                   
USE Variables_AMReX_Source, &
            ONLY :  iLeaf,                &
                    iTrunk
                   
                   
USE Variables_Tables, &
            ONLY :  Level_dx,                   &
                    Level_Ratios

USE Poseidon_Return_Routines_Module, &
            ONLY :  Calc_Var_at_Location

USE Return_Functions_FP , &
            ONLY :  Calc_Drv_At_Location_Type_B

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
                    CFA_ShortVars,                                  &
                    Kij_ShortVars

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
                    
                    
#ifdef POSEIDON_AMREX_FLAG
use amrex_base_module

USE amrex_box_module,   &
            ONLY :  amrex_box

USE amrex_boxarray_module, &
            ONLY :  amrex_boxarray

use amrex_fort_module, &
            ONLY :  amrex_spacedim

USE amrex_multifab_module,  &
            ONLY :  amrex_multifab,         &
                    amrex_multifab_build,   &
                    amrex_imultifab_build,  &
                    amrex_imultifab_destroy

USE Variables_AMReX_Core, &
            ONLY :  MF_Source,              &
                    AMReX_Num_Levels

USE Variables_AMReX_Source, &
            ONLY :  Source_PTR,             &
                    Mask_PTR,               &
                    iLeaf,                &
                    iTrunk

USE Poseidon_AMReX_MakeFineMask_Module, &
            ONLY :  AMReX_MakeFineMask

USE Poseidon_AMReX_Multilayer_Utilities_Module, &
            ONLY :  Find_Coarsest_Parent

#endif

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
LOGICAL,    INTENT(IN), OPTIONAL                            ::  U_Flag_Option
INTEGER,    INTENT(IN), OPTIONAL,   DIMENSION(1:5)          ::  U_Overide
LOGICAL,    INTENT(IN), OPTIONAL                            ::  Kij_Flag_Option
LOGICAL,    INTENT(IN), OPTIONAL                            ::  X_Flag_Option

INTEGER, DIMENSION(:), ALLOCATABLE                          ::  mFile_IDs
INTEGER                                                     ::  mNum_Files

INTEGER, DIMENSION(:), ALLOCATABLE                          ::  uFile_IDs
INTEGER                                                     ::  uNum_Files

INTEGER, DIMENSION(:), ALLOCATABLE                          ::  xFile_IDs
INTEGER                                                     ::  xNum_Files

INTEGER, DIMENSION(:), ALLOCATABLE                          ::  kFile_IDs
INTEGER                                                     ::  kNum_Files

INTEGER, DIMENSION(1:5)                                     ::  U_Flag_Used
INTEGER                                                     ::  Output_Locations_Mode



IF ( lPF_IO_Flags(iPF_IO_Write_Results) ) THEN


    IF ( PRESENT(U_Flag_Option) ) THEN
        IF( U_Flag_Option ) THEN

            IF ( PRESENT(U_Overide) ) THEN
                U_Flag_Used = U_Overide
            ELSE
                IF ( iPF_Core_Flags(iPF_Core_Method_Mode) == iPF_Core_Method_Newtonian ) THEN
                    U_Flag_Used = [1,0,0,0,0]
                ELSE
                    U_Flag_Used = Eq_Flags
                END IF
            END IF
        ELSE
            U_Flag_Used = [1,1,1,0,0]
        END IF
    ELSE
        U_Flag_Used = [1,1,1,0,0]
        
    END IF
    uNum_Files  = SUM(U_Flag_Used)     ! Main metric variables


    IF ( PRESENT(X_Flag_Option) ) THEN
        xNum_Files  = 0                ! Temporary
    ELSE
        xNum_Files  = 1
    END IF



    IF ( PRESENT(Kij_Flag_Option) ) THEN
        kNum_Files  = 0               ! Temporary
    ELSE
        kNum_Files  = 1
    END IF


    IF ( ANY( [uNum_Files,xNum_Files,kNum_Files].NE.0 ) ) THEN
        mNum_Files  = 4
    ELSE
        mNum_Files  = 0
    END IF


    IF ( mNum_Files .NE. 0 ) THEN


    
        IF ( PRESENT(Output_Locations_Flag) ) THEN
            Output_Locations_Mode = Output_Locations_Flag
        ELSE
            Output_Locations_Mode = 3  ! Default is radially nodal, angularly uniform.
        END IF

        Output_Locations_Mode = 2
        CALL Warning_Message('Output_Location_Mode Overridden in IO_Write_Final_Results.f90')

        ALLOCATE( uFile_IDs(1:uNum_Files) )
        ALLOCATE( xFile_IDs(1:xNum_Files) )
        ALLOCATE( kFile_IDs(1:kNum_Files) )
        ALLOCATE( mFile_IDs(1:mNum_Files) )


        CALL Create_Final_Results_Files(uNum_Files, uFile_IDs,      &
                                        xNum_Files, xFile_IDs,      &
                                        kNum_Files, kFile_IDs,      &
                                        mNum_Files, mFile_IDs,      &
                                        U_Flag_Used                 )




        !
        !   Write Results
        !
        !   Base Metric Variables

        
!        IF ( .TRUE. ) THEN
!            CALL Write_Final_Results_Samples( uNum_Files, uFile_IDs, U_Flag_Used, Output_Locations_Mode )
!        END IF
        CALL Write_Final_Results_All(   uNum_Files, uFile_IDs,      &
                                        xNum_Files, xFile_IDs,      &
                                        kNum_Files, kFile_IDs,      &
                                        mNum_Files, mFile_IDs,      &
                                        U_Flag_Used, Output_Locations_Mode )


        !
        !   Close Files
        !
        CALL Close_Final_Results_Files( uNum_Files, uFile_IDs,      &
                                        xNum_Files, xFile_IDs,      &
                                        kNum_Files, kFile_IDs,      &
                                        mNum_Files, mFile_IDs       )


        DEALLOCATE( uFile_IDs )
        DEALLOCATE( xFile_IDs )
        DEALLOCATE( kFile_IDs )
        DEALLOCATE( mFile_IDs )

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

INTEGER,                                            INTENT(IN)      ::  uNum_Files
INTEGER,                DIMENSION(1:uNum_Files),    INTENT(INOUT)   ::  uFile_IDs

INTEGER,                                            INTENT(IN)      ::  xNum_Files
INTEGER,                DIMENSION(1:xNum_Files),    INTENT(INOUT)   ::  xFile_IDs

INTEGER,                                            INTENT(IN)      ::  kNum_Files
INTEGER,                DIMENSION(1:kNum_Files),    INTENT(INOUT)   ::  kFile_IDs

INTEGER,                                            INTENT(IN)      ::  mNum_Files
INTEGER,                DIMENSION(1:mNum_Files),    INTENT(INOUT)   ::  mFile_IDs

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

DO i = 1,4
    CALL OPEN_NEW_FILE( mFilenames(i), mFile_IDs(i),200)
END DO

!   Base Metric Variables
DO i = 1,5
    IF ( U_Flag(i) == 1 ) THEN
        WRITE(uFilenames(i),116)                &
                Poseidon_Results_Dir,"Results_",TRIM(CFA_ShortVars(i)),"_",TRIM(File_Suffix),".out"

        CALL OPEN_NEW_FILE( uFilenames(i), uFile_IDs(i),205)
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

INTEGER,                DIMENSION(1:xNum_Files),    INTENT(INOUT)   ::  xFile_IDs
INTEGER,                                            INTENT(IN)      ::  xNum_Files

INTEGER,                DIMENSION(1:kNum_Files),    INTENT(INOUT)   ::  kFile_IDs
INTEGER,                                            INTENT(IN)      ::  kNum_Files

INTEGER,                DIMENSION(1:mNum_Files),    INTENT(INOUT)   ::  mFile_IDs
INTEGER,                                            INTENT(IN)      ::  mNum_Files

INTEGER                                                             ::  i

!   Dimension and Location Files
DO i = 1,mNum_Files
    CLOSE ( mFile_IDs(i) )
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
SUBROUTINE Write_Final_Results_All( uNum_Files, uFile_IDs,      &
                                    xNum_Files, xFile_IDs,      &
                                    kNum_Files, kFile_IDs,      &
                                    mNum_Files, mFile_IDs,      &
                                    CFA_Eq_Flags, OL_Flag )

INTEGER,                            INTENT(IN)              ::  uNum_Files
INTEGER,                            INTENT(IN)              ::  xNum_Files
INTEGER,                            INTENT(IN)              ::  kNum_Files
INTEGER,                            INTENT(IN)              ::  mNum_Files

INTEGER, DIMENSION(1:uNum_Files),    INTENT(IN)             ::  uFile_IDs
INTEGER, DIMENSION(1:xNum_Files),    INTENT(IN)             ::  xFile_IDs
INTEGER, DIMENSION(1:kNum_Files),    INTENT(IN)             ::  kFile_IDs
INTEGER, DIMENSION(1:mNum_Files),    INTENT(IN)             ::  mFile_IDs


INTEGER, DIMENSION(1:5),            INTENT(IN)              ::  CFA_Eq_Flags
INTEGER,                            INTENT(IN)              ::  OL_Flag

INTEGER                                                     ::  i, j, k

INTEGER                                                     ::  R_Dim
INTEGER                                                     ::  T_Dim
INTEGER                                                     ::  P_Dim

REAL(KIND = idp), DIMENSION(:,:,:,:),   ALLOCATABLE         ::  uVar_Holder
REAL(KIND = idp), DIMENSION(:,:,:,:),   ALLOCATABLE         ::  xVar_Holder
REAL(KIND = idp), DIMENSION(:,:,:,:,:), ALLOCATABLE         ::  xDrv_Holder
REAL(KIND = idp), DIMENSION(:,:,:,:),   ALLOCATABLE         ::  kVar_Holder

REAL(KIND = idp), DIMENSION(:),         ALLOCATABLE         ::  R_Holder
REAL(KIND = idp), DIMENSION(:),         ALLOCATABLE         ::  T_Holder
REAL(KIND = idp), DIMENSION(:),         ALLOCATABLE         ::  P_Holder

REAL(idp), DIMENSION(6)                                     ::  Units
REAL(idp)                                                   ::  CF

COMPLEX(idp)                                                ::  Tmp_Val
REAL(idp), DIMENSION(3)                                     ::  Tmp_Drv = 0.0_idp
REAL(idp), DIMENSION(3)                                     ::  Gamma

REAL(idp), DIMENSION(3,3,3)                                 ::  Christoffel
REAL(idp), DIMENSION(4)                                     ::  Reusable_Vals

Christoffel = 0.0_idp
Reusable_Vals=0.0_idp



CALL Create_Final_Results_Locations(OL_Flag,                        &
                                    R_Dim, T_Dim, P_Dim,            &
                                    R_Holder, T_Holder, P_Holder    )



ALLOCATE( uVar_Holder(1:P_Dim,1:T_Dim,1:R_Dim,5) )
ALLOCATE( xVar_Holder(1:P_Dim,1:T_Dim,1:R_Dim,3) )
ALLOCATE( xDrv_Holder(1:P_Dim,1:T_Dim,1:R_Dim,3,3) )
ALLOCATE( kVar_Holder(1:P_Dim,1:T_Dim,1:R_Dim,kNum_Files) )


xVar_Holder=0.0_idp
xDrv_Holder=0.0_idp

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


IF ( (xNum_Files .NE. 0 ) .OR. (kNum_Files .NE. 0 ) ) THEN

    DO k = 1,P_Dim
    DO j = 1,T_Dim
    DO i = 1,R_Dim
        xVar_Holder(k,j,i,1) = Calc_Var_at_Location(R_Holder(i),T_Holder(j),P_Holder(k), iU_X1, iVB_X)

        xVar_Holder(k,j,i,2) = 0.0_idp
        xVar_Holder(k,j,i,3) = 0.0_idp

    END DO ! i Loop
    END DO ! j Loop
    END DO ! k Loop

END IF


IF (kNum_Files .NE. 0 ) THEN

    DO k = 1,P_Dim
    DO j = 1,T_Dim
    DO i = 1,R_Dim

        CALL Calc_Drv_At_Location_Type_B(R_Holder(i),T_Holder(j),P_Holder(k), iU_X1, iVB_X, Tmp_Drv )
        xDrv_Holder(k,j,i,:,1) = Tmp_Drv

    END DO ! i Loop
    END DO ! j Loop
    END DO ! k Loop
    

    DO k = 1,P_Dim
    DO j = 1,T_Dim
    DO i = 1,R_Dim

        Gamma(1) = 1.0_idp
        Gamma(2) = 1.0_idp/(R_Holder(i)*R_Holder(i))
        Gamma(3) = Gamma(2) * 1.0_idp/( DSIN(T_Holder(j))*DSIN(T_Holder(j)) )

        Christoffel(1,2,2) = -R_Holder(i)
        Christoffel(1,3,3) = -R_Holder(i)*DSIN(T_Holder(j))*DSIN(T_Holder(j))

        Christoffel(2,1,2) = 1.0_idp/R_Holder(i)
        Christoffel(2,2,1) = 1.0_idp/R_Holder(i)
        Christoffel(2,3,3) = -DSIN(T_Holder(j))*DCOS(P_Holder(k))

        Christoffel(3,3,1) = 1.0_idp/R_Holder(i)
        Christoffel(3,1,3) = 1.0_idp/R_Holder(i)
        Christoffel(3,3,2) = 1.0_idp/DTAN(T_Holder(j))
        Christoffel(3,2,3) = 1.0_idp/DTAN(T_Holder(j))

        Reusable_Vals(1) = 2.0_idp/3.0_idp*(xDrv_Holder(k,j,i,1,1)+xDrv_Holder(k,j,i,2,2)+xDrv_Holder(k,j,i,3,3))
        Reusable_Vals(2) = 2.0_idp/3.0_idp*(Christoffel(1,1,1)+Christoffel(2,2,1)+Christoffel(3,3,1))
        Reusable_Vals(3) = 2.0_idp/3.0_idp*(Christoffel(1,1,2)+Christoffel(2,2,2)+Christoffel(3,3,2))
        Reusable_Vals(4) = 2.0_idp/3.0_idp*(Christoffel(1,1,3)+Christoffel(2,2,3)+Christoffel(3,3,3))

    
        Tmp_Val = Gamma(1)                                                                  &
                * ( 2.0_idp * xDrv_Holder(k,j,i,1,1) - Reusable_Vals(1)                     &
                  +(2.0_idp * Christoffel(1,1,1) - Reusable_Vals(2) )*xVar_Holder(k,j,i,1)  &
                  +(2.0_idp * Christoffel(1,1,2) - Reusable_Vals(3) )*xVar_Holder(k,j,i,2)  &
                  +(2.0_idp * Christoffel(1,1,3) - Reusable_Vals(4) )*xVar_Holder(k,j,i,3)  )

        kVar_Holder(k,j,i,1) = REAL( Tmp_Val/(Gamma(1)*Gamma(1)) )



    END DO ! i loop
    END DO ! j loop
    END DO ! k loop

END IF





! Write Output Location Files
WRITE(mFile_IDs(1),*)R_Dim, T_Dim, P_Dim
WRITE(mFile_IDs(2),*)R_Holder/Centimeter
WRITE(mFile_IDs(3),*)T_Holder
WRITE(mFile_IDs(4),*)P_Holder




Units = 1.0_idp
Units(3) = Shift_Units
! Write Output Value Files


DO i = 1,5
    IF ( CFA_Eq_Flags(i) == 1 ) THEN
        DO k = 1,P_Dim
        DO j = 1,T_Dim

            WRITE(uFile_IDs(i),*) uVar_Holder(k,j,:,i)/Units(i)

        END DO ! j Loop
        END DO ! k Loop
    END IF
END DO

DO k = 1,P_Dim
DO j = 1,T_Dim

    WRITE(xFile_IDs(1),*)xVar_Holder(k,j,:,1)/Units(3)

END DO ! j Loop
END DO ! k Loop



DO k = 1,P_Dim
DO j = 1,T_Dim

    WRITE(kFile_IDs(1),*)kVar_Holder(k,j,:,1)

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

REAL(idp),              DIMENSION(:), ALLOCATABLE,  INTENT(INOUT)   ::  R_Holder
REAL(idp),              DIMENSION(:), ALLOCATABLE,  INTENT(INOUT)   ::  T_Holder
REAL(idp),              DIMENSION(:), ALLOCATABLE,  INTENT(INOUT)   ::  P_Holder

REAL(idp)                                                   ::  DROT
REAL(idp)                                                   ::  dt
REAL(idp)                                                   ::  dp

REAL(idp),              DIMENSION(:), ALLOCATABLE           ::  Output_rc
REAL(idp),              DIMENSION(:), ALLOCATABLE           ::  Output_dr
REAL(idp),              DIMENSION(:), ALLOCATABLE           ::  Output_re

INTEGER                                                     ::  i, j, k

REAL(idp),              DIMENSION(:), ALLOCATABLE           ::  Node_x_Locs
!REAL(idp), DIMENSION(0:Degree)                              ::  Node_x_Locs


IF ( OL_Flag == 1 ) THEN

    R_Dim = WRITE_RESULTS_R_SAMPS
    T_Dim = WRITE_RESULTS_T_SAMPS
    P_Dim = WRITE_RESULTS_P_SAMPS

ELSE IF ( OL_Flag == 2 ) THEN

    R_Dim = Num_R_Elements*Degree
    T_Dim = Num_T_Elements
    P_Dim = Num_P_Elements

ELSE

    R_Dim = Num_R_Elements*Degree + 1
    T_Dim = Num_T_Elements
    P_Dim = Num_P_Elements

END IF



ALLOCATE( R_Holder(1:R_Dim) )
ALLOCATE( T_Holder(1:T_Dim) )
ALLOCATE( P_Holder(1:P_Dim) )




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


ELSE IF ( OL_Flag == 2 ) THEN

    ALLOCATE( Node_x_Locs(1:Degree))
    Node_x_Locs = Initialize_LG_Quadrature_Locations(DEGREE)


    
  
    DO i = 0,Num_R_Elements-1
    DO j = 1,Degree
        k = i*Degree + j
        DROT = 0.5_idp * (rlocs(i+1) - rlocs(i))

        R_Holder(k) = DROT * (Node_x_Locs(j)+1.0_idp) + rlocs(i)

    END DO
    END DO

    T_Holder = pi/2.0_idp
    P_Holder = pi

!    DO i = 0,Num_T_Elements-1
!        T_Holder(i+1) = 0.5_idp * (tlocs(i+1) + tlocs(i))
!    END DO
!
!    DO i = 0,Num_P_Elements-1
!        P_Holder(i+1) = 0.5_idp * (plocs(i+1) + plocs(i))
!    END DO

ELSE

    ALLOCATE( Node_x_Locs(0:Degree))
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

DEALLOCATE( Node_x_Locs )


END SUBROUTINE Create_Final_Results_Locations















 !+401+####################################################################!
!                                                                           !
!          Write_Final_Results                                              !
!                                                                           !
 !#########################################################################!
SUBROUTINE Write_Final_Results_AMReX( Output_Locations_Flag,      &
                                        U_Flag_Option,              &
                                        U_Overide,                  &
                                        Kij_Flag_Option,            &
                                        X_Flag_Option               )


INTEGER,    INTENT(IN), OPTIONAL                            ::  Output_Locations_Flag
LOGICAL,    INTENT(IN), OPTIONAL                            ::  U_Flag_Option
INTEGER,    INTENT(IN), OPTIONAL,   DIMENSION(1:5)          ::  U_Overide
LOGICAL,    INTENT(IN), OPTIONAL                            ::  Kij_Flag_Option
LOGICAL,    INTENT(IN), OPTIONAL                            ::  X_Flag_Option

INTEGER, DIMENSION(:), ALLOCATABLE                          ::  mFile_IDs
INTEGER                                                     ::  mNum_Files

INTEGER, DIMENSION(:), ALLOCATABLE                          ::  uFile_IDs
INTEGER                                                     ::  uNum_Files

INTEGER, DIMENSION(:), ALLOCATABLE                          ::  xFile_IDs
INTEGER                                                     ::  xNum_Files

INTEGER, DIMENSION(:), ALLOCATABLE                          ::  kFile_IDs
INTEGER                                                     ::  kNum_Files


INTEGER                                                     ::  i, Here
INTEGER                                                     ::  re, te, pe
INTEGER                                                     ::  rd, td, pd

INTEGER                                                     ::  R_Dim, T_Dim, P_Dim

REAL(idp), DIMENSION(3)                                     ::  Units

REAL(idp)                                                   ::  Quad_Span
REAL(idp), DIMENSION(1:Caller_NQ(1))                        ::  CUR_R_LOCS
REAL(idp), DIMENSION(1:Caller_NQ(1))                        ::  Cur_RX_Locs
REAL(idp), DIMENSION(1:Caller_NQ(2))                        ::  CUR_T_LOCS
REAL(idp), DIMENSION(1:Caller_NQ(2))                        ::  Cur_TX_Locs
REAL(idp), DIMENSION(1:Caller_NQ(3))                        ::  CUR_P_LOCS
REAL(idp), DIMENSION(1:Caller_NQ(3))                        ::  Cur_PX_Locs

INTEGER, DIMENSION(1:5)                                     ::  U_Flag_Used

REAL(idp)                                                   ::  DROT
REAL(idp)                                                   ::  DTOT
REAL(idp)                                                   ::  DPOT


TYPE(amrex_mfiter)                                          ::  mfi
TYPE(amrex_box)                                             ::  Box
TYPE(amrex_imultifab)                                       ::  Level_Mask
INTEGER, DIMENSION(3)                                       ::  iEL, iEU
INTEGER,    CONTIGUOUS, POINTER                             ::  Mask_PTR(:,:,:,:)
REAL(idp),  CONTIGUOUS, POINTER                             ::  Results_PTR(:,:,:,:)

INTEGER                                                     ::  Num_DOF
INTEGER                                                     ::  nComp
INTEGER                                                     ::  nLevels
INTEGER                                                     ::  lvl
INTEGER                                                     ::  MF_Results_nComps

INTEGER                                                     ::  iU_K11 = 6
INTEGER                                                     ::  iU_K12 = 7
INTEGER                                                     ::  iU_K13 = 8

TYPE(amrex_multifab),           ALLOCATABLE                 ::  MF_Results(:)
TYPE(amrex_boxarray),           ALLOCATABLE                 ::  BA_Results(:)
TYPE(amrex_distromap),          ALLOCATABLE                 ::  DM_Results(:)
TYPE(amrex_geometry),           ALLOCATABLE                 ::  GM_Results(:)

INTEGER                                                     ::  MF_Results_nVars    = 11
INTEGER                                                     ::  MF_Results_nGhosts  = 0




IF ( lPF_IO_Flags(iPF_IO_Write_Results) ) THEN


    Units = 1.0_idp
    Units(3) = Shift_Units

    Num_DOF = Caller_NQ(1)*Caller_NQ(2)*Caller_NQ(3)
    MF_Results_nComps = Num_DOF*MF_Results_nVars
    nLevels = AMReX_Num_Levels


    ALLOCATE( MF_Results(0:nLevels-1) )
    ALLOCATE( BA_Results(0:nLevels-1) )
    ALLOCATE( DM_Results(0:nLevels-1) )
    ALLOCATE( GM_Results(0:nLevels-1) )




    DO lvl = 0,nLevels-1
        CALL amrex_multifab_build(  MF_Results(lvl),            &
                                    MF_Source(Lvl)%BA,          &
                                    MF_Source(Lvl)%DM,          &
                                    MF_Results_nComps,          &
                                    MF_Results_nGhosts          )

        CALL MF_Results(lvl)%SetVal(0.0_idp)
    END DO



    


    IF ( PRESENT(U_Flag_Option) ) THEN
        IF( U_Flag_Option ) THEN

            IF ( PRESENT(U_Overide) ) THEN
                U_Flag_Used = U_Overide
            ELSE
                IF ( iPF_Core_Flags(iPF_Core_Method_Mode) == iPF_Core_Method_Newtonian ) THEN
                    U_Flag_Used = [1,0,0,0,0]
                ELSE
                    U_Flag_Used = Eq_Flags
                END IF
            END IF
        ELSE
            U_Flag_Used = [1,1,1,0,0]
        END IF
    ELSE
        U_Flag_Used = [1,1,1,0,0]
        
    END IF
    uNum_Files  = SUM(U_Flag_Used)     ! Main metric variables


    IF ( PRESENT(X_Flag_Option) ) THEN
        xNum_Files  = 0                ! Temporary
    ELSE
        xNum_Files  = 2
    END IF



    IF ( PRESENT(Kij_Flag_Option) ) THEN
        kNum_Files  = 0               ! Temporary
    ELSE
        kNum_Files  = 1
    END IF


    IF ( ANY( [uNum_Files,xNum_Files,kNum_Files].NE.0 ) ) THEN
        mNum_Files  = 4
    ELSE
        mNum_Files  = 0
    END IF



    IF ( mNum_Files .NE. 0 ) THEN


    
    
        ALLOCATE( uFile_IDs(1:uNum_Files) )
        ALLOCATE( xFile_IDs(1:xNum_Files) )
        ALLOCATE( kFile_IDs(1:kNum_Files) )
        ALLOCATE( mFile_IDs(1:mNum_Files) )


        CALL Create_Final_Results_Files(uNum_Files, uFile_IDs,      &
                                        xNum_Files, xFile_IDs,      &
                                        kNum_Files, kFile_IDs,      &
                                        mNum_Files, mFile_IDs,      &
                                        U_Flag_Used                 )



        ! Write Output Location Files
        R_Dim = Num_R_Elements*Degree + 1
        T_Dim = Num_T_Elements
        P_Dim = Num_P_Elements
        
        WRITE(mFile_IDs(1),*)R_Dim, T_Dim, P_Dim


        !
        !   Write Results
        !
        !   Base Metric Variables

        CALL Poseidon_Return_ALL( MF_Results )
        
        
        
        Quad_Span = Caller_xL(2) - Caller_xL(1)

        CUR_RX_LOCS = 2.0_idp * ( Caller_RQ_xlocs(:) - Caller_xL(1) )/Quad_Span - 1.0_idp
        CUR_TX_LOCS = 2.0_idp * ( Caller_TQ_xlocs(:) - Caller_xL(1) )/Quad_Span - 1.0_idp
        CUR_PX_LOCS = 2.0_idp * ( Caller_PQ_xlocs(:) - Caller_xL(1) )/Quad_Span - 1.0_idp

         
        DO lvl = nLevels-1,0,-1

            DROT = 0.5_idp * Level_dx(lvl,1)
            DTOT = 0.5_idp * Level_dx(lvl,2)
            DPOT = 0.5_idp * Level_dx(lvl,3)
            
            !
            !   MakeFineMask
            !
            IF ( lvl < nLevels-1 ) THEN
                CALL AMReX_MakeFineMask(  Level_Mask,               &
                                          MF_Results(lvl)%ba,        &
                                          MF_Results(lvl)%dm,        &
                                          MF_Results(lvl+1)%ba,      &
                                          iLeaf, iTrunk            )
            ELSE
                ! Create Level_Mask all equal to 1
                CALL amrex_imultifab_build( Level_Mask,             &
                                            MF_Results(lvl)%ba,      &
                                            MF_Results(lvl)%dm,      &
                                            1,                      &
                                            0                       )
                CALL Level_Mask%SetVal(iLeaf)
            END IF


            CALL amrex_mfiter_build(mfi, MF_Results(lvl), tiling = .true. )

            DO WHILE(mfi%next())

                Results_PTR => MF_Results(lvl)%dataPtr(mfi)
                Mask_PTR   => Level_Mask%dataPtr(mfi)

                Box = mfi%tilebox()
                nComp =  MF_Results(lvl)%ncomp()

                iEL = Box%lo
                iEU = Box%hi



                DO re = iEL(1),iEU(1)
                DO te = iEL(2),iEU(2)
                DO pe = iEL(3),iEU(3)

                    IF ( Mask_PTR(RE,TE,PE,1) == iLeaf ) THEN


                    Cur_R_Locs(:) = DROT * (CUR_RX_LOCS(:)+1.0_idp + 2.0_idp*re)
                    Cur_T_Locs(:) = DTOT * (CUR_TX_LOCS(:)+1.0_idp + 2.0_idp*te)
                    Cur_P_Locs(:) = DPOT * (CUR_PX_LOCS(:)+1.0_idp + 2.0_idp*pe)


                    DO rd = 1,Caller_NQ(1)
                        
                        WRITE(mFile_IDs(2),*)Cur_R_Locs(rd)/Centimeter
                    END DO
                    DO td = 1,Caller_NQ(2)
                        WRITE(mFile_IDs(3),*)Cur_T_Locs(td)
                    END DO
                    DO pd = 1,Caller_NQ(3)
                        WRITE(mFile_IDs(4),*)Cur_P_Locs(pd)
                    END DO




                    DO pd = 1,Caller_NQ(3)
                    DO td = 1,Caller_NQ(2)
                    DO rd = 1,Caller_NQ(1)


            
                        DO i = 1,3
                            Here = AMReX_nCOMP_Map( i, rd, td, pd, Caller_NQ )
                            WRITE(uFile_IDs(i),*) Results_PTR(re,te,pe,Here)/Units(i)
                        END DO

                        ! K_11
                        Here = AMReX_nCOMP_Map( iU_K11, rd, td, pd, Caller_NQ )
                        WRITE(kFile_IDs(1),*) Results_PTR(re,te,pe,Here)
                        
                        ! X1 Value
                        Here = AMReX_nCOMP_Map( iU_K12, rd, td, pd, Caller_NQ )
                        WRITE(xFile_IDs(1),*)Results_PTR(re,te,pe,Here)
                        
                        ! X1 Deriv
                        Here = AMReX_nCOMP_Map( iU_K13, rd, td, pd, Caller_NQ )
                        WRITE(xFile_IDs(2),*)Results_PTR(re,te,pe,Here)

                    END DO ! rd Loop
                    END DO ! td Loop
                    END DO ! pd Loop

                    END IF

                END DO ! pe
                END DO ! te
                END DO ! re

            END DO

            CALL amrex_mfiter_destroy(mfi)
            CALL amrex_imultifab_destroy( Level_Mask )

        END DO ! lvl



        !
        !   Close Files
        !
        CALL Close_Final_Results_Files( uNum_Files, uFile_IDs,      &
                                        xNum_Files, xFile_IDs,      &
                                        kNum_Files, kFile_IDs,      &
                                        mNum_Files, mFile_IDs       )


        DEALLOCATE( uFile_IDs )
        DEALLOCATE( xFile_IDs )
        DEALLOCATE( kFile_IDs )
        DEALLOCATE( mFile_IDs )

    

    END IF
    
    DEALLOCATE( MF_Results )
    DEALLOCATE( BA_Results )
    DEALLOCATE( DM_Results )
    DEALLOCATE( GM_Results )
    
END IF


END SUBROUTINE Write_Final_Results_AMReX





 !+202+########################################################!
!                                                               !
!       Poseidon_Return_AMReX_Extrinsic_Curvature               !
!                                                               !
 !#############################################################!
PURE FUNCTION AMReX_nCOMP_Map( iU, rd, td, pd, NQ )

INTEGER, INTENT(IN)                         ::  iU
INTEGER, INTENT(IN)                         ::  rd
INTEGER, INTENT(IN)                         ::  td
INTEGER, INTENT(IN)                         ::  pd
INTEGER, DIMENSION(3),  INTENT(IN)          ::  NQ

INTEGER                                     ::  AMReX_nCOMP_Map

INTEGER                                     ::  Here
INTEGER                                     ::  Num_QP

! CF = 1
! LF = 2
! S1 = 3
! S2 = 4
! S3 = 5
! K11 = 6
! K12 = 7
! K13 = 8
! K22 = 9
! K23 = 10
! K33 = 11

Here = (pd-1)*NQ(1)*NQ(2)   &
     + (td-1)*NQ(1)         &
     + rd
Num_QP = NQ(1)*NQ(2)*NQ(3)

AMReX_nCOMP_Map = (iU-1)*Num_QP + Here


END FUNCTION AMReX_nCOMP_Map



END MODULE IO_Write_Final_Results
