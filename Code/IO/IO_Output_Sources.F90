   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE IO_Output_Sources_Module                                              !##!
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
                    rlocs,                  &
                    tlocs,                  &
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
                    iPF_IO_Write_Sources

IMPLICIT NONE


CONTAINS



 !+101+################################################################!
!                                                                       !
!          Output_Poseidon_Sources_1D                                   !
!                                                                       !
 !#####################################################################!
SUBROUTINE Output_Poseidon_Sources_1D( Local_E, Local_S, Local_Si,                         &
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

IF ( lPF_IO_Flags(iPF_IO_Write_Sources) ) THEN
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

END SUBROUTINE Output_Poseidon_Sources_1D





 !+102+################################################################!
!                                                                       !
!          Output_Poseidon_Sources_3D                                   !
!                                                                       !
 !#####################################################################!
SUBROUTINE Output_Poseidon_Sources_3D( Local_E, Local_S, Local_Si,                         &
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
INTEGER                                                                 ::  iSF_Dim   = 4
INTEGER                                                                 ::  iSF_RLocs = 5
INTEGER                                                                 ::  iSF_TLocs = 6
INTEGER                                                                 ::  iSF_PLocs = 7
INTEGER                                                                 ::  iSF_Time  = 8



116 FORMAT (A,A,A,A)


IF ( lPF_IO_Flags(iPF_IO_Write_Sources) ) THEN
    Num_Files = 8

    ALLOCATE( Filenames(1:Num_Files) )
    ALLOCATE( File_IDs(1:Num_Files) )

    

    WRITE(Filenames(iSF_E),116) Poseidon_Sources_Dir,"Sources_E_",trim(File_Suffix),".out"
    WRITE(Filenames(iSF_S),116) Poseidon_Sources_Dir,"Sources_S_",trim(File_Suffix),".out"
    WRITE(Filenames(iSF_S1),116) Poseidon_Sources_Dir,"Sources_S1_",trim(File_Suffix),".out"
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

    

    WRITE(File_IDs(iSF_Dim),* )Local_RE_Dim, Local_TE_Dim, Local_PE_DIM
    WRITE(File_IDs(iSF_Dim),* )Local_RQ_Dim, Local_TQ_Dim, Local_PQ_DIM
!    WRITE(File_IDs(8),* )SELFSIM_T


    Delta_X = Right_Limit - Left_Limit
    


    DO pe = 0,Local_PE_Dim-1
    DO pd = 1,Local_PQ_Dim
    DO te = 0,Local_TE_Dim-1
    DO td = 1,Local_TQ_Dim
    DO re = 0,Local_RE_Dim-1

        DO rd = 1,Local_RQ_Dim

            tpd = (rd-1)*Local_TQ_Dim*Local_PQ_Dim      &
                + (td-1)*Local_PQ_Dim                   &
                + pd

            
            WRITE(File_IDs(iSF_E),*) Local_E(tpd,re,te,pe)/E_Units
            WRITE(File_IDs(iSF_S),*) Local_S(tpd,re,te,pe)/S_Units
            WRITE(File_IDs(iSF_S1),*) Local_Si(tpd,re,te,pe,1)/Si_Units

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
            WRITE(File_IDs(iSF_Rlocs),*) CUR_R_LOCS(rd)/Centimeter
        END DO
    END DO

    DO te = 0,Local_TE_Dim-1
        Dt_Over_Dx = (tlocs(te+1) - tlocs(te))/Delta_X
        CUR_T_LOCS(:) = Dt_Over_Dx * (Input_T_Quad(:)-Left_Limit) + tlocs(te)
        DO td = 1,Local_TQ_Dim
            WRITE(File_IDs(iSF_Tlocs),*) CUR_T_LOCS(td)
        END DO
    END DO

    DO pe = 0,Local_PE_Dim-1
    DO pd = 1,Local_PQ_Dim
        WRITE(File_IDs(iSF_Plocs),*) Input_P_Quad(pd)
    END DO
    END DO
    



    ! Close Files
    DO i = 1,Num_Files
        CLOSE( Unit = File_IDs(i))
    END DO

END IF




END SUBROUTINE Output_Poseidon_Sources_3D






 !+201+################################################################!
!                                                                       !
!          Output_Primatives                                            !
!                                                                       !
 !#####################################################################!
SUBROUTINE Output_Primatives( Density, Velocity, RadLocs, Num_Entries )

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






 !+202+################################################################!
!                                                                       !
!          Output_Yahil_Primatives                                      !
!                                                                       !
 !#####################################################################!
SUBROUTINE Output_Yahil_Primatives( Density, Velocity, Num_Entries )

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



END SUBROUTINE Output_Yahil_Primatives







END MODULE IO_Output_Sources_Module
