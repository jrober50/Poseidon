   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Driver_ReturnTest                                                     !##!
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
            ONLY :  idp

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag

USE Parameters_Variable_Indices, &
            ONLY :  iVB_X,                      &
                    iVB_S,                      &
                    iU_CF,                      &
                    iU_LF,                      &
                    iU_S1,                      &
                    iU_S2,                      &
                    iU_S3,                      &
                    iU_X1,                      &
                    iU_X2,                      &
                    iU_X3


USE amrex_base_module

USE amrex_box_module, &
            ONLY :  amrex_box

USE amrex_boxarray_module,  &
            ONLY :  amrex_boxarray,                 &
                    amrex_boxarray_build,           &
                    amrex_boxarray_destroy

USE amrex_distromap_module, &
            ONLY :  amrex_distromap,                &
                    amrex_distromap_build,          &
                    amrex_distromap_destroy

USE amrex_multifab_module,  &
            ONLY :  amrex_multifab,                 &
                    amrex_multifab_build

USE amrex_amrcore_module,   &
            ONLY :  amrex_amrcore_init,             &
                    amrex_init_virtual_functions,   &
                    amrex_init_from_scratch,        &
                    amrex_ref_ratio,                &
                    amrex_max_level

USE amrex_geometry_module,  &
            ONLY :  amrex_geometry


USE Poseidon_AMReX_MakeFineMask_Module, &
            ONLY :  AMReX_MakeFineMask

USE Functions_Quadrature, &
            ONLY :  Initialize_LG_Quadrature_Locations, &
                    Initialize_Trapezoid_Quadrature_Locations

USE Variables_AMReX_Source, &
            ONLY :  iLeaf,                &
                    iTrunk

USE Poseidon_Return_Routines_Module, &
            ONLY :  Poseidon_Return_Conformal_Factor,    &
                    Poseidon_Return_All

USE Variables_Quadrature, &
            ONLY :  NUM_R_QUAD_POINTS,          &
                    NUM_T_QUAD_POINTS,          &
                    NUM_P_QUAD_POINTS,          &
                    NUM_TP_QUAD_POINTS,         &
                    INT_R_LOCATIONS,            &
                    INT_T_LOCATIONS,            &
                    INT_P_LOCATIONS

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
                    File_Suffix

USE Poseidon_IO_Parameters, &
            ONLY :  Poseidon_Reports_Dir,                           &
                    Poseidon_IterReports_Dir,                       &
                    Poseidon_Objects_Dir,                           &
                    Poseidon_Mesh_Dir,                              &
                    Poseidon_LinSys_Dir,                            &
                    Poseidon_Results_Dir,                           &
                    Poseidon_Sources_Dir,                           &
                    CFA_ShortVars

USE Variables_Mesh, &
            ONLY :  Num_P_Elements,             &
                    rlocs,                      &
                    drlocs,                     &
                    tlocs,                      &
                    plocs

USE Maps_Domain, &
            ONLY :  Map_To_FEM_Node,            &
                    FEM_Elem_Map

USE Variables_Tables, &
            ONLY :  Ylm_CC_Values,              &
                    Ylm_Elem_CC_Values,         &
                    Lagrange_Poly_Table,        &
                    Lagpoly_MultiLayer_Table,   &
                    Level_dx,                   &
                    Level_Ratios


USE Poseidon_File_Routines_Module, &
            ONLY :  Open_New_File,                  &
                    Open_Existing_File

IMPLICIT NONE

INTEGER                 ::  iU_K11 = 6
INTEGER                 ::  iU_K12 = 7
INTEGER                 ::  iU_K13 = 8
INTEGER                 ::  iU_K22 = 9
INTEGER                 ::  iU_K23 = 10
INTEGER                 ::  iU_K33 = 11


CONTAINS



 !+101+####################################################!
!                                                           !
!          Return_Test	                                    !
!                                                           !
 !#########################################################!
SUBROUTINE Return_Test( nLevels, NQ, MF_Source )

INTEGER,                INTENT(IN)              ::  nLevels
INTEGER, DIMENSION(3),  INTENT(IN)              ::  NQ
TYPE(amrex_multifab),   INTENT(IN)              ::  MF_Source(0:nLevels-1)



INTEGER                                         ::  MF_Results_nVars    = 11
INTEGER                                         ::  MF_Results_nGhosts  = 0


INTEGER                                         ::  Num_Quad
INTEGER                                         ::  MF_Results_nComps


REAL(idp),  DIMENSION(NQ(1))                    ::  R_Quad
REAL(idp),  DIMENSION(NQ(2))                    ::  T_Quad
REAL(idp),  DIMENSION(NQ(3))                    ::  P_Quad
REAL(idp),  DIMENSION(2)                        ::  xL

TYPE(amrex_multifab),           ALLOCATABLE     ::  MF_Results(:)
TYPE(amrex_boxarray),           ALLOCATABLE     ::  BA_Results(:)
TYPE(amrex_distromap),          ALLOCATABLE     ::  DM_Results(:)
TYPE(amrex_geometry),           ALLOCATABLE     ::  GM_Results(:)

INTEGER                                         ::  lvl



ALLOCATE( MF_Results(0:nLevels-1) )
ALLOCATE( BA_Results(0:nLevels-1) )
ALLOCATE( DM_Results(0:nLevels-1) )
ALLOCATE( GM_Results(0:nLevels-1) )


Num_Quad = NQ(1)*NQ(2)*NQ(3)
MF_Results_nComps = Num_Quad*MF_Results_nVars

DO lvl = 0,nLevels-1

    CALL amrex_multifab_build(  MF_Results(lvl),            &
                                MF_Source(Lvl)%BA,          &
                                MF_Source(Lvl)%DM,          &
                                MF_Results_nComps,          &
                                MF_Results_nGhosts          )

    CALL MF_Results(lvl)%SetVal(0.0_idp)
END DO


xL(1) = -1.0_idp
xL(2) = +1.0_idp
R_Quad = Initialize_LG_Quadrature_Locations(NQ(1))
T_Quad = Initialize_LG_Quadrature_Locations(NQ(2))
P_Quad = Initialize_Trapezoid_Quadrature_Locations(NQ(3))

!CALL Poseidon_Return_Conformal_Factor(NQ,                   &
!                                      R_Quad,               &
!                                      T_Quad,               &
!                                      P_Quad,               &
!                                      xL(1),                &
!                                      xL(2),                &
!                                      nLevels,              &
!                                      MF_Results            )

CALL Poseidon_Return_Conformal_Factor( MF_Results )



!CALL Output_Variable( nLevels, NQ, MF_Results, iU_CF )


!CALL Poseidon_Return_ALL(   NQ,                   &
!                            R_Quad,               &
!                            T_Quad,               &
!                            P_Quad,               &
!                            xL(1),                &
!                            xL(2),                &
!                            nLevels,              &
!                            MF_Results            )



CALL Poseidon_Return_ALL( MF_Results )

 

CALL Output_All_Variables( nLevels, NQ, MF_Results )

!CALL Write_Final_Results_Kij( nLevels, NQ, MF_Results )


END SUBROUTINE Return_Test










 !+101+####################################################!
!                                                           !
!          Output_ConFactor                                 !
!                                                           !
 !#########################################################!
SUBROUTINE Output_Variable( nLevels, NQ, MF_Results, iU )

INTEGER,                INTENT(IN)              ::  nLevels
INTEGER, DIMENSION(3),  INTENT(IN)              ::  NQ
TYPE(amrex_multifab),   INTENT(IN)              ::  MF_Results(0:nLevels-1)
INTEGER,                INTENT(IN)              ::  iU

INTEGER                                         ::  lvl
TYPE(amrex_mfiter)                              ::  mfi
TYPE(amrex_box)                                 ::  Box
TYPE(amrex_imultifab)                           ::  Level_Mask
INTEGER                                         ::  nComp
INTEGER, DIMENSION(3)                           ::  iEL, iEU
INTEGER,    CONTIGUOUS, POINTER                 ::  Mask_PTR(:,:,:,:)
REAL(idp),  CONTIGUOUS, POINTER                 ::  Results_PTR(:,:,:,:)

INTEGER                                         ::  re, te, pe
INTEGER                                         ::  rd, td, pd
INTEGER                                         ::  Here
INTEGER                                         ::  Num_Quad

Num_Quad = NQ(1)*NQ(2)*NQ(3)


DO lvl = nLevels-1,0,-1

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

            DO pd = 1,NQ(3)
            DO td = 1,NQ(2)
            DO rd = 1,NQ(1)


                Here = (iU-1) * Num_Quad        &
                     + (pd-1)*NQ(1)*NQ(2)       &
                     + (td-1)*NQ(1)             &
                     + rd

                PRINT*,re,te,pe,rd,td,pd,iU,Results_PTR(re,te,pe,Here)

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


END SUBROUTINE Output_Variable




!+101+####################################################!
!                                                           !
!          Output_ConFactor                                 !
!                                                           !
 !#########################################################!
SUBROUTINE Output_ALL_Variables( nLevels, NQ, MF_Results )

INTEGER,                INTENT(IN)              ::  nLevels
INTEGER, DIMENSION(3),  INTENT(IN)              ::  NQ
TYPE(amrex_multifab),   INTENT(IN)              ::  MF_Results(0:nLevels-1)

INTEGER                                         ::  lvl
TYPE(amrex_mfiter)                              ::  mfi
TYPE(amrex_box)                                 ::  Box
TYPE(amrex_imultifab)                           ::  Level_Mask
INTEGER                                         ::  nComp
INTEGER, DIMENSION(3)                           ::  iEL, iEU
INTEGER,    CONTIGUOUS, POINTER                 ::  Mask_PTR(:,:,:,:)
REAL(idp),  CONTIGUOUS, POINTER                 ::  Results_PTR(:,:,:,:)

INTEGER                                         ::  re, te, pe
INTEGER                                         ::  rd, td, pd
INTEGER                                         ::  Here
INTEGER                                         ::  Num_Quad
INTEGER                                         ::  iU


Num_Quad = NQ(1)*NQ(2)*NQ(3)


DO lvl = nLevels-1,0,-1

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

            DO pd = 1,NQ(3)
            DO td = 1,NQ(2)
            DO rd = 1,NQ(1)




                PRINT*,re,te,pe,rd,td,pd
                DO iU = 1,11
                    Here = (iU-1) * Num_Quad        &
                         + (pd-1)*NQ(1)*NQ(2)       &
                         + (td-1)*NQ(1)             &
                         + rd
                    PRINT*,Results_PTR(re,te,pe,Here)
                END DO
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



END SUBROUTINE Output_ALL_Variables





 !+401+####################################################################!
!                                                                           !
!          Write_Final_Results                                              !
!                                                                           !
 !#########################################################################!
SUBROUTINE Write_Final_Results_Kij( nLevels, NQ, MF_Results )

INTEGER,                INTENT(IN)              ::  nLevels
INTEGER, DIMENSION(3),  INTENT(IN)              ::  NQ
TYPE(amrex_multifab),   INTENT(IN)              ::  MF_Results(0:nLevels-1)

INTEGER                                         ::  lvl
TYPE(amrex_mfiter)                              ::  mfi
TYPE(amrex_box)                                 ::  Box
TYPE(amrex_imultifab)                           ::  Level_Mask
INTEGER                                         ::  nComp
INTEGER, DIMENSION(3)                           ::  iEL, iEU
INTEGER,    CONTIGUOUS, POINTER                 ::  Mask_PTR(:,:,:,:)
REAL(idp),  CONTIGUOUS, POINTER                 ::  Results_PTR(:,:,:,:)

INTEGER                                         ::  re, te, pe
INTEGER                                         ::  rd, td, pd
INTEGER                                         ::  Num_Quad

INTEGER,    DIMENSION(3)                        ::  iEoff
REAL(idp), DIMENSION(3)                         ::  Gamma

CHARACTER(LEN = 100), DIMENSION(:), ALLOCATABLE             ::  Filenames
INTEGER, DIMENSION(:), ALLOCATABLE                          ::  File_IDs
INTEGER                                                     ::  Num_Files

REAL(idp), DIMENSION(1:NQ(1))                   ::  Cur_R_Locs
REAL(idp), DIMENSION(1:NQ(2))                   ::  Cur_T_Locs

REAL(idp)                                       ::  DROT
REAL(idp)                                       ::  DTOT

INTEGER                                         ::  FEM_Elem
INTEGER                                         ::  i

REAL(idp)                                       ::  Trace

116 FORMAT (A,A,A,A,A,A)

Num_Files = 1
ALLOCATE( Filenames(1:Num_Files) )
ALLOCATE( File_IDs(1:Num_Files) )

WRITE(Filenames(1),116) Poseidon_Results_Dir,"Results_Kij_",TRIM(File_Suffix),".out"

DO i = 1,Num_Files
    CALL OPEN_NEW_FILE( Filenames(i), File_IDs(i), 250 )
END DO


Num_Quad = NQ(1)*NQ(2)*NQ(3)


DO lvl = nLevels-1,0,-1

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

            IF ( amrex_spacedim == 1 ) THEN
                iEoff(2:3) = 0
            ELSEIF ( amrex_spacedim == 2) THEN
                iEoff(2)   = te
                iEoff(3)   = 0
            ELSEIF ( amrex_spacedim == 3 ) THEN
                iEoff(2) = te
                iEoff(3) = pe
            END IF


            FEM_Elem = FEM_Elem_Map(re,Lvl)

            DROT = drlocs(FEM_Elem)/2.0_idp
            DTOT = Level_dx(Lvl,2)/2.0_idp

            CUR_R_LOCS(:) = DROT * (Int_R_Locations(:) + 1.0_idp) + rlocs(FEM_Elem)
            CUR_T_LOCS(:) = DTOT * (Int_T_Locations(:) + 1.0_idp + iEOff(2)*2.0_idp)

            DO pd = 1,NQ(3)
            DO td = 1,NQ(2)
            DO rd = 1,NQ(1)

                Gamma(1) = Results_PTR(re,te,pe,AMReX_nCOMP_Map( iU_CF, rd, td, pd, NQ ))**4
                Gamma(2) = Gamma(1)/Cur_R_Locs(rd)**2
                Gamma(3) = Gamma(2)/DSIN(Cur_T_Locs(td))**2

                Trace = Results_PTR(re,te,pe,AMReX_nCOMP_Map( iU_K11, rd, td, pd, NQ ))/Gamma(1)    &
                      + Results_PTR(re,te,pe,AMReX_nCOMP_Map( iU_K11, rd, td, pd, NQ ))/Gamma(2)    &
                      + Results_PTR(re,te,pe,AMReX_nCOMP_Map( iU_K11, rd, td, pd, NQ ))/Gamma(3)
            
                PRINT*,Trace

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









END SUBROUTINE Write_Final_Results_Kij




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


END MODULE Driver_ReturnTest
