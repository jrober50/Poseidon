   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Maps_AMReX                                                            !##!
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

USE Variables_MPI, &
            ONLY :  myID_Poseidon

USE Timer_Routines_Module, &
            ONLY :  TimerStart,             &
                    TimerStop

USE Variables_AMReX_Source, &
            ONLY :  iTrunk,                 &
                    iLeaf

USE Poseidon_AMReX_MakeFineMask_Module,     &
            ONLY :  AMReX_MakeFineMask

USE Poseidon_AMReX_BoxArraySize_Module,     &
            ONLY :  AMReX_BoxArraySize

#ifdef POSEIDON_AMREX_FLAG
use amrex_base_module

USE amrex_fort_module, &
            ONLY :  amrex_spacedim

USE amrex_box_module,   &
            ONLY:   amrex_box

USE amrex_boxarray_module, &
            ONLY:   amrex_boxarray

USE amrex_distromap_module, &
            ONLY:   amrex_distromap,        &
                    amrex_distromap_build,  &
                    amrex_distromap_destroy

USE amrex_multifab_module,  &
            ONLY:   amrex_multifab,         &
                    amrex_multifab_build,   &
                    amrex_imultifab_build,  &
                    amrex_imultifab_destroy

USE Variables_AMReX_Core, &
            ONLY :  MF_Source
#endif

USE Variables_AMReX_Core, &
            ONLY :  AMReX_Num_Levels,       &
                    iNumLeafElements,       &
                    iLeafElementsPerLvl,    &
                    Findloc_Table,          &
                    FEM_Elem_Table,         &
                    Table_Offsets

USE Flags_Initial_Guess_Module, &
            ONLY :  lPF_IG_Flags,           &
                    iPF_IG_Flat_Guess

USE Flags_Initialization_Module, &
            ONLY :  lPF_Init_MTGV_Flags,    &
                    iPF_Init_MTGV_TransMat, &
                    lPF_Init_AMReX_Flags,   &
                    iPF_Init_AMReX_Maps


IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!                                                     				!
!                                                                               !
!###############################################################################!
PURE INTEGER FUNCTION AMReX_nCOMP_Map( iU, rd, td, pd, NQ )

INTEGER,                                        INTENT(IN)  ::  iU
INTEGER,                                        INTENT(IN)  ::  rd
INTEGER,                                        INTENT(IN)  ::  td
INTEGER,                                        INTENT(IN)  ::  pd
INTEGER,    DIMENSION(3),                       INTENT(IN)  ::  NQ


AMReX_nCOMP_Map = (iU-1)*NQ(1)*NQ(2)*NQ(3)          &
                + (pd-1)*NQ(1)*NQ(2)                &
                + (td-1)*NQ(1)                      &
                +  rd

END FUNCTION AMReX_nCOMP_Map









 !+101+############################################################!
!                                                                   !
!          Initialize_AMReX_Maps                                    !
!                                                                   !
 !#################################################################!
SUBROUTINE Initialize_AMReX_Maps()

INTEGER                         :: lvl



CALL Count_The_Leafs()

CALL Allocate_AMReX_Maps()

Table_Offsets(0) = 0
DO lvl = 1,AMReX_Num_Levels
    Table_Offsets(lvl) = SUM(iLeafElementsPerLvl(0:lvl-1))
END DO

CALL Create_FindLoc_Table()

CALL Create_FEM_Elem_Table()


lPF_Init_AMReX_Flags(iPF_Init_AMReX_Maps) = .TRUE.

END SUBROUTINE Initialize_AMReX_Maps



 !+102+############################################################!
!                                                                   !
!          Reinitialize_AMReX_Maps                                  !
!                                                                   !
 !#################################################################!
SUBROUTINE Reinitialize_AMReX_Maps()

INTEGER                         :: lvl

lPF_Init_AMReX_Flags(iPF_Init_AMReX_Maps) = .FALSE.

iNumLeafElements = 0
iLeafElementsPerLvl = 0

CALL Deallocate_AMReX_Maps()
CALL Initialize_AMReX_Maps()

lPF_Init_AMReX_Flags(iPF_Init_AMReX_Maps) = .TRUE.

END SUBROUTINE Reinitialize_AMReX_Maps





 !+101+############################################################!
!                                                                   !
!          Allocate_AMReX_Maps                                      !
!                                                                   !
 !#################################################################!
SUBROUTINE Allocate_AMReX_Maps()

ALLOCATE( FindLoc_Table(0:iNumLeafElements-1)  )
ALLOCATE( FEM_Elem_Table(0:iNumLeafElements-1) )
ALLOCATE( Table_Offsets(0:AMReX_Num_Levels)    )

END SUBROUTINE Allocate_AMReX_Maps






 !+101+############################################################!
!                                                                   !
!          Deallocate_AMReX_Maps                                    !
!                                                                   !
 !#################################################################!
SUBROUTINE Deallocate_AMReX_Maps()

DEALLOCATE( FindLoc_Table  )
DEALLOCATE( FEM_Elem_Table )
DEALLOCATE( Table_Offsets  )

END SUBROUTINE Deallocate_AMReX_Maps




 !+101+############################################################!
!                                                                   !
!          Reallocate_AMReX_Maps                                      !
!                                                                   !
 !#################################################################!
SUBROUTINE Reallocate_AMReX_Maps()

CALL Deallocate_AMReX_Maps()
CALL Allocate_AMReX_Maps

END SUBROUTINE Reallocate_AMReX_Maps




 !+101+############################################################!
!                                                                   !
!          Count_The_Leafs                                          !
!                                                                   !
 !#################################################################!
SUBROUTINE Count_The_Leafs()
#ifdef POSEIDON_AMREX_FLAG
INTEGER                                         ::  lvl
INTEGER                                         ::  numBoxes
INTEGER, ALLOCATABLE                            ::  mypmap(:)
INTEGER, CONTIGUOUS, POINTER                    ::  Mask(:,:,:,:)

TYPE(amrex_imultifab)                           ::  Level_Mask
TYPE(amrex_distromap)                           ::  DM
TYPE(amrex_mfiter)                              ::  mfi
TYPE(amrex_box)                                 ::  Box

INTEGER, DIMENSION(3)                           ::  ELo, EHi
INTEGER                                         ::  iLeaf  = 1
INTEGER                                         ::  iTrunk = 0

INTEGER, DIMENSION(3)                           ::  EOff
INTEGER, DIMENSION(1:3)                         ::  nGhost_Vec

IF ( amrex_spacedim == 1 ) THEN
    Eoff(2:3) = 1
ELSEIF ( amrex_spacedim == 2) THEN
    Eoff(2)   = 0
    Eoff(3)   = 1
ELSEIF ( amrex_spacedim == 3 ) THEN
    Eoff(2:3) = 0
END IF

nGhost_Vec = 0

iLeafElementsPerLvl = 0


DO lvl = 0,AMReX_Num_Levels-1

    ! Because AMReX will only store data owned by a process,
    ! using the Distrobution Mapping (DM) from MF_Source will
    ! only allow us to see the portion of the FineMask owned by a
    ! process. But we wish to view the whole domain.
    ! Therefore we must construct a DM that says the local process
    ! owns the entire domain's worth of data.
    ! mypmap is an array of integers with the id of the local process.
    numBoxes = AMReX_BoxArraySize(MF_Source(lvl)%BA)
    ALLOCATE(mypmap(0:numBoxes-1) )
    mypmap = myID_Poseidon
    CALL amrex_distromap_build(DM,mypmap)


    IF ( lvl < AMReX_Num_Levels-1 ) THEN
        CALL AMReX_MakeFineMask(  Level_Mask,               &
                                  MF_Source(lvl)%ba,        &
                                  DM,                       &
                                  nGhost_Vec,               &
                                  MF_Source(lvl+1)%ba,      &
                                  iLeaf, iTrunk  )
    ELSE
        ! Create Level_Mask all equal to 1
        CALL amrex_imultifab_build( Level_Mask,             &
                                    MF_Source(lvl)%ba,      &
                                    DM,      &
                                    1,                      &
                                    0                       )
        CALL Level_Mask%SetVal(iLeaf)
    END IF


    CALL amrex_mfiter_build(mfi, Level_Mask, tiling = .false. )
    DO WHILE(mfi%next())

        Mask => Level_Mask%dataPtr(mfi)
        Box = mfi%tilebox()

        ELo = Box%lo
        EHi = Box%hi

        IF (( ELo(2) == EOff(2)) .AND. (Elo(3) == EOff(3) )) THEN
            iLeafElementsPerLvl(lvl) = iLeafElementsPerLvl(lvl) + SUM(Mask(:,EOff(2),EOff(3),:))
        END IF
    END DO
    iNumLeafElements = Sum(iLeafElementsPerLvl)

    CALL amrex_mfiter_destroy(mfi)
    CALL amrex_imultifab_destroy(Level_Mask)
    CALL amrex_distromap_destroy(DM)
    DEALLOCATE( mypmap )
END DO


#endif
END SUBROUTINE Count_The_Leafs






 !+101+############################################################!
!                                                                   !
!          Create_Findloc_Table                                     !
!                                                                   !
 !#################################################################!
SUBROUTINE Create_Findloc_Table()
#ifdef POSEIDON_AMREX_FLAG

INTEGER                         :: lvl
INTEGER                         :: numBoxes

TYPE(amrex_mfiter)                              ::  mfi
TYPE(amrex_box)                                 ::  Box
TYPE(amrex_distromap)                           ::  DM
TYPE(amrex_imultifab)                           ::  Level_Mask
INTEGER, ALLOCATABLE                            ::  mypmap(:)

INTEGER                                         ::  Here
INTEGER                                         ::  Offset
INTEGER                                         ::  RE
INTEGER, DIMENSION(3)                           ::  ELo, EHi

INTEGER, CONTIGUOUS, POINTER                    :: Mask(:,:,:,:)


INTEGER, DIMENSION(3)                           ::  EOff
INTEGER, DIMENSION(1:3)                         ::  nGhost_Vec

nGhost_Vec = 0

IF ( amrex_spacedim == 1 ) THEN
    Eoff(2:3) = 1
ELSEIF ( amrex_spacedim == 2) THEN
    Eoff(2)   = 0
    Eoff(3)   = 1
ELSEIF ( amrex_spacedim == 3 ) THEN
    Eoff(2:3) = 0
END IF



DO lvl = 0,AMReX_Num_Levels-1

    ! Because AMReX will only store data owned by a process,
    ! using the Distrobution Mapping (DM) from MF_Source will
    ! only allow us to see the portion of the FineMask owned by a
    ! process. But we wish to view the whole domain.
    ! Therefore we must construct a DM that says the local process
    ! owns the entire domain's worth of data.
    ! mypmap is an array of integers with the id of the local process.
    numBoxes = AMReX_BoxArraySize(MF_Source(lvl)%BA)
    ALLOCATE(mypmap(0:numBoxes-1) )
    mypmap = myID_Poseidon
    CALL amrex_distromap_build(DM,mypmap)


    IF ( lvl < AMReX_Num_Levels-1 ) THEN
        CALL AMReX_MakeFineMask(  Level_Mask,               &
                                  MF_Source(lvl)%ba,        &
                                  DM,                       &
                                  nGhost_Vec,               &
                                  MF_Source(lvl+1)%ba,      &
                                  iLeaf, iTrunk  )
    ELSE
        ! Create Level_Mask all equal to 1
        CALL amrex_imultifab_build( Level_Mask,             &
                                    MF_Source(lvl)%ba,      &
                                    DM,      &
                                    1,                      &
                                    0                       )
        CALL Level_Mask%SetVal(iLeaf)
    END IF




    Offset = 0
    CALL amrex_mfiter_build(mfi, Level_Mask, tiling = .false. )
    DO WHILE(mfi%next())

        Mask => Level_Mask%dataPtr(mfi)

        Box = mfi%tilebox()

        ELo = Box%lo
        EHi = Box%hi

        IF (( ELo(2) == EOff(2)) .AND. (Elo(3) == EOff(3) )) THEN
        DO RE = ELo(1), EHi(1)
            
            Here = SUM(iLeafElementsPerLvl(0:lvl-1) ) + Offset

            IF ( SUM(Mask(RE,EOff(2),EOff(3),:)) == iLeaf ) THEN
                FindLoc_Table(Here) = RE
                Offset = Offset + 1
            END IF

        END DO
        END IF


    END DO

    CALL amrex_mfiter_destroy(mfi)
    CALL amrex_imultifab_destroy(Level_Mask)
    CALL amrex_distromap_destroy(DM)
    DEALLOCATE( mypmap )

END DO


#endif
END SUBROUTINE Create_Findloc_Table





 !+101+############################################################!
!                                                                   !
!          Create_FEM_Elem_Table                                    !
!                                                                   !
 !#################################################################!
SUBROUTINE Create_FEM_Elem_Table()
#ifdef POSEIDON_AMREX_FLAG


INTEGER, DIMENSION(0:iNumLeafElements-1)        ::  Tmp_Array
INTEGER                                         ::  lvl
INTEGER                                         ::  elem
INTEGER                                         ::  Here
INTEGER                                         ::  There

Tmp_Array = FindLoc_Table
FEM_Elem_Table = -1

! Modify Tmp_Arry.
! Multiply each row by 2^(AMReX_Max_Level - Cur_Level).
! This maps each element to its 'left-most element' on the finest mesh.
! Example: Element 2 on level 0, would contain elements 4 and 5 on level 1.
!           This mapping would return 4.
! This puts all elements "on the same level" and reveals their order.
DO lvl = 0,AMReX_Num_Levels-1
    Here  = Table_Offsets(lvl)
    There = Table_Offsets(lvl+1)-1
    Tmp_Array(Here:There) = Tmp_Array(Here:There)*2**(AMReX_Num_Levels-lvl-1)
END DO




! Fill in FEM_Elem_Table.
! No we search through Tmp_Array for the location with the largest value.
! The location represents the location in Findloc_Table that is associated with
! the 'right-most' element, i.e. the last element. That location in FEM_Elem_Table
! is then given the value iNumLeafElements-1.  That location in Tmp_Array is then set
! to -1.  The process repeats to find the element to the left of the first, and
! gives its location in FEM_Elem_Table iNumLeafElements-2. This process repeats
! until all values in FEM_Elem_Table are full.
DO elem = iNumLeafElements-1,0,-1
    Here = maxloc(Tmp_Array(:),dim=1)-1
    FEM_Elem_Table(Here) = elem             ! Since the arrays start their indexing at 0,
    Tmp_Array(Here)      = -1               ! and the Fortran standard is to start at 1,
                                            ! Here-1 must be used as the array index.
END DO



#endif
END SUBROUTINE Create_FEM_Elem_Table






END MODULE Maps_AMReX
