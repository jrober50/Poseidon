   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Initialization_AMReX                                                  !##!
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

USE Poseidon_Numbers_Module, &
            ONLY :  pi

USE Poseidon_Units_Module, &
            ONLY :  Set_Units

USE Poseidon_Parameters, &
            ONLY :  Domain_Dim,             &
                    Degree,                 &
                    L_Limit,                &
                    Method_Flag,            &
                    Verbose_Flag,           &
                    Convergence_Criteria,   &
                    Num_CFA_Vars,           &
                    Max_Iterations,         &
                    Poisson_Mode




USE Allocation_Core, &
            ONLY :  Allocate_Poseidon_CFA_Variables

USE Variables_MPI, &
            ONLY :  nProcs_Poseidon

USE Variables_Mesh, &
            ONLY :  Num_R_Elements,         &
                    Num_T_Elements,         &
                    Num_P_Elements,         &
                    R_Inner,                &
                    R_Outer

USE Variables_Quadrature, &
            ONLY :  Num_R_Quad_Points,      &
                    Num_T_Quad_Points,      &
                    Num_P_Quad_Points,      &
                    Num_TP_Quad_Points,     &
                    Num_Quad_DOF

USE Variables_Functions, &
            ONLY :  LM_Location

USE Initialization_XCFC, &
            ONLY :  Initialize_XCFC

USE Initialization_Quadrature, &
            ONLY :  Initialize_Quadrature

USE Initialization_Tables, &
            ONLY :  Initialize_Tables

USE Initialization_Derived, &
            ONLY :  Initialize_Derived_AMReX

USE Initialization_Subroutines, &
            ONLY :  Init_Fixed_Point_Params,        &
                    Init_IO_Params,                 &
                    Init_AMReX_Params


USE IO_Setup_Report_Module, &
            ONLY :  Output_Setup_Report

USE Timer_Routines_Module, &
            ONLY :  Init_Timers,                    &
                    Finalize_Timers,                &
                    TimerStart,                     &
                    TimerStop


USE Timer_Variables_Module, &
            ONLY :  Timer_Poisson_Matrix_Init,      &
                    Timer_Initialization_Core


USE Functions_Mapping, &
            ONLY :  CFA_3D_LM_Map




USE Variables_AMReX_Source, &
            ONLY :  Source_PTR,             &
                    Mask_PTR,               &
                    iCoarse,                &
                    iFine

USE Poseidon_AMReX_MakeFineMask_Module, &
            ONLY :  AMReX_MakeFineMask

USE Poseidon_AMReX_BoxArraySize_Module, &
            ONLY :  AMReX_BoxArraySize

USE Poseidon_AMReX_Multilayer_Utilities_Module, &
            ONLY :  Find_Coarsest_Parent

USE Variables_MPI, &
            ONLY :  myID_Poseidon



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
                    AMReX_Mode,             &
                    iNumLeafElements,       &
                    iLeafElementsPerLvl,    &
                    Findloc_Table,          &
                    FEM_Elem_Table,         &
                    Table_Offsets


IMPLICIT NONE


CONTAINS



 !+101+############################################################!
!                                                                   !
!          Initialize_Poseidon_with_AMReX                           !
!                                                                   !
 !#################################################################!
SUBROUTINE Initialize_Poseidon_with_AMReX(  FEM_Degree_Option,                  &
                                            L_Limit_Option,                     &
                                            Units_Option,                       &
                                            Domain_Edge_Option,                 &
                                            Coarse_NE_Option,                   &
                                            NQ_Option,                          &
                                            Max_Iterations_Option,              &
                                            Convergence_Criteria_Option,        &
                                            Anderson_M_Option,                  &
                                            CFA_Eq_Flags_Option,                &
                                            AMReX_Max_Levels_Option,            &
                                            AMReX_Max_Grid_Size_Option,         &
                                            AMReX_FEM_Refinement_Option,        &
                                            AMReX_Integral_Refinement_Option,   &
                                            Poisson_Mode_Option,                &
                                            Verbose_Option,                     &
                                            WriteAll_Option,                    &
                                            Print_Setup_Option,                 &
                                            Write_Setup_Option,                 &
                                            Print_Results_Option,               &
                                            Write_Results_Option,               &
                                            Print_Timetable_Option,             &
                                            Write_Timetable_Option,             &
                                            Write_Sources_Option,               &
                                            Suffix_Flag_Option,                 &
                                            Suffix_Tail_Option,                 &
                                            Frame_Option                        )

CHARACTER(LEN=1),        INTENT(IN), OPTIONAL               ::  Units_Option

INTEGER,                 INTENT(IN), OPTIONAL               ::  FEM_Degree_Option
INTEGER,                 INTENT(IN), OPTIONAL               ::  L_Limit_Option

REAL(idp), DIMENSION(2), INTENT(IN), OPTIONAL               ::  Domain_Edge_Option
INTEGER,   DIMENSION(3), INTENT(IN), OPTIONAL               ::  Coarse_NE_Option
INTEGER,   DIMENSION(3), INTENT(IN), OPTIONAL               ::  NQ_Option

INTEGER,                 INTENT(IN), OPTIONAL               ::  Max_Iterations_Option
REAL(idp),               INTENT(IN), OPTIONAL               ::  Convergence_Criteria_Option
INTEGER,                 INTENT(IN), OPTIONAL               ::  Anderson_M_Option
INTEGER,   DIMENSION(5), INTENT(IN), OPTIONAL               ::  CFA_EQ_Flags_Option


INTEGER,                 INTENT(IN), OPTIONAL               ::  AMReX_Max_Levels_Option
INTEGER,   DIMENSION(3), INTENT(IN), OPTIONAL               ::  AMReX_Max_Grid_Size_Option
INTEGER,                 INTENT(IN), OPTIONAL               ::  AMReX_FEM_Refinement_Option
INTEGER,                 INTENT(IN), OPTIONAL               ::  AMReX_Integral_Refinement_Option

LOGICAL,                 INTENT(IN), OPTIONAL               ::  Poisson_Mode_Option

LOGICAL,                 INTENT(IN), OPTIONAL               ::  Verbose_Option
LOGICAL,                 INTENT(IN), OPTIONAL               ::  WriteAll_Option
LOGICAL,                 INTENT(IN), OPTIONAL               ::  Print_Setup_Option
LOGICAL,                 INTENT(IN), OPTIONAL               ::  Write_Setup_Option
LOGICAL,                 INTENT(IN), OPTIONAL               ::  Print_Results_Option
LOGICAL,                 INTENT(IN), OPTIONAL               ::  Write_Results_Option
LOGICAL,                 INTENT(IN), OPTIONAL               ::  Print_Timetable_Option
LOGICAL,                 INTENT(IN), OPTIONAL               ::  Write_Timetable_Option
LOGICAL,                 INTENT(IN), OPTIONAL               ::  Write_Sources_Option

CHARACTER(LEN=10),       INTENT(IN), OPTIONAL               ::  Suffix_Flag_Option
CHARACTER(LEN=1),        INTENT(IN), OPTIONAL               ::  Suffix_Tail_Option
INTEGER,                 INTENT(IN), OPTIONAL               ::  Frame_Option


CALL Init_Timers
CALL TimerStart( Timer_Initialization_Core )





IF ( PRESENT( Verbose_Option ) ) THEN
    Verbose_Flag = Verbose_Option
ELSE
    Verbose_Flag = .FALSE.
END IF

IF ( Verbose_Flag .EQV. .TRUE. ) THEN
    PRINT*,"Initializing Poseidon with AMReX..."
END IF

AMReX_Mode = .TRUE.


IF ( PRESENT( Units_Option ) ) THEN
    CALL Set_Units(Units_Option)
ELSE
    CALL Set_Units("G")
END IF


IF ( PRESENT( FEM_Degree_Option ) ) THEN
    Degree = FEM_Degree_Option
ELSE
    Degree = 1
END IF


IF ( PRESENT( L_Limit_Option ) ) THEN
    L_Limit = L_Limit_Option
ELSE
    L_Limit = 0
END IF



IF ( .FALSE. ) THEN
!    Domain_Dim = Dimensions_Option
ELSE
    Domain_Dim = 3
END IF


IF ( PRESENT( Domain_Edge_Option ) ) THEN
    R_Inner = Domain_Edge_Option(1)
    R_Outer = Domain_Edge_Option(2)
ELSE
    R_Inner = 0.0_idp
    R_Outer = 1.0_idp
END IF


IF ( PRESENT( Coarse_NE_Option ) ) THEN
    Num_R_Elements = Coarse_NE_Option(1)
    Num_T_Elements = Coarse_NE_Option(2)
    Num_P_Elements = Coarse_NE_Option(3)
ELSE
    Num_R_Elements = 1
    Num_T_Elements = 1
    Num_P_Elements = 1
END IF



IF ( PRESENT( NQ_Option ) ) THEN
    Num_R_Quad_Points = NQ_Option(1)
    Num_T_Quad_Points = NQ_Option(2)
    Num_P_Quad_Points = NQ_Option(3)
ELSE
    Num_R_Quad_Points = 10
    Num_T_Quad_Points = 20
    Num_P_Quad_Points = 2*L_Limit + 1
END IF
Num_TP_Quad_Points = Num_T_Quad_Points*Num_P_Quad_Points
Num_Quad_DOF       = Num_R_Quad_Points*Num_TP_Quad_Points





CALL Init_IO_Params(WriteAll_Option,            &
                    Print_Setup_Option,         &
                    Write_Setup_Option,         &
                    Print_Results_Option,       &
                    Write_Results_Option,       &
                    Print_Timetable_Option,     &
                    Write_Timetable_Option,     &
                    Write_Sources_Option,       &
                    Suffix_Flag_Option,         &
                    Suffix_Tail_Option,         &
                    Frame_Option                )




CALL Init_AMReX_Params( AMReX_Max_Levels_Option,            &
                        AMReX_Max_Grid_Size_Option,         &
                        AMReX_FEM_Refinement_Option,        &
                        AMReX_Integral_Refinement_Option    )



IF ( PRESENT(Poisson_Mode_Option) ) THEN
    Poisson_Mode = Poisson_Mode_Option
ELSE
    Poisson_Mode = .FALSE.
END IF


IF ( Poisson_Mode ) THEN
    !=======================================================!
    !                                                       !
    !               Initialize Poisson Solver               !
    !                                                       !
    !=======================================================!


ELSE
    !=======================================================!
    !                                                       !
    !           Initialize CFA/XCFC Metric Solver           !
    !                                                       !
    !=======================================================!
    nProcs_Poseidon = 1
    Method_Flag = 3


    CALL Init_Fixed_Point_Params( Max_Iterations_Option,          &
                                  Convergence_Criteria_Option,    &
                                  Anderson_M_Option               )


    LM_Location => CFA_3D_LM_Map


    CALL Initialize_Derived_AMReX()
    CALL Allocate_Poseidon_CFA_Variables()
    CALL Initialize_Quadrature()
    CALL Initialize_Tables()


!    CALL Initialize_XCFC_with_AMReX(CFA_EQ_Flags_Option)



END IF ! Not Poisson Mode




IF ( Verbose_Flag ) THEN
    PRINT*,"Poseidon Initialization Complete"
END IF

IF ( Verbose_Flag ) THEN
    CALL Output_Setup_Report()
END IF

CALL TimerStop( Timer_Initialization_Core )




END SUBROUTINE Initialize_Poseidon_with_AMReX















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


END SUBROUTINE Initialize_AMReX_Maps






 !+101+############################################################!
!                                                                   !
!          Allocate_AMReX_Maps                                      !
!                                                                   !
 !#################################################################!
SUBROUTINE Allocate_AMReX_Maps()



IF ( .NOT. ALLOCATED(FindLoc_Table) ) THEN
    ALLOCATE( FindLoc_Table(0:iNumLeafElements-1) )
ELSE
    DEALLOCATE( FindLoc_Table )
    ALLOCATE( FindLoc_Table(0:iNumLeafElements-1) )
END IF


IF ( .NOT. ALLOCATED(FEM_Elem_Table) ) THEN
    ALLOCATE( FEM_Elem_Table(0:iNumLeafElements-1) )
ELSE
    DEALLOCATE( FEM_Elem_Table )
    ALLOCATE( FEM_Elem_Table(0:iNumLeafElements-1) )
END IF



IF ( .NOT. ALLOCATED(Table_Offsets) ) THEN
    ALLOCATE( Table_Offsets(0:AMReX_Num_Levels) )
ELSE
    DEALLOCATE( Table_Offsets )
    ALLOCATE( Table_Offsets(0:AMReX_Num_Levels) )
END IF


END SUBROUTINE Allocate_AMReX_Maps






 !+101+############################################################!
!                                                                   !
!          Deallocate_AMReX_Maps                                    !
!                                                                   !
 !#################################################################!
SUBROUTINE Deallocate_AMReX_Maps()

DEALLOCATE( FindLoc_Table )
DEALLOCATE( FEM_Elem_Table )

END SUBROUTINE Deallocate_AMReX_Maps






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


IF ( amrex_spacedim == 1 ) THEN
    Eoff(2:3) = 1
ELSEIF ( amrex_spacedim == 2) THEN
    Eoff(2)   = 0
    Eoff(3)   = 1
ELSEIF ( amrex_spacedim == 3 ) THEN
    Eoff(2:3) = 0
END IF

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
TYPE(amrex_boxarray)                            ::  BA
TYPE(amrex_box)                                 ::  Box
TYPE(amrex_distromap)                           ::  DM
TYPE(amrex_imultifab)                           ::  Level_Mask
INTEGER, ALLOCATABLE                            ::  mypmap(:)

INTEGER                                         ::  Here
INTEGER                                         ::  Offset
INTEGER                                         ::  RE, TE, PE
INTEGER, DIMENSION(3)                           ::  ELo, EHi

INTEGER, CONTIGUOUS, POINTER                    :: Mask(:,:,:,:)

INTEGER                                         :: iLeaf = 1
INTEGER                                         :: iTrunk   = 0

INTEGER, DIMENSION(3)                           ::  EOff


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

!PRINT*,"FindLoc_Table"
!PRINT*,FindLoc_Table

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
INTEGER                                         ::  Index
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

!PRINT*,"FEM_Elem_Table"
!PRINT*,FEM_Elem_Table


#endif
END SUBROUTINE Create_FEM_Elem_Table












END MODULE Initialization_AMReX
