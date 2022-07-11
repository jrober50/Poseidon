   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Poseidon_Interface_Initialization                                      !##!
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

USE Poseidon_Message_Routines_Module, &
            ONLY :  Init_Message

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
                    CFA_EQ_Flags,           &
                    Num_CFA_Eqs

USE Allocation_Sources, &
            ONLY :  Allocate_Poseidon_Source_Variables

USE XCFC_Source_Routine_Variables_Module, &
            ONLY :  Allocate_XCFC_Source_Routine_Variables

USE Variables_MPI, &
            ONLY :  nProcs_Poseidon,        &
                    myID_Poseidon,          &
                    ierr

USE Variables_Quadrature, &
            ONLY :  Num_R_Quad_Points,      &
                    Num_T_Quad_Points,      &
                    Num_P_Quad_Points,      &
                    Num_TP_Quad_Points,     &
                    Local_Quad_DOF,         &
                    Int_R_Locations,        &
                    Int_T_Locations,        &
                    Int_P_Locations,        &
                    Int_R_Weights,          &
                    Int_T_Weights,          &
                    Int_P_Weights,          &
                    Int_TP_Weights,         &
                    xLeftLimit,            &
                    xRightLimit

USE Initialization_XCFC, &
            ONLY :  Initialize_XCFC

USE Initialization_Poisson, &
            ONLY :  Initialize_Poisson

USE Initialization_Quadrature, &
            ONLY :  Initialize_Quadrature

USE Initialization_Tables, &
            ONLY :  Initialize_Tables

USE Initialization_Derived, &
            ONLY :  Initialize_Derived,         &
                    Initialize_Derived_AMReX

USE Initialization_Subroutines, &
            ONLY :  Init_Fixed_Point_Params,        &
                    Init_IO_Params,                 &
                    Init_MPI_Params,                &
                    Init_Quad_Params,               &
                    Init_Expansion_Params,          &
                    Init_Mesh_Params,               &
                    Initialize_Mesh,                &
                    Init_AMReX_Params,              &
                    Set_Caller_Data


USE Initialization_Subroutines_AMReX, &
            ONLY : Init_Parameters_From_AMReX_Input_File


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


USE Variables_AMReX_Source, &
            ONLY :  Source_PTR,             &
                    Mask_PTR,               &
                    iTrunk,                 &
                    iLeaf

USE Variables_Interface, &
            ONLY :  Caller_R_Units

USE Poseidon_AMReX_MakeFineMask_Module, &
            ONLY :  AMReX_MakeFineMask

USE Poseidon_AMReX_BoxArraySize_Module, &
            ONLY :  AMReX_BoxArraySize

USE Poseidon_AMReX_Multilayer_Utilities_Module, &
            ONLY :  Find_Coarsest_Parent


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

USE Functions_Translation_Matrix_Module, &
            ONLY :  Create_Translation_Matrix


USE Variables_Interface, &
            ONLY :  Caller_NQ,                  &
                    Caller_Quad_DOF,            &
                    Caller_xL,                  &
                    Caller_RQ_xlocs,            &
                    Caller_TQ_xlocs,            &
                    Caller_PQ_xlocs,            &
                    Caller_R_Units,             &
                    Translation_Matrix

USE Variables_AMReX_Core, &
            ONLY :  AMReX_Num_Levels,       &
                    iNumLeafElements,       &
                    iLeafElementsPerLvl,    &
                    Findloc_Table,          &
                    FEM_Elem_Table,         &
                    Table_Offsets

USE Flags_Core_Module, &
            ONLY :  lPF_Core_Flags,         &
                    iPF_Core_Poisson_Mode,  &
                    iPF_Core_CFA_Mode,      &
                    iPF_Core_XCFC_Mode,     &
                    iPF_Core_AMReX_Mode

USE Flags_Initial_Guess_Module, &
            ONLY :  lPF_IG_Flags,           &
                    iPF_IG_Flat_Guess

USE Flags_Initialization_Module, &
            ONLY :  lPF_Init_MTGV_Flags,    &
                    iPF_Init_MTGV_TransMat, &
                    lPF_Init_AMReX_Flags,   &
                    iPF_Init_AMReX_Maps

USE Flags_IO_Module, &
            ONLY :  lPF_IO_Flags,               &
                    iPF_IO_Write_Setup


IMPLICIT NONE


CONTAINS




 !+101+############################################################!
!                                                                   !
!          Initialize_Poseidon_with_AMReX                           !
!                                                                   !
 !#################################################################!
SUBROUTINE Initialize_Poseidon( Source_NE,                          &
                                Source_NQ,                          &
                                Source_xL,                          &
                                Source_RQ_xlocs,                    &
                                Source_TQ_xlocs,                    &
                                Source_PQ_xlocs,                    &
                                Source_Units,                       &
                                Source_Radial_Boundary_Units,       &
                                Domain_Edge_Option,                     &
                                Dimensions_Option,                      &
                                Method_Flag_Option,                     &
                                FEM_Degree_Option,                      &
                                L_Limit_Option,                         &
                                Source_R_Option,                        &
                                Source_T_Option,                        &
                                Source_P_Option,                        &
                                Source_DR_Option,                       &
                                Source_DT_Option,                       &
                                Source_DP_Option,                       &
                                Integration_NQ_Option,              &
                                Max_Iterations_Option,                  &
                                Convergence_Criteria_Option,        &
                                Anderson_M_Option,                      &
                                CFA_Eq_Flags_Option,                &
                                AMReX_FEM_Refinement_Option,        &
                                AMReX_Integral_Refinement_Option,   &
                                Poisson_Mode_Option,                &
                                Flat_Guess_Option,                  &
                                Verbose_Option,                     &
                                WriteAll_Option,                    &
                                Print_Setup_Option,                 &
                                Write_Setup_Option,                 &
                                Print_Results_Option,               &
                                Write_Results_Option,               &
                                Print_Timetable_Option,             &
                                Write_Timetable_Option,             &
                                Write_Sources_Option,               &
                                Print_Condition_Option,             &
                                Write_Condition_Option,             &
                                Suffix_Flag_Option,                 &
                                Suffix_Tail_Option,                 &
                                Frame_Option                        )



INTEGER,    DIMENSION(3),               INTENT(IN)              ::  Source_NQ
REAL(idp),  DIMENSION(2),               INTENT(IN)              ::  Source_xL
REAL(idp),  DIMENSION(Source_NQ(1)),    INTENT(IN)              ::  Source_RQ_xlocs
REAL(idp),  DIMENSION(Source_NQ(2)),    INTENT(IN)              ::  Source_TQ_xlocs
REAL(idp),  DIMENSION(Source_NQ(3)),    INTENT(IN)              ::  Source_PQ_xlocs
CHARACTER(LEN=1),                       INTENT(IN)              ::  Source_Units
CHARACTER(LEN=2),                       INTENT(IN)              ::  Source_Radial_Boundary_Units
REAL(idp), DIMENSION(:), INTENT(IN), OPTIONAL               ::  Source_R_Option
REAL(idp), DIMENSION(:), INTENT(IN), OPTIONAL               ::  Source_T_Option
REAL(idp), DIMENSION(:), INTENT(IN), OPTIONAL               ::  Source_P_Option
REAL(idp), DIMENSION(:), INTENT(IN), OPTIONAL               ::  Source_DR_Option
REAL(idp), DIMENSION(:), INTENT(IN), OPTIONAL               ::  Source_DT_Option
REAL(idp), DIMENSION(:), INTENT(IN), OPTIONAL               ::  Source_DP_Option


INTEGER,                 INTENT(IN), OPTIONAL               ::  Dimensions_Option
INTEGER,                 INTENT(IN), OPTIONAL               ::  FEM_Degree_Option
INTEGER,                 INTENT(IN), OPTIONAL               ::  L_Limit_Option

INTEGER,   DIMENSION(3), INTENT(IN), OPTIONAL               ::  Source_NE
REAL(idp), DIMENSION(2), INTENT(IN), OPTIONAL               ::  Domain_Edge_Option

INTEGER,   DIMENSION(3), INTENT(IN), OPTIONAL                   ::  Integration_NQ_Option
INTEGER,   DIMENSION(5), INTENT(IN), OPTIONAL                   ::  CFA_EQ_Flags_Option

INTEGER,                 INTENT(IN), OPTIONAL               ::  AMReX_FEM_Refinement_Option
INTEGER,                 INTENT(IN), OPTIONAL               ::  AMReX_Integral_Refinement_Option

LOGICAL,                 INTENT(IN), OPTIONAL               ::  Poisson_Mode_Option
LOGICAL,                 INTENT(IN), OPTIONAL               ::  Flat_Guess_Option

LOGICAL,                 INTENT(IN), OPTIONAL               ::  Verbose_Option
LOGICAL,                 INTENT(IN), OPTIONAL               ::  WriteAll_Option
LOGICAL,                 INTENT(IN), OPTIONAL               ::  Print_Setup_Option
LOGICAL,                 INTENT(IN), OPTIONAL               ::  Write_Setup_Option
LOGICAL,                 INTENT(IN), OPTIONAL               ::  Print_Results_Option
LOGICAL,                 INTENT(IN), OPTIONAL               ::  Write_Results_Option
LOGICAL,                 INTENT(IN), OPTIONAL               ::  Print_Timetable_Option
LOGICAL,                 INTENT(IN), OPTIONAL               ::  Write_Timetable_Option
LOGICAL,                 INTENT(IN), OPTIONAL               ::  Write_Sources_Option
LOGICAL,                 INTENT(IN), OPTIONAL               ::  Print_Condition_Option
LOGICAL,                 INTENT(IN), OPTIONAL               ::  Write_Condition_Option

CHARACTER(LEN=10),       INTENT(IN), OPTIONAL               ::  Suffix_Flag_Option
CHARACTER(LEN=1),        INTENT(IN), OPTIONAL               ::  Suffix_Tail_Option
INTEGER,                 INTENT(IN), OPTIONAL               ::  Frame_Option

INTEGER,                 INTENT(IN), OPTIONAL               ::  Max_Iterations_Option
REAL(idp),               INTENT(IN), OPTIONAL               ::  Convergence_Criteria_Option
INTEGER,                 INTENT(IN), OPTIONAL               ::  Anderson_M_Option

INTEGER,                 INTENT(IN), OPTIONAL               ::  Method_Flag_Option


CALL Init_Timers
CALL TimerStart( Timer_Initialization_Core )


IF ( PRESENT( Verbose_Option ) ) THEN
    Verbose_Flag = Verbose_Option
ELSE
    Verbose_Flag = .FALSE.
END IF
IF ( Verbose_Flag ) CALL Init_Message('Beginning Poseidon Core Initialization.')


CALL Set_Units(Source_Units)



IF ( PRESENT(Poisson_Mode_Option) ) THEN
    IF ( Poisson_Mode_Option ) THEN
        lPF_Core_Flags(iPF_Core_Poisson_Mode)   = .TRUE.
        lPF_Core_Flags(iPF_Core_CFA_Mode)       = .FALSE.
        lPF_Core_Flags(iPF_Core_XCFC_Mode)      = .FALSE.
    ELSE
        lPF_Core_Flags(iPF_Core_Poisson_Mode)   = .FALSE.
        lPF_Core_Flags(iPF_Core_CFA_Mode)       = .FALSE.
        lPF_Core_Flags(iPF_Core_XCFC_Mode)      = .TRUE.
    END IF
ELSE
    lPF_Core_Flags(iPF_Core_Poisson_Mode)   = .FALSE.
    lPF_Core_Flags(iPF_Core_CFA_Mode)       = .FALSE.
    lPF_Core_Flags(iPF_Core_XCFC_Mode)      = .TRUE.
END IF

IF ( PRESENT( Method_Flag_Option ) ) THEN
    Method_Flag = Method_Flag_Option
ELSE
    Method_Flag = 3
END IF


CALL Init_MPI_Params()

CALL Set_Caller_Data(   Source_NQ,                      &
                        Source_xL,                      &
                        Source_RQ_xlocs,                &
                        Source_TQ_xlocs,                &
                        Source_PQ_xlocs,                &
                        Source_Units,                   &
                        Source_Radial_Boundary_Units    )





#ifdef POSEIDON_AMREX_FLAG
    DOMAIN_DIM = amrex_spacedim
    lPF_Core_Flags(iPF_Core_AMReX_Mode) = .TRUE.
    CALL Init_Parameters_From_AMReX_Input_File()
#else
    lPF_Core_Flags(iPF_Core_AMReX_Mode) = .FALSE.

    IF ( PRESENT( Dimensions_Option ) ) THEN
        Domain_Dim = Dimensions_Option
    ELSE
        Domain_Dim = 3
    END IF

    CALL Init_Expansion_Params( FEM_Degree_Option, L_Limit_Option )

    CALL Initialize_Mesh(   Source_NE,                          &
                            Domain_Edge_Option,                 &
                            Source_Radial_Boundary_Units,       &
                            Source_R_Option,                    &
                            Source_T_Option,                    &
                            Source_P_Option,                    &
                            Source_DR_Option,                   &
                            Source_DT_Option,                   &
                            Source_DP_Option                    )

    CALL Init_Fixed_Point_Params(   Max_Iterations_Option,          &
                                    Convergence_Criteria_Option,    &
                                    Anderson_M_Option               )




#endif


CALL Init_IO_Params(    WriteAll_Option,                &
                        Print_Setup_Option,             &
                        Write_Setup_Option,             &
                        Print_Results_Option,           &
                        Write_Results_Option,           &
                        Print_Timetable_Option,         &
                        Write_Timetable_Option,         &
                        Write_Sources_Option,           &
                        Print_Condition_Option,     &
                        Write_Condition_Option,     &
                        Suffix_Flag_Option,             &
                        Suffix_Tail_Option,             &
                        Frame_Option                    )


CALL Init_Quad_Params( Integration_NQ_Option )
CALL Initialize_Quadrature()



Translation_Matrix = Create_Translation_Matrix( Caller_NQ,          &
                                                Caller_xL,          &
                                                Caller_RQ_xlocs,    &
                                                Caller_TQ_xlocs,    &
                                                Caller_PQ_xlocs,    &
                                                Caller_Quad_DOF,    &
                                                [Num_R_Quad_Points, Num_T_Quad_Points, Num_P_Quad_Points ],            &
                                                [xLeftLimit, xRightLimit ],            &
                                                Int_R_Locations,      &
                                                Int_T_Locations,      &
                                                Int_P_Locations,      &
                                                Local_Quad_DOF            )
lPF_Init_MTGV_Flags(iPF_Init_MTGV_TransMat) = .TRUE.



IF ( PRESENT(Flat_Guess_Option) ) THEN
    lPF_IG_Flags(iPF_IG_Flat_Guess) = Flat_Guess_Option
ELSE
    lPF_IG_Flags(iPF_IG_Flat_Guess) = .TRUE.
END IF







IF ( lPF_Core_Flags(iPF_Core_Poisson_Mode) ) THEN
    !=======================================================!
    !                                                       !
    !               Initialize Poisson Solver               !
    !                                                       !
    !=======================================================!
    CALL Initialize_Derived()

    CALL Allocate_Poseidon_Source_Variables()

    CALL Initialize_Tables()

    CALL Initialize_Poisson

ELSE
    !=======================================================!
    !                                                       !
    !           Initialize CFA/XCFC Metric Solver           !
    !                                                       !
    !=======================================================!
    nProcs_Poseidon = 1
    Method_Flag = 3


    IF ( PRESENT(CFA_Eq_Flags_Option) ) THEN
        CFA_EQ_Flags = CFA_Eq_Flags_Option
    ELSE
        CFA_EQ_Flags = [1,1,1,0,0]
    END IF

    NUM_CFA_Eqs = SUM(CFA_EQ_Flags)

#ifdef POSEIDON_AMREX_FLAG
    CALL Initialize_Derived_AMReX()
#else
    CALL Initialize_Derived()
#endif

    CALL Allocate_Poseidon_Source_Variables()
    CALL Allocate_XCFC_Source_Routine_Variables()
    CALL Initialize_Tables()

#ifndef POSEIDON_AMREX_FLAG
    CALL Initialize_XCFC()
#endif

END IF ! Not Poisson Mode



CALL TimerStop( Timer_Initialization_Core )





IF ( Verbose_Flag ) CALL Init_Message('Poseidon Initialization Core Complete.')


IF ( lPF_IO_Flags(iPF_IO_Write_Setup) ) THEN
    CALL Output_Setup_Report()
END IF





END SUBROUTINE Initialize_Poseidon
















END MODULE Poseidon_Interface_Initialization

