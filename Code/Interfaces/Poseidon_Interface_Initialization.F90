   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Poseidon_Interface_Initialization                                     !##!
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
            
USE Poseidon_Message_Routines_Module, &
            ONLY :  Warning_Message
            
USE Poseidon_Units_Module, &
            ONLY :  Set_Units

USE Poseidon_Parameters, &
            ONLY :  Domain_Dim,             &
                    Degree,                 &
                    L_Limit,                &
                    Method_Flag,            &
                    Verbose_Flag,           &
                    Convergence_Criteria,   &
                    Max_Iterations,         &
                    Eq_Flags,           &
                    Num_Eqs

USE Allocation_Sources, &
            ONLY :  Allocate_Poseidon_Source_Variables

USE XCFC_Source_Routine_Variables_Module, &
            ONLY :  Allocate_XCFC_Source_Routine_Variables

USE Variables_MPI, &
            ONLY :  nProcs_Poseidon,        &
                    myID_Poseidon,          &
                    MasterID_Poseidon,      &
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
            ONLY :  Initialize_Derived,             &
                    Initialize_Derived_AMReX_Part1

USE Initialization_Subroutines, &
            ONLY :  Init_Fixed_Point_Params,        &
                    Init_IO_Params,                 &
                    Init_MPI_Params,                &
                    Init_Quad_Params,               &
                    Init_Expansion_Params,          &
                    Init_Mesh_Params,               &
                    Initialize_Mesh,                &
                    Set_Caller_Data,                &
                    Set_Method_Flags


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

USE Variables_FP, &
            ONLY :  FP_Diagnostics_Flag


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

USE Initialization_FEM_Module, &
            ONLY :  Initialization_FEM

USE Variables_Interface, &
            ONLY :  Caller_NQ,                  &
                    Caller_Quad_DOF,            &
                    Caller_xL,                  &
                    Caller_RQ_xlocs,            &
                    Caller_TQ_xlocs,            &
                    Caller_PQ_xlocs,            &
                    Caller_R_Units,             &
                    Translation_Matrix

USE Flags_Core_Module, &
            ONLY :  iPF_Core_Flags,             &
                    iPF_Core_Unit_Mode,         &
                    iPF_Core_Method_Mode,       &
                    iPF_Core_AMReX_Mode,        &
                    iPF_Core_Units_CGS,         &
                    iPF_Core_Units_MKS,         &
                    iPF_Core_Units_Geometrized, &
                    iPF_Core_Units_Unitless,    &
                    iPF_Core_Method_Newtonian,  &
                    iPF_Core_Method_XCFC,       &
                    iPF_Core_AMReX_Off,         &
                    iPF_Core_AMReX_On

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

                    
USE Poseidon_Git_IO_Module, &
            ONLY :  Poseidon_Git_Output_Info



IMPLICIT NONE


CONTAINS




 !+101+############################################################!
!                                                                   !
!          Initialize_Poseidon                                      !
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
                                FEM_Degree_Option,                      &
                                L_Limit_Option,                         &
                                Source_R_Option,                        &
                                Source_T_Option,                        &
                                Source_P_Option,                        &
                                Source_DR_Option,                       &
                                Source_DT_Option,                       &
                                Source_DP_Option,                       &
                                Integration_NQ_Option,              &
                                Max_Iterations_Option,              &
                                Convergence_Criteria_Option,        &
                                Anderson_M_Option,                  &
                                Fixed_Point_Diagnostics_Option,     &
                                Eq_Flags_Option,                    &
                                AMReX_FEM_Refinement_Option,        &
                                AMReX_Integral_Refinement_Option,   &
                                Newtonian_Mode_Option,              &
                                CFA_Mode_Option,                    &
                                XCFC_Mode_Option,                   &
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
                                Write_FP_Diagnostics_Option,        &
                                Suffix_Flag_Option,                 &
                                Suffix_Param_Type_Option,           &
                                Suffix_Tail_Option,                 &
                                Suffix_Input_Option,                &
                                Frame_Option                        )



INTEGER,    DIMENSION(3),               INTENT(IN)          ::  Source_NQ
REAL(idp),  DIMENSION(2),               INTENT(IN)          ::  Source_xL
REAL(idp),  DIMENSION(Source_NQ(1)),    INTENT(IN)          ::  Source_RQ_xlocs
REAL(idp),  DIMENSION(Source_NQ(2)),    INTENT(IN)          ::  Source_TQ_xlocs
REAL(idp),  DIMENSION(Source_NQ(3)),    INTENT(IN)          ::  Source_PQ_xlocs
CHARACTER(LEN=1),                       INTENT(IN)          ::  Source_Units
CHARACTER(LEN=2),                       INTENT(IN)          ::  Source_Radial_Boundary_Units
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

INTEGER,   DIMENSION(3), INTENT(IN), OPTIONAL               ::  Integration_NQ_Option
INTEGER,   DIMENSION(5), INTENT(IN), OPTIONAL               ::  Eq_Flags_Option

INTEGER,                 INTENT(IN), OPTIONAL               ::  AMReX_FEM_Refinement_Option
INTEGER,                 INTENT(IN), OPTIONAL               ::  AMReX_Integral_Refinement_Option

LOGICAL,                 INTENT(IN), OPTIONAL               ::  Newtonian_Mode_Option
LOGICAL,                 INTENT(IN), OPTIONAL               ::  CFA_Mode_Option
LOGICAL,                 INTENT(IN), OPTIONAL               ::  XCFC_Mode_Option

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

LOGICAL,                 INTENT(IN), OPTIONAL               ::  Write_FP_Diagnostics_Option

CHARACTER(LEN=10),       INTENT(IN), OPTIONAL               ::  Suffix_Flag_Option
INTEGER,                 INTENT(IN), OPTIONAL               ::  Suffix_Param_Type_Option
CHARACTER(LEN=4),        INTENT(IN), OPTIONAL               ::  Suffix_Tail_Option
CHARACTER(LEN=*),        INTENT(IN), OPTIONAL               ::  Suffix_Input_Option
INTEGER,                 INTENT(IN), OPTIONAL               ::  Frame_Option

INTEGER,                 INTENT(IN), OPTIONAL               ::  Max_Iterations_Option
REAL(idp),               INTENT(IN), OPTIONAL               ::  Convergence_Criteria_Option
INTEGER,                 INTENT(IN), OPTIONAL               ::  Anderson_M_Option
LOGICAL,                 INTENT(IN), OPTIONAL               ::  Fixed_Point_Diagnostics_Option


CALL Init_Timers
CALL TimerStart( Timer_Initialization_Core )


IF ( PRESENT( Verbose_Option ) ) THEN
    Verbose_Flag = Verbose_Option
ELSE
    Verbose_Flag = .FALSE.
END IF
IF ( Verbose_Flag ) CALL Init_Message('Beginning Poseidon Core Initialization.')





IF ( Source_Units == "C" ) THEN
    iPF_Core_Flags(iPF_Core_Unit_Mode) = iPF_Core_Units_CGS
ELSE IF ( Source_Units == "S" ) THEN
    iPF_Core_Flags(iPF_Core_Unit_Mode) = iPF_Core_Units_MKS
ELSE IF ( Source_Units == "G" ) THEN
    iPF_Core_Flags(iPF_Core_Unit_Mode) = iPF_Core_Units_Geometrized
ELSE IF ( Source_Units == "U" ) THEN
    iPF_Core_Flags(iPF_Core_Unit_Mode) = iPF_Core_Units_Unitless
END IF
CALL Set_Units(Source_Units)



CALL Set_Method_Flags(  Newtonian_Mode_Option,      &
                        CFA_Mode_Option,            &
                        XCFC_Mode_Option            )



CALL Init_MPI_Params()

IF ( myID_Poseidon == MasterID_Poseidon ) THEN
    CALL Poseidon_Git_Output_Info()
END IF

CALL Set_Caller_Data(   Source_NQ,                      &
                        Source_xL,                      &
                        Source_RQ_xlocs,                &
                        Source_TQ_xlocs,                &
                        Source_PQ_xlocs,                &
                        Source_Units,                   &
                        Source_Radial_Boundary_Units    )





#ifdef POSEIDON_AMREX_FLAG

    DOMAIN_DIM = amrex_spacedim
    iPF_Core_Flags(iPF_Core_AMReX_Mode) = iPF_Core_AMReX_On
    CALL Init_Parameters_From_AMReX_Input_File()
    
    
#else
    iPF_Core_Flags(iPF_Core_AMReX_Mode) = iPF_Core_AMReX_Off

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


IF ( PRESENT(Fixed_Point_Diagnostics_Option) ) THEN
    FP_Diagnostics_Flag = Fixed_Point_Diagnostics_Option
ELSE
    FP_Diagnostics_Flag = .FALSE.
END IF




CALL Init_IO_Params(    WriteAll_Option,                &
                        Print_Setup_Option,             &
                        Write_Setup_Option,             &
                        Print_Results_Option,           &
                        Write_Results_Option,           &
                        Print_Timetable_Option,         &
                        Write_Timetable_Option,         &
                        Write_Sources_Option,           &
                        Print_Condition_Option,         &
                        Write_Condition_Option,         &
                        Write_FP_Diagnostics_Option,    &
                        Suffix_Flag_Option,             &
                        Suffix_Param_Type_Option,       &
                        Suffix_Tail_Option,             &
                        Suffix_Input_Option,            &
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


CALL Initialization_FEM()




IF ( iPF_Core_Flags(iPF_Core_Method_Mode) == iPF_Core_Method_Newtonian ) THEN
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
!    nProcs_Poseidon = 1
    Method_Flag = 3


    IF ( PRESENT(Eq_Flags_Option) ) THEN
        Eq_Flags = Eq_Flags_Option
    ELSE
        Eq_Flags = [1,1,1,0,0]
    END IF

    Num_Eqs = SUM(Eq_Flags)

#ifdef POSEIDON_AMREX_FLAG
    CALL Initialize_Derived_AMReX_Part1()
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

CALL Output_Setup_Report()





END SUBROUTINE Initialize_Poseidon
















END MODULE Poseidon_Interface_Initialization

