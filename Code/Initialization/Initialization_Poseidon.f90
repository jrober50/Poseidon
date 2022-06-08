   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Initialization_Poseidon                                                      !##!
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
                    Poisson_Mode,           &
                    CFA_EQ_Flags,           &
                    Num_CFA_Eqs

USE Variables_Functions, &
            ONLY :  LM_Location

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

USE Variables_Mesh, &
            ONLY :  Num_R_Elements,         &
                    Num_T_Elements,         &
                    Num_P_Elements,         &
                    Num_Loc_R_Elements,     &
                    Num_Loc_T_Elements,     &
                    Num_Loc_P_Elements,     &
                    rlocs,                  &
                    tlocs,                  &
                    plocs,                  &
                    drlocs,                 &
                    dtlocs,                 &
                    dplocs,                 &
                    R_Inner,                &
                    R_Outer,                &
                    R_Coarsen_Factor,       &
                    T_Coarsen_Factor,       &
                    P_Coarsen_Factor,       &
                    locs_set,               &
                    dlocs_set

USE Variables_Derived, &
            ONLY :  Prob_Dim,               &
                    Block_Prob_Dim,         &
                    SubShell_Prob_Dim,      &
                    Elem_Prob_Dim_Sqr,      &
                    Var_Dim,                &
                    Block_Var_Dim,          &
                    Num_Off_Diagonals,      &
                    ULM_Length,             &
                    Num_R_Nodes
                
USE Variables_MPI, &
            ONLY :  myID_Poseidon,          &
                    nProcs_Poseidon

USE Variables_AMReX_Core, &
            ONLY :  AMReX_Mode


USE Allocation_Sources, &
            ONLY :  Allocate_Poseidon_Source_Variables

USE XCFC_Source_Routine_Variables_Module, &
            ONLY :  Allocate_XCFC_Source_Routine_Variables

USE Allocation_Poisson, &
            ONLY :  Allocate_Poseidon_Poisson_Variables


USE Initialization_Subroutines, &
            ONLY :  Init_Fixed_Point_Params,        &
                    Init_Expansion_Params,          &
                    Init_Mesh_Params,               &
                    Init_IO_Params,                 &
                    Init_MPI_Params,                &
                    Init_Quad_Params,               &
                    Set_Caller_Data,                &
                    Initialize_Mesh

USE Initialization_Quadrature, &
            ONLY :  Initialize_Quadrature


USE Initialization_Tables, &
            ONLY :  Initialize_Tables

USE Initialization_Derived, &
            ONLY :  Initialize_Derived

USE Initialization_FP, &
            ONLY :  Initialize_FP,          &
                    Create_Eq_Maps

USE Initialization_XCFC, &
            ONLY :  Initialize_XCFC

USE Maps_Legacy, &
            ONLY :  CFA_3D_LM_Map

USE IO_Output_Mesh_Module, &
            ONLY :  Output_Mesh,            &
                    Output_Nodal_Mesh

USE IO_Setup_Report_Module, &
            ONLY :  Output_Setup_Report

USE Poisson_Matrix_Routines,    &
            ONLY :  Initialize_Stiffness_Matrix

USE Variables_Interface, &
            ONLY :  Caller_NQ,                  &
                    Caller_Quad_DOF,            &
                    Caller_xL,                  &
                    Caller_RQ_xlocs,            &
                    Caller_TQ_xlocs,            &
                    Caller_PQ_xlocs,            &
                    Caller_R_Units,             &
                    Translation_Matrix

USE Functions_Translation_Matrix_Module, &
            ONLY :  Create_Translation_Matrix


USE Timer_Routines_Module, &
            ONLY :  Init_Timers,                    &
                    Finalize_Timers,                &
                    TimerStart,                     &
                    TimerStop


USE Timer_Variables_Module, &
            ONLY :  Timer_Poisson_Matrix_Init,      &
                    Timer_Initialization_Core

USE Flags_IO_Module, &
            ONLY :  lPF_IO_Flags,                   &
                    iPF_IO_Write_Mesh

USE Flags_Initial_Guess_Module, &
            ONLY :  lPF_IG_Flags,                   &
                    iPF_IG_Flat_Guess

USE Flags_Initialization_Module, &
            ONLY :  lPF_Init_MTGV_Flags,    &
                    iPF_Init_MTGV_TransMat

IMPLICIT NONE


PUBLIC :: Initialize_Poseidon


CONTAINS

 !+101+####################################################################################!
!                                                                                           !
!       Initialize_Poseidon                                                                 !
!                                                                                           !
 !#########################################################################################!
SUBROUTINE Initialize_Poseidon( Dimensions_Option,                      &
                                FEM_Degree_Option,                      &
                                L_Limit_Option,                         &
                                Source_NE,                              &
                                Domain_Edge_Option,                     &
                                Source_NQ,                              &
                                Source_xL,                              &
                                Source_RQ_xlocs,                        &
                                Source_TQ_xlocs,                        &
                                Source_PQ_xlocs,                        &
                                Source_Units,                           &
                                Source_Radial_Boundary_Units,           &
                                Integration_NQ_Option,                  &
                                Source_R_Option,                        &
                                Source_T_Option,                        &
                                Source_P_Option,                        &
                                Source_DR_Option,                       &
                                Source_DT_Option,                       &
                                Source_DP_Option,                       &
                                Method_Flag_Option,                     &
                                CFA_Eq_Flags_Option,                    &
                                nProcs_Option,                          &
                                Max_Iterations_Option,                  &
                                Convergence_Criteria_Option,            &
                                Anderson_M_Option,                      &
                                Poisson_Mode_Option,                    &
                                Flat_Guess_Option,                      &
                                Verbose_Option,                         &
                                WriteAll_Option,                        &
                                Print_Setup_Option,                     &
                                Write_Setup_Option,                     &
                                Print_Results_Option,                   &
                                Write_Results_Option,                   &
                                Print_Timetable_Option,                 &
                                Write_Timetable_Option,                 &
                                Write_Sources_Option,                   &
                                Print_Condition_Option,                 &
                                Write_Condition_Option,                 &
                                Suffix_Flag_Option,                     &
                                Suffix_Tail_Option,                     &
                                Frame_Option                            )


INTEGER,    DIMENSION(3),               INTENT(IN)              ::  Source_NQ
REAL(idp),  DIMENSION(2),               INTENT(IN)              ::  Source_xL
REAL(idp),  DIMENSION(Source_NQ(1)),    INTENT(IN)              ::  Source_RQ_xlocs
REAL(idp),  DIMENSION(Source_NQ(2)),    INTENT(IN)              ::  Source_TQ_xlocs
REAL(idp),  DIMENSION(Source_NQ(3)),    INTENT(IN)              ::  Source_PQ_xlocs
CHARACTER(LEN=1),                       INTENT(IN)              ::  Source_Units
CHARACTER(LEN=2),                       INTENT(IN)              ::  Source_Radial_Boundary_Units

INTEGER,                 INTENT(IN), OPTIONAL               ::  FEM_Degree_Option
INTEGER,                 INTENT(IN), OPTIONAL               ::  L_Limit_Option

INTEGER,   DIMENSION(3), INTENT(IN), OPTIONAL               ::  Source_NE
REAL(idp), DIMENSION(2), INTENT(IN), OPTIONAL               ::  Domain_Edge_Option

INTEGER,                 INTENT(IN), OPTIONAL               ::  Method_Flag_Option

INTEGER,                 INTENT(IN), OPTIONAL               ::  Dimensions_Option
INTEGER,   DIMENSION(3), INTENT(IN), OPTIONAL               ::  Integration_NQ_Option


REAL(idp), DIMENSION(:), INTENT(IN), OPTIONAL               ::  Source_R_Option
REAL(idp), DIMENSION(:), INTENT(IN), OPTIONAL               ::  Source_T_Option
REAL(idp), DIMENSION(:), INTENT(IN), OPTIONAL               ::  Source_P_Option
REAL(idp), DIMENSION(:), INTENT(IN), OPTIONAL               ::  Source_DR_Option
REAL(idp), DIMENSION(:), INTENT(IN), OPTIONAL               ::  Source_DT_Option
REAL(idp), DIMENSION(:), INTENT(IN), OPTIONAL               ::  Source_DP_Option

INTEGER,   DIMENSION(5), INTENT(IN), OPTIONAL               ::  CFA_EQ_Flags_Option
INTEGER,                 INTENT(IN), OPTIONAL               ::  nProcs_Option

CHARACTER(LEN=10),       INTENT(IN), OPTIONAL               ::  Suffix_Flag_Option
CHARACTER(LEN=1),        INTENT(IN), OPTIONAL               ::  Suffix_Tail_Option
INTEGER,                 INTENT(IN), OPTIONAL               ::  Frame_Option

INTEGER,                 INTENT(IN), OPTIONAL               ::  Max_Iterations_Option

REAL(idp),               INTENT(IN), OPTIONAL               ::  Convergence_Criteria_Option
INTEGER,                 INTENT(IN), OPTIONAL               ::  Anderson_M_Option
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

CALL Init_Timers
CALL TimerStart( Timer_Initialization_Core )


IF ( PRESENT( Verbose_Option ) ) THEN
    Verbose_Flag = Verbose_Option
ELSE
    Verbose_Flag = .FALSE.
END IF
IF ( Verbose_Flag ) CALL Init_Message('Beginning Poseidon Core Initialization.')



CALL Set_Units(Source_Units)


IF ( PRESENT( Dimensions_Option ) ) THEN
    Domain_Dim = Dimensions_Option
ELSE
    Domain_Dim = 3
END IF

AMReX_Mode = .FALSE.

IF ( PRESENT(Poisson_Mode_Option) ) THEN
    Poisson_Mode = Poisson_Mode_Option
ELSE
    Poisson_Mode = .FALSE.
END IF




CALL Init_MPI_Params()

CALL Set_Caller_Data(   Source_NQ,                      &
                        Source_xL,                      &
                        Source_RQ_xlocs,                &
                        Source_TQ_xlocs,                &
                        Source_PQ_xlocs,                &
                        Source_Units,                   &
                        Source_Radial_Boundary_Units    )

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



CALL Init_IO_Params(WriteAll_Option,            &
                    Print_Setup_Option,         &
                    Write_Setup_Option,         &
                    Print_Results_Option,       &
                    Write_Results_Option,       &
                    Print_Timetable_Option,     &
                    Write_Timetable_Option,     &
                    Write_Sources_Option,       &
                    Print_Condition_Option,     &
                    Write_Condition_Option,     &
                    Suffix_Flag_Option,         &
                    Suffix_Tail_Option,         &
                    Frame_Option                )


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





IF ( Poisson_Mode ) THEN
    !=======================================================!
    !                                                       !
    !               Initialize Poisson Solver               !
    !                                                       !
    !=======================================================!


    LM_Location => CFA_3D_LM_Map
    CALL Initialize_Derived()
    CALL Initialize_Quadrature()
!    CALL Initialize_Tables()

    CALL Allocate_Poseidon_Poisson_Variables()



    CALL TimerStart( Timer_Poisson_Matrix_Init )
    CALL Initialize_Stiffness_Matrix()
    CALL TimerStop(  Timer_Poisson_Matrix_Init )




ELSE
    !=======================================================!
    !                                                       !
    !           Initialize CFA/XCFC Metric Solver           !
    !                                                       !
    !=======================================================!


    CALL Init_Fixed_Point_Params( Max_Iterations_Option,          &
                                  Convergence_Criteria_Option,    &
                                  Anderson_M_Option               )




    IF ( PRESENT(CFA_Eq_Flags_Option) ) THEN
        CFA_EQ_Flags = CFA_Eq_Flags_Option
    ELSE
        CFA_EQ_Flags = [1,1,1,0,0]
    END IF

    NUM_CFA_Eqs = SUM(CFA_EQ_Flags)
    CALL Create_Eq_Maps()

    IF ( PRESENT(nProcs_Option) ) THEN
        nProcs_Poseidon = nProcs_Option
    ELSE
        nProcs_Poseidon = 1
    END IF

    LM_Location => CFA_3D_LM_Map



    IF ( PRESENT( Method_Flag_Option ) ) THEN
        Method_Flag = Method_Flag_Option
    ELSE
        Method_Flag = 2
    END IF






    CALL Initialize_Derived()
    CALL Allocate_Poseidon_Source_Variables()
    CALL Allocate_XCFC_Source_Routine_Variables()
    CALL Initialize_Tables()


    CALL Initialize_XCFC()




END IF ! Not Poisson Mode

CALL TimerStop( Timer_Initialization_Core )




IF ( Verbose_Flag ) CALL Init_Message('Poseidon Initialization Core Complete.')


IF ( lPF_IO_Flags(iPF_IO_Write_Mesh) ) THEN
    CALL Output_Mesh( rlocs, NUM_R_ELEMENTS+1 )
    CALL Output_Nodal_Mesh( rlocs, NUM_R_ELEMENTS+1 )
END IF




CALL Output_Setup_Report()



END SUBROUTINE Initialize_Poseidon




































 !+201+####################################################################################!
!                                                                                           !
!       Initialize_Poseidon_From_File                                                       !
!                                                                                           !
 !#########################################################################################!
SUBROUTINE Initialize_Poseidon_From_File()




END SUBROUTINE Initialize_Poseidon_From_File










END MODULE Initialization_Poseidon
