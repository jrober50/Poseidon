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
                    Poisson_Mode,           &
                    CFA_EQ_Flags,           &
                    Num_CFA_Eqs

USE Variables_IO, &
            ONLY :  Report_Flags,           &
                    Report_IDs,             &
                    iRF_Setup,              &
                    iRF_Time,               &
                    iWF_Source,             &
                    iWF_Results,            &
                    Write_Flags,            &
                    File_Suffix

USE Variables_Functions, &
            ONLY :  LM_Location

USE Variables_Quadrature, &
            ONLY :  Num_R_Quad_Points,      &
                    Num_T_Quad_Points,      &
                    Num_P_Quad_Points,      &
                    Num_TP_Quad_Points,     &
                    Local_Quad_DOF

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
                    nProcs_Poseidon,        &
                    Num_Blocks_Per_Shell,   &
                    Num_R_Elems_Per_Block,  &
                    Num_T_Elems_Per_Block,  &
                    Num_P_Elems_Per_Block,  &
                    Num_R_Elems_Per_Shell,  &
                    Num_Shells,             &
                    Num_SubShells,          &
                    Num_SubShells_Per_Shell

USE Variables_AMReX_Core, &
            ONLY :  AMReX_Mode


USE Allocation_Core, &
            ONLY :  Allocate_Poseidon_CFA_Variables

USE Allocation_Poisson, &
            ONLY :  Allocate_Poseidon_Poisson_Variables

USE Allocation_Mesh, &
            ONLY :  Allocate_Mesh

USE Initialization_Mesh, &
            ONLY :  Initialize_Mesh

USE Initialization_Quadrature, &
            ONLY :  Initialize_Quadrature

USE Initialization_MPI, &
            ONLY :  Initialize_MPI

USE Initialization_Tables, &
            ONLY :  Initialize_Tables

USE Initialization_Derived, &
            ONLY :  Initialize_Derived

USE Initialization_FP, &
            ONLY :  Initialize_FP,          &
                    Create_Eq_Maps

USE Initialization_XCFC, &
            ONLY :  Initialize_XCFC

!USE Initialization_NR, &
!            ONLY :  Initialize_NR

USE Maps_Legacy, &
            ONLY :  CFA_3D_LM_Map

USE Poseidon_IO_Module, &
            ONLY :  Output_Mesh,            &
                    Output_Nodal_Mesh

USE IO_Setup_Report_Module, &
            ONLY :  Output_Setup_Report

USE Poisson_Matrix_Routines,    &
            ONLY :  Initialize_Stiffness_Matrix

USE Timer_Routines_Module, &
            ONLY :  Init_Timers,                    &
                    Finalize_Timers,                &
                    TimerStart,                     &
                    TimerStop


USE Timer_Variables_Module, &
            ONLY :  Timer_Poisson_Matrix_Init,      &
                    Timer_Initialization_Core


USE Initialization_Subroutines, &
            ONLY :  Init_Fixed_Point_Params,        &
                    Init_IO_Params

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
                                Units_Option,                           &
                                Domain_Edge_Option,                     &
                                NE_Option,                              &
                                NQ_Option,                              &
                                Coarsen_Option,                         &
                                r_Option, t_Option, p_Option,           &
                                dr_Option, dt_Option, dp_Option,        &
                                Method_Flag_Option,                     &
                                CFA_Eq_Flags_Option,                    &
                                nProcs_Option,                          &
                                Suffix_Flag_Option,                     &
                                Suffix_Tail_Option,                     &
                                Frame_Option,                           &
                                Max_Iterations_Option,                  &
                                Convergence_Criteria_Option,            &
                                Anderson_M_Option,                      &
                                Poisson_Mode_Option,                    &
                                AMReX_Mode_Option,                      &
                                Verbose_Option,                         &
                                WriteAll_Option,                        &
                                Print_Setup_Option,                     &
                                Write_Setup_Option,                     &
                                Print_Results_Option,                   &
                                Write_Results_Option,                   &
                                Print_Timetable_Option,                 &
                                Write_Timetable_Option,                 &
                                Write_Sources_Option                  )


CHARACTER(LEN=1),        INTENT(IN), OPTIONAL               ::  Units_Option

INTEGER,                 INTENT(IN), OPTIONAL               ::  FEM_Degree_Option
INTEGER,                 INTENT(IN), OPTIONAL               ::  L_Limit_Option

INTEGER,   DIMENSION(3), INTENT(IN), OPTIONAL               ::  NQ_Option
INTEGER,   DIMENSION(3), INTENT(IN), OPTIONAL               ::  NE_Option
REAL(idp), DIMENSION(2), INTENT(IN), OPTIONAL               ::  Domain_Edge_Option

INTEGER,                 INTENT(IN), OPTIONAL               ::  Method_Flag_Option

INTEGER,                 INTENT(IN), OPTIONAL               ::  Dimensions_Option

INTEGER,   DIMENSION(3), INTENT(IN), OPTIONAL               ::  Coarsen_Option

REAL(idp), DIMENSION(:), INTENT(IN), OPTIONAL               ::  r_Option
REAL(idp), DIMENSION(:), INTENT(IN), OPTIONAL               ::  t_Option
REAL(idp), DIMENSION(:), INTENT(IN), OPTIONAL               ::  p_Option
REAL(idp), DIMENSION(:), INTENT(IN), OPTIONAL               ::  dr_Option
REAL(idp), DIMENSION(:), INTENT(IN), OPTIONAL               ::  dt_Option
REAL(idp), DIMENSION(:), INTENT(IN), OPTIONAL               ::  dp_Option

INTEGER,   DIMENSION(5), INTENT(IN), OPTIONAL               ::  CFA_EQ_Flags_Option
INTEGER,                 INTENT(IN), OPTIONAL               ::  nProcs_Option

CHARACTER(LEN=10),       INTENT(IN), OPTIONAL               ::  Suffix_Flag_Option
CHARACTER(LEN=1),        INTENT(IN), OPTIONAL               ::  Suffix_Tail_Option
INTEGER,                 INTENT(IN), OPTIONAL               ::  Frame_Option

INTEGER,                 INTENT(IN), OPTIONAL               ::  Max_Iterations_Option

REAL(idp),               INTENT(IN), OPTIONAL               ::  Convergence_Criteria_Option
INTEGER,                 INTENT(IN), OPTIONAL               ::  Anderson_M_Option
LOGICAL,                 INTENT(IN), OPTIONAL               ::  Poisson_Mode_Option
LOGICAL,                 INTENT(IN), OPTIONAL               ::  AMReX_Mode_Option

LOGICAL,                 INTENT(IN), OPTIONAL               ::  Verbose_Option
LOGICAL,                 INTENT(IN), OPTIONAL               ::  WriteAll_Option
LOGICAL,                 INTENT(IN), OPTIONAL               ::  Print_Setup_Option
LOGICAL,                 INTENT(IN), OPTIONAL               ::  Write_Setup_Option
LOGICAL,                 INTENT(IN), OPTIONAL               ::  Print_Results_Option
LOGICAL,                 INTENT(IN), OPTIONAL               ::  Write_Results_Option
LOGICAL,                 INTENT(IN), OPTIONAL               ::  Print_Timetable_Option
LOGICAL,                 INTENT(IN), OPTIONAL               ::  Write_Timetable_Option
LOGICAL,                 INTENT(IN), OPTIONAL               ::  Write_Sources_Option


CALL Init_Timers
CALL TimerStart( Timer_Initialization_Core )


IF ( PRESENT( Verbose_Option ) ) THEN
    Verbose_Flag = Verbose_Option
ELSE
    Verbose_Flag = .FALSE.
END IF

IF ( Verbose_Flag ) THEN
    PRINT*,"Initializing Poseidon..."
END IF







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








IF ( PRESENT( Dimensions_Option ) ) THEN
    Domain_Dim = Dimensions_Option
ELSE
    Domain_Dim = 3
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
Local_Quad_DOF     = Num_R_Quad_Points*Num_TP_Quad_Points











!-------------------------------------------!
!                                           !
!                   Mesh                    !
!                                           !
!-------------------------------------------!
IF ( PRESENT( NE_Option ) ) THEN
    Num_R_Elements = NE_Option(1)
    Num_T_Elements = NE_Option(2)
    Num_P_Elements = NE_Option(3)
ELSE
    Num_R_Elements = 1
    Num_T_Elements = 1
    Num_P_Elements = 1
END IF



IF ( PRESENT( Coarsen_Option) ) THEN
    R_Coarsen_Factor = Coarsen_Option(1)
    T_Coarsen_Factor = Coarsen_Option(2)
    P_Coarsen_Factor = Coarsen_Option(3)
ELSE
    R_Coarsen_Factor = 1
    T_Coarsen_Factor = 1
    P_Coarsen_Factor = 1
END IF
Num_Loc_R_Elements = Num_R_Elements
Num_Loc_T_Elements = Num_R_Elements
Num_Loc_P_Elements = Num_R_Elements




CALL Allocate_Mesh()

IF ( PRESENT( r_Option ) ) THEN
    rlocs = r_Option
    locs_set(1) = .TRUE.
END IF
IF ( PRESENT( dr_Option ) ) THEN
    drlocs = dr_Option
    dlocs_set(1) = .TRUE.
END IF

IF ( PRESENT( t_Option ) ) THEN
    tlocs = t_Option
    locs_set(2) = .TRUE.
END IF

IF ( PRESENT( Domain_Edge_Option ) ) THEN
    R_Inner = Domain_Edge_Option(1)
    R_Outer = Domain_Edge_Option(2)
ELSE
    R_Inner = 0.0_idp
    R_Outer = 1.0_idp
END IF


IF ( PRESENT( dt_Option ) ) THEN
    dtlocs = dt_Option
    dlocs_set(2) = .TRUE.
ELSE
    dtlocs = pi/Num_T_Elements
    dlocs_set(2) = .TRUE.
END IF


IF ( PRESENT( p_Option ) ) THEN
    plocs = p_Option
    locs_set(3) = .TRUE.
END IF
IF ( PRESENT( dp_Option ) ) THEN
    dplocs = dp_Option
    dlocs_set(3) = .TRUE.
ELSE
    dplocs = 2.0_idp*pi/Num_P_Elements
    dlocs_set(3) = .TRUE.
END IF


CALL Initialize_Mesh( )


!
!   Verbose Options
!
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


    LM_Location => CFA_3D_LM_Map
    CALL Initialize_Derived()
    CALL Initialize_MPI()
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


    IF ( PRESENT(AMReX_Mode_Option) ) THEN
        AMReX_Mode = AMReX_Mode_Option
    ELSE
        AMReX_Mode = .FALSE.
    END IF



    LM_Location => CFA_3D_LM_Map



    IF ( PRESENT( Method_Flag_Option ) ) THEN
        Method_Flag = Method_Flag_Option
    ELSE
        Method_Flag = 2
    END IF






    CALL Initialize_Derived()
    CALL Initialize_MPI()
    CALL Allocate_Poseidon_CFA_Variables()
    CALL Initialize_Quadrature()

    CALL Initialize_Tables()


    IF ( Method_Flag == 1 ) THEN
        WRITE(*,'(A)')"The Newton-Raphson method is currently unavailable in Poseidon.  STOPING"
        STOP
!        CALL Initialize_NR(CFA_EQ_Flags_Option)
    ELSE IF ( Method_Flag == 2 ) THEN
        CALL Initialize_FP()
    ELSE IF ( Method_Flag == 3 ) THEN
        CALL Initialize_XCFC()
    END IF




END IF ! Not Poisson Mode






IF ( Verbose_Flag ) THEN
    PRINT*,"Poseidon Initialization Complete"
END IF


IF ( Write_Flags(7) == 1 ) THEN
    PRINT*,"Outputing Radial Mesh to file during initialization."
    CALL Output_Mesh( rlocs, NUM_R_ELEMENTS+1 )
    CALL Output_Nodal_Mesh( rlocs, NUM_R_ELEMENTS+1 )
END IF




CALL Output_Setup_Report()


CALL TimerStop( Timer_Initialization_Core )

END SUBROUTINE Initialize_Poseidon




































 !+201+####################################################################################!
!                                                                                           !
!       Initialize_Poseidon_From_File                                                       !
!                                                                                           !
 !#########################################################################################!
SUBROUTINE Initialize_Poseidon_From_File()




END SUBROUTINE Initialize_Poseidon_From_File










END MODULE Initialization_Poseidon
