   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Initialization_Subroutines                                            !##!
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
            ONLY :  Centimeter,     &
                    Meter,          &
                    Kilometer

USE Poseidon_Message_Routines_Module, &
            ONLY :  Init_Message

USE Poseidon_Parameters, &
            ONLY :  Degree,                         &
                    Degree_Default,                 &
                    L_Limit,                        &
                    L_Limit_Default,                &
                    Verbose_Flag,                   &
                    Convergence_Criteria,           &
                    Convergence_Criteria_Default,   &
                    Max_Iterations,                 &
                    Max_Iterations_Default

USE Variables_Mesh, &
            ONLY :  Num_R_Elements,         &
                    Num_T_Elements,         &
                    Num_P_Elements,         &
                    R_Inner,                &
                    R_Outer,                &
                    rlocs,                  &
                    tlocs,                  &
                    plocs,                  &
                    drlocs,                 &
                    dtlocs,                 &
                    dplocs

USE Variables_IO, &
            ONLY :  File_Suffix

USE Variables_FP, &
            ONLY :  FP_Anderson_M,          &
                    FP_Anderson_M_Default,  &
                    FP_Diagnostics_Flag

USE Variables_AMReX_Core, &
            ONLY :  AMReX_Max_Grid_Size,    &
                    AMReX_Max_Level,        &
                    AMReX_Num_Levels,       &
                    iFRL,                   &
                    iIRL,                   &
                    MF_Source_nComps

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


USE Variables_Interface, &
            ONLY :  Caller_Set,                 &
                    Caller_NQ,                  &
                    Caller_Quad_DOF,            &
                    Caller_xL,                  &
                    Caller_RQ_xlocs,            &
                    Caller_TQ_xlocs,            &
                    Caller_PQ_xlocs,            &
                    Caller_R_Units,             &
                    Translation_Matrix

USE Flags_Core_Module, &
            ONLY :  iPF_Core_Flags,             &
                    iPF_Core_Method_Mode,       &
                    iPF_Core_Method_Newtonian,  &
                    iPF_Core_Method_CFA,        &
                    iPF_Core_Method_XCFC,       &
                    iPF_Core_Method_Too_Many


USE Flags_IO_Module, &
            ONLY :  lPF_IO_Flags,               &
                    iPF_IO_Verbose,             &
                    iPF_IO_Print_Setup,         &
                    iPF_IO_Write_Setup,         &
                    iPF_IO_Print_Results,       &
                    iPF_IO_Write_Results,       &
                    iPF_IO_Print_Timetable,     &
                    iPF_IO_Write_Timetable,     &
                    iPF_IO_Write_Sources,       &
                    iPF_IO_Print_Cond,          &
                    iPF_IO_Write_Cond,          &
                    iPF_IO_Write_FP_Diagnostics

USE Flags_Initialization_Module, &
            ONLY :  lPF_Init_Flags,             &
                    lPF_Init_AMReX_Flags,       &
                    lPF_Init_Quad_Flags,        &
                    lPF_Init_Mesh_Flags,        &
                    iPF_Init_IO_Params,         &
                    iPF_Init_FP_Params,         &
                    iPF_Init_Caller_Vars,       &
                    iPF_Init_Quad_Params,       &
                    iPF_Init_MPI,               &
                    iPF_Init_Expansion_Params,  &
                    iPF_Init_AMReX_Params,      &
                    iPF_Init_Mesh_Params,       &
                    iPF_Init_Mesh_Init

USE Variables_MPI, &
            ONLY :  nProcs_Poseidon,        &
                    myID_Poseidon,          &
                    Poseidon_Comm_World,    &
                    ierr


USE Functions_Math, &
            ONLY :  Lagrange_Poly


USE Functions_Translation_Matrix_Module, &
            ONLY :  Create_Translation_Matrix

USE IO_Suffix_Module, &
            ONLY :  Create_Suffix

USE Maps_X_Space, &
            ONLY :  Map_To_X_Space

USE Allocation_Mesh, &
            ONLY :  Allocate_Mesh

USE MPI


IMPLICIT NONE


CONTAINS
 !+101+########################################################!
!                                                               !
!       Init_Expansion_Params                                   !
!                                                               !
 !#############################################################!
SUBROUTINE Init_Expansion_Params(   Degree_Option,   &
                                    L_Limit_Option   )

INTEGER,            INTENT(IN), OPTIONAL               ::  Degree_Option
INTEGER,            INTENT(IN), OPTIONAL               ::  L_Limit_Option

IF ( Verbose_Flag ) CALL Init_Message('Setting Expansion Parameters.')

IF ( PRESENT(Degree_Option) ) THEN
    Degree = Degree_Option
ELSE
    Degree = Degree_Default
END IF


IF ( PRESENT(L_Limit_Option) ) THEN
    L_Limit = L_Limit_Option
ELSE
    L_Limit = L_Limit_Default
END IF


lPF_Init_Flags(iPF_Init_Expansion_Params) = .TRUE.

END SUBROUTINE Init_Expansion_Params




 !+101+############################################################################!
!                                                                                   !
!       Init_IO_Params                                                              !
!                                                                                   !
 !#################################################################################!
SUBROUTINE Init_IO_Params(  WriteAll_Option,            &
                            Print_Setup_Option,         &
                            Write_Setup_Option,         &
                            Print_Results_Option,       &
                            Write_Results_Option,       &
                            Print_Timetable_Option,     &
                            Write_Timetable_Option,     &
                            Write_Sources_Option,       &
                            Print_Condition_Option,     &
                            Write_Condition_Option,     &
                            Write_FP_Diagnostics_Option,&
                            Suffix_Flag_Option,         &
                            Suffix_Param_Type_Option,   &
                            Suffix_Tail_Option,         &
                            Suffix_Input_Option,        &
                            Frame_Option                )


LOGICAL,            INTENT(IN), OPTIONAL               ::  WriteAll_Option
LOGICAL,            INTENT(IN), OPTIONAL               ::  Print_Setup_Option
LOGICAL,            INTENT(IN), OPTIONAL               ::  Write_Setup_Option
LOGICAL,            INTENT(IN), OPTIONAL               ::  Print_Results_Option
LOGICAL,            INTENT(IN), OPTIONAL               ::  Write_Results_Option
LOGICAL,            INTENT(IN), OPTIONAL               ::  Print_Timetable_Option
LOGICAL,            INTENT(IN), OPTIONAL               ::  Write_Timetable_Option
LOGICAL,            INTENT(IN), OPTIONAL               ::  Write_Sources_Option
LOGICAL,            INTENT(IN), OPTIONAL               ::  Print_Condition_Option
LOGICAL,            INTENT(IN), OPTIONAL               ::  Write_Condition_Option

LOGICAL,            INTENT(IN), OPTIONAL               ::  Write_FP_Diagnostics_Option
CHARACTER(LEN=10),  INTENT(IN), OPTIONAL               ::  Suffix_Flag_Option
CHARACTER(LEN=4),   INTENT(IN), OPTIONAL            ::  Suffix_Tail_Option
INTEGER,            INTENT(IN), OPTIONAL            ::  Frame_Option
INTEGER,            INTENT(IN), OPTIONAL            ::  Suffix_Param_Type_Option
CHARACTER(LEN=*),   INTENT(IN), OPTIONAL            ::  Suffix_Input_Option

CHARACTER(LEN=40)                                   ::  Suffix_Return


IF ( Verbose_Flag ) CALL Init_Message('Setting IO Parameters.')


IF ( Verbose_Flag ) THEN
    lPF_IO_Flags(iPF_IO_Verbose)        = .TRUE.
    lPF_IO_Flags(iPF_IO_Print_Setup)    = .TRUE.
    lPF_IO_Flags(iPF_IO_Print_Results)  = .TRUE.
    lPF_IO_Flags(iPF_IO_Print_Cond)     = .TRUE.
ELSE
    lPF_IO_Flags(iPF_IO_Verbose)        = .FALSE.
    lPF_IO_Flags(iPF_IO_Print_Setup)    = .FALSE.
    lPF_IO_Flags(iPF_IO_Print_Results)  = .FALSE.
    lPF_IO_Flags(iPF_IO_Print_Cond)     = .FALSE.
END IF



!
!   Setup
!
IF ( PRESENT(Print_Setup_Option) ) THEN
    IF ( Print_Setup_Option ) THEN
        lPF_IO_Flags(iPF_IO_Print_Setup) = .TRUE.
    ELSE
        lPF_IO_Flags(iPF_IO_Print_Setup) = .FALSE.
    END IF
END IF


IF ( PRESENT(Write_Setup_Option) ) THEN
    IF ( Write_Setup_Option ) THEN
        lPF_IO_Flags(iPF_IO_Write_Setup) = .TRUE.
    ELSE
        lPF_IO_Flags(iPF_IO_Write_Setup) = .FALSE.
    END IF
END IF



!
!   Results
!
IF ( PRESENT(Print_Results_Option) ) THEN
    IF ( Print_Results_Option ) THEN
        lPF_IO_Flags(iPF_IO_Print_Results) = .TRUE.
    ELSE
        lPF_IO_Flags(iPF_IO_Print_Results) = .FALSE.
    END IF
END IF


IF ( PRESENT(Write_Results_Option) ) THEN
    IF ( Write_Results_Option ) THEN
        lPF_IO_Flags(iPF_IO_Write_Results) = .TRUE.
    ELSE
        lPF_IO_Flags(iPF_IO_Write_Results) = .False.
    END IF
END IF





!
!   Timetable
!
IF ( PRESENT(Print_Timetable_Option) ) THEN
    IF ( Print_Timetable_Option ) THEN
        lPF_IO_Flags(iPF_IO_Print_Timetable) = .TRUE.
    ELSE
        lPF_IO_Flags(iPF_IO_Print_Timetable) = .FALSE.
    END IF
END IF



IF ( PRESENT(Write_Timetable_Option) ) THEN
    IF ( Write_Timetable_Option ) THEN
        lPF_IO_Flags(iPF_IO_Write_Timetable) = .TRUE.
    ELSE
        lPF_IO_Flags(iPF_IO_Write_Timetable) = .FALSE.
    END IF
END IF


!
!   Condition Number
!
IF ( PRESENT(Print_Condition_Option) ) THEN
    IF ( Print_Condition_Option ) THEN
        lPF_IO_Flags(iPF_IO_Print_Cond) = .TRUE.
    ELSE
        lPF_IO_Flags(iPF_IO_Print_Cond) = .FALSE.
    END IF
END IF



IF ( PRESENT(Write_Condition_Option) ) THEN
    IF ( Write_Condition_Option ) THEN
        lPF_IO_Flags(iPF_IO_Write_Cond) = .TRUE.
    ELSE
        lPF_IO_Flags(iPF_IO_Write_Cond) = .FALSE.
    END IF
END IF




!
!   Sources
!
IF ( PRESENT(Write_Sources_Option) ) THEN
    IF ( Write_Sources_Option ) THEN
        lPF_IO_Flags(iPF_IO_Write_Sources)  = .TRUE.
    ELSE
        lPF_IO_Flags(iPF_IO_Write_Sources)  = .FALSE.
    END IF
END IF




!
!   Fixed Point Diagnostics
!
IF ( PRESENT(Write_FP_Diagnostics_Option ) ) THEN
    IF ( Write_FP_Diagnostics_Option ) THEN
        lPF_IO_Flags(iPF_IO_Write_FP_Diagnostics) = .TRUE.
        FP_Diagnostics_Flag = .TRUE.
    ELSE
        lPF_IO_Flags(iPF_IO_Write_FP_Diagnostics) = .FALSE.
    END IF
END IF



!
!   Suffix
!
IF ( PRESENT(Suffix_Input_Option) ) THEN

    WRITE(File_Suffix,'(A)') TRIM(Suffix_Input_Option )

ELSE
    CALL Create_Suffix( Suffix_Return,              &
                        Suffix_Flag_Option,         &
                        Frame_Option,               &
                        Suffix_Param_Type_Option,   &
                        Suffix_Tail_Option          )
                    
    WRITE(File_Suffix,'(A)') TRIM(Suffix_Return)

END IF

lPF_Init_Flags(iPF_Init_IO_Params) = .TRUE.


END SUBROUTINE Init_IO_Params








 !+102+################################################!
!                                                       !
!          Init_Fixed_Point_Params                      !
!                                                       !
 !#####################################################!
SUBROUTINE Init_Fixed_Point_Params( Max_Iterations_Option,          &
                                    Convergence_Criteria_Option,    &
                                    Anderson_M_Option               )

INTEGER,                 INTENT(IN), OPTIONAL               ::  Max_Iterations_Option
REAL(idp),               INTENT(IN), OPTIONAL               ::  Convergence_Criteria_Option
INTEGER,                 INTENT(IN), OPTIONAL               ::  Anderson_M_Option


IF ( Verbose_Flag ) CALL Init_Message('Setting Fixed Point Parameters.')


IF ( PRESENT( Max_Iterations_Option ) ) THEN
    Max_Iterations = Max_Iterations_Option
ELSE
    Max_Iterations = Max_Iterations_Default
END IF


IF ( PRESENT( Convergence_Criteria_Option) ) THEN
    Convergence_Criteria = Convergence_Criteria_Option
ELSE
    Convergence_Criteria = Convergence_Criteria_Default
END IF


IF ( PRESENT(Anderson_M_Option) ) THEN
    FP_Anderson_M = Anderson_M_Option
ELSE
    FP_Anderson_M = FP_Anderson_M_Default
END IF


lPF_Init_Flags(iPF_Init_FP_Params) = .TRUE.


END SUBROUTINE Init_Fixed_Point_Params





 !+102+####################################################################!
!                                                                           !
!       Init_AMReX_Params                                                   !
!                                                                           !
 !#########################################################################!
SUBROUTINE Init_AMReX_Params(   nCells_Option,          &
                                Max_Level_Option,       &
                                Max_Grid_Size_Option    )

INTEGER,   DIMENSION(3), INTENT(IN), OPTIONAL       ::  nCells_Option
INTEGER,                 INTENT(IN), OPTIONAL       ::  Max_Level_Option
INTEGER,   DIMENSION(3), INTENT(IN), OPTIONAL       ::  Max_Grid_Size_Option


IF ( Verbose_Flag ) CALL Init_Message('Setting AMReX Parameters.')

IF ( PRESENT( Max_Level_Option ) ) THEN
    AMReX_Max_Level  = Max_Level_Option
ELSE
    AMReX_Max_Level  = 0
END IF
AMReX_Num_Levels = AMReX_Max_Level+1


IF ( PRESENT( Max_Level_Option ) ) THEN
    AMReX_Max_Grid_Size = Max_Grid_Size_Option
ELSE
    AMReX_Max_Grid_Size = 4
END IF



!IF ( PRESENT(AMReX_FEM_Refinement_Option) ) THEN
!    iFRL = AMReX_FEM_Refinement_Option
!ELSE
!    iFRL = AMReX_Max_Level
!END IF
!
!IF ( PRESENT(AMReX_Integral_Refinement_Option) ) THEN
!    iIRL = AMReX_Integral_Refinement_Option
!ELSE
!    iIRL = AMReX_Max_Level
!END IF

lPF_Init_AMReX_Flags(iPF_Init_AMReX_Params) = .TRUE.


END SUBROUTINE Init_AMReX_Params





 !+102+####################################################################!
!                                                                           !
!       Init_AMReX_Params                                                   !
!                                                                           !
 !#########################################################################!
SUBROUTINE Init_Mesh_Params(    NE_Option,          &
                                Inner_Edge_Option,  &
                                Outer_Edge_Option,  &
                                Radial_Units_Option )

INTEGER,    DIMENSION(3),   INTENT(IN),     OPTIONAL    ::  NE_Option
REAL(idp),  DIMENSION(3),   INTENT(IN),     OPTIONAL    ::  Inner_Edge_Option
REAL(idp),  DIMENSION(3),   INTENT(IN),     OPTIONAL    ::  Outer_Edge_Option
CHARACTER(LEN=2),           INTENT(IN),     OPTIONAL    ::  Radial_Units_Option

IF ( Verbose_Flag ) CALL Init_Message('Setting Mesh Parameters.')

IF ( PRESENT(NE_Option) ) THEN
    Num_R_Elements = NE_Option(1)
    Num_T_Elements = NE_Option(2)
    Num_P_Elements = NE_Option(3)
ELSE
    Num_R_Elements = 1
    Num_T_Elements = 1
    Num_P_Elements = 1
END IF


IF ( PRESENT(Inner_Edge_Option) ) THEN
    R_Inner = Inner_Edge_Option(1)
ELSE
    R_Inner = 0.0_idp
END IF

IF ( PRESENT(Outer_Edge_Option) ) THEN
    R_Outer = Outer_Edge_Option(1)
ELSE
    R_Outer = 0.0_idp
END IF



lPF_Init_Mesh_Flags(iPF_Init_Mesh_Params) = .TRUE.


END SUBROUTINE Init_Mesh_Params



 !+102+####################################################################!
!                                                                           !
!       Init_Quad_Params                                                   !
!                                                                           !
 !#########################################################################!
SUBROUTINE Init_Quad_Params(    NQ_Option   )

INTEGER,    DIMENSION(3),   INTENT(IN),     OPTIONAL    ::  NQ_Option

IF ( Verbose_Flag ) CALL Init_Message('Setting Quadrature Parameters.')

IF ( PRESENT( NQ_Option ) ) THEN
    Num_R_Quad_Points = NQ_Option(1)
    Num_T_Quad_Points = NQ_Option(2)
    Num_P_Quad_Points = NQ_Option(3)
ELSE
    Num_R_Quad_Points = 2*Degree + 2
    Num_T_Quad_Points = 1
    Num_P_Quad_Points = 2*L_Limit + 1
END IF
Num_TP_Quad_Points = Num_T_Quad_Points*Num_P_Quad_Points
Local_Quad_DOF     = Num_R_Quad_Points*Num_TP_Quad_Points


#ifdef POSEIDON_AMREX_FLAG
MF_Source_nComps   = 5*Local_Quad_DOF
#endif


lPF_Init_Quad_Flags(iPF_Init_Quad_Params) = .TRUE.


END SUBROUTINE Init_Quad_Params





 !+102+####################################################################!
!                                                                           !
!       Init_MPI_Params                                                   !
!                                                                           !
 !#########################################################################!
SUBROUTINE Init_MPI_Params(  )

IF ( Verbose_Flag ) CALL Init_Message('Setting MPI Parameters.')

CALL MPI_COMM_DUP(MPI_COMM_WORLD, Poseidon_Comm_World, ierr)
CALL MPI_COMM_SIZE(Poseidon_Comm_World, nProcs_Poseidon, ierr)
CALL MPI_COMM_RANK(Poseidon_Comm_World, myID_Poseidon, ierr)

lPF_Init_Flags(iPF_Init_MPI) = .TRUE.

END SUBROUTINE Init_MPI_Params




 !+102+####################################################################!
!                                                                           !
!       Set_Caller_Quadrature                                                   !
!                                                                           !
 !#########################################################################!
SUBROUTINE Set_Caller_Data( Source_NQ,          &
                            Source_xL,          &
                            Source_RQ_xlocs,    &
                            Source_TQ_xlocs,    &
                            Source_PQ_xlocs,    &
                            Source_Units,       &
                            Source_Radial_Boundary_Units )


INTEGER,    DIMENSION(3),               INTENT(IN)      ::  Source_NQ
REAL(idp),  DIMENSION(2),               INTENT(IN)      ::  Source_xL
REAL(idp),  DIMENSION(1:Source_NQ(1)),  INTENT(IN)      ::  Source_RQ_xlocs
REAL(idp),  DIMENSION(1:Source_NQ(2)),  INTENT(IN)      ::  Source_TQ_xlocs
REAL(idp),  DIMENSION(1:Source_NQ(3)),  INTENT(IN)      ::  Source_PQ_xlocs
CHARACTER(LEN=1),                       INTENT(IN)      ::  Source_Units
CHARACTER(LEN=2),                       INTENT(IN)      ::  Source_Radial_Boundary_Units




IF ( Verbose_Flag ) CALL Init_Message('Setting Caller Parameters.')

Caller_NQ = Source_NQ
Caller_xL = Source_xL

ALLOCATE( Caller_RQ_xlocs(1:Caller_NQ(1)) )
ALLOCATE( Caller_TQ_xlocs(1:Caller_NQ(2)) )
ALLOCATE( Caller_PQ_xlocs(1:Caller_NQ(3)) )

Caller_RQ_xlocs = Source_RQ_xlocs
Caller_TQ_xlocs = Source_TQ_xlocs
Caller_PQ_xlocs = Source_PQ_xlocs

Caller_Quad_DOF = Caller_NQ(1)*Caller_NQ(2)*Caller_NQ(3)

Caller_Set = .TRUE.

ALLOCATE( Translation_Matrix(1:Caller_Quad_DOF, 1:Local_Quad_DOF) )





IF ( Source_Radial_Boundary_Units == "cm" ) THEN
    Caller_R_Units = Centimeter
ELSE IF ( Source_Radial_Boundary_Units == " m" ) THEN
    Caller_R_Units = Meter
ELSE IF ( Source_Radial_Boundary_Units == "km" ) THEN
    Caller_R_Units = Kilometer
ELSE
    


    Caller_R_Units = Centimeter
END IF

lPF_Init_Flags(iPF_Init_Caller_Vars) = .TRUE.

END SUBROUTINE Set_Caller_Data




 !+101+################################################################################!
!                                                                                       !
!       Initialize_Mesh                                                                 !
!                                                                                       !
 !#####################################################################################!
SUBROUTINE Initialize_Mesh( NE_Option,                          &
                            Domain_Edge_Option,                 &
                            Radial_Boundary_Units_Option,       &
                            r_Option, t_Option, p_Option,       &
                            dr_Option, dt_Option, dp_Option     )

INTEGER,   DIMENSION(3), INTENT(IN), OPTIONAL               ::  NE_Option
CHARACTER(LEN=*),        INTENT(IN), OPTIONAL               ::  Radial_Boundary_Units_Option
REAL(idp), DIMENSION(2), INTENT(IN), OPTIONAL               ::  Domain_Edge_Option

REAL(idp), DIMENSION(:), INTENT(IN), OPTIONAL               ::  r_Option
REAL(idp), DIMENSION(:), INTENT(IN), OPTIONAL               ::  t_Option
REAL(idp), DIMENSION(:), INTENT(IN), OPTIONAL               ::  p_Option
REAL(idp), DIMENSION(:), INTENT(IN), OPTIONAL               ::  dr_Option
REAL(idp), DIMENSION(:), INTENT(IN), OPTIONAL               ::  dt_Option
REAL(idp), DIMENSION(:), INTENT(IN), OPTIONAL               ::  dp_Option


REAL(idp)                                                   ::  Caller_Units
LOGICAL,    DIMENSION(3)                                    ::  xSet
LOGICAL,    DIMENSION(3)                                    ::  dxSet
LOGICAL,    DIMENSION(3)                                    ::  Mesh_Set

IF ( Verbose_Flag ) CALL Init_Message('Initializing Mesh Variables.')



IF ( PRESENT( Radial_Boundary_Units_Option ) ) THEN
    IF ( Radial_Boundary_Units_Option == "cm" ) THEN
        Caller_Units = Centimeter
    ELSE IF ( Radial_Boundary_Units_Option == " m" ) THEN
        Caller_Units = Meter
    ELSE IF ( Radial_Boundary_Units_Option == "km" ) THEN
        Caller_Units = Kilometer
    ELSE
        Caller_Units = Centimeter
    END IF
ELSE
    Caller_Units = Centimeter
END IF


IF ( PRESENT( Domain_Edge_Option ) ) THEN
    R_Inner = Domain_Edge_Option(1)*Caller_Units
    R_Outer = Domain_Edge_Option(2)*Caller_Units
ELSE
    R_Inner = 0.0_idp*Caller_Units
    R_Outer = 1.0_idp*Caller_Units
END IF


CALL Init_Mesh_Params(  NE_Option,                      &
                        [R_Inner, 0.0_idp, 0.0_idp],    &
                        [R_Outer, pi, 2.0_idp*pi]       )



CALL Allocate_Mesh()


xSet  = .FALSE.
dxSet = .FALSE.

IF ( PRESENT( r_Option ) ) THEN
    rlocs = r_Option*Caller_Units
    xSet(1) = .TRUE.
END IF

IF ( PRESENT( dr_Option ) ) THEN
    drlocs = dr_Option*Caller_Units
    dxSet(1) = .TRUE.
END IF



IF ( PRESENT( t_Option ) ) THEN
    tlocs = t_Option
    xSet(2) = .TRUE.
END IF

IF ( PRESENT( dt_Option ) ) THEN
    dtlocs = dt_Option
    dxSet(2) = .TRUE.
END IF



IF ( PRESENT( p_Option ) ) THEN
    plocs = p_Option
    xSet(3) = .TRUE.
END IF

IF ( PRESENT( dp_Option ) ) THEN
    dplocs = dp_Option
    dxSet(3) = .TRUE.
END IF



Call Mesh_Check(xSet, dxSet, Mesh_Set)





IF (ALL(Mesh_Set) ) THEN
    lPF_Init_Mesh_Flags(iPF_Init_Mesh_Init) = .TRUE.
END IF


END SUBROUTINE Initialize_Mesh





SUBROUTINE Mesh_Check( xSet, dxSet, Mesh_Set )

LOGICAL,    DIMENSION(3),   INTENT(IN)      ::  xSet
LOGICAL,    DIMENSION(3),   INTENT(IN)      ::  dxSet
LOGICAL,    DIMENSION(3),   INTENT(OUT)     ::  Mesh_Set

INTEGER                                     ::  i


i = 1
IF ( (.NOT. xSet(i)) .AND. (.NOT. dxSet(i) ) ) THEN     ! Neither are Set
    drlocs = (R_Outer - R_Inner)/REAL(Num_R_Elements, KIND = idp)

    rlocs(0) = R_Inner
    DO i = 1,Num_R_Elements
        rlocs(i) = rlocs(i-1) + drlocs(i-1)
    END DO
ELSEIF ( (.NOT. xSet(i)) .AND. ( dxSet(i)) ) THEN       ! Only dx is set
    rlocs(0) = R_Inner
    DO i = 1,Num_R_Elements
        rlocs(i) = rlocs(i-1) + drlocs(i-1)
    END DO
ELSEIF ( ( xSet(i)) .AND. (.NOT. dxSet(i) ) ) THEN      ! Only x is set
    DO i = 0,Num_R_Elements-1
        drlocs(i) = rlocs(i+1)-rlocs(i)
    END DO
END IF
Mesh_Set(1) = .TRUE.




i = 2
IF ( (.NOT. xSet(i)) .AND. (.NOT. dxSet(i) ) ) THEN     ! Neither are Set
    dtlocs = ( pi )/REAL(Num_T_Elements, KIND = idp)

    tlocs(0) = 0.0_idp
    DO i = 1,Num_T_Elements
        tlocs(i) = tlocs(i-1) + dtlocs(i-1)
    END DO
ELSEIF ( (.NOT. xSet(i)) .AND. ( dxSet(i)) ) THEN       ! Only dx is set
    tlocs(0) = 0.0_idp
    DO i = 1,Num_T_Elements
        tlocs(i) = tlocs(i-1) + dtlocs(i-1)
    END DO
ELSEIF ( ( xSet(i)) .AND. (.NOT. dxSet(i) ) ) THEN      ! Only x is set
    DO i = 0,Num_T_Elements-1
        dtlocs(i) = tlocs(i+1)-tlocs(i)
    END DO
END IF
Mesh_Set(2) = .TRUE.


i = 3
IF ( (.NOT. xSet(i)) .AND. (.NOT. dxSet(i) ) ) THEN     ! Neither are Set
    dplocs = ( 2.0_idp*pi )/REAL(Num_P_Elements, KIND = idp)

    plocs(0) = 0.0_idp
    DO i = 1,Num_P_Elements
        plocs(i) = plocs(i-1) + dplocs(i-1)
    END DO
ELSEIF ( (.NOT. xSet(i)) .AND. ( dxSet(i)) ) THEN       ! Only dx is set
    plocs(0) = 0.0_idp
    DO i = 1,Num_P_Elements
        plocs(i) = plocs(i-1) + dplocs(i-1)
    END DO
ELSEIF ( ( xSet(i)) .AND. (.NOT. dxSet(i) ) ) THEN      ! Only x is set
    DO i = 0,Num_P_Elements-1
        dplocs(i) = plocs(i+1)-plocs(i)
    END DO
END IF
Mesh_Set(3) = .TRUE.


END SUBROUTINE Mesh_Check






 !+102+####################################################################!
!                                                                           !
!       Set_Method_Flags                                                   !
!                                                                           !
 !#########################################################################!
SUBROUTINE Set_Method_Flags(Newtonian_Mode_Option,      &
                            CFA_Mode_Option,            &
                            XCFC_Mode_Option            )

LOGICAL,                 INTENT(IN), OPTIONAL       ::  Newtonian_Mode_Option
LOGICAL,                 INTENT(IN), OPTIONAL       ::  CFA_Mode_Option
LOGICAL,                 INTENT(IN), OPTIONAL       ::  XCFC_Mode_Option

LOGICAL,    DIMENSION(3)                            ::  Roll_Call
INTEGER                                             ::  Attendence

LOGICAL,    DIMENSION(3)                            ::  Status
INTEGER                                             ::  Total


Roll_Call = (/  PRESENT(Newtonian_Mode_Option), &
                PRESENT(CFA_Mode_Option),       &
                PRESENT(XCFC_Mode_Option)       /)
Attendence = COUNT(Roll_Call)




IF ( Attendence == 0 ) THEN ! Nothing set => Go to Default => XCFC

    iPF_Core_Flags(iPF_Core_Method_Mode) = iPF_Core_Method_XCFC


ELSE



    Status = .FALSE.
    IF ( Roll_Call(1) ) THEN
        Status(1) = Newtonian_Mode_Option
    END IF

    IF ( Roll_Call(2) ) THEN
        Status(2) = CFA_Mode_Option
    END IF

    IF ( Roll_Call(3) ) THEN
        Status(3) = XCFC_Mode_Option
    END IF



    Total = COUNT(Status)



    IF ( Total .GE. 2 ) THEN
        iPF_Core_Flags(iPF_Core_Method_Mode) = iPF_Core_Method_Too_Many
    ELSE

        IF ( Status(1) ) THEN
            iPF_Core_Flags(iPF_Core_Method_Mode) = iPF_Core_Method_Newtonian
        END IF

        IF ( Status(2) ) THEN
            iPF_Core_Flags(iPF_Core_Method_Mode) = iPF_Core_Method_CFA
        END IF

        IF ( Status(3) ) THEN
            iPF_Core_Flags(iPF_Core_Method_Mode) = iPF_Core_Method_XCFC
        END IF

    END IF ! Total .GE. 2



END IF ! Attendence == 0




END SUBROUTINE Set_Method_Flags




END MODULE Initialization_Subroutines
