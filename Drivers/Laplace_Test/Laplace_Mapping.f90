   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
PROGRAM Laplace_Mapping                                                             !##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!

USE Poseidon_Kinds_Module, &
                ONLY :  idp

USE Initialization_Poseidon, &
                ONLY :  Initialize_Poseidon


USE Variables_MPI, &
                ONLY :  ierr

USE Variables_Functions, &
                ONLY :  Potential_Solution

USE Functions_Mesh, &
                ONLY :  Create_3D_Mesh

USE Poseidon_IO_Module, &
                ONLY :  Open_Run_Report_File,       &
                        Output_Run_Report,        &
                        Output_Final_Results

USE Poseidon_Main_Module, &
                ONLY :  Poseidon_CFA_Set_Uniform_Boundary_Conditions,   &
                        Poseidon_Close

USE FP_System_Solvers_Module, &
                ONLY :  Solve_FP_System

USE MPI


IMPLICIT NONE

!                                       !
!   Poseidon Initialization Variables   !
!                                       !


INTEGER                                                 ::  FEM_Degree_Input
INTEGER                                                 ::  L_Limit_Input

INTEGER                                                 ::  Dimension_Input

INTEGER                                                 ::  Mesh_Type
INTEGER                                                 ::  Solver_Type

CHARACTER(LEN = 1)                                      ::  Units_Input

INTEGER, DIMENSION(3)                                   ::  NE
INTEGER, DIMENSION(3)                                   ::  NQ
REAL(idp), DIMENSION(2)                                 ::  Domain_Edge

REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  dx_c, dy_c, dz_c
REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  x_e, y_e, z_e
REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  x_c, y_c, z_c

LOGICAL                                                 ::  Verbose

INTEGER,   DIMENSION(5)                                 ::  CFA_EQs
REAL(idp)                                               ::  Shift_Vector_BC
CHARACTER(LEN=1), DIMENSION(1:5)                        ::  INNER_BC_TYPES, OUTER_BC_TYPES
REAL(idp), DIMENSION(1:5)                               ::  INNER_BC_VALUES, OUTER_BC_VALUES

CHARACTER(LEN=10)                                       ::  Suffix_Flag_Input

INTEGER                                                 ::  RE_Index, RE_Index_Min, RE_Index_Max
INTEGER                                                 ::  Degree, Degree_Min, Degree_Max
INTEGER                                                 ::  L_Limit, L_Limit_Min, L_Limit_Max

INTEGER                                                 ::  Num_RE
INTEGER, DIMENSION(:), ALLOCATABLE                      ::  RE_Table

LOGICAL                                                 ::  Print_Results_Flag


CALL MPI_INIT(ierr)

!!                                       !!
!!   Poseidon Initialization Variables   !!
!!                                       !!

Units_Input         = "U"
Solver_Type         = 2
Suffix_Flag_Input   = "Params"


Dimension_Input     = 3

Mesh_Type           = 1
Domain_Edge(1)      = 1.0_idp   ! Inner Radius
Domain_Edge(2)      = 2.0_idp   ! Outer Radius

RE_Index_Min        = 2
RE_Index_Max        = 2

Num_RE = RE_Index_Max - RE_Index_Min + 1
ALLOCATE( RE_Table(1:Num_RE) )
RE_Table            = [8, 16, 32, 64, 128, 256, 512, 1024]

Degree_Min          = 1
Degree_Max          = 1

L_Limit_Min         = 1
L_Limit_Max         = 1

Verbose             = .TRUE.    !

CFA_Eqs = [0, 0, 1, 0, 0]

INNER_BC_TYPES = (/"N", "N","N","N","N"/)
OUTER_BC_TYPES = (/"D", "D","D","D","D"/)

Shift_Vector_BC = -1.0E2_idp
OUTER_BC_VALUES = (/0.0_idp, 0.0_idp, Shift_Vector_BC, 0.0_idp, 0.0_idp /)

NQ(1) = 10        ! Number of Radial Quadrature Points
NQ(2) = 5         ! Number of Theta Quadrature Points
NQ(3) = 5         ! Number of Phi Quadrature Points

Print_Results_Flag  = .TRUE.


DO RE_Index = RE_Index_Min, RE_Index_Max
    DO Degree = Degree_Min, Degree_Max
        DO L_Limit = L_Limit_Min, L_Limit_Max

            NE(1) = RE_Table(RE_Index)      ! Number of Radial Elements
            NE(2) = 1                       ! Number of Theta Elements
            NE(3) = 1                       ! Number of Phi Elements



            ALLOCATE( x_e(0:NE(1)), y_e(0:NE(2)), z_e(0:NE(3)) )
            ALLOCATE( x_c(1:NE(1)), y_c(1:NE(2)), z_c(1:NE(3)) )
            ALLOCATE( dx_c(1:NE(1)), dy_c(1:NE(2)), dz_c(1:NE(3)) )

            CALL Open_Run_Report_File()

            CALL Create_3D_Mesh( Mesh_Type,         &
                                 Domain_Edge(1),    &
                                 Domain_Edge(2),    &
                                 NE(1),             &
                                 NE(2),             &
                                 NE(3),             &
                                 x_e, x_c, dx_c,    &
                                 y_e, y_c, dy_c,    &
                                 z_e, z_c, dz_c     )



            CALL Initialize_Poseidon &
                (   Dimensions_Option       = Dimension_Input,      &
                    FEM_Degree_Option       = Degree,               &
                    L_Limit_Option          = L_Limit,              &
                    Units_Option            = Units_Input,          &
                    Domain_Edge_Option      = Domain_Edge,          &
                    NE_Option               = NE,                   &
                    NQ_Option               = NQ,                   &
                    r_Option                = x_e,                  &
                    t_Option                = y_e,                  &
                    p_Option                = z_e,                  &
            !        dr_Option               = dx_c,                 &
            !        dt_Option               = dy_c,                 &
            !        dp_Option               = dz_c                  &
                    Method_Flag_Option       = Solver_Type,         &
                    CFA_Eq_Flags_Option      = CFA_Eqs,             &
                    Suffix_Flag_Option       = Suffix_Flag_Input,   &
                    Verbose_Option           = Verbose,                       &
                    WriteAll_Option             = .FALSE.,                       &
                    Print_Setup_Option          = .TRUE.,                        &
                    Write_Setup_Option          = .TRUE.,                       &
                    Print_Results_Option        = Print_Results_Flag,            &
                    Write_Results_Option        = .TRUE.,                        &
                    Print_Timetable_Option      = .TRUE.,                       &
                    Write_Timetable_Option      = .TRUE.,                       &
                    Write_Sources_Option        = .TRUE.     )





            CALL Poseidon_CFA_Set_Uniform_Boundary_Conditions("I", INNER_BC_TYPES, INNER_BC_VALUES)
            CALL Poseidon_CFA_Set_Uniform_Boundary_Conditions("O", OUTER_BC_TYPES, OUTER_BC_VALUES)


            !IF ( .TRUE. )
            !    CALL Initialize_Flat_Space_Guess_Values
            !ELSE
            !    CALL Initialize_Calculated_Guess_Values
            !END IF



            Call Solve_FP_System()



            IF ( .TRUE. ) THEN
                CALL Output_Run_Report()
                CALL Output_Final_Results()
            END IF


            CALL Poseidon_Close()

            DEALLOCATE( x_e, y_e, z_e )
            DEALLOCATE( x_c, y_c, z_c )
            DEALLOCATE( dx_c, dy_c, dz_c )


        END DO ! L_Limit
    END DO ! Degree
END DO ! RE_Index


WRITE(*,'(//A18//)')"DING! You're Done!"

END PROGRAM Laplace_Mapping
