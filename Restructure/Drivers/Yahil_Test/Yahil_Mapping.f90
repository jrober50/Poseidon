   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
PROGRAM Beta_Mapping                                                                !##!
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

USE Variables_IO, &
                ONLY :  Write_Results_Flag

USE Variables_MPI, &
                ONLY :  ierr

USE Variables_FP, &
                ONLY :  FP_Coeff_Vector,        &
                        FP_Coeff_Vector_Beta,   &
                        FP_Source_Vector_Beta

USE Variables_Functions, &
                ONLY :  Potential_Solution

USE Functions_FP, &
                ONLY :  Laplace_Test_Sol

USE Functions_Mesh, &
                ONLY :  Create_3D_Mesh

USE Poseidon_IO_Module, &
                ONLY :  Open_Run_Report_File,       &
                        Output_Run_Report,        &
                        Output_Final_Results

USE Poseidon_Main_Module, &
                ONLY :  Poseidon_CFA_Set_Uniform_Boundary_Conditions,   &
                        Poseidon_Close

USE FP_Method_Module, &
                ONLY :  Solve_FP_System,        &
                        Solve_FP_System_Beta

USE IO_Print_Results, &
                ONLY :  Print_Results

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

CHARACTER(LEN=10)                                       ::  Suffix_Input

INTEGER                                                 ::  RE_Index, RE_Index_Min, RE_Index_Max
INTEGER                                                 ::  Degree_Min, Degree_Max
INTEGER                                                 ::  L_Limit_Min, L_Limit_Max

INTEGER                                                 ::  Num_RE
INTEGER, DIMENSION(:), ALLOCATABLE                      ::  RE_Table


CALL MPI_INIT(ierr)

!!                                       !!
!!   Poseidon Initialization Variables   !!
!!                                       !!

Units_Input         = "U"
Solver_Type         = 2
Suffix_Input        = "Params"


Dimension_Input     = 3

Mesh_Type           = 1
Domain_Edge(1)      = 1.0_idp   ! Inner Radius
Domain_Edge(2)      = 2.0_idp   ! Outer Radius

RE_Index_Min        = 1
RE_Index_Max        = 8

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




DO RE_Index = RE_Index_Min, RE_Index_Max
    DO FEM_Degree_Input = Degree_Min, Degree_Max
        DO L_Limit_Input = L_Limit_Min, L_Limit_Max

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
                (   Dimensions_Option       = Dimension_Input,  &
                    FEM_Degree_Option       = FEM_Degree_Input, &
                    L_Limit_Option          = L_Limit_Input,    &
                    Units_Option            = Units_Input,      &
                    Domain_Edge_Option      = Domain_Edge,      &
                    NE_Option               = NE,               &
                    NQ_Option               = NQ,               &
                    r_Option                = x_e,              &
                    t_Option                = y_e,              &
                    p_Option                = z_e,              &
            !        dr_Option               = dx_c,             &
            !        dt_Option               = dy_c,             &
            !        dp_Option               = dz_c              &
                    Suffix_Flag_Option       = Suffix_Input,    &
                    Solver_Type_Option       = Solver_Type,     &
                    CFA_Eq_Flags_Option      = CFA_Eqs,         &
                    Verbose_Option           = Verbose          )






            CALL Poseidon_CFA_Set_Uniform_Boundary_Conditions("I", INNER_BC_TYPES, INNER_BC_VALUES)
            CALL Poseidon_CFA_Set_Uniform_Boundary_Conditions("O", OUTER_BC_TYPES, OUTER_BC_VALUES)


            ! For the Laplace Test the source vector is zero
            FP_Source_Vector_Beta = 0.0_idp


            ! For this test, the initial guess is zero
            FP_Coeff_Vector       = 0.0_idp
            FP_Coeff_Vector_Beta  = 0.0_idp


            !Call Solve_FP_System()
            Call Solve_FP_System_Beta()

            IF (Verbose .EQV. .TRUE. ) THEN
                CALL Print_Results()
            END IF


            Write_Results_Flag = 1
            IF ( Write_Results_Flag == 1 ) THEN
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

END PROGRAM Beta_Mapping
