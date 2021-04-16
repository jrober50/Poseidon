   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
PROGRAM Yahil_Test                                                                  !##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!

USE Poseidon_Kinds_Module, &
                ONLY :  idp

USE Units_Module, &
                ONLY :  C_Square,       &
                        Set_Units,      &
                        Centimeter

USE Initialization_Poseidon, &
                ONLY :  Initialize_Poseidon

USE Variables_IO, &
                ONLY :  Write_Results_Flag

USE Variables_MPI, &
                ONLY :  ierr

USE Variables_FP, &
                ONLY :  FP_Coeff_Vector,            &
                        FP_Coeff_Vector_Beta,       &
                        FP_Source_Vector_Beta

USE Variables_Functions, &
                ONLY :  Potential_Solution

USE Functions_Mesh, &
                ONLY :  Create_3D_Mesh

USE Functions_Quadrature, &
                ONLY :  Initialize_LG_Quadrature_Locations

USE Functions_Mapping, &
                ONLY :  Map_From_X_Space

USE Poseidon_IO_Module, &
                ONLY :  Open_Run_Report_File,       &
                        Output_Final_Results

USE IO_Print_Results, &
                ONLY :  Print_Results

USE IO_Convergence_Output, &
                ONLY :  Output_Convergence_Data

USE Poseidon_Main_Module, &
                ONLY :  Poseidon_CFA_Set_Uniform_Boundary_Conditions,       &
                        Poseidon_Run,                                       &
                        Poseidon_Close

USE FP_Method_Module, &
                ONLY :  Solve_FP_System,        &
                        Solve_FP_System_Beta

USE FP_Initial_Guess_Module, &
                ONLY :  Init_FP_Guess_Flat,     &
                        Init_FP_Guess_Informed

USE SelfSimilar_Module, &
                ONLY :  Initialize_Yahil_Sources

USE Source_Input_Module, &
                ONLY :  Poseidon_Input_Sources

USE Poseidon_BC_Module, &
                ONLY :  Calc_Shift_BC_1D




USE MPI


IMPLICIT NONE

!                                       !
!   Poseidon Initialization Variables   !
!                                       !
INTEGER                                                 ::  Test_Number

INTEGER                                                 ::  FEM_Degree_Input
INTEGER                                                 ::  L_Limit_Input

INTEGER                                                 ::  Dimension_Input

INTEGER                                                 ::  Mesh_Type
INTEGER                                                 ::  Solver_Type

CHARACTER(LEN = 1)                                      ::  Units_Input

INTEGER, DIMENSION(3)                                   ::  NE
INTEGER, DIMENSION(3)                                   ::  NQ
REAL(idp), DIMENSION(2)                                 ::  Domain_Edge
INTEGER                                                 ::  Num_DOF


REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  dx_c, dy_c, dz_c
REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  x_e, y_e, z_e
REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  x_c, y_c, z_c

LOGICAL                                                 ::  Verbose

INTEGER,   DIMENSION(5)                                 ::  CFA_EQs

REAL(KIND = idp)                                        ::  Outer_Potential,    &
                                                            Inner_Potential,    &
                                                            Pot_to_Alphapsi,    &
                                                            Pot_to_Psi


REAL(idp)                                               ::  Shift_Vector_BC
CHARACTER(LEN=1), DIMENSION(1:5)                        ::  INNER_BC_TYPES, OUTER_BC_TYPES
REAL(idp), DIMENSION(1:5)                               ::  INNER_BC_VALUES, OUTER_BC_VALUES

CHARACTER(LEN=10)                                       ::  Suffix_Input

REAL(idp)                                               ::  Yahil_Time
REAL(idp)                                               ::  Kappa
REAL(idp)                                               ::  Gamma
REAL(idp)                                               ::  Ecc

REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  Input_R_Quad
REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  Input_T_Quad
REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  Input_P_Quad
REAL(idp)                                               ::  Left_Limit
REAL(idp)                                               ::  Right_Limit

REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              ::  Local_E
REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              ::  Local_S
REAL(idp), DIMENSION(:,:,:,:,:), ALLOCATABLE            ::  Local_Si

INTEGER                                                 ::  myID


CALL MPI_INIT(ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

!!                                       !!
!!   Poseidon Initialization Variables   !!
!!                                       !!

Units_Input         = "G"
Solver_Type         = 2
Test_Number         = 3

Yahil_Time          = 150.0_idp     ! ms
!Yahil_Time          = 0.01_idp     ! ms
Kappa               = 9.5390E+14
Gamma               = 1.30_idp
Ecc                 = 0.0_idp

FEM_Degree_Input    = 1
L_Limit_Input       = 0

Dimension_Input     = 3

Mesh_Type           = 2
Domain_Edge(1)      = 0.0_idp   ! Inner Radius (cm)
Domain_Edge(2)      = 1.0E+10   ! Outer Radius (cm)
NE(1)               = 64        ! Number of Radial Elements
NE(2)               = 1         ! Number of Theta Elements
NE(3)               = 1         ! Number of Phi Elements

NQ(1)               = 10        ! Number of Radial Quadrature Points
NQ(2)               = 10        ! Number of Theta Quadrature Points
NQ(3)               = 10        ! Number of Phi Quadrature Points

Verbose             = .FALSE.    !
Suffix_Input        = "Params"

CFA_Eqs = [1, 0, 0, 0, 0]


Num_DOF = NQ(1)*NQ(2)*NQ(3)

ALLOCATE( Input_R_Quad(1:NQ(1)) )
ALLOCATE( Input_T_Quad(1:NQ(2)) )
ALLOCATE( Input_P_Quad(1:NQ(3)) )

ALLOCATE( x_e(0:NE(1)), y_e(0:NE(2)), z_e(0:NE(3)) )
ALLOCATE( x_c(1:NE(1)), y_c(1:NE(2)), z_c(1:NE(3)) )
ALLOCATE( dx_c(1:NE(1)) )
ALLOCATE( dy_c(1:NE(2)) )
ALLOCATE( dz_c(1:NE(3)) )

ALLOCATE( Local_E(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1 )       )
ALLOCATE( Local_S(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1 )       )
ALLOCATE( Local_Si(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1, 1:3)  )


Input_R_Quad = Initialize_LG_Quadrature_Locations(NQ(1))
Input_T_Quad = Initialize_LG_Quadrature_Locations(NQ(2))
Input_P_Quad = Initialize_LG_Quadrature_Locations(NQ(3))

Left_Limit  = -0.50_idp
Right_Limit = +0.50_idp

Input_R_Quad = Map_From_X_Space(Left_Limit, Right_Limit, Input_R_Quad)
Input_T_Quad = Map_From_X_Space(Left_Limit, Right_Limit, Input_T_Quad)
Input_P_Quad = Map_From_X_Space(Left_Limit, Right_Limit, Input_P_Quad)




CALL Open_Run_Report_File()


CALL Set_Units(Units_Input)


CALL Create_3D_Mesh( Mesh_Type,         &
                     Domain_Edge(1),    &
                     Domain_Edge(2)*Centimeter,    &
                     NE(1),             &
                     NE(2),             &
                     NE(3),             &
                     x_e, x_c, dx_c,    &
                     y_e, y_c, dy_c,    &
                     z_e, z_c, dz_c     )

!PRINT*,Domain_Edge
!PRINT*,x_e
!
!STOP


CALL Initialize_Poseidon &
    (   Dimensions_Option       = Dimension_Input,  &
        FEM_Degree_Option       = FEM_Degree_Input, &
        L_Limit_Option          = L_Limit_Input,    &
        Units_Option            = Units_Input,      &
        Domain_Edge_Option      = Domain_Edge*Centimeter,      &
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


CALL Initialize_Yahil_Sources( Yahil_Time, Kappa, Gamma, Ecc,       &
                                NQ, Input_R_Quad, Input_T_Quad,     &
                                NE(1), NE(2), NE(3),                &
                                dx_c, x_e, y_e,                     &
                                Local_E, Local_S, Local_Si )

!PRINT*,"Local_E"
!PRINT*,Local_E
!PRINT*,"Local_S"
!PRINT*,Local_S
!PRINT*,"Local_Si"
!PRINT*,Local_Si

CALL Poseidon_Input_Sources(    myID, myID, myID,                           &
                                Local_E, Local_S, Local_Si,                 &
                                NE(1), NE(2), NE(3),                        &
                                NQ(1), NQ(2), NQ(3),                        &
                                Input_R_Quad, Input_T_Quad, Input_P_Quad,   &
                                Left_Limit, Right_Limit                     )





CALL Calc_Shift_BC_1D( Shift_Vector_BC,             &
                       NQ(1), NQ(2), NQ(3),         &
                       NE(1), NE(2), NE(3), 3,      &
                       Local_Si, x_e                )
!PRINT*,"Error with Shift Vector BC. Setting to 1.0"
!Shift_Vector_BC = 1.0_idp

! Calculate the Dirichlet Outer Boundary Value for the Lapse Function and Conformal Factor
Inner_Potential = Potential_Solution(x_e(0), 0.0_idp, 0.0_idp)
Outer_Potential = Potential_Solution(x_e(NE(1)), 0.0_idp, 0.0_idp)




!!PRINT*,"ALPHAPSI BC set to one"
!!Pot_to_AlphaPsi = 1.0_idp
Pot_to_Alphapsi = 1.0_idp + 0.5_idp*Outer_Potential/C_Square
Pot_to_Psi      = 1.0_idp - 0.5_idp*Outer_Potential/C_Square


INNER_BC_TYPES = (/"N", "N","N","N","N"/)
OUTER_BC_TYPES = (/"D", "D","D","D","D"/)


PRINT*,"BCs",Pot_to_Psi, Pot_to_Alphapsi, Shift_Vector_BC
INNER_BC_VALUES = (/0.0_idp, 0.0_idp, 0.0_idp,         0.0_idp, 0.0_idp /)
OUTER_BC_VALUES = (/Pot_to_Psi, Pot_to_Alphapsi, Shift_Vector_BC, 0.0_idp, 0.0_idp /)


CALL Poseidon_CFA_Set_Uniform_Boundary_Conditions("I", INNER_BC_TYPES, INNER_BC_VALUES)
CALL Poseidon_CFA_Set_Uniform_Boundary_Conditions("O", OUTER_BC_TYPES, OUTER_BC_VALUES)




CALL Init_FP_Guess_flat()

!CALL Init_FP_Guess_Informed(NE(1), NE(2), NE(3),    &
!                            NQ(1), NQ(2), NQ(3),    &
!                            x_e,                    &
!                            Local_Si                )





Call Poseidon_Run()






IF (Verbose .EQV. .TRUE. ) THEN
    CALL Print_Results()
END IF


Write_Results_Flag = 1
IF ( Write_Results_Flag == 1 ) THEN
    CALL Output_Convergence_Data()
    CALL Output_Final_Results()
END IF




CALL Poseidon_Close()


DEALLOCATE( x_e, y_e, z_e )
DEALLOCATE( x_c, y_c, z_c )
DEALLOCATE( dx_c, dy_c, dz_c )

DEALLOCATE( Local_E, Local_S, Local_Si )



END PROGRAM Yahil_Test

