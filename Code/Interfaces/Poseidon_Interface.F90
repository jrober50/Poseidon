   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_Interface                                                           !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains the interfacing subroutine that will handle action between         !##!
!##!    Poseidon and CHIMERA.                                                       !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!



USE Poseidon_Kinds_Module, &
                    ONLY : idp

USE Poseidon_Numbers_Module, &
                    ONLY : pi
USE Units_Module, &
            ONLY :  Set_Units,          &
                    C_Square


USE DRIVER_Parameters,  &
            ONLY :  DRIVER_R_ELEMS,                                     &
                    DRIVER_T_ELEMS,                                     &
                    DRIVER_P_ELEMS,                                     &
                    DRIVER_R_INPUT_NODES,                               &
                    DRIVER_T_INPUT_NODES,                               &
                    DRIVER_P_INPUT_NODES,                               &
                    nPROCS,                                             &
                    myID,                                               &
                    myID_Theta,                                         &
                    myID_Phi,                                           &
                    DRIVER_FIRST_GUESS_FLAG,                            &
                    DRIVER_SUBSEQUENT_GUESS_FLAG,                       &
                    DRIVER_FRAME

USE Variables_Functions, &
            ONLY :  Potential_Solution,             &
                    Shift_Solution



USE Poseidon_Parameters, &
            ONLY :  DOMAIN_DIM,             &
                    DEGREE,                 &
                    L_LIMIT,                &
                    Ratio_T_BNDLperBLCK,    &
                    Ratio_P_BNDLperBLCK,    &
                    Ratio_BNDLperBLCK


USE Variables_MPI, &
            ONLY :  Num_Block_Phi_Columns,  &
                    Num_Block_Theta_Rows,   &
                    ierr,                   &
                    POSEIDON_COMM_WORLD


USE Variables_Mesh, &
            ONLY :  NUM_R_ELEMENTS,         &
                    NUM_T_ELEMENTS,         &
                    NUM_P_ELEMENTS,         &
                    rlocs,                  &
                    tlocs,                  &
                    plocs


USE Poseidon_Main_Module, &
            ONLY :  Poseidon_Initialize,                                &
                    Poseidon_Run,                                       &
                    Poseidon_Close,                                     &
                    Poseidon_Set_Mesh,                                  &
                    Poseidon_CFA_Set_Uniform_Boundary_Conditions

USE Poseidon_Initialization_Module, &
            ONLY :  Poseidon_Initialize_From_File

USE FP_Initialization_Module, &
            ONLY :  Poseidon_FP_Init_From_File


USE Functions_Math, &
            ONLY :  Lagrange_Poly,                              &
                    Spherical_Harmonic

USE Functions_Mapping, &
            ONLY :  Map_To_X_Space,             &
                    Map_From_X_Space

USE Functions_Quadrature, &
            ONLY :  Initialize_LGL_Quadrature,                          &
                    Initialize_LGL_Quadrature_Locations,                &
                    Initialize_LG_Quadrature_Locations


USE Poseidon_IO_Module, &
            ONLY :  Clock_In,                                           &
                    OUTPUT_ITER_TIMETABLE,                              &
                    Write_Shift_1D

USE Poseidon_Parameter_Read_Module, &
            ONLY :  WRITE_CFA_COEFFICIENTS,                             &
                    READ_CFA_COEFFICIENTS 


USE Initial_Guess_Module, &
            ONLY :  Initialize_Flat_Space_Guess_Values,                 &
                    Initialize_Calculated_Guess_Values,                 &
                    Load_Initial_Guess_From_File,                    &
                    Calc_Shift_1D


USE Functions_Mesh, &
            ONLY :  Create_Uniform_1D_Mesh,                             &
                    Create_Logarithmic_1D_Mesh,                         &
                    Create_Split_1D_Mesh

USE Poseidon_BC_Module, &
            ONLY :  Calc_Shift_BC_1D


USE Source_Input_Module, &
            ONLY :   Poseidon_Input_Sources

USE Poseidon_Internal_Communication_Module, &
            ONLY :  Poseidon_Distribute_Solution





use mpi








IMPLICIT NONE

CONTAINS









SUBROUTINE Poseidon_CFA_3D(mode, modeb, units,                                          &
                            Solver_Mode, CFA_Eqs_Flag_Vector,                           &
                            imin, imax, nx, ij_ray_dim, ik_ray_dim,                     &
                            ny, nz, x_e, x_c, dx_c, y_e, y_c, dy_c, z_e, dz_c,          &
                            Num_Input_Nodes,                                            &
                            Input_R_Quad, Input_T_Quad, Input_P_Quad,                   &
                            Left_Limit, Right_Limit,                        &
                            Local_E, Local_S, Local_Si                                  )
!-----------------------------------------------------------------------
!
!
!     Potential values are calculated at zone interfaces and
!     stored in array POT(1:imax+1,1:jmax+1).
!
!     POT(i,1):      0     1     2     3  ...
!     zone:          |--1--|--2--|--3--|  ...
!
!
!     MODE = 1  ==>  Initialized the Poseidon. 
!                    Must be used in the first call to this subroutine.
!                    Else  MODE = 0
!
!     INPUT            rho_c(nx,ij_ray_dim,ik_ray_dim)    density [g cm^{-3}]
!                      x_e(nx+1)              radius (in cm) of left zone interface
!                      y_e(ny+1)              angle (in radians) of left zone interface
!                      z_e(nz+1)              angle (in radians) of left zone interface
!
!     OUTPUT           grav_x_c(ii,ij_ray_dim,ik_ray_dim)   zone-centered acceleration in x direction [dynes g^{-1} = cm s^{-2}]
!                      grav_y_c(ii,ij_ray_dim,ik_ray_dim)   zone-centered acceleration in y direction [dynes g^{-1} = cm s^{-2}]
!                      grav_z_c(ii,ij_ray_dim,ik_ray_dim)   zone-centered acceleration in z direction [dynes g^{-1} = cm s^{-2}]
!                      grav_pot_c(ii,ij_ray_dim,ik_ray_dim) zone-centered potential
!                      grav_x_e(ii,ij_ray_dim,ik_ray_dim)   zone-edged acceleration in x direction [dynes g^{-1} = cm s^{-2}]
!                      grav_y_e(ii,ij_ray_dim,ik_ray_dim)   zone-edged acceleration in y direction [dynes g^{-1} = cm s^{-2}]
!                      grav_z_e(ii,ij_ray_dim,ik_ray_dim)   zone-edged acceleration in z direction [dynes g^{-1} = cm s^{-2}]
!                      grav_pot_e(ii,ij_ray_dim,ik_ray_dim) zone-edged potential
!
!    Input arguments:
!  imin             : lower x-array index
!  imax             : upper x-array index
!  nx               : x-array extent
!  ij_ray_dim       : number of y-zones on a processor before swapping
!  ik_ray_dim       : number of z-zones on a processor before swapping
!  ny               : y-array extent
!  nz               : z-array extent
!  x_e              : x grid zone faces
!  x_c              ! x grid zone centers
!  dx_c             : x_e(i+1) - x_e(i)
!  y_e              : y grid zone left interfaces
!  y_c              ! y grid zone centers
!  dy_c             : y_e(i+1) - y_e(i)
!  z_e              : z grid zone left interfaces
!  dz_c             : z_e(i+1) - z_e(i)
!
!
!    Output arguments:
!
!-----------------------------------------------------------------------

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(IN)                                             :: mode
INTEGER, INTENT(IN)                                             :: modeb
INTEGER, INTENT(IN)                                             :: Solver_Mode
INTEGER, DIMENSION(1:5), INTENT(IN)                             :: CFA_EQs_Flag_Vector
CHARACTER(LEN = 1), INTENT(IN)                                  ::  Units

INTEGER, INTENT(IN)                                             :: imin       ! minimum x-array index
INTEGER, INTENT(IN)                                             :: imax       ! maximum x-array index
INTEGER, INTENT(IN)                                             :: nx         ! x-array extent
INTEGER, INTENT(IN)                                             :: ny         ! y-array extent
INTEGER, INTENT(IN)                                             :: nz         ! y-array extent
INTEGER, INTENT(IN)                                             :: ij_ray_dim ! number of y-zones on a processor before swapping
INTEGER, INTENT(IN)                                             :: ik_ray_dim ! number of z-zones on a processor before swapping

REAL(KIND = idp), INTENT(IN), DIMENSION(nx+1)                   :: x_e        ! x grid zone left interfaces [cm]
REAL(KIND = idp), INTENT(IN), DIMENSION(nx)                     :: x_c        ! x grid zone centers [cm]
REAL(KIND = idp), INTENT(IN), DIMENSION(nx)                     :: dx_c       ! x_e(i+1) - x_e(i) [cm]

REAL(KIND = idp), INTENT(IN), DIMENSION(-5:ny+7)                :: y_e        ! y grid zone left interfaces
REAL(KIND = idp), INTENT(IN), DIMENSION(-5:ny+6)                :: y_c        ! y grid zone centers
REAL(KIND = idp), INTENT(IN), DIMENSION(-5:ny+6)                :: dy_c       ! y_e(j+1) - y_e(j)

REAL(KIND = idp), INTENT(IN), DIMENSION(-5:nz+7)                :: z_e        ! z grid zone left interfaces
REAL(KIND = idp), INTENT(IN), DIMENSION(-5:nz+6)                :: dz_c       ! z_e(k+1) - z_e(k)

INTEGER, INTENT(IN), DIMENSION(1:3)                             ::  Num_Input_Nodes


REAL(KIND = idp), DIMENSION(1:Num_Input_Nodes(1)), INTENT(IN)   ::  Input_R_Quad
REAL(KIND = idp), DIMENSION(1:Num_Input_Nodes(2)), INTENT(IN)   ::  Input_T_Quad
REAL(KIND = idp), DIMENSION(1:Num_Input_Nodes(3)), INTENT(IN)   ::  Input_P_Quad


REAL(KIND = idp), INTENT(IN)                                    ::  Left_Limit,                 &
                                                                    Right_Limit

REAL(KIND = idp), INTENT(IN), DIMENSION(    1:Num_Input_Nodes(1)*Num_Input_Nodes(2)*Num_Input_Nodes(3),  &
                                            0:nx-1,                         &
                                            0:ij_ray_dim-1,                 &
                                            0:ik_ray_dim-1              )   ::  Local_E,       &
                                                                                Local_S

REAL(KIND = idp), INTENT(IN), DIMENSION(    1:Num_Input_Nodes(1)*Num_Input_Nodes(2)*Num_Input_Nodes(3),  &
                                            0:nx-1,                         &
                                            0:ij_ray_dim-1,                 &
                                            0:ik_ray_dim-1,                 &
                                            1:DOMAIN_DIM                )   ::  Local_Si




!                             !
!!                           !!
!!!     Local Variables     !!!
!!                           !!
!                             !
INTEGER                                                     ::  i
INTEGER                                                     ::  FEM_Degree, SH_Limit

REAL(KIND = idp), DIMENSION(1)                              ::  Input_Quad




INTEGER, DIMENSION(1:3)                                     ::  Num_Output_Nodes
REAL(KIND = idp), DIMENSION(1:2)                            ::  Output_Quad



REAL(KIND = idp)                                            ::  timea, timeb, tottime



CHARACTER(LEN=1), DIMENSION(1:5)                            ::  INNER_BC_TYPES, OUTER_BC_TYPES
REAL(KIND = idp), DIMENSION(1:5)                            ::  INNER_BC_VALUES, OUTER_BC_VALUES

REAL(KIND = idp)                                            ::  Outer_Potential,    &
                                                                Inner_Potential,    &
                                                                Pot_to_Alphapsi,    &
                                                                Pot_to_Psi



INTEGER                                                     ::  ierr

REAL(KIND = idp)                                            ::  ratio, q, pastloc


real(KIND = idp)                                            ::  Shift_Vector_BC

REAL(KIND = idp), DIMENSION(0:nx)                           ::  Shift_Vector




!                                         !
!!                                       !!
!!!         Initialize Poseidon         !!!
!!                                       !!
!                                         !
IF (( MODE == 0 ) .OR. ( MODE == 3 )) THEN
    ! Entered if this is the first call

    timea = MPI_Wtime()

    CALL Set_Units( Units )
    FEM_Degree = 1
    SH_Limit = 0
    IF ( Solver_Mode == 1 ) THEN
        CALL Poseidon_Initialize_From_File(   mode,                           & ! mode
                                    x_e(1),                         & ! Inner_Radius
                                    x_e(nx+1),                      & ! Outer_Radius
                                    nx,                             & ! NUM_R_ELEMENTS
                                    ny,                             & ! NUM_T_ELEMENTS
                                    nz,                             & ! NUM_P_ELEMENTS
                                    nx,                             & ! NUM_LOC_R_ELEMENTS
                                    ij_ray_dim,                     & ! NUM_LOC_T_ELEMENTS
                                    ik_ray_dim,                     & ! NUM_LOC_P_ELEMENTS
                                    dx_c,                           & ! Delta_R_Vector
                                    dy_c(1:ny),                     & ! Delta_T_Vector
                                    dz_c(1:nz)                      ) ! Delta_P_Vector)





!    CALL Poseidon_Initialize &
!         ( Units                  = "G",                  &
!           Dimensions             = 1,                    &
!           FEM_Degree_Input       = 1,                    &
!           L_Limit_Input          = 0,                    &
!           Inner_Radius           = 0.0_idp,              &
!           Outer_Radius           = 1.0E8_idp,            &
!           R_Elements_Input       = DRIVER_R_ELEMS,       &
!           T_Elements_Input       = DRIVER_T_ELEMS,       &
!           P_Elements_Input       = DRIVER_P_ELEMS,       &
!           Local_R_Elements_Input = DRIVER_R_ELEMS,       &
!           Local_T_Elements_Input = DRIVER_T_ELEMS,       &
!           Local_P_Elements_Input = DRIVER_P_ELEMS,       &
!           Num_R_Quad_Input       = DRIVER_R_INPUT_NODES, &
!           Num_T_Quad_Input       = DRIVER_T_INPUT_NODES, &
!           Num_P_Quad_Input       = DRIVER_P_INPUT_NODES, &
!           Input_Delta_R_Vector   = dx_c )

    ELSEIF ( Solver_Mode == 2 ) THEN

                CALL Poseidon_FP_Init_From_File( mode,                  & ! mode
                                                CFA_EQs_Flag_Vector,    & ! CFA_EQ_Flags
                                                x_e(1),                 & ! Inner_Radius
                                                x_e(nx+1),              & ! Outer_Radius
                                                nx,                     & ! NUM_R_ELEMENTS
                                                ny,                     & ! NUM_T_ELEMENTS
                                                nz,                     & ! NUM_P_ELEMENTS
                                                nx,                     & ! NUM_LOC_R_ELEMENTS
                                                ij_ray_dim,             & ! NUM_LOC_T_ELEMENTS
                                                ik_ray_dim,             & ! NUM_LOC_P_ELEMENTS
                                                dx_c,                   & ! Delta_R_Vector
                                                dy_c(1:ny),             & ! Delta_T_Vector
                                                dz_c(1:nz)              ) ! Delta_P_Vector)


    END IF




    timeb = MPI_Wtime()

!    PRINT*,"               Initialize Time :", timeb-timea, myID
    CALL Clock_In(timeb-timea, 1)


END IF



PRINT*,"After Poseidon_Init"


! (Global_TE/Ray_TE) * (1/NUM_THETA_BLOCKS) = NUM_CHIMERA_THETA_RAYS/NUM_THETA_BLOCKS
Ratio_T_BNDLperBLCK = ny/( NUM_BLOCK_THETA_ROWS*ij_ray_dim )
Ratio_P_BNDLperBLCK = nz/( NUM_BLOCK_PHI_COLUMNS*ik_ray_dim )

Ratio_BNDLperBLCK = Ratio_T_BNDLperBLCK * Ratio_P_BNDLperBLCK


!
!   When the Radial Mesh changes the stiffness matrix needs to be recomputed.
!   This call sets the new mesh and generates the new stiffness matrix. 
!   This doesn't need to be called if Poseidon was just initialized above
!   as the mesh is set and matrix generated as part of the initialization.
!
IF ( ( MODE .NE. 0) .AND. ( MODEb == 1 ) )  THEN

    CALL Poseidon_Set_Mesh(dx_c)

END IF





!                                         !
!!                                       !!
!!!       Insert Source Variables       !!!
!!                                       !!
!                                         !
timea = MPI_Wtime()


CALL Poseidon_Input_Sources( myID, myID_Theta, myID_Phi,                                &
                             Local_E, Local_S, Local_Si,                                &
                             nx, ij_ray_dim, ik_ray_dim,                                &
                             Num_Input_Nodes(1),Num_Input_Nodes(2),Num_Input_Nodes(3),  &
                             Input_R_Quad, Input_T_Quad, Input_P_Quad,                  &
                             Left_Limit, Right_Limit                                    )



timeb = MPI_Wtime()
CALL Clock_In(timeb-timea, 2)




PRINT*,"After Poseidon_Input_Sources"



IF ( POSEIDON_COMM_WORLD .NE. MPI_COMM_NULL ) THEN


    !                                         !
    !!                                       !!
    !!!    Calculate Boundary Conditions    !!!
    !!                                       !!
    !                                         !
    timea = MPI_Wtime()

    
    ! Calculate the Dirichlet Outer Boundary Value for the Shift Vector
    CALL Calc_Shift_BC_1D( Shift_Vector_BC,                                                 &
                           NUM_Input_Nodes(1), NUM_Input_Nodes(2), NUM_Input_Nodes(3),      &
                           NUM_R_ELEMENTS, ij_ray_dim, ik_ray_dim, DOMAIN_DIM,              &
                           Local_Si, rlocs          )
!    Shift_Vector_BC = 0.0_idp


    ! Calculate the Dirichlet Outer Boundary Value for the Lapse Function and Conformal Factor
    Inner_Potential = Potential_Solution(rlocs(0), 0.0_idp, 0.0_idp)
    Outer_Potential = Potential_Solution(rlocs(NUM_R_ELEMENTS), 0.0_idp, 0.0_idp)




    !PRINT*,"ALPHAPSI BC set to one"
    !Pot_to_AlphaPsi = 1.0_idp
    Pot_to_Alphapsi = 1.0_idp + 0.5_idp*Outer_Potential/C_Square
    Pot_to_Psi      = 1.0_idp - 0.5_idp*Outer_Potential/C_Square



    !                                         !
    !!                                       !!
    !!!       Set Boundary Conditions       !!!
    !!                                       !!
    !                                         !

    ! Set BC type, N = Neumann, D = Dirichlet
    INNER_BC_TYPES = (/"N", "N","N","N","N"/)
    OUTER_BC_TYPES = (/"D", "D","D","D","D"/)




    ! Set BC Values
    INNER_BC_VALUES = (/0.0_idp, 0.0_idp, 0.0_idp, 0.0_idp, 0.0_idp /)
    OUTER_BC_VALUES = (/Pot_to_Psi, Pot_to_AlphaPsi, Shift_Vector_BC, 0.0_idp, 0.0_idp /)
!    OUTER_BC_VALUES = (/Pot_to_Psi, Pot_to_AlphaPsi, 0.0_idp, 0.0_idp, 0.0_idp /)
!    OUTER_BC_VALUES = (/1.0_idp, 1.0_idp, 0.0_idp, 0.0_idp, 0.0_idp /)

!    Shift_Vector_BC = -1.0E2_idp
!    OUTER_BC_VALUES = (/0.0_idp, 0.0_idp, Shift_Vector_BC, 0.0_idp, 0.0_idp /)


!    PRINT*,"FRAME = ",DRIVER_FRAME," Shift BC = ",Shift_Vector_BC
    
    ! Commit BC types and values to Poseidon
    CALL Poseidon_CFA_Set_Uniform_Boundary_Conditions("I", INNER_BC_TYPES, INNER_BC_VALUES)
    CALL Poseidon_CFA_Set_Uniform_Boundary_Conditions("O", OUTER_BC_TYPES, OUTER_BC_VALUES)


    timeb = MPI_Wtime()
    CALL Clock_In(timeb-timea, 3)

    Print*,"After BC_Set"


    !                                         !
    !!                                       !!
    !!!          Set Initial Guess          !!!
    !!                                       !!
    !                                         !
    IF (( MODE == 0 ) .OR. ( MODE == 3 )) THEN
        ! Entered only on first Frame

        SELECT CASE ( DRIVER_FIRST_GUESS_FLAG )
            CASE( 1 )   ! Flat Space Guess
                CALL Initialize_Flat_Space_Guess_Values
            CASE( 2 )   ! Calculated Guess ( Tests Only )
                CALL Initialize_Calculated_Guess_Values
            CASE( 3 )   ! Load From File
                CALL Load_Initial_Guess_From_File
        END SELECT
    ELSE
        ! Entered on all subsequent frames

        SELECT CASE ( DRIVER_SUBSEQUENT_GUESS_FLAG )
            CASE( 1 )   ! Flat Space Guess
                CALL Initialize_Flat_Space_Guess_Values
            CASE( 2 )   ! Calculated Guess ( Tests Only )
                CALL Initialize_Calculated_Guess_Values
            CASE( 3 )   ! Load From File

            
            CASE( 4 )   ! Use Solution from Last Frame
                ! No Action Required

        END SELECT
    END IF

    PRINT*,"After Guess"
    !                                         !
    !!                                       !!
    !!!             Run Poseidon            !!!
    !!                                       !!
    !                                         !
    ! Run Poseidon. Upon Completion of this subroutine, the coefficients
    ! for the solution vector will be known.
    CALL Poseidon_Run()
    PRINT*,"After Run"

END IF

timeb = MPI_Wtime()
CALL Poseidon_Distribute_Solution()
timea = MPI_Wtime()

CALL Clock_In(timea-timeb, 19)





        !                                         !
        !!                                       !!
        !!!            Output Results           !!!
        !!                                       !!
        !                                         !
!IF ( myID == 0 ) THEN
!    
!    CALL WRITE_CFA_COEFFICIENTS()
!END IF






END SUBROUTINE Poseidon_CFA_3D




























 !+202+############################################################################!
!                                                                                   !
!                    Poseidon_Shutdown                                              !
!                                                                                   !
 !#################################################################################!
SUBROUTINE Poseidon_Shutdown()

CALL Poseidon_Close()

END SUBROUTINE Poseidon_Shutdown






END MODULE Poseidon_Interface
