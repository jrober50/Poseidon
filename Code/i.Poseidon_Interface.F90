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



USE Poseidon_Constants_Module, &
            ONLY :  idp, pi,speed_of_light


USE Units_Module, &
            ONLY :  Set_Units, Grav_Constant_G


USE CHIMERA_Parameters,  &
            ONLY :  Analytic_Solution,      &
                    Shift_Solution,         &
                    Enclosed_Mass,          &
                    CHIMERA_R_LOCS,         &
                    nPROCS,                 &
                    myID,                   &
                    myID_Theta,             &
                    myID_Phi,               &
                    Ratio_T_BNDLperBLCK,    &
                    Ratio_P_BNDLperBLCK,    &
                    Ratio_BNDLperBLCK,      &
                    NUM_ENTRIES,            &
                    SELFSIM_R_VALS,         &
                    SELFSIM_T


USE Poseidon_Parameters, &
            ONLY :  DOMAIN_DIM,             &
                    DEGREE,                 &
                    L_LIMIT,                &
                    NUM_CFA_VARS,           &
                    DATA_DIST_MODE,         &
                    NUM_R_ELEMS_PER_SHELL,  &
                    NUM_SHELLS,             &
                    NUM_BLOCKS_PER_SHELL,   &
                    NUM_BLOCK_THETA_ROWS,   &
                    NUM_BLOCK_PHI_COLUMNS,  &
                    nPROCS_POSEIDON,        &
                    STF_MAPPING_FLAG,       &
                    NUM_R_ELEMS_PER_BLOCK,  &
                    NUM_T_ELEMS_PER_BLOCK,  &
                    NUM_P_ELEMS_PER_BLOCK,  &
                    NUM_R_QUAD_POINTS,      &
                    NUM_T_QUAD_POINTS,      &
                    NUM_P_QUAD_POINTS,      &
                    R_COARSEN_FACTOR,       &
                    T_COARSEN_FACTOR,       &
                    P_COARSEN_FACTOR,       &
                    SOL_DIST_SCHEME


USE Global_Variables_And_Parameters, &
            ONLY :  ierr,                                               &
                    myID_Poseidon,                                      &
                    myID_Shell,                                         &
                    myID_PETSc,                                         &
                    myID_Dist,                                          &
                    myShell,                                            &
                    FirstCall_Flag,                                     &
                    Stiffness_Matrix_Initialized_Flag,                  &
                    Matrix_Cholesky_Factorized_Flag,                    &
                    rlocs,                                              &
                    NUM_R_ELEMENTS,                                     &
                    NUM_T_ELEMENTS,                                     &
                    NUM_P_ELEMENTS,                                     &
                    R_INNER, R_OUTER,                                   &
                    Block_Source_E,                                     &
                    Block_Source_S,                                     &
                    Block_Source_Si,                                    &
                    POSEIDON_COMM_WORLD,                                &
                    POSEIDON_COMM_SHELL,                                &
                    POSEIDON_COMM_DIST,                                 &
                    Prob_Dim,                                           &
                    Block_Prob_Dim,                                     &
                    LM_Length,                                          &
                    ULM_Length,                                         &
                    Local_Length,                                       &
                    Coefficient_Vector,                                 &
                    rlocs, tlocs, plocs


USE Poseidon_Main_Module, &
            ONLY :  Poseidon_Initialize,                                &
                    Poseidon_Run,                                       &
                    Poseidon_Close,                                     &
                    Poseidon_Set_Mesh,                                  &
                    Poseidon_CFA_Set_Uniform_Boundary_Conditions




USE CFA_3D_Master_Build_Module,     &
            ONLY : Calc_3D_Values_At_Location


USE Additional_Functions_Module, &
            ONLY :  Lagrange_Poly,                                      &
                    Spherical_Harmonic,                                 &
                    Map_To_X_Space, Map_From_X_Space,                   &
                    Initialize_LGL_Quadrature,                          &
                    Initialize_LGL_Quadrature_Locations,                &
                    Initialize_LG_Quadrature_Locations,                 &
                    MVMULT_FULL


USE IO_Functions_Module, &
            ONLY :  Clock_In,                &
                    PRINT_TIMETABLE

USE Poseidon_Parameter_Read_Module, &
            ONLY :  WRITE_CFA_COEFFICIENTS,  &
                    READ_CFA_COEFFICIENTS 



USE Jacobian_Internal_Functions_Module, &
            ONLY :  Calc_Shift_BC_1D,                 &
                    Calc_Shift_1D,                    &
                    Write_Shift_1D,                   &
                    Initialize_Special_Guess_Values


USE Mesh_Module, &
            ONLY :  Create_Uniform_1D_Mesh,             &
                    Create_Logarithmic_1D_Mesh,         &
                    Create_Split_1D_Mesh



USE Poseidon_Internal_Communication_Module,             &
            ONLY :  Poseidon_CrsThrd_Block_Share,       &
                    Poseidon_Distribute_Solution


use mpi








IMPLICIT NONE

CONTAINS









SUBROUTINE Poseidon_CFA_3D(mode, modeb, imin, imax, nx, ij_ray_dim, ik_ray_dim,         &
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





INTEGER                                                     ::  h,i,j,k,l,my_j, my_k






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



REAL(KIND = idp)                                            ::  deltar, r,          &
                                                                theta, phi,         &
                                                                Analytic_Val,       &
                                                                Solver_Val,         &
                                                                Solver_Valb,        & 
                                                                Error_Val


REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Output_re,          &
                                                                Output_rc,          &
                                                                Output_dr



INTEGER                                                     ::  NUM_SAMPLES


REAL(KIND = idp)                                            ::  Return_Psi,         &
                                                                Return_AlphaPsi,    &
                                                                Return_Beta1,       &
                                                                Return_Beta2,       &
                                                                Return_Beta3


REAL(KIND = idp)                                            ::  csqr
INTEGER                                                     ::  ierr

REAL(KIND = idp)                                            ::  ratio, q, pastloc


real(KIND = idp)                                            ::  Shift_Vector_BC

REAL(KIND = idp), DIMENSION(0:nx)                           ::  Shift_Vector


LOGICAL                        ::  Output_to_File
CHARACTER(LEN = 19)            ::  filenamea,filenameb
INTEGER                        ::  file_ida, file_idb


Output_to_File = .TRUE. 

110 FORMAT (11X,A1,16X,A18,9X,A13,10X,A18,10X,A10)                     !!! Output Header

111 FORMAT (11X,A1,18X,A13,10X,A18,10X,A11,14X,A11,14X,A11)             !!! Output Header for Results file
112 FORMAT (11X,A1,16X,A18,9x,A14)                                     !!! Output Header for Analytic Solution file

113 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15)               !!! Output
114 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15)    !!! Output for Results file
115 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15)                                     !!! Output for Analytic Solution file


!                                         !
!!                                       !!
!!!         Initialize Poseidon         !!!
!!                                       !!
!                                         !

IF (MODE == 1) THEN


    !                                         !
    !!                                       !!
    !!!         Initialize Poseidon         !!!
    !!                                       !!
    !                                         !

    timea = MPI_Wtime()

    CALL Set_Units( "C" )



    CALL Poseidon_Initialize(   x_e(1),                         & ! Inner_Radius
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

    timeb = MPI_Wtime()

!    PRINT*,"               Initialize Time :", timeb-timea, myID
    CALL Clock_In(timeb-timea, 1)





END IF



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
IF ( ( MODE .NE. 1) .AND. ( MODEb == 1 ) )  THEN
    PRINT*,"here"
    CALL Poseidon_Set_Mesh(dx_c)

END IF





!                                         !
!!                                       !!
!!!       Insert Source Variables       !!!
!!                                       !!
!                                         !



timea = MPI_Wtime()




IF (myID == 0) THEN
     PRINT*,"Before Pseidon_CrsThrd_Block_Share",myID
END IF

CALL Poseidon_CrsThrd_Block_Share(  myID, myID_Theta, myID_Phi,                              &
                                   Local_E, Local_S, Local_Si,                               &
                                   nx, ij_ray_dim, ik_ray_dim,                               &
                                   Num_Input_Nodes(1),Num_Input_Nodes(2),Num_Input_Nodes(2), &
                                   Input_R_Quad, Input_T_Quad, Input_P_Quad,                 &
                                   Left_Limit, Right_Limit,                                  &
                                   nx, ny, nz,                                               &
                                   dx_c, dy_c(1:ny), dz_c(1:nz),                             &
                                   Block_Source_E, Block_Source_S, Block_Source_Si           )




timeb = MPI_Wtime()
CALL Clock_In(timeb-timea, 2)









IF ( POSEIDON_COMM_WORLD .NE. MPI_COMM_NULL ) THEN


    !                                         !
    !!                                       !!
    !!!       Set Boundary Conditions       !!!
    !!                                       !!
    !                                         !

    CALL Initialize_Special_Guess_Values()


!    CALL Calc_Shift_1D( Shift_Vector,                                    &
!                           NUM_Input_Nodes(1), NUM_Input_Nodes(2), NUM_Input_Nodes(3),     &
!                           NUM_R_ELEMENTS, ij_ray_dim, ik_ray_dim, DOMAIN_DIM,             &
!                           Local_Si, rlocs          )


!    CALL Write_Shift_1D( Shift_Vector,                                                     &
!                         NUM_R_ELEMENTS,                                                   &
!                         rlocs,                                                            &
!                         SELFSIM_T                  )

    CALL Calc_Shift_BC_1D( Shift_Vector_BC,                       &
                           NUM_Input_Nodes(1), NUM_Input_Nodes(2), NUM_Input_Nodes(3),     &
                           NUM_R_ELEMENTS, ij_ray_dim, ik_ray_dim, DOMAIN_DIM,             &
                           Local_Si, rlocs          ) 





    timea = MPI_Wtime()


    Inner_Potential = Analytic_Solution(rlocs(0), 0.0_idp, 0.0_idp)
    Outer_Potential = Analytic_Solution(rlocs(NUM_R_ELEMENTS), 0.0_idp, 0.0_idp)
    
!    Inner_Potential = 0.0_idp
!    Outer_Potential = Grav_Constant_G *Enclosed_Mass(NUM_R_ELEMENTS)            &
!                    / CHIMERA_R_LOCS(NUM_R_ELEMENTS)



    csqr = Speed_of_Light*Speed_of_Light

    Pot_to_Alphapsi = 1.0_idp + 0.5_idp*Outer_Potential/csqr
    Pot_to_Psi = 1.0_idp - 0.5_idp*Outer_Potential/csqr


    INNER_BC_TYPES = (/"N", "N","N","N","N"/)

    OUTER_BC_TYPES = (/"D", "D","D","D","D"/)


    Pot_to_Alphapsi = 1.0_idp + 0.5_idp*Inner_Potential/csqr


    INNER_BC_VALUES = (/0.0_idp, 0.0_idp, 0.0_idp, 0.0_idp, 0.0_idp /)

    Pot_to_Alphapsi = 1.0_idp + 0.5_idp*Outer_Potential/csqr

    OUTER_BC_VALUES = (/Pot_to_Psi, Pot_to_AlphaPsi, Shift_Vector_BC, 0.0_idp, 0.0_idp /)


    CALL Poseidon_CFA_Set_Uniform_Boundary_Conditions("I", INNER_BC_TYPES, INNER_BC_VALUES)

    CALL Poseidon_CFA_Set_Uniform_Boundary_Conditions("O", OUTER_BC_TYPES, OUTER_BC_VALUES)


    timeb = MPI_Wtime()
    CALL Clock_In(timeb-timea, 3)







        !                                         !
        !!                                       !!
        !!!             Run Poseidon            !!!
        !!                                       !!
        !                                         !
    CALL Poseidon_Run()


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
IF (( .TRUE. ) .AND. (myID == 0)) THEN
    NUM_SAMPLES = 1000

    ALLOCATE( Output_re(0:NUM_SAMPLES) )
    ALLOCATE( Output_rc(1:NUM_SAMPLES) )
    ALLOCATE( Output_dr(1:NUM_SAMPLES) )


    ! Open Results File
    file_ida = 42
    filenamea = "OUTPUT/Results.out"
    OPEN( Unit = file_ida, file = filenamea )
    WRITE(file_ida,111)"r","Psi Potential","AlphaPsi Potential","Beta1 Value","Beta2 Value","Beta3 Value"

    ! Open Solution File
    file_idb = 43
    filenameb = "OUTPUT/Solution.out"
    OPEN( Unit = file_idb, file = filenameb )
    WRITE(file_idb,112)"r","Analytic Potential","Beta1 Solution"


    ! Set Output locations
    theta = pi/2.0_idp
    phi = pi/2.0_idp

    CALL Create_Logarithmic_1D_Mesh( R_INNER, R_OUTER, NUM_SAMPLES,        &
                                    output_re, output_rc, output_dr                 )



    DO i = 1,NUM_SAMPLES

        CALL Calc_3D_Values_At_Location( output_rc(i), theta, phi,                   &
                                         Return_Psi, Return_AlphaPsi,                &
                                         Return_Beta1, Return_Beta2, Return_Beta3    )

        Analytic_Val = Analytic_Solution(output_rc(i),theta,phi)
        Solver_Val = 2.0_idp*csqr*(1.0_idp - Return_Psi)
        Solver_Valb = 2.0_idp*csqr*(Return_AlphaPsi - 1.0_idp)
        Error_Val = ABS((Analytic_Val - Solver_Val)/Analytic_Val)


       WRITE(42,114) output_rc(i), Return_Psi, Return_AlphaPsi, Return_Beta1,Return_Beta2,Return_Beta3
       WRITE(43,115) output_rc(i),Analytic_Val, Shift_Solution(output_rc(i),rlocs,NUM_R_ELEMENTS)



    END DO


    ! Close Files
    CLOSE( Unit = file_ida)
    CLOSE( Unit = file_idb)


END IF




        !                                         !
        !!                                       !!
        !!!         Output Coefficients         !!!
        !!                                       !!
        !                                         !
IF ( myID == 0) THEN

    CALL WRITE_CFA_COEFFICIENTS()

    CALL READ_CFA_COEFFICIENTS()

END IF





IF ( myID_PETSC == 0) THEN 

   CALL Print_TimeTable(0)

END IF




END SUBROUTINE Poseidon_CFA_3D










SUBROUTINE Poseidon_Newt_3D(mode, modeb, imin, imax, nx, ij_ray_dim, ik_ray_dim,         &
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





INTEGER                                                     ::  h,i,j,k,l,my_j, my_k






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



REAL(KIND = idp)                                            ::  deltar, r,          &
                                                                theta, phi,         &
                                                                Analytic_Val,       &
                                                                Solver_Val,         &
                                                                Solver_Valb,        &
                                                                Error_Val


REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Output_re,          &
                                                                Output_rc,          &
                                                                Output_dr



INTEGER                                                     ::  NUM_SAMPLES


REAL(KIND = idp)                                            ::  Return_Psi,         &
                                                                Return_AlphaPsi,    &
                                                                Return_Beta1,       &
                                                                Return_Beta2,       &
                                                                Return_Beta3


REAL(KIND = idp)                                            ::  csqr
INTEGER                                                     ::  ierr

REAL(KIND = idp)                                            ::  ratio, q, pastloc


real(KIND = idp)                                            ::  Shift_Vector_BC

REAL(KIND = idp), DIMENSION(0:nx)                           ::  Shift_Vector


LOGICAL                        ::  Output_to_File
CHARACTER(LEN = 19)            ::  filenamea,filenameb
INTEGER                        ::  file_ida, file_idb


Output_to_File = .TRUE.

110 FORMAT (11X,A1,16X,A18,9X,A13,10X,A18,10X,A10)                     !!! Output Header

111 FORMAT (11X,A1,18X,A13,10X,A18,10X,A11,14X,A11,14X,A11)             !!! Output Header for Results file
112 FORMAT (11X,A1,16X,A18,9x,A14)                                     !!! Output Header for Analytic Solution file

113 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15)               !!! Output
114 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15)    !!! Output for Results file
115 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15)                                     !!! Output for Analytic Solution file


!                                         !
!!                                       !!
!!!         Initialize Poseidon         !!!
!!                                       !!
!                                         !

IF (MODE == 1) THEN


    !                                         !
    !!                                       !!
    !!!         Initialize Poseidon         !!!
    !!                                       !!
    !                                         !

    timea = MPI_Wtime()

    CALL Set_Units( "C" )



    CALL Poseidon_Initialize(   x_e(1),                         & ! Inner_Radius
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

    timeb = MPI_Wtime()

!    PRINT*,"               Initialize Time :", timeb-timea, myID
    CALL Clock_In(timeb-timea, 1)





END IF










IF ( POSEIDON_COMM_WORLD .NE. MPI_COMM_NULL ) THEN


    !                                         !
    !!                                       !!
    !!!       Set Boundary Conditions       !!!
    !!                                       !!
    !                                         !

    CALL Initialize_Special_Guess_Values()


!    CALL Calc_Shift_1D( Shift_Vector,                                    &
!                           NUM_Input_Nodes(1), NUM_Input_Nodes(2), NUM_Input_Nodes(3),     &
!                           NUM_R_ELEMENTS, ij_ray_dim, ik_ray_dim, DOMAIN_DIM,             &
!                           Local_Si, rlocs          )


!    CALL Write_Shift_1D( Shift_Vector,                                                     &
!                         NUM_R_ELEMENTS,                                                   &
!                         rlocs,                                                            &
!                         SELFSIM_T                  )

    CALL Calc_Shift_BC_1D( Shift_Vector_BC,                       &
                           NUM_Input_Nodes(1), NUM_Input_Nodes(2), NUM_Input_Nodes(3),     &
                           NUM_R_ELEMENTS, ij_ray_dim, ik_ray_dim, DOMAIN_DIM,             &
                           Local_Si, rlocs          )





    timea = MPI_Wtime()


    Inner_Potential = Analytic_Solution(rlocs(0), 0.0_idp, 0.0_idp)
    Outer_Potential = Analytic_Solution(rlocs(NUM_R_ELEMENTS), 0.0_idp, 0.0_idp)

!    Inner_Potential = 0.0_idp
!    Outer_Potential = Grav_Constant_G *Enclosed_Mass(NUM_R_ELEMENTS)            &
!                    / CHIMERA_R_LOCS(NUM_R_ELEMENTS)



    csqr = Speed_of_Light*Speed_of_Light

    Pot_to_Alphapsi = 1.0_idp + 0.5_idp*Outer_Potential/csqr
    Pot_to_Psi = 1.0_idp - 0.5_idp*Outer_Potential/csqr


    INNER_BC_TYPES = (/"N", "N","N","N","N"/)

    OUTER_BC_TYPES = (/"D", "D","D","D","D"/)


    Pot_to_Alphapsi = 1.0_idp + 0.5_idp*Inner_Potential/csqr


    INNER_BC_VALUES = (/0.0_idp, 0.0_idp, 0.0_idp, 0.0_idp, 0.0_idp /)

    Pot_to_Alphapsi = 1.0_idp + 0.5_idp*Outer_Potential/csqr

    OUTER_BC_VALUES = (/Pot_to_Psi, Pot_to_AlphaPsi, Shift_Vector_BC, 0.0_idp, 0.0_idp /)


    CALL Poseidon_CFA_Set_Uniform_Boundary_Conditions("I", INNER_BC_TYPES, INNER_BC_VALUES)

    CALL Poseidon_CFA_Set_Uniform_Boundary_Conditions("O", OUTER_BC_TYPES, OUTER_BC_VALUES)


    timeb = MPI_Wtime()
    CALL Clock_In(timeb-timea, 3)







        !                                         !
        !!                                       !!
        !!!             Run Poseidon            !!!
        !!                                       !!
        !                                         !
    CALL Poseidon_Run()


END IF





        !                                         !
        !!                                       !!
        !!!         Output Coefficients         !!!
        !!                                       !!
        !                                         !
IF ( myID == 0) THEN

    CALL WRITE_CFA_COEFFICIENTS()

    CALL READ_CFA_COEFFICIENTS()

END IF





IF ( myID_PETSC == 0) THEN

   CALL Print_TimeTable(0)

END IF




END SUBROUTINE Poseidon_Newt_3D































END MODULE Poseidon_Interface
