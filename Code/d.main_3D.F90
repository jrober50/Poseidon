   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
PROGRAM Poseidon_CFA                                                                !##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!    This program is an example of how to run the Newtonian side of Poseidon.    !##!
!##!                                                                                !##!
!##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!

!*D*================================!
!                                   !
!           Dependencies            !
!                                   !
!===================================!





USE Poseidon_Constants_Module, &
            ONLY :  idp, pi

USE CHIMERA_Parameters, &
            ONLY :  CHIMERA_R_ELEMS,        &
                    CHIMERA_C_ELEMS,        &
                    CHIMERA_T_ELEMS,        &
                    CHIMERA_P_ELEMS,        &
                    CHIMERA_R_INPUT_NODES,  &
                    CHIMERA_T_INPUT_NODES,  &
                    CHIMERA_P_INPUT_NODES,  &
                    CHIMERA_LEFT_LIMIT,     &
                    CHIMERA_RIGHT_LIMIT,    &
                    CHIMERA_INNER_RADIUS,   &
                    CHIMERA_CORE_RADIUS,    &
                    CHIMERA_OUTER_RADIUS,   &
                    CHIMERA_MESH_TYPE,      &
                    CHIMERA_R_LOCS,         &
                    CHIMERA_T_LOCS,         &
                    CHIMERA_P_LOCS,         &
                    CHIMERA_Delta_R,        &
                    CHIMERA_Delta_T,        &
                    CHIMERA_Delta_P,        &
                    ENCLOSED_MASS,          &
                    CHIMERA_Potential,      &
                    CHIMERA_E,              &
                    CHIMERA_S,              &
                    CHIMERA_Si,             &
                    CHIMERA_DIMENSION,      &
                    CHIMERA_PROCS,          &
                    CHIMERA_y_PROCS,        &
                    CHIMERA_z_PROCS,        &
                    MPI_COMM_XY,            &
                    MPI_COMM_XZ,            &
                    MPI_COMM_GRID,          &
                    nPROCS,                 &
                    ij_ray_dim,             &
                    ik_ray_dim,             &
                    myID,                   &
                    myID_theta,             &
                    myID_phi,               &
                    ngrid,                  &
                    SELFSIM_T,              &
                    SELFSIM_KAPPA,          &
                    SELFSIM_GAMMA

USE CHIMERA_Params_Read_Module, &
            ONLY :  Unpack_CHIMERA_Parameters

USE SelfSimilar_Module, &
            ONLY :  UNPACK_SELF_SIMILAR

USE CHIMERA_HDF5_Module, &
            ONLY :  Load_CHIMERA_HDF5


USE Mesh_Module, &
            ONLY :  Create_3D_Mesh


USE Units_Module, &
            ONLY :  Set_Units, Grav_Constant_G, Speed_of_Light




USE Test_Functions_Module, &
            ONLY :  Poseidon_Initialize_CFA_Test_Prob_Par_CHIMERA






USE Poseidon_Interface, &
            ONLY :  Poseidon_CFA_3D


USE Additional_Functions_Module, &
            ONLY :  Map_From_X_Space,                       &
                    Initialize_LG_Quadrature_Locations




USE MPI

IMPLICIT NONE











                !*I*============================================!
                !                                               !
                !            Variable Initialization            !
                !                                               !
                !===============================================!

INTEGER                                                     ::   i, j, re      ! DO Loop Counter Variables







!                                       !
!   Poseidon Initialization Variables   !
!                                       !

REAL(KIND = idp)                                            ::  Inner_Radius,           &
                                                                Core_Radius,            &
                                                                Outer_Radius


INTEGER                                                     ::  nx,     &
                                                                nc,     &
                                                                ny,     &
                                                                nz

REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  dx_c
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  dy_c, tmp_dy_c
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  dz_c, tmp_dz_c
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  x_e, x_c, y_e, y_c, z_e, z_c






!                               !
!   Boundary Value Variables    !
!                               !

!REAL(KIND = idp)                                            :: Enclosed_Mass




!                       !
!   Output Variables    !
!                       !

INTEGER                                                     ::  NUM_SAMPLES
REAL(KIND = idp)                                            ::  potential


REAL(KIND = idp)                                            ::  Output_Left_Limit,  &
                                                                Output_Right_Limit

INTEGER, DIMENSION(1:3)                                     ::  Num_Output_Nodes

REAL(KIND = idp), DIMENSION(:,:,:,:), ALLOCATABLE           ::  Potential_Output



REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Output_R_Quad
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Output_T_Quad
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Output_P_Quad




!                                                       !
!   Other Variables Used to Initialize or Store Values  !
!                                                       !

REAL(KIND = idp)                                            ::  deltar, r, theta, phi
INTEGER                                                     ::  Num_Out_DOF

LOGICAL                                                     ::  RUN_POSEIDON_FLAG = .TRUE.


INTEGER                                                     ::  Problem_Dimension


INTEGER                                                     :: mode, modeb, imin, imax




INTEGER                                                     :: Num_Procs
INTEGER                                                     :: Num_y_Procs, Num_z_Procs

INTEGER                                                     :: n_y_procs, n_z_procs
INTEGER                                                     :: j_block_col, k_block_row

INTEGER                                                     ::  CHIMERA_COMM_WORLD,     &
                                                                CHIMERA_GROUP_WORLD,    &
                                                                MPI_GROUP_WORLD

INTEGER, DIMENSION(:), ALLOCATABLE                          ::  Workers_Array
INTEGER                                                     ::  myID_CHIMERA


INTEGER                                                     ::  ierr

!                           !
!   Source Input Variables  !
!                           !
REAL(KIND = idp)                                            ::  t, kappa, gamma


REAL(KIND = idp)                                            ::  Left_Limit,             &
                                                                Right_Limit

INTEGER                                                     ::  Num_DOF
INTEGER, DIMENSION(1:3)                                     ::  Num_Input_Nodes

REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Input_R_Quad
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Input_T_Quad
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Input_P_Quad
REAL(KIND = idp), DIMENSION(:,:,:,:), ALLOCATABLE           ::  Rho

REAL(KIND = idp), DIMENSION(:,:,:,:), ALLOCATABLE           ::  Local_E
REAL(KIND = idp), DIMENSION(:,:,:,:), ALLOCATABLE           ::  Local_S
REAL(KIND = idp), DIMENSION(:,:,:,:,:), ALLOCATABLE         ::  Local_Si


INTEGER                       :: Test_Num
INTEGER                       :: Mesh_Type


INTEGER                       :: file_number=00002
INTEGER                       :: stride=1

LOGICAL                       :: frame_flag=.TRUE.

CHARACTER(len=12)              :: path='Data3/Frames'

CALL MPI_INIT(ierr)

CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nPROCS, ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)



!!                                      !!
!!      CHIMERA Setup Parameters        !!
!!                                      !!

CALL Unpack_CHIMERA_Parameters()





ALLOCATE( CHIMERA_R_LOCS(0:CHIMERA_R_ELEMS+1), CHIMERA_T_LOCS(0:CHIMERA_T_ELEMS+1), CHIMERA_P_LOCS(0:CHIMERA_P_ELEMS+1))
ALLOCATE( CHIMERA_Delta_R(0:CHIMERA_R_ELEMS), CHIMERA_Delta_T(0:CHIMERA_T_ELEMS), CHIMERA_DELTA_P(0:CHIMERA_P_ELEMS))
ALLOCATE( Enclosed_Mass(0:CHIMERA_R_ELEMS+1) )
ALLOCATE( CHIMERA_Potential(0:CHIMERA_R_ELEMS) )
ALLOCATE( CHIMERA_E(1:CHIMERA_R_ELEMS,1:CHIMERA_T_ELEMS, 1:CHIMERA_P_ELEMS) )
ALLOCATE( CHIMERA_S(1:CHIMERA_R_ELEMS,1:CHIMERA_T_ELEMS, 1:CHIMERA_P_ELEMS) )
ALLOCATE( CHIMERA_Si(1:CHIMERA_R_ELEMS,1:CHIMERA_T_ELEMS, 1:CHIMERA_P_ELEMS, 1:3) )

Left_Limit = CHIMERA_LEFT_LIMIT
Right_Limit = CHIMERA_RIGHT_LIMIT

Num_Input_Nodes(1) = CHIMERA_R_INPUT_NODES
Num_Input_Nodes(2) = CHIMERA_T_INPUT_NODES
Num_Input_Nodes(3) = CHIMERA_P_INPUT_NODES
Num_DOF = Num_Input_Nodes(1)*Num_Input_Nodes(2)*Num_Input_Nodes(3)

ALLOCATE(   Input_R_Quad(1:Num_Input_Nodes(1)),         &
            Input_T_Quad(1:Num_Input_Nodes(2)),         &
            Input_P_Quad(1:Num_Input_Nodes(3))          )







!RUN_POSEIDON_FLAG = .FALSE.

CALL Set_Units( "C" )





!!                                       !!
!!   Poseidon Initialization Variables   !!
!!                                       !!


Num_Procs = CHIMERA_PROCS
Num_y_Procs = CHIMERA_y_PROCS
Num_z_Procs = CHIMERA_z_PROCS

Problem_Dimension = CHIMERA_DIMENSION


!   Set the Inner Radius of the Problem !
Inner_Radius = CHIMERA_INNER_RADIUS

!   Set the Core Radius of the Problem !
Core_Radius = CHIMERA_CORE_RADIUS

!   Set the Outer Radius of the Problem !
Outer_Radius = CHIMERA_OUTER_RADIUS !* 15 * 1000 * 100





!   Set the Number of Radial Elements   !
!nx = 4
nx = CHIMERA_R_ELEMS
nc = CHIMERA_C_ELEMS


!   Set the Number of Theta Elements    !
!ny = 2*Num_Procs                  !   For a 1-Dimensional Simulation Set equal to 1
ny = CHIMERA_T_ELEMS


!   Set the Number of Phi Elements      !
!nz = 2*Num_Procs                  !   For a 1 or 2-Dimensional Simulations Set equal to 1
nz = CHIMERA_P_ELEMS





MESH_TYPE = CHIMERA_MESH_TYPE
TEST_NUM = 3

NUM_SAMPLES = 20



PRINT*,"Test Number ",TEST_NUM
PRINT*,"Mesh Type   ",Mesh_Type





!                                       !
!       Initialize the MPI World        !
!                                       !


!PRINT*,"BEFORE MPI_INIT"


!CALL MPI_INIT(ierr)
!CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nPROCS, ierr)
!CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

!IF( myID == 0) THEN

!PRINT*,"BEFORE CREATE CHIMERA_COMM_WORLD",NUM_PROCS, nPROCS
!PRINT*,"   "
!PRINT*,"  "
!END IF


IF ( NUM_PROCS == nPROCS ) THEN

    CALL MPI_COMM_DUP(MPI_COMM_WORLD, CHIMERA_COMM_WORLD, ierr)
    CALL MPI_COMM_RANK(CHIMERA_COMM_WORLD, myID_Chimera, ierr)

ELSE IF ( NUM_PROCS < nPROCs ) THEN



    CALL MPI_Comm_group(MPI_COMM_WORLD, MPI_GROUP_WORLD, ierr)


    ALLOCATE(   Workers_Array(0:NUM_PROCS-1)      )
    Workers_Array = (/(i, i=0,NUM_PROCS-1, 1)/)

    CALL MPI_Group_incl(    MPI_GROUP_WORLD,            &
                            NUM_PROCS,                  &
                            Workers_Array,              &
                            CHIMERA_GROUP_WORLD,        &
                            ierr                        )


    CALL MPI_Comm_create_group( MPI_COMM_WORLD,             &
                                CHIMERA_GROUP_WORLD,        &
                                0,                          &
                                CHIMERA_COMM_WORLD,         &
                                ierr                        )

    IF ( MPI_COMM_NULL .NE. CHIMERA_COMM_WORLD ) THEN

        CALL MPI_COMM_RANK(CHIMERA_COMM_WORLD, myID_Chimera, ierr)
    ELSE

        myID_CHIMERA = -1

    END IF


ELSE IF ( NUM_PROCS > nPROCS ) THEN

    PRINT*," CHIMERA expects more processes than in MPI_COMM_WORLD "
    PRINT*," Expected : ",NUM_PROCS," Found : ", nPROCS

    RUN_POSEIDON_FLAG = .FALSE.

END IF

IF( myID == 0) THEN

    PRINT*,"  "
    PRINT*,"AFTER CREATE CHIMERA_COMM_WORLD"
    PRINT*,"  "

END IF
CALL SLEEP(1)
CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)


!IF (nPROCS .NE. Num_Procs ) THEN
!    IF (myID == 0 ) THEN
!       PRINT*,"Not the expected number of CHIMERA processes, may cause errors"
!       PRINT*," Expected : ",Num_Procs," Found : ",nPROCS
!       PRINT*,"Checked in d.main_3D.F90, line 360"
!    END IF
!    RUN_POSEIDON_FLAG = .FALSE.
!
!END IF

IF ( Num_y_Procs * Num_Z_procs .NE. NUM_PROCS ) THEN
    IF (myID == 0) THEN
       PRINT*,"Grid will not divide evenly onto processes"
       PRINT*,"Checked in d.main_3D.F90, line 370"
    END IF
    RUN_POSEIDON_FLAG = .FALSE.

END IF






!                                   !
!   Create Space for Local Rho Data !
!                                   !
IF ( Problem_Dimension .EQ. 1 ) THEN

    ! 1-D Set Up

    ij_ray_dim = 1
    ik_ray_dim = 1

    myid_theta = 0
    myid_phi = 0



ELSE IF ( Problem_Dimension .EQ. 2 ) THEN



    !   2-D Set up

    ij_ray_dim = ny/Num_y_Procs
    ik_ray_dim = nz

    myid_theta = myid_CHIMERA
    myid_phi = 0



ELSE IF ( Problem_Dimension .EQ. 3 ) THEN

    !   3-D Set up

    ij_ray_dim = ny/Num_y_Procs
    ik_ray_dim = nz/Num_z_Procs

    j_block_col = MOD(myid_CHIMERA, Num_y_PROCS)
    k_block_row = myid_CHIMERA/Num_y_PROCS


    CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, k_block_row, myid_CHIMERA, MPI_COMM_XY, ierr)
    CALL MPI_COMM_RANK(MPI_COMM_XY, myid_theta, ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_XY, num_y_procs, ierr)


    CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, j_block_col, myid_CHIMERA, MPI_COMM_XZ, ierr)
    CALL MPI_COMM_RANK(MPI_COMM_XZ, myid_phi, ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_XZ, num_z_procs, ierr)


END IF


!  Create Communicators if Ying-Yang Grid is used
!  Not complete
IF ( ngrid == 2 ) THEN




ELSE

    MPI_COMM_GRID = MPI_COMM_WORLD

END IF



IF( myID == -130) THEN

     PRINT*,"ny",ny,"nz",nz
     PRINT*,"Num_Y_PROCS",Num_y_Procs,"Num_z_Procs",Num_z_Procs
     PRINT*,"ij_ray_dim",ij_ray_dim,"ik_ray_dim",ik_ray_dim
END IF










!
!   RECREATE CHIMERA DOMAIN DECOMPOSITION
!

ALLOCATE( x_e(0:nx), x_c(1:nx), y_e(-5:ny+7), y_c(-5:ny+6), z_e(-5:ny+7), z_c(-5:ny+6) )
ALLOCATE(dx_c(1:nx))
ALLOCATE(dy_c(-5:ny+6), tmp_dy_c(1:ny))
ALLOCATE(dz_c(-5:nz+6), tmp_dz_c(1:nz))

!                                               !
!   For this example, we build uniform meshes.  !
!                                               !

IF ( TEST_NUM == 2 ) THEN

    PRINT*,"Before Load_CHIMERA_HDF5",myID_chimera
    CALL Load_CHIMERA_HDF5(file_number, stride, path, frame_flag)
    PRINT*,"After  Load_CHIMERA_HDF5",myID_chimera

    dx_c(0:nx) = CHIMERA_Delta_R(0:nx)
    dy_c(0:ny) = CHIMERA_Delta_T(0:ny)
    dz_c(0:nz) = CHIMERA_Delta_P(0:nz)


    x_e(0:nx) = CHIMERA_R_LOCS(0:nx)
    y_e(0:ny) = CHIMERA_T_LOCS(0:ny)
    z_e(0:nz) = CHIMERA_P_LOCS(0:nz)

ELSE

    CALL Create_3D_Mesh( Mesh_Type,                                    &
                         Inner_Radius, Core_Radius, Outer_Radius,      &
                         nx, nc, ny, nz,                               &
                         x_e, x_c, dx_c,                               &
                         y_e(0:ny), y_c(1:ny), dy_c(1:ny),             &
                         z_e(0:nz), z_c(1:nz), dz_c(1:nz)              )

    OPEN( UNIT = 42, file = 'OUTPUT/rmesh.out')
    WRITE(42,*) x_e
    CLOSE( UNIT = 42)

END IF





!
!   RECREATE CHIMERA LOCAL DATA
!


!                                                                       !
!   For this example the quadrature points in each dimension are the    !
!   Lobatto quadrature points.                                          !
!                                                                       !
Input_R_Quad = Initialize_LG_Quadrature_Locations(Num_Input_Nodes(1))
Input_T_Quad = Initialize_LG_Quadrature_Locations(Num_Input_Nodes(2))
Input_P_Quad = Initialize_LG_Quadrature_Locations(Num_Input_Nodes(3))




!                                                                       !
!   As the Lobatto points are defined in [-1,1] space, but we want to   !
!   provide source input on a [-.5, .5] space, we need to map the       !
!   the locations from one space to the other.                          !
!                                                                       !

!   Function - Map_From_X_Space - This Function maps Input_R_Quad (Real Double Value/Vector in [-1,1] space) to
!                                   a space defined by the limits (Real Double Values)
Input_R_Quad = Map_From_X_Space(Left_Limit, Right_Limit, Input_R_Quad)
Input_T_Quad = Map_From_X_Space(Left_Limit, Right_Limit, Input_T_Quad)
Input_P_Quad = Map_From_X_Space(Left_Limit, Right_Limit, Input_P_Quad)




ALLOCATE(   Local_E(1:Num_Input_Nodes(1)*Num_Input_Nodes(2)*Num_Input_Nodes(3),       &
                    0:nx-1, 0:ij_ray_dim-1, 0:ik_ray_dim-1    ),                  &
            Local_S(1:Num_Input_Nodes(1)*Num_Input_Nodes(2)*Num_Input_Nodes(3),       &
                    0:nx-1, 0:ij_ray_dim-1, 0:ik_ray_dim-1     ),                 &
            Local_Si(1:Num_Input_Nodes(1)*Num_Input_Nodes(2)*Num_Input_Nodes(3),      &
                    0:nx-1, 0:ij_ray_dim-1, 0:ik_ray_dim-1, 1:3)                  )





CALL Poseidon_Initialize_CFA_Test_Prob_Par_CHIMERA(Test_Num, nx, ny, nz,                &
                                                    nx, ij_ray_dim, ik_ray_dim,         &
                                                    Num_Input_Nodes,                    &
                                                    Input_R_Quad,                       &
                                                    dx_c, x_e,                          &
                                                    Local_E, Local_S, Local_Si          )





!DO i = 0,nPROCS
!
!    if (myid == i ) THEN
!
!        PRINT*,myid,"Local_E"
!        PRINT*,Local_E
!        PRINT*,"++++++++++++++++"
!
!    END IF
!
!    CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
!
!!    PRINT*,"after barrier",myid
!END DO







!
!   Make Poseidon Calls
!


!RUN_POSEIDON_FLAG = .FALSE.

IF ( RUN_POSEIDON_FLAG .EQV. .TRUE. ) THEN


    IF ( myid == 0 ) THEN

        PRINT*,"----------  CHIMERA Variables ----------"
        PRINT*,"        CHIMERA_Y_nPROCS = ",CHIMERA_y_PROCS
        PRINT*,"        CHIMERA_z_nPROCS = ",CHIMERA_z_PROCS

    END IF

!   PRINT*,"Left & Right Limits",Left_Limit, Right_Limit,CHIMERA_LEFT_LIMIT,CHIMERA_RIGHT_LIMIT


    mode = 1
    modeb = 0
    imin = 1
    imax = nx
    CALL Poseidon_CFA_3D(   mode, modeb, imin, imax, nx,                    &
                            ij_ray_dim, ik_ray_dim, ny, nz,                 &
                            x_e, x_c, dx_c, y_e, y_c, dy_c, z_e, dz_c,      &
                            Num_Input_Nodes,                                &
                            Input_R_Quad, Input_T_Quad, Input_P_Quad,       &
                            Left_Limit, Right_Limit,                        &
                            Local_E, Local_S, Local_Si                      )




END IF









CALL MPI_Finalize(ierr)








END PROGRAM Poseidon_CFA
