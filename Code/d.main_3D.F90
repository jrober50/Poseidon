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

USE DRIVER_Parameters, &
            ONLY :  DRIVER_R_ELEMS,        &
                    DRIVER_C_ELEMS,        &
                    DRIVER_T_ELEMS,        &
                    DRIVER_P_ELEMS,        &
                    DRIVER_R_INPUT_NODES,  &
                    DRIVER_T_INPUT_NODES,  &
                    DRIVER_P_INPUT_NODES,  &
                    DRIVER_LEFT_LIMIT,     &
                    DRIVER_RIGHT_LIMIT,    &
                    DRIVER_INNER_RADIUS,   &
                    DRIVER_CORE_RADIUS,    &
                    DRIVER_OUTER_RADIUS,   &
                    DRIVER_MESH_TYPE,      &
                    DRIVER_R_LOCS,         &
                    DRIVER_T_LOCS,         &
                    DRIVER_P_LOCS,         &
                    DRIVER_Delta_R,        &
                    DRIVER_Delta_T,        &
                    DRIVER_Delta_P,        &
                    ENCLOSED_MASS,          &
                    DRIVER_POTENTIAL,      &
                    DRIVER_SHIFT_VAL,      &
                    DRIVER_LAPSE_VAL,      &
                    DRIVER_E,              &
                    DRIVER_S,              &
                    DRIVER_Si,             &
                    DRIVER_TEST_NUMBER,    &
                    DRIVER_FRAME,          &
                    DRIVER_START_FRAME,    &
                    DRIVER_END_FRAME,      &
                    DRIVER_TOTAL_FRAMES,    &
                    DRIVER_DIMENSION,      &
                    DRIVER_PROCS,          &
                    DRIVER_y_PROCS,        &
                    DRIVER_z_PROCS,        &
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
                    SELFSIM_START_T,        &
                    SELFSIM_END_T,          &
                    SELFSIM_NUM_FRAMES,     &
                    SELFSIM_KAPPA,          &
                    SELFSIM_GAMMA,          &
                    CHIMERA_START_FRAME,    &
                    CHIMERA_END_FRAME,      &
                    Iteration_History,      &
                    DRIVER_RUN_TIME_TABLE,  &
                    DRIVER_FRAME_TIME_TABLE

USE Driver_Params_Read_Module, &
            ONLY :  Unpack_Driver_Parameters


USE CHIMERA_HDF5_Module, &
            ONLY :  Load_CHIMERA_HDF5


USE Mesh_Module, &
            ONLY :  Create_3D_Mesh


USE Units_Module, &
            ONLY :  Set_Units, Grav_Constant_G, Speed_of_Light



USE Driver_Test_Functions_Module, &
            ONLY :  Poseidon_Initialize_CFA_Test_Problem_CHIMERA


USE Driver_IO_Functions_Module, &
            ONLY :  Open_Run_Report_File,                           &
                    Output_Run_Report,                              &
                    Close_Run_Report_File,                          &
                    Open_Frame_Report_File,                         &
                    Output_Frame_Report,                            &
                    Close_Frame_Report_File,                        &
                    Output_Iteration_History



USE Poseidon_Interface, &
            ONLY :  Poseidon_CFA_3D,                                &
                    Poseidon_Shutdown

USE Poseidon_Info_Module, &
            ONLY :  PQ_Iterations_Max,                              &
                    PQ_Iterations_Used


USE Driver_Additional_Functions_Module, &
            ONLY :  Map_From_X_Space,                       &
                    Initialize_LG_Quadrature_Locations

USE Poseidon_Parameters, &
            ONLY :  Poseidon_Frame


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

INTEGER                                                     ::  DRIVER_COMM_WORLD,     &
                                                                DRIVER_GROUP_WORLD,    &
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


INTEGER                         :: Test_Num
INTEGER                         :: Mesh_Type
INTEGER                         :: FIRST_FRAME_FLAG
INTEGER                         :: CHIMERA_FRAME
REAL(KIND = idp)                :: SELFSIM_TIME
REAL(KIND = idp)                :: SELFSIM_DELTA_T

INTEGER                         :: file_number=00001
INTEGER                         :: stride=1

LOGICAL                         :: frame_flag=.TRUE.

CHARACTER(len=12)               :: path='Data3/Frames'

CALL MPI_INIT(ierr)

CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nPROCS, ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)



!!                                      !!
!!      CHIMERA Setup Parameters        !!
!!                                      !!

CALL Unpack_DRIVER_Parameters()




ALLOCATE( DRIVER_R_LOCS(0:DRIVER_R_ELEMS),  &
          DRIVER_T_LOCS(0:DRIVER_T_ELEMS),  &
          DRIVER_P_LOCS(0:DRIVER_P_ELEMS)   )

ALLOCATE( DRIVER_Delta_R(0:DRIVER_R_ELEMS-1),   &
          DRIVER_Delta_T(0:DRIVER_T_ELEMS-1),   &
          DRIVER_DELTA_P(0:DRIVER_P_ELEMS-1)    )

ALLOCATE( Enclosed_Mass(0:DRIVER_R_ELEMS+1) )
ALLOCATE( DRIVER_Potential(0:DRIVER_R_ELEMS) )
ALLOCATE( DRIVER_SHIFT_VAL(0:DRIVER_R_ELEMS) )



Num_Input_Nodes(1) = DRIVER_R_INPUT_NODES
Num_Input_Nodes(2) = DRIVER_T_INPUT_NODES
Num_Input_Nodes(3) = DRIVER_P_INPUT_NODES
Num_DOF = Num_Input_Nodes(1)*Num_Input_Nodes(2)*Num_Input_Nodes(3)


ALLOCATE( DRIVER_E(1:Num_DOF, 0:DRIVER_R_ELEMS-1,0:DRIVER_T_ELEMS-1,0:DRIVER_P_ELEMS-1) )
ALLOCATE( DRIVER_S(1:Num_DOF, 0:DRIVER_R_ELEMS-1,0:DRIVER_T_ELEMS-1,0:DRIVER_P_ELEMS-1) )
ALLOCATE( DRIVER_Si(1:Num_DOF, 0:DRIVER_R_ELEMS-1,0:DRIVER_T_ELEMS-1,0:DRIVER_P_ELEMS-1, 1:3) )


ALLOCATE(   Input_R_Quad(1:Num_Input_Nodes(1)),         &
            Input_T_Quad(1:Num_Input_Nodes(2)),         &
            Input_P_Quad(1:Num_Input_Nodes(3))          )

ALLOCATE( Driver_Frame_Time_Table(1:25) )
ALLOCATE( Driver_Run_Time_Table(1:25) )




!RUN_POSEIDON_FLAG = .FALSE.

CALL Set_Units( "C" )





!!                                       !!
!!   Poseidon Initialization Variables   !!
!!                                       !!
Left_Limit          = DRIVER_LEFT_LIMIT
Right_Limit         = DRIVER_RIGHT_LIMIT
Num_Procs           = DRIVER_PROCS
Num_y_Procs         = DRIVER_y_PROCS
Num_z_Procs         = DRIVER_z_PROCS
Problem_Dimension   = DRIVER_DIMENSION
Inner_Radius        = DRIVER_INNER_RADIUS
Core_Radius         = DRIVER_CORE_RADIUS
Outer_Radius        = DRIVER_OUTER_RADIUS
nx                  = DRIVER_R_ELEMS
nc                  = DRIVER_C_ELEMS
ny                  = DRIVER_T_ELEMS
nz                  = DRIVER_P_ELEMS
MESH_TYPE           = DRIVER_MESH_TYPE
TEST_NUM            = DRIVER_TEST_NUMBER






!                                       !
!       Initialize the MPI World        !
!                                       !

IF ( NUM_PROCS == nPROCS ) THEN

    CALL MPI_COMM_DUP(MPI_COMM_WORLD, DRIVER_COMM_WORLD, ierr)
    CALL MPI_COMM_RANK(DRIVER_COMM_WORLD, myID_Chimera, ierr)

ELSE IF ( NUM_PROCS < nPROCs ) THEN



    CALL MPI_Comm_group(MPI_COMM_WORLD, MPI_GROUP_WORLD, ierr)


    ALLOCATE(   Workers_Array(0:NUM_PROCS-1)      )
    Workers_Array = (/(i, i=0,NUM_PROCS-1, 1)/)

    CALL MPI_Group_incl(    MPI_GROUP_WORLD,            &
                            NUM_PROCS,                  &
                            Workers_Array,              &
                            DRIVER_GROUP_WORLD,        &
                            ierr                        )


    CALL MPI_Comm_create_group( MPI_COMM_WORLD,             &
                                DRIVER_GROUP_WORLD,        &
                                0,                          &
                                DRIVER_COMM_WORLD,         &
                                ierr                        )

    IF ( MPI_COMM_NULL .NE. DRIVER_COMM_WORLD ) THEN

        CALL MPI_COMM_RANK(DRIVER_COMM_WORLD, myID_Chimera, ierr)
    ELSE

        ! If process doesn't belong give bogus id.
        myID_CHIMERA = -1

    END IF


ELSE IF ( NUM_PROCS > nPROCS ) THEN

    PRINT*," CHIMERA expects more processes than in MPI_COMM_WORLD "
    PRINT*," Expected : ",NUM_PROCS," Found : ", nPROCS

    RUN_POSEIDON_FLAG = .FALSE.

END IF



! Check if there are enough processes to fill process grid
IF ( Num_y_Procs * Num_Z_procs .NE. NUM_PROCS ) THEN
    IF (myID == 0) THEN
       PRINT*,"Grid will not divide evenly onto processes"
       PRINT*,"Checked in d.main_3D.F90, line 445"
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
    PRINT*,"In Driver, myID_Theta",myID_Theta

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











!
!   RECREATE CHIMERA DOMAIN DECOMPOSITION
!

ALLOCATE( x_e(0:nx), x_c(1:nx), y_e(-5:ny+7), y_c(-5:ny+6), z_e(-5:ny+7), z_c(-5:ny+6) )
ALLOCATE( dx_c(1:nx))
ALLOCATE( dy_c(-5:ny+6), tmp_dy_c(1:ny))
ALLOCATE( dz_c(-5:nz+6), tmp_dz_c(1:nz))

ALLOCATE(   Local_E(1:Num_DOF, 0:nx-1, 0:ij_ray_dim-1, 0:ik_ray_dim-1    ),             &
            Local_S(1:Num_DOF, 0:nx-1, 0:ij_ray_dim-1, 0:ik_ray_dim-1     ),            &
            Local_Si(1:Num_DOF, 0:nx-1, 0:ij_ray_dim-1, 0:ik_ray_dim-1, 1:3)            )

ALLOCATE( DRIVER_LAPSE_VAL(0:nx-1, 0:ij_ray_dim-1, 0:ik_ray_dim-1) )






!
!   Set Up Frame Iterations
!
IF ( TEST_NUM == 2 ) THEN

    DRIVER_START_FRAME = 1
    DRIVER_END_FRAME = CHIMERA_END_FRAME - CHIMERA_START_FRAME + 1

ELSE IF ( TEST_NUM == 3 ) THEN
    DRIVER_START_FRAME = 1
    DRIVER_END_FRAME = SELFSIM_NUM_FRAMES
    IF ( SELFSIM_NUM_FRAMES .NE. 1 ) THEN
        SELFSIM_DELTA_T = (SELFSIM_END_T - SELFSIM_START_T)/(SELFSIM_NUM_FRAMES-1)
    ELSE
        SELFSIM_DELTA_T = SELFSIM_START_T
    END IF

END IF
FIRST_FRAME_FLAG = 1


IF ( myID == 0 ) THEN
    CALL Open_Run_Report_File()
END IF

DRIVER_TOTAL_FRAMES = DRIVER_END_FRAME - DRIVER_START_FRAME + 1
ALLOCATE( Iteration_History(1:DRIVER_TOTAL_FRAMES) )








DO DRIVER_FRAME = DRIVER_START_FRAME,DRIVER_END_FRAME
    Poseidon_Frame = Driver_Frame
    CALL OPEN_FRAME_REPORT_FILE(DRIVER_FRAME)


    WRITE(*,'(A16,I3.3,A4,I3.3)')"Current Frame = ",DRIVER_FRAME," of ",DRIVER_END_FRAME
    IF ( TEST_NUM == 3 ) THEN
        !  Yahil Self-Similar Test, Calculate Time Parameter
        SELFSIM_T = (DRIVER_FRAME - 1 )*SELFSIM_DELTA_T + SELFSIM_START_T
        WRITE(*,'(A13,ES17.10,A4)')"Yahil Time = ",SELFSIM_T," ms."
    END IF





    IF ( TEST_NUM == 2 ) THEN
        !  CHIMERA Test, Calculate Current CHIMERA Frame

        CHIMERA_FRAME = CHIMERA_START_FRAME + DRIVER_FRAME - 1
        WRITE(*,'(A28,I3.3)')"CHIMERA HDF5 Frame number = ",CHIMERA_FRAME


        CALL Load_CHIMERA_HDF5(CHIMERA_FRAME, stride, path, frame_flag)

        dx_c(1:nx) = DRIVER_Delta_R(0:nx-1)
        dy_c(1:ny) = DRIVER_Delta_T(0:ny-1)
        dz_c(1:nz) = DRIVER_Delta_P(0:nz-1)


        x_e(0:nx) = DRIVER_R_LOCS(0:nx)
        y_e(0:ny) = DRIVER_T_LOCS(0:ny)
        z_e(0:nz) = DRIVER_P_LOCS(0:nz)



    ELSE

        ! Other Test, Create Mesh
        CALL Create_3D_Mesh( Mesh_Type,                                    &
                             Inner_Radius, Core_Radius, Outer_Radius,      &
                             nx, nc, ny, nz,                               &
                             x_e, x_c, dx_c,                               &
                             y_e(0:ny), y_c(1:ny), dy_c(1:ny),             &
                             z_e(0:nz), z_c(1:nz), dz_c(1:nz)              )


    END IF




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











    !
    !   Initialize Source Variables Based on Test Problem
    !
    CALL Poseidon_Initialize_CFA_Test_Problem_CHIMERA(  Test_Num, nx, ny, nz,               &
                                                        nx, ij_ray_dim, ik_ray_dim,         &
                                                        Num_Input_Nodes,                    &
                                                        Input_R_Quad, Input_T_Quad,         &
                                                        dx_c, x_e, y_e(0:ny),               &
                                                        Local_E, Local_S, Local_Si          )


    !
    !   Call Poseidon Interface
    !

    modeb = 0
    imin = 1
    imax = nx
    !   CHOOSE OPERATIONAL MODE
    !
    !   Mode = 0    :   Initalize, Run
    !   Mode = 1    :              Run
    !   Mode = 2    :              Run, Close
    !   Mode = 3    :   Initalize, Run, Close


    IF ( FIRST_FRAME_FLAG == 1 ) THEN
        mode = 0
        modeb = 0
    ELSE IF ( DRIVER_FRAME == DRIVER_END_FRAME) THEN
        mode = 2
        modeb = 1
    ELSE IF ( ( DRIVER_FRAME == DRIVER_END_FRAME) .AND. ( FIRST_FRAME_FLAG == 1 ) ) THEN
        mode = 3
        modeb = 0
    ELSE
        mode = 1
        modeb = 1
    END IF

    CALL Poseidon_CFA_3D(   mode, modeb, imin, imax, nx,                    &
                            ij_ray_dim, ik_ray_dim, ny, nz,                 &
                            x_e, x_c, dx_c, y_e, y_c, dy_c, z_e, dz_c,      &
                            Num_Input_Nodes,                                &
                            Input_R_Quad, Input_T_Quad, Input_P_Quad,       &
                            Left_Limit, Right_Limit,                        &
                            Local_E, Local_S, Local_Si                      )




    CALL PQ_ITERATIONS_USED( Iteration_History(DRIVER_FRAME) )

    CALL OUTPUT_FRAME_REPORT(DRIVER_FRAME)
    CALL CLOSE_FRAME_REPORT_FILE()
    FIRST_FRAME_FLAG = 0
END DO ! Frame Loop


IF ( myID == 0 ) THEN
    CALL OUTPUT_RUN_REPORT()
    CALL OUTPUT_ITERATION_HISTORY()
    CALL CLOSE_RUN_REPORT_FILE()
END IF
CALL Poseidon_Shutdown()
CALL MPI_Finalize(ierr)


WRITE(*,'(//A18//)')"DING! You're Done!"





END PROGRAM Poseidon_CFA

