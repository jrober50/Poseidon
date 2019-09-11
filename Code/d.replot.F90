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

USE Poseidon_Parameter_Read_Module, &
            ONLY :  UNPACK_POSEIDON_PARAMETERS,             &
                    WRITE_CFA_COEFFICIENTS,                 &
                    READ_CFA_COEFFICIENTS
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
                    CHIMERA_POTENTIAL,      &
                    CHIMERA_SHIFT_VAL,      &
                    CHIMERA_LAPSE_VAL,      &
                    CHIMERA_E,              &
                    CHIMERA_S,              &
                    CHIMERA_Si,             &
                    CHIMERA_TEST_NUMBER,    &
                    CHIMERA_FRAME,          &
                    CHIMERA_START_FRAME,    &
                    CHIMERA_MAX_FRAME,      &
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


IMPLICIT NONE







CALL Unpack_CHIMERA_Parameters()


Left_Limit          = CHIMERA_LEFT_LIMIT
Right_Limit         = CHIMERA_RIGHT_LIMIT
Num_Procs           = CHIMERA_PROCS
Num_y_Procs         = CHIMERA_y_PROCS
Num_z_Procs         = CHIMERA_z_PROCS
Problem_Dimension   = CHIMERA_DIMENSION
Inner_Radius        = CHIMERA_INNER_RADIUS
Core_Radius         = CHIMERA_CORE_RADIUS
Outer_Radius        = CHIMERA_OUTER_RADIUS
nx                  = CHIMERA_R_ELEMS
nc                  = CHIMERA_C_ELEMS
ny                  = CHIMERA_T_ELEMS
nz                  = CHIMERA_P_ELEMS
MESH_TYPE           = CHIMERA_MESH_TYPE
TEST_NUM            = CHIMERA_TEST_NUMBER



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


END IF




ALLOCATE( x_e(0:nx), x_c(1:nx), y_e(-5:ny+7), y_c(-5:ny+6), z_e(-5:ny+7), z_c(-5:ny+6) )
ALLOCATE(dx_c(1:nx))
ALLOCATE(dy_c(-5:ny+6), tmp_dy_c(1:ny))
ALLOCATE(dz_c(-5:nz+6), tmp_dz_c(1:nz))


IF ( TEST_NUM == 2 ) THEN

    WRITE(*,'(A16,I3.3,A4,I3.3)')"Current Frame = ",CHIMERA_FRAME," of ",CHIMERA_MAX_FRAME
    CALL Load_CHIMERA_HDF5(CHIMERA_FRAME, stride, path, frame_flag)

    dx_c(1:nx) = CHIMERA_Delta_R(0:nx-1)
    dy_c(1:ny) = CHIMERA_Delta_T(0:ny-1)
    dz_c(1:nz) = CHIMERA_Delta_P(0:nz-1)


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

    OPEN( UNIT = 42, file = 'OUTPUT/R_Mesh.out')
    WRITE(42,*) x_e
    CLOSE( UNIT = 42)

END IF




CALL UNPACK_POSEIDON_PARAMETERS()

Mode = 1


CALL Poseidon_Initialize(   mode,                           & ! mode
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
                            dz_c(1:nz)                      )

CALL READ_CFA_COEFFICIENTS()









END PROGRAM Poseidon_CFA


