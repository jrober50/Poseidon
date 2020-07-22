   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE CHIMERA_HDF5_Module                                                          !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
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
USE Poseidon_Constants_Module, &
                        ONLY :  idp,pi,eps
USE Units_Module, &
                        ONLY :  Grav_Constant_G,        &
                                Speed_of_Light,         &
                                C_Square,               &
                                GR_Source_Scalar

USE Driver_Parameters, &
                        ONLY :  DRIVER_R_ELEMS,        &
                                DRIVER_T_ELEMS,        &
                                DRIVER_P_ELEMS,        &
                                DRIVER_R_LOCS,         &
                                DRIVER_T_LOCS,         &
                                DRIVER_P_LOCS,         &
                                DRIVER_R_INPUT_NODES,  &
                                DRIVER_T_INPUT_NODES,  &
                                DRIVER_P_INPUT_NODES,  &
                                DRIVER_Delta_R,        &
                                DRIVER_Delta_T,        &
                                DRIVER_Delta_P,        &
                                Enclosed_Mass,          &
                                DRIVER_POTENTIAL,      &
                                DRIVER_LAPSE_VAL,      &
                                SOURCE_OUTPUT_FLAG,    &
                                DRIVER_E,              &
                                DRIVER_S,              &
                                DRIVER_Si,             &
                                DRIVER_T_LOCS_LOCAL,   &
                                DRIVER_P_LOCS_LOCAL,   &
                                myID,                   &
                                myID_theta,             &
                                myID_phi,               &
                                ij_ray_dim,             &
                                ik_ray_dim,             &
                                MPI_COMM_GRID,          &
                                DRIVER_PROCS,          &
                                DRIVER_y_PROCS,        &
                                DRIVER_z_PROCS

USE CHIMERA_TEST_FUNCS_Module,  &
                        ONLY :  INIT_CHIMERA_POTENTIAL, &
                                INIT_CHIMERA_SHIFT_VAL

USE HDF5_IO_Module
USE MPI
USE HDF5

IMPLICIT NONE

INTEGER                                 :: my_grid = 1
INTEGER                                 :: nlog=6


INTEGER                                 :: imin, imax
INTEGER                                 :: jmin, jmax
INTEGER                                 :: kmin, kmax

INTEGER                                 :: i_lower, i_upper
INTEGER                                 :: j_lower, j_upper
INTEGER                                 :: k_lower, k_upper


INTEGER                                 :: nz_hyperslabs
INTEGER, PRIVATE                        :: j_offset
INTEGER, PRIVATE                        :: k_offset
INTEGER, PRIVATE                        :: nz_hyperslab_width
INTEGER, PRIVATE                        :: nz_hyperslab_procwidth
INTEGER, PRIVATE                        :: my_hyperslab_group
INTEGER, PRIVATE                        :: k_hyperslab_offset
INTEGER, PRIVATE                        :: mpi_comm_per_hyperslab_group
INTEGER, PRIVATE                        :: myid_hyperslab_group

LOGICAL                                 :: hyperslab_group_master

CHARACTER(len=22)                       :: suffix
CHARACTER(len= 8), PARAMETER            :: prefix = 'chimera_'

INTEGER(HID_T), PUBLIC                  :: file_id
INTEGER(HID_T), PUBLIC                  :: group_id

INTEGER(HSIZE_T), DIMENSION(1)          :: datasize1d
INTEGER(HSIZE_T), DIMENSION(2)          :: datasize2d
INTEGER(HSIZE_T), DIMENSION(3)          :: datasize3d
INTEGER(HSIZE_T), DIMENSION(4)          :: datasize4d
INTEGER(HSIZE_T), DIMENSION(5)          :: datasize5d

INTEGER(HSIZE_T), DIMENSION(2)          :: slab_offset2d
INTEGER(HSIZE_T), DIMENSION(3)          :: slab_offset3d
INTEGER(HSIZE_T), DIMENSION(4)          :: slab_offset4d
INTEGER(HSIZE_T), DIMENSION(5)          :: slab_offset5d


REAL(KIND = idp), DIMENSION(:,:,:), ALLOCATABLE     ::  Density
REAL(KIND = idp), DIMENSION(:,:,:), ALLOCATABLE     ::  Pressure
REAL(KIND = idp), DIMENSION(:,:,:), ALLOCATABLE     ::  vel_r
REAL(KIND = idp), DIMENSION(:,:,:), ALLOCATABLE     ::  vel_t
REAL(KIND = idp), DIMENSION(:,:,:), ALLOCATABLE     ::  vel_p
REAL(KIND = idp), DIMENSION(:,:,:), ALLOCATABLE     ::  int_e
REAL(KIND = idp), DIMENSION(:,:,:), ALLOCATABLE     ::  Lapse_Val
REAL(KIND = idp), DIMENSION(1)                                  :: time


CONTAINS
!+201+###########################################################################!
!                                                                                !
!               LOAD_CHIMERA_HDF5                                                !
!                                                                                !
!################################################################################!
SUBROUTINE LOAD_CHIMERA_HDF5(file_number, stride, path, frame_flag )

INTEGER, INTENT(IN)                                 :: file_number
INTEGER, INTENT(IN)                                 :: stride

LOGICAL, INTENT(IN)                                 :: frame_flag

CHARACTER(len=*), INTENT(IN)                        :: path

INTEGER                                             :: here
INTEGER                                             :: i,j,k
INTEGER                                             :: pd, td, rd

REAL(KIND = idp)                                    :: Specific_Enthalpy
REAL(KIND = idp)                                    :: V_Square
REAL(KIND = idp)                                    :: Reusable_Value

INTEGER, DIMENSION(1:3)                             ::  Num_Nodes

REAL(KIND = idp), DIMENSION(:,:,:,:), ALLOCATABLE   ::  HR_Density

REAL(KIND = idp), DIMENSION(:), ALLOCATABLE         ::  Quad_Locs
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE         ::  Quad_Weights

REAL(KIND = idp)                                    ::  R_Center
REAL(KIND = idp)                                    ::  deltar_overtwo
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE         ::  CUR_R_LOCS

INTEGER :: error







imin = 1
imax = DRIVER_R_ELEMS

jmin = 1
jmax = DRIVER_T_ELEMS

kmin = 1
kmax = DRIVER_P_ELEMS

i_lower = imin
i_upper = imax

j_lower = jmin
j_upper = jmax

k_lower = kmin
k_upper = kmax

IF ( DRIVER_T_ELEMS > 1) j_lower = 1
IF ( DRIVER_T_ELEMS > 1) j_upper = ij_ray_dim
IF ( DRIVER_P_ELEMS > 1) k_lower = 1
IF ( DRIVER_P_ELEMS > 1) k_upper = ik_ray_dim


ALLOCATE( HR_DENSITY(DRIVER_R_INPUT_NODES,k_lower:k_upper,j_lower:j_upper,i_lower:i_upper))
ALLOCATE( Quad_Locs(1:DRIVER_R_INPUT_NODES) )
ALLOCATE( Quad_Weights(1:DRIVER_R_INPUT_NODES) )
ALLOCATE( Cur_R_Locs(1:DRIVER_R_INPUT_NODES) )



!  Allocate space for temporary storage !
ALLOCATE( Density(k_lower:k_upper,j_lower:j_upper,i_lower:i_upper) )
ALLOCATE( Pressure(k_lower:k_upper,j_lower:j_upper,i_lower:i_upper) )
ALLOCATE( vel_r(k_lower:k_upper,j_lower:j_upper,i_lower:i_upper) ) 
ALLOCATE( vel_t(k_lower:k_upper,j_lower:j_upper,i_lower:i_upper) )
ALLOCATE( vel_p(k_lower:k_upper,j_lower:j_upper,i_lower:i_upper) )
ALLOCATE( int_e(k_lower:k_upper,j_lower:j_upper,i_lower:i_upper) )
ALLOCATE( Lapse_Val(k_lower:k_upper,j_lower:j_upper,i_lower:i_upper) )


CALL Initialize_LG_Quadrature(DRIVER_R_INPUT_NODES, Quad_Locs, Quad_Weights)


!  Open CHIMERA HDF5 File  !
CALL OPEN_CHIMERA_HDF5( file_number, path, frame_flag )

!  Read Relevant data into arrays !
CALL READ_CHIMERA_HDF5( file_number, stride )

!  Close CHIMERA HDF5 File  !
CALL CLOSE_CHIMERA_HDF5()





Density(:,:,DRIVER_R_ELEMS) = Density(:,:,DRIVER_R_ELEMS-1)


!DO i = i_lower, i_upper-1
!
!    r_center = (DRIVER_R_LOCS(i)+DRIVER_R_LOCS(i+1))/2.0_idp
!    deltar_overtwo = (DRIVER_R_LOCS(i+1)-DRIVER_R_LOCS(i))/2.0_idp
!    CUR_R_LOCS(:) = deltar_overtwo * (Quad_Locs(:) + 1.0_idp ) + DRIVER_R_LOCS(i)
!
!    DO rd = 1,DRIVER_R_INPUT_NODES
!
!        IF ( CUR_R_LOCS(rd) < r_center  ) THEN
!
!
!        ELSE
!
!
!        END IF
!    END DO  ! rd Loop
!END DO  ! i Loop







!  Convert CHIMERA data into CFA Input !
DO i = i_lower, i_upper
    DO j = j_lower, j_upper
        DO k = k_lower, k_upper

            Specific_Enthalpy = C_Square + int_e(k,j,i) + Pressure(k,j,i)/Density(k,j,i)
            V_Square = vel_r(k,j,i)*vel_r(k,j,i) + vel_t(k,j,i)*vel_t(k,j,i) + vel_p(k,j,i)*vel_p(k,j,i)

            Reusable_Value = (Density(k,j,i)*Specific_Enthalpy)/(1.0_idp - V_Square/C_Square)

            DO pd = 1, DRIVER_P_INPUT_NODES
                DO td = 1, DRIVER_T_INPUT_NODES
                    DO rd = 1, DRIVER_R_INPUT_NODES



                        here = (pd-1) * DRIVER_R_INPUT_NODES*DRIVER_T_INPUT_NODES             &
                             + (td-1) * DRIVER_R_INPUT_NODES                                   &
                             + rd

                        DRIVER_E(here,i-1,j-1,k-1) = Reusable_Value - Pressure(k,j,i)
                        DRIVER_S(here,i-1,j-1,k-1) = Reusable_Value*V_Square/C_Square + 3.0_idp * Pressure(k,j,i)
                        DRIVER_Si(here,i-1,j-1,k-1,1) = Reusable_Value*vel_r(k,j,i)/C_Square
                        DRIVER_Si(here,i-1,j-1,k-1,2) = Reusable_Value*vel_t(k,j,i)/C_Square
                        DRIVER_Si(here,i-1,j-1,k-1,3) = Reusable_Value*vel_p(k,j,i)/C_Square

                        DRIVER_LAPSE_VAL(i-1,j-1,k-1) = Lapse_Val(k,j,i)

                    END DO
                END DO
            END DO

        END DO
    END DO
END DO







CALL INIT_CHIMERA_POTENTIAL( k_lower, k_upper, j_lower, j_upper, i_lower, i_upper, Density )

Num_Nodes = (/ DRIVER_R_INPUT_NODES, DRIVER_T_INPUT_NODES, DRIVER_P_INPUT_NODES /)
CALL INIT_CHIMERA_SHIFT_VAL( Num_Nodes,                                         &
                             DRIVER_R_ELEMS,DRIVER_T_ELEMS,DRIVER_P_ELEMS,   &
                             DRIVER_Si,                                        &
                             DRIVER_R_LOCS                                     &
                            )




CALL WRITE_CHIMERA_OUTPUT( file_number )



!  Release Temporary Variables !
DEALLOCATE( Density )
DEALLOCATE( Pressure )
DEALLOCATE( vel_r )
DEALLOCATE( vel_t )
DEALLOCATE( vel_p )
DEALLOCATE( int_e )
DEALLOCATE( Lapse_Val )



END SUBROUTINE LOAD_CHIMERA_HDF5



!+201+###########################################################################!
!                                                                                !
!               OPEN_CHIMERA_HDF5                                                !
!                                                                                !
!################################################################################!
SUBROUTINE OPEN_CHIMERA_HDF5(file_number, path, frame_flag)


!ik_ray_dim

INTEGER, INTENT(IN)                     :: file_number
LOGICAL, INTENT(IN)                     :: frame_flag

CHARACTER(len=*), INTENT(IN)            :: path



INTEGER                                 :: error
INTEGER                                 :: iproc
INTEGER                                 :: nproc_per_hyperslab_group
INTEGER                                 :: mpi_group_per_hyperslab_group
INTEGER                                 :: mpi_world_group

INTEGER, DIMENSION(1)                   :: slab_data
INTEGER, DIMENSION(:), ALLOCATABLE      :: hyperslab_group_member








IF (myID == 0 ) THEN

    ! Create File name !
    IF ( frame_flag ) THEN
        WRITE(suffix,fmt='(i5.5,a6,i1,a1,i2.2,a3)') file_number,'_grid_',my_grid,'_',1,'.h5    '
    ELSE
        WRITE(suffix, fmt='(i9.9,a6,i1,a1,i2.2,a3)') file_number,'_grid_',my_grid,'_',1,'.h5'
    END IF

    ! Turn on HDF5 !
    CALL h5open_f(error)

    ! Open File !
    CALL h5fopen_f(TRIM(path)//'/'//prefix//TRIM(suffix), H5F_ACC_RDONLY_F, file_id, error)

    ! Open 'mesh' group !
    CALL h5gopen_f(file_id, '/mesh', group_id, error)

    ! Read nz_hyperslabs from '/mesh' !
    datasize1d(1) = 1
    CALL read_1d_slab('nz_hyperslabs',slab_data, group_id, datasize1d)
    nz_hyperslabs = slab_data(1)

    ! Close group !
    CALL h5gclose_f(group_id, error)

    ! Close File !
    CALL h5fclose_f(file_id, error)

    ! Close HDF5 !
    CALL h5close_f(error)

END IF ! myID == 0

CALL MPI_BCAST( nz_hyperslabs, 1, MPI_INTEGER, 0, MPI_COMM_GRID, error )


j_offset = 0
k_offset = myID_phi * ik_ray_dim

!-----------------------------------------------------------------------
! Initialize hyperslabs and communicators
!-----------------------------------------------------------------------

nz_hyperslab_width = DRIVER_P_ELEMS / nz_hyperslabs
IF( MOD( DRIVER_z_PROCS, nz_hyperslabs ) /= 0 ) THEN
    PRINT*, 'DRIVER_z_PROCS should be evenly divisible by nz_hyperslabs'
    PRINT*, 'DRIVER_z_PROCS = ', DRIVER_z_PROCS, ' nz_hyperslabs = ', nz_hyperslabs
    CALL MPI_ABORT( MPI_COMM_WORLD, 666, error)
END IF

nz_hyperslab_procwidth = DRIVER_z_PROCS / nz_hyperslabs
IF( nz_hyperslabs > DRIVER_z_PROCS ) THEN
    PRINT*, 'nz_hyperslabs cannot be more than DRIVER_z_PROCS'
    PRINT*, 'DRIVER_z_PROCS = ', DRIVER_z_PROCS, ' nz_hyperslabs = ', nz_hyperslabs
    CALL MPI_ABORT( MPI_COMM_WORLD, 666, error)
END IF

!-- Set the k index in the hyperslab group
k_hyperslab_offset = MOD( k_offset, nz_hyperslab_width )

!-- Create MPI communicator for each group of writers
my_hyperslab_group = myID_phi / nz_hyperslab_procwidth + 1
nproc_per_hyperslab_group = DRIVER_z_PROCS / nz_hyperslabs

ALLOCATE(hyperslab_group_member(nproc_per_hyperslab_group))

hyperslab_group_member &
&       = (/(iproc, iproc = (my_hyperslab_group-1)*nproc_per_hyperslab_group, &
                                my_hyperslab_group*nproc_per_hyperslab_group-1)/)

!PRINT*,"my_hyperslab_group",my_hyperslab_group, myID, myID_phi,nz_hyperslab_procwidth
!PRINT*,"nproc_per_hyperslab_group",nproc_per_hyperslab_group
!IF(myID == 0) THEN
!   PRINT*,"hyperslab_Group_member",hyperslab_group_member
!END IF

!CALL MPI_COMM_GROUP(MPI_COMM_GRID, mpi_world_group, error)

!CALL MPI_GROUP_INCL(mpi_world_group, nproc_per_hyperslab_group, &
!&      hyperslab_group_member, mpi_group_per_hyperslab_group, error)

!CALL MPI_COMM_CREATE(MPI_COMM_GRID, mpi_group_per_hyperslab_group, &
!&      mpi_comm_per_hyperslab_group, error)

!PRINT*,"HERE error",error, myID, my_hyperslab_group, nproc_per_hyperslab_group



!CALL MPI_GROUP_FREE(mpi_group_per_hyperslab_group, error)
!DEALLOCATE(hyperslab_group_member)

!-- Set the "master" (ie. first processor) of each hyperslab group
!hyperslab_group_master = .false.
!PRINT*,"HERE",myID
!CALL MPI_COMM_RANK(mpi_comm_per_hyperslab_group, myid_hyperslab_group, error)
!PRINT*,"THERE",myID

!IF( myid_hyperslab_group == 0 ) hyperslab_group_master = .true.

!-----------------------------------------------------------------------
! Open HDF5 file for reading
!-----------------------------------------------------------------------
IF ( frame_flag ) THEN
  WRITE(suffix, fmt='(i5.5,a6,i1,a1,i2.2,a3)') file_number,'_grid_',my_grid,'_',my_hyperslab_group,'.h5    '
ELSE
  WRITE(suffix, fmt='(i9.9,a6,i1,a1,i2.2,a3)') file_number,'_grid_',my_grid,'_',my_hyperslab_group,'.h5'
END IF

CALL h5open_f(error)
CALL h5fopen_f(TRIM(path)//'/'//prefix//TRIM(suffix), H5F_ACC_RDONLY_F, file_id, error)

IF ( error /= 0 ) THEN
  PRINT*, '***ERROR in trying to open ', TRIM(path)//'/'//prefix//suffix
  WRITE(nlog, '(a50)') '***ERROR in trying to open ', TRIM(path)//'/'//prefix//suffix
  CALL MPI_ABORT(MPI_COMM_WORLD,666, error)
END IF






END SUBROUTINE OPEN_CHIMERA_HDF5






!+202+###########################################################################!
!                                                                                !
!              CLOSE_CHIMERA_HDF5                                                !
!                                                                                !
!################################################################################!
SUBROUTINE CLOSE_CHIMERA_HDF5()

INTEGER                                 :: error




! Close file !
CALL h5fclose_f(file_id, error)
IF ( error /= 0 ) THEN
    PRINT*, '***ERROR in trying to close file with id = ', file_id
    WRITE(nlog, '(a50)') '***ERROR in trying to close file with id = ', file_id
    CALL MPI_ABORT(MPI_COMM_WORLD,666, error)
END IF


! Close HDF5 !
CALL h5close_f(error)
IF ( error /= 0 ) THEN
    PRINT*, '***ERROR in trying to close HDF5 '
    WRITE(nlog, '(a50)') '***ERROR in trying to close HDF5'
    CALL MPI_ABORT(MPI_COMM_WORLD,666, error)
END IF




END SUBROUTINE CLOSE_CHIMERA_HDF5




















!+301+###########################################################################!
!                                                                                !
!              READ_CHIMERA_HDF5                                                 !
!                                                                                !
!################################################################################!
SUBROUTINE READ_CHIMERA_HDF5( file_number, stride )

INTEGER, INTENT(in)                                 :: file_number
INTEGER, INTENT(in)                                 :: stride



INTEGER                                             :: error
INTEGER                                             :: i,j,k ! space loop counter

INTEGER, DIMENSION(3)                               :: array_dimensions

INTEGER, DIMENSION(2)                               :: radial_index_bound
INTEGER, DIMENSION(2)                               :: theta_index_bound
INTEGER, DIMENSION(2)                               :: phi_index_bound



!REAL(KIND = idp), DIMENSION(0:DRIVER_R_ELEMS)                      :: x_ef
!REAL(KIND = idp), DIMENSION(0:DRIVER_T_ELEMS)                      :: y_ef
!REAL(KIND = idp), DIMENSION(0:DRIVER_P_ELEMS)                      :: z_ef
!
!REAL(KIND = idp), DIMENSION(0:DRIVER_R_ELEMS-1)                    :: dx_cf
!REAL(KIND = idp), DIMENSION(0:DRIVER_T_ELEMS-1)                    :: dy_cf
!REAL(KIND = idp), DIMENSION(0:DRIVER_P_ELEMS-1)                    :: dz_cf

REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                        :: x_ef
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                        :: y_ef
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                        :: z_ef

REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                        :: dx_cf
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                        :: dy_cf
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                        :: dz_cf



REAL(KIND = idp), DIMENSION(DRIVER_T_ELEMS,DRIVER_P_ELEMS)    :: ongrid_mask

REAL(KIND = idp)                                                :: r2dr     ! 1/3 (r_outer^3 - r_inner^3)
REAL(KIND = idp), DIMENSION(DRIVER_R_ELEMS)                    :: send_buff, recv_buff

REAL(KIND = idp)                                                :: lw_max, l_x, l_y, l_z, rsinth  ! widths of zones
REAL(KIND = idp)                                                :: min_omega



! Open 'mesh' group !
CALL h5gopen_f(file_id, '/mesh', group_id, error)


! Read and Check Array Dimensions !
! Spatial Elements !
datasize1d(1) = 3
CALL read_1d_slab('array_dimensions', array_dimensions, group_id, datasize1d)

!IF ( array_dimensions(1) /= DRIVER_R_ELEMS ) THEN
!    PRINT *, '*** ERROR: non-matching value of nx:', DRIVER_R_ELEMS, '/=', array_dimensions(1)
!    CALL MPI_ABORT(MPI_COMM_WORLD,101,error)
!END IF
!
!IF ( array_dimensions(2) /= DRIVER_T_ELEMS ) THEN
!    PRINT *, '*** ERROR: non-matching value of ny:', DRIVER_T_ELEMS, '/=', array_dimensions(2)
!    CALL MPI_ABORT(MPI_COMM_WORLD,102,error)
!END IF
!
!IF ( array_dimensions(3) /= DRIVER_P_ELEMS ) THEN
!    PRINT *, '*** ERROR: non-matching value of nz:', DRIVER_P_ELEMS, '/=', array_dimensions(3)
!    CALL MPI_ABORT(MPI_COMM_WORLD,103,error)
!END IF


IF ( array_dimensions(1) < DRIVER_R_ELEMS ) THEN
    PRINT *, '*** ERROR: Incompatable value of nx, too large:', DRIVER_R_ELEMS, '>', array_dimensions(1)
    CALL MPI_ABORT(MPI_COMM_WORLD,101,error)
END IF

IF ( array_dimensions(2) < DRIVER_T_ELEMS ) THEN
    PRINT *, '*** ERROR: Incompatable value of ny, too large:', DRIVER_T_ELEMS, '>', array_dimensions(2)
    CALL MPI_ABORT(MPI_COMM_WORLD,102,error)
END IF

IF ( array_dimensions(3) < DRIVER_P_ELEMS ) THEN
    PRINT *, '*** ERROR: Incompatable value of nz, too large:', DRIVER_P_ELEMS, '>', array_dimensions(3)
    CALL MPI_ABORT(MPI_COMM_WORLD,103,error)
END IF


ALLOCATE( x_ef(0:array_dimensions(1)), dx_cf(0:array_dimensions(1)-1) )
ALLOCATE( y_ef(0:array_dimensions(2)), dy_cf(0:array_dimensions(2)-1) )
ALLOCATE( z_ef(0:array_dimensions(3)), dz_cf(0:array_dimensions(3)-1) )


datasize1d(1) = 1
CALL HDF5_READ('time',time, group_id, datasize1d, error )


datasize1d(1) = array_dimensions(1)
CALL read_1d_slab('x_ef', x_ef, group_id, datasize1d)
CALL read_1d_slab('dx_cf',dx_cf, group_id, datasize1d)

datasize1d(1) = array_dimensions(2)
CALL read_1d_slab('y_ef', y_ef, group_id, datasize1d)
CALL read_1d_slab('dy_cf',dy_cf, group_id, datasize1d)

datasize1d(1) = array_dimensions(2)
CALL read_1d_slab('z_ef', z_ef, group_id, datasize1d)
CALL read_1d_slab('dz_cf',dz_cf, group_id, datasize1d)


DRIVER_R_LOCS = x_ef(0:DRIVER_R_ELEMS)
DRIVER_T_LOCS = y_ef(0:DRIVER_T_ELEMS)
DRIVER_P_LOCS = z_ef(0:DRIVER_P_ELEMS)


DRIVER_Delta_R = dx_cf(0:DRIVER_R_ELEMS-1)
DRIVER_Delta_T = dy_cf(0:DRIVER_T_ELEMS-1)
DRIVER_Delta_P = dz_cf(0:DRIVER_P_ELEMS-1)



! Close Mesh Group !
CALL h5gclose_f(group_id, error)








! Read Density Information !

! Open fluid Group !
CALL h5gopen_f(file_id, '/fluid', group_id, error)

datasize3d = (/DRIVER_R_ELEMS, ij_ray_dim, ik_ray_dim/)
slab_offset3d = (/0,j_offset,k_hyperslab_offset/)


! Read Density !
CALL read_scalar_HDF(Density, 'rho_c', group_id, datasize3d, slab_offset3d )

! Read Pressure !
CALL read_scalar_HDF(Pressure,'press', group_id, datasize3d, slab_offset3d )

! Read Radial Velocity !
CALL read_scalar_HDF(vel_r,'u_c', group_id, datasize3d, slab_offset3d )

! Read Theta Velocity !
CALL read_scalar_HDF(vel_t,'v_c', group_id, datasize3d, slab_offset3d )


! Read Phi Velocity !
CALL read_scalar_HDF(vel_p,'w_c', group_id, datasize3d, slab_offset3d )

! Read Internal Energy  !
CALL read_scalar_HDF(int_e,'e_int', group_id, datasize3d, slab_offset3d )

! Close fluid group
CALL h5gclose_f(group_id, error)




!! Open fluid Group !
CALL h5gopen_f(file_id, '/metadata', group_id, error)

datasize3d = (/DRIVER_R_ELEMS, ij_ray_dim, ik_ray_dim/)
slab_offset3d = (/0,j_offset,k_hyperslab_offset/)


! Read Lapse Function !
CALL read_scalar_HDF(Lapse_Val,'agr_e',group_id, datasize3d, slab_offset3d )


! Close fluid group
CALL h5gclose_f(group_id, error)







END SUBROUTINE READ_CHIMERA_HDF5










!+401+###########################################################################!
!                                                                                !
!              READ_SCALAR_HDF                                                   !
!                                                                                !
!################################################################################!
SUBROUTINE READ_SCALAR_HDF( var_out, var_name, group_id, datasize, slab_offset )


CHARACTER(*), INTENT(IN)                            :: var_name

INTEGER(HID_T)                                      :: group_id
INTEGER(HSIZE_T), dimension(3), INTENT(IN)          :: datasize
INTEGER(HSIZE_T), dimension(3), INTENT(IN)          :: slab_offset

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)     :: var_out



INTEGER                                                             :: ik, ij, ix
REAL(KIND = idp), DIMENSION(datasize(1),datasize(2),datasize(3))    :: var_in


CALL read_ray_hyperslab( var_name, var_in, group_id, datasize, slab_offset )


DO ij = 1,datasize(2)
    DO ik = 1,datasize(3)
        DO ix = imin, imax
            var_out(ik,ij,ix) = var_in(ix,ij,ik)
        END DO
    END DO
END DO


END SUBROUTINE READ_SCALAR_HDF







!+501+###########################################################################!
!                                                                                !
!              WRITE_CHIMERA_OUTPUT                                              !
!                                                                                !
!################################################################################!
SUBROUTINE WRITE_CHIMERA_OUTPUT( file_number )

INTEGER, INTENT(IN)                                 ::  file_number

INTEGER                                             ::  File_IDa = 72
INTEGER                                             ::  File_IDb = 73
INTEGER                                             ::  File_IDc = 74
INTEGER                                             ::  File_IDd = 75
INTEGER                                             ::  File_IDe = 76
INTEGER                                             ::  File_IDf = 77
INTEGER                                             ::  File_IDg = 78
INTEGER                                             ::  File_IDh = 79

CHARACTER(LEN = :), ALLOCATABLE                     ::  file_namea
CHARACTER(LEN = :), ALLOCATABLE                     ::  file_nameb
CHARACTER(LEN = :), ALLOCATABLE                     ::  file_namec
CHARACTER(LEN = :), ALLOCATABLE                     ::  file_named
CHARACTER(LEN = :), ALLOCATABLE                     ::  file_namee
CHARACTER(LEN = :), ALLOCATABLE                     ::  file_namef
CHARACTER(LEN = :), ALLOCATABLE                     ::  file_nameg
CHARACTER(LEN = :), ALLOCATABLE                     ::  file_nameh

LOGICAL                                             ::  exist
INTEGER                                             ::  i, j, k

110 FORMAT (11X,A1,24X,A5)
111 FORMAT (11X,A1,24X,A1,22X,A2,22X,A1,22X,A3)
112 FORMAT (11X,A1,23X,A3,17X,A8)

114 FORMAT (E22.15,3x,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15)
115 FORMAT (E22.15,3X,ES22.15,3X,ES22.15)


IF ( SOURCE_OUTPUT_FLAG == 1 ) THEN

    ALLOCATE( character(len=50) :: file_namea )
    ALLOCATE( character(len=48) :: file_nameb )
    ALLOCATE( character(len=49) :: file_namec )
    ALLOCATE( character(len=46) :: file_named )
    ALLOCATE( character(len=46) :: file_namee )
    ALLOCATE( character(len=46) :: file_namef )
    ALLOCATE( character(len=50) :: file_nameg )
    ALLOCATE( character(len=50) :: file_nameh )

    WRITE(file_namea,'(A37,I5.5,A4)')'OUTPUT/CHIMERA_RESULTS/CHIMERA_Lapse_',file_number,'.out'
    WRITE(file_nameb,'(A39,I5.5,A4)')'OUTPUT/CHIMERA_RESULTS/CHIMERA_Density_',file_number,'.out'
    WRITE(file_namec,'(A40,I5.5,A4)')'OUTPUT/CHIMERA_RESULTS/CHIMERA_Velocity_',file_number,'.out'
    WRITE(file_named,'(A37,I5.5,A4)')'OUTPUT/CHIMERA_RESULTS/CHIMERA_rlocs_',file_number,'.out'
    WRITE(file_namee,'(A37,I5.5,A4)')'OUTPUT/CHIMERA_RESULTS/CHIMERA_tlocs_',file_number,'.out'
    WRITE(file_namef,'(A37,I5.5,A4)')'OUTPUT/CHIMERA_RESULTS/CHIMERA_plocs_',file_number,'.out'
    WRITE(file_nameg,'(A41,I5.5,A4)')'OUTPUT/CHIMERA_RESULTS/CHIMERA_potential_',file_number,'.out'


    OPEN(unit = file_ida,file = file_namea)
    OPEN(unit = file_idb,file = file_nameb)
    OPEN(unit = file_idc,file = file_namec)
    OPEN(unit = file_idd,file = file_named)
    OPEN(unit = file_ide,file = file_namee)
    OPEN(unit = file_idf,file = file_namef)
    OPEN(unit = file_idg,file = file_nameg)


    ! Output CHIMERA rlocs
    WRITE(file_IDd,*)1.0_idp
    DO i = 0,DRIVER_R_ELEMS-1
        WRITE(file_IDd,*) (DRIVER_R_LOCS(i+1)+DRIVER_R_LOCS(i))/2.0_idp
    END DO

    ! Output CHIMERA tlocs
    DO i = 0,ij_ray_dim-1
        WRITE(file_IDe,*) (DRIVER_T_LOCS(i+1)+DRIVER_T_LOCS(i))/2.0_idp
    END DO

    ! Output CHIMERA tlocs
    DO i = 0,ik_ray_dim-1
        WRITE(file_IDf,*) (DRIVER_P_LOCS(i+1)+DRIVER_P_LOCS(i))/2.0_idp
    END DO


    ! Output CHIMERA Density
    DO k = 1,ik_ray_dim
        DO j = 1, ij_ray_dim

            WRITE(file_idb,*)Density(k,j,1),Density(k,j,:)

        END DO ! j Loop
    END DO ! k Loop

    ! Output CHIMERA Velocity
    DO k = 1,ik_ray_dim
        DO j = 1, ij_ray_dim

            WRITE(file_idc,*)Vel_R(k,j,1),Vel_R(k,j,:)

        END DO ! j Loop
    END DO ! k Loop


    ! Output CHIMERA Potential
    WRITE(file_IDg,*)DRIVER_Potential


    WRITE(File_IDa,*)DRIVER_Lapse_val



    CLOSE(unit = file_ida)
    CLOSE(unit = file_idb)
    CLOSE(unit = file_idc)
    CLOSE(unit = file_idd)
    CLOSE(unit = file_ide)
    CLOSE(unit = file_idf)
    CLOSE(unit = file_IDg)

END IF






!WRITE(file_nameh,'(A)') 'OUTPUT/Poseidon_Results/Frame_Time.out'
!inquire(file = file_nameh, exist=exist)
!
!IF (( exist ) .AND. (file_number .NE. 1 )) THEN
!    OPEN(unit = file_idh,file = file_nameh,status="old",position="append", action="write" )
!    WRITE(file_idh,'(I5.5,E22.15)') file_number, time
!    CLOSE(unit=file_idh)
!ELSE
!    OPEN(unit = file_idh,file = file_nameh)
!    WRITE(file_idh,'(I5.5,E22.15)') file_number, time
!    CLOSE(unit=file_idh)
!END IF

END SUBROUTINE WRITE_CHIMERA_OUTPUT












!+503+##################################################################!
!                                                                       !
!   Initialize_LG_Quadrature - Calculate the Legendre-Gauss Quadrature  !
!                              node locations and weights.              !
!                                                                       !
!#######################################################################!
SUBROUTINE Initialize_LG_Quadrature(Ord, xloc, weights)

INTEGER, INTENT(IN)                                     ::  Ord
REAL(KIND = idp), INTENT(INOUT), DIMENSION(1:Ord)       ::  xloc, weights

INTEGER                                                 :: i, j, m
REAL(KIND = idp)                                        :: p1, p2, p3, pp, z, z1


m = (Ord + 1)/2

DO i = 1,m

    z = cos(pi * (i-0.25_idp)/(Ord + 0.5_idp))
    z1 = 42.0_idp

    DO WHILE ( ABS(z - z1) .GT. eps)
        p1 = 1.0_idp
        p2 = 0.0_idp

        DO j = 1, Ord

            p3 = p2
            p2 = p1
            p1 = ((2.0_idp*j - 1.0_idp)*z*p2 - (j - 1.0_idp)*p3)/j

        END DO


        pp = Ord*(z*p1 - p2)/(z*z-1.0_idp)
        z1 = z
        z = z1 - p1/pp


    END DO

    xloc(i) = -z
    xloc(Ord-i+1) = +z

    weights(i) = 2.0_idp/((1.0_idp - z*z)*pp*pp)
    weights(Ord-i+1) = weights(i)

END DO

!PRINT*,"!!!!!!!!!!!!!!!"
!PRINT*,xloc
!PRINT*," "
!PRINT*,weights
!PRINT*,"!!!!!!!!!!!!!!!"


END SUBROUTINE Initialize_LG_Quadrature


END MODULE CHIMERA_HDF5_Module
