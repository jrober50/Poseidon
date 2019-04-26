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
                        ONLY :  Grav_Constant_G,    &
                                C_Square,           &
                                GR_Source_Scalar

USE CHIMERA_Parameters, &
                        ONLY :  CHIMERA_R_ELEMS,    &
                                CHIMERA_T_ELEMS,    &
                                CHIMERA_P_ELEMS,    &
                                CHIMERA_R_LOCS,     &
                                CHIMERA_T_LOCS,     &
                                CHIMERA_P_LOCS,     &
                                CHIMERA_Delta_R,    &
                                CHIMERA_Delta_T,    &
                                CHIMERA_Delta_P,    &
                                Enclosed_Mass,      &
                                CHIMERA_Potential,  &
                                CHIMERA_E,          &
                                CHIMERA_S,          &
                                CHIMERA_Si,         &
                                CHIMERA_T_LOCS_LOCAL,&
                                CHIMERA_P_LOCS_LOCAL,&
                                myID,               &
                                myID_theta,         &
                                myID_phi,           &
                                ij_ray_dim,         &
                                ik_ray_dim,         &
                                MPI_COMM_GRID,      &
                                CHIMERA_PROCS,      &
                                CHIMERA_y_PROCS,    &
                                CHIMERA_z_PROCS

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



CONTAINS
!+201+###########################################################################!
!                                                                                !
!               LOAD_CHIMERA_HDF5                                                !
!                                                                                !
!################################################################################!
SUBROUTINE LOAD_CHIMERA_HDF5(file_number, stride, path, frame_flag )

INTEGER, INTENT(IN)           :: file_number
INTEGER, INTENT(IN)           :: stride

LOGICAL, INTENT(IN)           :: frame_flag

CHARACTER(len=*), INTENT(IN)  :: path

INTEGER, DIMENSION(:), ALLOCATABLE      :: here
INTEGER                       :: i,j,k

REAL(KIND = idp), DIMENSION(:,:,:), ALLOCATABLE  :: Specific_Enthalpy
REAL(KIND = idp), DIMENSION(:,:,:), ALLOCATABLE  :: V_Square
REAL(KIND = idp), DIMENSION(:,:,:), ALLOCATABLE  :: Reusable_Value
REAL(KIND = idp), DIMENSION(:,:,:), ALLOCATABLE  :: Lorentz_Factor_Sqr


REAL(KIND = idp), DIMENSION(:,:,:), ALLOCATABLE  :: E
REAL(KIND = idp), DIMENSION(:,:,:,:), ALLOCATABLE  :: Si
REAL(KIND = idp), DIMENSION(:,:,:), ALLOCATABLE  :: S

INTEGER :: error

imin = 1
imax = CHIMERA_R_ELEMS

jmin = 1
jmax = CHIMERA_T_ELEMS

kmin = 1
kmax = CHIMERA_P_ELEMS

i_lower = imin
i_upper = imax

j_lower = jmin
j_upper = jmax

k_lower = kmin
k_upper = kmax

IF ( CHIMERA_T_ELEMS > 1) j_lower = 0
IF ( CHIMERA_T_ELEMS > 1) j_upper = ij_ray_dim + 1
IF ( CHIMERA_P_ELEMS > 1) k_lower = 0
IF ( CHIMERA_P_ELEMS > 1) k_upper = ik_ray_dim + 1




!  Allocate space for temporary storage !
ALLOCATE( Density(k_lower:k_upper,j_lower:j_upper,i_lower:i_upper) )
ALLOCATE( Pressure(k_lower:k_upper,j_lower:j_upper,i_lower:i_upper) )
ALLOCATE( vel_r(k_lower:k_upper,j_lower:j_upper,i_lower:i_upper) ) 
ALLOCATE( vel_t(k_lower:k_upper,j_lower:j_upper,i_lower:i_upper) )
ALLOCATE( vel_p(k_lower:k_upper,j_lower:j_upper,i_lower:i_upper) )
ALLOCATE( int_e(k_lower:k_upper,j_lower:j_upper,i_lower:i_upper) )

ALLOCATE( Specific_Enthalpy(k_lower:k_upper,j_lower:j_upper,i_lower:i_upper) )
ALLOCATE( V_Square(k_lower:k_upper,j_lower:j_upper,i_lower:i_upper) )
ALLOCATE( Lorentz_Factor_Sqr(k_lower:k_upper,j_lower:j_upper,i_lower:i_upper) )
ALLOCATE( Reusable_Value(k_lower:k_upper,j_lower:j_upper,i_lower:i_upper) )

ALLOCATE( E(k_lower:k_upper,j_lower:j_upper,i_lower:i_upper) )
ALLOCATE( Si(1:3,k_lower:k_upper,j_lower:j_upper,i_lower:i_upper) )
ALLOCATE( S(k_lower:k_upper,j_lower:j_upper,i_lower:i_upper) )

ALLOCATE( here(k_lower:k_upper) )

PRINT*,"Before OPEN_CHIMERA_HDF5",myID
!  Open CHIMERA HDF5 File  !
CALL OPEN_CHIMERA_HDF5( file_number, path, frame_flag )


CALL MPI_BARRIER(MPI_COMM_WORLD, error)
PRINT*,"Before READ_CHIMERA_HDF5",myID
!  Read Relevant data into arrays !
CALL READ_CHIMERA_HDF5( file_number, stride )


CALL MPI_BARRIER(MPI_COMM_WORLD, error)
PRINT*,"Before CLOSE_CHIMERA_HDF5",myID
!  Close CHIMERA HDF5 File  !
CALL CLOSE_CHIMERA_HDF5()










Density(:,:,CHIMERA_R_ELEMS) = Density(:,:,CHIMERA_R_ELEMS-1)



!  Convert CHIMERA data into CFA Input !

!Enclosed_Mass = (4.0_idp/3.0_idp)*pi* ( r_{i+1)^3 - r_{i}^3 ) * rho_i
Enclosed_Mass(0) = (4.0_idp/3.0_idp)*pi*DENSITY(k_lower,j_lower,1)     &
                                       *CHIMERA_R_LOCS(0)**3    

DO i = 1,CHIMERA_R_ELEMS

    Enclosed_Mass(i) = Enclosed_Mass(i-1)                                 &
                       + (4.0_idp/3.0_idp)*pi*DENSITY(k_lower,j_lower,i)  &
                                          *( CHIMERA_R_LOCS(i)**3       &
                                           - CHIMERA_R_LOCS(i-1)**3)

END DO 








CHIMERA_Potential(CHIMERA_R_ELEMS) = -Grav_Constant_G * Enclosed_Mass(CHIMERA_R_ELEMS)       &
                                                      /CHIMERA_R_LOCS(CHIMERA_R_ELEMS)


DO i = CHIMERA_R_ELEMS-1,2, -1

     CHIMERA_Potential(i) = CHIMERA_Potential(i+1) - Grav_Constant_G                         &
                                                   * Enclosed_Mass(i)                        &
                                                   / (CHIMERA_R_LOCS(i)*CHIMERA_R_LOCS(i))   &
                                                   * CHIMERA_Delta_R(i)

END DO
!CHIMERA_Potential(1) = CHIMERA_Potential(2)/2.0_idp

CHIMERA_Potential(1) = CHIMERA_Potential(2) - 3*Grav_Constant_G*Enclosed_Mass(2)     &
                                              /(2*CHIMERA_R_LOCS(2))






int_e = 0.0_idp
Pressure = 0.0_idp

Specific_Enthalpy = C_Square + (int_e + Pressure)/Density

V_Square = vel_r*vel_r + vel_t*vel_t + vel_p*vel_p
Lorentz_Factor_Sqr = 1.0_idp/(1.0_idp - V_Square/C_Square)

Reusable_Value = (GR_Source_Scalar*Density*Specific_Enthalpy) / (1.0_idp - V_Square/C_Square)

IF (ANY(DENSITY .eq. 0.0_idp) ) THEN

    DO i = i_lower, i_upper

        IF ( DENSITY(k_lower, j_lower, i) == 0.0_idp) THEN

            Reusable_Value(:,:,i) = 0.0_idp

        END IF
    END DO

END IF


!CHIMERA_E = Reusable_Value - GR_Source_Scalar*Pressure
!CHIMERA_Si(1,k_lower:k_upper,j_lower:j_upper,i_lower:i_upper) = Reusable_Value(k_lower:k_upper,j_lower:j_upper,i_lower:i_upper)   &
!                                                              * vel_r(k_lower:k_upper,j_lower:j_upper,i_lower:i_upper)
!CHIMERA_Si(2,:,:,:) = Reusable_Value*vel_t
!CHIMERA_Si(3,:,:,:) = Reusable_Value*vel_p
!CHIMERA_S = Reusable_Value*V_Square - 3.0_idp*GR_Source_Scalar*Pressure



CHIMERA_S = 0.0_idp
CHIMERA_Si = 0.0_idp


DO k = k_lower, k_upper
   DO j = j_lower, j_upper
      DO i = i_lower, i_upper

!            CHIMERA_E(i,j,k) = GR_Source_Scalar*Density(k,j,i)*C_Square

           CHIMERA_E(i,j,k) = Reusable_Value(k,j,i) - GR_SOURCE_Scalar*Pressure(k,j,i)
!           CHIMERA_S(i,j,k) = Reusable_Value(k,j,i)*V_Square(k,j,i) - 3.0_idp * GR_SOURCE_Scalar*Pressure(k,j,i)
!           CHIMERA_Si(i,j,k,1) = Reusable_Value(k,j,i)*vel_r(k,j,i)
!           CHIMERA_Si(i,j,k,2) = Reusable_Value(k,j,i)*vel_t(k,j,i)
!           CHIMERA_Si(i,j,k,2) = Reusable_Value(k,j,i)*vel_p(k,j,i)

      END DO
   END DO
END DO





IF (myID == -130) THEN

   PRINT*," "
   PRINT*,"Density"
   PRINT*,"i_upper",i_upper
   PRINT*," "
   DO i = i_lower, i_upper
      DO j = j_lower, j_upper
!          PRINT*,i,CHIMERA_R_LOCS(i),CHIMERA_E(i,j_lower,k_lower)
          PRINT*,i,CHIMERA_R_LOCS(i),Density(k_lower,j_lower,i), int_e(k_lower, j_lower, i), Pressure(k_lower, j_lower, i)
!          PRINT*,i,CHIMERA_R_LOCS(i), Density(k_lower, j_lower, i), CHIMERA_E(i,j_lower, k_lower)
!          PRINT*,i,CHIMERA_R_LOCS(i), Enclosed_Mass(i), CHIMERA_Potential(i)
!          PRINT*,i,CHIMERA_E(i,j_lower, k_lower),Reusable_Value(k_lower,j_lower,i), Density(k_lower, j_lower,i)
      END DO
   END DO
   PRINT*," "
   PRINT*," "


END IF





!  Release Temporary Variables !
DEALLOCATE( Density )
DEALLOCATE( Pressure )
DEALLOCATE( vel_r )
DEALLOCATE( vel_t )
DEALLOCATE( vel_p )
DEALLOCATE( int_e )

DEALLOCATE( Specific_Enthalpy )
DEALLOCATE( V_Square )
!DEALLOCATE( Lorentz_Factor_Sqr )
DEALLOCATE( Reusable_Value )

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

    PRINT*," End of myID==0"

END IF ! myID == 0

CALL MPI_BCAST( nz_hyperslabs, 1, MPI_INTEGER, 0, MPI_COMM_GRID, error )


j_offset = 0
k_offset = myID_phi * ik_ray_dim

!-----------------------------------------------------------------------
! Initialize hyperslabs and communicators
!-----------------------------------------------------------------------

nz_hyperslab_width = CHIMERA_P_ELEMS / nz_hyperslabs
IF( MOD( CHIMERA_z_PROCS, nz_hyperslabs ) /= 0 ) THEN
    PRINT*, 'CHIMERA_z_PROCS should be evenly divisible by nz_hyperslabs'
    PRINT*, 'CHIMERA_z_PROCS = ', CHIMERA_z_PROCS, ' nz_hyperslabs = ', nz_hyperslabs
    CALL MPI_ABORT( MPI_COMM_WORLD, 666, error)
END IF

nz_hyperslab_procwidth = CHIMERA_z_PROCS / nz_hyperslabs
IF( nz_hyperslabs > CHIMERA_z_PROCS ) THEN
    PRINT*, 'nz_hyperslabs cannot be more than CHIMERA_z_PROCS'
    PRINT*, 'CHIMERA_z_PROCS = ', CHIMERA_z_PROCS, ' nz_hyperslabs = ', nz_hyperslabs
    CALL MPI_ABORT( MPI_COMM_WORLD, 666, error)
END IF

!-- Set the k index in the hyperslab group
k_hyperslab_offset = MOD( k_offset, nz_hyperslab_width )

!-- Create MPI communicator for each group of writers
my_hyperslab_group = myID_phi / nz_hyperslab_procwidth + 1
nproc_per_hyperslab_group = CHIMERA_z_PROCS / nz_hyperslabs

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
INTEGER                                 :: iproc
INTEGER                                 :: nproc_per_hyperslab_group
INTEGER                                 :: mpi_world_group
INTEGER                                 :: mpi_group_per_hyperslab_group

INTEGER, DIMENSION(1)                   :: slab_data

INTEGER, DIMENSION(:), ALLOCATABLE      :: hyperslab_group_member



! Close file !
CALL h5fclose_f(file_id, error)

! Close HDF5 !
CALL h5close_f(error)


IF ( error /= 0 ) THEN
    PRINT*, '***ERROR in trying to close file '
    WRITE(nlog, '(a50)') '***ERROR in trying to close file'
    CALL MPI_ABORT(MPI_COMM_WORLD,666, error)
END IF

! Destroy hyperslab comminicator so it may be recreated for a new frame

!CALL MPI_COMM_FREE(mpi_comm_per_hyperslab_group, error)

IF ( error /= 0 ) THEN
    PRINT*, '***ERROR in trying to close file '
    WRITE(nlog, '(a50)') '***ERROR in trying to close file'
    CALL MPI_ABORT(MPI_COMM_WORLD,666, error)
END IF

RETURN
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



REAL(KIND = idp), DIMENSION(0:CHIMERA_R_ELEMS)                  :: x_ef
REAL(KIND = idp), DIMENSION(0:CHIMERA_T_ELEMS)                  :: y_ef
REAL(KIND = idp), DIMENSION(0:CHIMERA_P_ELEMS)                  :: z_ef

REAL(KIND = idp), DIMENSION(0:CHIMERA_R_ELEMS)                    :: dx_cf
REAL(KIND = idp), DIMENSION(0:CHIMERA_T_ELEMS)                    :: dy_cf
REAL(KIND = idp), DIMENSION(0:CHIMERA_P_ELEMS)                    :: dz_cf

REAL(KIND = idp), DIMENSION(1)                                  :: time_in

REAL(KIND = idp), DIMENSION(CHIMERA_T_ELEMS,CHIMERA_P_ELEMS)    :: ongrid_mask

REAL(KIND = idp)                                                :: r2dr     ! 1/3 (r_outer^3 - r_inner^3)
REAL(KIND = idp), DIMENSION(CHIMERA_R_ELEMS)                    :: send_buff, recv_buff

REAL(KIND = idp)                                                :: lw_max, l_x, l_y, l_z, rsinth  ! widths of zones
REAL(KIND = idp)                                                :: min_omega




! Open 'mesh' group !
CALL h5gopen_f(file_id, '/mesh', group_id, error)


! Read and Check Array Dimensions !
! Spatial Elements !
datasize1d(1) = 3
CALL read_1d_slab('array_dimensions', array_dimensions, group_id, datasize1d)

IF ( array_dimensions(1) /= CHIMERA_R_ELEMS ) THEN
    PRINT *, '*** ERROR: non-matching value of nx:', CHIMERA_R_ELEMS, '/=', array_dimensions(1)
    CALL MPI_ABORT(MPI_COMM_WORLD,101,error)
END IF

IF ( array_dimensions(2) /= CHIMERA_T_ELEMS ) THEN
    PRINT *, '*** ERROR: non-matching value of ny:', CHIMERA_T_ELEMS, '/=', array_dimensions(2)
    CALL MPI_ABORT(MPI_COMM_WORLD,102,error)
END IF

IF ( array_dimensions(3) /= CHIMERA_P_ELEMS ) THEN
    PRINT *, '*** ERROR: non-matching value of nz:', CHIMERA_P_ELEMS, '/=', array_dimensions(3)
    CALL MPI_ABORT(MPI_COMM_WORLD,103,error)
END IF





! Read Element Edge Locations !
datasize1d(1) = CHIMERA_R_ELEMS+1
CALL read_1d_slab('x_ef', x_ef(1:CHIMERA_R_ELEMS+1), group_id, datasize1d)
CALL read_1d_slab('dx_cf',dx_cf(0:CHIMERA_R_ELEMS), group_id, datasize1d)

datasize1d(1) = CHIMERA_T_ELEMS+1
CALL read_1d_slab('y_ef', y_ef(1:CHIMERA_T_ELEMS+1), group_id, datasize1d)
CALL read_1d_slab('dy_cf',dy_cf(0:CHIMERA_R_ELEMS), group_id, datasize1d)

datasize1d(1) = CHIMERA_P_ELEMS+1
CALL read_1d_slab('z_ef', z_ef(1:CHIMERA_P_ELEMS+1), group_id, datasize1d)
CALL read_1d_slab('dz_cf',dz_cf(0:CHIMERA_R_ELEMS), group_id, datasize1d)


CHIMERA_R_LOCS = x_ef
CHIMERA_T_LOCS = y_ef
CHIMERA_P_LOCS = z_ef

CHIMERA_Delta_R = dx_cf

CHIMERA_Delta_T = dy_cf
CHIMERA_Delta_P = dz_cf

!CHIMERA_T_LOCS_LOCAL = y_ef(j_offset+1:j_offset+ij_ray_dim+1)
!CHIMERA_P_LOCS_LOCAL = z_ef(k_offset+1:k_offset+ik_ray_dim+1)




! Close Mesh Group !
CALL h5gclose_f(group_id, error)






! Read Density Information !

! Open fluid Group !
CALL h5gopen_f(file_id, '/fluid', group_id, error)




datasize3d = (/CHIMERA_R_ELEMS, ij_ray_dim, ik_ray_dim/)
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



! Read Lapse Function !
!   datasize1d(1) = imax-imin+1
!   CALL read_1d_slab('agr_c', lapse(imin:imax), group_id, datasize1d )


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
