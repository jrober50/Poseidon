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
USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Poseidon_Numbers_Module, &
            ONLY :  pi, eps

USE Poseidon_Units_Module, &
            ONLY :  C_Square

USE Quadrature_Mapping_Functions, &
            ONLY : Quad_Map

USE Variables_IO, &
            ONLY :  Write_Flags

USE CHIMERA_TEST_FUNCS_Module,  &
            ONLY :  INIT_CHIMERA_POTENTIAL, &
                    INIT_CHIMERA_SHIFT_VAL

USE HDF5_IO_Module
USE MPI
USE HDF5

IMPLICIT NONE

INTEGER                                 :: my_grid = 1
INTEGER                                 :: nlog=6

INTEGER, DIMENSION(1:3)                             ::  Num_Nodes
INTEGER, DIMENSION(1:3)                             ::  Num_Elems

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
REAL(KIND = idp), DIMENSION(1)                      ::  time

REAL(KIND = idp), DIMENSION(:), ALLOCATABLE         :: x_ef
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE         :: y_ef
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE         :: z_ef

REAL(KIND = idp), DIMENSION(:), ALLOCATABLE         :: dx_cf
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE         :: dy_cf
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE         :: dz_cf

INTEGER                                 :: ij_ray_dim
INTEGER                                 :: ik_ray_dim
INTEGER                                 :: myid
INTEGER                                 :: myid_phi
INTEGER                                 :: zProcs



CONTAINS
!+201+###########################################################################!
!                                                                                !
!               LOAD_CHIMERA_HDF5                                                !
!                                                                                !
!################################################################################!
SUBROUTINE LOAD_CHIMERA_HDF5(file_number, stride, path, frame_flag, &
                             Num_Quad,                          &
                             R_Quad, T_Quad, P_Quad,            &
                             Left_Limit, Right_Limit,           &
                             Output_Elem,                       &
                             R_Locs, T_Locs, P_Locs,            &
                             Output_E, Output_S, Output_Si      )

INTEGER, INTENT(IN)                                 :: file_number
INTEGER, INTENT(IN)                                 :: stride

LOGICAL, INTENT(IN)                                 :: frame_flag

CHARACTER(len=*), INTENT(IN)                        :: path

INTEGER,    INTENT(IN), DIMENSION( 1:3 )                                    :: Num_Quad
REAL(idp),  INTENT(IN), DIMENSION( 1:Num_Quad(1) )                          :: R_Quad
REAL(idp),  INTENT(IN), DIMENSION( 1:Num_Quad(2) )                          :: T_Quad
REAL(idp),  INTENT(IN), DIMENSION( 1:Num_Quad(3) )                          :: P_Quad
REAL(idp),  INTENT(IN)                                                      :: Left_Limit
REAL(idp),  INTENT(IN)                                                      :: Right_Limit

INTEGER,    INTENT(OUT), DIMENSION( 1:3 )                                    :: Output_Elem
REAL(idp),  INTENT(OUT), DIMENSION( 1:Num_Elems(1)+1 )                        :: R_Locs
REAL(idp),  INTENT(OUT), DIMENSION( 1:Num_Elems(2)+1 )                        :: T_Locs
REAL(idp),  INTENT(OUT), DIMENSION( 1:Num_Elems(3)+1 )                        :: P_Locs

REAL(idp),  INTENT(OUT),DIMENSION( 1:Num_Quad(1)*Num_Quad(2)*Num_Quad(3),   &
                                   0:Num_Elems(1)-1,                         &
                                   0:Num_Elems(2)-1,                         &
                                   0:Num_Elems(3)-1  )                       ::  Output_E,   &
                                                                                Output_S

REAL(idp), INTENT(OUT), DIMENSION( 1:Num_Quad(1)*Num_Quad(2)*Num_Quad(3),   &
                                   0:Num_Elems(1)-1,                         &
                                   0:Num_Elems(2)-1,                         &
                                   0:Num_Elems(3)-1,                         &
                                   1:3                  )                   ::  Output_Si


INTEGER                                             :: here
INTEGER                                             :: i,j,k
INTEGER                                             :: pd, td, rd

REAL(KIND = idp)                                    :: Specific_Enthalpy
REAL(KIND = idp)                                    :: V_Square
REAL(KIND = idp)                                    :: Reusable_Value



REAL(KIND = idp), DIMENSION(:,:,:,:), ALLOCATABLE   ::  HR_Density

REAL(KIND = idp), DIMENSION(:), ALLOCATABLE         ::  Quad_Locs
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE         ::  Quad_Weights

REAL(KIND = idp)                                    ::  R_Center
REAL(KIND = idp)                                    ::  deltar_overtwo
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE         ::  CUR_R_LOCS

INTEGER :: error










!  Open CHIMERA HDF5 File  !
CALL OPEN_CHIMERA_HDF5( file_number, path, frame_flag )


!  Read Relevant data into arrays !
CALL READ_CHIMERA_Base( Num_Elems, file_number, stride )




imin = 1
imax = Num_Elems(1)

jmin = 1
jmax = Num_Elems(2)

kmin = 1
kmax = Num_Elems(3)

i_lower = imin
i_upper = imax

j_lower = jmin
j_upper = jmax

k_lower = kmin
k_upper = kmax

ij_ray_dim = Num_Elems(2)
ik_ray_dim = Num_Elems(3)

IF ( Num_Elems(2) > 1) j_lower = 1
IF ( Num_Elems(2) > 1) j_upper = ij_ray_dim
IF ( Num_Elems(3) > 1) k_lower = 1
IF ( Num_Elems(3) > 1) k_upper = ik_ray_dim


ALLOCATE( HR_DENSITY(Num_Quad(1),k_lower:k_upper,j_lower:j_upper,i_lower:i_upper))
ALLOCATE( Quad_Locs(1:Num_Quad(1)) )
ALLOCATE( Quad_Weights(1:Num_Quad(1)) )
ALLOCATE( Cur_R_Locs(1:Num_Quad(1)) )



!  Allocate space for temporary storage !
ALLOCATE( Density(k_lower:k_upper,j_lower:j_upper,i_lower:i_upper) )
ALLOCATE( Pressure(k_lower:k_upper,j_lower:j_upper,i_lower:i_upper) )
ALLOCATE( vel_r(k_lower:k_upper,j_lower:j_upper,i_lower:i_upper) )
ALLOCATE( vel_t(k_lower:k_upper,j_lower:j_upper,i_lower:i_upper) )
ALLOCATE( vel_p(k_lower:k_upper,j_lower:j_upper,i_lower:i_upper) )
ALLOCATE( int_e(k_lower:k_upper,j_lower:j_upper,i_lower:i_upper) )
ALLOCATE( Lapse_Val(k_lower:k_upper,j_lower:j_upper,i_lower:i_upper) )


ALLOCATE( x_ef(0:Num_Elems(1)), dx_cf(0:Num_Elems(1)-1) )
ALLOCATE( y_ef(0:Num_Elems(2)), dy_cf(0:Num_Elems(2)-1) )
ALLOCATE( z_ef(0:Num_Elems(3)), dz_cf(0:Num_Elems(3)-1) )


CALL Initialize_LG_Quadrature(Num_Quad(1), Quad_Locs, Quad_Weights)










!  Read Relevant data into arrays !
CALL READ_CHIMERA_HDF5( file_number, stride )

!  Close CHIMERA HDF5 File  !
CALL CLOSE_CHIMERA_HDF5()








Density(:,:,Num_Elems(1)) = Density(:,:,Num_Elems(1)-1)




!  Convert CHIMERA data into CFA Input !
DO i = i_lower, i_upper
DO j = j_lower, j_upper
DO k = k_lower, k_upper

    Specific_Enthalpy = C_Square + int_e(k,j,i) + Pressure(k,j,i)/Density(k,j,i)
    V_Square = vel_r(k,j,i)*vel_r(k,j,i) + vel_t(k,j,i)*vel_t(k,j,i) + vel_p(k,j,i)*vel_p(k,j,i)

    Reusable_Value = (Density(k,j,i)*Specific_Enthalpy)/(1.0_idp - V_Square/C_Square)

    DO rd = 1, Num_Quad(1)
    DO td = 1, Num_Quad(2)
    DO pd = 1, Num_Quad(3)

        Here = Quad_Map(rd,td,pd,Num_Quad)

        Output_E(here,i-1,j-1,k-1) = Reusable_Value - Pressure(k,j,i)
        Output_S(here,i-1,j-1,k-1) = Reusable_Value*V_Square/C_Square + 3.0_idp * Pressure(k,j,i)
        Output_Si(here,i-1,j-1,k-1,1) = Reusable_Value*vel_r(k,j,i)/C_Square
        Output_Si(here,i-1,j-1,k-1,2) = Reusable_Value*vel_t(k,j,i)/C_Square
        Output_Si(here,i-1,j-1,k-1,3) = Reusable_Value*vel_p(k,j,i)/C_Square

    END DO
    END DO
    END DO

END DO
END DO
END DO







CALL INIT_CHIMERA_POTENTIAL( k_lower, k_upper, j_lower, j_upper, i_lower, i_upper, Density )

Num_Nodes = (/ Num_Quad(1), Num_Quad(2), Num_Quad(3) /)




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

!CALL MPI_BCAST( nz_hyperslabs, 1, MPI_INTEGER, 0, MPI_COMM_GRID, error )


j_offset = 0
k_offset = myID_phi * ik_ray_dim

!-----------------------------------------------------------------------
! Initialize hyperslabs and communicators
!-----------------------------------------------------------------------

nz_hyperslab_width = Num_Elems(3)/ nz_hyperslabs
IF( MOD( zProcs, nz_hyperslabs ) /= 0 ) THEN
    PRINT*, 'zProcs should be evenly divisible by nz_hyperslabs'
    PRINT*, 'zProcs = ', zProcs, ' nz_hyperslabs = ', nz_hyperslabs
    CALL MPI_ABORT( MPI_COMM_WORLD, 666, error)
END IF

nz_hyperslab_procwidth = zProcs / nz_hyperslabs
IF( nz_hyperslabs > zProcs ) THEN
    PRINT*, 'nz_hyperslabs cannot be more than zProcs'
    PRINT*, 'zProcs = ', zProcs, ' nz_hyperslabs = ', nz_hyperslabs
    CALL MPI_ABORT( MPI_COMM_WORLD, 666, error)
END IF

!-- Set the k index in the hyperslab group
k_hyperslab_offset = MOD( k_offset, nz_hyperslab_width )

!-- Create MPI communicator for each group of writers
my_hyperslab_group = myID_phi / nz_hyperslab_procwidth + 1
nproc_per_hyperslab_group = zProcs / nz_hyperslabs

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
!              READ_CHIMERA_Base                                                 !
!                                                                                !
!################################################################################!
SUBROUTINE READ_CHIMERA_Base( Num_Elems, file_number, stride )

INTEGER, INTENT(OUT), DIMENSION(1:3)                :: Num_Elems
INTEGER, INTENT(in)                                 :: file_number
INTEGER, INTENT(in)                                 :: stride

INTEGER                                             :: error

! Open 'mesh' group !
CALL h5gopen_f(file_id, '/mesh', group_id, error)


! Read and Check Array Dimensions !
! Spatial Elements !
datasize1d(1) = 3
CALL read_1d_slab('array_dimensions', Num_Elems, group_id, datasize1d)
! Close Mesh Group !
CALL h5gclose_f(group_id, error)




END SUBROUTINE READ_CHIMERA_Base










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


REAL(KIND = idp), DIMENSION(Num_Elems(2),Num_Elems(3))     :: ongrid_mask

REAL(KIND = idp)                                            :: r2dr     ! 1/3 (r_outer^3 - r_inner^3)
REAL(KIND = idp), DIMENSION(Num_Elems(1))                   :: send_buff, recv_buff

REAL(KIND = idp)                                            :: lw_max, l_x, l_y, l_z, rsinth  ! widths of zones
REAL(KIND = idp)                                            :: min_omega



! Open 'mesh' group !
CALL h5gopen_f(file_id, '/mesh', group_id, error)



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




! Close Mesh Group !
CALL h5gclose_f(group_id, error)








! Read Density Information !

! Open fluid Group !
CALL h5gopen_f(file_id, '/fluid', group_id, error)

datasize3d = (/num_elems(1), ij_ray_dim, ik_ray_dim/)
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

datasize3d = (/Num_Elems(1), ij_ray_dim, ik_ray_dim/)
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
