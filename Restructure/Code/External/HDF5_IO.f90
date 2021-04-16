   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE HDF5_IO_Module                                                               !##!
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
            ONLY : idp

USE HDF5

IMPLICIT NONE


PUBLIC HDF5_OPENFILE
PUBLIC HDF5_CLOSEFILE
PUBLIC HDF5_OPENGROUP
PUBLIC HDF5_CLOSEGROUP
PUBLIC HDF5_READ

PRIVATE HDF5_READ_DOUBLE_1D
PRIVATE HDF5_READ_DOUBLE_2D
PRIVATE HDF5_READ_DOUBLE_3D
PRIVATE HDF5_READ_INTEGER_1D
PRIVATE HDF5_READ_INTEGER_2D
PRIVATE HDF5_READ_INTEGER_3D

PUBLIC  READ_1D_SLAB
PRIVATE READ_1D_SLAB_INT
PRIVATE READ_1D_SLAB_DOUBLE

PUBLIC  READ_RAY_HYPERSLAB
PRIVATE READ_RAY_HYPERSLAB_DBL_2D
PRIVATE READ_RAY_HYPERSLAB_DBL_3D
PRIVATE READ_RAY_HYPERSLAB_DBL_4D
PRIVATE READ_RAY_HYPERSLAB_DBL_5D
PRIVATE READ_RAY_HYPERSLAB_INT_3D


INTERFACE HDF5_READ
  MODULE PROCEDURE HDF5_READ_DOUBLE_1D
  MODULE PROCEDURE HDF5_READ_DOUBLE_2D
  MODULE PROCEDURE HDF5_READ_DOUBLE_3D
  MODULE PROCEDURE HDF5_READ_INTEGER_1D
  MODULE PROCEDURE HDF5_READ_INTEGER_2D
  MODULE PROCEDURE HDF5_READ_INTEGER_3D
END INTERFACE HDF5_READ


INTERFACE READ_1D_SLAB
  MODULE PROCEDURE READ_1D_SLAB_INT
  MODULE PROCEDURE READ_1D_SLAB_DOUBLE
END INTERFACE READ_1D_SLAB

INTERFACE READ_RAY_HYPERSLAB
  MODULE PROCEDURE READ_RAY_HYPERSLAB_DBL_2D
  MODULE PROCEDURE READ_RAY_HYPERSLAB_DBL_3D
  MODULE PROCEDURE READ_RAY_HYPERSLAB_DBL_4D
  MODULE PROCEDURE READ_RAY_HYPERSLAB_DBL_5D
  MODULE PROCEDURE READ_RAY_HYPERSLAB_INT_3D
END INTERFACE READ_RAY_HYPERSLAB


CONTAINS
!+101+###########################################################################!
!                                                                                !
!                  HDF5_OPENFILE                                                     !
!                                                                                !
!################################################################################!
SUBROUTINE HDF5_OPENFILE( FileName, NewFlag, file_id, hdferr )

CHARACTER(LEN = *), INTENT(IN)                        ::  FileName
LOGICAL, INTENT(IN)                                   ::  NewFlag
INTEGER(HID_T), INTENT(OUT)                           ::  file_id
INTEGER, INTENT(INOUT)                                ::  hdferr

IF ( NewFlag ) THEN

     CALL h5fcreate_f( TRIM( FileName ), H5F_ACC_TRUNC_F, file_id, hdferr )

ELSE

     CALL h5fopen_f( TRIM( FileName ), H5F_ACC_RDONLY_F, file_id, hdferr )

END IF 

END SUBROUTINE HDF5_OPENFILE





!+102+###########################################################################!
!                                                                                !
!                  HDF5_CLOSEFILE                                                !
!                                                                                !
!################################################################################!
SUBROUTINE HDF5_CLOSEFILE( file_id, hdferr )

INTEGER(HID_T), INTENT(IN)                            ::  file_id
INTEGER, INTENT(INOUT)                                ::  hdferr

CALL h5fclose_f( file_id, hdferr )


END SUBROUTINE HDF5_CLOSEFILE












!+103+###########################################################################!
!                                                                                !
!                  HDF5_OPENGROUP                                                !
!                                                                                !
!################################################################################!
SUBROUTINE HDF5_OPENGROUP( GroupName, NewFlag, file_id, group_id, hdferr )

CHARACTER(LEN = *), INTENT(IN)                        ::  GroupName
LOGICAL, INTENT(IN)                                   ::  NewFlag
INTEGER(HID_T), INTENT(IN)                            ::  file_id
INTEGER(HID_T), INTENT(OUT)                           ::  group_id
INTEGER, INTENT(INOUT)                                ::  hdferr

IF ( NewFlag ) THEN

     CALL h5gcreate_f( file_id, TRIM( GroupName ), group_id, hdferr )

ELSE

     CALL h5gopen_f( file_id, TRIM( GroupName ), group_id, hdferr )

END IF



END SUBROUTINE HDF5_OPENGROUP



!+104+###########################################################################!
!                                                                                !
!                  HDF5_CLOSEGROUP                                               !
!                                                                                !
!################################################################################!
SUBROUTINE HDF5_CLOSEGROUP( group_id, hdferr )

INTEGER(HID_T), INTENT(IN)                            ::  group_id
INTEGER, INTENT(INOUT)                                ::  hdferr

CALL h5gclose_f( group_id, hdferr )


END SUBROUTINE HDF5_CLOSEGROUP








!+201+###########################################################################!
!                                                                                !
!                  HDF5_READ_DOUBLE_1D                                          !
!                                                                                !
!################################################################################!
SUBROUTINE HDF5_READ_DOUBLE_1D( Name, Values, group_id, datasize, hdferr )

CHARACTER(LEN = *), INTENT(IN)                        ::  Name
REAL(KIND = idp), DIMENSION(:), INTENT(OUT)           ::  Values
INTEGER(HID_T), INTENT(IN)                            ::  group_id
INTEGER(HSIZE_T), DIMENSION(1), INTENT(IN)            ::  datasize
INTEGER, INTENT(INOUT)                                ::  hdferr

INTEGER(HID_T)                                        ::  dataset_id


CALL h5dopen_f( group_id, name, dataset_id, hdferr)
CALL h5dread_f( dataset_id, H5T_NATIVE_DOUBLE, Values, datasize, hdferr)
CALL h5dclose_f( dataset_id, hdferr ) 


END SUBROUTINE HDF5_READ_DOUBLE_1D







!+202+###########################################################################!
!                                                                                !
!                  HDF5_READ_DOUBLE_2D                                           !
!                                                                                !
!################################################################################!
SUBROUTINE HDF5_READ_DOUBLE_2D( Name, Values, group_id, datasize, hdferr )

CHARACTER(LEN = *), INTENT(IN)                        ::  Name
REAL(KIND = idp), DIMENSION(:,:), INTENT(OUT)         ::  Values
INTEGER(HID_T), INTENT(IN)                            ::  group_id
INTEGER(HSIZE_T), DIMENSION(2), INTENT(IN)            ::  datasize
INTEGER, INTENT(INOUT)                                ::  hdferr

INTEGER(HID_T)                                        ::  dataset_id


CALL h5dopen_f( group_id, name, dataset_id, hdferr)
CALL h5dread_f( dataset_id, H5T_NATIVE_DOUBLE, Values, datasize, hdferr)
CALL h5dclose_f( dataset_id, hdferr )


END SUBROUTINE HDF5_READ_DOUBLE_2D







!+203+###########################################################################!
!                                                                                !
!                  HDF5_READ_DOUBLE_3D                                           !
!                                                                                !
!################################################################################!
SUBROUTINE HDF5_READ_DOUBLE_3D( Name, Values, group_id, datasize, hdferr )

CHARACTER(LEN = *), INTENT(IN)                        ::  Name
REAL(KIND = idp), DIMENSION(:,:,:), INTENT(OUT)       ::  Values
INTEGER(HID_T), INTENT(IN)                            ::  group_id
INTEGER(HSIZE_T), DIMENSION(3), INTENT(IN)            ::  datasize
INTEGER, INTENT(INOUT)                                ::  hdferr

INTEGER(HID_T)                                        ::  dataset_id


CALL h5dopen_f( group_id, name, dataset_id, hdferr)
CALL h5dread_f( dataset_id, H5T_NATIVE_DOUBLE, Values, datasize, hdferr)
CALL h5dclose_f( dataset_id, hdferr )


END SUBROUTINE HDF5_READ_DOUBLE_3D








!+204+###########################################################################!
!                                                                                !
!                  HDF5_READ_INTEGER_1D                                          !
!                                                                                !
!################################################################################!
SUBROUTINE HDF5_READ_INTEGER_1D( Name, Values, group_id, datasize, hdferr )

CHARACTER(LEN = *), INTENT(IN)                        ::  Name
INTEGER, DIMENSION(:), INTENT(OUT)                    ::  Values
INTEGER(HID_T), INTENT(IN)                            ::  group_id
INTEGER(HSIZE_T), DIMENSION(1), INTENT(IN)            ::  datasize
INTEGER, INTENT(INOUT)                                ::  hdferr

INTEGER(HID_T)                                        ::  dataset_id


CALL h5dopen_f( group_id, name, dataset_id, hdferr)
CALL h5dread_f( dataset_id, H5T_NATIVE_INTEGER, Values, datasize, hdferr)
CALL h5dclose_f( dataset_id, hdferr )


END SUBROUTINE HDF5_READ_INTEGER_1D



!+205+###########################################################################!
!                                                                                !
!                  HDF5_READ_INTEGER_2D                                          !
!                                                                                !
!################################################################################!
SUBROUTINE HDF5_READ_INTEGER_2D( Name, Values, group_id, datasize, hdferr )

CHARACTER(LEN = *), INTENT(IN)                        ::  Name
INTEGER, DIMENSION(:,:), INTENT(OUT)                  ::  Values
INTEGER(HID_T), INTENT(IN)                            ::  group_id
INTEGER(HSIZE_T), DIMENSION(2), INTENT(IN)            ::  datasize
INTEGER, INTENT(INOUT)                                ::  hdferr

INTEGER(HID_T)                                        ::  dataset_id


CALL h5dopen_f( group_id, name, dataset_id, hdferr)
CALL h5dread_f( dataset_id, H5T_NATIVE_INTEGER, Values, datasize, hdferr)
CALL h5dclose_f( dataset_id, hdferr )


END SUBROUTINE HDF5_READ_INTEGER_2D




!+206+###########################################################################!
!                                                                                !
!                  HDF5_READ_INTEGER_3D                                          !
!                                                                                !
!################################################################################!
SUBROUTINE HDF5_READ_INTEGER_3D( Name, Values, group_id, datasize, hdferr )

CHARACTER(LEN = *), INTENT(IN)                        ::  Name
INTEGER, DIMENSION(:,:,:), INTENT(OUT)                ::  Values
INTEGER(HID_T), INTENT(IN)                            ::  group_id
INTEGER(HSIZE_T), DIMENSION(3), INTENT(IN)            ::  datasize
INTEGER, INTENT(INOUT)                                ::  hdferr

INTEGER(HID_T)                                        ::  dataset_id


CALL h5dopen_f( group_id, name, dataset_id, hdferr)
CALL h5dread_f( dataset_id, H5T_NATIVE_INTEGER, Values, datasize, hdferr)
CALL h5dclose_f( dataset_id, hdferr )


END SUBROUTINE HDF5_READ_INTEGER_3D









!+301+###########################################################################!
!                                                                                !
!                  READ_1D_SLAB_INT ( Same as HDF5_READ_DOUBLE_1D )              !
!                                                                                !
!################################################################################!
SUBROUTINE READ_1D_SLAB_INT(name, value, group_id, datasize )

CHARACTER(*), INTENT(IN)                        :: name
INTEGER(HID_T), INTENT(IN)                      :: group_id
INTEGER(HSIZE_T), DIMENSION(1), INTENT(IN)      :: datasize

INTEGER, DIMENSION(:), INTENT(OUT)              :: value

INTEGER(HID_T)                                  :: dataset_id
INTEGER                                         :: error

!  Open Dataset  !
CALL h5dopen_f(group_id, name, dataset_id, error)

!  Read Dataset  !
CALL h5dread_f(dataset_id, H5T_NATIVE_INTEGER, value, datasize, error)

! Close Dataset  !
CALL h5dclose_f(dataset_id, error)


END SUBROUTINE READ_1D_SLAB_INT





!+301+###########################################################################!
!                                                                                !
!                  READ_1D_SLAB_DOUBLE ( Same as HDF5_READ_DOUBLE_1D )           !
!                                                                                !
!################################################################################!
SUBROUTINE READ_1D_SLAB_DOUBLE(name, value, group_id, datasize )

CHARACTER(*), INTENT(IN)                        :: name
INTEGER(HID_T), INTENT(IN)                      :: group_id
INTEGER(HSIZE_T), DIMENSION(1), INTENT(IN)      :: datasize

REAL(KIND = idp), DIMENSION(:), INTENT(OUT)     :: value

INTEGER(HID_T)                                  :: dataset_id
INTEGER                                         :: error

!  Open Dataset  !
CALL h5dopen_f(group_id, name, dataset_id, error)

!  Read Dataset  !
CALL h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, value, datasize, error)

! Close Dataset  !
CALL h5dclose_f(dataset_id, error)


END SUBROUTINE READ_1D_SLAB_DOUBLE













!+401+###########################################################################!
!                                                                                !
!                  READ_RAY_HYPERSLAB_DBL_2D                                     !
!                                                                                !
!################################################################################!
SUBROUTINE READ_RAY_HYPERSLAB_DBL_2D( name, value, group_id, datasize, slab_offset)

CHARACTER(*), INTENT(IN)                                 :: name
INTEGER(HID_T), INTENT(IN)                               :: group_id
INTEGER(HSIZE_T), DIMENSION(2), INTENT(IN)               :: datasize
INTEGER(HSIZE_T), DIMENSION(2), INTENT(IN)               :: slab_offset
REAL(KIND = idp), DIMENSION(:,:), INTENT(OUT)            :: value


INTEGER(HID_T)                                           :: filespace
INTEGER(HID_T)                                           :: memspace
INTEGER(HID_T)                                           :: dataset_id
INTEGER(HSIZE_T), DIMENSION(2)                           :: null_offset
INTEGER                                                  :: error

null_offset = 0

CALL h5dopen_f(group_id, name, dataset_id, error)
CALL h5dget_space_f(dataset_id, filespace, error)
CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, slab_offset, datasize, error)
CALL h5screate_simple_f(2, datasize, memspace, error)
CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, null_offset, datasize, error)
CALL h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, value, datasize,     &
               error, file_space_id=filespace, mem_space_id=memspace)
CALL h5sclose_f(filespace, error)
CALL h5sclose_f(memspace, error)
CALL h5dclose_f(dataset_id, error)

END SUBROUTINE READ_RAY_HYPERSLAB_DBL_2D



!+402+###########################################################################!
!                                                                                !
!                  READ_RAY_HYPERSLAB_DBL_3D                                     !
!                                                                                !
!################################################################################!
SUBROUTINE READ_RAY_HYPERSLAB_DBL_3D( name, value, group_id, datasize, slab_offset)

CHARACTER(*), INTENT(IN)                                 :: name
INTEGER(HID_T), INTENT(IN)                               :: group_id
INTEGER(HSIZE_T), DIMENSION(3), INTENT(IN)               :: datasize
INTEGER(HSIZE_T), DIMENSION(3), INTENT(IN)               :: slab_offset
REAL(KIND = idp), DIMENSION(:,:,:), INTENT(OUT)          :: value


INTEGER(HID_T)                                           :: filespace
INTEGER(HID_T)                                           :: memspace
INTEGER(HID_T)                                           :: dataset_id
INTEGER(HSIZE_T), DIMENSION(3)                           :: null_offset
INTEGER                                                  :: error

null_offset = 0

CALL h5dopen_f(group_id, name, dataset_id, error)
CALL h5dget_space_f(dataset_id, filespace, error)
CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, slab_offset, datasize, error)
CALL h5screate_simple_f(3, datasize, memspace, error)
CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, null_offset, datasize, error)
CALL h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, value, datasize,     &
               error, file_space_id=filespace, mem_space_id=memspace)
CALL h5sclose_f(filespace, error)
CALL h5sclose_f(memspace, error)
CALL h5dclose_f(dataset_id, error)

END SUBROUTINE READ_RAY_HYPERSLAB_DBL_3D


!+403+###########################################################################!
!                                                                                !
!                  READ_RAY_HYPERSLAB_DBL_4D                                     !
!                                                                                !
!################################################################################!
SUBROUTINE READ_RAY_HYPERSLAB_DBL_4D( name, value, group_id, datasize, slab_offset)

CHARACTER(*), INTENT(IN)                                 :: name
INTEGER(HID_T), INTENT(IN)                               :: group_id
INTEGER(HSIZE_T), DIMENSION(4), INTENT(IN)               :: datasize
INTEGER(HSIZE_T), DIMENSION(4), INTENT(IN)               :: slab_offset
REAL(KIND = idp), DIMENSION(:,:,:,:), INTENT(OUT)        :: value


INTEGER(HID_T)                                           :: filespace
INTEGER(HID_T)                                           :: memspace
INTEGER(HID_T)                                           :: dataset_id
INTEGER(HSIZE_T), DIMENSION(4)                           :: null_offset
INTEGER                                                  :: error

null_offset = 0

CALL h5dopen_f(group_id, name, dataset_id, error)
CALL h5dget_space_f(dataset_id, filespace, error)
CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, slab_offset, datasize, error)
CALL h5screate_simple_f(4, datasize, memspace, error)
CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, null_offset, datasize, error)
CALL h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, value, datasize,     &
               error, file_space_id=filespace, mem_space_id=memspace)
CALL h5sclose_f(filespace, error)
CALL h5sclose_f(memspace, error)
CALL h5dclose_f(dataset_id, error)

END SUBROUTINE READ_RAY_HYPERSLAB_DBL_4D




!+404+###########################################################################!
!                                                                                !
!                  READ_RAY_HYPERSLAB_DBL_5D                                     !
!                                                                                !
!################################################################################!
SUBROUTINE READ_RAY_HYPERSLAB_DBL_5D( name, value, group_id, datasize, slab_offset)

CHARACTER(*), INTENT(IN)                                 :: name
INTEGER(HID_T), INTENT(IN)                               :: group_id
INTEGER(HSIZE_T), DIMENSION(5), INTENT(IN)               :: datasize
INTEGER(HSIZE_T), DIMENSION(5), INTENT(IN)               :: slab_offset
REAL(KIND = idp), DIMENSION(:,:,:,:,:), INTENT(OUT)      :: value


INTEGER(HID_T)                                           :: filespace
INTEGER(HID_T)                                           :: memspace
INTEGER(HID_T)                                           :: dataset_id
INTEGER(HSIZE_T), DIMENSION(5)                           :: null_offset
INTEGER                                                  :: error

null_offset = 0

CALL h5dopen_f(group_id, name, dataset_id, error)
CALL h5dget_space_f(dataset_id, filespace, error)
CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, slab_offset, datasize, error)
CALL h5screate_simple_f(5, datasize, memspace, error)
CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, null_offset, datasize, error)
CALL h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, value, datasize,     &
               error, file_space_id=filespace, mem_space_id=memspace)
CALL h5sclose_f(filespace, error)
CALL h5sclose_f(memspace, error)
CALL h5dclose_f(dataset_id, error)

END SUBROUTINE READ_RAY_HYPERSLAB_DBL_5D








!+403+###########################################################################!
!                                                                                !
!                  READ_RAY_HYPERSLAB_INT_3D                                     !
!                                                                                !
!################################################################################!
SUBROUTINE READ_RAY_HYPERSLAB_INT_3D( name, value, group_id, datasize, slab_offset)

CHARACTER(*), INTENT(IN)                                 :: name
INTEGER(HID_T), INTENT(IN)                               :: group_id
INTEGER(HSIZE_T), DIMENSION(3), INTENT(IN)               :: datasize
INTEGER(HSIZE_T), DIMENSION(3), INTENT(IN)               :: slab_offset
INTEGER, DIMENSION(:,:,:), INTENT(OUT)                   :: value


INTEGER(HID_T)                                           :: filespace
INTEGER(HID_T)                                           :: memspace
INTEGER(HID_T)                                           :: dataset_id
INTEGER(HSIZE_T), DIMENSION(3)                           :: null_offset
INTEGER                                                  :: error

null_offset = 0

CALL h5dopen_f(group_id, name, dataset_id, error)
CALL h5dget_space_f(dataset_id, filespace, error)
CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, slab_offset, datasize, error)
CALL h5screate_simple_f(3, datasize, memspace, error)
CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, null_offset, datasize, error)
CALL h5dread_f(dataset_id, H5T_NATIVE_INTEGER, value, datasize,     &
               error, file_space_id=filespace, mem_space_id=memspace)
CALL h5sclose_f(filespace, error)
CALL h5sclose_f(memspace, error)
CALL h5dclose_f(dataset_id, error)

END SUBROUTINE READ_RAY_HYPERSLAB_INT_3D








END MODULE HDF5_IO_Module
