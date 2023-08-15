   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Variables_AMReX_Core                                                         !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!


USE Poseidon_Kinds_Module,  &
            ONLY : idp

#ifdef POSEIDON_AMREX_FLAG
USE amrex_multifab_module,  &
            ONLY :  amrex_multifab

USE amrex_boxarray_module,  &
            ONLY :  amrex_boxarray

USE amrex_distromap_module, &
            ONLY :  amrex_distromap

USE amrex_geometry_module,  &
            ONLY :  amrex_geometry
#endif


IMPLICIT NONE






LOGICAL, PUBLIC                                 ::  AMReX_Tiling = .TRUE.

LOGICAL, PUBLIC                                 ::  AMReX_Regrid_Flag = .TRUE.

INTEGER, PUBLIC                                 ::  AMReX_MaxLevel      = -1
INTEGER, PUBLIC                                 ::  AMReX_Num_Levels    = -1
INTEGER, PUBLIC                                 ::  AMReX_Old_Levels    = -1
INTEGER, PUBLIC, DIMENSION(3)                   ::  AMReX_Max_Grid_Size = -1

INTEGER, PUBLIC                                 ::  AMReX_Tag_Style     = 3

INTEGER, PUBLIC                                 ::  iFRL    ! FEM Refinement Level
INTEGER, PUBLIC                                 ::  iFEME   ! Num FEM Elements

INTEGER, PUBLIC                                 ::  iIRL    ! Integral Refinement Level
INTEGER, PUBLIC                                 ::  iINTE   ! Num Integral Elements


INTEGER, PUBLIC                                 ::  AMReX_Dims

INTEGER, PUBLIC                                 ::  iNumLeafElements
INTEGER, PUBLIC, DIMENSION(:), ALLOCATABLE      ::  iLeafElementsPerLvl
INTEGER, PUBLIC, DIMENSION(:), ALLOCATABLE      ::  FindLoc_Table
INTEGER, PUBLIC, DIMENSION(:), ALLOCATABLE      ::  FEM_Elem_Table

INTEGER, PUBLIC, DIMENSION(:), ALLOCATABLE      ::  Table_Offsets


INTEGER, PUBLIC                                 ::  MF_Source_nComps


REAL(idp), CONTIGUOUS, POINTER                  :: Source_PTR(:,:,:,:)
INTEGER,   CONTIGUOUS, POINTER                  :: Mask_PTR(:,:,:,:)
INTEGER,   CONTIGUOUS, POINTER                  :: Ghost_PTR(:,:,:,:)


#ifdef POSEIDON_AMREX_FLAG

TYPE(amrex_multifab),  PUBLIC,  ALLOCATABLE     ::  MF_Source(:)
TYPE(amrex_boxarray),  PUBLIC,  ALLOCATABLE     ::  BA_Source(:)
TYPE(amrex_distromap), PUBLIC,  ALLOCATABLE     ::  DM_Source(:)
TYPE(amrex_geometry),  PUBLIC,  ALLOCATABLE     ::  GM_Source(:)

#endif




END MODULE Variables_AMReX_Core





