   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Variables_AMReX_Multifabs                                                    !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!


USE Poseidon_Kinds_Module,  &
            ONLY :  idp

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


#ifdef POSEIDON_AMREX_FLAG

REAL(idp)           , ALLOCATABLE           :: xL(:), xR(:)
INTEGER                                     :: coord_sys
INTEGER             , ALLOCATABLE           :: nCells(:)
INTEGER                                     :: nLevels
INTEGER                                     :: MaxLevel
INTEGER                                     :: MaxGridSizeX1
INTEGER                                     :: MaxGridSizeX2
INTEGER                                     :: MaxGridSizeX3
INTEGER                                     :: MaxGridSizeX(3)
INTEGER                                     :: BlockingFactorX1
INTEGER                                     :: BlockingFactorX2
INTEGER                                     :: BlockingFactorX3
LOGICAL                                     :: UseTiling


TYPE(amrex_multifab),  ALLOCATABLE                  ::  MF_Source(:)
TYPE(amrex_boxarray),  ALLOCATABLE                  ::  BA_Source(:)
TYPE(amrex_distromap), ALLOCATABLE                  ::  DM_Source(:)
TYPE(amrex_geometry),  ALLOCATABLE                  ::  Geom_Source(:)
#endif



END MODULE Variables_AMReX_Multifabs




