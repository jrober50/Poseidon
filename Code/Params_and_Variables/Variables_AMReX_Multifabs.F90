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
            ONLY : idp

#ifdef POSEIDON_AMREX_FLAG
USE amrex_multifab_module,  &
            ONLY: amrex_multifab

USE amrex_boxarray_module,  &
            ONLY: amrex_boxarray

USE amrex_distromap_module, &
            ONLY: amrex_distromap
#endif

IMPLICIT NONE


#ifdef POSEIDON_AMREX_FLAG
TYPE(amrex_multifab), ALLOCATABLE                   ::  MF_Source(:)
TYPE(amrex_boxarray), ALLOCATABLE                   ::  BA_Source(:)
TYPE(amrex_distromap), ALLOCATABLE                  ::  DM_Source(:)
#endif



END MODULE Variables_AMReX_Multifabs




