   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Poseidon_AMReX_BoxArraySize_Module                                    !##!
!##!                                                                         !##!
!##!_________________________________________________________________________!##!
!##!                                                                         !##!
!##!                                                                         !##!
!##!=========================================================================!##!
!##!                                                                         !##!
!##!    Contains:                                                            !##!
!##!                                                                         !##!
!##!                                                                         !##!
!###############################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !##########################################################################!


!*D*================================!
!                                   !
!           Dependencies            !
!                                   !
!===================================!
USE iso_c_binding

#ifdef POSEIDON_AMREX_FLAG
use amrex_base_module

USE amrex_box_module,       ONLY: &
  amrex_box
USE amrex_boxarray_module,  ONLY: &
  amrex_boxarray
USE amrex_distromap_module, ONLY: &
  amrex_distromap
USE amrex_multifab_module,  ONLY: &
  amrex_multifab, &
  amrex_imultifab
#endif




IMPLICIT NONE


INTERFACE
!+101+##########################################################################!
!                                                                               !
!                                                     				!
!                                                                               !
!###############################################################################!
    PURE FUNCTION poseidon_amrex_fi_boxarraysize( BA ) BIND(c)

    IMPORT
    IMPLICIT NONE
    type(c_ptr), VALUE, INTENT(IN)      :: BA
    INTEGER(amrex_long)                 :: poseidon_amrex_fi_boxarraysize


    END FUNCTION poseidon_amrex_fi_boxarraysize


END INTERFACE


CONTAINS

#ifdef POSEIDON_AMREX_FLAG
!+101+##########################################################################!
!                                                                               !
!                                                                     !
!                                                                               !
!###############################################################################!
PURE FUNCTION AMReX_BoxArraySize(  BA ) result( Size )
    class(amrex_boxarray),  INTENT(IN)      :: BA
    INTEGER(AMReX_long)                     :: Size

    Size = poseidon_amrex_fi_boxarraysize( BA%p )

END FUNCTION AMReX_BoxArraySize

#else

!+101+##########################################################################!
!                                                                               !
!                                                                     !
!                                                                               !
!###############################################################################!
PURE FUNCTION AMReX_BoxArraySize(  )
INTEGER         :: AMReX_BoxArraySize

AMReX_BoxArraySize = 1

END FUNCTION AMReX_BoxArraySize

#endif




END MODULE Poseidon_AMReX_BoxArraySize_Module
