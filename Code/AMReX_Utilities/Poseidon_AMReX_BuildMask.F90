   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Poseidon_AMReX_BuildMask_Module                                       !##!
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

USE amrex_box_module, &
            ONLY:   amrex_box

USE amrex_boxarray_module, &
            ONLY:   amrex_boxarray

USE amrex_distromap_module, &
            ONLY:   amrex_distromap

USE amrex_multifab_module, &
            ONLY:   amrex_imultifab
#endif


IMPLICIT NONE


INTERFACE

    !+101+##########################################################################!
    !                                                                               !
    !                                                                                   !
    !                                                                               !
    !###############################################################################!
    SUBROUTINE amrex_fi_buildmask(  Mask,                   &
                                    Geom,                   &
                                    C_Covered,              &
                                    C_NotCovered,           &
                                    C_PhysBnd,              &
                                    C_Interior              ) BIND(c)

        IMPORT
        IMPLICIT NONE
        type(c_ptr)                                     :: Mask
        type(c_ptr), VALUE                              :: Geom

        INTEGER(c_int), VALUE                           :: C_Covered
        INTEGER(c_int), VALUE                           :: C_NotCovered
        INTEGER(c_int), VALUE                           :: C_PhysBnd
        INTEGER(c_int), VALUE                           :: C_Interior

    END SUBROUTINE amrex_fi_buildmask

END INTERFACE

CONTAINS

#ifdef POSEIDON_AMREX_FLAG
!+101+##########################################################################!
!                                                                               !
!                                                                     !
!                                                                               !
!###############################################################################!
SUBROUTINE AMReX_BuildMask( Mask,                   &
                            Geom,                   &
                            C_Covered,              &
                            C_NotCovered,           &
                            C_PhysBnd,              &
                            C_Interior              )

    type(amrex_imultifab),  INTENT(INOUT)           :: Mask
    type(amrex_geometry),   INTENT(IN)              :: Geom

    INTEGER(c_int), INTENT(IN)                      :: C_Covered
    INTEGER(c_int), INTENT(IN)                      :: C_NotCovered
    INTEGER(c_int), INTENT(IN)                      :: C_PhysBnd
    INTEGER(c_int), INTENT(IN)                      :: C_Interior

    
    Mask%owner  = .TRUE.
    Mask%nc     = 1
!    Mask%ng     = 0
    CALL amrex_fi_buildmask(    Mask%p,             &
                                Geom%p,             &
                                C_Covered,          &
                                C_NotCovered,       &
                                C_PhysBnd,          &
                                C_Interior          )



END SUBROUTINE AMReX_BuildMask



#else




!+101+##########################################################################!
!                                                                               !
!                                                                     !
!                                                                               !
!###############################################################################!
SUBROUTINE AMReX_BuildMask(   )


END SUBROUTINE AMReX_BuildMask


#endif



END MODULE Poseidon_AMReX_BuildMask_Module

