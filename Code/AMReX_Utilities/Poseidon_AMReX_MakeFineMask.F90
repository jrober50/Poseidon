   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Poseidon_AMReX_MakeFineMask_Module                                    !##!
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
    SUBROUTINE amrex_fi_makefinemask( Mask,                   &
                                      Coarse_BA,              &
                                      Coarse_DM,              &
                                      nGhost_Vec,             &
                                      Periodicity,            &
                                      Fine_BA,                &
                                      C_Coarse, C_Fine ) BIND(c)

        IMPORT
        IMPLICIT NONE
        type(c_ptr)                                     :: Mask
        type(c_ptr), VALUE                              :: Coarse_BA
        type(c_ptr), VALUE                              :: Coarse_DM
        INTEGER(c_int), DIMENSION(1:3)                  :: nGhost_Vec
        INTEGER(c_int), DIMENSION(1:3)                  :: Periodicity
        type(c_ptr), VALUE                              :: Fine_BA

        INTEGER(c_int), VALUE                           :: C_Coarse
        INTEGER(c_int), VALUE                           :: C_Fine

    END SUBROUTINE amrex_fi_makefinemask

END INTERFACE

CONTAINS

#ifdef POSEIDON_AMREX_FLAG
!+101+##########################################################################!
!                                                                               !
!                                                                     !
!                                                                               !
!###############################################################################!
SUBROUTINE AMReX_MakeFineMask(  Mask,                   &
                                Coarse_BA,              &
                                Coarse_DM,              &
                                nGhost_Vec,             &
                                Fine_BA,                &
                                C_Coarse,               &
                                C_Fine                  )

    type(amrex_imultifab), INTENT(INOUT)            :: Mask
    type(amrex_boxarray),  INTENT(IN)               :: Coarse_BA
    type(amrex_distromap), INTENT(IN)               :: Coarse_DM
    type(amrex_boxarray),  INTENT(IN)               :: Fine_BA

    INTEGER(c_int), DIMENSION(1:3), INTENT(IN)      :: nGhost_Vec
    INTEGER(c_int), INTENT(IN)                      :: C_Coarse
    INTEGER(c_int), INTENT(IN)                      :: C_Fine

    INTEGER(c_int), DIMENSION(1:3)                         :: Periodicity

    Periodicity = 0

    Mask%owner  = .TRUE.
    Mask%nc     = 1
!    Mask%ng     = 0
    CALL amrex_fi_makefinemask( Mask%p,                 &
                                Coarse_BA%p,            &
                                Coarse_DM%p,            &
                                nGhost_Vec,             &
                                Periodicity,            &
                                Fine_BA%p,              &
                                C_Coarse, C_Fine     )



    print*,"In AMReX_MakeFineMask",Mask%ng
END SUBROUTINE AMReX_MakeFineMask



#else




!+101+##########################################################################!
!                                                                               !
!                                                                     !
!                                                                               !
!###############################################################################!
SUBROUTINE AMReX_MakeFineMask(   )


END SUBROUTINE AMReX_MakeFineMask


#endif



END MODULE Poseidon_AMReX_MakeFineMask_Module
