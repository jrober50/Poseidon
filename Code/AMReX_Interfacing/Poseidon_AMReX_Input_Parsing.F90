   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Poseidon_AMReX_Input_Parsing_Module                                   !##!
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
USE Poseidon_Kinds_Module, &
               ONLY :  idp

#ifdef POSEIDON_AMREX_FLAG
use amrex_base_module

USE amrex_box_module,       ONLY: &
  amrex_box
USE amrex_boxarray_module,  ONLY: &
  amrex_boxarray,         &
  amrex_boxarray_build,   &
  amrex_boxarray_destroy
USE amrex_distromap_module, ONLY: &
  amrex_distromap,       &
  amrex_distromap_build, &
  amrex_distromap_destroy
USE amrex_multifab_module,  ONLY: &
  amrex_multifab, &
  amrex_multifab_build, &
  amrex_mfiter, &
  amrex_mfiter_build, &
  amrex_mfiter_destroy
USE amrex_parmparse_module, ONLY: &
  amrex_parmparse,       &
  amrex_parmparse_build, &
  amrex_parmparse_destroy

USE Variables_AMReX_Multifabs, &
        ONLY :  xL,                 &
                xR,                 &
                coord_sys,          &
                nCells,             &
                nLevels,            &
                MaxLevel,           &
                MaxGridSizeX1,      &
                MaxGridSizeX2,      &
                MaxGridSizeX3,      &
                MaxGridSizeX,       &
                BlockingFactorX1,   &
                BlockingFactorX2,   &
                BlockingFactorX3,   &
                UseTiling


#endif

IMPLICIT NONE



CONTAINS


!+301+###########################################################################!
!                                                                                !
!                  Init_AMReX_Parameters                                             !
!                                                                                !
!################################################################################!
SUBROUTINE Init_AMReX_Parameters()


#ifdef POSEIDON_AMREX_FLAG
TYPE(amrex_parmparse)                       :: PP

ALLOCATE( xL(3), xR(3) )
ALLOCATE( nCells(3) )


CALL amrex_parmparse_build( PP, 'geometry' )
    CALL PP % get   ( 'coord_sys', coord_sys )
    CALL PP % getarr( 'prob_lo'  , xL        )
    CALL PP % getarr( 'prob_hi'  , xR        )
CALL amrex_parmparse_destroy( PP )

MaxGridSizeX1    = 1
MaxGridSizeX2    = 1
MaxGridSizeX3    = 1
BlockingFactorX1 = 1
BlockingFactorX2 = 1
BlockingFactorX3 = 1
UseTiling        = .FALSE.
CALL amrex_parmparse_build( PP, 'amr' )
  CALL PP % getarr( 'n_cell'           , nCells           )
  CALL PP % query ( 'max_grid_size_x'  , MaxGridSizeX1    )
  CALL PP % query ( 'max_grid_size_y'  , MaxGridSizeX2    )
  CALL PP % query ( 'max_grid_size_z'  , MaxGridSizeX3    )
  CALL PP % query ( 'blocking_factor_x', BlockingFactorX1 )
  CALL PP % query ( 'blocking_factor_y', BlockingFactorX2 )
  CALL PP % query ( 'blocking_factor_z', BlockingFactorX3 )
  CALL PP % get   ( 'max_level'        , MaxLevel         )
  CALL PP % query ( 'UseTiling'        , UseTiling        )
CALL amrex_parmparse_destroy( PP )

MaxGridSizeX = [ MaxGridSizeX1, MaxGridSizeX2, MaxGridSizeX3 ]
nLevels = MaxLevel+1
#endif


END SUBROUTINE Init_AMReX_Parameters





END MODULE Poseidon_AMReX_Input_Parsing_Module
