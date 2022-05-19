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



USE Variables_AMReX_Core, &
            ONLY :  AMReX_Max_Grid_Size,    &
                    AMReX_Max_Level,        &
                    AMReX_Num_Levels,       &
                    AMReX_Tiling


USE Poseidon_Parameters, &
            ONLY :  Degree,                         &
                    L_Limit,                        &
                    Verbose_Flag,                   &
                    Convergence_Criteria,           &
                    Convergence_Criteria_Default,   &
                    Max_Iterations,                 &
                    Max_Iterations_Default

USE Variables_FP, &
            ONLY :  FP_Anderson_M,          &
                    FP_Anderson_M_Default


USE Variables_Mesh, &
            ONLY :  Num_R_Elements,         &
                    Num_T_Elements,         &
                    Num_P_Elements,         &
                    R_Inner,                &
                    R_Outer


#ifdef POSEIDON_AMREX_FLAG
use amrex_base_module
use amrex_fort_module,      ONLY: &
  amrex_spacedim
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
USE amrex_bc_types_module, ONLY: &
  amrex_bc_foextrap, &
  amrex_bc_bogus

USE Variables_Driver_AMReX, &
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
                    UseTiling,          &
                    StepNo,             &
                    t_new,              &
                    t_old,              &
                    dt,                 &
                    lo_Bc,              &
                    hi_bc

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








!+301+###########################################################################!
!                                                                                !
!                  Init_AMReX_Parameters                                             !
!                                                                                !
!################################################################################!
SUBROUTINE Init_AMReX_Parameters_From_Input_File()


#ifdef POSEIDON_AMREX_FLAG
TYPE(amrex_parmparse)                       :: PP

REAL(idp),  ALLOCATABLE                     :: xL_In(:), xR_In(:)
INTEGER,    ALLOCATABLE                     :: nCells_In(:)

ALLOCATE( xL_In(3), xR_In(3) )
ALLOCATE( nCells_In(3) )


IF ( Verbose_Flag ) THEN
    WRITE(*,'(A,I2.2)')"-Initializing Poseidon variables from inputs file."
END IF



CALL amrex_parmparse_build( PP, 'geometry' )
    CALL PP % getarr( 'prob_lo'  , xL_In   )
    CALL PP % getarr( 'prob_hi'  , xR_In   )
CALL amrex_parmparse_destroy( PP )


R_Inner = xL_In(1)
R_Outer = xR_In(1)



MaxGridSizeX1    = 1
MaxGridSizeX2    = 1
MaxGridSizeX3    = 1
BlockingFactorX1 = 1
BlockingFactorX2 = 1
BlockingFactorX3 = 1
UseTiling        = .FALSE.
CALL amrex_parmparse_build( PP, 'amr' )
  CALL PP % getarr( 'n_cell'           , nCells_In              )
  CALL PP % query ( 'max_grid_size_x'  , AMReX_Max_Grid_Size(1) )
  CALL PP % query ( 'max_grid_size_y'  , AMReX_Max_Grid_Size(2) )
  CALL PP % query ( 'max_grid_size_z'  , AMReX_Max_Grid_Size(3) )
  CALL PP % get   ( 'max_level'        , AMReX_Max_Level        )
  CALL PP % query ( 'UseTiling'        , AMReX_Tiling           )
CALL amrex_parmparse_destroy( PP )


Num_R_Elements = nCells_In(1)
Num_T_Elements = nCells_In(2)
Num_P_Elements = nCells_In(3)

AMReX_Num_Levels = AMReX_Max_Level+1


Degree = 1
L_Limit = 0
Max_Iterations = 10
FP_Anderson_M = FP_Anderson_M_Default
Max_Iterations = Max_Iterations_Default
Convergence_Criteria = Convergence_Criteria_Default
CALL amrex_parmparse_build( PP, 'poseidon' )
  CALL PP % query ( 'fem_degree'        , Degree )
  CALL PP % query ( 'l_limit'           , L_Limit               )
  CALL PP % query ( 'max_fp_iters'      , Max_Iterations        )
  CALL PP % query ( 'anderson_m'        , FP_Anderson_M         )
  CALL PP % query ( 'converge_criteria' , Convergence_Criteria  )
CALL amrex_parmparse_destroy( PP )



#endif




END SUBROUTINE Init_AMReX_Parameters_From_Input_File



END MODULE Poseidon_AMReX_Input_Parsing_Module
