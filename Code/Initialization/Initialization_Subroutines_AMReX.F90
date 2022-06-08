   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Initialization_Subroutines_AMReX                                      !##!
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

USE Poseidon_Message_Routines_Module, &
            ONLY :  Init_Message

USE Initialization_Subroutines, &
            ONLY :  Init_Expansion_Params,      &
                    Init_Fixed_Point_Params,    &
                    Init_Mesh_Params,           &
                    Init_AMReX_Params

USE Variables_AMReX_Core, &
            ONLY :  AMReX_Max_Grid_Size,    &
                    AMReX_Max_Level,        &
                    AMReX_Num_Levels,       &
                    AMReX_Tiling


USE Poseidon_Parameters, &
            ONLY :  Degree_Default,                 &
                    L_Limit_Default,                &
                    Verbose_Flag,                   &
                    Convergence_Criteria_Default,   &
                    Max_Iterations_Default

USE Variables_FP, &
            ONLY :  FP_Anderson_M_Default

USE Variables_Mesh, &
            ONLY :  Num_R_Elements,         &
                    Num_T_Elements,         &
                    Num_P_Elements,         &
                    R_Inner,                &
                    R_Outer

USE Variables_Interface, &
            ONLY :  Caller_R_Units


#ifdef POSEIDON_AMREX_FLAG
USE amrex_base_module

USE amrex_parmparse_module, &
            ONLY :  amrex_parmparse,            &
                    amrex_parmparse_build,      &
                    amrex_parmparse_destroy
#endif

IMPLICIT NONE


CONTAINS



!+301+###########################################################################!
!                                                                                !
!                  Init_AMReX_Parameters                                             !
!                                                                                !
!################################################################################!
SUBROUTINE Init_Parameters_From_AMReX_Input_File()


#ifdef POSEIDON_AMREX_FLAG

TYPE(amrex_parmparse)                       :: PP

REAL(idp),  DIMENSION(:),   ALLOCATABLE     ::  xL_In
REAL(idp),  DIMENSION(:),   ALLOCATABLE     ::  xR_In
INTEGER,    DIMENSION(:),   ALLOCATABLE     ::  nCells_In
INTEGER,    DIMENSION(:),   ALLOCATABLE     ::  Max_Grid_Size_In
INTEGER                                     ::  Max_Level_In
LOGICAL                                     ::  Tiling_In

INTEGER                                     ::  Degree_In
INTEGER                                     ::  L_Limit_In

INTEGER                                     ::  FP_Anderson_M_In
INTEGER                                     ::  Max_Iterations_In
REAL(idp)                                   ::  Convergence_Criteria_In


ALLOCATE( xL_In(3) )
ALLOCATE( xR_In(3) )
ALLOCATE( nCells_In(3) )
ALLOCATE( Max_Grid_Size_In(3) )


IF ( Verbose_Flag ) &
    CALL Init_Message('Initializing Poseidon variables from inputs file.')

xL_In = 1
xR_In = 1
CALL amrex_parmparse_build( PP, 'geometry' )
    CALL PP % getarr( 'prob_lo'  , xL_In   )
    CALL PP % getarr( 'prob_hi'  , xR_In   )
CALL amrex_parmparse_destroy( PP )


nCells_In               = 1
Max_Grid_Size_In        = 1
Max_Level_In            = 1
Tiling_In               = .FALSE.
CALL amrex_parmparse_build( PP, 'amr' )
    CALL PP % getarr( 'n_cell'           , nCells_In            )
    CALL PP % query ( 'max_grid_size_x'  , Max_Grid_Size_In(1)  )
    CALL PP % query ( 'max_grid_size_y'  , Max_Grid_Size_In(2)  )
    CALL PP % query ( 'max_grid_size_z'  , Max_Grid_Size_In(3)  )
    CALL PP % get   ( 'max_level'        , Max_Level_In         )
    CALL PP % query ( 'UseTiling'        , Tiling_In            )
CALL amrex_parmparse_destroy( PP )



Degree_In = Degree_Default
L_Limit_In = L_Limit_Default
Max_Iterations_In = 10
FP_Anderson_M_In = FP_Anderson_M_Default
Max_Iterations_In = Max_Iterations_Default
Convergence_Criteria_In = Convergence_Criteria_Default
CALL amrex_parmparse_build( PP, 'poseidon' )
    CALL PP % query ( 'fem_degree'        , Degree_In )
    CALL PP % query ( 'l_limit'           , L_Limit_In               )
    CALL PP % query ( 'max_fp_iters'      , Max_Iterations_In        )
    CALL PP % query ( 'anderson_m'        , FP_Anderson_M_In         )
    CALL PP % query ( 'converge_criteria' , Convergence_Criteria_In  )
CALL amrex_parmparse_destroy( PP )






CALL Init_Expansion_Params( Degree_In, L_Limit_In )

CALL Init_Fixed_Point_Params(   Max_Iterations_In,          &
                                Convergence_Criteria_In,    &
                                FP_Anderson_M_In            )

CALL Init_AMReX_Params( nCells_In,          &
                        Max_Level_In,       &
                        Max_Grid_Size_In    )


xL_In(1) = xL_In(1)*Caller_R_Units
xR_In(1) = xR_In(1)*Caller_R_Units

CALL Init_Mesh_Params(  nCells_In,  &
                        xL_In,      &
                        xR_In       )





#endif



END SUBROUTINE Init_Parameters_From_AMReX_Input_File











END MODULE Initialization_Subroutines_AMReX
