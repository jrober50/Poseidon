   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Allocation_Tables                                                            !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
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

USE Poseidon_Message_Routines_Module, &
            ONLY :  Init_Message

USE Poseidon_Parameters, &
            ONLY :  Degree,                 &
                    L_Limit,                &
                    Max_Iterations,         &
                    Verbose_Flag

USE Variables_Derived, &
            ONLY :  LM_Length

USE Variables_Mesh, &
            ONLY :  Num_R_Elements,         &
                    Num_T_Elements,         &
                    Num_P_Elements

USE Variables_Quadrature, &
            ONLY :  Num_R_Quad_Points,      &
                    Num_T_Quad_Points,      &
                    Num_P_Quad_Points,      &
                    Num_TP_Quad_Points

USE Variables_Tables, &
            ONLY :  Ylm_Table_Block,            &
                    Ylm_Values,                 &
                    Ylm_dt_Values,              &
                    Ylm_dp_Values,              &
                    Ylm_CC_Values,              &
                    Ylm_CC_dt_Values,           &
                    Ylm_CC_dp_Values,           &
                    Lagrange_Poly_Table,        &
                    LPT_LPT,                    &
                    M_Values,                   &
                    Ylm_Norm_Table,             &
                    Ylm_Sqrt_Table,             &
                    rBT_NormedLegendre,         &
!                    rBT_NormedLegendre_dt,  &
!                    rBT_NormedLegendre_CC,  &
                    Ylm_Elem_Values,            &
                    Ylm_Elem_dt_Values,         &
                    Ylm_Elem_dp_Values,         &
                    Ylm_Elem_CC_Values,         &
                    Level_dx,                   &
                    Level_Ratios,               &
                    LagPoly_MultiLayer_Table,   &
                    LagPoly_Num_Tables

USE Variables_AMReX_Core, &
            ONLY :  AMReX_Max_Grid_Size,        &
                    AMReX_Num_Levels

USE Variables_IO, &
            ONLY :  Frame_Update_Table,         &
                    Frame_Residual_Table,       &
                    Iteration_Histogram

USE Flags_Initialization_Module, &
            ONLY :  lPF_Init_Tables_Flags,    &
                    iPF_Init_Tables_Alloc


IMPLICIT NONE


CONTAINS




!+101+##########################################################################!
!                                                                               !
!                            Allocate_Mesh                                      !
!                                                                               !
!###############################################################################!
SUBROUTINE Allocate_Tables()

IF ( Verbose_Flag ) CALL Init_Message('Allocating Table Variables.')

#ifdef POSEIDON_AMREX_FLAG

ALLOCATE( Ylm_Norm_Table( -L_Limit:L_Limit, 0:L_Limit ) )
ALLOCATE( Ylm_Sqrt_Table( -L_Limit:L_Limit, 0:L_Limit ) )


ALLOCATE( rBT_NormedLegendre(   -L_Limit:L_Limit,           &
                                -1:L_Limit,                 &
                                1:Num_T_Quad_Points,        &
                                0:AMReX_Max_Grid_Size(2)-1 )   )

!ALLOCATE( rBT_NormedLegendre_dt(-L_Limit:L_Limit,           &
!                                -1:L_Limit,                 &
!                                1:Num_T_Quad_Points,        &
!                                0:MaxGridSizeX2-1       )   )
!
!ALLOCATE( rBT_NormedLegendre_CC(-L_Limit:L_Limit,           &
!                                -1:L_Limit,                 &
!                                1:Num_T_Quad_Points,        &
!                                0:MaxGridSizeX2-1       )   )


ALLOCATE( Ylm_Elem_Values(      1:LM_Length,                &
                                1:Num_TP_Quad_Points    )   )
ALLOCATE( Ylm_Elem_dt_Values(   1:LM_Length,                &
                                1:Num_TP_Quad_Points    )   )
ALLOCATE( Ylm_Elem_dp_Values(   1:LM_Length,                &
                                1:Num_TP_Quad_Points    )   )
ALLOCATE( Ylm_Elem_CC_Values(   1:Num_TP_Quad_Points,       &
                                1:LM_Length            )   )


ALLOCATE( LagPoly_MultiLayer_Table( 0:DEGREE,               &
                                    1:NUM_R_QUAD_POINTS,    &
                                    0:1,                    &
                                    0:LagPoly_Num_Tables-1  )   )

ALLOCATE( Level_dx( 0:AMReX_Num_Levels-1, 3 ) )
ALLOCATE( Level_Ratios(0:AMReX_Num_Levels) )


#else


ALLOCATE( Ylm_Values(       1:LM_Length,                    &
                            1:NUM_TP_QUAD_POINTS,           &
                            0:Num_T_Elements-1,      &
                            0:Num_P_Elements-1)      )

ALLOCATE( Ylm_dt_Values(    1:LM_Length,                    &
                            1:NUM_TP_QUAD_POINTS,           &
                            0:Num_T_Elements-1,      &
                            0:Num_P_Elements-1)      )

ALLOCATE( Ylm_dp_Values(    1:LM_Length,                    &
                            1:NUM_TP_QUAD_POINTS,           &
                            0:Num_T_Elements-1,      &
                            0:Num_P_Elements-1)      )

ALLOCATE( Ylm_CC_Values(    1:NUM_TP_QUAD_POINTS,           &
                            1:LM_Length,                    &
                            0:Num_T_Elements-1,      &
                            0:Num_P_Elements-1)      )

ALLOCATE( Ylm_CC_dt_Values( 1:NUM_TP_QUAD_POINTS,           &
                            1:LM_Length,                    &
                            0:Num_T_Elements-1,      &
                            0:Num_P_Elements-1)      )

ALLOCATE( Ylm_CC_dp_Values( 1:NUM_TP_QUAD_POINTS,           &
                            1:LM_Length,                    &
                            0:Num_T_Elements-1,      &
                            0:Num_P_Elements-1)      )


#endif

ALLOCATE( M_VALUES(0:L_LIMIT) )
ALLOCATE( Lagrange_Poly_Table(0:DEGREE, 1:NUM_R_QUAD_POINTS, 0:2)   )
ALLOCATE( LPT_LPT( 1:NUM_R_QUAD_POINTS,0:DEGREE,0:DEGREE,0:1,0:2)       )



ALLOCATE( Frame_Update_Table(1:MAX_ITERATIONS,1:5) )
ALLOCATE( Frame_Residual_Table(1:3, 1:MAX_ITERATIONS,1:5) )
ALLOCATE( Iteration_Histogram(1:MAX_ITERATIONS) )

Frame_Update_Table = 0.0_idp
Frame_Residual_Table = 0.0_idp
Iteration_Histogram = 0

lPF_Init_Tables_Flags(iPF_Init_Tables_Alloc) = .TRUE.

END SUBROUTINE Allocate_Tables









!+102+##########################################################################!
!                                                                               !
!                           Deallocate_Mesh                                     !
!                                                                               !
!###############################################################################!
SUBROUTINE Deallocate_Tables()

#ifdef POSEIDON_AMREX_FLAG

DEALLOCATE( Ylm_Norm_Table )
DEALLOCATE( Ylm_Sqrt_Table )
DEALLOCATE( rBT_NormedLegendre )

DEALLOCATE( Ylm_Elem_Values )
DEALLOCATE( Ylm_Elem_dt_Values )
DEALLOCATE( Ylm_Elem_dp_Values )
DEALLOCATE( Ylm_Elem_CC_Values )


DEALLOCATE( LagPoly_MultiLayer_Table )

DEALLOCATE( Level_dx )
DEALLOCATE( Level_Ratios )

#else

DEALLOCATE( Ylm_Values )
DEALLOCATE( Ylm_dt_Values )
DEALLOCATE( Ylm_dp_Values )

DEALLOCATE( Ylm_CC_Values )
DEALLOCATE( Ylm_CC_dt_Values )
DEALLOCATE( Ylm_CC_dp_Values )

#endif

DEALLOCATE( M_Values )


DEALLOCATE( Lagrange_Poly_Table )
DEALLOCATE( LPT_LPT )

DEALLOCATE( Frame_Update_Table )
DEALLOCATE( Frame_Residual_Table )
DEALLOCATE( ITERATION_HISTOGRAM )


END SUBROUTINE Deallocate_Tables





END MODULE Allocation_Tables


