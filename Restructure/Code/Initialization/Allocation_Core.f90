   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Allocation_Core                                                              !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains the functions and subroutines associated with the stiffness matrix !##!
!##!        involved in the linear system formed by the expansions. This includes   !##!
!##!        functions that build the matrix in various storage formats.             !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   Allocate_Poseidon_Variables                                         !##!
!##!    +102+   Deallocate_Poseidon_Variables                                       !##!
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

USE Poseidon_Parameters, &
            ONLY :  DEGREE,                 &
                    L_LIMIT,                &
                    Max_Iterations

USE Variables_Quadrature, &
            ONLY :  Num_R_Quad_Points,      &
                    Num_T_Quad_Points,      &
                    Num_P_Quad_Points,      &
                    Num_TP_Quad_Points,     &
                    Local_Node_Locations

USE Variables_MPI, &
            ONLY :  Num_R_Elems_Per_Block,  &
                    Num_T_Elems_Per_Block,  &
                    Num_P_Elems_Per_Block

Use Variables_Derived, &
            ONLY :  Num_R_Nodes,            &
                    Prob_Dim,               &
                    Block_Prob_Dim,         &
                    Elem_Prob_Dim,          &
                    Elem_Prob_Dim_Sqr,      &
                    SubShell_Prob_Dim,      &
                    Num_Off_Diagonals,      &
                    LM_Length

USE Variables_IO, &
            ONLY :  Frame_Time_Table,       &
                    Frame_Update_Table,     &
                    Frame_Residual_Table,   &
                    Run_Time_Table,         &
                    Num_Timer_Calls,        &
                    Iteration_Histogram,    &
                    Iter_Time_Table

USE Variables_NR, &
            ONLY :  Update_Vector,          &
                    Coefficient_Vector,     &
                    Block_RHS_Vector,       &
                    Block_STF_Mat,          &
                    Block_Elem_STF_MatVec

USE Variables_Source, &
            ONLY :  Block_Source_E,             &
                    Block_Source_S,             &
                    Block_Source_Si





IMPLICIT NONE


CONTAINS



!+101+###########################################################################!
!                                                                                !
!                            Allocate_Poseidon_Variables                         !
!                                                                                !
!################################################################################!
SUBROUTINE Allocate_Poseidon_CFA_Variables()





ALLOCATE(Block_Source_E(    1:NUM_R_QUAD_POINTS,        &
                            1:NUM_T_QUAD_POINTS,        &
                            1:NUM_P_QUAD_POINTS,        &
                            0:NUM_R_ELEMS_PER_BLOCK-1,  &
                            0:NUM_T_ELEMS_PER_BLOCK-1,  &
                            0:NUM_P_ELEMS_PER_BLOCK-1   )   )

ALLOCATE(Block_Source_S(    1:NUM_R_QUAD_POINTS,        &
                            1:NUM_T_QUAD_POINTS,        &
                            1:NUM_P_QUAD_POINTS,        &
                            0:NUM_R_ELEMS_PER_BLOCK-1,  &
                            0:NUM_T_ELEMS_PER_BLOCK-1,  &
                            0:NUM_P_ELEMS_PER_BLOCK-1   )   )

ALLOCATE(Block_Source_Si(   1:NUM_R_QUAD_POINTS,        &
                            1:NUM_T_QUAD_POINTS,        &
                            1:NUM_P_QUAD_POINTS,        &
                            0:NUM_R_ELEMS_PER_BLOCK-1,  &
                            0:NUM_T_ELEMS_PER_BLOCK-1,  &
                            0:NUM_P_ELEMS_PER_BLOCK-1,  &
                            1:3          )   )



ALLOCATE( ITER_TIME_TABLE(1:NUM_TIMER_CALLS) )
ALLOCATE( FRAME_TIME_TABLE(1:NUM_TIMER_CALLS) )
ALLOCATE( RUN_TIME_TABLE(1:NUM_TIMER_CALLS) )

ITER_TIME_TABLE = 0.0_idp
FRAME_TIME_TABLE = 0.0_idp
RUN_TIME_TABLE = 0.0_idp

ALLOCATE( Frame_Update_Table(1:MAX_ITERATIONS) )
ALLOCATE( Frame_Residual_Table(1:MAX_ITERATIONS) )
ALLOCATE( Iteration_Histogram(1:MAX_ITERATIONS) )
Iteration_Histogram = 0

END SUBROUTINE Allocate_Poseidon_CFA_Variables











!+102+###########################################################################!
!                                                                                !
!                           Deallocate_Poseidon_Variables                        !
!                                                                                !
!################################################################################!
SUBROUTINE Deallocate_Poseidon_CFA_Variables()




DEALLOCATE( Block_Source_E )
DEALLOCATE( Block_Source_S )
DEALLOCATE( Block_Source_Si )



DEALLOCATE( ITER_TIME_TABLE )
DEALLOCATE( Frame_TIME_TABLE )
DEALLOCATE( RUN_TIME_TABLE  )

DEALLOCATE( Frame_Update_Table )
DEALLOCATE( Frame_Residual_Table )
DEALLOCATE( ITERATION_HISTOGRAM )



END SUBROUTINE Deallocate_Poseidon_CFA_Variables





END MODULE Allocation_Core
