   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Allocate_Variables_Module                                                    !##!
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
USE Poseidon_Constants_Module, &
            ONLY : idp

USE Poseidon_Parameters, &
            ONLY :  DOMAIN_DIM,             &
                    DEGREE,                 &
                    L_LIMIT,                &
                    NUM_CFA_VARS,           &
                    DATA_DIST_MODE,         &
                    NUM_R_ELEMS_PER_SHELL,  &
                    NUM_SHELLS,             &
                    NUM_BLOCKS_PER_SHELL,   &
                    nPROCS_POSEIDON,        &
                    STF_MAPPING_FLAG,       &
                    NUM_R_ELEMS_PER_BLOCK,  &
                    NUM_T_ELEMS_PER_BLOCK,  &
                    NUM_P_ELEMS_PER_BLOCK,  &
                    NUM_R_QUAD_POINTS,      &
                    NUM_T_QUAD_POINTS,      &
                    NUM_P_QUAD_POINTS,      &
                    MAX_ITERATIONS


USE Poseidon_Variables_Module, &
            ONLY :  NUM_R_ELEMENTS,             &
                    NUM_T_ELEMENTS,             &
                    NUM_P_ELEMENTS,             &
                    NUM_TP_QUAD_POINTS,         &
                    NUM_R_NODES,                &
                    PROB_DIM,                   &
                    Block_PROB_DIM,             &
                    SUBSHELL_PROB_DIM,          &
                    PHYSICS_TYPE,               &
                    NUM_OFF_DIAGONALS,          &
                    rlocs,                      &
                    tlocs,                      &
                    plocs,                      &
                    RHS_Vector,                 &
                    Coefficient_Vector,         &
                    Source_Term_Coefficients,   &
                    Block_Source_E,             &
                    Block_Source_S,             &
                    Block_Source_Si,            &
                    Update_Vector,              &
                    LOCAL_NODE_LOCATIONS,       &
                    INT_R_LOCATIONS,            &
                    INT_R_WEIGHTS,              &
                    INT_T_LOCATIONS,            &
                    INT_T_WEIGHTS,              &
                    INT_P_LOCATIONS,            &
                    INT_P_WEIGHTS,              &
                    INT_TP_WEIGHTS,             &
                    Ylm_Table_Block,            &
                    Ylm_Values,                 &
                    Ylm_dt_Values,              &
                    Ylm_dp_Values,              &
                    Ylm_CC_Values,          &
                    Ylm_CC_dt_Values,       &
                    Ylm_CC_dp_Values,       &
                    LM_Length,                  &
                    M_VALUES,                   &
                    Lagrange_Poly_Table,        &
                    LPT_LPT,                    &
                    Block_RHS_Vector,           &
                    Block_STF_Mat,              &
                    BLOCK_ELEM_STF_MATVEC,      &
                    Elem_PROB_DIM,              &
                    Elem_PROB_DIM_SQR,          &
                    Block_Prob_Dim,             &
                    Iter_Time_Table,            &
                    Frame_Time_Table,           &
                    Frame_Convergence_Table,    &
                    Iteration_Histogram,        &
                    Run_Time_Table,             &
                    Num_Timer_Calls,            &
                    STF_NNZ,                    &
                    STF_MAT,                    &
                    STF_ELEM_VAL,               &
                    STF_COL_PTR,                &
                    STF_ROW_IND






IMPLICIT NONE


CONTAINS



!+101+###########################################################################!
!                                                                                !
!                            Allocate_Poseidon_Variables                         !
!                                                                                !
!################################################################################!
SUBROUTINE Allocate_Poseidon_CFA_Variables()



ALLOCATE(rlocs(0:NUM_R_ELEMENTS))
ALLOCATE(tlocs(0:NUM_T_ELEMENTS))
ALLOCATE(plocs(0:NUM_P_ELEMENTS))


!ALLOCATE( RHS_Vector( 0:PROB_DIM-1 ) )
ALLOCATE( Update_Vector(0:PROB_DIM-1 ) )
ALLOCATE( Coefficient_Vector(0:PROB_DIM-1 ) )


ALLOCATE( Block_RHS_Vector( 0:Block_PROB_DIM-1) )
ALLOCATE( Block_STF_Mat( 0:2*NUM_OFF_DIAGONALS, 0:SUBSHELL_PROB_DIM-1) )


ALLOCATE( BLOCK_ELEM_STF_MATVEC(0:ELEM_PROB_DIM_SQR-1 ,0:NUM_R_ELEMS_PER_BLOCK-1) )



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






ALLOCATE( LOCAL_NODE_LOCATIONS(0:DEGREE) )



ALLOCATE( M_VALUES(0:L_LIMIT) )



!ALLOCATE( Ylm_Table_Block(  -L_LIMIT:L_LIMIT, -1:L_LIMIT,               &
!                            1:NUM_T_QUAD_POINTS, 1:NUM_P_QUAD_POINTS,   &
!                            0:NUM_T_ELEMS_PER_BLOCK-1, 0:NUM_P_ELEMS_PER_BLOCK-1)     )



ALLOCATE( Ylm_Values(       0:LM_Length-1,                                          &
                            1:NUM_TP_QUAD_POINTS,                                   &
                            0:NUM_T_ELEMS_PER_BLOCK-1, 0:NUM_P_ELEMS_PER_BLOCK-1)   )

ALLOCATE( Ylm_dt_Values(    0:LM_Length-1,                                          &
                            1:NUM_TP_QUAD_POINTS,                                   &
                            0:NUM_T_ELEMS_PER_BLOCK-1, 0:NUM_P_ELEMS_PER_BLOCK-1)   )

ALLOCATE( Ylm_dp_Values(    0:LM_Length-1,                                          &
                            1:NUM_TP_QUAD_POINTS,                                   &
                            0:NUM_T_ELEMS_PER_BLOCK-1, 0:NUM_P_ELEMS_PER_BLOCK-1)   )

ALLOCATE( Ylm_CC_Values(    1:NUM_TP_QUAD_POINTS,                                   &
                            0:LM_Length-1,                                          &
                            0:NUM_T_ELEMS_PER_BLOCK-1, 0:NUM_P_ELEMS_PER_BLOCK-1)   )

ALLOCATE( Ylm_CC_dt_Values( 1:NUM_TP_QUAD_POINTS,                                   &
                            0:LM_Length-1,                                          &
                            0:NUM_T_ELEMS_PER_BLOCK-1, 0:NUM_P_ELEMS_PER_BLOCK-1)   )

ALLOCATE( Ylm_CC_dp_Values( 1:NUM_TP_QUAD_POINTS,                                   &
                            0:LM_Length-1,                                          &
                            0:NUM_T_ELEMS_PER_BLOCK-1, 0:NUM_P_ELEMS_PER_BLOCK-1)   )


ALLOCATE( Lagrange_Poly_Table(0:DEGREE, 1:NUM_R_QUAD_POINTS, 0:2)   )
ALLOCATE( LPT_LPT( 1:NUM_R_QUAD_POINTS,0:DEGREE,0:DEGREE,0:1,0:2)       )

ALLOCATE(INT_R_LOCATIONS(1:NUM_R_QUAD_POINTS), INT_R_WEIGHTS(1:NUM_R_QUAD_POINTS))
ALLOCATE(INT_T_LOCATIONS(1:NUM_T_QUAD_POINTS), INT_T_WEIGHTS(1:NUM_T_QUAD_POINTS))
ALLOCATE(INT_P_LOCATIONS(1:NUM_P_QUAD_POINTS), INT_P_WEIGHTS(1:NUM_P_QUAD_POINTS))

ALLOCATE( INT_TP_WEIGHTS(1:NUM_TP_QUAD_POINTS) )


ALLOCATE( ITER_TIME_TABLE(1:NUM_TIMER_CALLS) )
ALLOCATE( FRAME_TIME_TABLE(1:NUM_TIMER_CALLS) )
ALLOCATE( RUN_TIME_TABLE(1:NUM_TIMER_CALLS) )

ITER_TIME_TABLE = 0.0_idp
FRAME_TIME_TABLE = 0.0_idp
RUN_TIME_TABLE = 0.0_idp

ALLOCATE( FRAME_CONVERGENCE_TABLE(1:MAX_ITERATIONS))
ALLOCATE( Iteration_Histogram(1:MAX_ITERATIONS) )
Iteration_Histogram = 0

END SUBROUTINE Allocate_Poseidon_CFA_Variables











!+102+###########################################################################!
!                                                                                !
!                           Deallocate_Poseidon_Variables                        !
!                                                                                !
!################################################################################!
SUBROUTINE Deallocate_Poseidon_CFA_Variables()


DEALLOCATE(rlocs)
DEALLOCATE(tlocs)
DEALLOCATE(plocs)

!DEALLOCATE( RHS_Vector )
DEALLOCATE( Update_Vector )
DEALLOCATE( Coefficient_Vector )


DEALLOCATE( Block_RHS_Vector )
DEALLOCATE( Block_STF_Mat )
DEALLOCATE( Block_ELEM_STF_MATVEC )


DEALLOCATE( Lagrange_Poly_Table )
DEALLOCATE( LPT_LPT )

DEALLOCATE(INT_R_LOCATIONS, INT_R_WEIGHTS)
DEALLOCATE(INT_T_LOCATIONS, INT_T_WEIGHTS)
DEALLOCATE(INT_P_LOCATIONS, INT_P_WEIGHTS)

DEALLOCATE( LOCAL_NODE_LOCATIONS )
DEALLOCATE( M_VALUES )

DEALLOCATE( Block_Source_E )
DEALLOCATE( Block_Source_S )
DEALLOCATE( Block_Source_Si )

!DEALLOCATE( Ylm_Table_Block )
DEALLOCATE( Ylm_Values )
DEALLOCATE( Ylm_CC_Values )
DEALLOCATE( Ylm_dt_Values )
DEALLOCATE( Ylm_dp_Values)



DEALLOCATE( ITER_TIME_TABLE )
DEALLOCATE( RUN_TIME_TABLE  )



DEALLOCATE( ITERATION_HISTOGRAM )



END SUBROUTINE Deallocate_Poseidon_CFA_Variables













!+201+###########################################################################!
!                                                                                !
!                     Allocate_Poseidon_Newtonian_Variables                      !
!                                                                                !
!################################################################################!
SUBROUTINE Allocate_Poseidon_Newtonian_Variables()



ALLOCATE(rlocs(0:NUM_R_ELEMENTS))
ALLOCATE(tlocs(0:NUM_T_ELEMENTS))
ALLOCATE(plocs(0:NUM_P_ELEMENTS))



ALLOCATE( RHS_Vector( 0:PROB_DIM-1 ) )
ALLOCATE( Update_Vector(0:PROB_DIM-1 ) )
ALLOCATE( Coefficient_Vector(0:PROB_DIM-1 ) )











ALLOCATE(Source_Term_Coefficients(  1:NUM_P_QUAD_POINTS,        &
                                    1:NUM_T_QUAD_POINTS,        &
                                    1:NUM_R_QUAD_POINTS,        &
                                    0:NUM_P_ELEMENTS-1,         &
                                    0:NUM_T_ELEMENTS-1,         &
                                    0:NUM_R_ELEMENTS-1     )    )





ALLOCATE( LOCAL_NODE_LOCATIONS(0:DEGREE) )

ALLOCATE( M_VALUES(0:L_LIMIT) )


ALLOCATE( Ylm_Values(   0:LM_Length,                                            &
                        1:NUM_TP_QUAD_POINTS,                                   &
                        0:NUM_T_ELEMS_PER_BLOCK-1, 0:NUM_P_ELEMS_PER_BLOCK-1)   )

ALLOCATE( Ylm_CC_Values(    0:LM_Length,                                            &
                            1:NUM_TP_QUAD_POINTS,                                   &
                            0:NUM_T_ELEMS_PER_BLOCK-1, 0:NUM_P_ELEMS_PER_BLOCK-1)   )


ALLOCATE( Lagrange_Poly_Table(0:DEGREE, 1:NUM_R_QUAD_POINTS, 0:2)   )
ALLOCATE( LPT_LPT( 1:NUM_R_QUAD_POINTS,0:DEGREE,0:DEGREE,0:1,0:2)       )

ALLOCATE(INT_R_LOCATIONS(1:NUM_R_QUAD_POINTS), INT_R_WEIGHTS(1:NUM_R_QUAD_POINTS))
ALLOCATE(INT_T_LOCATIONS(1:NUM_T_QUAD_POINTS), INT_T_WEIGHTS(1:NUM_T_QUAD_POINTS))
ALLOCATE(INT_P_LOCATIONS(1:NUM_P_QUAD_POINTS), INT_P_WEIGHTS(1:NUM_P_QUAD_POINTS))





STF_NNZ = NUM_R_ELEMENTS*(DEGREE + 1)*(DEGREE + 1) - NUM_R_ELEMENTS + 1

!!! Allocate CCS Matrix Arrays !!!
ALLOCATE(STF_ELEM_VAL(0:STF_NNZ-1,0:L_LIMIT), STF_ROW_IND(0:STF_NNZ-1), STF_COL_PTR(0:NUM_R_NODES))







ALLOCATE( ITER_TIME_TABLE(1:Num_Timer_Calls) )
ITER_TIME_TABLE = 0.0_idp


END SUBROUTINE Allocate_Poseidon_Newtonian_Variables










!+202+###########################################################################!
!                                                                                !
!                    Deallocate_Poseidon_Newtonian_Variables                     !
!                                                                                !
!################################################################################!
SUBROUTINE Deallocate_Poseidon_Newtonian_Variables()


DEALLOCATE(rlocs)
DEALLOCATE(tlocs)
DEALLOCATE(plocs)






DEALLOCATE( Lagrange_Poly_Table )


DEALLOCATE(INT_R_LOCATIONS, INT_R_WEIGHTS)
DEALLOCATE(INT_T_LOCATIONS, INT_T_WEIGHTS)
DEALLOCATE(INT_P_LOCATIONS, INT_P_WEIGHTS)

DEALLOCATE( Source_Term_Coefficients )



DEALLOCATE( Coefficient_Vector )


DEALLOCATE(RHS_Vector)



END SUBROUTINE Deallocate_Poseidon_Newtonian_Variables






END MODULE Allocate_Variables_Module
