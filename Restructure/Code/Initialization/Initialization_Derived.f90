   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Initialization_Derived                                                       !##!
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

USE Poseidon_Parameters, &
                    ONLY :  Domain_Dim,             &
                            Degree,                 &
                            L_Limit,                &
                            Num_CFA_Vars


USE Variables_Mesh, &
                    ONLY :  Num_R_Elements


USE Variables_MPI, &
                    ONLY :  Num_R_Elems_Per_Block,  &
                            Num_R_Elems_Per_SubShell

USE Variables_Derived, &
                    ONLY :  Num_R_Nodes,            &
                            Block_Num_R_Nodes,      &
                            SubShell_Num_R_Nodes,   &
                            Var_Dim,                &
                            Elem_Var_Dim,           &
                            Block_Var_Dim,          &
                            SubShell_Var_Dim,       &
                            LM_Length,              &
                            ULM_Length,             &
                            Prob_Dim,               &
                            Elem_Prob_Dim,          &
                            Elem_Prob_Dim_Sqr,      &
                            Block_Prob_Dim,         &
                            SubShell_Prob_Dim,      &
                            Num_Off_Diagonals,      &
                            Beta_Prob_Dim,          &
                            Beta_Elem_Prob_Dim
                            
USE Variables_Functions, &
                    ONLY :  Matrix_Location,            &
                            LM_Location

USE Functions_Mapping, &
                    ONLY :  CFA_1D_Matrix_Map,      &
                            CFA_2D_Matrix_Map,      &
                            CFA_3D_Matrix_Map,      &
                            CFA_1D_LM_Map,          &
                            CFA_2D_LM_Map,          &
                            CFA_3D_LM_Map


IMPLICIT NONE

CONTAINS
 !+101+################################################################################!
!                                                                                       !
!       Initialize_Mesh                                                                 !
!                                                                                       !
 !#####################################################################################!
SUBROUTINE Initialize_Derived()



IF ( DOMAIN_DIM == 1 ) THEN

    LM_LENGTH = 1
    Matrix_Location => CFA_1D_Matrix_Map
    LM_Location => CFA_1D_LM_Map

ELSE IF ( DOMAIN_DIM == 2 ) THEN

    LM_LENGTH = L_LIMIT + 1
    Matrix_Location => CFA_2D_Matrix_Map
    LM_Location => CFA_2D_LM_Map

ELSE IF ( DOMAIN_DIM == 3 ) THEN

    LM_LENGTH = (L_LIMIT + 1)*(L_LIMIT + 1)
    Matrix_Location => CFA_3D_Matrix_Map
    LM_Location => CFA_3D_LM_Map

END IF

NUM_R_NODES             = DEGREE*NUM_R_ELEMENTS + 1
BLOCK_NUM_R_NODES       = DEGREE*NUM_R_ELEMS_PER_BLOCK + 1
SUBSHELL_NUM_R_NODES    = DEGREE*NUM_R_ELEMS_PER_SUBSHELL + 1


!
!   VAR_DIM - Length of vector required to hold coefficients for 1 variable.
!
VAR_DIM             = LM_LENGTH*NUM_R_NODES
ELEM_VAR_DIM        = LM_LENGTH*(DEGREE + 1)
BLOCK_VAR_DIM       = LM_LENGTH*BLOCK_NUM_R_NODES
SUBSHELL_VAR_DIM    = LM_LENGTH*SUBSHELL_NUM_R_NODES


!
!   PROB_DIM - Length of vector required to hold coefficients for all variables
!
ULM_LENGTH          = NUM_CFA_VARS*LM_LENGTH
PROB_DIM            = NUM_CFA_VARS*VAR_DIM
ELEM_PROB_DIM       = NUM_CFA_VARS*ELEM_VAR_DIM
ELEM_PROB_DIM_SQR   = ELEM_PROB_DIM*ELEM_PROB_DIM
BLOCK_PROB_DIM      = NUM_CFA_VARS*BLOCK_VAR_DIM
SUBSHELL_PROB_DIM   = NUM_CFA_VARS*SUBSHELL_VAR_DIM



Beta_Prob_Dim       = 3*Var_Dim
Beta_Elem_Prob_Dim  = 3*Elem_Var_Dim


NUM_OFF_DIAGONALS = ULM_LENGTH*(DEGREE + 1) - 1



END SUBROUTINE Initialize_Derived





END MODULE Initialization_Derived


