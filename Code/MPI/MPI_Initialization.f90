   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Initialization_MPI                                                           !##!
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
                        Verbose_Flag,           &
                        Num_CFA_Vars

USE Variables_Mesh, &
                ONLY :  Num_R_Elements,             &
                        Num_T_Elements,             &
                        Num_P_Elements
USE Variables_MPI, &
                ONLY :  Num_R_Elems_Per_Shell,      &
                        Num_R_Elems_Per_SubShell,   &
                        Num_R_Elems_Per_Block,      &
                        Num_T_Elems_Per_Block,      &
                        Num_P_Elems_Per_Block,      &
                        Ratio_BNDLperBLCK,          &
                        Ratio_T_BNDLperBLCK,        &
                        Ratio_P_BNDLperBLCK

USE Variables_Derived, &
                ONLY :  Block_Num_R_Nodes,      &
                        Block_Var_Dim,          &
                        Block_Prob_Dim,         &
                        SubShell_Num_R_Nodes,   &
                        SubShell_Var_Dim,       &
                        SubShell_Prob_Dim,      &
                        LM_Length

USE Functions_MPI, &
                ONLY : Create_Poseidon_Communicators

IMPLICIT NONE

CONTAINS

 !+101+################################################################################!
!                                                                                       !
!       Initialize_Quad                                                                 !
!                                                                                       !
 !#####################################################################################!
SUBROUTINE Initialize_MPI()


IF ( Verbose_Flag ) THEN
    PRINT*,"-Initializing MPI Variables. "
END IF


Num_R_Elems_Per_Shell = Num_R_Elements
Num_R_Elems_Per_SubShell = Num_R_Elements
Num_R_Elems_Per_Block = Num_R_Elements
Num_T_Elems_Per_Block = Num_T_Elements
Num_P_Elems_Per_Block = Num_P_Elements

BLOCK_NUM_R_NODES   = DEGREE*NUM_R_ELEMS_PER_BLOCK + 1
BLOCK_VAR_DIM       = LM_LENGTH*BLOCK_NUM_R_NODES
BLOCK_PROB_DIM      = NUM_CFA_VARS*BLOCK_VAR_DIM


SUBSHELL_NUM_R_NODES    = DEGREE*NUM_R_ELEMS_PER_SUBSHELL + 1
SUBSHELL_VAR_DIM        = LM_LENGTH*SUBSHELL_NUM_R_NODES
SUBSHELL_PROB_DIM       = NUM_CFA_VARS*SUBSHELL_VAR_DIM


Ratio_T_BNDLperBLCK = 1
Ratio_P_BNDLperBLCK = 1
Ratio_BNDLperBLCK = Ratio_T_BNDLperBLCK * Ratio_P_BNDLperBLCK

CALL Create_Poseidon_Communicators()


END SUBROUTINE Initialize_MPI





END MODULE Initialization_MPI

