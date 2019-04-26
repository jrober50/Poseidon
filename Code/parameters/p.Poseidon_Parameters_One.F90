   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_Parameters                                                          !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains the top level subroutines needed to inialize, run, and close       !##!
!##!        Poseidon.                                                               !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!


USE Poseidon_Constants_Module, &
            ONLY :  idp, pi, fdp


IMPLICIT NONE




INTEGER, PARAMETER           ::  DOMAIN_DIM    = 3
INTEGER, PARAMETER           ::  DEGREE        = 1
INTEGER, PARAMETER           ::  L_LIMIT       = 0

INTEGER, PARAMETER           ::  nPROCS_POSEIDON     = 1

INTEGER, PARAMETER           ::  NUM_CFA_VARS  = 5  ! 2+DOMAIN_DIM

INTEGER, PARAMETER           ::  DATA_DIST_MODE          = 4
INTEGER, PARAMETER           ::  SOL_DIST_SCHEME         = 2
INTEGER, PARAMETER           ::  STF_MAPPING_FLAG        = 2


INTEGER, PARAMETER           ::  NUM_SHELLS              = 1
INTEGER, PARAMETER           ::  NUM_SUBSHELLS_PER_SHELL = 1

INTEGER, PARAMETER           ::  NUM_BLOCKS_PER_SHELL    = 1 ! \
INTEGER, PARAMETER           ::  NUM_BLOCK_THETA_ROWS    = 1 ! - NUM_BLOCKS_PER_SHELL = NUM_BLOCK_THETA_ROWS
INTEGER, PARAMETER           ::  NUM_BLOCK_PHI_COLUMNS   = 1 ! /                      * NUM_BLOCK_PHI_COLUMNS




INTEGER, PARAMETER           ::  RAD_ELEMS_PER_SHELL     = 722 ! NUM_SHELS*RAD_ELEMS_PER_SHELL = CHIMERA_R_ELEM
INTEGER, PARAMETER           ::  NUM_R_ELEMS_PER_BLOCK   = RAD_ELEMS_PER_SHELL
INTEGER, PARAMETER           ::  NUM_T_ELEMS_PER_BLOCK   = 1  ! NUM_BLOCK_THETA_ROW*NUM_T_ELEMS_PER_BLOCK = CHIMERA_T_ELEM
INTEGER, PARAMETER           ::  NUM_P_ELEMS_PER_BLOCK   = 1  ! NUM_BLOCK_PHI_COLUMNS*NUM_P_ELEMS_PER_BLOCK = CHIMERA_P_ELEM

INTEGER, PARAMETER           ::  NUM_R_ELEMS_PER_SUBSHELL   = 722 ! NUM_R_ELEMS_PER_SUBSHELL = NUM_R_ELEMS_PER_BLOCK/NUM_SUBSHELLS_PER_SHELL


INTEGER, PARAMETER           ::  NUM_R_QUAD_POINTS = 3
INTEGER, PARAMETER           ::  NUM_T_QUAD_POINTS = 1
INTEGER, PARAMETER           ::  NUM_P_QUAD_POINTS = 1






END MODULE Poseidon_Parameters
