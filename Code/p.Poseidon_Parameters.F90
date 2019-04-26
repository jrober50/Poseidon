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




INTEGER           ::  DOMAIN_DIM
INTEGER           ::  DEGREE
INTEGER           ::  L_LIMIT

INTEGER           ::  nPROCS_POSEIDON

INTEGER, PARAMETER           ::  NUM_CFA_VARS = 5

INTEGER, PARAMETER           ::  DATA_DIST_MODE = 4
INTEGER, PARAMETER           ::  SOL_DIST_SCHEME = 2
INTEGER, PARAMETER           ::  STF_MAPPING_FLAG = 2

INTEGER           ::  NUM_SHELLS
INTEGER           ::  NUM_SUBSHELLS_PER_SHELL


INTEGER           ::  NUM_BLOCKS_PER_SHELL
INTEGER           ::  NUM_BLOCK_THETA_ROWS
INTEGER           ::  NUM_BLOCK_PHI_COLUMNS


INTEGER           ::  NUM_R_ELEMS_PER_SHELL     ! NUM_SHELS*RAD_ELEMS_PER_SHELL = CHIMERA_R_ELEM
INTEGER           ::  NUM_R_ELEMS_PER_SUBSHELL  ! NUM_R_ELEMS_PER_SUBSHELL = NUM_R_ELEMS_PER_BLOCK/NUM_SUBSHELLS_PER_SHELL


INTEGER           ::  NUM_R_ELEMS_PER_BLOCK
INTEGER           ::  NUM_T_ELEMS_PER_BLOCK     ! NUM_BLOCK_THETA_ROW*NUM_T_ELEMS_PER_BLOCK = CHIMERA_T_ELEM
INTEGER           ::  NUM_P_ELEMS_PER_BLOCK     ! NUM_BLOCK_PHI_COLUMNS*NUM_P_ELEMS_PER_BLOCK = CHIMERA_P_ELEM


INTEGER           ::  NUM_R_QUAD_POINTS
INTEGER           ::  NUM_T_QUAD_POINTS
INTEGER           ::  NUM_P_QUAD_POINTS

INTEGER           ::  NUM_QUAD_DOF

INTEGER           ::  R_COARSEN_FACTOR
INTEGER           ::  T_COARSEN_FACTOR
INTEGER           ::  P_COARSEN_FACTOR


INTEGER           ::  NUM_BLOCKS            
INTEGER           ::  NUM_SUBSHELLS

INTEGER           ::  MAX_ITERATIONS
REAL(KIND =idp)   ::  CONVERGENCE_CRITERIA











END MODULE Poseidon_Parameters
