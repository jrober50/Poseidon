   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_Initialization_Module                                                   !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains the top level subroutines needed to inialize, run, and close       !##!
!##!        Poseidon.                                                               !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!    +101+   Poseidon_Initialize                                                 !##!
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
USE Poseidon_Constants_Module, &
            ONLY :  idp, pi, fdp

USE Units_Module, &
            ONLY :  Grav_Constant_G,            &
                    C_Square

USE Poseidon_Parameters, &
            ONLY :  DOMAIN_DIM,                 &
                    DEGREE,                     &
                    L_LIMIT,                    &
                    DATA_DIST_MODE,             &
                    NUM_CFA_VARS,               &
                    NUM_R_ELEMS_PER_SHELL,      &
                    NUM_R_ELEMS_PER_SUBSHELL,   &
                    NUM_SHELLS,                 &
                    NUM_SUBSHELLS,              &
                    NUM_SUBSHELLS_PER_SHELL,    &
                    NUM_BLOCKS,                 &
                    NUM_BLOCKS_PER_SHELL,       &
                    NUM_BLOCK_THETA_ROWS,       &
                    NUM_BLOCK_PHI_COLUMNS,      &
                    nPROCS_POSEIDON,            &
                    STF_MAPPING_FLAG,           &
                    NUM_R_ELEMS_PER_BLOCK,      &
                    NUM_T_ELEMS_PER_BLOCK,      &
                    NUM_P_ELEMS_PER_BLOCK,      &
                    Ratio_T_BNDLperBLCK,        &
                    Ratio_P_BNDLperBLCK,        &
                    Ratio_BNDLperBLCK,          &
                    R_COARSEN_FACTOR,           &
                    T_COARSEN_FACTOR,           &
                    P_COARSEN_FACTOR,           &
                    NUM_R_QUAD_POINTS,          &
                    NUM_T_QUAD_POINTS,          &
                    NUM_P_QUAD_POINTS,          &
                    CUR_ITERATION,              &
                    MAX_ITERATIONS,             &
                    CONVERGENCE_CRITERIA,       &
                    CONVERGENCE_FLAG,           &
                    ITER_REPORT_NUM_SAMPLES,    &
                    OUTPUT_SETUP_TABLE_FLAG,    &
                    WRITE_REPORT_FLAG,          &
                    RUN_REPORT_FILE_ID,         &
                    POSEIDON_INITIALIZED_FLAG



USE Poseidon_Variables_Module, &
            ONLY :  R_INNER, R_OUTER,                               &
                    NUM_R_ELEMENTS, NUM_T_ELEMENTS, NUM_P_ELEMENTS, &
                    NUM_LOC_R_ELEMENTS,         &
                    NUM_LOC_T_ELEMENTS,         &
                    NUM_LOC_P_ELEMENTS,         &
                    NUM_TP_QUAD_POINTS,                             &
                    NUM_R_NODES,                                    &
                    BLOCK_NUM_R_NODES,                              &
                    SUBSHELL_NUM_R_NODES,                           &
                    rlocs,                                          &
                    tlocs,                                          &
                    plocs,                                          &
                    drlocs,                                         &
                    dtlocs,                                         &
                    dplocs,                                         &
                    RHS_Vector,                                     &
                    Coefficient_Vector,                             &
                    Source_Term_Coefficients,                       &
                    Test_Space_Allocated_Flag,                      &
                    Stiffness_Matrix_Initialized_Flag,              &
                    FirstCall_Flag,                                 &
                    RADIAL_MESH_SET_FLAG,                           &
                    THETA_MESH_SET_FLAG,                            &
                    PHI_MESH_SET_FLAG,                              &
                    INT_R_LOCATIONS,                                &
                    INT_R_WEIGHTS,                                  &
                    INT_T_LOCATIONS,                                &
                    INT_T_WEIGHTS,                                  &
                    INT_P_LOCATIONS,                                &
                    INT_P_WEIGHTS,                                  &
                    INT_TP_WEIGHTS,                                 &
                    LOCAL_NODE_LOCATIONS,                           &
                    ierr,                                           &
                    myID_Poseidon,                                  &
                    PHYSICS_TYPE,                                   &
                    VAR_DIM,                                        &
                    ELEM_VAR_DIM,                                   &
                    BLOCK_VAR_DIM,                                  &
                    SUBSHELL_VAR_DIM,                               &
                    PROB_DIM,                                       &
                    ELEM_PROB_DIM,                                  &
                    ELEM_PROB_DIM_SQR,                              &
                    BLOCK_PROB_DIM,                                 &
                    SUBSHELL_PROB_DIM,                              &
                    INNER_CFA_BC_VALUES,                            &
                    OUTER_CFA_BC_VALUES,                            &
                    INNER_CFA_BC_TYPE,                              &
                    OUTER_CFA_BC_TYPE,                              &
                    NUM_OFF_DIAGONALS,                              &
                    LM_LENGTH,                                      &
                    ULM_LENGTH,                                     &
                    M_VALUES,                                       &
                    Matrix_Location,                                &
                    LM_Location,                                    &
                    POSEIDON_COMM_WORLD,                            &
                    RUN_TIME_TABLE



USE Poseidon_Quadrature_Module, &
            ONLY :  Initialize_LG_Quadrature


USE Poseidon_Mesh_Module, &
            ONLY :  Generate_Defined_Coarse_Mesh


USE Poseidon_Tables_Module, &
            ONLY :  Initialize_Ylm_Tables,                          &
                    Initialize_Lagrange_Poly_Tables

USE Poseidon_Allocation_Module, &
            ONLY :  Allocate_Poseidon_CFA_Variables,                &
                    Deallocate_Poseidon_CFA_Variables


USE Poseidon_Calculate_Results_Module, &
            ONLY :  Calc_3D_Values_At_Location

USE Poseidon_MPI_Module, &
            ONLY :  CREATE_POSEIDON_COMMUNICATORS

USE Poseidon_Parameter_Read_Module, &
            ONLY :  UNPACK_POSEIDON_PARAMETERS

USE Poseidon_Mapping_Functions_Module,  &
            ONLY :  CFA_3D_Matrix_Map,                              &
                    CFA_2D_Matrix_Map,                              &
                    CFA_1D_Matrix_Map,                              &
                    CFA_1D_LM_Map,                                  &
                    CFA_2D_LM_Map,                                  &
                    CFA_3D_LM_Map


USE mpi




IMPLICIT NONE




                    !*F&S*==========================================!
                    !                                               !
                    !           Functions & Subroutines             !
                    !                                               !
                    !===============================================!
CONTAINS










 !+101+####################################################################################!
!                                                                                           !
!       Poseidon_Initialize_From_File                                                       !
!                                                                                           !
!===========================================================================================!
!                                                                                           !
!   Sets code parameters, allocates space, and initializes functions and variables needed   !
!   to run Poseidon.                                                                        !
!                                                                                           !
!-------------------------------------------------------------------------------------------!
!                                                                                           !
!       Input Variables     :                                                               !
!                                                                                           !
!       FEM_Degree_Input        -       Integer, Order of Finite Element Method solver to   !
!                                               be performed.                               !
!                                                                                           !
!       L_Limit_Input           -       Integer, Limit on the l value of the Spherical      !
!                                               Harmonic spectral decomposition.            !
!                                                                                           !
!       Inner_Radius            -       Real, Inner radius of the computational domain.     !
!                                                                                           !
!       Outer_Radius            -       Real, Outer radius of the computational domain.     !
!                                                                                           !
!       R_Elements_Input        -       Integer, Number of radial elements.                 !
!                                                                                           !
!       T_Elements_Input        -       Integer, Number of theta elements.                  !
!                                                                                           !
!       P_Elements_Input        -       Integer, Number of phi elements.                    !
!                                                                                           !
!       Input_Delta_R_Vector    -       Optional Real Vector, Dimension(1:R_Elements_Input) !
!                                       Each value corresponds to the radial length of the  !
!                                       each radial shell of elements starting from         !
!                                       Inner_Radius.                                       !
!                                                                                           !
!       Input_Delta_T_Vector    -       Optional Real Vector, Dimension(1:T_Elements_Input) !
!                                       Each value correspons to the angular width of each  !
!                                       theta wedge of elements starting from 0.            !
!                                                                                           !
!       Input_Delta_P_Vector    -       Optional Real Vector, Dimension(1:P_Elements_Input) !
!                                       Each value correspons to the angular width of each  !
!                                       phi wedge of elements starting from 0.              !
!                                                                                           !
 !#########################################################################################!
SUBROUTINE Poseidon_Initialize_From_File( mode,                                                         &
                                Inner_Radius, Outer_Radius,                                             &
                                R_Elements_Input, T_Elements_Input, P_Elements_Input,                   &
                                Local_R_Elements_Input, Local_T_Elements_Input, Local_P_Elements_Input, &
                                Input_Delta_R_Vector, Input_Delta_T_Vector, Input_Delta_P_Vector        )



                                         !                          !
                                        !!      Input Variables     !!
                                         !

INTEGER, INTENT(IN)                                                             ::  mode,                   &
                                                                                    R_Elements_Input,       &
                                                                                    T_Elements_Input,       &
                                                                                    P_Elements_Input,       &
                                                                                    Local_R_Elements_Input, &
                                                                                    Local_T_Elements_Input, &
                                                                                    Local_P_Elements_Input

REAL(KIND = idp), INTENT(IN)                                                    ::  Inner_Radius,           &
                                                                                    Outer_Radius


REAL(KIND = idp), DIMENSION(1:R_Elements_Input),    OPTIONAL,   INTENT(IN)      ::  Input_Delta_R_Vector
REAL(KIND = idp), DIMENSION(1:T_Elements_Input),    OPTIONAL,   INTENT(IN)      ::  Input_Delta_T_Vector
REAL(KIND = idp), DIMENSION(1:P_Elements_Input),    OPTIONAL,   INTENT(IN)      ::  Input_Delta_P_Vector





                                         !                              !
                                        !!     Subroutine Variables     !!
                                         !                              !

REAL(KIND = idp), DIMENSION(1:R_Elements_Input)                                 ::  Delta_R_Vector
REAL(KIND = idp), DIMENSION(1:T_Elements_Input)                                 ::  Delta_T_Vector
REAL(KIND = idp), DIMENSION(1:P_Elements_Input)                                 ::  Delta_P_Vector


INTEGER                                                                         ::  td, pd
INTEGER                                                                         ::  nPROCS, tmpID
INTEGER                                                                         ::  l



REAL(KIND = idp)                                                                :: timea, timeb, timec



CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nPROCS,ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, tmpID, ierr)


 !                                      !
!!  Load Poseidon Parameters From File  !!
 !                                      !
CALL UNPACK_POSEIDON_PARAMETERS()



 !                                          !
!!  Set Global Variables to Input Values    !!
 !                                          !
R_INNER = Inner_Radius
R_OUTER = Outer_Radius


NUM_R_ELEMENTS = R_Elements_Input/R_Coarsen_Factor
NUM_T_ELEMENTS = T_Elements_Input/T_Coarsen_Factor
NUM_P_ELEMENTS = P_Elements_Input/P_Coarsen_Factor

NUM_LOC_R_ELEMENTS = Local_R_Elements_Input
NUM_LOC_T_ELEMENTS = Local_T_Elements_Input
NUM_LOC_P_ELEMENTS = Local_P_Elements_Input

!
!   *_NUM_R_NODES - Number of total radial nodes in total, per block, per subshell
!
NUM_R_NODES             = DEGREE*NUM_R_ELEMENTS + 1
BLOCK_NUM_R_NODES       = DEGREE*NUM_R_ELEMS_PER_BLOCK + 1
SUBSHELL_NUM_R_NODES    = DEGREE*NUM_R_ELEMS_PER_SUBSHELL + 1


!
!   NUM_TP_QUAD_POINTS = number of angular quadrature points per radial point
!
NUM_TP_QUAD_POINTS = NUM_T_QUAD_POINTS*NUM_P_QUAD_POINTS


CALL CHECK_SETUP(R_Elements_Input, T_Elements_Input, P_Elements_Input, TmpID )







!
!   Associate the Correct Map Functions, and Set Spherical Harmonic Length
!
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
IF ( PHYSICS_TYPE == "NEW" ) THEN

    ULM_LENGTH = LM_LENGTH
    PROB_DIM = VAR_DIM
    ELEM_PROB_DIM = ELEM_VAR_DIM
    ELEM_PROB_DIM_SQR = ELEM_VAR_DIM*ELEM_VAR_DIM
    BLOCK_PROB_DIM = BLOCK_VAR_DIM
    SUBSHELL_PROB_DIM = SUBSHELL_VAR_DIM

ELSE IF ( PHYSICS_TYPE == "CFA" ) THEN

    ULM_LENGTH          = NUM_CFA_VARS*LM_LENGTH
    PROB_DIM            = NUM_CFA_VARS*VAR_DIM
    ELEM_PROB_DIM       = NUM_CFA_VARS*ELEM_VAR_DIM
    ELEM_PROB_DIM_SQR   = ELEM_PROB_DIM*ELEM_PROB_DIM
    BLOCK_PROB_DIM      = NUM_CFA_VARS*BLOCK_VAR_DIM
    SUBSHELL_PROB_DIM   = NUM_CFA_VARS*SUBSHELL_VAR_DIM

END IF


NUM_OFF_DIAGONALS = ULM_LENGTH*(DEGREE + 1) - 1




                                 !                                      !
                                !!      Allocate Space for Poseidon     !!
                                 !                                      !
CALL Allocate_Poseidon_CFA_Variables()




                                 !                                          !
                                !!      Initialize Reusable Quadratures     !!
                                 !                                          !

CALL Initialize_LG_Quadrature(NUM_R_QUAD_POINTS, INT_R_LOCATIONS, INT_R_WEIGHTS)
CALL Initialize_LG_Quadrature(NUM_T_QUAD_POINTS, INT_T_LOCATIONS, INT_T_WEIGHTS)
CALL Initialize_LG_Quadrature(NUM_P_QUAD_POINTS, INT_P_LOCATIONS, INT_P_WEIGHTS)

IF ( NUM_T_QUAD_POINTS == 1 ) THEN
    ! 1D Cheat !
    INT_T_WEIGHTS = 4.0_idp/pi
END IF


DO td = 1,NUM_T_QUAD_POINTS
    DO pd = 1,NUM_P_QUAD_POINTS
        INT_TP_WEIGHTS( (td-1)*NUM_P_QUAD_POINTS+pd ) = INT_T_WEIGHTS(td)*INT_P_WEIGHTS(pd)
    END DO
END DO



                                 !                                          !
                                !!   Initialize Lagrange Polynomial Tables  !!
                                 !                                          !
CALL Initialize_Lagrange_Poly_Tables()








                                 !                                          !
                                !!      Set Initial Mesh (Optional)         !!
                                 !                                          !

CALL Generate_Defined_Coarse_Mesh(R_Elements_Input, NUM_R_ELEMENTS, R_Coarsen_Factor,     &
                                  R_INNER, Input_Delta_R_Vector, rlocs, drlocs)
RADIAL_MESH_SET_FLAG = .TRUE.


CALL Generate_Defined_Coarse_Mesh(T_Elements_Input, NUM_T_ELEMENTS, T_Coarsen_Factor,     &
                                  0.0_idp, Input_Delta_T_Vector, tlocs, dtlocs)
THETA_MESH_SET_FLAG = .TRUE.


CALL Generate_Defined_Coarse_Mesh(P_Elements_Input, NUM_P_ELEMENTS, P_Coarsen_Factor,     &
                                  0.0_idp, Input_Delta_P_Vector, plocs, dplocs)
PHI_MESH_SET_FLAG = .TRUE.





IF ( DOMAIN_DIM == 1 ) THEN
     M_VALUES = 0
ELSE IF ( DOMAIN_DIM == 2 ) THEN
     M_VALUES = 0
ELSE IF ( DOMAIN_DIM == 3 ) THEN
     M_VALUES = (/(l,l=0,L_LIMIT,1)/)
END IF




Ratio_T_BNDLperBLCK = NUM_T_ELEMENTS/( NUM_BLOCK_THETA_ROWS*NUM_LOC_T_ELEMENTS )
Ratio_P_BNDLperBLCK = NUM_P_ELEMENTS/( NUM_BLOCK_PHI_COLUMNS*NUM_LOC_P_ELEMENTS )
Ratio_BNDLperBLCK = Ratio_T_BNDLperBLCK * Ratio_P_BNDLperBLCK

CALL CREATE_POSEIDON_COMMUNICATORS( DATA_DIST_MODE )

CALL Initialize_Ylm_Tables()



CALL OUTPUT_SETUP_TABLE( nPROCS, R_Elements_Input, T_Elements_Input, P_Elements_Input )




END SUBROUTINE Poseidon_Initialize_From_File




!+201+######################################################################################!
!                                                                                           !
!       Check_Setup                                                                         !
!                                                                                           !
!###########################################################################################!
SUBROUTINE CHECK_SETUP( R_Elements_Input, T_Elements_Input, P_Elements_Input, TmpID )

INTEGER, INTENT(IN)                                         ::  R_Elements_Input,       &
                                                                T_Elements_Input,       &
                                                                P_Elements_Input,       &
                                                                TmpID





IF (NUM_R_ELEMS_PER_SHELL*NUM_SHELLS .NE. NUM_R_ELEMENTS ) THEN

    IF ( tmpID == 0) THEN
          PRINT*,"Domain not properly broken into shells.  Number of radial elements does not match."
          PRINT*, " # of Elements Per Shell : ",NUM_R_ELEMS_PER_SHELL ,             &
              " # Shells : ",NUM_SHELLS ,                                         &
              " # Global Radial Elements",NUM_R_ELEMENTS
          PRINT*,"Checked in z.Poseidon_Initialization_Module.F90, line 840 "
    END IF
    CALL MPI_FINALIZE(ierr)
    STOP


END IF


IF ( NUM_BLOCKS .NE. NUM_BLOCKS_PER_SHELL*NUM_SHELLS ) THEN


    IF ( tmpID == 0) THEN
       PRINT*,"Domain not properly decomposed.  Expected number of global blocks does not match decomposition."
       PRINT*, " # Total Blocks : ",NUM_BLOCKS,                        &
            " # of Blocks Per Shell : ",NUM_BLOCKS_PER_SHELL,       &
            " # of Shells",NUM_SHELLS
       PRINT*,"Checked in z.Poseidon_Initialization_Module.F90, line 857 "
    END IF
    CALL MPI_FINALIZE(ierr)
    STOP

END IF


IF ( NUM_BLOCKS_PER_SHELL .NE. NUM_BLOCK_THETA_ROWS*NUM_BLOCK_PHI_COLUMNS ) THEN


    IF ( tmpID == 0) THEN
       PRINT*,"Domain not properly decomposed.  Expected number of shell-wise blocks does not match decomposition."
       PRINT*, " # of Blocks Per Shell : ",NUM_BLOCKS_PER_SHELL,       &
            " # of Theta rows",NUM_BLOCK_THETA_ROWS,                &
            " # of Phi columns",NUM_BLOCK_PHI_COLUMNS
       PRINT*,"Checked in z.Poseidon_Initialization_Module.F90, line 873 "
    END IF

    CALL MPI_FINALIZE(ierr)
    STOP

END IF




IF ( NUM_BLOCKS .NE. nPROCS_POSEIDON ) THEN


    IF ( tmpID == 0) THEN
       PRINT*,"Poseidon's number of processes does not match the number of blocks"
       PRINT*, " # of Blocks : ",NUM_BLOCKS,                        &
       " # of Poseidon Processes expected : ",nPROCs_POSEIDON
       PRINT*,"Checked in z.Poseidon_Initialization_Module.F90, line 891 "
    END IF
    CALL MPI_FINALIZE(ierr)
    STOP

END IF




IF ( NUM_R_ELEMS_PER_BLOCK.NE. NUM_R_ELEMS_PER_SHELL ) THEN

    IF ( tmpID == 0) THEN
       PRINT*,"The number of elements per block does not match."
       PRINT*, " # of Elements per Block : ",NUM_R_ELEMS_PER_BLOCK
       PRINT*, " # of Radial Elements Per Shell : ",NUM_R_ELEMS_PER_SHELL
       PRINT*,"Checked in z.Poseidon_Initialization_Module.F90, line 907 "
    END IF
    CALL MPI_FINALIZE(ierr)
    STOP

END IF






IF ( R_COARSEN_FACTOR*NUM_R_ELEMS_PER_SHELL*NUM_SHELLS .NE. R_Elements_Input) THEN

     IF (tmpID == 0) THEN

        PRINT*,"Poseidon's Radial Decomposition is incompatible with Input"
        PRINT*,"Expected Global Radial Elements",R_COARSEN_FACTOR*NUM_R_ELEMS_PER_SHELL*NUM_SHELLS
        PRINT*,"Input Global Radial Elements",R_Elements_Input
        PRINT*,"Checked in Z.Poseidon_Initialization_Module.F90, Poseidon_Initialize"

     END IF
     CALL MPI_FINALIZE(ierr)
     STOP

END IF

IF ( T_COARSEN_FACTOR*NUM_T_ELEMS_PER_BLOCK*NUM_BLOCK_THETA_ROWS .NE. T_Elements_Input) THEN

     IF (tmpID == 0) THEN

        PRINT*,"Poseidon's Phi Decomposition is incompatible with Input"
        PRINT*,"Expected Global Theta Elements",T_COARSEN_FACTOR*NUM_T_ELEMS_PER_BLOCK*NUM_BLOCK_THETA_ROWS
        PRINT*,"Input Global Theta Elements",T_Elements_Input
        PRINT*,"Checked in Z.Poseidon_Initialization_Module.F90, Poseidon_Initialize"

     END IF
     CALL MPI_FINALIZE(ierr)
     STOP

END IF


IF ( P_COARSEN_FACTOR*NUM_P_ELEMS_PER_BLOCK*NUM_BLOCK_PHI_COLUMNS .NE. P_Elements_Input) THEN

     IF (tmpID == 0) THEN

        PRINT*,"Poseidon's Phi Decomposition is incompatible with Input"
        PRINT*,"Expected Global Phi Elements",P_COARSEN_FACTOR*NUM_P_ELEMS_PER_BLOCK*NUM_BLOCK_PHI_COLUMNS
        PRINT*,"Input Global Phi Elements",P_Elements_Input
        PRINT*,"Checked in Z.Poseidon_Initialization_Module.F90, Poseidon_Initialize"

     END IF
     CALL MPI_FINALIZE(ierr)
     STOP

END IF


END SUBROUTINE CHECK_SETUP












!+301+######################################################################################!
!                                                                                           !
!       OUTPUT_SETUP_TABLE                                                                  !
!                                                                                           !
!###########################################################################################!
SUBROUTINE OUTPUT_SETUP_TABLE( nPROCS, R_ELEMENTS_INPUT, T_ELEMENTS_INPUT, P_ELEMENTS_INPUT )

INTEGER, INTENT(IN)                                                 :: nPROCS
INTEGER, INTENT(IN)                                                 :: R_ELEMENTS_INPUT
INTEGER, INTENT(IN)                                                 :: T_ELEMENTS_INPUT
INTEGER, INTENT(IN)                                                 :: P_ELEMENTS_INPUT
INTEGER                                                                         ::  NUM_JACOBIAN_ELEMENTS
INTEGER                                                                         ::  NONZEROS
INTEGER                                                                         ::  BLOCK_NONZEROS
REAL(KIND = idp)                                                                ::  SPARSITY



1401 FORMAT('------------- POSEIDON PARAMETERS --------------'/)
1402 FORMAT('         ACTIVE DIMENSIONS = ',I12.1)
1403 FORMAT('                    DEGREE = ',I12.1)
1404 FORMAT('                   L_LIMIT = ',I12.1)
1405 FORMAT('              NUM_CFA_VARS = ',I12.1)
1406 FORMAT('                    nPROCS = ',I12.1)
1407 FORMAT('           nPROCS_POSEIDON = ',I12.1)
1408 FORMAT('    SOURCE RADIAL ELEMENTS = ',I12.1)
1409 FORMAT('     SOURCE THETA ELEMENTS = ',I12.1)
1410 FORMAT('       SOURCE PHI ELEMENTS = ',I12.1)
1411 FORMAT('  POSEIDON RADIAL ELEMENTS = ',I12.1)
1412 FORMAT('   POSEIDON THETA ELEMENTS = ',I12.1)
1413 FORMAT('     POSEIDON PHI ELEMENTS = ',I12.1)
1414 FORMAT('          R_COARSEN_FACTOR = ',I12.1)
1415 FORMAT('          T_COARSEN_FACTOR = ',I12.1)
1416 FORMAT('          P_COARSEN_FACTOR = ',I12.1)
1417 FORMAT('               NUM_R_NODES = ',I12.1)
1418 FORMAT('          NUMBER OF SHELLS = ',I12.1)
1419 FORMAT('       NUMBER OF SUBSHELLS = ',I12.1)
1420 FORMAT('       SUBSHELLS PER SHELL = ',I12.1)
1421 FORMAT('      NUM BLOCKS PER SHELL = ',I12.1)
1422 FORMAT('     NUM_R_ELEMS_PER_SHELL = ',I12.1)
1423 FORMAT('     NUM_T_ELEMS_PER_BLOCK = ',I12.1)
1424 FORMAT('     NUM_P_ELEMS_PER_BLOCK = ',I12.1)
1425 FORMAT('                   VAR_DIM = ',I12.1)
1426 FORMAT('             BLOCK_VAR_DIM = ',I12.1)
1427 FORMAT('                  PROB_DIM = ',I12.1)
1428 FORMAT('            BLOCK_PROB_DIM = ',I12.1)
1429 FORMAT('         SUBSHELL_PROB_DIM = ',I12.1)
1430 FORMAT('         NUM_OFF_DIAGONALS = ',I12.1)
1431 FORMAT('     NUM JACOBIAN ELEMENTS = ',I12.1)
1432 FORMAT(' NUM NNZ JACOBIAN ELEMENTS = ',I12.1)
1433 FORMAT('        JACOBIAN SPARSITY  = ',F12.6)
1434 FORMAT('         NUM NNZ PER BLOCK = ',I12.1)
1435 FORMAT('              Inner Radius = ',ES20.12E3)
1436 FORMAT('              Outer Radius = ',ES20.12E3)
1437 FORMAT('      # Radial Quad Points = ',I12.1)
1438 FORMAT('       # Theta Quad Points = ',I12.1)
1439 FORMAT('         # Phi Quad Points = ',I12.1)
1440 FORMAT('        Maximum Iterations = ',I12.1)
1441 FORMAT('      Convergence Criteria = ',ES20.12E3)

IF ( OUTPUT_SETUP_TABLE_FLAG == 1 ) THEN

    NUM_JACOBIAN_ELEMENTS = PROB_DIM*PROB_DIM

    NONZEROS = NUM_R_ELEMENTS*ELEM_PROB_DIM_SQR          &
             - (NUM_R_ELEMENTS-1)*ULM_LENGTH*ULM_LENGTH

    BLOCK_NONZEROS = NUM_R_ELEMS_PER_BLOCK*ELEM_PROB_DIM_SQR          &
                   - (NUM_R_ELEMS_PER_BLOCK-1)*ULM_LENGTH*ULM_LENGTH


    SPARSITY = REAL(NONZEROS,KIND =idp)/REAL(NUM_JACOBIAN_ELEMENTS, KIND = idp)


    IF  ( myID_Poseidon == 0 )  THEN

        WRITE(*,1401)
        WRITE(*,1402)DOMAIN_DIM
        WRITE(*,1403)DEGREE
        WRITE(*,1404)L_LIMIT
        WRITE(*,1405)NUM_CFA_VARS
        WRITE(*,1406)nPROCS
        WRITE(*,1407)nPROCS_POSEIDON
        WRITE(*,1408)R_ELEMENTS_INPUT
        WRITE(*,1409)T_ELEMENTS_INPUT
        WRITE(*,1410)P_ELEMENTS_INPUT
        WRITE(*,1411)NUM_R_ELEMENTS
        WRITE(*,1412)NUM_T_ELEMENTS
        WRITE(*,1413)NUM_P_ELEMENTS
        WRITE(*,1414)R_COARSEN_FACTOR
        WRITE(*,1415)T_COARSEN_FACTOR
        WRITE(*,1416)P_COARSEN_FACTOR
        WRITE(*,1417)NUM_R_NODES
        WRITE(*,1418)NUM_SHELLS
        WRITE(*,1419)NUM_SUBSHELLS
        WRITE(*,1420)NUM_SUBSHELLS_PER_SHELL
        WRITE(*,1421)NUM_BLOCKS_PER_SHELL
        WRITE(*,1422)NUM_R_ELEMS_PER_SHELL
        WRITE(*,1423)NUM_T_ELEMS_PER_BLOCK
        WRITE(*,1424)NUM_P_ELEMS_PER_BLOCK
        WRITE(*,1425)VAR_DIM
        WRITE(*,1426)BLOCK_VAR_DIM
        WRITE(*,1427)PROB_DIM
        WRITE(*,1428)BLOCK_PROB_DIM
        WRITE(*,1429)SUBSHELL_PROB_DIM
        WRITE(*,1430)NUM_OFF_DIAGONALS
        WRITE(*,1431)NUM_JACOBIAN_ELEMENTS
        WRITE(*,1432)NONZEROS
        WRITE(*,1433)SPARSITY
        WRITE(*,1434)BLOCK_NONZEROS
        WRITE(*,1435)R_INNER
        WRITE(*,1436)R_OUTER
        WRITE(*,1437)NUM_R_QUAD_POINTS
        WRITE(*,1438)NUM_T_QUAD_POINTS
        WRITE(*,1439)NUM_P_QUAD_POINTS
        WRITE(*,1440)MAX_ITERATIONS
        WRITE(*,1441)CONVERGENCE_CRITERIA

    END IF
    IF ( myID_Poseidon == 0 ) THEN

        WRITE(RUN_REPORT_FILE_ID,1401)
        WRITE(RUN_REPORT_FILE_ID,1402)DOMAIN_DIM
        WRITE(RUN_REPORT_FILE_ID,1403)DEGREE
        WRITE(RUN_REPORT_FILE_ID,1404)L_LIMIT
        WRITE(RUN_REPORT_FILE_ID,1405)NUM_CFA_VARS
        WRITE(RUN_REPORT_FILE_ID,1406)nPROCS
        WRITE(RUN_REPORT_FILE_ID,1407)nPROCS_POSEIDON
        WRITE(RUN_REPORT_FILE_ID,1408)R_ELEMENTS_INPUT
        WRITE(RUN_REPORT_FILE_ID,1409)T_ELEMENTS_INPUT
        WRITE(RUN_REPORT_FILE_ID,1410)P_ELEMENTS_INPUT
        WRITE(RUN_REPORT_FILE_ID,1411)NUM_R_ELEMENTS
        WRITE(RUN_REPORT_FILE_ID,1412)NUM_T_ELEMENTS
        WRITE(RUN_REPORT_FILE_ID,1413)NUM_P_ELEMENTS
        WRITE(RUN_REPORT_FILE_ID,1414)R_COARSEN_FACTOR
        WRITE(RUN_REPORT_FILE_ID,1415)T_COARSEN_FACTOR
        WRITE(RUN_REPORT_FILE_ID,1416)P_COARSEN_FACTOR
        WRITE(RUN_REPORT_FILE_ID,1417)NUM_R_NODES
        WRITE(RUN_REPORT_FILE_ID,1418)NUM_SHELLS
        WRITE(RUN_REPORT_FILE_ID,1419)NUM_SUBSHELLS
        WRITE(RUN_REPORT_FILE_ID,1420)NUM_SUBSHELLS_PER_SHELL
        WRITE(RUN_REPORT_FILE_ID,1421)NUM_BLOCKS_PER_SHELL
        WRITE(RUN_REPORT_FILE_ID,1422)NUM_R_ELEMS_PER_SHELL
        WRITE(RUN_REPORT_FILE_ID,1423)NUM_T_ELEMS_PER_BLOCK
        WRITE(RUN_REPORT_FILE_ID,1424)NUM_P_ELEMS_PER_BLOCK
        WRITE(RUN_REPORT_FILE_ID,1425)VAR_DIM
        WRITE(RUN_REPORT_FILE_ID,1426)BLOCK_VAR_DIM
        WRITE(RUN_REPORT_FILE_ID,1427)PROB_DIM
        WRITE(RUN_REPORT_FILE_ID,1428)BLOCK_PROB_DIM
        WRITE(RUN_REPORT_FILE_ID,1429)SUBSHELL_PROB_DIM
        WRITE(RUN_REPORT_FILE_ID,1430)NUM_OFF_DIAGONALS
        WRITE(RUN_REPORT_FILE_ID,1431)NUM_JACOBIAN_ELEMENTS
        WRITE(RUN_REPORT_FILE_ID,1432)NONZEROS
        WRITE(RUN_REPORT_FILE_ID,1433)SPARSITY
        WRITE(RUN_REPORT_FILE_ID,1434)BLOCK_NONZEROS
        WRITE(RUN_REPORT_FILE_ID,1435)R_INNER
        WRITE(RUN_REPORT_FILE_ID,1436)R_OUTER
        WRITE(RUN_REPORT_FILE_ID,1437)NUM_R_QUAD_POINTS
        WRITE(RUN_REPORT_FILE_ID,1438)NUM_T_QUAD_POINTS
        WRITE(RUN_REPORT_FILE_ID,1439)NUM_P_QUAD_POINTS
        WRITE(RUN_REPORT_FILE_ID,1440)MAX_ITERATIONS
        WRITE(RUN_REPORT_FILE_ID,1441)CONVERGENCE_CRITERIA
        WRITE(RUN_REPORT_FILE_ID,'(/ / / / / / / /)')


    END IF
END IF


END SUBROUTINE OUTPUT_SETUP_TABLE





 !+101+####################################################################################!
!                                                                                           !
!       Poseidon_Initialize                                                                 !
!                                                                                           !
!===========================================================================================!
!                                                                                           !
!   Sets code parameters, allocates space, and initializes functions and variables needed   !
!   to run Poseidon.                                                                        !
!                                                                                           !
!-------------------------------------------------------------------------------------------!
!                                                                                           !
!       Input Variables     :                                                               !
!                                                                                           !
!       FEM_Degree_Input        -       Integer, Order of Finite Element Method solver to   !
!                                               be performed.                               !
!                                                                                           !
!       L_Limit_Input           -       Integer, Limit on the l value of the Spherical      !
!                                               Harmonic spectral decomposition.            !
!                                                                                           !
!       Inner_Radius            -       Real, Inner radius of the computational domain.     !
!                                                                                           !
!       Outer_Radius            -       Real, Outer radius of the computational domain.     !
!                                                                                           !
!       R_Elements_Input        -       Integer, Number of radial elements.                 !
!                                                                                           !
!       T_Elements_Input        -       Integer, Number of theta elements.                  !
!                                                                                           !
!       P_Elements_Input        -       Integer, Number of phi elements.                    !
!                                                                                           !
!       Input_Delta_R_Vector    -       Optional Real Vector, Dimension(1:R_Elements_Input) !
!                                       Each value corresponds to the radial length of the  !
!                                       each radial shell of elements starting from         !
!                                       Inner_Radius.                                       !
!                                                                                           !
!       Input_Delta_T_Vector    -       Optional Real Vector, Dimension(1:T_Elements_Input) !
!                                       Each value correspons to the angular width of each  !
!                                       theta wedge of elements starting from 0.            !
!                                                                                           !
!       Input_Delta_P_Vector    -       Optional Real Vector, Dimension(1:P_Elements_Input) !
!                                       Each value correspons to the angular width of each  !
!                                       phi wedge of elements starting from 0.              !
!                                                                                           !
 !#########################################################################################!
SUBROUTINE Poseidon_Initialize_1D( FEM_Degree, SH_Limit,                                                   &
                                Inner_Radius, Outer_Radius,                                             &
                                R_Elements_Input, T_Elements_Input, P_Elements_Input,                   &
                                Local_R_Elements_Input, Local_T_Elements_Input, Local_P_Elements_Input, &
                                Num_R_Quad_Input, Num_T_Quad_Input, Num_P_Quad_Input,                   &
                                Input_Delta_R_Vector, Input_Delta_T_Vector, Input_Delta_P_Vector        )



                                         !                          !
                                        !!      Input Variables     !!
                                         !

INTEGER, INTENT(IN)                                                             ::  FEM_Degree,             &
                                                                                    SH_Limit,               &
                                                                                    R_Elements_Input,       &
                                                                                    T_Elements_Input,       &
                                                                                    P_Elements_Input,       &
                                                                                    Local_R_Elements_Input, &
                                                                                    Local_T_Elements_Input, &
                                                                                    Local_P_Elements_Input, &
                                                                                    Num_R_Quad_Input,       &
                                                                                    Num_T_Quad_Input,       &
                                                                                    Num_P_Quad_Input

REAL(KIND = idp), INTENT(IN)                                                    ::  Inner_Radius,           &
                                                                                    Outer_Radius


REAL(KIND = idp), DIMENSION(1:R_Elements_Input),    OPTIONAL,   INTENT(IN)      ::  Input_Delta_R_Vector
REAL(KIND = idp), DIMENSION(1:T_Elements_Input),    OPTIONAL,   INTENT(IN)      ::  Input_Delta_T_Vector
REAL(KIND = idp), DIMENSION(1:P_Elements_Input),    OPTIONAL,   INTENT(IN)      ::  Input_Delta_P_Vector





                                         !                              !
                                        !!     Subroutine Variables     !!
                                         !                              !

REAL(KIND = idp), DIMENSION(1:R_Elements_Input)                                 ::  Delta_R_Vector
REAL(KIND = idp), DIMENSION(1:T_Elements_Input)                                 ::  Delta_T_Vector
REAL(KIND = idp), DIMENSION(1:P_Elements_Input)                                 ::  Delta_P_Vector


INTEGER                                                                         ::  td, pd
INTEGER                                                                         ::  nPROCS, tmpID
INTEGER                                                                         ::  l



REAL(KIND = idp)                                                                :: timea, timeb, timec


nPROCS = 1
nPROCs_POSEIDON = 1
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nPROCS,ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, tmpID, ierr)


DEGREE = FEM_Degree
L_LIMIT = SH_Limit
DOMAIN_DIM = 3



 !                                          !
!!  Set Global Variables to Input Values    !!
 !                                          !
R_INNER = Inner_Radius
R_OUTER = Outer_Radius



NUM_R_ELEMENTS = R_Elements_Input
NUM_T_ELEMENTS = T_Elements_Input
NUM_P_ELEMENTS = P_Elements_Input



NUM_LOC_R_ELEMENTS = Local_R_Elements_Input
NUM_LOC_T_ELEMENTS = Local_T_Elements_Input
NUM_LOC_P_ELEMENTS = Local_P_Elements_Input

NUM_SHELLS                  = 1
NUM_SUBSHELLS               = 1
NUM_SUBSHELLS_PER_SHELL     = 1
NUM_BLOCKS                  = 1
NUM_BLOCKS_PER_SHELL        = 1
NUM_BLOCK_THETA_ROWS        = 1
NUM_BLOCK_PHI_COLUMNS       = 1

NUM_R_ELEMS_PER_BLOCK       = NUM_R_ELEMENTS
NUM_T_ELEMS_PER_BLOCK       = 1
NUM_P_ELEMS_PER_BLOCK       = 1

NUM_R_ELEMS_PER_SHELL       = NUM_R_ELEMENTS
NUM_R_ELEMS_PER_SUBSHELL    = NUM_R_ELEMENTS


!
!   *_NUM_R_NODES - Number of total radial nodes in total, per block, per subshell
!
NUM_R_NODES             = DEGREE*NUM_R_ELEMENTS + 1
BLOCK_NUM_R_NODES       = DEGREE*NUM_R_ELEMS_PER_BLOCK + 1
SUBSHELL_NUM_R_NODES    = DEGREE*NUM_R_ELEMS_PER_SUBSHELL + 1



NUM_R_QUAD_POINTS = NUM_R_Quad_Input
NUM_T_QUAD_POINTS = NUM_T_Quad_Input
NUM_P_QUAD_POINTS = NUM_P_Quad_Input
!
!   NUM_TP_QUAD_POINTS = number of angular quadrature points per radial point
!
NUM_TP_QUAD_POINTS = NUM_T_QUAD_POINTS*NUM_P_QUAD_POINTS

CALL CHECK_SETUP(R_Elements_Input, T_Elements_Input, P_Elements_Input, TmpID )


!
!   Associate the Correct Map Functions, and Set Spherical Harmonic Length
!
LM_LENGTH = 1
Matrix_Location => CFA_3D_Matrix_Map
LM_Location => CFA_3D_LM_Map




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


NUM_OFF_DIAGONALS = ULM_LENGTH*(DEGREE + 1) - 1





                                 !                                      !
                                !!      Allocate Space for Poseidon     !!
                                 !                                      !
CALL Allocate_Poseidon_CFA_Variables()




                                 !                                          !
                                !!      Initialize Reusable Quadratures     !!
                                 !                                          !

CALL Initialize_LG_Quadrature(NUM_R_QUAD_POINTS, INT_R_LOCATIONS, INT_R_WEIGHTS)
CALL Initialize_LG_Quadrature(NUM_T_QUAD_POINTS, INT_T_LOCATIONS, INT_T_WEIGHTS)
CALL Initialize_LG_Quadrature(NUM_P_QUAD_POINTS, INT_P_LOCATIONS, INT_P_WEIGHTS)

IF ( NUM_T_QUAD_POINTS == 1 ) THEN
    ! 1D Cheat !
    INT_T_WEIGHTS = 4.0_idp/pi
END IF


DO td = 1,NUM_T_QUAD_POINTS
    DO pd = 1,NUM_P_QUAD_POINTS
        INT_TP_WEIGHTS( (td-1)*NUM_P_QUAD_POINTS+pd ) = INT_T_WEIGHTS(td)*INT_P_WEIGHTS(pd)
    END DO
END DO


                                 !                                          !
                                !!   Initialize Lagrange Polynomial Tables  !!
                                 !                                          !
CALL Initialize_Lagrange_Poly_Tables()


                                 !                                          !
                                !!      Set Initial Mesh (Optional)         !!
                                 !                                          !

CALL Generate_Defined_Coarse_Mesh(R_Elements_Input, NUM_R_ELEMENTS, R_Coarsen_Factor,     &
                                  R_INNER, Input_Delta_R_Vector, rlocs, drlocs)

RADIAL_MESH_SET_FLAG = .TRUE.


CALL Generate_Defined_Coarse_Mesh(T_Elements_Input, NUM_T_ELEMENTS, T_Coarsen_Factor,     &
                                  0.0_idp, (/ pi /), tlocs, dtlocs)
THETA_MESH_SET_FLAG = .TRUE.


CALL Generate_Defined_Coarse_Mesh(P_Elements_Input, NUM_P_ELEMENTS, P_Coarsen_Factor,     &
                                  0.0_idp, (/ 2.0_idp*pi /), plocs, dplocs)
PHI_MESH_SET_FLAG = .TRUE.





IF ( DOMAIN_DIM == 1 ) THEN
     M_VALUES = 0
ELSE IF ( DOMAIN_DIM == 2 ) THEN
     M_VALUES = 0
ELSE IF ( DOMAIN_DIM == 3 ) THEN
     M_VALUES = (/(l,l=0,L_LIMIT,1)/)
END IF


Ratio_T_BNDLperBLCK = NUM_T_ELEMENTS/( NUM_BLOCK_THETA_ROWS*NUM_LOC_T_ELEMENTS )
Ratio_P_BNDLperBLCK = NUM_P_ELEMENTS/( NUM_BLOCK_PHI_COLUMNS*NUM_LOC_P_ELEMENTS )
Ratio_BNDLperBLCK = Ratio_T_BNDLperBLCK * Ratio_P_BNDLperBLCK

CALL CREATE_POSEIDON_COMMUNICATORS( DATA_DIST_MODE )

CALL Initialize_Ylm_Tables()



!CALL OUTPUT_SETUP_TABLE( nPROCS, R_Elements_Input, T_Elements_Input, P_Elements_Input )



END SUBROUTINE Poseidon_Initialize_1D





 !+101+####################################################################################!
!                                                                                           !
!       Poseidon_Initialize                                                                 !
!                                                                                           !
!===========================================================================================!
!                                                                                           !
!   Sets code parameters, allocates space, and initializes functions and variables needed   !
!   to run Poseidon.                                                                        !
!                                                                                           !
!-------------------------------------------------------------------------------------------!
!                                                                                           !
!       Input Variables     :                                                               !
!                                                                                           !
!       FEM_Degree_Input        -       Integer, Order of Finite Element Method solver to   !
!                                               be performed.                               !
!                                                                                           !
!       L_Limit_Input           -       Integer, Limit on the l value of the Spherical      !
!                                               Harmonic spectral decomposition.            !
!                                                                                           !
!       Inner_Radius            -       Real, Inner radius of the computational domain.     !
!                                                                                           !
!       Outer_Radius            -       Real, Outer radius of the computational domain.     !
!                                                                                           !
!       R_Elements_Input        -       Integer, Number of radial elements.                 !
!                                                                                           !
!       T_Elements_Input        -       Integer, Number of theta elements.                  !
!                                                                                           !
!       P_Elements_Input        -       Integer, Number of phi elements.                    !
!                                                                                           !
!       Input_Delta_R_Vector    -       Optional Real Vector, Dimension(1:R_Elements_Input) !
!                                       Each value corresponds to the radial length of the  !
!                                       each radial shell of elements starting from         !
!                                       Inner_Radius.                                       !
!                                                                                           !
!       Input_Delta_T_Vector    -       Optional Real Vector, Dimension(1:T_Elements_Input) !
!                                       Each value correspons to the angular width of each  !
!                                       theta wedge of elements starting from 0.            !
!                                                                                           !
!       Input_Delta_P_Vector    -       Optional Real Vector, Dimension(1:P_Elements_Input) !
!                                       Each value correspons to the angular width of each  !
!                                       phi wedge of elements starting from 0.              !
!                                                                                           !
 !#########################################################################################!
SUBROUTINE Poseidon_Initialize_3D( FEM_Degree, SH_Limit,                                                   &
                                Inner_Radius, Outer_Radius,                                             &
                                R_Elements_Input, T_Elements_Input, P_Elements_Input,                   &
                                Local_R_Elements_Input, Local_T_Elements_Input, Local_P_Elements_Input, &
                                Num_R_Quad_Input, Num_T_Quad_Input, Num_P_Quad_Input,                   &
                                Input_Delta_R_Vector, Input_Delta_T_Vector, Input_Delta_P_Vector        )



                                         !                          !
                                        !!      Input Variables     !!
                                         !

INTEGER, INTENT(IN)                                                             ::  FEM_Degree,             &
                                                                                    SH_Limit,               &
                                                                                    R_Elements_Input,       &
                                                                                    T_Elements_Input,       &
                                                                                    P_Elements_Input,       &
                                                                                    Local_R_Elements_Input, &
                                                                                    Local_T_Elements_Input, &
                                                                                    Local_P_Elements_Input, &
                                                                                    Num_R_Quad_Input,       &
                                                                                    Num_T_Quad_Input,       &
                                                                                    Num_P_Quad_Input

REAL(KIND = idp), INTENT(IN)                                                    ::  Inner_Radius,           &
                                                                                    Outer_Radius


REAL(KIND = idp), DIMENSION(1:R_Elements_Input),    OPTIONAL,   INTENT(IN)      ::  Input_Delta_R_Vector
REAL(KIND = idp), DIMENSION(1:T_Elements_Input),    OPTIONAL,   INTENT(IN)      ::  Input_Delta_T_Vector
REAL(KIND = idp), DIMENSION(1:P_Elements_Input),    OPTIONAL,   INTENT(IN)      ::  Input_Delta_P_Vector





                                         !                              !
                                        !!     Subroutine Variables     !!
                                         !                              !

REAL(KIND = idp), DIMENSION(1:R_Elements_Input)                                 ::  Delta_R_Vector
REAL(KIND = idp), DIMENSION(1:T_Elements_Input)                                 ::  Delta_T_Vector
REAL(KIND = idp), DIMENSION(1:P_Elements_Input)                                 ::  Delta_P_Vector


INTEGER                                                                         ::  td, pd
INTEGER                                                                         ::  nPROCS, tmpID
INTEGER                                                                         ::  l



REAL(KIND = idp)                                                                :: timea, timeb, timec



CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nPROCS,ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, tmpID, ierr)


DEGREE = FEM_Degree
L_LIMIT = SH_Limit


 !                                          !
!!  Set Global Variables to Input Values    !!
 !                                          !
R_INNER = Inner_Radius
R_OUTER = Outer_Radius




NUM_R_ELEMENTS = R_Elements_Input
NUM_T_ELEMENTS = T_Elements_Input
NUM_P_ELEMENTS = P_Elements_Input

NUM_LOC_R_ELEMENTS = Local_R_Elements_Input
NUM_LOC_T_ELEMENTS = Local_T_Elements_Input
NUM_LOC_P_ELEMENTS = Local_P_Elements_Input

!
!   *_NUM_R_NODES - Number of total radial nodes in total, per block, per subshell
!
NUM_R_NODES             = DEGREE*NUM_R_ELEMENTS + 1
BLOCK_NUM_R_NODES       = DEGREE*NUM_R_ELEMS_PER_BLOCK + 1
SUBSHELL_NUM_R_NODES    = DEGREE*NUM_R_ELEMS_PER_SUBSHELL + 1



NUM_R_QUAD_POINTS = NUM_R_Quad_Input
NUM_T_QUAD_POINTS = NUM_T_Quad_Input
NUM_P_QUAD_POINTS = NUM_P_Quad_Input
!
!   NUM_TP_QUAD_POINTS = number of angular quadrature points per radial point
!
NUM_TP_QUAD_POINTS = NUM_T_QUAD_POINTS*NUM_P_QUAD_POINTS


CALL CHECK_SETUP(R_Elements_Input, T_Elements_Input, P_Elements_Input, TmpID )




!
!   Associate the Correct Map Functions, and Set Spherical Harmonic Length
!
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
IF ( PHYSICS_TYPE == "NEW" ) THEN

    ULM_LENGTH = LM_LENGTH
    PROB_DIM = VAR_DIM
    ELEM_PROB_DIM = ELEM_VAR_DIM
    ELEM_PROB_DIM_SQR = ELEM_VAR_DIM*ELEM_VAR_DIM
    BLOCK_PROB_DIM = BLOCK_VAR_DIM
    SUBSHELL_PROB_DIM = SUBSHELL_VAR_DIM

ELSE IF ( PHYSICS_TYPE == "CFA" ) THEN

    ULM_LENGTH          = NUM_CFA_VARS*LM_LENGTH
    PROB_DIM            = NUM_CFA_VARS*VAR_DIM
    ELEM_PROB_DIM       = NUM_CFA_VARS*ELEM_VAR_DIM
    ELEM_PROB_DIM_SQR   = ELEM_PROB_DIM*ELEM_PROB_DIM
    BLOCK_PROB_DIM      = NUM_CFA_VARS*BLOCK_VAR_DIM
    SUBSHELL_PROB_DIM   = NUM_CFA_VARS*SUBSHELL_VAR_DIM

END IF


NUM_OFF_DIAGONALS = ULM_LENGTH*(DEGREE + 1) - 1






                                 !                                      !
                                !!      Allocate Space for Poseidon     !!
                                 !                                      !
CALL Allocate_Poseidon_CFA_Variables()








                                 !                                          !
                                !!      Initialize Reusable Quadratures     !!
                                 !                                          !

CALL Initialize_LG_Quadrature(NUM_R_QUAD_POINTS, INT_R_LOCATIONS, INT_R_WEIGHTS)
CALL Initialize_LG_Quadrature(NUM_T_QUAD_POINTS, INT_T_LOCATIONS, INT_T_WEIGHTS)
CALL Initialize_LG_Quadrature(NUM_P_QUAD_POINTS, INT_P_LOCATIONS, INT_P_WEIGHTS)

IF ( NUM_T_QUAD_POINTS == 1 ) THEN
    ! 1D Cheat !
    INT_T_WEIGHTS = 4.0_idp/pi
END IF


DO td = 1,NUM_T_QUAD_POINTS
    DO pd = 1,NUM_P_QUAD_POINTS
        INT_TP_WEIGHTS( (td-1)*NUM_P_QUAD_POINTS+pd ) = INT_T_WEIGHTS(td)*INT_P_WEIGHTS(pd)
    END DO
END DO



                                 !                                          !
                                !!   Initialize Lagrange Polynomial Tables  !!
                                 !                                          !
CALL Initialize_Lagrange_Poly_Tables()








                                 !                                          !
                                !!      Set Initial Mesh (Optional)         !!
                                 !                                          !

CALL Generate_Defined_Coarse_Mesh(R_Elements_Input, NUM_R_ELEMENTS, R_Coarsen_Factor,     &
                                  R_INNER, Input_Delta_R_Vector, rlocs, drlocs)

RADIAL_MESH_SET_FLAG = .TRUE.


CALL Generate_Defined_Coarse_Mesh(T_Elements_Input, NUM_T_ELEMENTS, T_Coarsen_Factor,     &
                                  0.0_idp, Input_Delta_T_Vector, tlocs, dtlocs)
THETA_MESH_SET_FLAG = .TRUE.


CALL Generate_Defined_Coarse_Mesh(P_Elements_Input, NUM_P_ELEMENTS, P_Coarsen_Factor,     &
                                  0.0_idp, Input_Delta_P_Vector, plocs, dplocs)
PHI_MESH_SET_FLAG = .TRUE.





IF ( DOMAIN_DIM == 1 ) THEN
     M_VALUES = 0
ELSE IF ( DOMAIN_DIM == 2 ) THEN
     M_VALUES = 0
ELSE IF ( DOMAIN_DIM == 3 ) THEN
     M_VALUES = (/(l,l=0,L_LIMIT,1)/)
END IF



Ratio_T_BNDLperBLCK = NUM_T_ELEMENTS/( NUM_BLOCK_THETA_ROWS*NUM_LOC_T_ELEMENTS )
Ratio_P_BNDLperBLCK = NUM_P_ELEMENTS/( NUM_BLOCK_PHI_COLUMNS*NUM_LOC_P_ELEMENTS )
Ratio_BNDLperBLCK = Ratio_T_BNDLperBLCK * Ratio_P_BNDLperBLCK


CALL CREATE_POSEIDON_COMMUNICATORS( DATA_DIST_MODE )

CALL Initialize_Ylm_Tables()



CALL OUTPUT_SETUP_TABLE( nPROCS, R_Elements_Input, T_Elements_Input, P_Elements_Input )



END SUBROUTINE Poseidon_Initialize_3D





END MODULE Poseidon_Initialization_Module
