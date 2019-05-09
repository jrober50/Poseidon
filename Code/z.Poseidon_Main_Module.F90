   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_Main_Module                                                         !##!
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
!##!    +102+   Poseidon_Run                                                        !##!
!##!    +103+   Poseidon_Close                                                      !##!
!##!    +104+   Poseidon_Set_Mesh                                                   !##!
!##!    +105+   Poseidon_CFA_Set_Boundary_Condtions                                 !##!
!##!                                                                                !##!
!##!    +201+   Poseidon_Readiness_Check                                            !##!
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
                    WRITE_REPORT_FLAG,          &
                    RUN_REPORT_FILE_ID



USE Global_Variables_And_Parameters, &
            ONLY :  R_INNER, R_OUTER,                               &
                    NUM_R_ELEMENTS, NUM_T_ELEMENTS, NUM_P_ELEMENTS, &
                    NUM_R_NODES,                                    &
                    BLOCK_NUM_R_NODES,                              &
                    SUBSHELL_NUM_R_NODES,                           &
                    rlocs, tlocs, plocs,                            &
                    RHS_Vector,                                     &
                    Coefficient_Vector,                             &
                    Source_Term_Coefficients,                       &
                    Test_Space_Allocated_Flag,                      &
                    Stiffness_Matrix_Initialized_Flag,              &
                    FirstCall_Flag,                                 &
                    INNER_BC_SET_FLAG, OUTER_BC_SET_FLAG,           &
                    INNER_BC_TYPE, OUTER_BC_TYPE,                   &
                    INNER_DIR_BC_INPUT, INNER_NEU_BC_INPUT,         &
                    OUTER_DIR_BC_INPUT, OUTER_NEU_BC_INPUT,         &
                    INNER_UNIFORM_DIR_BC_FLAG,                      &
                    OUTER_UNIFORM_DIR_BC_FLAG,                      &
                    RADIAL_MESH_SET_FLAG,                           &
                    THETA_MESH_SET_FLAG,                            &
                    PHI_MESH_SET_FLAG,                              &
                    INT_R_LOCATIONS,                                &
                    INT_R_WEIGHTS,                                  &
                    INT_T_LOCATIONS,                                &
                    INT_T_WEIGHTS,                                  &
                    INT_P_LOCATIONS,                                &
                    INT_P_WEIGHTS,                                  &
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
                    NONZEROS,                                       &
                    BLOCK_NONZEROS,                                 &
                    RUN_TIME_TABLE






USE Additional_Functions_Module, &
            ONLY :  Spherical_Harmonic,                             &
                    Map_To_X_Space, Map_From_X_Space,               &
                    Initialize_LG_Quadrature,                       &
                    Initialize_LG_Quadrature_Locations,             &
                    CFA_3D_Matrix_Map,                              &
                    CFA_2D_Matrix_Map,                              &
                    CFA_1D_Matrix_Map,                              &
                    CFA_1D_LM_Map,                                  &
                    CFA_2D_LM_Map,                                  &
                    CFA_3D_LM_Map,                                  &
                    Generate_Defined_Mesh,                          &
                    Generate_Defined_Coarse_Mesh


USE Allocate_Variables_Module, &
            ONLY :  Allocate_Poseidon_CFA_Variables,                &
                    Deallocate_Poseidon_CFA_Variables


USE Jacobian_Internal_Functions_Module,  &
            ONLY :  Initialize_Guess_Values,                        &
                    Initialize_Ylm_Table,                           &
                    Initialize_Ylm_Tables,                          &
                    Initialize_Lagrange_Poly_Tables

USE CFA_Newton_Raphson_Module, &
            ONLY :  CFA_Newton_Raphson

USE CFA_3D_Master_Build_Module, &
                        ONLY :  Calc_3D_Values_At_Location

USE Poseidon_MPI_Module, &
            ONLY :  CREATE_POSEIDON_COMMUNICATORS

USE Poseidon_Parameter_Read_Module, &
            ONLY :  UNPACK_POSEIDON_PARAMETERS


USE IO_Functions_Module, &
            ONLY :  OPEN_RUN_REPORT_FILE,                           &
                    CLOSE_RUN_REPORT_FILE

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

SUBROUTINE Poseidon_Initialize( Inner_Radius, Outer_Radius,                                             &
                                R_Elements_Input, T_Elements_Input, P_Elements_Input,                   &
                                Local_R_Elements_Input, Local_T_Elements_Input, Local_P_Elements_Input, &
                                Input_Delta_R_Vector, Input_Delta_T_Vector, Input_Delta_P_Vector        )



                                         !                          !
                                        !!      Input Variables     !!
                                         !                          

INTEGER, INTENT(IN)                                                             ::  R_Elements_Input,       &
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





INTEGER                                                                         ::  tmpID, nPROCS
INTEGER                                                                         ::  l 
INTEGER(8) :: TMPA, TMPB, TMPC

INTEGER    ::  NON_ZEROS_B

REAL(KIND = idp)            :: timea, timeb, timec



CALL MPI_COMM_RANK(MPI_COMM_WORLD, tmpID, ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nPROCS, ierr)



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

!
!   *_NUM_R_NODES - Number of total radial nodes in total, per block, per subshell
!
NUM_R_NODES             = DEGREE*NUM_R_ELEMENTS + 1
BLOCK_NUM_R_NODES       = DEGREE*NUM_R_ELEMS_PER_BLOCK + 1
SUBSHELL_NUM_R_NODES    = DEGREE*NUM_R_ELEMS_PER_SUBSHELL + 1






IF (NUM_R_ELEMS_PER_SHELL*NUM_SHELLS .NE. NUM_R_ELEMENTS ) THEN

    IF ( tmpID == 0) THEN
          PRINT*,"Domain not properly broken into shells.  Number of radial elements does not match."
          PRINT*, " # of Elements Per Shell : ",NUM_R_ELEMS_PER_SHELL ,             &
              " # Shells : ",NUM_SHELLS ,                                         &
              " # Global Radial Elements",NUM_R_ELEMENTS
          PRINT*,"Checked in z.Poseidon_Main_Modules.F90, line 330 "
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
       PRINT*,"Checked in z.Poseidon_Main_Modules.F90, line 344 "
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
       PRINT*,"Checked in z.Poseidon_Main_Modules.F90, line 344 "
    END IF

    CALL MPI_FINALIZE(ierr)
    STOP

END IF




IF ( NUM_BLOCKS .NE. nPROCS_POSEIDON ) THEN


    IF ( tmpID == 0) THEN
       PRINT*,"Poseidon's number of processes does not match the number of blocks"
       PRINT*, " # of Blocks : ",NUM_BLOCKS,                        &
       " # of Poseidon Processes expected : ",nPROCs_POSEIDON
       PRINT*,"Checked in z.Poseidon_Main_Modules.F90, line 344 "
    END IF
    CALL MPI_FINALIZE(ierr)
    STOP

END IF




IF ( NUM_R_ELEMS_PER_BLOCK.NE. NUM_R_ELEMS_PER_SHELL ) THEN

    IF ( tmpID == 0) THEN
       PRINT*,"The number of elements per block does not match."
       PRINT*, " # of Elements per Block : ",NUM_R_ELEMS_PER_BLOCK
       PRINT*, " # of Radial Elements Per Shell : ",NUM_R_ELEMS_PER_SHELL
       PRINT*,"Checked in z.Poseidon_Main_Modules.F90, line 393 "
    END IF
    CALL MPI_FINALIZE(ierr)
    STOP

END IF






IF ( R_COARSEN_FACTOR*NUM_R_ELEMS_PER_SHELL*NUM_SHELLS .NE. R_Elements_Input) THEN

     IF (tmpID == 0) THEN 

        PRINT*,"Poseidon's Radial Decomposition is incompatible with Input"
        PRINT*,"Expected Global Radial Elements",R_COARSEN_FACTOR*NUM_R_ELEMS_PER_SHELL*NUM_SHELLS
        PRINT*,"Input Global Radial Elements",R_Elements_Input
        PRINT*,"Checked in Z.Poseidon_Main_Module.F90, Poseidon_Initialize"

     END IF 
     CALL MPI_FINALIZE(ierr)
     STOP

END IF 

IF ( T_COARSEN_FACTOR*NUM_T_ELEMS_PER_BLOCK*NUM_BLOCK_THETA_ROWS .NE. T_Elements_Input) THEN

     IF (tmpID == 0) THEN

        PRINT*,"Poseidon's Phi Decomposition is incompatible with Input"
        PRINT*,"Expected Global Theta Elements",T_COARSEN_FACTOR*NUM_T_ELEMS_PER_BLOCK*NUM_BLOCK_THETA_ROWS
        PRINT*,"Input Global Theta Elements",T_Elements_Input
        PRINT*,"Checked in Z.Poseidon_Main_Module.F90, Poseidon_Initialize"

     END IF
     CALL MPI_FINALIZE(ierr)
     STOP

END IF


IF ( P_COARSEN_FACTOR*NUM_P_ELEMS_PER_BLOCK*NUM_BLOCK_PHI_COLUMNS .NE. P_Elements_Input) THEN

     IF (tmpID == 0) THEN

        PRINT*,"Poseidon's Phi Decomposition is incompatible with Input"
        PRINT*,"Expected Global Phi Elements",P_COARSEN_FACTOR*NUM_P_ELEMS_PER_BLOCK*NUM_BLOCK_PHI_COLUMNS
        PRINT*,"Input Global Phi Elements",P_Elements_Input
        PRINT*,"Checked in Z.Poseidon_Main_Module.F90, Poseidon_Initialize"

     END IF 
     CALL MPI_FINALIZE(ierr)
     STOP

END IF 








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

TMPA = 3*NUM_OFF_DIAGONALS + 1
TMPB = PROB_DIM
TMPC = TMPA*TMPB



NUM_OFF_DIAGONALS = ULM_LENGTH*(DEGREE + 1) - 1

TMPB = NUM_OFF_DIAGONALS
BLOCK_NONZEROS = (2*NUM_OFF_DIAGONALS + 1)* BLOCK_PROB_DIM


NONZEROS = NUM_R_ELEMENTS*ELEM_PROB_DIM_SQR          &
         - (NUM_R_ELEMENTS-1)*ULM_LENGTH*ULM_LENGTH

BLOCK_NONZEROS = NUM_R_ELEMS_PER_BLOCK*ELEM_PROB_DIM_SQR          &
               - (NUM_R_ELEMS_PER_BLOCK-1)*ULM_LENGTH*ULM_LENGTH







IF ( tmpID == 0 ) THEN

    PRINT*,"------------ POSEIDON PARAMETERS -------------"
    PRINT*,"       ACTIVE DIMENSIONS = ",DOMAIN_DIM
    PRINT*,"                  DEGREE = ",DEGREE
    PRINT*,"                 L_LIMIT = ",L_LIMIT
    PRINT*,"            NUM_CFA_VARS = ",NUM_CFA_VARS
    PRINT*,"                  nPROCS = ",nPROCS
    PRINT*,"         nPROCS_POSEIDON = ",nPROCS_POSEIDON
    PRINT*," CHIMERA RADIAL ELEMENTS = ",R_ELEMENTS_INPUT
    PRINT*,"  CHIMERA THETA ELEMENTS = ",T_ELEMENTS_INPUT
    PRINT*,"    CHIMERA PHI ELEMENTS = ",P_ELEMENTS_INPUT
    PRINT*,"POSEIDON RADIAL ELEMENTS = ",NUM_R_ELEMENTS
    PRINT*," POSEIDON THETA ELEMENTS = ",NUM_T_ELEMENTS
    PRINT*,"   POSEIDON PHI ELEMENTS = ",NUM_P_ELEMENTS
    PRINT*,"        R_COARSEN_FACTOR = ",R_COARSEN_FACTOR
    PRINT*,"        T_COARSEN_FACTOR = ",T_COARSEN_FACTOR
    PRINT*,"        P_COARSEN_FACTOR = ",P_COARSEN_FACTOR
    PRINT*,"             NUM_R_NODES = ",NUM_R_NODES
    PRINT*,"        NUMBER OF SHELLS = ",NUM_SHELLS
    PRINT*,"     NUMBER OF SUBSHELLS = ",NUM_SUBSHELLS
    PRINT*,"     SUBSHELLS PER SHELL = ",NUM_SUBSHELLS_PER_SHELL
    PRINT*,"    NUM BLOCKS PER SHELL = ",NUM_BLOCKS_PER_SHELL
    PRINT*,"   NUM_R_ELEMS_PER_SHELL = ",NUM_R_ELEMS_PER_SHELL
    PRINT*,"   NUM_T_ELEMS_PER_BLOCK = ",NUM_T_ELEMS_PER_BLOCK
    PRINT*,"   NUM_P_ELEMS_PER_BLOCK = ",NUM_P_ELEMS_PER_BLOCK
    PRINT*,"                 VAR_DIM = ",VAR_DIM
    PRINT*,"           BLOCK_VAR_DIM = ",BLOCK_VAR_DIM
    PRINT*,"                PROB_DIM = ",PROB_DIM
    PRINT*,"          BLOCK_PROB_DIM = ",BLOCK_PROB_DIM
    PRINT*,"       SUBSHELL_PROB_DIM = ",SUBSHELL_PROB_DIM
    PRINT*,"       NUM_OFF_DIAGONALS = ",NUM_OFF_DIAGONALS
    PRINT*,"        NUM STF ELEMENTS = ",TMPC
    PRINT*,"    NUM NNZ STF ELEMENTS = ",NONZEROS
    PRINT*,"       NUM NNZ PER BLOCK = ",BLOCK_NONZEROS
    PRINT*,"            Inner Radius = ",R_INNER
    PRINT*,"            Outer Radius = ",R_OUTER
    PRINT*,"    # Radial Quad Points = ",NUM_R_QUAD_POINTS
    PRINT*,"     # Theta Quad Points = ",NUM_T_QUAD_POINTS
    PRINT*,"       # Phi Quad Points = ",NUM_P_QUAD_POINTS
    PRINT*,"      Maximum Iterations = ",MAX_ITERATIONS
    PRINT*,"    Convergence Criteria = ",CONVERGENCE_CRITERIA
END IF





                                 !                                      !
                                !!      Allocate Space for Poseidon     !!
                                 !                                      !






!!!!  Allocate Data Space !!!!
CALL Allocate_Poseidon_CFA_Variables()





CALL Initialize_LG_Quadrature(NUM_R_QUAD_POINTS, INT_R_LOCATIONS, INT_R_WEIGHTS)
CALL Initialize_LG_Quadrature(NUM_T_QUAD_POINTS, INT_T_LOCATIONS, INT_T_WEIGHTS)
CALL Initialize_LG_Quadrature(NUM_P_QUAD_POINTS, INT_P_LOCATIONS, INT_P_WEIGHTS)






CALL Initialize_Lagrange_Poly_Tables()



                     !                                          !
                    !!      Set Initial Mesh (Optional)         !!
                     !                                          !




CALL Generate_Defined_Coarse_Mesh(R_Elements_Input, NUM_R_ELEMENTS, R_Coarsen_Factor,     &
                                  R_INNER, Input_Delta_R_Vector, rlocs)
RADIAL_MESH_SET_FLAG = .TRUE.


CALL Generate_Defined_Coarse_Mesh(T_Elements_Input, NUM_T_ELEMENTS, T_Coarsen_Factor,     &
                                  0.0_idp, Input_Delta_T_Vector, tlocs)
THETA_MESH_SET_FLAG = .TRUE.


CALL Generate_Defined_Coarse_Mesh(P_Elements_Input, NUM_P_ELEMENTS, P_Coarsen_Factor,     &
                                  0.0_idp, Input_Delta_P_Vector, plocs)
PHI_MESH_SET_FLAG = .TRUE.






IF ( DOMAIN_DIM == 1 ) THEN
     M_VALUES = 0
ELSE IF ( DOMAIN_DIM == 2 ) THEN
     M_VALUES = 0
ELSE IF ( DOMAIN_DIM == 3 ) THEN
     M_VALUES = (/(l,l=0,L_LIMIT,1)/)
END IF



CALL CREATE_POSEIDON_COMMUNICATORS( DATA_DIST_MODE )


CALL Initialize_Ylm_Tables()


END SUBROUTINE Poseidon_Initialize















!+102+#######################################################################################!
!                                                                                           !
!       Poseidon_Run - Calculates Source Vector and Solves for solution coefficients        !
!                                                                                           !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Run()

LOGICAL                                             ::  Readiness_Flag


IF ( myID_Poseidon == 0 ) THEN
    CALL OPEN_RUN_REPORT_FILE()
END IF



!CALL Poseidon_Readiness_Check(Readiness_Flag)
Readiness_Flag = .TRUE.

IF ( Readiness_Flag ) THEN

    CALL CFA_Newton_Raphson()

ELSE

    PRINT*, "ERROR IN POSEIDON : There was an error in setting up Poseidon, therefore it did not run."

END IF

    IF ( myID_Poseidon == 0 ) THEN
        CALL OUTPUT_RUN_REPORT(Cur_Iteration, myID_Poseidon)
        CALL CLOSE_RUN_REPORT_FILE()
    END IF



END SUBROUTINE Poseidon_Run
















!+103+######################################################################################!
!                                                                                           !
!       Poseidon_Close - Deallocate all Poseidon variables.                                 !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Close()



!!!!  Deallocate Data Space !!!!
CALL Deallocate_Poseidon_CFA_Variables

DEALLOCATE(LOCAL_NODE_LOCATIONS)


Stiffness_Matrix_Initialized_Flag = .FALSE.
Test_Space_Allocated_Flag = .FALSE.

FirstCall_Flag = .TRUE.


RADIAL_MESH_SET_FLAG = .FALSE.
THETA_MESH_SET_FLAG = .FALSE.
PHI_MESH_SET_FLAG = .FALSE.

INNER_BC_SET_FLAG = .FALSE.
OUTER_BC_SET_FLAG = .FALSE.


END SUBROUTINE Poseidon_Close
















 !+104+####################################################################################!
!                                                                                           !
!      Poseidon_Set_Mesh                                                                    !
!                                                                                           !
!-------------------------------------------------------------------------------------------!
!                                                                                           !
!   This subroutine sets the locations of the element end locations using vectors giving    !
!   the length of elements in each dimension.                                               !
!                                                                                           !
!-------------------------------------------------------------------------------------------!
!                                                                                           !
!   Input Variables     :      * = R, T, P                                                  !
!                                                                                           !
!   Input_Delta_*_Vector    -   Optional, Real Vector, Dimension(1:NUM_*_ELEMENTS)          !
!                               Gives the length of elements in the * dimension.            !
!                                                                                           !
 !#########################################################################################!
SUBROUTINE Poseidon_Set_Mesh(Input_Delta_R_Vector, Input_Delta_T_Vector, Input_Delta_P_Vector)



REAL(KIND = idp), DIMENSION(1:NUM_R_ELEMENTS),    OPTIONAL,   INTENT(IN)      ::  Input_Delta_R_Vector
REAL(KIND = idp), DIMENSION(1:NUM_T_ELEMENTS),    OPTIONAL,   INTENT(IN)      ::  Input_Delta_T_Vector
REAL(KIND = idp), DIMENSION(1:NUM_P_ELEMENTS),    OPTIONAL,   INTENT(IN)      ::  Input_Delta_P_Vector

REAL(KIND = idp), DIMENSION(1:NUM_R_ELEMENTS)                                  ::  TMP_Delta_R_Vector
REAL(KIND = idp), DIMENSION(1:NUM_T_ELEMENTS)                                  ::  TMP_Delta_T_Vector
REAL(KIND = idp), DIMENSION(1:NUM_P_ELEMENTS)                                  ::  TMP_Delta_P_Vector



  !                                                                                                 !
 !!     If Input_Delta_R_Vector is present then user wishes to define a non-uniform radial mesh.    !!
  !                                                                                                 !
IF (    PRESENT(Input_Delta_R_Vector)   ) THEN




    CALL Generate_Defined_Mesh(NUM_R_ELEMENTS, R_INNER, Input_Delta_R_Vector, rlocs)
    RADIAL_MESH_SET_FLAG = .TRUE.

     !                                                                              !
    !!  If the stiffness matrix has been initialized then this call is resetting    !!
    !!  the radial mesh, which means the stiffness matrix needs to be regenerated.  !!
     !                                                                              !
    IF ( Stiffness_Matrix_Initialized_Flag .EQV. .TRUE. ) THEN


        IF ( PHYSICS_TYPE == 'NEW' ) THEN


            !
            !   The Newtonian stiffness matrix is indepenent of input variables
            !   and therefore we can go ahead and generate it.
            !
        !   CALL Generate_Newton_Stiffness_Matrix()

        END IF



    END IF




END IF





 !                                                                                                 !
!!     If Input_Delta_R_Vector is present then user wishes to define a non-uniform radial mesh.    !!
 !                                                                                                 !
IF (    PRESENT(Input_Delta_T_Vector)   ) THEN


    CALL Generate_Defined_Mesh(NUM_T_ELEMENTS, 0.0_idp, Input_Delta_T_Vector, tlocs)
    THETA_MESH_SET_FLAG = .TRUE.


ELSE IF ( ( THETA_MESH_SET_FLAG .EQV. .FALSE. ) .AND. ( .NOT. PRESENT(Input_Delta_T_Vector) ) ) THEN


    TMP_Delta_T_Vector = (pi)/ REAL(NUM_T_ELEMENTS)

    CALL Generate_Defined_Mesh(NUM_T_ELEMENTS, 0.0_idp, (/ pi /), tlocs)
    THETA_MESH_SET_FLAG = .TRUE.



END IF





 !                                                                                                 !
!!     If Input_Delta_R_Vector is present then user wishes to define a non-uniform radial mesh.    !!
 !                                                                                                 !
IF (    PRESENT(Input_Delta_P_Vector)   ) THEN


    CALL Generate_Defined_Mesh(NUM_P_ELEMENTS, 0.0_idp, Input_Delta_P_Vector, plocs)
    PHI_MESH_SET_FLAG = .TRUE.


ELSE IF ( ( PHI_MESH_SET_FLAG .EQV. .FALSE. ) .AND. ( .NOT. PRESENT(Input_Delta_P_Vector) ) ) THEN

    TMP_Delta_P_Vector = (2*pi) / REAL(Num_P_Elements )


    CALL Generate_Defined_Mesh(NUM_P_ELEMENTS, 0.0_idp, TMP_Delta_P_Vector, plocs)
    PHI_MESH_SET_FLAG = .TRUE.

END IF



END SUBROUTINE Poseidon_Set_Mesh










 !+105+####################################################################################!
!                                                                                           !
!      Poseidon_CFA_Set_Boundary_Condtion                                                       !
!                                                                                           !
!-------------------------------------------------------------------------------------------!
!                                                                                           !
!   Set boundary condition values for the system.  This needs to be done once, but can be   !
!   repeatedly used if the boundary conditions are changed.                                 !
!                                                                                           !
!-------------------------------------------------------------------------------------------!
!                                                                                           !
!   Input Variables     :                                                                   !
!                                                                                           !
!           BC_Location_Input       -   "I" for Inner Boundary                              !
!                                       "O" for Outer Boundary                              !
!                                                                                           !
!           BC_Type_Input           -   "D" for Dirichlet Boundary Condition                !
!                                       "N" for Neumann Boundary Condition                  !
!                                                                                           !
!           BC_Value_Input          -   For a Dirichlet Boundary Condition specify the      !
!                                           Newtonian potential at the boundary.            !
!                                       For a Neumann Boundary Condition specify the        !
!                                           radial derivative of the potential at the       !
!                                           boundary.                                       !
!                                                                                           !
 !#########################################################################################!
SUBROUTINE Poseidon_CFA_Set_Uniform_Boundary_Conditions(BC_Location_Input, BC_Type_Input, BC_Value_Input)



CHARACTER(LEN = 1), INTENT(IN)                                  ::  BC_Location_Input
CHARACTER(LEN = 1), DIMENSION(1:NUM_CFA_VARS), INTENT(IN)       ::  BC_Type_Input

REAL(KIND = idp), DIMENSION(1:NUM_CFA_VARS), INTENT(IN)         ::  BC_Value_Input



INTEGER                         ::  ui



IF (    BC_Location_Input == "I"    ) THEN

    INNER_CFA_BC_TYPE(1:NUM_CFA_VARS) = BC_Type_Input(1:NUM_CFA_VARS)
    INNER_CFA_BC_VALUES(1:NUM_CFA_VARS) = BC_Value_Input(1:NUM_CFA_VARS)




ELSE IF (    BC_Location_Input == "O"    ) THEN

    OUTER_CFA_BC_TYPE(1:NUM_CFA_VARS) = BC_Type_Input(1:NUM_CFA_VARS)
    OUTER_CFA_BC_VALUES(1:NUM_CFA_VARS) = BC_Value_Input(1:NUM_CFA_VARS)


END IF







END SUBROUTINE Poseidon_CFA_Set_Uniform_Boundary_Conditions

















!+201+######################################################################################!
!                                                                                           !
!       Poseidon_Readiness_Check                                                            !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Readiness_Check(Readiness_Flag)


LOGICAL,                                INTENT(INOUT)       :: Readiness_Flag
LOGICAL                                                     :: Error_Flag


Error_Flag = .FALSE.

                                    !               !
                                    !   Mesh Check  !
                                    !               !
IF ( .NOT. RADIAL_MESH_SET_FLAG ) THEN

    PRINT*,"ERROR IN POSEIDON : The radial mesh was not set before 'Poseidon_Run' was called."
    Error_Flag = .TRUE.

END IF

IF ( .NOT. THETA_MESH_SET_FLAG ) THEN

    PRINT*,"ERROR IN POSEIDON : The theta mesh was not set before 'Poseidon_Run' was called."
    Error_Flag = .TRUE.

END IF

IF ( .NOT. PHI_MESH_SET_FLAG ) THEN

    PRINT*,"ERROR IN POSEIDON : The phi mesh was not set before 'Poseidon_Run' was called."
    Error_Flag = .TRUE.

END IF





                            !                               !
                            !   Boundary Conditions Check   !
                            !                               !
IF ( .NOT. INNER_BC_SET_FLAG ) THEN

    PRINT*,"ERROR IN POSEIDON : The inner boundary condition was not set before 'Poseidon_Run' was called."
    Error_Flag = .TRUE.

END IF


IF ( .NOT. OUTER_BC_SET_FLAG ) THEN

    PRINT*,"ERROR IN POSEIDON : Ther outer boundary condition was not set before 'Poseidon_Run' was called."
    Error_Flag = .TRUE.

END IF



Readiness_Flag = .NOT. Error_Flag





END SUBROUTINE Poseidon_Readiness_Check








!+501+##########################################################################!
!                                                                               !
!                   OUTPUT_ITERATION_REPORT                                     !
!                                                                               !
!###############################################################################!
SUBROUTINE OUTPUT_RUN_REPORT(Iter, Rank)

INTEGER, INTENT(IN)                 :: Iter, Rank

INTEGER                                         ::  FILE_ID
INTEGER                                         ::  i
REAL(KIND = idp)                                ::  r, theta, phi, deltar
REAL(KIND = idp)                                ::  Analytic_Val, Solver_Val
REAL(KIND = idp)                                ::  Return_Psi, Return_AlphaPsi
REAL(KIND = idp)                                ::  Return_Beta1, Return_Beta2, Return_Beta3
REAL(KIND = idp)                                ::  PsiPot_Val, AlphaPsiPot_Val

120 FORMAT (A61)
121 FORMAT (A1)
122 FORMAT (A41,I2.2)
123 FORMAT (A38,ES22.15)

109 FORMAT (A,I2.2,A,I2.2)
110 FORMAT (11X,A1,18X,A13,10X,A18,10X,A11,14X,A11,14X,A11)
111 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15)
112 FORMAT (A43,I2.2,A2,I2.2,A4)

FILE_ID = RUN_REPORT_FILE_ID



! Write Title to File
IF (( WRITE_REPORT_FLAG == 2) .OR. (WRITE_REPORT_FLAG == 3) ) THEN

    WRITE(FILE_ID,'(A)')"                                Average Timing Results"
    WRITE(FILE_ID,'(A)')"            ============================================================="
    WRITE(FILE_ID,'(A)')" "
    WRITE(FILE_ID,123)"                    Initialize Time : ",RUN_TIME_TABLE(1)
    WRITE(FILE_ID,123)" Input/Communicate Source Data Time : ",RUN_TIME_TABLE(2)
    WRITE(FILE_ID,123)"     Input Boundary Conditions Time : ",RUN_TIME_TABLE(3)
    WRITE(FILE_ID,123)"        CFA_3D_Apply_BCs_Part1 Time : ",RUN_TIME_TABLE(4)
    WRITE(FILE_ID,120)"-------------------------------------------------------------"
    WRITE(FILE_ID,123)" ||     Calc_3D_Current_Values Time : ",RUN_TIME_TABLE(5)
    WRITE(FILE_ID,123)" ||    CREATE_3D_SubJcbn_Terms Time : ",RUN_TIME_TABLE(6)
    WRITE(FILE_ID,123)" ||       CREATE_3D_RHS_VECTOR Time : ",RUN_TIME_TABLE(7)
    WRITE(FILE_ID,123)"\  /     CREATE_3D_JCBN_MATRIX Time : ",RUN_TIME_TABLE(8)
    WRITE(FILE_ID,120)"-\/ ---------------------------------------------------------"
    WRITE(FILE_ID,123)"CREATE_3D_NONLAPLACIAN_STF_MAT Time : ",RUN_TIME_TABLE(9)
    WRITE(FILE_ID,123)"REDUCE_3D_NONLAPLACIAN_STF_MAT Time : ",RUN_TIME_TABLE(10)
    WRITE(FILE_ID,123)"FINISH_3D_NONLAPLACIAN_STF_MAT Time : ",RUN_TIME_TABLE(11)
    WRITE(FILE_ID,123)"          FINISH_3D_RHS_VECTOR Time : ",RUN_TIME_TABLE(12)
    WRITE(FILE_ID,123)"        CFA_3D_Apply_BCs_Part2 Time : ",RUN_TIME_TABLE(13)
    WRITE(FILE_ID,123)"                    CFA_Solver Time : ",RUN_TIME_TABLE(14)
    WRITE(FILE_ID,123)"        CFA_Coefficient_Update Time : ",RUN_TIME_TABLE(15)
    WRITE(FILE_ID,123)"   CFA_Coefficient_Share_PETSc Time : ",RUN_TIME_TABLE(16)
    WRITE(FILE_ID,123)"         CFA_Convergence_Check Time : ",RUN_TIME_TABLE(17)
    WRITE(FILE_ID,123)"               Total Iteration Time : ",RUN_TIME_TABLE(18)
    WRITE(FILE_ID,123)"             Poseidon_Dist_Sol Time : ",RUN_TIME_TABLE(19)
    WRITE(FILE_ID,120)"============================================================="
    WRITE(FILE_ID,121)" "
    WRITE(FILE_ID,121)" "
    WRITE(FILE_ID,121)" "
    WRITE(FILE_ID,121)" "

END IF









WRITE(FILE_ID,'(A)')"                                 Convergence Results"
WRITE(FILE_ID,'(A)')"            ============================================================="
WRITE(FILE_ID,'(A)')""
IF ( Cur_Iteration == 2 ) THEN
    WRITE(FILE_ID,'(A,I2.2,A)')"The Newton-Raphson solver exited after ",Cur_Iteration-1," Iteration."
ELSE
    WRITE(FILE_ID,'(A,I2.2,A)')"The Newton-Raphson solver exited after ",Cur_Iteration-1," Iterations."
END IF


IF ( CONVERGENCE_FLAG == 1 ) THEN

    WRITE(FILE_ID,'(A)')"The solution converged within the tolerance set."

ELSE IF ( CONVERGENCE_FLAG == 2 ) THEN

    WRITE(FILE_ID,'(A)')"The solution did not converge within the maximum number of iterations allowed."
    WRITE(FILE_ID,'(A,I2.2)')"The maximum number of iterations allowed is ",Max_Iterations

END IF
WRITE(FILE_ID,'(A)')" "
WRITE(FILE_ID,'(A)')" "





















! Write Results Table Header to Screen
IF (( WRITE_REPORT_FLAG == 1) .OR. (WRITE_REPORT_FLAG == 3) ) THEN
    IF ( Rank == 0 ) THEN
        PRINT*,"++++++++++++++++++++++++++ myID,",Rank,"Sample Run Results ++++++++++++++++++++++++++"
        WRITE(*,110)"r","Psi Potential","AlphaPsi Potential","Beta Value1","Beta Value2","Beta Value3"
    END IF
END IF

! Write Results Table Header to File
IF (( WRITE_REPORT_FLAG == 2) .OR. (WRITE_REPORT_FLAG == 3) ) THEN
    WRITE(FILE_ID,'(A)')"                                 Sample of Final Results"
    WRITE(FILE_ID,'(A)')"            ============================================================="
    WRITE(FILE_ID,'(A)')" "
    WRITE(FILE_ID,110)"r","Psi Potential","AlphaPsi Potential","Beta Value1","Beta Value2","Beta Value3"
END IF


deltar = ( R_OUTER - R_INNER )/ REAL(ITER_REPORT_NUM_SAMPLES, KIND = idp)
DO i = 0,ITER_REPORT_NUM_SAMPLES

    r = i*deltar + R_INNER
    theta = pi/2.0_idp
    phi = pi/2.0_idp


    CALL Calc_3D_Values_At_Location( r, theta, phi,                              &
                                    Return_Psi, Return_AlphaPsi,                &
                                    Return_Beta1, Return_Beta2, Return_Beta3    )




    ! AlphaPsi_to_Pot   =   2*C_Square*(AlphaPsi - 1)
    ! Psi_to_Pot        =   2*C_Square*(1 - Psi)

    ! Calculate Conformal Factor value from Newtonian Potential
    PsiPot_Val = 2.0_idp*C_Square*(1.0_idp - Return_Psi)

    ! Calculate the product of the Conformal Factor and Lapse Function from Newtonian Potential
    AlphaPsiPot_Val = 2.0_idp*C_Square*(Return_AlphaPsi - 1.0_idp)


    ! Write Results to Screen
    IF (( WRITE_REPORT_FLAG == 1) .OR. (WRITE_REPORT_FLAG == 3) ) THEN
        IF ( Rank == 0 ) THEN
            WRITE(*,111) r,PsiPot_Val,AlphaPsiPot_Val,Return_Beta1,Return_Beta2,Return_Beta3
        END IF
    END IF

    ! Write Results to File
    IF (( WRITE_REPORT_FLAG == 2) .OR. (WRITE_REPORT_FLAG == 3) ) THEN
        WRITE(FILE_ID,111) r,PsiPot_Val,AlphaPsiPot_Val,Return_Beta1,Return_Beta2,Return_Beta3
    END IF

END DO






END SUBROUTINE OUTPUT_RUN_REPORT








END MODULE Poseidon_Main_Module
