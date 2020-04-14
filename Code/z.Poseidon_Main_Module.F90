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



USE Poseidon_Variables_Module, &
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
                    LM_LENGTH,                                      &
                    ULM_LENGTH,                                     &
                    M_VALUES,                                       &
                    Matrix_Location,                                &
                    LM_Location,                                    &
                    POSEIDON_COMM_WORLD,                            &
                    RUN_TIME_TABLE






USE Poseidon_Additional_Functions_Module, &
            ONLY :  Spherical_Harmonic,                             &
                    Map_To_X_Space, Map_From_X_Space,               &
                    Initialize_LG_Quadrature,                       &
                    Initialize_LG_Quadrature_Locations,             &
                    Generate_Defined_Mesh,                          &
                    Generate_Defined_Coarse_Mesh


USE Allocate_Variables_Module, &
            ONLY :  Allocate_Poseidon_CFA_Variables,                &
                    Deallocate_Poseidon_CFA_Variables


USE Jacobian_Internal_Functions_Module,  &
            ONLY :  Initialize_Guess_Values,                        &
                    Initialize_Ylm_Tables,                          &
                    Initialize_Lagrange_Poly_Tables

USE CFA_Newton_Raphson_3D_Module, &
            ONLY :  CFA_Newton_Raphson_3D

USE Poseidon_MPI_Module, &
            ONLY :  CREATE_POSEIDON_COMMUNICATORS

USE Poseidon_Parameter_Read_Module, &
            ONLY :  UNPACK_POSEIDON_PARAMETERS

USE Poseidon_Additional_Functions_Module, &
                        ONLY :  Calc_3D_Values_At_Location


USE mpi



IMPLICIT NONE




                    !*F&S*==========================================!
                    !                                               !
                    !           Functions & Subroutines             !
                    !                                               !
                    !===============================================!
CONTAINS















!+102+#######################################################################################!
!                                                                                           !
!       Poseidon_Run - Calculates Source Vector and Solves for solution coefficients        !
!                                                                                           !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Run()

LOGICAL                                             ::  Readiness_Flag






!CALL Poseidon_Readiness_Check(Readiness_Flag)
Readiness_Flag = .TRUE.

IF ( Readiness_Flag ) THEN

    CALL CFA_Newton_Raphson_3D()
ELSE

    PRINT*, "ERROR IN POSEIDON : There was an error in setting up Poseidon, therefore it did not run."

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
















END MODULE Poseidon_Main_Module
