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
            ONLY :  idp, pi

USE Units_Module, &
            ONLY :  Set_Units

USE Poseidon_Parameters, &
            ONLY :  DEGREE,                     &
                    L_LIMIT,                    &
                    NUM_CFA_VARS,               &
                    POSEIDON_INITIALIZED_FLAG,  &
                    NUM_R_QUAD_POINTS,          &
                    NUM_T_QUAD_POINTS,          &
                    NUM_P_QUAD_POINTS,          &
                    Poseidon_Frame,             &
                    SOLVER_TYPE_FLAG


USE Poseidon_Variables_Module, &
            ONLY :  R_INNER, R_OUTER,                               &
                    NUM_R_ELEMENTS, NUM_T_ELEMENTS, NUM_P_ELEMENTS, &
                    rlocs, tlocs, plocs,                            &
                    Test_Space_Allocated_Flag,                      &
                    Stiffness_Matrix_Initialized_Flag,              &
                    FirstCall_Flag,                                 &
                    INNER_BC_SET_FLAG,                              &
                    OUTER_BC_SET_FLAG,                              &
                    RADIAL_MESH_SET_FLAG,                           &
                    THETA_MESH_SET_FLAG,                            &
                    PHI_MESH_SET_FLAG,                              &
                    PHYSICS_TYPE,                                   &
                    VAR_DIM,                                        &
                    PROB_DIM,                                       &
                    INNER_CFA_BC_VALUES,                            &
                    OUTER_CFA_BC_VALUES,                            &
                    INNER_CFA_BC_TYPE,                              &
                    OUTER_CFA_BC_TYPE,                              &
                    LM_LENGTH,                                      &
                    Block_Source_E,             &
                    Block_Source_S,             &
                    Block_Source_Si,            &
                    INT_R_LOCATIONS,            &
                    INT_T_LOCATIONS,            &
                    INT_P_LOCATIONS

USE Poseidon_Initialization_Module,                                 &
            ONLY :  Poseidon_Initialize_1D,                         &
                    Poseidon_Initialize_3D


USE Poseidon_Mesh_Module, &
            ONLY :  Generate_Defined_Mesh
                    

USE Poseidon_Allocation_Module, &
            ONLY :  Deallocate_Poseidon_CFA_Variables


USE CFA_Newton_Raphson_3D_Module, &
            ONLY :  CFA_Newton_Raphson_3D


USE CFA_GMRES_Module, &
            ONLY :  CFA_GMRES

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
!       Poseidon_Initialize                                                                 !
!                                                                                           !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Initialize( Units,                                                                  &
                                Dimensions,                                                             &
                                FEM_Degree_Input, L_Limit_Input,                                        &
                                Inner_Radius, Outer_Radius,                                             &
                                R_Elements_Input, T_Elements_Input, P_Elements_Input,                   &
                                Local_R_Elements_Input, Local_T_Elements_Input, Local_P_Elements_Input, &
                                Num_R_Quad_Input, Num_T_Quad_Input, Num_P_Quad_Input,                   &
                                Input_Delta_R_Vector, Input_Delta_T_Vector, Input_Delta_P_Vector        )



                                         !                          !
                                        !!      Input Variables     !!
                                         !

INTEGER, INTENT(IN)                                                             ::  Dimensions,             &
                                                                                    FEM_Degree_Input,       &
                                                                                    L_Limit_Input,          &
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
CHARACTER(LEN = 1), INTENT(IN)                                                  ::  Units

REAL(KIND = idp), DIMENSION(1:R_Elements_Input),    OPTIONAL,   INTENT(IN)      ::  Input_Delta_R_Vector
REAL(KIND = idp), DIMENSION(1:T_Elements_Input),    OPTIONAL,   INTENT(IN)      ::  Input_Delta_T_Vector
REAL(KIND = idp), DIMENSION(1:P_Elements_Input),    OPTIONAL,   INTENT(IN)      ::  Input_Delta_P_Vector



CALL Set_Units(Units)

IF ( POSEIDON_INITIALIZED_FLAG .EQV. .FALSE. ) THEN


    IF ( Dimensions == 1 ) THEN
        CALL Poseidon_Initialize_1D( FEM_Degree_Input, L_Limit_Input,                                        &
                                     Inner_Radius, Outer_Radius,                                             &
                                     R_Elements_Input, T_Elements_Input, P_Elements_Input,                   &
                                     Local_R_Elements_Input, Local_T_Elements_Input, Local_P_Elements_Input, &
                                     Num_R_Quad_Input, Num_T_Quad_Input, Num_P_Quad_Input,                   &
                                     Input_Delta_R_Vector, Input_Delta_T_Vector, Input_Delta_P_Vector        )

    ELSE

        CALL Poseidon_Initialize_3D( FEM_Degree_Input, L_Limit_Input,                                        &
                                     Inner_Radius, Outer_Radius,                                             &
                                     R_Elements_Input, T_Elements_Input, P_Elements_Input,                   &
                                     Local_R_Elements_Input, Local_T_Elements_Input, Local_P_Elements_Input, &
                                     Num_R_Quad_Input, Num_T_Quad_Input, Num_P_Quad_Input,                   &
                                     Input_Delta_R_Vector, Input_Delta_T_Vector, Input_Delta_P_Vector        )

    END IF

    POSEIDON_INITIALIZED_FLAG = .TRUE.

END IF

END SUBROUTINE Poseidon_Initialize











!+102+#######################################################################################!
!                                                                                           !
!       Poseidon_Run - Calculates Source Vector and Solves for solution coefficients        !
!                                                                                           !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Run()

LOGICAL                                             ::  Readiness_Flag
INTEGER, DIMENSION(1:5)                             ::  Eq_Flag_Array



PRINT*,"Quad",NUM_R_QUAD_POINTS,NUM_T_QUAD_POINTS,NUM_P_QUAD_POINTS

!CALL Poseidon_Readiness_Check(Readiness_Flag)
Readiness_Flag = .TRUE.

IF ( Readiness_Flag ) THEN
    IF ( SOLVER_TYPE_FLAG == 1 ) THEN

        CALL CFA_Newton_Raphson_3D()
    
    ELSE IF ( SOLVER_TYPE_FLAG == 2 ) THEN

        Eq_Flag_Array = [ 1, 1, 1, 1, 1]
        CALL CFA_GMRES(Eq_Flag_Array)

    ELSE

        PRINT*,"ERROR IN POSEIDON : Solver Type Flag has invalid value. "

    END IF

ELSE

    PRINT*, "ERROR IN POSEIDON : There was an error in setting up Poseidon, therefore it did not run."

END IF
Poseidon_Frame = Poseidon_Frame + 1




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
