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
USE Poseidon_Kinds_Module, &
                    ONLY : idp

USE Poseidon_Numbers_Module, &
                    ONLY : pi


USE Poseidon_Parameters, &
            ONLY :  Num_CFA_Vars,               &
                    Poseidon_Frame,             &
                    Poisson_Mode,               &
                    Method_Flag


USE Variables_Mesh, &
            ONLY :  R_Inner,                    &
                    Num_R_Elements,             &
                    Num_T_Elements,             &
                    Num_P_Elements,             &
                    rlocs,                      &
                    tlocs,                      &
                    plocs,                      &
                    RADIAL_MESH_SET_FLAG,       &
                    THETA_MESH_SET_FLAG,        &
                    PHI_MESH_SET_FLAG
            
USE Variables_BC, &
            ONLY :  INNER_CFA_BC_VALUES,            &
                    OUTER_CFA_BC_VALUES,            &
                    INNER_CFA_BC_TYPE,              &
                    OUTER_CFA_BC_TYPE,              &
                    INNER_Poisson_BC_VALUE,         &
                    OUTER_Poisson_BC_VALUE,         &
                    INNER_Poisson_BC_TYPE,          &
                    OUTER_Poisson_BC_TYPE

USE Variables_Interface, &
            ONLY :  Caller_Set,                     &
                    Caller_RQ_xlocs,                &
                    Caller_TQ_xlocs,                &
                    Caller_PQ_xlocs,                &
                    Translation_Matrix

USE Functions_Mesh, &
            ONLY :  Generate_Defined_Mesh
                    

USE Allocation_Sources, &
            ONLY :  Deallocate_Poseidon_Source_Variables

USE Allocation_Poisson, &
            ONLY :  Deallocate_Poseidon_Poisson_Variables


USE XCFC_Source_Routine_Variables_Module, &
            ONLY :  Deallocate_XCFC_Source_Routine_Variables

USE Poisson_Main_Module, &
            ONLY :  Poisson_Solve

USE FP_AndersonM_Module, &
            ONLY : Fixed_Point_AndersonM

USE XCFC_Method_Module, &
            ONLY : XCFC_Method

USE Allocation_Mesh, &
            ONLY : Deallocate_Mesh

USE Allocation_Quadrature, &
            ONLY : Deallocate_Quadrature

USE Allocation_Tables, &
            ONLY : Deallocate_Tables

!USE Allocation_NR, &
!            ONLY : Deallocate_NR

USE Allocation_FP, &
            ONLY : Deallocate_FP

USE Allocation_XCFC_Linear_Systems, &
            ONLY : Deallocate_XCFC_Linear_Systems

USE Allocation_SelfSimilar, &
            ONLY : Deallocate_SelfSim

USE Timer_Routines_Module, &
            ONLY : Finalize_Timers

USE IO_Print_Results, &
            ONLY :  Print_Results

USE Flags_Main_Module, &
            ONLY : Poseidon_Clear_All_Flags

USE Flags_Boundary_Conditions_Module, &
            ONLY :  lPF_BC_Flags,           &
                    iPF_BC_Outer_Set,       &
                    iPF_BC_Inner_Set

USE Flags_IO_Module, &
            ONLY :  lPF_IO_Flags,           &
                    iPF_IO_Print_Results

USE Flags_Check_Routines, &
            ONLY :  Poseidon_Run_Check

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


!Readiness_Flag = .TRUE.
Readiness_Flag = Poseidon_Run_Check()


IF ( Readiness_Flag ) THEN

    IF ( Poisson_Mode .eqv. .TRUE. ) THEN

        CALL Poisson_Solve()
    
    ELSE IF ( Method_Flag == 1 ) THEN
        WRITE(*,'(A)')"The Newton-Raphson method is not currently available in Poseidon. STOPING"
        STOP

!        CALL CFA_Newton_Raphson_3D()

    ELSE IF ( Method_Flag == 2 ) THEN

        Call Fixed_Point_AndersonM()

    ELSE IF ( Method_Flag == 3 ) THEN
        CALL XCFC_Method()
    ELSE

        PRINT*,"ERROR IN POSEIDON : Solver Type Flag has invalid value. "
        STOP

    END IF

ELSE

    PRINT*, "ERROR IN POSEIDON : There was an error in setting up Poseidon, therefore it did not run."
    STOP

END IF
Poseidon_Frame = Poseidon_Frame + 1



IF ( lPF_IO_Flags(iPF_IO_Print_Results) ) THEN
    Call Print_Results()
END IF



END SUBROUTINE Poseidon_Run
















!+103+######################################################################################!
!                                                                                           !
!       Poseidon_Close - Deallocate all Poseidon variables.                                 !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Close()


CALL Deallocate_Mesh()

IF ( Caller_Set ) THEN
    DEALLOCATE( Caller_RQ_xlocs )
    DEALLOCATE( Caller_TQ_xlocs )
    DEALLOCATE( Caller_PQ_xlocs )
    DEALLOCATE( Translation_Matrix )
    Caller_Set = .FALSE.
END IF


IF ( Poisson_Mode ) THEN
    Call Deallocate_Poseidon_Poisson_Variables

ELSE
    !!!!  Deallocate Data Space !!!!
    CALL Deallocate_Poseidon_Source_Variables
    Call Deallocate_XCFC_Source_Routine_Variables

    CALL Deallocate_Quadrature()
    CALL Deallocate_Tables()

    CALL Deallocate_SelfSim()

    IF ( Method_Flag == 1 ) THEN
        WRITE(*,'(A)')"The Newton-Raphson method is not currently available in Poseidon. STOPING"
        STOP
!        CALL Deallocate_NR
    ELSE IF ( Method_Flag == 2 ) THEN
        CALL Deallocate_FP
    ELSE IF ( Method_Flag == 3 ) THEN
        CALL Deallocate_XCFC_Linear_Systems
    END IF

END IF


CALL Finalize_Timers()
CALL Poseidon_Clear_All_Flags()


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

REAL(KIND = idp), DIMENSION(1:NUM_T_ELEMENTS)                                  ::  TMP_Delta_T_Vector
REAL(KIND = idp), DIMENSION(1:NUM_P_ELEMENTS)                                  ::  TMP_Delta_P_Vector



  !                                                                                                 !
 !!     If Input_Delta_R_Vector is present then user wishes to define a non-uniform radial mesh.    !!
  !                                                                                                 !
IF (    PRESENT(Input_Delta_R_Vector)   ) THEN




    CALL Generate_Defined_Mesh(NUM_R_ELEMENTS, R_INNER, Input_Delta_R_Vector, rlocs)
    RADIAL_MESH_SET_FLAG = .TRUE.




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



IF (    BC_Location_Input == "I"    ) THEN

    INNER_CFA_BC_TYPE(1:NUM_CFA_VARS) = BC_Type_Input(1:NUM_CFA_VARS)
    INNER_CFA_BC_VALUES(1:NUM_CFA_VARS) = BC_Value_Input(1:NUM_CFA_VARS)




ELSE IF (    BC_Location_Input == "O"    ) THEN

    OUTER_CFA_BC_TYPE(1:NUM_CFA_VARS) = BC_Type_Input(1:NUM_CFA_VARS)
    OUTER_CFA_BC_VALUES(1:NUM_CFA_VARS) = BC_Value_Input(1:NUM_CFA_VARS)


END IF



END SUBROUTINE Poseidon_CFA_Set_Uniform_Boundary_Conditions









 !+105+####################################################################################!
!                                                                                           !
!      Poseidon_Poisson_Set_Boundary_Condtion                                               !
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
SUBROUTINE Poseidon_Poisson_Set_Uniform_Boundary_Conditions(BC_Location_Input,      &
                                                            BC_Type_Input,          &
                                                            BC_Value_Input          )



CHARACTER(LEN = 1), INTENT(IN)          ::  BC_Location_Input
CHARACTER(LEN = 1), INTENT(IN)          ::  BC_Type_Input
REAL(KIND = idp),   INTENT(IN)          ::  BC_Value_Input



IF (    BC_Location_Input == "I"    ) THEN

    INNER_Poisson_BC_TYPE  = BC_Type_Input
    INNER_Poisson_BC_VALUE = BC_Value_Input

ELSE IF (    BC_Location_Input == "O"    ) THEN

    OUTER_Poisson_BC_TYPE  = BC_Type_Input
    OUTER_Poisson_BC_VALUE = BC_Value_Input


END IF







END SUBROUTINE Poseidon_Poisson_Set_Uniform_Boundary_Conditions








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
IF ( .NOT. lPF_BC_Flags(iPF_BC_Inner_Set) ) THEN

    PRINT*,"ERROR IN POSEIDON : The inner boundary condition was not set before 'Poseidon_Run' was called."
    Error_Flag = .TRUE.

END IF


IF ( .NOT. lPF_BC_Flags(iPF_BC_Outer_Set) ) THEN

    PRINT*,"ERROR IN POSEIDON : Ther outer boundary condition was not set before 'Poseidon_Run' was called."
    Error_Flag = .TRUE.

END IF



Readiness_Flag = .NOT. Error_Flag





END SUBROUTINE Poseidon_Readiness_Check
















END MODULE Poseidon_Main_Module
