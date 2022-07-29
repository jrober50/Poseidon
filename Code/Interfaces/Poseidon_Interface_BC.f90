   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Poseidon_Interface_Boundary_Conditions                                !##!
!##!                                                                         !##!
!##!_________________________________________________________________________!##!
!##!                                                                         !##!
!##!                                                                         !##!
!##!=========================================================================!##!
!##!                                                                         !##!
!##!    Contains:                                                            !##!
!##!                                                                         !##!
!##!                                                                         !##!
!###############################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !##########################################################################!


!*D*================================!
!                                   !
!           Dependencies            !
!                                   !
!===================================!
USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Poseidon_Parameters, &
            ONLY :  Num_Vars

USE Variables_BC, &
            ONLY :  INNER_CFA_BC_VALUES,            &
                    OUTER_CFA_BC_VALUES,            &
                    INNER_CFA_BC_TYPE,              &
                    OUTER_CFA_BC_TYPE,              &
                    INNER_Poisson_BC_VALUE,         &
                    OUTER_Poisson_BC_VALUE,         &
                    INNER_Poisson_BC_TYPE,          &
                    OUTER_Poisson_BC_TYPE

IMPLICIT NONE



INTERFACE Poseidon_Set_Uniform_Boundary_Conditions
    MODULE PROCEDURE Poseidon_Set_Uniform_BC_XCFC
    MODULE PROCEDURE Poseidon_Set_Uniform_BC_Poisson
END INTERFACE Poseidon_Set_Uniform_Boundary_Conditions


CONTAINS




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
 SUBROUTINE Poseidon_Set_Uniform_BC_XCFC(   BC_Location_Input,      &
                                            BC_Type_Input,          &
                                            BC_Value_Input          )



 CHARACTER(LEN = 1),                                INTENT(IN)      ::  BC_Location_Input
 CHARACTER(LEN = 1),    DIMENSION(1:Num_Vars),  INTENT(IN)      ::  BC_Type_Input

 REAL(idp),             DIMENSION(1:Num_Vars),  INTENT(IN)      ::  BC_Value_Input



 IF (    BC_Location_Input == "I"    ) THEN

     INNER_CFA_BC_TYPE(1:Num_Vars) = BC_Type_Input(1:Num_Vars)
     INNER_CFA_BC_VALUES(1:Num_Vars) = BC_Value_Input(1:Num_Vars)

 ELSE IF (    BC_Location_Input == "O"    ) THEN

     OUTER_CFA_BC_TYPE(1:Num_Vars) = BC_Type_Input(1:Num_Vars)
     OUTER_CFA_BC_VALUES(1:Num_Vars) = BC_Value_Input(1:Num_Vars)

 END IF


 END SUBROUTINE Poseidon_Set_Uniform_BC_XCFC









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
 SUBROUTINE Poseidon_Set_Uniform_BC_Poisson( BC_Location_Input,      &
                                             BC_Type_Input,          &
                                             BC_Value_Input          )



 CHARACTER(LEN = 1), INTENT(IN)          ::  BC_Location_Input
 CHARACTER(LEN = 1), INTENT(IN)          ::  BC_Type_Input
 REAL(KIND = idp),   INTENT(IN)          ::  BC_Value_Input


 IF (    BC_Location_Input == "I"    ) THEN

    INNER_CFA_BC_TYPE(1)   = BC_Type_Input
    INNER_CFA_BC_VALUES(1) = BC_Value_Input

 ELSE IF (    BC_Location_Input == "O"    ) THEN

    OUTER_CFA_BC_TYPE(1)   = BC_Type_Input
    OUTER_CFA_BC_VALUES(1) = BC_Value_Input

 END IF


 END SUBROUTINE Poseidon_Set_Uniform_BC_Poisson









END MODULE Poseidon_Interface_Boundary_Conditions
