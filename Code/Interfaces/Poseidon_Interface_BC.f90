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
            ONLY :  Inner_BC_Values,            &
                    Outer_BC_Values,            &
                    Inner_BC_Type,              &
                    Outer_BC_Type

IMPLICIT NONE



INTERFACE Poseidon_Set_Uniform_Boundary_Conditions
    MODULE PROCEDURE Poseidon_Set_Uniform_BC_XCFC
    MODULE PROCEDURE Poseidon_Set_Uniform_BC_Newtonian
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

     Inner_BC_Type(1:Num_Vars) = BC_Type_Input(1:Num_Vars)
     Inner_BC_Values(1:Num_Vars) = BC_Value_Input(1:Num_Vars)

 ELSE IF (    BC_Location_Input == "O"    ) THEN

     Outer_BC_Type(1:Num_Vars) = BC_Type_Input(1:Num_Vars)
     Outer_BC_Values(1:Num_Vars) = BC_Value_Input(1:Num_Vars)

 END IF


 END SUBROUTINE Poseidon_Set_Uniform_BC_XCFC









  !+105+####################################################################################!
 !                                                                                           !
 !      Poseidon_Newtonian_Set_Boundary_Condtion                                               !
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
 SUBROUTINE Poseidon_Set_Uniform_BC_Newtonian( BC_Location_Input,      &
                                               BC_Type_Input,          &
                                               BC_Value_Input          )



 CHARACTER(LEN = 1), INTENT(IN)          ::  BC_Location_Input
 CHARACTER(LEN = 1), INTENT(IN)          ::  BC_Type_Input
 REAL(KIND = idp),   INTENT(IN)          ::  BC_Value_Input


 IF (    BC_Location_Input == "I"    ) THEN

    Inner_BC_Type   = BC_Type_Input
    Inner_BC_Values = BC_Value_Input

 ELSE IF (    BC_Location_Input == "O"    ) THEN

    Outer_BC_Type   = BC_Type_Input
    Outer_BC_Values = BC_Value_Input

 END IF


 END SUBROUTINE Poseidon_Set_Uniform_BC_Newtonian









END MODULE Poseidon_Interface_Boundary_Conditions
