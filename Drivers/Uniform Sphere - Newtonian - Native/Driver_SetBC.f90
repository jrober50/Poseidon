   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Driver_SetBC_Module                                                   !##!
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

USE Poseidon_Numbers_Module, &
            ONLY :  pi

USE Poseidon_Units_Module, &
            ONLY :  Grav_Constant_G,        &
                    Centimeter,             &
                    Second

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag

USE Poseidon_Message_Routines_Module, &
            ONLY :  Driver_Init_Message

USE Variables_Mesh, &
            ONLY :  Num_R_Elements,         &
                    rlocs

USE External_UST_Solution_Module, &
            ONLY :  UST_Solution

USE Poseidon_Interface_Boundary_Conditions, &
            ONLY :  Poseidon_Set_Uniform_Boundary_Conditions

USE Variables_External, &
            ONLY :  UST_Rhoo,          &
                    UST_Star_Radius


IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!   Driver_SetBC                                                                !
!                                                                               !
!###############################################################################!
SUBROUTINE Driver_SetBC(Surface_RE,         &
                        rho_o               )

INTEGER,    INTENT(INOUT)                       ::  Surface_RE
REAL(idp),  INTENT(INOUT)                       ::  rho_o

CHARACTER(LEN=1)                                ::  Inner_BC_Type
CHARACTER(LEN=1)                                ::  Outer_BC_Type
REAL(idp)                                       ::  Inner_BC_Value
REAL(idp)                                       ::  Outer_BC_Value

IF ( Verbose_Flag ) CALL Driver_Init_Message('Calculating boundary conditions.')


Inner_BC_Type = "N"
Outer_BC_Type = "D"

Inner_BC_Value = 0.0_idp


IF ( Surface_RE .LE. Num_R_Elements-1 ) THEN
    ! Surface is Inside the Domain  !
    Outer_BC_Value = -2.0_idp*pi*Grav_Constant_G        &
                   * rho_o                              &
                   * (3.0_idp*rlocs(Surface_RE)**2      &
                        - rlocs(Num_R_Elements)**2  )   &
                   / 3.0_idp

ELSE
    ! Surface is Outside the Domain !
    Outer_BC_Value = -4.0_idp*pi*Grav_Constant_G         &
                   * rho_o                              &
                   * rlocs(Surface_RE)**3               &
                   /(3.0_idp*rlocs(Num_R_Elements))

END IF


Outer_BC_Value = UST_Solution(rlocs(Num_R_Elements))


CALL Poseidon_Set_Uniform_Boundary_Conditions("I", Inner_BC_Type, Inner_BC_Value)
CALL Poseidon_Set_Uniform_Boundary_Conditions("O", Outer_BC_Type, Outer_BC_Value)




END SUBROUTINE Driver_SetBC




END MODULE Driver_SetBC_Module


