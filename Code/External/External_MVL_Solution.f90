   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE External_MVL_Solution_Module                                          !##!
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
            ONLY : idp

USE Poseidon_Message_Routines_Module, &
            ONLY :  Warning_Message

USE Poseidon_Units_Module, &
            ONLY :  Centimeter

USE Variables_Mesh, &
            ONLY :  R_Inner,        &
                    R_Outer

USE Variables_External, &
            ONLY :  MVL_Boundary_Value,     &
                    MVL_C1,                 &
                    MVL_C2

IMPLICIT NONE


CONTAINS



 !+101+####################################################!
!                                                           !
!          Set_MVL_Test_Params                              !
!                                                           !
 !#########################################################!
SUBROUTINE Set_MVL_Test_Params( Boundary_Value_In )

REAL(idp),  INTENT(IN)          ::  Boundary_Value_In

IF ( R_Inner .LE. 0.0_idp) THEN

    CALL Warning_Message("Inner boundary location is incompatible with this test. Must be greater that zero.")
    STOP

ELSE

    MVL_Boundary_Value = Boundary_Value_In

    MVL_C1 = (2.0_idp*R_Outer*R_Outer*MVL_Boundary_Value/(R_Inner*R_Inner*R_Inner)) &
            /(1.0_idp + 2.0_idp*R_Outer*R_Outer*R_Outer/(R_Inner*R_Inner*R_Inner) )
    MVL_C2 = (R_Outer*R_Outer*MVL_Boundary_Value)       &
            /(1.0_idp + 2.0_idp*R_Outer*R_Outer*R_Outer/(R_Inner*R_Inner*R_Inner) )

    PRINT*,"MVL_C1"
    PRINT*,MVL_C1
    PRINT*,"MVL_C2"
    PRINT*,MVL_C2
    PRINT*,"Centimeter",Centimeter
END IF

END SUBROUTINE Set_MVL_Test_Params







 !+201+####################################################!
!                                                           !
!        MVL_Solution                                       !
!                                                           !
 !#########################################################!
FUNCTION MVL_Solution(r)


REAL(idp), INTENT(IN)       ::  r

REAL(idp)                   ::  r_wUnits
REAL(idp)                   ::  MVL_Solution


r_wUnits = r/Centimeter

MVL_Solution = MVL_C1*r_wUnits + MVL_C2/(r_wUnits*r_wUnits)

END FUNCTION MVL_Solution






END MODULE External_MVL_Solution_Module

