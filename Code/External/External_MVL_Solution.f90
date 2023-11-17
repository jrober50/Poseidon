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

USE Variables_BC, &
            ONLY :  Outer_BC_Values,        &
                    Inner_BC_Values

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
SUBROUTINE Set_MVL_Test_Params( )

REAL(idp)               :: alpha
REAL(idp)               :: gamma

REAL(idp)               :: a
REAL(idp)               :: b

REAL(idp)               :: C1_Numer, C2_Numer
REAL(idp)               :: Common_Denom

IF ( R_Inner .LE. -1.0_idp) THEN

    CALL Warning_Message("Inner boundary location is incompatible with this test. Must be greater that zero.")
    STOP

ELSE

    MVL_Boundary_Value = Outer_BC_Values(3)
    a = R_Inner
    b = R_Outer
    
    alpha = Inner_BC_Values(3)
    gamma = Outer_BC_Values(3)

!    MVL_C1 = (2.0_idp*R_Outer*R_Outer*MVL_Boundary_Value/(R_Inner*R_Inner*R_Inner)) &
!            /(1.0_idp + 2.0_idp*R_Outer*R_Outer*R_Outer/(R_Inner*R_Inner*R_Inner) )
!    MVL_C2 = (R_Outer*R_Outer*MVL_Boundary_Value)       &
!            /(1.0_idp + 2.0_idp*R_Outer*R_Outer*R_Outer/(R_Inner*R_Inner*R_Inner) )

    PRINT*,"a    :",a
    print*,"b    :",b
    print*,"Alpha: ",alpha
    print*,"Gamma: ",gamma

    C1_Numer = alpha + 2.0_idp*gamma*b*b/(a*a*a)
    C2_Numer = b*b*gamma - alpha*b*b*b
    Common_Denom = 1 + 2*b*b*b/(a*a*a)
    
!    MVL_C1 = C1_Numer/Common_Denom
!    MVL_C2 = C2_Numer/Common_Denom

    MVL_C1 = -1.0E+2_idp
    MVL_C2 = -5.0E+1_idp

!    MVL_C1 = -1.0E+2_idp
!    MVL_C2 = 0.0_idp

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

