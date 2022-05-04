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

USE Poseidon_Units_Module, &
        ONLY :  C_Square

USE Variables_Functions, &
        ONLY :  Potential_Solution

USE Poseidon_Main_Module, &
        ONLY :  Poseidon_CFA_Set_Uniform_Boundary_Conditions


IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!   Driver_SetBC                                                                !
!                                                                               !
!###############################################################################!
SUBROUTINE Driver_SetBC( NE, x_e, i )

INTEGER,   DIMENSION(3),       INTENT(IN)       ::  NE
REAL(idp), DIMENSION(0:NE(1)), INTENT(IN)       ::  x_E
INTEGER,                       INTENT(IN)       ::  i

REAL(idp)                                       ::  Psi_BC
REAL(idp)                                       ::  AlphaPsi_BC
REAL(idp)                                       ::  Shift_Vector_BC
CHARACTER(LEN=1), DIMENSION(1:5)                ::  INNER_BC_TYPES, OUTER_BC_TYPES
REAL(idp), DIMENSION(1:5)                       ::  INNER_BC_VALUES, OUTER_BC_VALUES



Shift_Vector_BC = -1.0E2_idp


!WRITE(*,'(A,ES12.4E3)') "Shift Vector BC : ",Shift_Vector_BC




INNER_BC_TYPES = (/"N", "N","N","N","N"/)
OUTER_BC_TYPES = (/"D", "D","D","D","D"/)

INNER_BC_VALUES = (/0.0_idp, 0.0_idp, 0.0_idp, 0.0_idp, 0.0_idp /)
OUTER_BC_VALUES = (/1.0_idp, 1.0_idp, Shift_Vector_BC, 0.0_idp, 0.0_idp /)

CALL Poseidon_CFA_Set_Uniform_Boundary_Conditions("I", INNER_BC_TYPES, INNER_BC_VALUES)
CALL Poseidon_CFA_Set_Uniform_Boundary_Conditions("O", OUTER_BC_TYPES, OUTER_BC_VALUES)


END SUBROUTINE Driver_SetBC








FUNCTION Calc_C1(a,b,alpha,gamma)

REAL(idp), INTENT(IN)           ::  a
REAL(idp), INTENT(IN)           ::  b
REAL(idp), INTENT(IN)           ::  alpha
REAL(idp), INTENT(IN)           ::  gamma

REAL(idp)                       ::  Calc_C1

Calc_C1 = ( alpha + 2*b*b*gamma/(a*a*a))/(1+ 2*b*b*b/(a*a*a))

END FUNCTION Calc_C1



FUNCTION Calc_C2(a,b,alpha,gamma)

REAL(idp), INTENT(IN)           ::  a
REAL(idp), INTENT(IN)           ::  b
REAL(idp), INTENT(IN)           ::  alpha
REAL(idp), INTENT(IN)           ::  gamma

REAL(idp)                       ::  Calc_C2

Calc_C2 = ( b*b*gamma - alpha*b*b*b)/(1+ 2*b*b*b/(a*a*a))

END FUNCTION Calc_C2




FUNCTION MVL_Solution(r, C1, C2)

REAL(idp), INTENT(IN)           ::  r
REAL(idp), INTENT(IN)           ::  C1
REAL(idp), INTENT(IN)           ::  C2

REAL(idp)                       ::  MVL_Solution

MVL_Solution = C1*r + C2/(r*r)

END FUNCTION MVL_Solution









END MODULE Driver_SetBC_Module

