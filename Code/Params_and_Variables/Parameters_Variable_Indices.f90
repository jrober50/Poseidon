   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Parameters_Variable_Indices                                                 !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!

IMPLICIT NONE

!
!   Solution Variable, U, Indices
!
INTEGER, PUBLIC, PARAMETER          ::  iU_CF = 1 ! Conformal Factor
INTEGER, PUBLIC, PARAMETER          ::  iU_LF = 2 ! Lapse Function
INTEGER, PUBLIC, PARAMETER          ::  iU_S1 = 3 ! Shift Vector, Radial Component
INTEGER, PUBLIC, PARAMETER          ::  iU_S2 = 4 ! Shift Vector, Theta Component
INTEGER, PUBLIC, PARAMETER          ::  iU_S3 = 5 ! Shift Vector, Phi Component
INTEGER, PUBLIC, PARAMETER          ::  iU_X1 = 6 ! X Vector, Radial Component
INTEGER, PUBLIC, PARAMETER          ::  iU_X2 = 7 ! X Vector, Theta Component
INTEGER, PUBLIC, PARAMETER          ::  iU_X3 = 8 ! X Vector, Phi Component

!
!   Vector Type B, VB, Flag Indices
!
INTEGER, PUBLIC, PARAMETER          ::  iVB_S = 1
INTEGER, PUBLIC, PARAMETER          ::  iVB_X = 2


!
!   Source Varaible, S, Indices
!
INTEGER, PUBLIC, PARAMETER          ::  iS_E  = 1
INTEGER, PUBLIC, PARAMETER          ::  iS_S  = 5
INTEGER, PUBLIC, PARAMETER          ::  iS_S1 = 2
INTEGER, PUBLIC, PARAMETER          ::  iS_S2 = 3
INTEGER, PUBLIC, PARAMETER          ::  iS_S3 = 4




!
!   Quadrature Variable, Q, Indices
!
!INTEGER, PUBLIC, PARAMETER          ::  iQ_R   = 1
!INTEGER, PUBLIC, PARAMETER          ::  iQ_T   = 2
!INTEGER, PUBLIC, PARAMETER          ::  iQ_P   = 3
!INTEGER, PUBLIC, PARAMETER          ::  iQ_TP  = 4
!INTEGER, PUBLIC, PARAMETER          ::  iQ_All = 5


END MODULE Parameters_Variable_Indices
