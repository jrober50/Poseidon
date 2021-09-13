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

INTEGER, PUBLIC, PARAMETER          ::  iU_CF = 1 ! Conformal Factor
INTEGER, PUBLIC, PARAMETER          ::  iU_LF = 2 ! Lapse Function
INTEGER, PUBLIC, PARAMETER          ::  iU_S1 = 3 ! Shift Vector, Radial Component
INTEGER, PUBLIC, PARAMETER          ::  iU_S2 = 4 ! Shift Vector, Theta Component
INTEGER, PUBLIC, PARAMETER          ::  iU_S3 = 5 ! Shift Vector, Phi Component
INTEGER, PUBLIC, PARAMETER          ::  iU_X1 = 6 ! X Vector, Radial Component
INTEGER, PUBLIC, PARAMETER          ::  iU_X2 = 7 ! X Vector, Theta Component
INTEGER, PUBLIC, PARAMETER          ::  iU_X3 = 8 ! X Vector, Phi Component

INTEGER, PUBLIC, PARAMETER          ::  iVB_S = 1
INTEGER, PUBLIC, PARAMETER          ::  iVB_X = 2


INTEGER, PARAMETER                  ::  iS_E  = 1
INTEGER, PARAMETER                  ::  iS_S  = 2
INTEGER, PARAMETER                  ::  iS_S1 = 3
INTEGER, PARAMETER                  ::  iS_S2 = 4
INTEGER, PARAMETER                  ::  iS_S3 = 5

END MODULE Parameters_Variable_Indices
