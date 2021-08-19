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

INTEGER, PUBLIC, PARAMETER          ::  iCF = 1 ! Conformal Factor
INTEGER, PUBLIC, PARAMETER          ::  iLF = 2 ! Lapse Function
INTEGER, PUBLIC, PARAMETER          ::  iS1 = 3 ! Shift Vector, Radial Component
INTEGER, PUBLIC, PARAMETER          ::  iS2 = 4 ! Shift Vector, Theta Component
INTEGER, PUBLIC, PARAMETER          ::  iS3 = 5 ! Shift Vector, Phi Component
INTEGER, PUBLIC, PARAMETER          ::  iX1 = 6 ! X Vector, Radial Component
INTEGER, PUBLIC, PARAMETER          ::  iX2 = 7 ! X Vector, Theta Component
INTEGER, PUBLIC, PARAMETER          ::  iX3 = 8 ! X Vector, Phi Component


END MODULE Parameters_Variable_Indices
