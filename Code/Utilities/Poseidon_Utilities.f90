   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Poseidon_Utilities_Module                                             !##!
!##!                                                                         !##!
!##!_________________________________________________________________________!##!
!##!                                                                         !##!
!##!                                                                         !##!
!##!=========================================================================!##!
!##!                                                                         !##!
!##!    Contains:                                                            !##!
!##!                                                                         !##!
!##!        +101+                   Poseidon_Calc_ADM_Mass                   !##!
!##!        +102+                   Poseidon_Calc_ADM_Mass_Parts             !##!
!##!                                                                         !##!
!##!        +201+                   Poseidon_Calc_Komar_Mass                 !##!
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

USE ADM_Mass_Module, &
            ONLY : Calc_ADM_Mass

USE ADM_Mass_In_Parts_Module, &
            ONLY : Calc_ADM_Mass_In_Parts

USE Komar_Mass_Module, &
            ONLY : Calc_Komar_Mass

USE Timer_Routines_Module, &
            ONLY :  TimerStart,     &
                    TimerStop

USE Timer_Variables_Module, &
            ONLY :  Timer_Core_Utilities


IMPLICIT NONE


CONTAINS



!+101+##################################################################!
!                                                                       !
!          Poseidon_Calc_ADM_Mass                     				    !
!                                                                       !
!#######################################################################!
SUBROUTINE Poseidon_Calc_ADM_Mass( ADM_Mass )

REAL(idp), INTENT(INOUT)                           ::  ADM_Mass

CALL TimerStart( Timer_Core_Utilities )
CALL Calc_ADM_Mass( ADM_Mass )
CALL TimerStop( Timer_Core_Utilities )

END SUBROUTINE Poseidon_Calc_ADM_Mass


!+102+##################################################################!
!                                                                       !
!          Poseidon_Calc_ADM_Mass_Parts                                 !
!                                                                       !
!#######################################################################!
SUBROUTINE Poseidon_Calc_ADM_Mass_Parts( ADM_Mass, ADM_Phys, ADM_Curve )

REAL(idp), INTENT(INOUT)                           ::  ADM_Mass
REAL(idp), INTENT(INOUT)                           ::  ADM_Phys
REAL(idp), INTENT(INOUT)                           ::  ADM_Curve

CALL TimerStart( Timer_Core_Utilities )
CALL Calc_ADM_Mass_In_Parts( ADM_Mass, ADM_Phys, ADM_Curve )
CALL TimerStop( Timer_Core_Utilities )

END SUBROUTINE Poseidon_Calc_ADM_Mass_Parts





!+201+##################################################################!
!                                                                       !
!          Poseidon_Calc_Komar_Mass                                     !
!                                                                       !
!#######################################################################!
SUBROUTINE Poseidon_Calc_Komar_Mass( Komar_Mass )

REAL(idp), INTENT(INOUT)                           ::  Komar_Mass

CALL TimerStart( Timer_Core_Utilities )
CALL Calc_Komar_Mass( Komar_Mass )
CALL TimerStop( Timer_Core_Utilities )

END SUBROUTINE Poseidon_Calc_Komar_Mass



END MODULE Poseidon_Utilities_Module
