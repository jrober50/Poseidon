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


IMPLICIT NONE


CONTAINS



!+101+##################################################################!
!                                                                       !
!          Poseidon_Calc_ADM_Mass                     				    !
!                                                                       !
!#######################################################################!
SUBROUTINE Poseidon_Calc_ADM_Mass( ADM_Mass )

REAL(idp), INTENT(INOUT)                           ::  ADM_Mass

CALL Calc_ADM_Mass( ADM_Mass )

END SUBROUTINE Poseidon_Calc_ADM_Mass


!+101+##################################################################!
!                                                                       !
!          Poseidon_Calc_ADM_Mass_Parts                                 !
!                                                                       !
!#######################################################################!
SUBROUTINE Poseidon_Calc_ADM_Mass_Parts( ADM_Mass, ADM_Phys, ADM_Curve )

REAL(idp), INTENT(INOUT)                           ::  ADM_Mass
REAL(idp), INTENT(INOUT)                           ::  ADM_Phys
REAL(idp), INTENT(INOUT)                           ::  ADM_Curve

CALL Calc_ADM_Mass_In_Parts( ADM_Mass, ADM_Phys, ADM_Curve )

END SUBROUTINE Poseidon_Calc_ADM_Mass_Parts



END MODULE Poseidon_Utilities_Module
