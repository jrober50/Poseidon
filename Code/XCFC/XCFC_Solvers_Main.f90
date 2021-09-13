   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE XCFC_Solvers_Main_Module                                                     !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!


USE Poseidon_Kinds_Module, &
        ONLY :  idp

USE Poseidon_Parameters, &
        ONLY :  Verbose_Flag

USE Parameters_Variable_Indices, &
        ONLY :  iVB_X,                      &
                iVB_S,                      &
                iU_CF,                      &
                iU_LF,                      &
                iU_S1,                      &
                iU_S2,                      &
                iU_S3,                      &
                iU_X1,                      &
                iU_X2,                      &
                iU_X3

USE XCFC_Source_Vector_TypeB_Module, &
        ONLY :  XCFC_Calc_Source_Vector_TypeB

USE XCFC_System_Solvers_Module, &
        ONLY :  XCFC_Solve_System_TypeB

USE XCFC_Fixed_Point_Module, &
        ONLY :  XCFC_Fixed_Point


IMPLICIT NONE

CONTAINS

!+201+##########################################################################!
!                                                                               !
!                       Solve_X_System                                             !
!                                                                               !
!###############################################################################!
SUBROUTINE XCFC_X_Solve()

INTEGER, DIMENSION(3)                                   ::  iU
INTEGER                                                 ::  iVB

IF ( Verbose_Flag ) THEN
    PRINT*,"Begining X system. "
END IF

iU = [iU_X1, iU_X2, iU_X3]
iVB = iVB_X

CALL XCFC_Calc_Source_Vector_TypeB( iU, iVB )
CALL XCFC_Solve_System_TypeB( iU )

END SUBROUTINE XCFC_X_Solve




!+201+##########################################################################!
!                                                                               !
!                       Solve_Shift_System                                      !
!                                                                               !
!###############################################################################!
SUBROUTINE XCFC_Shift_Solve()

INTEGER, DIMENSION(3)                                   ::  iU
INTEGER                                                 ::  iVB

IF ( Verbose_Flag ) THEN
    PRINT*,"Begining Shift system. "
END IF

iU = [iU_S1, iU_S2, iU_S3]
iVB = iVB_S

CALL XCFC_Calc_Source_Vector_TypeB( iU, iVB )
CALL XCFC_Solve_System_TypeB( iU )

END SUBROUTINE XCFC_Shift_Solve






!+101+##########################################################################!
!                                                                               !
!                       XCFC_ConFactor_Solve                                  !
!                                                                               !
!###############################################################################!
SUBROUTINE XCFC_ConFactor_Solve()


IF ( Verbose_Flag ) THEN
    PRINT*,"Begining Conformal Factor system. "
END IF

CALL XCFC_Fixed_Point(iU_CF)


END SUBROUTINE XCFC_ConFactor_Solve




!+101+##########################################################################!
!                                                                               !
!                       XCFC_ConFactor_Solve                                  !
!                                                                               !
!###############################################################################!
SUBROUTINE XCFC_Lapse_Solve()


IF ( Verbose_Flag ) THEN
    PRINT*,"Begining Lapse Function system. "
END IF

CALL XCFC_Fixed_Point(iU_LF)


END SUBROUTINE XCFC_Lapse_Solve




END MODULE XCFC_Solvers_Main_Module

