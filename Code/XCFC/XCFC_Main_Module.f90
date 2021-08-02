   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE XCFC_Main_Module                                                             !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains the top level subroutines needed to inialize, run, and close       !##!
!##!        Poseidon.                                                               !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!

!*D*================================!
!                                   !
!           Dependencies            !
!                                   !
!===================================!
USE Poseidon_Kinds_Module, &
            ONLY : idp

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag


USE Variables_FP,  &
            ONLY :  FP_Coeff_Vector,            &
                    FP_Coeff_Vector_Orig,       &
                    CFA_EQ_Flags

USE Variables_IO, &
            ONLY :  Write_Flags,                &
                    Print_Flags

USE IO_Print_Results, &
            ONLY :  Print_Results

USE Poseidon_IO_Module, &
            ONLY :  Output_Final_Results

USE XCFC_X_Module ,  &
            ONLY :  XCFC_X_Solve

USE XCFC_ConFactor_Module, &
            ONLY :  XCFC_ConFactor_Solve

USE XCFC_Lapse_Module, &
            ONLY :  XCFC_Lapse_Solve

USE XCFC_Shift_Module, &
            ONLY :  XCFC_Shift_Solve

USE XCFC_Source_Vector_Module, &
            ONLY :  Allocate_XCFC_Source_Variables, &
                    Deallocate_XCFC_Source_Variables


IMPLICIT NONE




CONTAINS


!+101+##########################################################################!
!                                                                               !
!                       XCFC_Method                                             !
!                                                                               !
!###############################################################################!
SUBROUTINE XCFC_Method()

INTEGER :: i, ConFact_Loop_Number
LOGICAL                                                 :: PR = .FALSE.



ConFact_Loop_Number = 1

CALL Allocate_XCFC_Source_Variables()


IF ( Verbose_Flag ) THEN
    PRINT*,"Begining XCFC Fixed Point Iterative Solve."
END IF

IF ( (Write_Flags(5) == 1) .OR. (Write_Flags(5) == 3) ) THEN
    PR = .TRUE.
    WRITE(*,'(A)')"Initial Guess Values"
    CALL Print_Results()
    PRINT*," "
END IF

! Save Origional FP_Coeff_Vector
FP_Coeff_Vector_Orig = FP_Coeff_Vector


! Solve for X
CALL XCFC_X_Solve()



! Solve for Conformal Factor
DO i = 1,ConFact_Loop_Number
    FP_Coeff_Vector_Orig( :, :, 1 ) = FP_Coeff_Vector( :, :, 1 )
    IF ( CFA_Eq_Flags(1) == 1 ) THEN
        CALL XCFC_ConFactor_Solve()
    END IF
END DO

! Define S*Block_Source_Si(rd, td, pd, re, te, pe, ui)

! Solve for Lapse Function
IF ( CFA_Eq_Flags(2) == 1 ) THEN
    CALL XCFC_Lapse_Solve()
END IF

! Solve for Shift Vector

IF ( ANY(CFA_Eq_Flags(3:5) == 1) ) THEN
    CALL XCFC_Shift_Solve()
END IF

CALL Deallocate_XCFC_Source_Variables()



CALL Output_Final_Results()

END SUBROUTINE XCFC_Method






















!+101+##########################################################################!
!                                                                               !
!                       XCFC_Method                                             !
!                                                                               !
!###############################################################################!
SUBROUTINE XCFC_Method_Part1()

INTEGER :: i, ConFact_Loop_Number
LOGICAL                                                 :: PR = .FALSE.



ConFact_Loop_Number = 1

CALL Allocate_XCFC_Source_Variables()


IF ( Verbose_Flag ) THEN
    PRINT*,"Begining XCFC Metric Solve (Part 1 of 2)."
END IF

IF ( (Write_Flags(5) == 1) .OR. (Write_Flags(5) == 3) ) THEN
    PR = .TRUE.
    WRITE(*,'(A)')"Initial Guess Values"
    CALL Print_Results()
    PRINT*," "
END IF


! Save Origional FP_Coeff_Vector
FP_Coeff_Vector_Orig = FP_Coeff_Vector



! Solve for X
CALL XCFC_X_Solve()



! Solve for Conformal Factor
DO i = 1,ConFact_Loop_Number
    FP_Coeff_Vector_Orig( :, :, 1 ) = FP_Coeff_Vector( :, :, 1 )
    IF ( CFA_Eq_Flags(1) == 1 ) THEN
        CALL XCFC_ConFactor_Solve()
    END IF
END DO



END SUBROUTINE XCFC_Method_Part1


















!+101+##########################################################################!
!                                                                               !
!                       XCFC_Method                                             !
!                                                                               !
!###############################################################################!
SUBROUTINE XCFC_Method_Part2()


IF ( Verbose_Flag ) THEN
    PRINT*,"Begining XCFC Metric Solve (Part 2 of 2)."
END IF


! Solve for Lapse Function
IF ( CFA_Eq_Flags(2) == 1 ) THEN
    CALL XCFC_Lapse_Solve()
END IF


! Solve for Shift Vector
IF ( ANY(CFA_Eq_Flags(3:5) == 1) ) THEN
    CALL XCFC_Shift_Solve()
END IF




CALL Deallocate_XCFC_Source_Variables()



CALL Output_Final_Results()


END SUBROUTINE XCFC_Method_Part2













END MODULE XCFC_Main_Module
