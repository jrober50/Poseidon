   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE XCFC_Method_Module                                                           !##!
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

USE Parameters_Variable_Indices, &
            ONLY :  iU_CF,                      &
                    iU_LF,                      &
                    iU_S1,                      &
                    iU_S2,                      &
                    iU_S3,                      &
                    iU_X1,                      &
                    iU_X2,                      &
                    iU_X3,                      &
                    iVB_S,                      &
                    iVB_X

USE Variables_MPI, &
            ONLY :  myID_Poseidon,              &
                    MasterID_Poseidon,          &
                    Poseidon_Comm_World,        &
                    nPROCS_Poseidon

USE Variables_Derived, &
            ONLY :  Num_R_Nodes,                &
                    LM_Length

USE Variables_FP,  &
            ONLY :  FP_Coeff_Vector_A,          &
                    FP_Coeff_Vector_B,          &
                    CFA_EQ_Flags

USE Variables_IO, &
            ONLY :  Write_Flags,                &
                    Print_Flags

USE IO_Print_Results, &
            ONLY :  Print_Results

USE Poseidon_IO_Module, &
            ONLY :  Output_Final_Results

USE XCFC_Solvers_Main_Module ,  &
            ONLY :  XCFC_X_Solve,               &
                    XCFC_ConFactor_Solve,       &
                    XCFC_Lapse_Solve,           &
                    XCFC_Shift_Solve

USE XCFC_Source_Variables_Module, &
            ONLY :  Allocate_XCFC_Source_Variables, &
                    Deallocate_XCFC_Source_Variables

USE Poseidon_MPI_Utilities_Module, &
            ONLY :  STOP_MPI,               &
                    MPI_Master_Print,       &
                    MPI_All_Print

USE MPI

IMPLICIT NONE




CONTAINS


!+101+##########################################################################!
!                                                                               !
!                       XCFC_Method                                             !
!                                                                               !
!###############################################################################!
SUBROUTINE XCFC_Method()

LOGICAL                     :: PR = .FALSE.


CALL Allocate_XCFC_Source_Variables()


IF ( Verbose_Flag ) THEN
    PRINT*,"Begining XCFC Fixed Point Iterative Solve."
END IF



CALL Output_Initial_Guess(PR)



! Solve for X
CALL XCFC_X_Solve()



! Solve for Conformal Factor
IF ( CFA_Eq_Flags(iU_CF) == 1 ) THEN
    CALL XCFC_ConFactor_Solve()
END IF




! Solve for Lapse Function
IF ( CFA_Eq_Flags(iU_LF) == 1 ) THEN
    CALL XCFC_Lapse_Solve()
END IF




! Solve for Shift Vector
IF ( ANY(CFA_Eq_Flags(iU_S1:iU_S3) == 1) ) THEN
    CALL XCFC_Shift_Solve()
END IF




CALL Deallocate_XCFC_Source_Variables()


IF (myID_Poseidon == MasterID_Poseidon ) THEN
    CALL Output_Final_Results()
END IF


END SUBROUTINE XCFC_Method












!+101+##########################################################################!
!                                                                               !
!                       XCFC_Method                                             !
!                                                                               !
!###############################################################################!
SUBROUTINE XCFC_Method_New()

!INTEGER                         :: ierr
LOGICAL                         :: PR = .FALSE.




CALL Allocate_XCFC_Source_Variables()


IF ( Verbose_Flag ) THEN
    PRINT*,"Begining XCFC Fixed Point Iterative Solve."
END IF


CALL Output_Initial_Guess(PR)


CALL XCFC_X_Solve()

!PRINT*,"STOPing in XCFC_Method"
!CALL MPI_Finalize(ierr)
!STOP


CALL XCFC_ConFactor_Solve()

CALL XCFC_Lapse_Solve()

CALL XCFC_Shift_Solve()


CALL Deallocate_XCFC_Source_Variables()


CALL Output_Final_Results()

END SUBROUTINE XCFC_Method_New









!+101+##########################################################################!
!                                                                               !
!                       XCFC_Method                                             !
!                                                                               !
!###############################################################################!
SUBROUTINE XCFC_Method_Part1()

LOGICAL                             :: PR = .FALSE.




CALL Allocate_XCFC_Source_Variables()


IF ( Verbose_Flag ) THEN
    PRINT*,"Begining XCFC Metric Solve (Part 1 of 2)."
END IF


CALL Output_Initial_Guess(PR)


! Solve for X
CALL XCFC_X_Solve()



! Solve for Conformal Factor
IF ( CFA_Eq_Flags(1) == 1 ) THEN
    CALL XCFC_ConFactor_Solve()
END IF

CALL Deallocate_XCFC_Source_Variables()

END SUBROUTINE XCFC_Method_Part1


















!+101+##########################################################################!
!                                                                               !
!                       XCFC_Method                                             !
!                                                                               !
!###############################################################################!
SUBROUTINE XCFC_Method_Part2()

CALL Allocate_XCFC_Source_Variables()

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








!+301+##########################################################################!
!                                                                               !
!         Output_Initial_Guess                                                  !
!                                                                               !
!###############################################################################!
SUBROUTINE Output_Initial_Guess( PR )

LOGICAL, INTENT(INOUT)             :: PR

IF ( myID_Poseidon == 1 ) THEN
IF ( (Write_Flags(5) == 1) .OR. (Write_Flags(5) == 3) ) THEN
    PR = .TRUE.
    WRITE(*,'(A)')"Initial Guess Values"
    CALL Print_Results()
    PRINT*," "
END IF
END IF

END SUBROUTINE OUTPUT_Initial_Guess





END MODULE XCFC_Method_Module
