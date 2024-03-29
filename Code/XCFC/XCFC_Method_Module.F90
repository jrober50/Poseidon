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
            ONLY :  idp

USE Poseidon_Message_Routines_Module, &
            ONLY :  Run_Message

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag,               &
                    Eq_Flags

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

USE Variables_IO, &
            ONLY :  Write_Flags,                &
                    Print_Flags

USE IO_Print_Results, &
            ONLY :  Print_Results


USE XCFC_Solvers_Main_Module ,  &
            ONLY :  XCFC_X_Solve,               &
                    XCFC_ConFactor_Solve,       &
                    XCFC_Lapse_Solve,           &
                    XCFC_Shift_Solve


USE Poseidon_MPI_Utilities_Module, &
            ONLY :  STOP_MPI,               &
                    MPI_Master_Print,       &
                    MPI_All_Print

USE Flags_IO_Module, &
            ONLY :  lPF_IO_Flags,           &
                    iPF_IO_Print_Results



#ifdef POSEIDON_MEMORY_FLAG
USE Poseidon_Memory_Routines, &
            ONLY :  Poseidon_Mark_Memory

USE Memory_Variables_Module, &
            ONLY :  Memory_Method_Start,     &
                    Memory_Method_Before_CF, &
                    Memory_Method_Before_LF, &
                    Memory_Method_Before_SV, &
                    Memory_Method_End
#endif

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



IF ( Verbose_Flag ) CALL Run_Message('Beginning XCFC System Solve.')


CALL Output_Initial_Guess(PR)


#ifdef POSEIDON_MEMORY_FLAG
CALL Poseidon_Mark_Memory(Memory_Method_Start)
PRINT*,"Before XV Solve                     : ",Memory_Method_Start
#endif

! Solve for X
CALL XCFC_X_Solve()



#ifdef POSEIDON_MEMORY_FLAG
CALL Poseidon_Mark_Memory(Memory_Method_Before_CF)
PRINT*,"Before CF Solve                     : ",Memory_Method_Before_CF

#endif

! Solve for Conformal Factor
IF ( Eq_Flags(iU_CF) == 1 ) THEN
    CALL XCFC_ConFactor_Solve()
END IF


#ifdef POSEIDON_MEMORY_FLAG
CALL Poseidon_Mark_Memory(Memory_Method_Before_LF)
PRINT*,"Before LF Solve                     : ",Memory_Method_Before_LF
#endif

! Solve for Lapse Function
IF ( Eq_Flags(iU_LF) == 1 ) THEN
    CALL XCFC_Lapse_Solve()
END IF


#ifdef POSEIDON_MEMORY_FLAG
CALL Poseidon_Mark_Memory(Memory_Method_Before_SV)
PRINT*,"Before SV Solve                     : ",Memory_Method_Before_SV
#endif

! Solve for Shift Vector
IF ( ANY(Eq_Flags(iU_S1:iU_S3) == 1) ) THEN
    CALL XCFC_Shift_Solve()
END IF


#ifdef POSEIDON_MEMORY_FLAG
CALL Poseidon_Mark_Memory(Memory_Method_End)
PRINT*,"After SV Solve                      : ",Memory_Method_End
#endif



END SUBROUTINE XCFC_Method












!+101+##########################################################################!
!                                                                               !
!                       XCFC_Method                                             !
!                                                                               !
!###############################################################################!
SUBROUTINE XCFC_Method_New()

!INTEGER                         :: ierr
LOGICAL                         :: PR = .FALSE.





IF ( Verbose_Flag ) CALL Run_Message('Beginning XCFC System Solve.')



CALL Output_Initial_Guess(PR)


CALL XCFC_X_Solve()

!PRINT*,"STOPing in XCFC_Method"
!CALL MPI_Finalize(ierr)
!STOP


CALL XCFC_ConFactor_Solve()

CALL XCFC_Lapse_Solve()

CALL XCFC_Shift_Solve()



END SUBROUTINE XCFC_Method_New









!+101+##########################################################################!
!                                                                               !
!                       XCFC_Method                                             !
!                                                                               !
!###############################################################################!
SUBROUTINE XCFC_Method_Part1()

LOGICAL                             :: PR = .FALSE.


IF ( Verbose_Flag ) CALL Run_Message('Begining XCFC Metric Solve (Part 1 of 2).')


CALL Output_Initial_Guess(PR)


! Solve for X
CALL XCFC_X_Solve()



! Solve for Conformal Factor
IF ( Eq_Flags(1) == 1 ) THEN
    CALL XCFC_ConFactor_Solve()
END IF



!IF ( lPF_IO_Flags(IPF_IO_Print_Results) ) THEN
!    CALL Print_Results()
!END IF



END SUBROUTINE XCFC_Method_Part1


















!+101+##########################################################################!
!                                                                               !
!                       XCFC_Method                                             !
!                                                                               !
!###############################################################################!
SUBROUTINE XCFC_Method_Part2()


IF ( Verbose_Flag ) CALL Run_Message('Begining XCFC Metric Solve (Part 2 of 2).')


! Solve for Lapse Function
IF ( Eq_Flags(2) == 1 ) THEN
    CALL XCFC_Lapse_Solve()
END IF


! Solve for Shift Vector
IF ( ANY(Eq_Flags(3:5) == 1) ) THEN
    CALL XCFC_Shift_Solve()
END IF


!IF ( lPF_IO_Flags(IPF_IO_Print_Results) ) THEN
!    CALL Print_Results()
!END IF





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
