   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_XCFC_Interface_Module                                               !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains the top level subroutines needed to run the Poseidon in XCFC mode. !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   Poseidon_XCFC_Run_Part1                                             !##!
!##!    +102+   Poseidon_XCFC_Run_Part2                                             !##!
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
USE Variables_MPI, &
            ONLY :  myID_Poseidon,              &
                    MasterID_Poseidon,          &
                    nPROCS_Poseidon,            &
                    Poseidon_Comm_World

USE XCFC_Method_Module, &
            ONLY :  XCFC_Method_Part1,      &
                    XCFC_Method_Part2

USE IO_Print_Results, &
            ONLY :  Print_Results

USE IO_Write_Final_Results, &
            ONLY :  Write_Final_Results

USE Flags_IO_Module, &
            ONLY :  lPF_IO_Flags,           &
                    iPF_IO_Print_Results,   &
                    iPF_IO_Write_Results

USE MPI


CONTAINS


!+101+######################################################################################!
!                                                                                           !
!       Poseidon_XCFC_Run_Part1 - Performs the steps to calculate the XCFC Conformal        !
!                                 Factor.  Upon completion of this call, the routine        !
!                                 Poseidon_Return_ConFactor() can be used to acquire        !
!                                 values of the conformal factor at desired locations.      !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_XCFC_Run_Part1()

INTEGER         :: i, ierr

CALL XCFC_Method_Part1()

IF ( lPF_IO_Flags(iPF_IO_Print_Results) ) THEN
    DO i = 0,nPROCs_Poseidon
        IF (myID_Poseidon == MasterID_Poseidon) THEN
            PRINT*,"myID_Poseidon :",myID_Poseidon
            Call Print_Results()
        END IF
        CALL MPI_Barrier(Poseidon_Comm_World,ierr)
    END DO

END IF

IF ( lPF_IO_Flags(iPF_IO_Write_Results) ) THEN
IF ( myID_Poseidon == MasterID_Poseidon ) THEN
    Call Write_Final_Results()
END IF
END IF


END SUBROUTINE Poseidon_XCFC_Run_Part1





!+102+######################################################################################!
!                                                                                           !
!       Poseidon_XCFC_Run_Part2 - Performs the steps to calculate the XCFC Lapse Function,  !
!                                 and Shift Vector.  Upon completion of this call, the      !
!                                 routine, Poseidon_Return_Lapse() and                      !
!                                 Poseidon_Return_Shift() can be used to acquire            !
!                                 values of the those variables at desired locations.       !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_XCFC_Run_Part2()

INTEGER         :: i, ierr

CALL XCFC_Method_Part2()

IF ( lPF_IO_Flags(iPF_IO_Print_Results) ) THEN
    DO i = 0,nPROCs_Poseidon
        IF (myID_Poseidon == MasterID_Poseidon) THEN
            PRINT*,"myID_Poseidon :",myID_Poseidon
            Call Print_Results()
        END IF
        CALL MPI_Barrier(Poseidon_Comm_World,ierr)
    END DO
END IF

IF ( lPF_IO_Flags(iPF_IO_Write_Results) ) THEN
IF ( myID_Poseidon == MasterID_Poseidon ) THEN
    Call Write_Final_Results()
END IF
END IF


END SUBROUTINE Poseidon_XCFC_Run_Part2







END MODULE Poseidon_XCFC_Interface_Module
