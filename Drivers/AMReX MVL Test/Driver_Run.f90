   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Driver_Run_Module                                                     !##!
!##!                                                                         !##!
!##!_________________________________________________________________________!##!
!##!                                                                         !##!
!##!                                                                         !##!
!##!=========================================================================!##!
!##!                                                                         !##!
!##!    Contains:                                                            !##!
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
USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag

USE Variables_MPI, &
            ONLY :  myID_Poseidon,              &
                    MasterID_Poseidon

USE Poseidon_IO_Module, &
            ONLY :  Output_Final_Results

USE XCFC_Solvers_Main_Module ,  &
            ONLY :  XCFC_X_Solve

USE XCFC_Source_Variables_Module, &
            ONLY :  Allocate_XCFC_Source_Variables, &
                    Deallocate_XCFC_Source_Variables

IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!                                                     				!
!                                                                               !
!###############################################################################!
SUBROUTINE Driver_Run()

CALL Allocate_XCFC_Source_Variables()


IF ( Verbose_Flag ) THEN
    PRINT*,"Begining Single Modified Vector Laplacian Solve."
END IF

! Solve for X
CALL XCFC_X_Solve()



CALL Deallocate_XCFC_Source_Variables()


IF (myID_Poseidon == MasterID_Poseidon ) THEN
    CALL Output_Final_Results()
END IF



END SUBROUTINE Driver_Run


END MODULE Driver_Run_Module
