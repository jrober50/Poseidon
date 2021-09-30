   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Poseidon_MPI_Utilities_Module                                         !##!
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
USE Variables_MPI, &
        ONLY :  myID_Poseidon,      &
                nPROCS_Poseidon,    &
                MasterID_Poseidon,  &
                Poseidon_Comm_World

USE MPI


IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!          STOP_MPI                                                				!
!                                                                               !
!###############################################################################!
SUBROUTINE STOP_MPI( ierr )

INTEGER, INTENT(INOUT)                 :: ierr

INTEGER                             :: errcode = 13

CALL MPI_Abort(Poseidon_Comm_World, errcode, ierr)
STOP

END SUBROUTINE STOP_MPI



!+101+##########################################################################!
!                                                                               !
!          STOP_MPI                                                                !
!                                                                               !
!###############################################################################!
SUBROUTINE MPI_Master_Print( str )

CHARACTER(Len=*), INTENT(IN)                 :: Str

IF ( myID_Poseidon == MasterID_Poseidon ) THEN
    WRITE(*,'(A,I2.2,A,A)')"myID : ",myID_Poseidon," - ",Str
END IF

END SUBROUTINE MPI_Master_Print




!                                                                               !
!          STOP_MPI                                                                !
!                                                                               !
!###############################################################################!
SUBROUTINE MPI_All_Print( str )

CHARACTER(Len=*), INTENT(IN)                 :: Str

WRITE(*,'(A,I2.2,A,A)')"myID : ",myID_Poseidon," - ",Str

END SUBROUTINE MPI_All_Print



END MODULE Poseidon_MPI_Utilities_Module
