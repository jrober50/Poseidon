   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Poseidon_Bailout_Module                                               !##!
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


IMPLICIT NONE


CONTAINS

 !+101+####################################################!
!                                                           !
!       Poseidon_Bailout   	                                !
!                                                           !
 !#########################################################!
SUBROUTINE Poseidon_Bailout( Message )

CHARACTER(LEN = *),     INTENT(IN),     OPTIONAL    ::  Message

LOGICAL                                             ::  MPI_Flag
INTEGER                                             ::  iErr

IF ( MPI_Initialize(MPI_Flag,iErr)) THEN
    CALL MPI_Finalize(iErr)
END IF


STOP Message


END SUBROUTINE Poseidon_Bailout











END MODULE Poseidon_Bailout_Module