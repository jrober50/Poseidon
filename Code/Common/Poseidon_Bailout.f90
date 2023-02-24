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

USE MPI

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


CALL MPI_Initialized(MPI_Flag,iErr)
IF (MPI_Flag) CALL MPI_Finalize(iErr)



STOP Message


END SUBROUTINE Poseidon_Bailout











END MODULE Poseidon_Bailout_Module
