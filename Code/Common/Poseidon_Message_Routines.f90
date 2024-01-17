   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Poseidon_Message_Routines_Module                                      !##!
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
PRIVATE

PUBLIC ::   Init_Message
PUBLIC ::   Run_Message
PUBLIC ::   Warning_Message
PUBLIC ::   Driver_Init_Message


CONTAINS



 !+101+####################################################!
!                                                           !
!       Initialization_Message   	                        !
!                                                           !
 !#########################################################!
SUBROUTINE Init_Message( Message )

CHARACTER(LEN=*), INTENT(IN)    :: Message

WRITE(*,'(A,A)')'Poseidon Initialization : ',Message

END SUBROUTINE Init_Message

 !+201+####################################################!
!                                                           !
!       Run_Message                                         !
!                                                           !
 !#########################################################!
SUBROUTINE Run_Message( Message )

CHARACTER(LEN=*), INTENT(IN)    :: Message

WRITE(*,'(A,A)')'Poseidon Running : ',Message

END SUBROUTINE Run_Message


 !+301+####################################################!
!                                                           !
!       Warning_Message                                     !
!                                                           !
 !#########################################################!
SUBROUTINE Warning_Message( Message )

CHARACTER(LEN=*), INTENT(IN)    :: Message

WRITE(*,'(A,A)')'POSEIDON WARNING : ',Message

END SUBROUTINE Warning_Message



 !+301+####################################################!
!                                                           !
!       Driver_Init_Message                                 !
!                                                           !
 !#########################################################!
SUBROUTINE Driver_Init_Message( Message )

CHARACTER(LEN=*), INTENT(IN)    :: Message

WRITE(*,'(A,A)')'Poseidon Driver Initialization : ',Message

END SUBROUTINE Driver_Init_Message



END MODULE Poseidon_Message_Routines_Module
