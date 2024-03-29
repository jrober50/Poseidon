   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Flags_Check_Routines                                                  !##!
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


USE Flags_Initialization_Module, &
            ONLY :  Poseidon_Initialization_Check

IMPLICIT NONE


CONTAINS



!+101+####################################################!
!                                                           !
!          Poseidon_Run_Check                               !
!                                                           !
!#########################################################!
LOGICAL FUNCTION Poseidon_Run_Check()

LOGICAL, DIMENSION(1:4)             ::  Flags

Flags = [   Poseidon_Initialization_Check(),        &
            Poseidon_Source_Input_Check(),          &
            Poseidon_Boundary_Condition_Check(),    &
            Poseidon_Guess_Check()                  ]

IF ( ALL(Flags) ) THEN
    Poseidon_Run_Check = .TRUE.
    IF ( Verbose_Flag ) THEN
        WRITE(*,'(A)') "Poseidon Run Check :  Passed."
    END IF
ELSE
    Poseidon_Run_Check = .FALSE.
    STOP "Poseidon Run Check : Failed."
END IF


END FUNCTION Poseidon_Run_Check







 !+103+####################################################!
!                                                           !
!          Poseidon_Source_Input_Check                      !
!                                                           !
 !#########################################################!
LOGICAL FUNCTION Poseidon_Source_Input_Check()

IF ( .TRUE. ) THEN
    Poseidon_Source_Input_Check = .TRUE.
ELSE
    Poseidon_Source_Input_Check = .FALSE.
END IF

END FUNCTION Poseidon_Source_Input_Check



 !+103+####################################################!
!                                                           !
!          Poseidon_Boundary_Condition_Check                !
!                                                           !
 !#########################################################!
LOGICAL FUNCTION Poseidon_Boundary_Condition_Check()

IF ( .TRUE. ) THEN
    Poseidon_Boundary_Condition_Check = .TRUE.
ELSE
    Poseidon_Boundary_Condition_Check = .FALSE.
END IF

END FUNCTION Poseidon_Boundary_Condition_Check




 !+103+####################################################!
!                                                           !
!          Poseidon_Boundary_Condition_Check                !
!                                                           !
 !#########################################################!
LOGICAL FUNCTION Poseidon_Guess_Check()

IF ( .TRUE. ) THEN
    Poseidon_Guess_Check = .TRUE.
ELSE
    Poseidon_Guess_Check = .FALSE.
END IF

END FUNCTION Poseidon_Guess_Check




END MODULE Flags_Check_Routines
