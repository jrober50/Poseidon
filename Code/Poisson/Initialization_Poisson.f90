   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Initialization_Poisson                                                       !##!
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
!##!    +101+   Poseidon_Initialize                                                 !##!
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
            ONLY :  Init_Message

USE Poseidon_Parameters, &
            ONLY :  Degree,                 &
                    Verbose_Flag

USE Variables_Mesh, &
            ONLY :  Num_R_Elements

USE Variables_Matrices, &
            ONLY :  Laplace_NNZ

USE Allocation_Poisson_Linear_System, &
            ONLY :  Allocate_Poisson_Linear_System

USE Matrix_Initialization_Module, &
            ONLY :  Initialize_Poisson_Matrix

USE Timer_Routines_Module, &
            ONLY :  TimerStart,                     &
                    TimerStop

USE Timer_Variables_Module, &
            ONLY :  Timer_Poisson_Initialization,      &
                    Timer_Poisson_Matrix_Init

USE Flags_Initialization_Module, &
            ONLY :  lPF_Init_Flags,     &
                    iPF_Init_Method_Vars
IMPLICIT NONE




                    !*F&S*==========================================!
                    !                                               !
                    !           Functions & Subroutines             !
                    !                                               !
                    !===============================================!
CONTAINS










 !+101+####################################################################################!
!                                                                                           !
!       Initialize_XCFC                                                          !
!                                                                                           !
!===========================================================================================!
!                                                                                           !
 !#########################################################################################!
SUBROUTINE Initialize_Poisson( )


IF ( .NOT. lPF_Init_Flags(iPF_Init_Method_Vars) ) THEN


    IF ( Verbose_Flag ) CALL Init_Message('Initializing XCFC system variables.')

    CALL TimerStart( Timer_Poisson_Initialization )


    Laplace_NNZ = NUM_R_ELEMENTS*(DEGREE + 1)*(DEGREE + 1) - NUM_R_ELEMENTS + 1

    CALL Allocate_Poisson_Linear_System()

    CALL TimerStart( Timer_Poisson_Matrix_Init )
    CALL Initialize_Poisson_Matrix()
    CALL TimerStop( Timer_Poisson_Matrix_Init )

    CALL TimerStop( Timer_Poisson_Initialization )

    lPF_Init_Flags(iPF_Init_Method_Vars) = .TRUE.

END IF

END SUBROUTINE Initialize_Poisson








END MODULE Initialization_Poisson


