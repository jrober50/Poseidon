   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Initialization_XCFC                                                          !##!
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

USE Variables_Derived, &
            ONLY :  iVB_Elem_Prob_Dim

USE Variables_Mesh, &
            ONLY :  Num_R_Elements
USE Variables_Matrices, &
            ONLY :  Laplace_NNZ,                &
                    iMB_Diagonals,             &
                    iMB_Bandwidth

USE Allocation_XCFC_Linear_Systems, &
            ONLY :  Allocate_XCFC_Linear_Systems

USE Matrix_Initialization_Module, &
            ONLY :  Initialize_XCFC_Matrices

USE Timer_Routines_Module, &
            ONLY :  TimerStart,                     &
                    TimerStop

USE Timer_Variables_Module, &
            ONLY :  Timer_Initialization,      &
                    Timer_Matrix_Init

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










 !+101+################################################################!
!                                                                       !
!          Initialize_XCFC                                              !
!                                                                       !
!=======================================================================!
!                                                                       !
 !#####################################################################!
SUBROUTINE Initialize_XCFC( )


IF ( .NOT. lPF_Init_Flags(iPF_Init_Method_Vars) ) THEN

    IF ( Verbose_Flag ) CALL Init_Message('Initializing XCFC system variables.')

    CALL TimerStart( Timer_Initialization )



    iMB_Diagonals = iVB_Elem_Prob_Dim-1
    iMB_Bandwidth = 2*iMB_Diagonals+1

    Laplace_NNZ = NUM_R_ELEMENTS*(DEGREE + 1)*(DEGREE + 1) - NUM_R_ELEMENTS + 1


    ! Allocate Arrays
    CALL Allocate_XCFC_Linear_Systems()


    ! Construct Matrices
    CALL Initialize_XCFC_Matrices()

    !lPF_Init_Flag(iPF_Init_XCFC) = .TRUE.
    CALL TimerStop( Timer_Initialization )

    lPF_Init_Flags(iPF_Init_Method_Vars) = .TRUE.
END IF

END SUBROUTINE Initialize_XCFC








END MODULE Initialization_XCFC

