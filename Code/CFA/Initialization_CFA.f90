   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Initialization_CFA                                                            !##!
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
            ONLY :  idp, fdp

USE Poseidon_Numbers_Module, &
            ONLY :  pi


USE Poseidon_Parameters, &
            ONLY :  Degree,                 &
                    L_Limit,                &
                    Verbose_Flag

USE Variables_Derived, &
            ONLY :  LM_Length,              &
                    iVB_Elem_Prob_Dim

USE Variables_Mesh, &
            ONLY :  Num_R_Elements

USE Variables_Functions, &
            ONLY :  Calc_3D_Values_At_Location,    &
                    Calc_1D_CFA_Values

USE Variables_Matrices, &
            ONLY :  Laplace_NNZ,                &
                    iMB_Diagonals,             &
                    iMB_Bandwidth

USE Allocation_CFA_Linear_Systems, &
            ONLY :  Allocate_CFA_Linear_Systems


USE Return_Functions_FP,   &
            ONLY :  Calc_FP_Values_At_Location,  &
                    Calc_1D_CFA_Values_FP

USE Matrix_Initialization_Module, &
            ONLY :  Initialize_XCFC_Matrices

USE Timer_Routines_Module, &
            ONLY :  TimerStart,                     &
                    TimerStop

USE Timer_Variables_Module, &
            ONLY :  Timer_Initialization,      &
                    Timer_Matrix_Init






IMPLICIT NONE




                    !*F&S*==========================================!
                    !                                               !
                    !           Functions & Subroutines             !
                    !                                               !
                    !===============================================!
CONTAINS










 !+101+####################################################################################!
!                                                                                           !
!       Initialize_CFA                                                          !
!                                                                                           !
!===========================================================================================!
!                                                                                           !
 !#########################################################################################!
SUBROUTINE Initialize_CFA( CFA_EQ_Flags_Input )

INTEGER, DIMENSION(5), INTENT(IN), OPTIONAL             ::  CFA_EQ_Flags_Input



IF ( Verbose_Flag ) THEN
    PRINT*,"-Initializing Fixed Point Method variables. "
END IF
CALL TimerStart( Timer_Initialization )


Laplace_NNZ = NUM_R_ELEMENTS*(DEGREE + 1)*(DEGREE + 1) - NUM_R_ELEMENTS + 1
iMB_Diagonals = iVB_Elem_Prob_Dim
iMB_Bandwidth = 2*iMB_Diagonals+1


CALL Allocate_CFA_Linear_Systems()

CALL TimerStart( Timer_Matrix_Init )
CALL Initialize_XCFC_Matrices()
CALL TimerStop( Timer_Matrix_Init )


Calc_3D_Values_At_Location  => Calc_FP_Values_At_Location
Calc_1D_CFA_Values          => Calc_1D_CFA_Values_FP




CALL TimerStop( Timer_Initialization )


END SUBROUTINE Initialize_CFA











END MODULE Initialization_CFA

