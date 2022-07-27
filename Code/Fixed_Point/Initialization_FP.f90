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
            ONLY :  Domain_Dim,             &
                    DEGREE,                 &
                    L_LIMIT,                &
                    Num_CFA_Eqs,            &
                    CFA_Eq_Flags,           &
                    Verbose_Flag

USE Variables_Derived, &
            ONLY :  LM_Length,              &
                    Beta_Elem_Prob_Dim

USE Variables_Mesh, &
            ONLY :  Num_R_Elements

USE Variables_Functions, &
            ONLY :  Calc_3D_Values_At_Location,    &
                    Calc_1D_CFA_Values

USE Variables_Matrices, &
            ONLY :  Laplace_NNZ,                &
                    Beta_Diagonals,             &
                    Beta_Bandwidth

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
            ONLY :  Timer_FP_Initialization,      &
                    Timer_FP_Matrix_Init






IMPLICIT NONE




                    !*F&S*==========================================!
                    !                                               !
                    !           Functions & Subroutines             !
                    !                                               !
                    !===============================================!
CONTAINS










 !+101+####################################################################################!
!                                                                                           !
!       Initialize_FP                                                          !
!                                                                                           !
!===========================================================================================!
!                                                                                           !
 !#########################################################################################!
SUBROUTINE Initialize_CFA( CFA_EQ_Flags_Input )

INTEGER, DIMENSION(5), INTENT(IN), OPTIONAL             ::  CFA_EQ_Flags_Input



IF ( Verbose_Flag ) THEN
    PRINT*,"-Initializing Fixed Point Method variables. "
END IF
CALL TimerStart( Timer_FP_Initialization )


Laplace_NNZ = NUM_R_ELEMENTS*(DEGREE + 1)*(DEGREE + 1) - NUM_R_ELEMENTS + 1
Beta_Diagonals = Beta_Elem_Prob_Dim
Beta_Bandwidth = 2*Beta_Diagonals+1


CALL Allocate_CFA_Linear_Systems()

CALL TimerStart( Timer_FP_Matrix_Init )
CALL Initialize_XCFC_Matrices()
CALL TimerStop( Timer_FP_Matrix_Init )


Calc_3D_Values_At_Location  => Calc_FP_Values_At_Location
Calc_1D_CFA_Values          => Calc_1D_CFA_Values_FP




CALL TimerStop( Timer_FP_Initialization )


END SUBROUTINE Initialize_CFA











END MODULE Initialization_CFA

