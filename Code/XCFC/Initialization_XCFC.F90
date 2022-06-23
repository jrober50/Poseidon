   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Initialization_XCFC                                                            !##!
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

USE Poseidon_Message_Routines_Module, &
            ONLY :  Init_Message


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

USE Variables_FP, &
            ONLY :  CFA_EQ_Map,                 &
                    CFA_Var_Map,                &
                    CFA_Mat_Map

USE Variables_Matrices, &
            ONLY :  Laplace_NNZ,                &
                    Beta_Diagonals,             &
                    Beta_Bandwidth,             &
                    Num_Matrices

USE Allocation_XCFC_Linear_Systems, &
            ONLY :  Allocate_XCFC_Linear_Systems

USE Return_Functions_FP,   &
            ONLY :  Calc_FP_Values_At_Location,  &
                    Calc_1D_CFA_Values_FP

USE Matrix_Initialization_Module, &
            ONLY :  Initialize_XCFC_Matrices

USE Timer_Routines_Module, &
            ONLY :  TimerStart,                     &
                    TimerStop

USE Timer_Variables_Module, &
            ONLY :  Timer_XCFC_Initialization,      &
                    Timer_XCFC_Matrix_Init

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
SUBROUTINE Initialize_XCFC( )


IF ( .NOT. lPF_Init_Flags(iPF_Init_Method_Vars) ) THEN

    IF ( Verbose_Flag ) CALL Init_Message('Initializing XCFC system variables.')

    CALL TimerStart( Timer_XCFC_Initialization )



    Beta_Diagonals = Beta_Elem_Prob_Dim-1
    Beta_Bandwidth = 2*Beta_Diagonals+1


    Laplace_NNZ = NUM_R_ELEMENTS*(DEGREE + 1)*(DEGREE + 1) - NUM_R_ELEMENTS + 1

    CALL Allocate_XCFC_Linear_Systems()

    CALL TimerStart( Timer_XCFC_Matrix_Init )
    CALL Initialize_XCFC_Matrices()
    CALL TimerStop( Timer_XCFC_Matrix_Init )



    Calc_3D_Values_At_Location  => Calc_FP_Values_At_Location
    Calc_1D_CFA_Values          => Calc_1D_CFA_Values_FP



    !lPF_Init_Flag(iPF_Init_XCFC) = .TRUE.
    CALL TimerStop( Timer_XCFC_Initialization )

    lPF_Init_Flags(iPF_Init_Method_Vars) = .TRUE.
END IF

END SUBROUTINE Initialize_XCFC






 !+102+####################################################################################!
!                                                                                           !
!       Create_Eq_Maps                                                          !
!                                                                                           !
 !#########################################################################################!
SUBROUTINE Create_Eq_Maps()

INTEGER                                                 ::  i, j


CFA_EQ_Map = -1
j = 1
DO i = 1,5
    IF ( CFA_EQ_Flags(i) == 1 ) THEN
        CFA_EQ_Map(j) = i
        j = j+1
    END IF
END DO


CFA_Var_Map = -1
j = 1
DO i = 1,5
    IF ( CFA_EQ_Flags(i) == 1 ) THEN
        CFA_Var_Map(i) = j
        j = j+1
    END IF
END DO




!
!   Calculate the number of matrices to be created and stored.
!
! The Psi and AlphaPsi equations can share the same Lapace Matrix.
! Each of the Shift Vector componenets will need its own matrix.
Num_Matrices = 0
IF (CFA_EQ_Flags(1) == 1)  THEN
    Num_Matrices = Num_Matrices + 1
END IF
IF (CFA_EQ_Flags(2) == 1)  THEN
    Num_Matrices = Num_Matrices + 1
END IF



IF ( ( CFA_EQ_Flags(1) == 1 ) .OR. ( CFA_EQ_Flags(2) == 1 ) ) THEN
    IF ( CFA_EQ_Flags(1) == 1 ) THEN
        CFA_Mat_Map(1) = 1
    END IF
    IF ( CFA_EQ_Flags(2) == 1 ) THEN
        CFA_Mat_Map(2) = 1
    END IF
    j = 2
    DO i = 3,5
        IF ( CFA_EQ_Flags(i) == 1 ) THEN
            CFA_Mat_Map(i) = j
            j = j + 1
        END IF
    END DO
ELSE

    j = 1
    DO i = 3,5
        IF ( CFA_EQ_Flags(i) == 1 ) THEN
            CFA_Mat_Map(i) = j
            j = j + 1
        END IF
    END DO
END IF


!lPF_Init_Flag(iPF_Init_Eq_Maps) = .TRUE.


END SUBROUTINE Create_Eq_Maps








END MODULE Initialization_XCFC

