   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Source_Input_XCFC_Native_Module                                      	     !##!
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
USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag

USE Poseidon_Message_Routines_Module, &
            ONLY :  Run_Message

USE Parameters_Variable_Indices, &
            ONLY :  iS_E,               &
                    iS_S,               &
                    iS_S1,              &
                    iS_S2,              &
                    iS_S3

USE Variables_Interface, &
            ONLY :  Caller_Set,                     &
                    Caller_NQ,                      &
                    Caller_Quad_DOF,                &
                    Caller_xL,                      &
                    Caller_RQ_xlocs,                &
                    Caller_TQ_xlocs,                &
                    Caller_PQ_xlocs,                &
                    Translation_Matrix


USE Variables_Mesh, &
            ONLY :  NUM_R_ELEMENTS,     &
                    NUM_T_ELEMENTS,     &
                    NUM_P_ELEMENTS,     &
                    drlocs,             &
                    dtlocs,             &
                    dplocs

USE Variables_Source, &
            ONLY :  Block_Source_E,     &
                    Block_Source_S,     &
                    Block_Source_Si

USE Variables_Quadrature, &
            ONLY :  Num_R_Quad_Points,      &
                    Num_T_Quad_Points,      &
                    Num_P_Quad_Points,      &
                    Num_TP_Quad_Points,     &
                    Local_Quad_DOF,         &
                    xLeftLimit,             &
                    xRightLimit,            &
                    Int_R_Locations


USE Functions_Translation_Matrix_Module, &
            ONLY :  Create_Translation_Matrix

USE Maps_Quadrature, &
            ONLY :  Quad_Map

USE Timer_Routines_Module, &
            ONLY :  TimerStart,     &
                    TimerSTop

USE Timer_Variables_Module, &
            ONLY :  Timer_GR_SourceInput,           &
                    Timer_GR_SourceInput_PartA




IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!                           Poseidon_Input_Sources                              !
!                                                                               !
!###############################################################################!
SUBROUTINE Poseidon_Input_Sources_Part1_Native( Input_E,            &
                                                Input_Si,           &
                                                Input_NE,           &
                                                Input_NQ,           &
                                                Input_R_Quad,       &
                                                Input_T_Quad,       &
                                                Input_P_Quad,       &
                                                Input_xL            )


INTEGER,    INTENT(IN), DIMENSION(3)                                ::  Input_NE
INTEGER,    INTENT(IN), DIMENSION(3)                                ::  Input_NQ

REAL(idp),  INTENT(IN), DIMENSION(  1:Input_NQ(1)*Input_NQ(2)*Input_NQ(3),  &
                                    0:Input_NE(1)-1,                        &
                                    0:Input_NE(2)-1,                        &
                                    0:Input_NE(3)-1  )                      ::  Input_E

REAL(idp),  INTENT(IN), DIMENSION(  1:Input_NQ(1)*Input_NQ(2)*Input_NQ(3),  &
                                    0:Input_NE(1)-1,                        &
                                    0:Input_NE(2)-1,                        &
                                    0:Input_NE(3)-1,                        &
                                    1:3         )                           ::  Input_Si



REAL(idp),  INTENT(IN), DIMENSION( 1:Input_NQ(1) )                  ::  Input_R_Quad
REAL(idp),  INTENT(IN), DIMENSION( 1:Input_NQ(2) )                  ::  Input_T_Quad
REAL(idp),  INTENT(IN), DIMENSION( 1:Input_NQ(3) )                  ::  Input_P_Quad
REAL(idp),  INTENT(IN), DIMENSION(2)                                ::  Input_xL




INTEGER                                             ::  Local_Here


REAL(idp), DIMENSION(:,:), ALLOCATABLE              ::  TransMat
INTEGER                                             ::  Their_DOF

INTEGER                                             ::  re, te, pe


IF ( Verbose_Flag ) CALL Run_Message('Receiving Part 1 XCFC Sources. Container : Fortran Array.')

CALL TimerStart(Timer_GR_SourceInput)
CALL TimerStart(Timer_GR_SourceInput_PartA)



! Define Interpolation Matrix
Their_DOF = Input_NQ(1)*Input_NQ(2)*Input_NQ(3)


ALLOCATE( TransMat(1:Their_DOF, 1:Local_Quad_DOF) )

TransMat =  Create_Translation_Matrix(  Input_NQ,          &
                                        Input_xL,          &
                                        Input_R_Quad,    &
                                        Input_T_Quad,    &
                                        Input_P_Quad,    &
                                        Their_DOF,         &
                                        [Num_R_Quad_Points, Num_T_Quad_Points, Num_P_Quad_Points ],            &
                                        [xLeftLimit, xRightLimit ],            &
                                        Int_R_Locations,      &
                                        Int_R_Locations,      &
                                        Int_R_Locations,      &
                                        Local_Quad_DOF            )


DO pe = 0,Input_NE(3)-1
DO te = 0,Input_NE(2)-1
DO re = 0,Input_NE(1)-1

DO Local_Here = 1,Local_Quad_DOF

    Block_Source_E(Local_Here,re,te,pe) = DOT_PRODUCT( TransMat(:,Local_Here),      &
                                                       Input_E(:,re,te,pe)          )

    Block_Source_Si(Local_Here,re,te,pe,1) = DOT_PRODUCT( TransMat(:,Local_Here),   &
                                                          Input_Si(:,re,te,pe,1)    )


    Block_Source_Si(Local_Here,re,te,pe,2) = DOT_PRODUCT( TransMat(:,Local_Here),   &
                                                          Input_Si(:,re,te,pe,2)    )

    Block_Source_Si(Local_Here,re,te,pe,3) = DOT_PRODUCT( TransMat(:,Local_Here),   &
                                                          Input_Si(:,re,te,pe,3)    )


END DO  ! Local_Here



END DO ! re Loop
END DO ! te Loop
END DO ! pe Loop


CALL TimerStop(Timer_GR_SourceInput)
CALL TimerStop(Timer_GR_SourceInput_PartA)




END SUBROUTINE Poseidon_Input_Sources_Part1_Native





!+101+##########################################################################!
!                                                                               !
!                           Poseidon_Input_Sources                              !
!                                                                               !
!###############################################################################!
SUBROUTINE Poseidon_Input_Sources_Part1_Native_Caller( Input_E, Input_Si  )



REAL(idp),  INTENT(IN), DIMENSION(  1:Caller_Quad_DOF,              &
                                    0:Num_R_Elements-1,             &
                                    0:Num_T_Elements-1,             &
                                    0:Num_P_Elements-1  )           ::  Input_E

REAL(idp),  INTENT(IN), DIMENSION(  1:Caller_Quad_DOF,              &
                                    0:Num_R_Elements-1,             &
                                    0:Num_T_Elements-1,             &
                                    0:Num_P_Elements-1,             &
                                    1:3         )                   ::  Input_Si





INTEGER                                                             ::  Local_Here
INTEGER                                                             ::  re, te, pe

IF ( Verbose_Flag ) CALL Run_Message('Receiving Part 1 XCFC Sources. Container : Fortran Array.')
CALL TimerStart(Timer_GR_SourceInput)
CALL TimerStart(Timer_GR_SourceInput_PartA)



DO pe = 0,Num_P_Elements-1
DO te = 0,Num_T_Elements-1
DO re = 0,Num_R_Elements-1

DO Local_Here = 1,Local_Quad_DOF

    Block_Source_E(Local_Here,re,te,pe) = DOT_PRODUCT( Translation_Matrix(:,Local_Here),    &
                                                       Input_E(:,re,te,pe)                  )

    
    Block_Source_Si(Local_Here,re,te,pe,1) = DOT_PRODUCT( Translation_Matrix(:,Local_Here), &
                                                          Input_Si(:,re,te,pe,1)            )


    Block_Source_Si(Local_Here,re,te,pe,2) = DOT_PRODUCT( Translation_Matrix(:,Local_Here), &
                                                          Input_Si(:,re,te,pe,2)            )

    Block_Source_Si(Local_Here,re,te,pe,3) = DOT_PRODUCT( Translation_Matrix(:,Local_Here), &
                                                          Input_Si(:,re,te,pe,3)            )

END DO  ! Local_Here

END DO ! re Loop
END DO ! te Loop
END DO ! pe Loop


CALL TimerStop(Timer_GR_SourceInput)
CALL TimerStop(Timer_GR_SourceInput_PartA)


END SUBROUTINE Poseidon_Input_Sources_Part1_Native_Caller




!+101+##########################################################################!
!                                                                               !
!                           Poseidon_Input_Sources                              !
!                                                                               !
!###############################################################################!
SUBROUTINE Poseidon_Input_Sources_XCFC_Native(  Input_E,          &
                                                Input_Si,         &
                                                Input_S,          &
                                                Input_NE,         &
                                                Input_NQ,         &
                                                Input_R_Quad,     &
                                                Input_T_Quad,     &
                                                Input_P_Quad,     &
                                                Input_xL          )

INTEGER,    INTENT(IN), DIMENSION(1:3)                                ::  Input_NE
INTEGER,    INTENT(IN), DIMENSION(1:3)                                ::  Input_NQ

REAL(idp),  INTENT(IN), DIMENSION(  1:Input_NQ(1)*Input_NQ(2)*Input_NQ(3),  &
                                    0:Input_NE(1)-1,                        &
                                    0:Input_NE(2)-1,                        &
                                    0:Input_NE(3)-1  )                      ::  Input_E

REAL(idp),  INTENT(IN), DIMENSION(  1:Input_NQ(1)*Input_NQ(2)*Input_NQ(3),  &
                                    0:Input_NE(1)-1,                        &
                                    0:Input_NE(2)-1,                        &
                                    0:Input_NE(3)-1,                        &
                                    1:3         )                           ::  Input_Si

REAL(idp),  INTENT(IN), DIMENSION(  1:Input_NQ(1)*Input_NQ(2)*Input_NQ(3),  &
                                    0:Input_NE(1)-1,                        &
                                    0:Input_NE(2)-1,                        &
                                    0:Input_NE(3)-1  )                      ::  Input_S


REAL(idp),  INTENT(IN), DIMENSION( 1:Input_NQ(1) )                  ::  Input_R_Quad
REAL(idp),  INTENT(IN), DIMENSION( 1:Input_NQ(2) )                  ::  Input_T_Quad
REAL(idp),  INTENT(IN), DIMENSION( 1:Input_NQ(3) )                  ::  Input_P_Quad
REAL(idp),  INTENT(IN), DIMENSION(2)                                ::  Input_xL


INTEGER                                             ::  Local_Here


REAL(idp), DIMENSION(:,:), ALLOCATABLE              ::  TransMat
INTEGER                                             ::  Their_DOF

INTEGER                                             ::  re, te, pe
INTEGER                                             ::  rd, td, pd


IF ( Verbose_Flag ) CALL Run_Message('Receiving XCFC Sources. Container : Fortran Array.')
CALL TimerStart(Timer_GR_SourceInput)
CALL TimerStart(Timer_GR_SourceInput_PartA)



! Define Interpolation Matrix
Their_DOF = Input_NQ(1)*Input_NQ(2)*Input_NQ(3)


ALLOCATE(TransMat(1:Their_DOF, 1:Local_Quad_DOF))

TransMat =  Create_Translation_Matrix(  Input_NQ,          &
                                        Input_xL,          &
                                        Input_R_Quad,    &
                                        Input_T_Quad,    &
                                        Input_P_Quad,    &
                                        Their_DOF,         &
                                        [Num_R_Quad_Points, Num_T_Quad_Points, Num_P_Quad_Points ],            &
                                        [xLeftLimit, xRightLimit ],            &
                                        Int_R_Locations,      &
                                        Int_R_Locations,      &
                                        Int_R_Locations,      &
                                        Local_Quad_DOF            )


DO pe = 0,Input_NE(3)-1
DO te = 0,Input_NE(2)-1
DO re = 0,Input_NE(1)-1

DO Local_Here = 1,Local_Quad_DOF

    Block_Source_E(Local_Here,re,te,pe) = DOT_PRODUCT( TransMat(Local_Here,:),      &
                                                       Input_E(:,re,te,pe)          )


    Block_Source_Si(Local_Here,re,te,pe,1) = DOT_PRODUCT( TransMat(:,Local_Here),   &
                                                          Input_Si(:,re,te,pe,1)    )


    Block_Source_Si(Local_Here,re,te,pe,2) = DOT_PRODUCT( TransMat(:,Local_Here),   &
                                                          Input_Si(:,re,te,pe,2)    )

    Block_Source_Si(Local_Here,re,te,pe,3) = DOT_PRODUCT( TransMat(:,Local_Here),   &
                                                          Input_Si(:,re,te,pe,3)    )

    Block_Source_S(Local_Here,re,te,pe) = DOT_PRODUCT( TransMat(:,Local_Here),      &
                                                       Input_S(:,re,te,pe)          )

END DO  ! Local_Here

END DO ! re Loop
END DO ! te Loop
END DO ! pe Loop


CALL TimerStop(Timer_GR_SourceInput)
CALL TimerStop(Timer_GR_SourceInput_PartA)



END SUBROUTINE Poseidon_Input_Sources_XCFC_Native








!+101+##########################################################################!
!                                                                               !
!                           Poseidon_Input_Sources                              !
!                                                                               !
!###############################################################################!
SUBROUTINE Poseidon_Input_Sources_XCFC_Native_Caller( Input_E, Input_Si, Input_S  )



REAL(idp),  INTENT(IN), DIMENSION(  1:Caller_Quad_DOF,              &
                                    0:Num_R_Elements-1,             &
                                    0:Num_T_Elements-1,             &
                                    0:Num_P_Elements-1  )           ::  Input_E

REAL(idp),  INTENT(IN), DIMENSION(  1:Caller_Quad_DOF,              &
                                    0:Num_R_Elements-1,             &
                                    0:Num_T_Elements-1,             &
                                    0:Num_P_Elements-1,             &
                                    1:3         )                   ::  Input_Si

REAL(idp),  INTENT(IN), DIMENSION(  1:Caller_Quad_DOF,              &
                                    0:Num_R_Elements-1,             &
                                    0:Num_T_Elements-1,             &
                                    0:Num_P_Elements-1  )           ::  Input_S



INTEGER                                                             ::  Local_Here
INTEGER                                                             ::  re, te, pe

IF ( Verbose_Flag ) CALL Run_Message('Receiving XCFC Sources. Container : Fortran Array.')
CALL TimerStart(Timer_GR_SourceInput)



DO pe = 0,Num_P_Elements-1
DO te = 0,Num_T_Elements-1
DO re = 0,Num_R_Elements-1

!PRINT*,Input_E(:,re,te,pe)
!PRINT*,Input_Si(:,re,te,pe,1)

DO Local_Here = 1,Local_Quad_DOF

    Block_Source_E(Local_Here,re,te,pe) = DOT_PRODUCT( Translation_Matrix(:,Local_Here),    &
                                                       Input_E(:,re,te,pe)                  )

    Block_Source_Si(Local_Here,re,te,pe,1) = DOT_PRODUCT( Translation_Matrix(:,Local_Here), &
                                                          Input_Si(:,re,te,pe,1)            )


    Block_Source_Si(Local_Here,re,te,pe,2) = DOT_PRODUCT( Translation_Matrix(:,Local_Here), &
                                                          Input_Si(:,re,te,pe,2)            )

    Block_Source_Si(Local_Here,re,te,pe,3) = DOT_PRODUCT( Translation_Matrix(:,Local_Here), &
                                                          Input_Si(:,re,te,pe,3)            )

    Block_Source_S(Local_Here,re,te,pe) = DOT_PRODUCT( Translation_Matrix(:,Local_Here),    &
                                                       Input_S(:,re,te,pe)                  )

END DO  ! Local_Here

END DO ! re Loop
END DO ! te Loop
END DO ! pe Loop


CALL TimerStop(Timer_GR_SourceInput)


END SUBROUTINE Poseidon_Input_Sources_XCFC_Native_Caller











END MODULE Source_Input_XCFC_Native_Module
