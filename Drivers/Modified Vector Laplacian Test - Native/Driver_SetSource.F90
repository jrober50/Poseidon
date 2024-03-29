   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Driver_SetSource_Module                                              !##!
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
            ONLY :  Driver_Init_Message

USE Poseidon_Interface_Source_Input, &
            ONLY :  Poseidon_Input_Sources

USE Timer_Routines_Module, &
            ONLY :  TimerStart,     &
                    TimerSTop

USE Timer_Variables_Module, &
            ONLY :  Timer_Driver_SetSource_InitTest,        &
                    Timer_Driver_SetSource_SetSource,       &
                    Timer_Driver_SetSource_Scale

USE IO_Output_Sources_Module, &
            ONLY :  Output_Poseidon_Sources_3D

USE Flags_IO_Module, &
            ONLY :  lPF_IO_Flags,       &
                    iPF_IO_Write_Sources

IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!     Driver_SetSource                                                			!
!                                                                               !
!###############################################################################!
SUBROUTINE Driver_SetSource( NE, NQ,                    &
                             dx_c, x_e, y_e,            &
                             R_Quad, T_Quad, P_Quad,    &
                             LeftLimit, RightLimit,   &
                             myID                       )

INTEGER, INTENT(IN), DIMENSION(3)                       ::  NE
INTEGER, INTENT(IN), DIMENSION(3)                       ::  NQ

REAL(idp), INTENT(IN), DIMENSION(1:NE(1))               ::  dx_c
REAL(idp), INTENT(IN), DIMENSION(0:NE(1))               ::  x_e
REAL(idp), INTENT(IN), DIMENSION(0:NE(2))               ::  y_e

REAL(idp), INTENT(IN), DIMENSION(1:NQ(1))               ::  R_Quad
REAL(idp), INTENT(IN), DIMENSION(1:NQ(2))               ::  T_Quad
REAL(idp), INTENT(IN), DIMENSION(1:NQ(3))               ::  P_Quad

REAL(idp), INTENT(IN)                                   ::  LeftLimit
REAL(idp), INTENT(IN)                                   ::  RightLimit

INTEGER,   INTENT(IN)                                   ::  myID

INTEGER                                                 ::  Num_DOF
INTEGER                                                 ::  Here

INTEGER                                                 ::  re, te, pe
INTEGER                                                 ::  rd, td, pd
REAL(idp), DIMENSION(1:NQ(1))                           ::  Cur_R_Locs

REAL(idp)                                               ::  Psi_Holder
REAL(idp)                                               ::  Psi_Power

REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              ::  Local_E
REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              ::  Local_S
REAL(idp), DIMENSION(:,:,:,:,:), ALLOCATABLE            ::  Local_Si


IF ( Verbose_Flag ) CALL Driver_Init_Message('Initializing the Modified Vector Laplacian test source.')


Num_DOF = NQ(1)*NQ(2)*NQ(3)

ALLOCATE( Local_E(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1 )       )
ALLOCATE( Local_S(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1 )       )
ALLOCATE( Local_Si(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1, 1:3)  )



Local_E  = 0.0_idp
Local_S  = 0.0_idp
Local_Si = 0.0_idp

CALL TimerStart( Timer_Driver_SetSource_SetSource )




CALL Poseidon_Input_Sources(Local_E,                &
                            Local_Si,               &
                            Local_S,                &
                            NE,                     &
                            NQ,                     &
                            R_Quad,                 &
                            T_Quad,                 &
                            P_Quad,                 &
                            [LeftLimit, RightLimit] )



IF ( lPF_IO_Flags(iPF_IO_Write_Sources) ) THEN

    CALL Output_Poseidon_Sources_3D( Local_E, Local_S, Local_Si,     &
                                     NE(1), NE(2), NE(3),            &
                                     NQ(1), NQ(2), NQ(3),            &
                                     R_Quad, T_Quad, P_Quad,         &
                                     LeftLimit, RightLimit         )

 
END IF

CALL TimerStop( Timer_Driver_SetSource_SetSource )

DEALLOCATE( Local_E, Local_S, Local_Si )

END SUBROUTINE Driver_SetSource









END MODULE Driver_SetSource_Module
