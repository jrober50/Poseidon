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

USE External_UST_Solution_Module, &
                ONLY :  Set_UST_Test_Params

USE Poseidon_Interface_Source_Input, &
                ONLY :  Poseidon_Input_Sources

USE Timer_Routines_Module, &
            ONLY :  TimerStart,     &
                    TimerSTop

USE Timer_Variables_Module, &
            ONLY :  Timer_Driver_SetSource_InitTest,        &
                    Timer_Driver_SetSource_SetSource,       &
                    Timer_Driver_SetSource_Scale



IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!     Driver_SetSource                                                			!
!                                                                               !
!###############################################################################!
SUBROUTINE Driver_SetSource( NE, NQ,                    &
                             Surface_RE,                &
                             rho_o,                     &
                             dx_c, x_e,                 &
                             R_Quad, T_Quad, P_Quad,    &
                             LeftLimit, RightLimit      )

INTEGER,    INTENT(IN), DIMENSION(3)                    ::  NE
INTEGER,    INTENT(IN), DIMENSION(3)                    ::  NQ

INTEGER,    INTENT(INOUT)                               ::  Surface_RE
REAL(idp),  INTENT(INOUT)                               ::  rho_o

REAL(idp),  INTENT(IN), DIMENSION(1:NE(1))              ::  dx_c
REAL(idp),  INTENT(IN), DIMENSION(0:NE(1))              ::  x_e

REAL(idp),  INTENT(IN), DIMENSION(1:NQ(1))              ::  R_Quad
REAL(idp),  INTENT(IN), DIMENSION(1:NQ(2))              ::  T_Quad
REAL(idp),  INTENT(IN), DIMENSION(1:NQ(3))              ::  P_Quad

REAL(idp),  INTENT(IN)                                  ::  LeftLimit
REAL(idp),  INTENT(IN)                                  ::  RightLimit



INTEGER                                                 ::  re
INTEGER                                                 ::  Num_DOF

REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              ::  Density


Num_DOF = NQ(1)*NQ(2)*NQ(3)

ALLOCATE( Density(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1 )  )


CALL TimerStart( Timer_Driver_SetSource_InitTest )


CALL Set_UST_Test_Params(   rho_o,               &
                            x_e(Surface_RE)      )


Density = 0.0_idp
DO re = 0,NE(1)-1

    IF ( re .LE. Surface_RE-1 ) THEN
        Density(:,re,:,:) = rho_o
    END IF
    
END DO ! re

CALL TimerStop( Timer_Driver_SetSource_InitTest )







CALL TimerStart( Timer_Driver_SetSource_SetSource )

CALL Poseidon_Input_Sources(    Density,                &
                                NE,                     &
                                NQ,                     &
                                R_Quad,                 &
                                T_Quad,                 &
                                P_Quad,                 &
                                [LeftLimit, RightLimit] )

CALL TimerStop( Timer_Driver_SetSource_SetSource )




END SUBROUTINE Driver_SetSource









END MODULE Driver_SetSource_Module
