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

USE Poseidon_Units_Module, &
            ONLY :  C_Square,   &
                    Centimeter

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag

USE Poseidon_Message_Routines_Module, &
            ONLY :  Driver_Init_Message

USE External_Yahil_Profile_Module, &
            ONLY :  Initialize_Yahil_Sources

USE Poseidon_Interface_Source_Input, &
            ONLY :  Poseidon_Input_Sources

USE Variables_Functions, &
            ONLY :  Potential_Solution

USE Timer_Routines_Module, &
            ONLY :  TimerStart,     &
                    TimerSTop

USE Timer_Variables_Module, &
            ONLY :  Timer_Driver_SetSource_InitTest,        &
                    Timer_Driver_SetSource_SetSource,       &
                    Timer_Driver_SetSource_Scale

USE Maps_Quadrature, &
            ONLY :  Quad_Map

USE IO_Output_Sources_Module, &
            ONLY :  Output_Poseidon_Sources_3D

USE Flags_IO_Module, &
            ONLY :  lPF_IO_Flags,                       &
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
                             Left_Limit, Right_Limit,   &
                             Solver_Type,               &
                             myID,                      &
                             Yahil_Params               )

INTEGER, INTENT(IN), DIMENSION(3)                       ::  NE
INTEGER, INTENT(IN), DIMENSION(3)                       ::  NQ

REAL(idp), INTENT(IN), DIMENSION(1:NE(1))               ::  dx_c
REAL(idp), INTENT(IN), DIMENSION(0:NE(1))               ::  x_e
REAL(idp), INTENT(IN), DIMENSION(0:NE(2))               ::  y_e

REAL(idp), INTENT(IN), DIMENSION(1:NQ(1))               ::  R_Quad
REAL(idp), INTENT(IN), DIMENSION(1:NQ(2))               ::  T_Quad
REAL(idp), INTENT(IN), DIMENSION(1:NQ(3))               ::  P_Quad

REAL(idp), INTENT(IN)                                   ::  Left_Limit
REAL(idp), INTENT(IN)                                   ::  Right_Limit

INTEGER,   INTENT(IN)                                   ::  Solver_Type
INTEGER,   INTENT(IN)                                   ::  myID
REAL(idp), INTENT(IN), DIMENSION(4)                     ::  Yahil_Params

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


Num_DOF = NQ(1)*NQ(2)*NQ(3)

ALLOCATE( Local_E(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1 )       )
ALLOCATE( Local_S(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1 )       )
ALLOCATE( Local_Si(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1, 1:3)  )



CALL TimerStart( Timer_Driver_SetSource_InitTest )

CALL Initialize_Yahil_Sources( Yahil_Params(1),                 &
                                Yahil_Params(2),                &
                                Yahil_Params(3),                &
                                Yahil_Params(4),                &
                                NQ, R_Quad, T_Quad,             &
                                NE(1), NE(2), NE(3),            &
                                dx_c*Centimeter,                &
                                x_e*Centimeter,                 &
                                y_e,                            &
                                Local_E, Local_S, Local_Si      )





CALL TimerStop( Timer_Driver_SetSource_InitTest )

CALL TimerStart( Timer_Driver_SetSource_Scale )

IF ( Solver_Type == 3 ) THEN


    DO pe = 1,NE(3)
    DO te = 1,NE(2)
    DO re = 1,NE(1)


    Cur_R_Locs = dx_c(re)*(R_Quad(:) - Left_Limit) + x_e(re-1)
    Cur_R_Locs = Cur_R_Locs*Centimeter

    DO rd = 1,NQ(1)

    Psi_Holder = 1.0_idp    &
               - 0.5_idp*Potential_Solution(cur_r_locs(rd), 0.0_idp, 0.0_idp)/C_Square

    Psi_Power  = Psi_Holder**6

    DO pd = 1,NQ(3)
    DO td = 1,NQ(2)

        here = Quad_Map(rd,td,pd)

        Local_E(Here,re-1,te-1,pe-1) = Local_E(Here,re-1,te-1,pe-1)*Psi_Power
        Local_S(Here,re-1,te-1,pe-1) = Local_S(Here,re-1,te-1,pe-1)*Psi_Power
        Local_Si(Here,re-1,te-1,pe-1,1:3) = Local_Si(Here,re-1,te-1,pe-1,1:3)*Psi_Power
        

    END DO  ! pd
    END DO  ! td
    END DO  ! rd

    END DO ! re
    END DO ! te
    END DO ! pe



 END IF

CALL TimerStop( Timer_Driver_SetSource_Scale )




CALL TimerStart( Timer_Driver_SetSource_SetSource )


!CALL Poseidon_Input_Sources(    myID, myID, myID,               &
!                                Local_E, Local_S, Local_Si,     &
!                                NE(1), NE(2), NE(3),            &
!                                NQ(1), NQ(2), NQ(3),            &
!                                R_Quad, T_Quad, P_Quad,         &
!                                Left_Limit, Right_Limit         )


!CALL Poseidon_Input_Sources( Local_E,                    &
!                             Local_Si,                   &
!                             Local_S,                    &
!                             NE,                         &
!                             NQ,                         &
!                             R_Quad,                     &
!                             T_Quad,                     &
!                             P_Quad,                     &
!                             [Left_Limit, Right_Limit]   )

CALL Poseidon_Input_Sources( Local_E,                    &
                            Local_Si,                   &
                            Local_S                     )


IF ( lPF_IO_Flags(iPF_IO_Write_Sources) ) THEN

    CALL Output_Poseidon_Sources_3D( Local_E, Local_S, Local_Si,     &
                                     NE(1), NE(2), NE(3),            &
                                     NQ(1), NQ(2), NQ(3),            &
                                     R_Quad, T_Quad, P_Quad,         &
                                     Left_Limit, Right_Limit         )

 
END IF


CALL TimerStop( Timer_Driver_SetSource_SetSource )

DEALLOCATE( Local_E, Local_S, Local_Si )

END SUBROUTINE Driver_SetSource









END MODULE Driver_SetSource_Module
