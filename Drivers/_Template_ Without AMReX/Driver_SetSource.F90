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

USE SelfSimilar_Module, &
                ONLY :  Initialize_Yahil_Sources

USE Source_Input_Module, &
                ONLY :  Poseidon_Input_Sources





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
                             LeftLimit, RightLimit,     &
                             Yahil_Params,              &
                             Solver_Type                )

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

REAL(idp), INTENT(IN), DIMENSION(3)                     ::  Yahil_Params
INTEGER,   INTENT(IN)                                   ::  Solver_Type

REAL(idp)                                               ::  Psi_Holder


INTEGER                                                 ::  nLevels
INTEGER                                                 ::  nVars_Source
INTEGER                                                 ::  nghost
INTEGER, DIMENSION(3)                                   ::  NE_Lower
INTEGER, DIMENSION(3)                                   ::  NE_Upper

INTEGER                                                 ::  re, te, pe
INTEGER                                                 ::  rq, tq, pq
INTEGER                                                 ::  here, lvl
INTEGER                                                 ::  there, var


INTEGER                                                 ::  Num_DOF


REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              ::  Local_E
REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              ::  Local_S
REAL(idp), DIMENSION(:,:,:,:,:), ALLOCATABLE            ::  Local_Si


Num_DOF = NQ(1)*NQ(2)*NQ(3)

ALLOCATE( Local_E(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1 )       )
ALLOCATE( Local_S(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1 )       )
ALLOCATE( Local_Si(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1, 1:3)  )


CALL Initialize_Yahil_Sources(  Yahil_Params(1),                &
                                Yahil_Params(2),                &
                                Yahil_Params(3),                &
                                0.0_idp,                        &
                                NQ, R_Quad, T_Quad,             &
                                NE(1), NE(2), NE(3),            &
                                dx_c, x_e, y_e,                 &
                                Local_E, Local_S, Local_Si      )



IF ( Verbose_Flag ) THEN
    WRITE(*,'(A)')"-Inputing Sources"
END IF

!CALL Poseidon_Input_Sources(myID_Poseidon,                      &
!                            myID_Poseidon,                      &
!                            myID_Poseidon,                      &
!                            Local_E, Local_S, Local_Si,         &
!                            NE(1), NE(2), NE(3),                &
!                            NQ(1), NQ(2), NQ(3),                &
!                            R_Quad, T_Quad, P_Quad,             &
!                            LeftLimit, RightLimit               )







IF ( .FALSE. ) THEN

    CALL Output_Poseidon_Sources_3D( Local_E, Local_S, Local_Si,    &
                                     NE(1), NE(2), NE(3),           &
                                     NQ(1), NQ(2), NQ(3),           &
                                     R_Quad, T_Quad, P_Quad,        &
                                     LeftLimit, RightLimit          )

END IF




END SUBROUTINE Driver_SetSource









END MODULE Driver_SetSource_Module
