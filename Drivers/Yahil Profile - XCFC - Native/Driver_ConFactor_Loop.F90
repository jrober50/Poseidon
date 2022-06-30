   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Driver_ConFactor_Loop_Module                                          !##!
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

USE Poseidon_Numbers_Module, &
            ONLY :  pi

USE Poseidon_Units_Module, &
            ONLY :  Centimeter

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag

USE Poseidon_Message_Routines_Module, &
            ONLY :  Driver_Init_Message

USE Poseidon_Interface_Source_Input, &
            ONLY :  Poseidon_Input_Sources

USE External_Yahil_Profile_Module, &
            ONLY :  Initialize_Yahil_Sources

USE Poseidon_Interface_Run, &
            ONLY :  Poseidon_Run

USE Poseidon_Interface_Return_Routines, &
            ONLY :  Poseidon_Return_Conformal_Factor

USE IO_Print_Results, &
            ONLY :  Print_Results


IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!     Driver_ConFactor_Loop                                                            !
!                                                                               !
!###############################################################################!
SUBROUTINE Driver_ConFactor_Loop(   NE, NQ,                    &
                                    dx_c, x_e, y_e,            &
                                    R_Quad, T_Quad, P_Quad,    &
                                    Left_Limit, Right_Limit,   &
                                    Solver_Type,               &
                                    myID,                      &
                                    Yahil_Params,              &
                                    Tolerance,                 &
                                    Iter_Max                    )

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

REAL(idp), INTENT(IN)                                   ::  Tolerance

INTEGER,   INTENT(IN)                                   ::  Iter_Max

INTEGER                                                 ::  Num_DOF
INTEGER                                                 ::  Here

INTEGER                                                 ::  re, te, pe
INTEGER                                                 ::  rq, tq, pq

REAL(idp)                                               ::  Psi_Holder
REAL(idp)                                               ::  Psi_Power


REAL(idp)                                               ::  Max_Difference

REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              ::  Difference
REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              ::  ConFactor_New
REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              ::  ConFactor_Old
REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              ::  Local_E
REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              ::  Local_S
REAL(idp), DIMENSION(:,:,:,:,:), ALLOCATABLE            ::  Local_Si

REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              ::  Scaled_E
REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              ::  Scaled_S
REAL(idp), DIMENSION(:,:,:,:,:), ALLOCATABLE            ::  Scaled_Si

REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  Cur_R_Locs

LOGICAL                                                 ::  Flag

INTEGER                                                 ::  Iter

Num_DOF = NQ(1)*NQ(2)*NQ(3)

ALLOCATE( Difference(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1 )    )
ALLOCATE( ConFactor_New(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1 ) )
ALLOCATE( ConFactor_Old(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1 ) )
ALLOCATE( Local_E(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1 )       )
ALLOCATE( Local_S(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1 )       )
ALLOCATE( Local_Si(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1, 1:3)  )

ALLOCATE( Scaled_E(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1 )       )
ALLOCATE( Scaled_S(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1 )       )
ALLOCATE( Scaled_Si(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1, 1:3)  )

ALLOCATE( Cur_R_Locs(1:NQ(1)) )



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


IF ( Verbose_Flag ) CALL Driver_Init_Message('Starting Yahil Conformal Factor Loop.')






CALL Poseidon_Return_Conformal_Factor( NE, NQ,                 &
                                        R_Quad,                 &
                                        T_Quad,                 &
                                        P_Quad,                 &
                                        Left_Limit,              &
                                        Right_Limit,             &
                                        ConFactor_Old       )




Iter = 0
Flag = .TRUE.
DO WHILE ( Flag )

    DO re = 0,NE(1)-1
    DO te = 0,NE(2)-1
    DO pe = 0,NE(3)-1
    DO pq = 1,NQ(3)
    DO tq = 1,NQ(2)
    DO rq = 1,NQ(1)


        Here = (pq-1)*NQ(1)*NQ(2)       &
             + (tq-1)*NQ(1)             &
             + rq

        Scaled_E(Here,re,te,pe) = ConFactor_Old(Here,re,te,pe)**6 * Local_E(Here,re,te,pe)
        Scaled_Si(Here,re,te,pe,1) = ConFactor_Old(Here,re,te,pe)**6 * Local_Si(Here,re,te,pe,1)
        Scaled_Si(Here,re,te,pe,2) = ConFactor_Old(Here,re,te,pe)**6 * Local_Si(Here,re,te,pe,2)
        Scaled_Si(Here,re,te,pe,3) = ConFactor_Old(Here,re,te,pe)**6 * Local_Si(Here,re,te,pe,3)
        Scaled_S(Here,re,te,pe) = ConFactor_Old(Here,re,te,pe)**6 * Local_S(Here,re,te,pe)

        PRINT*,Scaled_E(Here,re,te,pe),ConFactor_Old(Here,re,te,pe)**6, Local_E(Here,re,te,pe)
    END DO ! rq
    END DO ! tq
    END DO ! pq
    END DO ! pe
    END DO ! te
    END DO ! re
    






    IF ( Verbose_Flag ) THEN
        WRITE(*,'(A)')"-Inputing Sources"
    END IF

    CALL Poseidon_Input_Sources(Scaled_E,                &
                                Scaled_Si,               &
                                Scaled_S,                &
                                NE,                     &
                                NQ,                     &
                                R_Quad,                 &
                                T_Quad,                 &
                                P_Quad,                 &
                                [Left_Limit, Right_Limit] )


    Call Poseidon_Run()


    CALL Poseidon_Return_Conformal_Factor( NE, NQ,                 &
                                            R_Quad,                 &
                                            T_Quad,                 &
                                            P_Quad,                 &
                                            Left_Limit,              &
                                            Right_Limit,             &
                                            ConFactor_New        )

    Difference = abs(ConFactor_New - ConFactor_Old)
    Max_Difference = maxval(Difference)
    ConFactor_Old = ConFactor_New

    IF ( Verbose_Flag ) THEN
        WRITE(*,'(A,I3.3)') "End of ConFactor Iter ",Iter
        PRINT*,"Max_Difference",Max_Difference, Tolerance
    END IF


    Iter = Iter + 1
    IF ( Iter > Iter_Max ) THEN
        Flag = .FALSE.
    ELSEIF ( Max_Difference .LE. Tolerance) THEN
        Flag = .FALSE.
    ELSE
        Flag = .TRUE.
    END IF

END DO ! While Iter < Iter_Max





DEALLOCATE( Cur_R_Locs )


END SUBROUTINE Driver_ConFactor_Loop









END MODULE Driver_ConFactor_Loop_Module

