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
            ONLY :  Centimeter, gram, C_Square

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag

USE Poseidon_Message_Routines_Module, &
            ONLY :  Driver_Init_Message

USE Poseidon_IO_Parameters, &
            ONLY :  Poseidon_Results_Dir

USE Poseidon_Interface_Source_Input, &
            ONLY :  Poseidon_Input_Sources

USE Variables_IO, &
            ONLY :  File_Suffix

USE Variables_External, &
            ONLY :  HCT_Rhoo,       &
                    HCT_Star_Radius

USE Poseidon_File_Routines_Module, &
            ONLY :  Open_New_File

USE Driver_SetSource_Module, &
            ONLY :  Driver_SetSource

USE Poseidon_Interface_Run, &
            ONLY :  Poseidon_Run

USE Poseidon_Interface_Return_Routines, &
            ONLY :  Poseidon_Return_Conformal_Factor

USE Maps_Quadrature, &
            ONLY :  Quad_Map

USE IO_Print_Results, &
            ONLY :  Print_Results

USE Driver_SetGuess_Module, &
            ONLY :  Driver_SetGuess

IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!     Driver_SetSource                                                            !
!                                                                               !
!###############################################################################!
SUBROUTINE Driver_ConFactor_Loop( NE, NQ,                    &
                                  dx_c, x_e, y_e,            &
                                  R_Quad, T_Quad, P_Quad,    &
                                  LeftLimit, RightLimit,     &
                                  myID,                      &
                                  Alpha,                     &
                                  Star_Radius,               &
                                  Tolerance,                 &
                                  Iter_Max                   )

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

REAL(idp), INTENT(IN)                                   ::  Alpha
REAL(idp), INTENT(IN)                                   ::  Star_Radius

REAL(idp), INTENT(IN)                                   ::  Tolerance

INTEGER,   INTENT(IN)                                   ::  Iter_Max



INTEGER                                                 ::  re, te, pe
INTEGER                                                 ::  rq, tq, pq

INTEGER                                                 ::  Num_DOF
INTEGER                                                 ::  Iter
INTEGER                                                 ::  Here


REAL(idp)                                               ::  Max_Difference

REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              ::  Difference
REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              ::  ConFactor_New
REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              ::  ConFactor_Old
REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              ::  Local_E
REAL(idp), DIMENSION(:,:,:,:), ALLOCATABLE              ::  Local_S
REAL(idp), DIMENSION(:,:,:,:,:), ALLOCATABLE            ::  Local_Si

REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  Cur_R_Locs

INTEGER                                                 ::  HCT_Fileid
CHARACTER(LEN = 100)                                    ::  HCT_Filename


REAL(idp)                                               ::  Beta
REAL(idp)                                               ::  C
REAL(idp)                                               ::  rho_o
REAL(idp)                                               ::  uaR
REAL(idp)                                               ::  fofalpha

LOGICAL                                                 ::  Flag

Num_DOF = NQ(1)*NQ(2)*NQ(3)

ALLOCATE( Difference(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1 )    )
ALLOCATE( ConFactor_New(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1 ) )
ALLOCATE( ConFactor_Old(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1 ) )
ALLOCATE( Local_E(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1 )       )
ALLOCATE( Local_S(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1 )       )
ALLOCATE( Local_Si(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1, 1:3)  )

ALLOCATE( Cur_R_Locs(1:NQ(1)) )



IF ( Verbose_Flag ) CALL Driver_Init_Message('Initializing the Hamiltonian constraint test source.')




CALL Poseidon_Return_Conformal_Factor( NE, NQ,                 &
                                        R_Quad,                 &
                                        T_Quad,                 &
                                        P_Quad,                 &
                                        LeftLimit,              &
                                        RightLimit,             &
                                        ConFactor_New        )




Iter = 0
Flag = .TRUE.
DO WHILE ( Flag )


    DO re = 1,NE(1)
    DO te = 1,NE(2)
    DO pe = 1,NE(3)
        Cur_R_Locs(:) = dx_c(re)*(R_Quad(:) + LeftLimit)+  x_e(re)

        DO pq = 1,NQ(3)
        DO tq = 1,NQ(2)
        DO rq = 1,NQ(1)


            Here = Quad_Map(rq, tq, pq, NQ(1), NQ(2),NQ(3))


            IF ( cur_r_locs(rq) .LE. HCT_Star_Radius ) THEN
                Local_E(Here, re-1, te-1, pe-1) = HCT_Rhoo*ConFactor_New(Here,re-1,te-1,pe-1)**6
!                Local_E(Here, re-1, te-1, pe-1) = HCT_Rhoo
!                PRINT*,HCT_Rhoo*ConFactor_New(Here,re-1,te-1,pe-1)**6,      &
!                        HCT_Rhoo,ConFactor_New(Here,re-1,te-1,pe-1)**6
            ELSE
                Local_E(Here, re-1, te-1, pe-1) = 0.0_idp
            END IF

        END DO ! rq
        END DO ! tq
        END DO ! pq
    END DO ! pe
    END DO ! te
    END DO ! re

    Local_S  = 0.0_idp
    Local_Si = 0.0_idp








    CALL Poseidon_Input_Sources(Local_E,                &
                                Local_Si,               &
                                Local_S,                &
                                NE,                     &
                                NQ,                     &
                                R_Quad,                 &
                                T_Quad,                 &
                                P_Quad,                 &
                                [LeftLimit, RightLimit] )


    Call Poseidon_Run()


    CALL Poseidon_Return_Conformal_Factor( NE, NQ,                 &
                                            R_Quad,                 &
                                            T_Quad,                 &
                                            P_Quad,                 &
                                            LeftLimit,              &
                                            RightLimit,             &
                                            ConFactor_New        )

    Difference = abs(ConFactor_New - ConFactor_Old)
    Max_Difference = maxval(Difference)
    ConFactor_Old = ConFactor_New

    IF ( .TRUE. ) THEN
!    IF ( Verbose_Flag ) THEN
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



IF ( .FALSE. ) THEN

    CALL Output_Poseidon_Sources_3D( Local_E, Local_S, Local_Si,    &
                                     NE(1), NE(2), NE(3),           &
                                     NQ(1), NQ(2), NQ(3),           &
                                     R_Quad, T_Quad, P_Quad,        &
                                     LeftLimit, RightLimit          )

END IF



DEALLOCATE( Cur_R_Locs )


END SUBROUTINE Driver_ConFactor_Loop









END MODULE Driver_ConFactor_Loop_Module

