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

USE Poseidon_Numbers_Module, &
            ONLY :  pi

USE Poseidon_Units_Module, &
            ONLY :  Centimeter

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag

USE Poseidon_IO_Parameters, &
            ONLY :  Poseidon_Results_Dir

USE SelfSimilar_Module, &
            ONLY :  Initialize_Yahil_Sources

USE Source_Input_Module, &
            ONLY :  Poseidon_Input_Sources

USE Variables_IO, &
            ONLY :  File_Suffix

USE Poseidon_File_Routines_Module, &
            ONLY :  Open_New_File


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
                             myID,                      &
                             Alpha,                     &
                             Star_Radius                )

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

REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  Cur_R_Locs

INTEGER                                                 ::  HCT_Fileid
CHARACTER(LEN = 100)                                    ::  HCT_Filename


REAL(idp)                                               ::  fofalpha
REAL(idp)                                               ::  rho_o



Num_DOF = NQ(1)*NQ(2)*NQ(3)


ALLOCATE( Local_E(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1 )       )
ALLOCATE( Local_S(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1 )       )
ALLOCATE( Local_Si(1:Num_DOF, 0:NE(1)-1, 0:NE(2)-1, 0:NE(3)-1, 1:3)  )

ALLOCATE( Cur_R_Locs(1:NQ(1)) )




IF ( Verbose_Flag ) THEN
    WRITE(*,'(A)')"-Creating Sources."
END IF

HCT_Fileid = 4242

WRITE(HCT_Filename,'(A,A,A,A)') Poseidon_Results_Dir,"HCT_Params_",TRIM(File_Suffix),".out"
Call Open_New_File(HCT_Filename,HCT_Fileid)
WRITE(HCT_Fileid,'(2ES22.15)') Alpha, Star_Radius/Centimeter
CLOSE( HCT_Fileid)




fofalpha            =  Alpha**5/(1.0_idp+Alpha*Alpha)**3

rho_o               =  (3.0_idp/(2.0_idp*pi*Star_Radius*Star_Radius) )*fofalpha*fofalpha



DO re = 1,NE(1)
DO te = 1,NE(2)
DO pe = 1,NE(3)
    Cur_R_Locs(:) = dx_c(re)*(R_Quad(:) + LeftLimit)+  x_e(re)

    DO pq = 1,NQ(3)
    DO tq = 1,NQ(2)
    DO rq = 1,NQ(1)


        Here = (pq-1)*NQ(1)*NQ(2)       &
             + (tq-1)*NQ(1)             &
             + rq

        IF ( cur_r_locs(rq) .LE. Star_Radius ) THEN

            Local_E(Here, re-1, te-1, pe-1) = rho_o
!            Print*,Local_E(rq, re-1, 0, 0)
        ELSE
            Local_E(Here, re-1, te-1, pe-1) = 0.0_idp
        END IF
!            PRINT*,Cur_r_locs(rq), STar_radius,Local_E(rq, re-1, 0, 0)

    END DO ! rq
    END DO ! tq
    END DO ! pq
END DO ! pe
END DO ! te
END DO ! re

Local_S  = 0.0_idp
Local_Si = 0.0_idp




IF ( Verbose_Flag ) THEN
    WRITE(*,'(A)')"-Inputing Sources"
END IF

CALL Poseidon_Input_Sources(myID,                           &
                            myID,                           &
                            myID,                           &
                            Local_E, Local_S, Local_Si,     &
                            NE(1), NE(2), NE(3),            &
                            NQ(1), NQ(2), NQ(3),            &
                            R_Quad, T_Quad, P_Quad,         &
                            LeftLimit, RightLimit           )




IF ( .FALSE. ) THEN

    CALL Output_Poseidon_Sources_3D( Local_E, Local_S, Local_Si,    &
                                     NE(1), NE(2), NE(3),           &
                                     NQ(1), NQ(2), NQ(3),           &
                                     R_Quad, T_Quad, P_Quad,        &
                                     LeftLimit, RightLimit          )

END IF



DEALLOCATE( Cur_R_Locs )


END SUBROUTINE Driver_SetSource









END MODULE Driver_SetSource_Module
