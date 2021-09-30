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

use amrex_base_module

USE amrex_box_module,       ONLY: &
  amrex_box
USE amrex_boxarray_module,  ONLY: &
  amrex_boxarray,         &
  amrex_boxarray_build,   &
  amrex_boxarray_destroy
USE amrex_distromap_module, ONLY: &
  amrex_distromap,       &
  amrex_distromap_build, &
  amrex_distromap_destroy
USE amrex_multifab_module,  ONLY: &
  amrex_multifab, &
  amrex_multifab_build




USE Poseidon_Kinds_Module, &
                ONLY :  idp

USE Poseidon_Parameters, &
                ONLY :  Verbose_Flag

USE Units_Module, &
                ONLY :  C_Square

USE Variables_MPI, &
                ONLY :  myID_Poseidon,      &
                        nPROCS_Poseidon,    &
                        Poseidon_Comm_World

USE Variables_Functions, &
                ONLY :  Potential_Solution

USE SelfSimilar_Module, &
                ONLY :  Initialize_Yahil_Sources

USE Source_Input_Module, &
                ONLY :  Poseidon_Input_Sources

USE Source_Input_AMReX, &
                ONLY :  Poseidon_Input_Sources_AMREX

USE Poseidon_AMReX_Utilities_Module, &
                ONLY :  PackSources_AMReX

USE Poseidon_MPI_Utilities_Module, &
                ONLY :  STOP_MPI,               &
                        MPI_Master_Print,       &
                        MPI_All_Print

USE Variables_AMReX_Multifabs, &
            ONLY :  MaxGridSizeX,       &
                    MF_Source,          &
                    BA_Source,          &
                    DM_Source,          &
                    Geom_Source

USE Variables_AMReX_Source, &
            ONLY :  Source_PTR

USE MPI

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


!TYPE(amrex_box)                                         ::  Domain
TYPE(amrex_box)                                         ::  Box


INTEGER                                                 ::  nLevels
INTEGER                                                 ::  nVars_Source
INTEGER                                                 ::  nghost
INTEGER, DIMENSION(3)                                   ::  NE_Lower
INTEGER, DIMENSION(3)                                   ::  NE_Upper

INTEGER                                                 ::  re, te, pe
INTEGER                                                 ::  rq, tq, pq
INTEGER                                                 ::  here, lvl
INTEGER                                                 ::  there, var

INTEGER, DIMENSION(3)                               :: ELo, EHi
TYPE(amrex_mfiter)                                  :: mfi
REAL(idp), CONTIGUOUS, POINTER                      :: p(:,:,:,:)


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
                                Local_E, Local_S, Local_Si                   )



IF ( Solver_Type == 3 ) THEN

    IF ( Verbose_Flag ) THEN
        WRITE(*,'(A)')"-Scaling Sources"
    END IF

    DO pe = 1,NE(3)
    DO te = 1,NE(2)
    DO re = 1,NE(1)
    DO pq = 1,NQ(3)
    DO tq = 1,NQ(2)
    DO rq = 1,NQ(1)

        here = (pq-1)*NQ(2)*NQ(1)   &
             + (tq-1)*NQ(1)                &
             + rq
        Psi_Holder = 1.0_idp    &
                - 0.5_idp*Potential_Solution(x_e(re-1), 0.0_idp, 0.0_idp)/C_Square

        Local_E(Here,re-1,te-1,pe-1) = Local_E(Here,re-1,te-1,pe-1)*Psi_Holder**6
        Local_S(Here,re-1,te-1,pe-1) = Local_S(Here,re-1,te-1,pe-1)*Psi_Holder**6
        Local_Si(Here,re-1,te-1,pe-1,1:3) = Local_Si(Here,re-1,te-1,pe-1,1:3)*Psi_Holder**6

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
END IF



IF ( Verbose_Flag ) THEN
    WRITE(*,'(A)')"-Creating AMReX Variables"
END IF
nghost = 0
NE_Lower = [ 1,1,1 ]
NE_Upper = [ NE(1),NE(2),NE(3) ]
nVars_Source = 5
nLevels      = 1

!ALLOCATE( MF_Source(0:nLevels-1) )
!ALLOCATE( Geom_Source(0:nLevels-1) )
!ALLOCATE( DM_Source(0:nLevels-1) )
!ALLOCATE( BA_Source(0:nLevels-1) )

Box = amrex_box(NE_Lower, NE_upper)

DO lvl = 0,nLevels-1
    CALL amrex_boxarray_build(BA_Source(lvl), Box)
END DO
DO lvl = 0,nLevels-1
    CALL BA_Source(lvl)%maxSize( MaxGridSizeX )
END DO

DO lvl = 0,nLevels-1
    CALL amrex_geometry_build( Geom_Source(lvl), Box )
    CALL amrex_distromap_build(DM_Source(lvl),BA_Source(lvl))
END DO

DO lvl = 0,nLevels-1
    CALL amrex_multifab_build(  MF_Source(lvl),        &
                                BA_Source(lvl),          &
                                DM_Source(lvl),         &
                                Num_DOF*nVars_Source,   &
                                nghost                  )
    CALL MF_Source(lvl) % SetVal(0.0_idp)
END DO ! lvl


!Print*,"Calling SetSource_Parallel",myID_Poseidon
CALL SetSource_Parallel( Num_DOF, nVars_Source, nLevels, &
                         NE_Lower, NE_Upper,             &
                         Local_E, Local_S, Local_Si,     &
                         MF_Source(:)                   )







DO lvl = 0,nLevels-1
    CALL amrex_mfiter_build(mfi, MF_Source(lvl), tiling = .false. )
    DO WHILE(mfi%next())

        CALL MPI_All_Print("In mfiter loop B")

        Source_Ptr => MF_Source(lvl)%dataPtr(mfi)
        Box = mfi%tilebox()

        ELo = Box%lo
        EHi = Box%hi


!        DO PE = ELo(3), EHi(3)
!        DO TE = ELo(2), EHi(2)
!        DO RE = ELo(1), EHi(1)
!
!
!            DO Var = 3,5
!                Here  = (Var-1)*Num_DOF+1
!                There = Var*Num_DOF
!                PRINT*,"myID ",myID_Poseidon," - ",Source_PTR(RE,TE,PE,Here:There)
!            END DO
!
!
!
!        END DO
!        END DO
!        END DO

    END DO
    CALL amrex_mfiter_destroy(mfi)
END DO ! lvl





IF ( Verbose_Flag ) THEN
    WRITE(*,'(A)')"-Inputing Sources",myID_Poseidon
END IF
IF ( .TRUE. ) THEN

    CALL Poseidon_Input_Sources_AMReX( MF_Source,                  &
                                       nLevels,                     &
                                       nVars_Source,                &
                                       NE, NQ,                      &
                                       R_Quad, T_Quad, P_Quad,      &
                                       LeftLimit, RightLimit        )


ELSE
    CALL Poseidon_Input_Sources(myID_Poseidon,                      &
                                myID_Poseidon,                      &
                                myID_Poseidon,                      &
                                Local_E, Local_S, Local_Si,         &
                                NE(1), NE(2), NE(3),                &
                                NQ(1), NQ(2), NQ(3),                &
                                R_Quad, T_Quad, P_Quad,             &
                                LeftLimit, RightLimit               )
END IF







IF ( .FALSE. ) THEN

    CALL Output_Poseidon_Sources_3D( Local_E, Local_S, Local_Si,    &
                                     NE(1), NE(2), NE(3),           &
                                     NQ(1), NQ(2), NQ(3),           &
                                     R_Quad, T_Quad, P_Quad,        &
                                     LeftLimit, RightLimit          )

END IF




END SUBROUTINE Driver_SetSource






!+201+##########################################################################!
!                                                                               !
!     SetSource_Parallel                                                            !
!                                                                               !
!###############################################################################!
SUBROUTINE SetSource_Parallel(  nDOF, nVars, nLevels,          &
                                NE_Lower, NE_Upper,             &
                                Local_E, Local_S, Local_Si,     &
                                MF_Sources                      )

INTEGER,    INTENT(IN)                              ::  nDOF, nVars, nLevels
INTEGER,    INTENT(IN), DIMENSION(3)                ::  NE_Lower, NE_Upper
REAL(idp),  INTENT(IN), DIMENSION( nDOF,        &
                                NE_Upper(1),    &
                                NE_Upper(2),    &
                                NE_Upper(3)     )   ::  Local_E

REAL(idp),  INTENT(IN), DIMENSION( nDOF,        &
                                NE_Upper(1),    &
                                NE_Upper(2),    &
                                NE_Upper(3)     )   ::  Local_S

REAL(idp),  INTENT(IN), DIMENSION( nDOF,        &
                                NE_Upper(1),    &
                                NE_Upper(2),    &
                                NE_Upper(3),    &
                                1:3             )   ::  Local_Si


TYPE(amrex_multifab), INTENT(INOUT)                 ::  MF_Sources(0:nLevels-1)

REAL(idp), CONTIGUOUS, POINTER                      :: p(:,:,:,:)

INTEGER                                             :: PE, TE, RE, Var
INTEGER, DIMENSION(3)                               :: ELo, EHi

INTEGER                                             :: Here, There, lvl

TYPE(amrex_mfiter)                                  :: mfi
TYPE(amrex_box)                                     :: Box

INTEGER                                             :: ierr, process

PRINT*,"In SetSource_Parallel ", myID_Poseidon




DO lvl = 0,nLevels-1
    CALL amrex_mfiter_build(mfi, MF_Sources(lvl), tiling = .false. )
    DO WHILE(mfi%next())

        CALL MPI_All_Print("In mfiter loop")

        p => MF_Sources(lvl)%dataPtr(mfi)
        Box = mfi%tilebox()

        ELo = Box%lo
        EHi = Box%hi


        DO PE = ELo(3), EHi(3)
        DO TE = ELo(2), EHi(2)
        DO RE = ELo(1), EHi(1)

            Here  = 1
            There = nDOF
            p(RE,TE,PE, Here:There ) = Local_E(:,RE,TE,PE)
            
            Here  = 1*nDOF+1
            There = 2*nDOF
            p(RE,TE,PE, Here:There ) = Local_S(:,RE,TE,PE)

            DO Var = 3,5
                Here  = (Var-1)*nDOF+1
                There = Var*nDOF
                p(RE,TE,PE,Here:There) = Local_Si(:,RE,TE,PE,Var-2)
            END DO

        

        END DO
        END DO
        END DO

    END DO
    CALL amrex_mfiter_destroy(mfi)
END DO ! lvl




!CALL MPI_Barrier(Poseidon_Comm_World, ierr)
!CALL MPI_Master_Print("At the end of SetSource_Parallel")
!CALL STOP_MPI(ierr)


END SUBROUTINE SetSource_Parallel


END MODULE Driver_SetSource_Module
