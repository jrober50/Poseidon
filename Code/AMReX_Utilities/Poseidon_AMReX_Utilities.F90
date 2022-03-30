  !##################################################################!
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################!
!##!                                                                !##!
!##!                                                                !##!
MODULE Poseidon_AMReX_Utilities_Module                              !##!
!##!                                                                !##!
!##!________________________________________________________________!##!
!##!                                                                !##!
!##!    Contains:                                                   !##!
!##!                                                                !##!
!##!    +101+       Poseidon2AMReX_ES                               !##!
!##!    +102+       Poseidon2AMReX_Si                               !##!
!##!                                                                !##!
!##!    +201+       AMReX2Poseidon_ES                               !##!
!##!    +202+       AMReX2Poseidon_Si                               !##!
!##!    +203+       AMReX2Poseidon_ES_Long                          !##!
!##!    +204+       AMReX2Poseidon_Si_Long                          !##!
!##!                                                                !##!
!######################################################################!
!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  !##################################################################!


USE Poseidon_Kinds_Module, &
               ONLY :  idp

#ifdef POSEIDON_AMREX_FLAG
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
  amrex_multifab_build, &
  amrex_mfiter, &
  amrex_mfiter_build, &
  amrex_mfiter_destroy
USE amrex_parmparse_module, ONLY: &
  amrex_parmparse,       &
  amrex_parmparse_build, &
  amrex_parmparse_destroy
#endif

IMPLICIT NONE





#ifdef POSEIDON_AMREX_FLAG
INTERFACE Poseidon2AMReX
  MODULE PROCEDURE Poseidon2AMReX_ES
  MODULE PROCEDURE Poseidon2AMReX_Si
END INTERFACE Poseidon2AMReX

INTERFACE AMReX2Poseidon
  MODULE PROCEDURE AMReX2Poseidon_ES
  MODULE PROCEDURE AMReX2Poseidon_Si
  MODULE PROCEDURE AMReX2Poseidon_ES_Long
  MODULE PROCEDURE AMReX2Poseidon_Si_Long
END INTERFACE AMReX2Poseidon
#endif






CONTAINS






#ifdef POSEIDON_AMREX_FLAG
!+101+##########################################################################!
!                                                                               !
!                  Poseidon2AMReX_ES                                            !
!                                                                               !
!###############################################################################!
SUBROUTINE Poseidon2AMReX_ES(  VarSize,                 &
                                LowerNE, UpperNE,       &
                                PVar, AVar              )

INTEGER,  INTENT(IN)                    :: VarSize

INTEGER, DIMENSION(3), INTENT(IN)       :: LowerNE
INTEGER, DIMENSION(3), INTENT(IN)       :: UpperNE

REAL(idp), DIMENSION(VarSize,UpperNE(1),UpperNE(2),UpperNE(3)), INTENT(IN)      :: PVar

TYPE(amrex_multifab), INTENT(INOUT)             ::  AVar

REAL(idp), CONTIGUOUS, POINTER                  :: p(:,:,:,:)

INTEGER     :: PE, TE, RE

type(amrex_mfiter) :: mfi


! Transfer Source_Si !

!$OMP parallel private(mfi,PE,TE,RE,Var,p)
CALL amrex_mfiter_build(mfi, AVar, tiling = .true. )
DO WHILE(mfi%next())

    p => AVar%dataPtr(mfi)

    DO PE = LowerNE(3), UpperNE(3)
    DO TE = LowerNE(2), UpperNE(2)
    DO RE = LowerNE(1), UpperNE(1)

    p(RE,TE,PE,1:VarSize)          &
        = PVar(1:VarSize,RE,TE,PE)

    END DO
    END DO
    END DO

END DO
CALL amrex_mfiter_destroy(mfi)
!$OMP end parallel


END SUBROUTINE Poseidon2AMReX_ES








!+102+##########################################################################!
!                                                                               !
!                  Poseidon2AMReX_Si                                            !
!                                                                               !
!###############################################################################!
SUBROUTINE Poseidon2AMReX_Si(  VarSize, NumVars,       &
                                LowerNE, UpperNE,       &
                                PVar, AVar              )

INTEGER,  INTENT(IN)                    :: VarSize
INTEGER,  INTENT(IN)                    :: NumVars

INTEGER, DIMENSION(3), INTENT(IN)       :: LowerNE
INTEGER, DIMENSION(3), INTENT(IN)       :: UpperNE

REAL(idp), DIMENSION(VarSize,UpperNE(1),UpperNE(2),UpperNE(3),NumVars), INTENT(IN)    :: PVar

TYPE(amrex_multifab), INTENT(INOUT)             ::  AVar

REAL(idp), CONTIGUOUS, POINTER                  :: p(:,:,:,:)

INTEGER     :: PE, TE, RE, Var

type(amrex_mfiter) :: mfi


! Transfer Source_Si !

!$OMP parallel private(mfi,PE,TE,RE,Var,p)
CALL amrex_mfiter_build(mfi, AVar, tiling = .true. )
DO WHILE(mfi%next())

    p => AVar%dataPtr(mfi)

    DO PE = LowerNE(3), UpperNE(3)
    DO TE = LowerNE(2), UpperNE(2)
    DO RE = LowerNE(1), UpperNE(1)

    DO Var = 1,NumVars
        p(RE,TE,PE,(Var-1)*VarSize+1:Var*VarSize)          &
            = PVar(1:VarSize,RE,TE,PE,Var)
    END DO

    END DO
    END DO
    END DO

END DO
CALL amrex_mfiter_destroy(mfi)
!$OMP end parallel


END SUBROUTINE Poseidon2AMReX_Si










!+201+###########################################################################!
!                                                                                !
!                  AMReX2Poseidon_ES                                             !
!                                                                                !
!################################################################################!
SUBROUTINE AMReX2Poseidon_ES( VarSize,                &
                              LowerNE, UpperNE,       &
                              PVar, AVar              )

INTEGER,  INTENT(IN)                    :: VarSize

INTEGER, DIMENSION(3), INTENT(IN)       :: LowerNE
INTEGER, DIMENSION(3), INTENT(IN)       :: UpperNE

REAL(idp), DIMENSION(VarSize,UpperNE(1),UpperNE(2),UpperNE(3)), INTENT(OUT)      :: PVar

TYPE(amrex_multifab), INTENT(IN)                ::  AVar

REAL(idp), CONTIGUOUS, POINTER                  :: p(:,:,:,:)

INTEGER                                         :: RE, TE, PE

type(amrex_mfiter) :: mfi


!$OMP parallel private(mfi,PE,TE,RE,Var,p)
CALL amrex_mfiter_build(mfi, AVar, tiling = .false. )
DO WHILE(mfi%next())

    p => AVar%dataPtr(mfi)

    DO PE = LowerNE(3), UpperNE(3)
    DO TE = LowerNE(2), UpperNE(2)
    DO RE = LowerNE(1), UpperNE(1)
        PVar(1:VarSize,RE,TE,PE)                &
            = p(RE,TE,PE,1:VarSize)
    END DO
    END DO
    END DO

END DO
CALL amrex_mfiter_destroy(mfi)
!$OMP end parallel


END SUBROUTINE AMReX2Poseidon_ES



!+202+###########################################################################!
!                                                                                !
!                  AMReX2Poseidon_Si                                             !
!                                                                                !
!################################################################################!
SUBROUTINE AMReX2Poseidon_Si( VarSize, NumVars,       &
                              LowerNE, UpperNE,       &
                              PVar, AVar              )

INTEGER,  INTENT(IN)                    :: VarSize
INTEGER,  INTENT(IN)                    :: NumVars

INTEGER, DIMENSION(3), INTENT(IN)       :: LowerNE
INTEGER, DIMENSION(3), INTENT(IN)       :: UpperNE

REAL(idp), DIMENSION(VarSize,UpperNE(1),UpperNE(2),UpperNE(3),NumVars), INTENT(OUT)      :: PVar

TYPE(amrex_multifab), INTENT(IN)                ::  AVar

REAL(idp), CONTIGUOUS, POINTER                  :: p(:,:,:,:)

INTEGER                                         :: RE, TE, PE, var

type(amrex_mfiter) :: mfi


!$OMP parallel private(mfi,PE,TE,RE,Var,p)
CALL amrex_mfiter_build(mfi, AVar, tiling = .false. )
DO WHILE(mfi%next())

    p => AVar%dataPtr(mfi)

    DO PE = LowerNE(3), UpperNE(3)
    DO TE = LowerNE(2), UpperNE(2)
    DO RE = LowerNE(1), UpperNE(1)
    DO var = 1,NumVars
        PVar(1:VarSize,RE,TE,PE,var)                &
            = p(RE,TE,PE,(var-1)*VarSize+1:var*VarSize)
    END DO
    END DO
    END DO
    END DO

END DO
CALL amrex_mfiter_destroy(mfi)
!$OMP end parallel


END SUBROUTINE AMReX2Poseidon_Si



!+204+###########################################################################!
!                                                                                !
!                  AMReX2Poseidon_ES_Long                                       !
!                                                                                !
!################################################################################!
SUBROUTINE AMReX2Poseidon_ES_Long( VarSize,                &
                                    LowerNE, UpperNE,       &
                                    PVar, AVar              )

INTEGER, DIMENSION(3), INTENT(IN)       :: VarSize

INTEGER, DIMENSION(3), INTENT(IN)       :: LowerNE
INTEGER, DIMENSION(3), INTENT(IN)       :: UpperNE

REAL(idp), DIMENSION(VarSize(1),VarSize(2),VarSize(3),          &
                     UpperNE(1),UpperNE(2),UpperNE(3)), INTENT(OUT)      :: PVar

TYPE(amrex_multifab), INTENT(IN)                ::  AVar

REAL(idp), CONTIGUOUS, POINTER                  :: p(:,:,:,:)

INTEGER                                         :: RE, TE, PE
INTEGER                                         :: RD, TD, PD, Here

type(amrex_mfiter) :: mfi


!$OMP parallel private(mfi,PE,TE,RE,Var,p)
CALL amrex_mfiter_build(mfi, AVar, tiling = .false. )
DO WHILE(mfi%next())

    p => AVar%dataPtr(mfi)

    DO PE = LowerNE(3), UpperNE(3)
    DO TE = LowerNE(2), UpperNE(2)
    DO RE = LowerNE(1), UpperNE(1)
    DO PD = 1,VarSize(3)
    DO TD = 1,VarSize(2)
    DO RD = 1,VarSize(1)

        Here = RD + (TD-1)*VarSize(1) + (PD-1)*VarSize(1)*VarSize(2)

        PVar(RD,TD,PD,RE,TE,PE) = p(RE,TE,PE,Here)

    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

END DO
CALL amrex_mfiter_destroy(mfi)
!$OMP end parallel


END SUBROUTINE AMReX2Poseidon_ES_Long





!+204+###########################################################################!
!                                                                                !
!                  AMReX2Poseidon_Si_Long                                       !
!                                                                                !
!################################################################################!
SUBROUTINE AMReX2Poseidon_Si_Long( VarSize, NumVars,       &
                                    LowerNE, UpperNE,       &
                                    PVar, AVar              )

INTEGER, DIMENSION(3), INTENT(IN)       :: VarSize
INTEGER,  INTENT(IN)                    :: NumVars

INTEGER, DIMENSION(3), INTENT(IN)       :: LowerNE
INTEGER, DIMENSION(3), INTENT(IN)       :: UpperNE

REAL(idp), DIMENSION(VarSize(1),VarSize(2),VarSize(3),          &
                    UpperNE(1),UpperNE(2),UpperNE(3),NumVars), INTENT(OUT)      :: PVar

TYPE(amrex_multifab), INTENT(IN)                ::  AVar

REAL(idp), CONTIGUOUS, POINTER                  :: p(:,:,:,:)

INTEGER                                         :: RE, TE, PE, Var
INTEGER                                         :: RD, TD, PD, Here
INTEGER                                         :: Num_DOF

type(amrex_mfiter) :: mfi


Num_DOF = VarSize(1)*VarSize(2)*VarSize(3)

!$OMP parallel private(mfi,PE,TE,RE,Var,p)
CALL amrex_mfiter_build(mfi, AVar, tiling = .false. )
DO WHILE(mfi%next())

    p => AVar%dataPtr(mfi)

    DO PE = LowerNE(3), UpperNE(3)
    DO TE = LowerNE(2), UpperNE(2)
    DO RE = LowerNE(1), UpperNE(1)
    DO PD = 1,VarSize(3)
    DO TD = 1,VarSize(2)
    DO RD = 1,VarSize(1)
    DO Var = 1,NumVars

        Here = RD + (TD-1)*VarSize(1) + (PD-1)*VarSize(1)*VarSize(2)
        PVar(RD,TD,PD,RE,TE,PE,var)                &
            = p(RE,TE,PE,(Var-1)*Num_DOF+Here)
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO
    END DO

END DO
CALL amrex_mfiter_destroy(mfi)
!$OMP end parallel


END SUBROUTINE AMReX2Poseidon_Si_Long




!+204+###########################################################################!
!                                                                                !
!                  PackageSources_AMReX                                             !
!                                                                                !
!################################################################################!
SUBROUTINE PackSources_AMReX(   nDOF, nVars, nLevels,          &
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


DO lvl = 0,nLevels-1
    CALL amrex_mfiter_build(mfi, MF_Sources(lvl), tiling = .true. )
    DO WHILE(mfi%next())

        p => MF_Sources(lvl)%dataPtr(mfi)
        Box = mfi%tilebox()

        ELo = Box%lo
        EHi = Box%hi

        DO PE = ELo(3), EHi(3)
        DO TE = ELo(3), EHi(3)
        DO RE = ELo(3), EHi(3)

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

END SUBROUTINE PackSources_AMReX







!+204+##########################################################################!
!                                                                               !
!                  UnpackSources_AMReX                                          !
!                                                                               !
!###############################################################################!
SUBROUTINE UnpackSources_AMReX( nDOF, nLevels,                  &
                                NE_Lower, NE_Upper,             &
                                Local_E, Local_S, Local_Si,     &
                                MF_Sources                      )

INTEGER, INTENT(IN)                                     ::  nDOF, nLevels
INTEGER, INTENT(IN),  DIMENSION(3)                      ::  NE_Lower, NE_Upper
REAL(idp), INTENT(OUT), DIMENSION(  nDOF,           &
                                    NE_Upper(1),    &
                                    NE_Upper(2),    &
                                    NE_Upper(3)     )   ::  Local_E

REAL(idp), INTENT(OUT), DIMENSION(  nDOF,           &
                                    NE_Upper(1),    &
                                    NE_Upper(2),    &
                                    NE_Upper(3)     )   ::  Local_S

REAL(idp), INTENT(OUT), DIMENSION(  nDOF,           &
                                    NE_Upper(1),    &
                                    NE_Upper(2),    &
                                    NE_Upper(3),    &
                                    1:3             )   ::  Local_Si


TYPE(amrex_multifab), INTENT(IN)                   ::   MF_Sources(0:nLevels-1)

REAL(idp), CONTIGUOUS, POINTER                      :: p(:,:,:,:)

INTEGER                                             :: PE, TE, RE, Var
INTEGER, DIMENSION(3)                               :: ELo, EHi
 
INTEGER                                             :: Here, There, lvl

TYPE(amrex_mfiter)                                  :: mfi
TYPE(amrex_box)                                     :: Box


DO lvl = 0,nLevels-1
    CALL amrex_mfiter_build(mfi, MF_Sources(lvl), tiling = .true. )
    DO WHILE(mfi%next())

        p => MF_Sources(lvl)%dataPtr(mfi)
        Box = mfi%tilebox()

        ELo = Box%lo
        EHi = Box%hi

        DO PE = ELo(3), EHi(3)
        DO TE = ELo(3), EHi(3)
        DO RE = ELo(3), EHi(3)

            Here  = 0
            There = nDOF-1
            Local_E(:,RE,TE,PE) = p(RE,TE,PE, Here:There )

            Here  = 1*nDOF
            There = 2*nDOF-1
            Local_S(:,RE,TE,PE) = p(RE,TE,PE, Here:There )

            DO Var = 3,5
                Here  = (Var-1)*nDOF
                There = Var*nDOF - 1
                 Local_Si(:,RE,TE,PE,Var-2) = p(RE,TE,PE,Here:There)
            END DO

        END DO
        END DO
        END DO

    END DO
    CALL amrex_mfiter_destroy(mfi)
END DO ! lvl

END SUBROUTINE UnpackSources_AMReX





!!+204+###########################################################################!
!!                                                                                !
!!                  Poseidon2Thornado                                             !
!!                                                                                !
!!################################################################################!
!SUBROUTINE Poseidon2Thornado()
!
!TYPE(amrex_multifab), ALLOCATABLE           :: PVars(:)
!TYPE(amrex_multifab), ALLOCATABLE           :: tVars(:)
!
!TYPE(amrex_box)                                         ::  Box
!TYPE(amrex_boxarray)                                    ::  BoxArray
!TYPE(amrex_distromap)                                   ::  DistroMap
!
!
!INTEGER                                     :: nDOF    = 12
!INTEGER                                     :: nLevels = 1
!INTEGER                                     :: nPVars  = 5
!INTEGER                                     :: ntVars  = 13
!INTEGER                                     :: nGhost  = 0
!
!INTEGER                                     :: nGhost_Copy = 0
!
!INTEGER                                     :: Copy_Length
!INTEGER                                     :: SRC_Start
!INTEGER                                     :: DST_Start
!
!INTEGER, DIMENSION(3)                       :: NE_Lower = (/1,1,1/)
!INTEGER, DIMENSION(3)                       :: NE_Upper = (/10,10,10/)
!
!ALLOCATE( PVars(0:nLevels-1) )
!ALLOCATE( tVars(0:nLevels-1) )
!
!
!
!Box = amrex_box(NE_Lower, NE)
!CALL amrex_boxarray_build(BoxArray, Box)
!CALL amrex_distromap_build(DistroMap,BoxArray)
!
!CALL amrex_multifab_build(PVars,BoxArray,DistroMap,nPVars*nDOF, nghost)
!CALL amrex_multifab_build(tVars,BoxArray,DistroMap,ntVars*nDOF, nghost)
!
!! PVars(1 = Conformal Factor, 2 = Lapse Function, 3,4,5 = Shift Vector )
!! tVars(13 = Conformal Factor, 9 = Lapse Function, 10,11,12 = Shift Vector)
!
!Copy_Length = nDOF
!SRC_Start   = 1 !(1-1)*nDOF+1
!DST_Start   = (13-1)*nDOF+1  ! Depends on index of first component
!
!tVars%COPY(PVars, SRC_Start, DST_Start, Copy_Length, nGhost_Copy)   ! Copy Conformal Factor
!
!
!
!Copy_Length = nDOF
!SRC_Start   = (2-1)*nDOF+1
!DST_Start   = (9-1)*nDOF+1  ! Depends on index of first component
!
!tVars%COPY(PVars, SRC_Start, DST_Start, Copy_Length, nGhost_Copy)   ! Copy Lapse Function
!
!
!Copy_Length = 3*nDOF
!SRC_Start   = (3-1)*nDOF+1
!DST_Start   = (10-1)*nDOF+1  ! Depends on index of first component
!
!tVars%COPY(PVars, SRC_Start, DST_Start, Copy_Length, nGhost_Copy)   ! Copy Shift Vector Components together
!
!END SUBROUTINE Poseidon2Thornado



#endif










END MODULE Poseidon_AMReX_Utilities_Module
