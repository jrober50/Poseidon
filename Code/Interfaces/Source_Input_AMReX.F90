  !################################################################################!
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Source_Input_AMReX                                                           !##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  !################################################################################!


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
#endif

USE Poseidon_Internal_Communication_Module, &
            ONLY :  Poseidon_CFA_Block_Share

USE Poseidon_Parameters, &
            ONLY :  DOMAIN_DIM

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

USE Poseidon_IO_Module, &
            ONLY :  OUTPUT_POSEIDON_SOURCES_1D



USE Variables_AMReX_Core,   &
            ONLY :  AMReX_Levels

#ifdef POSEIDON_AMREX_FLAG
USE Variables_AMReX_Multifabs,  &
            ONLY :  MF_Source,          &
                    BA_Source,          &
                    DM_Source

USE Poseidon_AMReX_Utilities_Module,    &
            ONLY :  AMReX2Poseidon,     &
                    UnpackSources_AMReX
#endif

use mpi



IMPLICIT NONE

CONTAINS


#ifdef POSEIDON_AMREX_FLAG
!+101+##########################################################################!
!                                                                               !
!                           Poseidon_Input_Sources                              !
!                                                                               !
!###############################################################################!
SUBROUTINE Poseidon_Input_Sources_AMREX( MF_SRC_Input,              &
                                         nLevels_Input,             &
                                         nVars_Input,               &
                                         NE, NQ,                    &
                                         R_Quad, T_Quad, P_Quad,    &
                                         LeftLimit, RightLimit      )

TYPE(amrex_multifab), INTENT(IN)                    ::  MF_SRC_Input(0:nLevels_Input-1)
INTEGER, INTENT(IN)                                 ::  nLevels_Input
INTEGER, INTENT(IN)                                 ::  nVars_Input

INTEGER, INTENT(IN), DIMENSION(3)                   ::  NE
INTEGER, INTENT(IN), DIMENSION(3)                   ::  NQ

REAL(idp), INTENT(IN), DIMENSION(NQ(1))             ::  R_Quad
REAL(idp), INTENT(IN), DIMENSION(NQ(2))             ::  T_Quad
REAL(idp), INTENT(IN), DIMENSION(NQ(3))             ::  P_Quad

REAL(idp), INTENT(IN)                               ::  LeftLimit
REAL(idp), INTENT(IN)                               ::  RightLimit

INTEGER                                             ::  RE, TE, PE
INTEGER, DIMENSION(3)                               ::  ELo, EHi
INTEGER                                             ::  lvl
INTEGER                                             ::  nghost = 0
INTEGER                                             ::  Num_DOF

INTEGER                                             :: nComp
REAL(idp), CONTIGUOUS, POINTER                      :: SRC(:,:,:,:)
REAL(idp), CONTIGUOUS, POINTER                      :: DST(:,:,:,:)

TYPE(amrex_mfiter)                                  :: mfi
TYPE(amrex_box)                                     :: Box

Num_DOF = NQ(1)*NQ(2)*NQ(3)


DO lvl = 0,nLevels_Input-1

    BA_Source(lvl) = MF_SRC_Input(lvl)%ba
    DM_Source(lvl) = MF_SRC_Input(lvl)%dm
    nComp = MF_SRC_Input(lvl)%ncomp()
    PRINT*,"nComp ",nComp," Num_DOF*nVars_Input ",Num_DOF*nVars_Input

    CALL amrex_multifab_build(  MF_Source(lvl),         &
                                BA_Source(lvl),         &
                                DM_Source(lvl),         &
                                Num_DOF*nVars_Input,   &
                                nghost                  )
END DO ! lvl Loop


DO lvl = 0,nLevels_Input-1
    CALL amrex_mfiter_build(mfi, MF_Source(lvl), tiling = .false. )
    DO WHILE(mfi%next())

        SRC => MF_SRC_Input(lvl)%dataPtr(mfi)
        DST => MF_Source(lvl)%dataPTR(mfi)
        Box = mfi%tilebox()

        ELo = Box%lo
        EHi = Box%hi

!        PRINT*,"ELo ",Elo
!        PRINT*,"EHi ",EHi
        DO PE = ELo(3), EHi(3)
        DO TE = ELo(2), EHi(2)
        DO RE = ELo(1), EHi(1)


            DST(RE,TE,PE,:) = SRC(RE,TE,PE,:)


        END DO
        END DO
        END DO

    END DO
    CALL amrex_mfiter_destroy(mfi)
END DO ! lvl




END SUBROUTINE Poseidon_Input_Sources_AMREX














!+101+##########################################################################!
!                                                                               !
!                           Poseidon_Input_Sources                              !
!                                                                               !
!###############################################################################!
SUBROUTINE Poseidon_XCFC_Input_Sources_AMREX(   MF_Sources, AMReX_Levels,   &
                                                NE_Lower, NE_Upper,         &
                                                NQ, Si_Dim                  )



TYPE(amrex_multifab), INTENT(IN)                                ::  MF_Sources(0:AMReX_Levels-1)
INTEGER             , INTENT(IN)                                ::  AMReX_Levels


INTEGER, DIMENSION(3), INTENT(IN)                               ::  NE_Lower
INTEGER, DIMENSION(3), INTENT(IN)                               ::  NE_Upper

INTEGER, DIMENSION(3), INTENT(IN)                               ::  NQ
INTEGER              , INTENT(IN)                               ::  Si_Dim

INTEGER                                                         ::  Num_DOF

Num_DOF = NQ(1)*NQ(2)*NQ(3)


CALL UnpackSources_AMReX( Num_DOF, AMReX_Levels,                &
                          NE_Lower, NE_Upper,                   &
                          Block_Source_E, Block_Source_S, Block_Source_Si,  &
                          MF_Sources                            )


END SUBROUTINE Poseidon_XCFC_Input_Sources_AMREX



#else


SUBROUTINE Poseidon_Input_Sources_AMREX( )
END SUBROUTINE Poseidon_Input_Sources_AMREX

SUBROUTINE Poseidon_XCFC_Input_Sources_AMREX( )
END SUBROUTINE Poseidon_XCFC_Input_Sources_AMREX


#endif




END MODULE Source_Input_AMReX
