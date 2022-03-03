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
            ONLY :  DOMAIN_DIM,         &
                    Poseidon_Remesh_Flag

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




#ifdef POSEIDON_AMREX_FLAG
USE Variables_AMReX_Core, &
            ONLY :  MF_Source,          &
                    AMReX_Num_Levels

USE Poseidon_AMReX_Utilities_Module,    &
            ONLY :  AMReX2Poseidon,     &
                    UnpackSources_AMReX


USE Initialization_XCFC_with_AMReX_Module, &
            ONLY :  Initialization_XCFC_with_AMReX


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
SUBROUTINE Poseidon_Input_Sources_AMREX( MF_SRC_Input, MF_Src_nComps, Num_Levels )

TYPE(amrex_multifab),   INTENT(IN)                  ::  MF_SRC_Input(0:Num_Levels-1)
INTEGER,                INTENT(IN)                  ::  MF_Src_nComps
INTEGER,                INTENT(IN)                  ::  Num_Levels

INTEGER                                             ::  RE, TE, PE
INTEGER, DIMENSION(3)                               ::  ELo, EHi
INTEGER                                             ::  level
INTEGER                                             ::  nghost = 0
INTEGER                                             ::  Num_DOF

INTEGER                                             :: nComp
REAL(idp), CONTIGUOUS, POINTER                      :: SRC(:,:,:,:)
REAL(idp), CONTIGUOUS, POINTER                      :: DST(:,:,:,:)

TYPE(amrex_mfiter)                                  :: mfi
TYPE(amrex_box)                                     :: Box



DO level = 0,AMReX_Num_Levels-1

    CALL amrex_multifab_build(  MF_Source(level),           &
                                MF_Src_Input(Level)%BA,     &
                                MF_Src_Input(Level)%DM,     &
                                MF_Src_nComps, 1                        )

    PRINT*,"To Do : Interpolate between Source quadrature points."
    MF_Source(Level)=MF_Src_Input(Level)
END DO




IF ( Poseidon_Remesh_Flag ) THEN

    Call Initialization_XCFC_with_AMReX()

END IF





END SUBROUTINE Poseidon_Input_Sources_AMREX














!!+101+##########################################################################!
!!                                                                               !
!!                           Poseidon_Input_Sources                              !
!!                                                                               !
!!###############################################################################!
!SUBROUTINE Poseidon_Input_XCFC_Sources_AMREX(   MF_Sources, nLevels_Input,   &
!                                                NE_Lower, NE_Upper,         &
!                                                NQ, Si_Dim                  )
!
!
!
!TYPE(amrex_multifab), INTENT(IN)        ::  MF_Sources(0:nLevels_Input-1)
!
!
!END SUBROUTINE Poseidon_Input_XCFC_Sources_AMREX





#else


SUBROUTINE Poseidon_Input_Sources_AMREX( )
END SUBROUTINE Poseidon_Input_Sources_AMREX

SUBROUTINE Poseidon_XCFC_Input_Sources_AMREX( )
END SUBROUTINE Poseidon_XCFC_Input_Sources_AMREX


#endif




END MODULE Source_Input_AMReX
