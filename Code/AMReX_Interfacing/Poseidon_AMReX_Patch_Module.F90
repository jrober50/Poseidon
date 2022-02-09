   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Poseidon_AMReX_Patch_Module                                           !##!
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

USE ISO_C_BINDING

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
  amrex_multifab_build
USE amrex_amr_module, ONLY: &
  amrex_geom, &
  amrex_ref_ratio,  &
  amrex_interp_cell_cons
USE amrex_fillpatch_module, ONLY: &
  amrex_fillpatch, &
  amrex_fillcoarsepatch


USE Variables_Driver_AMReX,  &
            ONLY :  MF_Driver_Source,      &
                    GM_Driver_Source,      &
                    MF_Src_nComps,  &
                    MF_Src_nGhost,  &
                    lo_bc,          &
                    hi_bc,          &
                    t_new,          &
                    t_old,          &
                    UseTiling

#endif

USE Poseidon_Kinds_Module, &
            ONLY :  idp


IMPLICIT NONE


CONTAINS


#ifdef POSEIDON_AMREX_FLAG
!+101+##########################################################################!
!                                                                               !
!                                                     				!
!                                                                               !
!###############################################################################!
SUBROUTINE FillPatch( Level, Time, Src )

INTEGER,                INTENT(in)          ::  Level
REAL(idp),              INTENT(in)          ::  Time
TYPE(amrex_multifab),   INTENT(inout)       ::  Src

INTEGER, PARAMETER                          ::  sComp = 1
INTEGER, PARAMETER                          ::  dComp = 1
INTEGER                                     ::  nComp

!INTEGER :: lo_bc(amrex_spacedim,MF_Src_nComps) = amrex_bc_int_dir
!INTEGER :: hi_bc(amrex_spacedim,MF_Src_nComps) = amrex_bc_int_dir

nComp = MF_Driver_Source(Level)%nComp()


IF ( Level == 0 ) THEN

    CALL amrex_fillpatch( Src,  t_old(Level),  MF_Driver_Source(Level),    &
                                t_new(Level),  MF_Driver_Source(Level),    &
                                GM_Driver_Source(Level), FillPhysicalBC,   &
                                Time,                               &
                                sComp,                              &
                                dComp,                              &
                                nComp                               )

ELSE


    CALL amrex_fillpatch( Src,  t_old(Level-1),  MF_Driver_Source(Level-1),    &
                                t_new(Level-1),  MF_Driver_Source(Level-1),    &
                                GM_Driver_Source(Level-1), FillPhysicalBC,     &
                                t_old(Level),  MF_Driver_Source(Level),        &
                                t_new(Level),  MF_Driver_Source(Level),        &
                                GM_Driver_Source(Level), FillPhysicalBC,       &
                                Time,                                   &
                                sComp,                                  &
                                dComp,                                  &
                                nComp,                                  &
                                amrex_ref_ratio(Level-1),               &
                                amrex_interp_cell_cons,                 &
                                lo_bc,                                  &
                                hi_bc                                   )


END IF



END SUBROUTINE FillPatch





!+102+##########################################################################!
!                                                                               !
!                                                                     !
!                                                                               !
!###############################################################################!
SUBROUTINE FillCoarsePatch( Level, Time, Src )

INTEGER,                INTENT(in)          ::  Level
REAL(idp),              INTENT(in)          ::  Time
TYPE(amrex_multifab),   INTENT(inout)       ::  Src

INTEGER, PARAMETER                          ::  sComp = 1
INTEGER, PARAMETER                          ::  dComp = 1
INTEGER                                     ::  nComp

!INTEGER :: lo_bc(amrex_spacedim,MF_Src_nComps) = amrex_bc_int_dir
!INTEGER :: hi_bc(amrex_spacedim,MF_Src_nComps) = amrex_bc_int_dir

nComp = MF_Driver_Source(Level)%nComp()

CALL amrex_fillcoarsepatch(Src, t_old(Level-1),  MF_Driver_Source(Level-1),    &
                                t_new(Level-1),  MF_Driver_Source(Level-1),    &
                                amrex_geom(Level-1), FillPhysicalBC,    &
                                amrex_geom(Level), FillPhysicalBC,      &
                                Time,                                   &
                                sComp,                                  &
                                dComp,                                  &
                                nComp,                                  &
                                amrex_ref_ratio(Level-1),               &
                                amrex_interp_cell_cons,                 &
                                lo_bc,                                  &
                                hi_bc                                   )

END SUBROUTINE FillCoarsePatch













!+201+##########################################################################!
!                                                                               !
!                                                                     !
!                                                                               !
!###############################################################################!
SUBROUTINE FillPhysicalBC( pMF, sComp, nComp, Time, pGEOM ) BIND(c)

  ! --- No INTENT here because amrex source code doesn't have it ---

TYPE(c_ptr),    VALUE :: pMF, pGEOM
INTEGER(c_int), VALUE :: sComp, nComp
REAL(idp),       VALUE :: Time

TYPE(amrex_geometry) :: GEOM
TYPE(amrex_multifab) :: MF
TYPE(amrex_mfiter)   :: MFI
INTEGER              :: pLo(4), pHi(4)
REAL(idp), CONTIGUOUS, POINTER :: p(:,:,:,:)

IF( .NOT. amrex_is_all_periodic() )THEN

GEOM = pGEOM
MF   = pMF

!$OMP PARALLEL PRIVATE(MFI,p,pLo,pHi)
CALL amrex_mfiter_build( MFI, MF, tiling = UseTiling )

DO WHILE( MFI % next() )

  p => MF % DataPtr( MFI )

  ! Check if part of this box is outside the domain
  IF( .NOT. GEOM % DOMAIN % CONTAINS( p ) )THEN

    pLo = LBOUND( p )
    pHi = UBOUND( p )

    CALL amrex_filcc &
           ( p, pLo, pHi, &
             GEOM % DOMAIN % lo, GEOM % DOMAIN % hi, &
             GEOM % dX, &
             GEOM % get_physical_location( pLo ), &
             lo_bc, hi_bc )

  END IF

END DO
!$OMP END PARALLEL

CALL amrex_mfiter_destroy( MFI )

END IF

END SUBROUTINE FillPhysicalBC


#else




SUBROUTINE FillPatch()
END SUBROUTINE FillPatch


SUBROUTINE FillCoarsePatch()
END SUBROUTINE FillCoarsePatch


SUBROUTINE FillPhysicalBC()
END SUBROUTINE FillPhysicalBC


#endif

END MODULE Poseidon_AMReX_Patch_Module
