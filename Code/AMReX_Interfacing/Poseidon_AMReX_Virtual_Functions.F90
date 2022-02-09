   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Poseidon_AMReX_Virtual_Functions_Module                               !##!
!##!                                                                         !##!
!##!_________________________________________________________________________!##!
!##!                                                                         !##!
!##!                                                                         !##!
!##!=========================================================================!##!
!##!                                                                         !##!
!##!    Contains:                                                            !##!
!##!                                                                         !##!
!##!    +101+   VF_Make_New_Level_From_Scratch                               !##!
!##!    +102+   VF_Make_New_Level_From_Coarse                                !##!
!##!    +103+   VF_Remake_Level                                              !##!
!##!    +104+   VF_Clear_Level                                               !##!
!##!    +105+   VF_Error_Estimate                                            !##!
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
USE amrex_tagbox_module, ONLY: &
  amrex_tagboxarray

USE Variables_Driver_AMReX,  &
ONLY :  MF_Driver_Source,      &
        MF_Src_nComps,  &
        MF_Src_nGhost
#endif


USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Variables_AMReX_Core, &
            ONLY :  AMReX_Tiling

USE Poseidon_AMReX_Tagging_Module, &
            ONLY : Tag_Source_Elements

USE Poseidon_AMReX_Patch_Module, &
            ONLY :  FillCoarsePatch,    &
                    FillPatch

USE Poseidon_Init_Data_On_Level_Module, &
            ONLY :  Poseidon_Init_Data_On_Level

IMPLICIT NONE



CONTAINS


#ifdef POSEIDON_AMREX_FLAG
 !+101+################################################################!
!                                                                       !
!          VF_Make_New_Level_From_Scratch                               !
!                                                                       !
 !#####################################################################!
SUBROUTINE VF_Make_New_Level_From_Scratch(lev, time, pba, pdm) bind(c)

INTEGER,            INTENT(IN), VALUE           ::  lev
REAL(amrex_real),   INTENT(IN), VALUE           ::  time
TYPE(c_ptr),        INTENT(IN), VALUE           ::  pba, pdm

TYPE(amrex_boxarray)                            :: BA
TYPE(amrex_distromap)                           :: DM
TYPE(amrex_mfiter)                              :: mfi
TYPE(amrex_box)                                 :: BX
REAL(amrex_real),   CONTIGUOUS, POINTER         :: p(:,:,:,:)


BA = pba
DM = pdm

!t_new(lev) = time
!t_old(lev) = time - 1.e200_amrex_real
CALL VF_Clear_Level(lev)

!PRINT*,MF_Src_nComps, MF_Src_nGhost

CALL amrex_multifab_build( MF_Driver_Source(lev), BA, DM, MF_Src_nComps, MF_Src_nGhost )

! Fill New Level with Source Data !
CALL amrex_mfiter_build(mfi, MF_Driver_Source(lev))


DO WHILE (mfi%next())
    BX = mfi%tilebox()
    p  => MF_Driver_Source(lev)%dataptr(mfi)

    CALL Poseidon_Init_Data_On_Level(lev,           &
                                     BX%lo, BX%hi,  &
                                     BX%lo, BX%hi,  &
                                     MF_Src_nComps, &
                                     p              )
END DO

CALL amrex_mfiter_destroy(mfi)



END SUBROUTINE VF_Make_New_Level_From_Scratch






 !+102+################################################################!
!                                                                       !
!          VF_Make_New_Level_From_Coarse                                !
!                                                                       !
 !#####################################################################!
SUBROUTINE VF_Make_New_Level_From_Coarse(lev, time, pba, pdm) bind(c)

INTEGER,            INTENT(IN), VALUE           ::  lev
REAL(amrex_real),   INTENT(IN), VALUE           ::  time
TYPE(c_ptr),        INTENT(IN), VALUE           ::  pba, pdm

TYPE(amrex_boxarray)                            :: BA
TYPE(amrex_distromap)                           :: DM
!TYPE(amrex_mfiter)                              :: mfi
!TYPE(amrex_box)                                 :: BX
!REAL(amrex_real),   CONTIGUOUS, POINTER         :: p(:,:,:,:)

BA = pba
DM = pdm

!t_new(lev) = time
!t_old(lev) = time - 1.e200_amrex_real

CALL VF_Clear_Level(lev)

CALL amrex_multifab_build( MF_Driver_Source(lev), BA, DM, MF_Src_nComps, MF_Src_nGhost )

CALL FillCoarsePatch(lev, Time, MF_Driver_Source(lev))


END SUBROUTINE VF_Make_New_Level_From_Coarse






 !+103+################################################################!
!                                                                       !
!          VF_Remake_Level                                              !
!                                                                       !
 !#####################################################################!
SUBROUTINE VF_Remake_Level(lev, time, pba, pdm) bind(c)

INTEGER,            INTENT(IN), VALUE           ::  lev
REAL(amrex_real),   INTENT(IN), VALUE           ::  time
TYPE(c_ptr),        INTENT(IN), VALUE           ::  pba, pdm

TYPE(amrex_boxarray)                            :: BA
TYPE(amrex_distromap)                           :: DM
TYPE(amrex_multifab)                            :: New_MF_Driver_Source

BA = pba
DM = pdm

CALL amrex_multifab_build(New_MF_Driver_Source, BA, DM, MF_Src_nComps, 0)

CALL FillPatch(lev, time, New_MF_Driver_Source)

CALL VF_Clear_Level(lev)

!t_new(lev) = time
!t_old(lev) = time - 1.e200_amrex_real

CALL amrex_multifab_build(MF_Driver_Source(lev), BA, DM, MF_Src_nComps, MF_Src_nGhost )

CALL MF_Driver_Source(lev)%copy(New_MF_Driver_Source, 1, 1, MF_Src_nComps, 0)

CALL amrex_multifab_destroy(New_MF_Driver_Source)


END  SUBROUTINE VF_Remake_Level








 !+104+################################################################!
!                                                                       !
!          VF_Clear_Level                                               !
!                                                                       !
 !#####################################################################!
SUBROUTINE VF_Clear_Level(lev) bind(c)

INTEGER,            INTENT(IN), VALUE           ::  lev

CALL amrex_multifab_destroy(MF_Driver_Source(lev))

END SUBROUTINE VF_Clear_Level








 !+105+################################################################!
!                                                                       !
!          VF_Error_Estimate                                            !
!                                                                       !
 !#####################################################################!
SUBROUTINE VF_Error_Estimate(lev, cp, t, SetTag, ClearTag ) bind(c)

INTEGER,                INTENT(IN), VALUE       ::  lev
TYPE(c_ptr),            INTENT(IN), VALUE       ::  cp
REAL(amrex_real),       INTENT(IN), VALUE       ::  t
CHARACTER(KIND=c_char), INTENT(IN), VALUE       ::  SetTag, ClearTag

TYPE(amrex_tagboxarray)                         ::  Tag
TYPE(amrex_mfiter)                              ::  mfi
TYPE(amrex_box)                                 ::  BX
REAL(idp),              CONTIGUOUS, POINTER     ::  Src(:,:,:,:)
CHARACTER(KIND=c_char), CONTIGUOUS, POINTER     ::  TagArr(:,:,:,:)


Tag = cp

CALL amrex_mfiter_build( mfi, MF_Driver_Source(lev), Tiling = AMReX_Tiling)
DO WHILE( mfi%next() )

    BX = mfi%tilebox()

    Src     => MF_Driver_Source(lev)%dataptr(mfi)
    TagArr  => Tag%dataptr(mfi)
    CALL Tag_Source_Elements( lev,                              &
                              BX%lo, BX%hi,                     &
                              LBOUND(Src), UBOUND(Src),         &
                              LBOUND(TagArr), UBOUND(TagArr),   &
                              MF_Src_nComps,                    &
                              Src,                              &
                              SetTag, ClearTag,                 &
                              TagArr                            )

END DO

CALL amrex_mfiter_destroy(mfi)




END SUBROUTINE VF_Error_Estimate


#else

SUBROUTINE VF_Make_New_Level_From_Scratch
END SUBROUTINE VF_Make_New_Level_From_Scratch


SUBROUTINE VF_Make_New_Level_From_Coarse
END SUBROUTINE VF_Make_New_Level_From_Coarse

SUBROUTINE VF_Remake_Level
END SUBROUTINE VF_Remake_Level

SUBROUTINE VF_Clear_Level
END SUBROUTINE VF_Clear_Level

SUBROUTINE VF_Error_Estimate
END SUBROUTINE VF_Error_Estimate

#endif











END MODULE Poseidon_AMReX_Virtual_Functions_Module
