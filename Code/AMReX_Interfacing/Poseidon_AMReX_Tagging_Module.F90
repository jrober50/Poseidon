   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Poseidon_AMReX_Tagging_Module                                         !##!
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
#endif



USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Variables_Mesh, &
            ONLY :  NUM_R_ELEMENTS,             &
                    NUM_T_ELEMENTS,             &
                    NUM_P_ELEMENTS

IMPLICIT NONE


CONTAINS


 !+101+################################################################!
!                                                                       !
!          Tag_Source_Elements                                          !
!                                                                       !
 !#####################################################################!
SUBROUTINE Tag_Source_Elements( Level,              &
                                BLo, BHi,           & ! Current Box Bounds
                                SLo, SHi,           & ! Src Bounds
                                TLo, THi,           & ! Tag Bounds
                                nComps,             &
                                Src,                & ! Pointer to Source Data
                                SetTag, ClearTag,   &
                                Tag                 ) ! Pointer to Tag Data

INTEGER,                INTENT(IN)      ::  Level
INTEGER,                INTENT(IN)      ::  BLo(3), BHi(3)
INTEGER,                INTENT(IN)      ::  SLo(3), SHi(3)
INTEGER,                INTENT(IN)      ::  TLo(4), THi(4)
INTEGER,                INTENT(IN)      ::  nComps
REAL(idp),              INTENT(IN)      ::  Src(SLo(1):SHi(1),  &
                                                SLo(2):SHi(2),  &
                                                SLo(3):SHi(3),  &
                                                nComps          )

CHARACTER(KIND=c_char), INTENT(IN)      ::  SetTag
CHARACTER(KIND=c_char), INTENT(IN)      ::  ClearTag

CHARACTER(KIND=c_char), INTENT(INOUT)   ::  Tag(TLo(1):THi(1),  &
                                                TLo(2):THi(2),  &
                                                TLo(3):THi(3),  &
                                                TLo(4):THi(4)   )

INTEGER                                 ::  Third
INTEGER                                 ::  Fourth
INTEGER                                 ::  Half
INTEGER                                 ::  UFourth
INTEGER                                 ::  re, te, pe


Third   = Num_R_Elements
Half    = Num_R_Elements/2**(Level+1)
Fourth  = Half/2
UFourth = Num_R_Elements - Fourth


!PRINT*,Num_R_Elements, 2**(Level+1),Half

DO pe = BLo(3),BHi(3)
DO te = BLo(2),BHi(2)
DO re = BLo(1),BHi(1)
    

    
    IF ( re < Half ) THEN
!        PRINT*,"Refining",re
        Tag(re,te,pe,1) = SetTag
    ELSE
        Tag(re,te,pe,1) = ClearTag
    END IF

!    IF ( re < Fourth ) THEN
!        Tag(re,te,pe,1) = SetTag
!    ELSE IF ( re > UFourth ) THEN
!        Tag(re,te,pe,1) = SetTag
!    ELSE
!        Tag(re,te,pe,1) = ClearTag
!    END IF



END DO ! re
END DO ! te
END DO ! pe




END SUBROUTINE Tag_Source_Elements



END MODULE Poseidon_AMReX_Tagging_Module
