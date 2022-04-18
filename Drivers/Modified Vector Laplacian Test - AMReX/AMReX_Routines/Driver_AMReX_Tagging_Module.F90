   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Driver_AMReX_Tagging_Module                                         !##!
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

USE Variables_External, &
            ONLY :  Central_E


USE Variables_AMReX_Core, &
            ONLY :  AMReX_Max_Level


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

INTEGER                                 ::  Tag_Style = 1


IF ( Tag_Style == 1 ) THEN

    CALL Tag_By_Half( Level,              &
                      BLo, BHi,           & ! Current Box Bounds
                      SLo, SHi,           & ! Src Bounds
                      TLo, THi,           & ! Tag Bounds
                      nComps,             &
                      Src,                & ! Pointer to Source Data
                      SetTag, ClearTag,   &
                      Tag                 )



ELSE IF ( Tag_Style == 2 ) THEN

    CALL Tag_By_Source( Level,              &
                        BLo, BHi,           & ! Current Box Bounds
                        SLo, SHi,           & ! Src Bounds
                        TLo, THi,           & ! Tag Bounds
                        nComps,             &
                        Src,                & ! Pointer to Source Data
                        SetTag, ClearTag,   &
                        Tag                 )

END IF





END SUBROUTINE Tag_Source_Elements











 !+101+################################################################!
!                                                                       !
!          Tag_Source_Elements                                          !
!                                                                       !
 !#####################################################################!
SUBROUTINE Tag_By_Half( Level,              &
                        BLo, BHi,           & ! Current Box Bounds
                        SLo, SHi,           & ! Src Bounds
                        TLo, THi,           & ! Tag Bounds
                        nComps,             &
                        Src,                & ! Pointer to Source Data
                        SetTag, ClearTag,   &
                        Tag                 )


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
INTEGER                                 ::  TwoThirds
INTEGER                                 ::  Fourth
INTEGER                                 ::  Half
INTEGER                                 ::  UFourth
INTEGER                                 ::  re, te, pe


Third     = Num_R_Elements/3.0_idp
TwoThirds = 2.0_idp*Third
Half      = Num_R_Elements/2
Fourth    = Half/2
UFourth   = Num_R_Elements - Fourth


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




END SUBROUTINE Tag_By_Half









!+101+################################################################!
!                                                                       !
!          Tag_Source_Elements                                          !
!                                                                       !
 !#####################################################################!
SUBROUTINE Tag_By_Source( Level,              &
                          BLo, BHi,           & ! Current Box Bounds
                          SLo, SHi,           & ! Src Bounds
                          TLo, THi,           & ! Tag Bounds
                          nComps,             &
                          Src,                & ! Pointer to Source Data
                          SetTag, ClearTag,   &
                          Tag                 )


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


INTEGER                                 ::  re, te, pe

INTEGER                                 ::  Num_DOF
REAL(idp)                               ::  E_Level_Threshold
REAL(idp)                               ::  E_Element_Max

REAL(idp), DIMENSION(0:10)               ::  Half_Decade_Table


Half_Decade_Table = [ 1,5,10,50,100,500,1000,5000,10000,50000,100000 ]


Num_DOF = nComps/5

!E_Level_Threshold = Central_E/10**(AMReX_Max_Level-Level)

E_Level_Threshold = Central_E/Half_Decade_Table(Level)

!PRINT*,"Level ",level," Threshold ",E_Level_Threshold," Central E ",Central_E

DO pe = BLo(3),BHi(3)
DO te = BLo(2),BHi(2)
DO re = BLo(1),BHi(1)
    
    E_Element_Max = MAXVAL( Src(re,te,pe,0:Num_DOF-1))
    
!    PRINT*,"E_Max ",E_Element_Max," Threshold ",E_Level_Threshold
    IF ( E_Element_Max > E_Level_Threshold ) THEN
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




END SUBROUTINE Tag_By_Source




END MODULE Driver_AMReX_Tagging_Module
