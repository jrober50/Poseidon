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
            ONLY :  MacLaurin_SemiMinor,        &
                    MacLaurin_SemiMajor,        &
                    MacLaurin_RE,               &
                    MacLaurin_SphereType


USE Variables_Quadrature, &
            ONLY :  INT_R_LOCATIONS,            &
                    Int_T_Locations,            &
                    Int_P_Locations,            &
                    NUM_R_QUAD_POINTS,          &
                    NUM_T_QUAD_POINTS,          &
                    NUM_P_QUAD_POINTS,          &
                    Local_Quad_DOF

USE Variables_AMReX_Core, &
            ONLY :  AMReX_Max_Level

USE Variables_Tables,   &
            ONLY :  Level_DX

USE Maps_Quadrature, &
            ONLY :  Quad_Map


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

INTEGER                                 ::  Tag_Style = 2


! Be careful with Tagging. You have to tag in radial shells for now.




IF ( Tag_Style == 1 ) THEN

    CALL Tag_By_Intersection( Level,              &
                              BLo, BHi,           & ! Current Box Bounds
                              SLo, SHi,           & ! Src Bounds
                              TLo, THi,           & ! Tag Bounds
                              nComps,             &
                              Src,                & ! Pointer to Source Data
                              SetTag, ClearTag,   &
                              Tag                 )



ELSE IF ( Tag_Style == 2 ) THEN

    CALL Tag_To_Core( Level,              &
                        BLo, BHi,           & ! Current Box Bounds
                        SLo, SHi,           & ! Src Bounds
                        TLo, THi,           & ! Tag Bounds
                        nComps,             &
                        Src,                & ! Pointer to Source Data
                        SetTag, ClearTag,   &
                        Tag                 )



ELSE IF ( Tag_Style == 3 ) THEN

    CALL Tag_By_Bracket( Level,              &
                        BLo, BHi,           & ! Current Box Bounds
                        SLo, SHi,           & ! Src Bounds
                        TLo, THi,           & ! Tag Bounds
                        nComps,             &
                        Src,                & ! Pointer to Source Data
                        SetTag, ClearTag,   &
                        Tag                 ) ! Pointer to Tag Data


END IF






END SUBROUTINE Tag_Source_Elements








 !+101+################################################################!
!                                                                       !
!          Tag_Source_Elements                                          !
!                                                                       !
 !#####################################################################!
SUBROUTINE Tag_By_Intersection( Level,              &
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

INTEGER                                 ::  Tag_Style = 2

INTEGER                                 ::  re, te, pe



REAL(idp)                               :: SemiMajor_Axis
REAL(idp)                               :: SemiMinor_Axis

REAL(idp)                               :: DROT, DTOT, DPOT
REAL(idp)                               :: A, B, C
REAL(idp)                               :: AA, BB, CC, AB

REAL(idp)                               :: R_Hi
REAL(idp)                               :: R_Low

REAL(idp)                               :: Radius_L
REAL(idp)                               :: Radius_R

REAL(idp)                               :: Theta_L
REAL(idp)                               :: Theta_R

REAL(idp)                               :: SinSqr_L
REAL(idp)                               :: SinSqr_R
REAL(idp)                               :: CosSqr_L
REAL(idp)                               :: CosSqr_R

LOGICAL                                 :: Cond1
LOGICAL                                 :: Cond2


DROT = Level_dx(Level,1)/2.0_idp
DTOT = Level_dx(Level,2)/2.0_idp
DPOT = Level_dx(Level,3)/2.0_idp


SemiMajor_Axis = MacLaurin_SemiMajor
SemiMinor_Axis = MacLaurin_SemiMinor

IF ( MacLaurin_SphereType == 'P') THEN
!    Spheroid_Type_Flag  = 2
!    Spheroid_Name       = 'Prolate'

    A = SemiMajor_Axis
    B = SemiMinor_Axis
    C = B
    
ELSE
!    Spheroid_Type_Flag  = 1
!    Spheroid_Name       = 'Oblate '

    A = SemiMajor_Axis
    B = A
    C = SemiMinor_Axis

END IF


AA = A*A
BB = B*B
CC = C*C


AB = A*B



DO pe = BLo(3),BHi(3)
DO te = BLo(2),BHi(2)
DO re = BLo(1),BHi(1)
    

    ! The lower and upper radii of the current element
    R_Low   = 2.0_idp*DROT*re
    R_Hi    = 2.0_idp*DROT*(re+1)

    ! The left and right theta values of the current element.
    Theta_L = 2.0_idp*DTOT*te
    Theta_R = 2.0_idp*DTOT*(te+1)

    SinSqr_L = DSIN(Theta_L)*DSIN(Theta_L)
    SinSqr_R = DSIN(Theta_R)*DSIN(Theta_R)

    CosSqr_L = DCOS(Theta_L)*DCOS(Theta_L)
    CosSqr_R = DCOS(Theta_R)*DCOS(Theta_R)


    ! Calculate the radius of the spheroid at the left and right theta values
    Radius_L = AB /sqrt(BB*SinSqr_L + AA*CosSqr_L)
    Radius_R = AB /sqrt(BB*SinSqr_R + AA*CosSqr_R)


    ! Does the spheroid radius intersect the left side of the element?
    IF ( (R_Low .LE. Radius_L) .and. (Radius_L .LE. R_Hi ) ) THEN
        Cond1 = .TRUE.
    ELSE
        Cond1 = .FALSE.
    END IF

    ! Does the spheroid radius intersect the right side of the element?
    IF ( (R_Low .LE. Radius_R) .and. (Radius_R .LE. R_Hi ) ) THEN
        Cond2 = .TRUE.
    ELSE
        Cond2 = .FALSE.
    END IF


    ! If either of the above are true, the element is tagged.
    IF ( Cond1 .or. Cond2 ) THEN
!        PRINT*,"Refining",level,re,te
        Tag(re,te,pe,1) = SetTag
    ELSE
        Tag(re,te,pe,1) = ClearTag
    END IF




END DO ! re
END DO ! te
END DO ! pe





END SUBROUTINE Tag_By_Intersection








 !+101+################################################################!
!                                                                       !
!          Tag_Source_Elements                                          !
!                                                                       !
 !#####################################################################!
SUBROUTINE Tag_To_Core( Level,              &
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

INTEGER                                 ::  Tag_Style = 2

INTEGER                                 ::  re, te, pe
INTEGER                                 ::  rd, td, pd
INTEGER                                 ::  Here

REAL(idp)                               :: SemiMajor_Axis
REAL(idp)                               :: SemiMinor_Axis

REAL(idp)                               :: DROT, DTOT, DPOT

REAL(idp)                                   ::  A,  B,  C
REAL(idp)                                   ::  AA, BB, CC, AB

REAL(idp), DIMENSION(1:Num_R_Quad_Points)   ::  Cur_R_Locs
REAL(idp), DIMENSION(1:Num_T_Quad_Points)   ::  Cur_T_Locs
REAL(idp), DIMENSION(1:Num_P_Quad_Points)   ::  Cur_P_Locs

REAL(idp), DIMENSION(1:Num_R_Quad_Points)   ::  RSqr

REAL(idp), DIMENSION(1:Num_T_Quad_Points)   ::  SinSqr_T
REAL(idp), DIMENSION(1:Num_T_Quad_Points)   ::  CosSqr_T

REAL(idp), DIMENSION(1:Num_P_Quad_Points)   ::  SinSqr_P
REAL(idp), DIMENSION(1:Num_P_Quad_Points)   ::  CosSqr_P


REAL(idp)                               :: Value

LOGICAL                                 :: Cond1
LOGICAL                                 :: Cond2
LOGICAL, DIMENSION(:,:,:), ALLOCATABLE  :: Values


ALLOCATE( Values(1:Num_R_Quad_Points, 1:Num_T_Quad_Points, 1:Num_P_Quad_Points) )

DROT = Level_dx(Level,1)/2.0_idp
DTOT = Level_dx(Level,2)/2.0_idp
DPOT = Level_dx(Level,3)/2.0_idp


SemiMajor_Axis = MacLaurin_SemiMajor
SemiMinor_Axis = MacLaurin_SemiMinor

IF ( MacLaurin_SphereType == 'P') THEN
!    Spheroid_Type_Flag  = 2
!    Spheroid_Name       = 'Prolate'

    A = SemiMajor_Axis
    B = SemiMinor_Axis
    C = B
    
ELSE
!    Spheroid_Type_Flag  = 1
!    Spheroid_Name       = 'Oblate '

    A = SemiMajor_Axis
    B = A
    C = SemiMinor_Axis

END IF




DO pe = BLo(3),BHi(3)
DO te = BLo(2),BHi(2)
DO re = BLo(1),BHi(1)
    



    IF ( 2.0_idp*DROT*(RE+1) .LE. SemiMajor_Axis ) THEN
!        PRINT*,"Refining",level,re,te
        Tag(re,te,pe,1) = SetTag
    ELSE
        Tag(re,te,pe,1) = ClearTag
    END IF


END DO ! re
END DO ! te
END DO ! pe





END SUBROUTINE Tag_to_Core






 !+101+################################################################!
!                                                                       !
!          Tag_Source_Elements                                          !
!                                                                       !
 !#####################################################################!
SUBROUTINE Tag_By_Bracket( Level,              &
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

INTEGER                                 ::  Tag_Style = 2

INTEGER                                 ::  re, te, pe
INTEGER                                 ::  rd, td, pd
INTEGER                                 ::  Here

REAL(idp)                               ::  SemiMajor_Axis
REAL(idp)                               ::  SemiMinor_Axis

REAL(idp)                               ::  DROT

REAL(idp)                               ::  Inner_Edge
REAL(idp)                               ::  Outer_Edge

DROT = Level_dx(Level,1)/2.0_idp

SemiMajor_Axis = MacLaurin_SemiMajor
SemiMinor_Axis = MacLaurin_SemiMinor


DO pe = BLo(3),BHi(3)
DO te = BLo(2),BHi(2)
DO re = BLo(1),BHi(1)
    
    Inner_Edge = 2.0_idp*DROT*RE
    Outer_Edge = 2.0_idp*DROT*(RE+1)


    IF (( Outer_Edge .LE. SemiMajor_Axis )          &
                    .AND.                           &
        (Inner_Edge .GE. SemiMinor_Axis  ) )   THEN
!        PRINT*,"Refining",level,re,te
        Tag(re,te,pe,1) = SetTag
    ELSE
        Tag(re,te,pe,1) = ClearTag
    END IF


END DO ! re
END DO ! te
END DO ! pe





END SUBROUTINE Tag_By_Bracket



END MODULE Driver_AMReX_Tagging_Module
