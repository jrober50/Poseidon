   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Driver_Init_Data_On_Level_Module                                    !##!
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

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag

USE Poseidon_Units_Module, &
            ONLY :  Grav_Constant_G,    &
                    Speed_of_Light,     &
                    C_Square,           &
                    GR_Source_Scalar,   &
                    Centimeter,         &
                    Second,             &
                    Millisecond,         &
                    Erg,                &
                    Gram

USE Variables_External, &
            ONLY :  MacLaurin_SemiMinor,    &
                    MacLaurin_SemiMajor,    &
                    MacLaurin_Ecc,          &
                    MacLaurin_SphereType,   &
                    MacLaurin_Rho



USE Variables_Quadrature, &
            ONLY :  INT_R_LOCATIONS,            &
                    Int_T_Locations,            &
                    Int_P_Locations,            &
                    NUM_R_QUAD_POINTS,          &
                    NUM_T_QUAD_POINTS,          &
                    NUM_P_QUAD_POINTS,          &
                    Num_Quad_DOF


USE Variables_Tables,   &
            ONLY :  Level_DX

USE Timer_Routines_Module, &
            ONLY :  TimerStart,                     &
                    TimerStop

USE Timer_Variables_Module, &
            ONLY :  Timer_Core_Init_Test_Problem

USE Variables_Functions, &
            ONLY :  Potential_Solution

USE Quadrature_Mapping_Functions, &
            ONLY :  Quad_Map

IMPLICIT NONE


CONTAINS




 !+101+################################################################!
!                                                                       !
!          Poseidon_Init_Data_On_Level                                  !
!                                                                       !
 !#####################################################################!
SUBROUTINE Poseidon_Init_Data_On_Level( Level,          &
                                        BLo, BHi,       &
                                        SLo, SHi,       &
                                        nComps,         &
                                        Src             )

INTEGER,                INTENT(IN)          ::  Level
INTEGER,                INTENT(IN)          ::  BLo(3), BHi(3)
INTEGER,                INTENT(IN)          ::  SLo(3), SHi(3)
INTEGER,                INTENT(IN)          ::  nComps
REAL(idp),              INTENT(INOUT)       ::  Src(SLo(1):SHi(1),  &
                                                    SLo(2):SHi(2),  &
                                                    SLo(3):SHi(3),  &
                                                    nComps          )

INTEGER                                     ::  Here
INTEGER                                     ::  re, te, pe
INTEGER                                     ::  rd, td, pd

REAL(idp)                                   ::  DROT
REAL(idp)                                   ::  DTOT
REAL(idp)                                   ::  DPOT

REAL(idp)                                   ::  SemiMajor_Axis
REAL(idp)                                   ::  SemiMinor_Axis

REAL(idp)                                   ::  A,  B,  C
REAL(idp)                                   ::  AA, BB, CC

REAL(idp)                                   ::  Density
REAL(idp)                                   ::  Pressure
REAL(idp)                                   ::  Energy
REAL(idp)                                   ::  Spec_Ent
REAL(idp)                                   ::  Value

REAL(idp), DIMENSION(1:Num_R_Quad_Points)   ::  Cur_R_Locs
REAL(idp), DIMENSION(1:Num_T_Quad_Points)   ::  Cur_T_Locs
REAL(idp), DIMENSION(1:Num_P_Quad_Points)   ::  Cur_P_Locs

REAL(idp), DIMENSION(1:Num_R_Quad_Points)   ::  RSqr

REAL(idp), DIMENSION(1:Num_T_Quad_Points)   ::  SinSqr_T
REAL(idp), DIMENSION(1:Num_T_Quad_Points)   ::  CosSqr_T

REAL(idp), DIMENSION(1:Num_P_Quad_Points)   ::  SinSqr_P
REAL(idp), DIMENSION(1:Num_P_Quad_Points)   ::  CosSqr_P
CHARACTER(LEN=7)                            :: Spheroid_Name


CALL TimerStart( Timer_Core_Init_Test_Problem )

Src = 0.0_idp





SemiMajor_Axis = MacLaurin_SemiMajor
SemiMinor_Axis = MacLaurin_SemiMinor

IF ( MacLaurin_SphereType == 'P') THEN
!    Spheroid_Type_Flag  = 2
    Spheroid_Name       = 'Prolate'

    A = SemiMajor_Axis
    B = SemiMinor_Axis
    C = B
    
ELSE
!    Spheroid_Type_Flag  = 1
    Spheroid_Name       = 'Oblate '

    A = SemiMajor_Axis
    B = A
    C = SemiMinor_Axis

END IF


AA = A*A
BB = B*B
CC = C*C

Density  = MacLaurin_Rho
Pressure = 0.0_idp
Energy   = 0.0_idp


!Pressure = kappa * Density**Gamma
!Energy   = Pressure / (Gamma - 1.0_dip )
Spec_Ent = C_Square + (Energy + Pressure)/Density



DROT = Level_dx(Level,1)/2.0_idp
DTOT = Level_dx(Level,2)/2.0_idp
DPOT = Level_dx(Level,3)/2.0_idp







DO pe = BLo(3),BHi(3)
DO te = BLo(2),BHi(2)
DO re = BLo(1),BHi(1)

!    PRINT*,DROT*2.0_idp*re
    Cur_R_locs(:) = DROT * (Int_R_Locations(:) + 1.0_idp + re*2.0_idp)
    Cur_T_locs(:) = DTOT * (Int_T_Locations(:) + 1.0_idp + te*2.0_idp)
    Cur_P_locs(:) = DPOT * (Int_P_Locations(:) + 1.0_idp + pe*2.0_idp)


    RSqr(:)     = Cur_R_Locs(:) * Cur_R_Locs(:)
    SinSqr_t(:) = SIN( Cur_T_Locs(:) ) * SIN( CUR_T_LOCS(:) )
    CosSqr_t(:) = COS( Cur_T_Locs(:) ) * COS( CUR_T_LOCS(:) )
    SinSqr_p(:) = SIN( Cur_P_Locs(:) ) * SIN( CUR_P_LOCS(:) )
    CosSqr_p(:) = COS( Cur_P_Locs(:) ) * COS( CUR_P_LOCS(:) )


    DO rd = 1,Num_R_Quad_Points
    DO td = 1,Num_T_Quad_Points
    DO pd = 1,Num_P_Quad_Points

        Here = Quad_Map(rd, td, pd )


        Value = ( rsqr(rd) * cossqr_p(pd) * sinsqr_t(td) ) / AA      &
              + ( rsqr(rd) * sinsqr_p(pd) * sinsqr_t(td) ) / BB      &
              + ( rsqr(rd) * cossqr_t(td) ) / CC


        IF ( Value .LE. 1.0_idp ) THEN
            Src(re,te,pe,Here) = Density * Spec_Ent - Pressure
        END IF

!        PRINT*,re,te,pe,Here,Src(re,te,pe,Here)

    END DO ! pd
    END DO ! td
    END DO ! rd

END DO ! re
END DO ! te
END DO ! pe


CALL TimerStop( Timer_Core_Init_Test_Problem )


END SUBROUTINE Poseidon_Init_Data_On_Level











END MODULE Driver_Init_Data_On_Level_Module
