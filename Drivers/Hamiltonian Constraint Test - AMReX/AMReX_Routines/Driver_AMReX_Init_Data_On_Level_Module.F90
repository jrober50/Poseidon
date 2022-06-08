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

USE Poseidon_Numbers_Module, &
            ONLY :  pi

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


USE Variables_Source, &
            ONLY :  Block_Source_E,             &
                    Block_Source_S,             &
                    Block_Source_Si

USE Variables_Quadrature, &
            ONLY :  INT_R_LOCATIONS,            &
                    NUM_R_QUAD_POINTS,          &
                    NUM_T_QUAD_POINTS,          &
                    NUM_P_QUAD_POINTS,          &
                    Local_Quad_DOF


USE Variables_External, &
            ONLY :  HCT_Alpha,                  &
                    HCT_Star_Radius

USE Poseidon_AMReX_Multilayer_Utilities_Module, &
            ONLY :  Find_Coarsest_Parent

USE Variables_Tables, &
            ONLY :  Level_dx

USE Functions_Math, &
            ONLY :  Lagrange_Poly

USE Timer_Routines_Module, &
            ONLY :  TimerStart,                     &
                    TimerStop

USE Timer_VAriables_Module, &
            ONLY :  Timer_Core_Init_Test_Problem

USE Variables_Functions, &
            ONLY :  Potential_Solution

USE Maps_Quadrature, &
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



INTEGER                                     ::  re, te, pe
INTEGER                                     ::  rd, td, pd
INTEGER                                     ::  Num_DOF


REAL(idp), DIMENSION(1:Num_R_Quad_Points)   ::  cur_r_locs
REAL(idp)                                   ::  DROT

REAL(idp)                                   ::  fofalpha
REAL(idp)                                   ::  rho_o
INTEGER                                     ::  Here


CALL TimerStart( Timer_Core_Init_Test_Problem )



fofalpha            =  HCT_Alpha**5/(1.0_idp+HCT_Alpha*HCT_Alpha)**3
rho_o               =  (3.0_idp/(2.0_idp*pi*HCT_Star_Radius*HCT_Star_Radius) )*fofalpha*fofalpha


DROT = Level_dx(Level,1)/2.0_idp
Num_DOF = nComps/5

Src = 0.0_idp

DO pe = BLo(3),BHi(3)
DO te = BLo(2),BHi(2)

DO re = BLo(1),BHi(1)

    Cur_r_locs(:) = DROT * (Int_R_Locations(:) + 1.0_idp + re*2.0_idp)
    

    DO pd = 1,Num_P_Quad_Points
    DO td = 1,Num_T_Quad_Points
    DO rd = 1,Num_R_Quad_Points


        Here = (pd-1)*Num_R_Quad_Points*Num_T_Quad_Points       &
             + (td-1)*Num_R_Quad_Points                         &
             + rd

        IF ( cur_r_locs(rd) .LE. HCT_Star_Radius ) THEN
            Src(re,te,pe,Here) = rho_o
        ELSE
            Src(re,te,pe,Here) = 0.0_idp
        END IF

    END DO ! rq
    END DO ! tq
    END DO ! pq

END DO ! re
END DO ! te
END DO ! pe









CALL TimerStop( Timer_Core_Init_Test_Problem )


END SUBROUTINE Poseidon_Init_Data_On_Level











END MODULE Driver_Init_Data_On_Level_Module
