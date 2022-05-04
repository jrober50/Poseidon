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


USE Variables_Source, &
            ONLY :  Block_Source_E,             &
                    Block_Source_S,             &
                    Block_Source_Si

USE Variables_Quadrature, &
            ONLY :  INT_R_LOCATIONS,            &
                    NUM_R_QUAD_POINTS,          &
                    NUM_T_QUAD_POINTS,          &
                    NUM_P_QUAD_POINTS,          &
                    Num_Quad_DOF


USE Variables_Yahil, &
            ONLY :  SelfSim_T,              &
                    SelfSim_Kappa,          &
                    SelfSim_Gamma

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



INTEGER                                     ::  re, te, pe
INTEGER                                     ::  rd, td, pd
INTEGER                                     ::  Num_DOF


REAL(idp), DIMENSION(1:Num_R_Quad_Points)   ::  cur_r_locs
REAL(idp)                                   ::  DROT



CALL TimerStart( Timer_Core_Init_Test_Problem )



fofalpha            =  Alpha**5/(1.0_idp+Alpha*Alpha)**3
rho_o               =  (3.0_idp/(2.0_idp*pi*Star_Radius*Star_Radius) )*fofalpha*fofalpha


DROT = Level_dx(Level,1)/2.0_idp
Num_DOF = nComps/5

DO pe = BLo(3),BHi(3)
DO te = BLo(2),BHi(2)

line_min = 1
DO re = BLo(1),BHi(1)

    Cur_r_locs(:) = DROT * (Int_R_Locations(:) + 1.0_idp + re*2.0_idp)
    

    DO pd = 1,NQ(3)
    DO td = 1,NQ(2)
    DO rd = 1,NQ(1)


        Here = (pd-1)*NQ(1)*NQ(2)       &
             + (td-1)*NQ(1)             &
             + rd

        IF ( cur_r_locs(rq) .LE. Star_Radius ) THEN
            Local_E(re,te,pe,0*Num_DOF+Here) = rho_o
        ELSE
            Local_E(Here, re-1, te-1, pe-1) = 0.0_idp
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
