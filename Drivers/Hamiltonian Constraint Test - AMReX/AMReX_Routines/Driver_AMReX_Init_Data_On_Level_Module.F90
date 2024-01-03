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
            ONLY :  Centimeter
            
USE Parameters_Variable_Indices, &
            ONLY :  iU_CF

USE Driver_Variables, &
            ONLY :  Driver_NQ,          &
                    Driver_RQ_xLocs,    &
                    Driver_xL
                    
USE Variables_Driver_AMReX, &
            ONLY :  xL,                 &
                    xR,                 &
                    MaxLevel,           &
                    nCells

USE Variables_Quadrature, &
            ONLY :  INT_R_LOCATIONS,            &
                    NUM_R_QUAD_POINTS,          &
                    NUM_T_QUAD_POINTS,          &
                    NUM_P_QUAD_POINTS,          &
                    Local_Quad_DOF


USE Variables_External, &
            ONLY :  HCT_Alpha,              &
                    HCT_Star_Radius,        &
                    HCT_Rhoo

USE Timer_Routines_Module, &
            ONLY :  TimerStart,                     &
                    TimerStop

USE Timer_Variables_Module, &
            ONLY :  Timer_Core_Init_Test_Problem

USE Flags_Initial_Guess_Module, &
            ONLY :  lPF_IG_Flags,           &
                    iPF_IG_Set
                    
USE Poseidon_Return_Routines_Module, &
            ONLY :  Calc_Var_at_Location


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


REAL(idp)                                   ::  xwidth

REAL(idp), DIMENSION(1:Driver_NQ(1))        ::  cur_r_locs
REAL(idp)                                   ::  DROT

REAL(idp)                                   ::  Psi
REAL(idp)                                   ::  rho_o
INTEGER                                     ::  Here



CALL TimerStart( Timer_Core_Init_Test_Problem )

rho_o    =  HCT_Rhoo

xwidth   = Driver_xL(2)-Driver_xL(1)


DROT     = (((xR(1)-xL(1))*Centimeter)/(nCells(1)*2.0_idp**Level ))/xwidth
Num_DOF  = nComps/5


Src = 0.0_idp
DO pe = BLo(3),BHi(3)
DO te = BLo(2),BHi(2)

DO re = BLo(1),BHi(1)

    Cur_R_Locs(:) = DROT * (Driver_RQ_xlocs(:) - Driver_xL(1) + re*xwidth)
    

    DO pd = 1,Driver_NQ(3)
    DO td = 1,Driver_NQ(2)
    DO rd = 1,Driver_NQ(1)


        Here = (pd-1)*Num_R_Quad_Points*Num_T_Quad_Points       &
             + (td-1)*Num_R_Quad_Points                         &
             + rd

        
        IF ( Cur_R_Locs(rd) .LE. HCT_Star_Radius ) THEN
            IF ( lPF_IG_Flags(iPF_IG_Set) ) THEN
                Psi = Calc_Var_at_Location(Cur_R_Locs(rd),0.0_idp,0.0_idp, iU_CF)
                Src(re,te,pe,Here) = Psi**6 * rho_o
            ELSE
                Src(re,te,pe,Here) = rho_o
            END IF
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
