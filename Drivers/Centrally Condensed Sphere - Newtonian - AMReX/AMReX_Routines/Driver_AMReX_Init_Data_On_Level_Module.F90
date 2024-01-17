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

USE amrex_box_module, &
            ONLY :  amrex_box
            
USE amrex_boxarray_module, &
            ONLY :  amrex_boxarray,         &
                    amrex_boxarray_build,   &
                    amrex_boxarray_destroy
                    
USE amrex_distromap_module, &
            ONLY :  amrex_distromap,        &
                    amrex_distromap_build,  &
                    amrex_distromap_destroy
                    
USE amrex_multifab_module, &
            ONLY :  amrex_multifab,         &
                    amrex_multifab_build
#endif

USE Poseidon_Kinds_Module, &
            ONLY :  idp
            
USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag

USE Poseidon_Units_Module, &
            ONLY :  Centimeter,         &
                    Density_Units,      &
                    C_Square
                    
USE Parameters_Variable_Indices, &
            ONLY :  iS_E

USE Variables_External, &
            ONLY :  CCS_SurfaceRadius,              &
                    CCS_CoreRadius,                 &
                    CCS_CoreDensity


USE Variables_Quadrature, &
            ONLY :  INT_R_LOCATIONS,            &
                    NUM_R_QUAD_POINTS,          &
                    NUM_T_QUAD_POINTS,          &
                    NUM_P_QUAD_POINTS

USE Timer_Routines_Module, &
            ONLY :  TimerStart,                     &
                    TimerStop

USE Timer_VAriables_Module, &
            ONLY :  Timer_Core_Init_Test_Problem

USE Maps_Quadrature, &
            ONLY :  Quad_Map_Long_Array

USE Driver_Variables, &
            ONLY :  Driver_NQ,          &
                    Driver_RQ_xLocs,    &
                    Driver_xL
                    
USE Variables_Driver_AMReX, &
            ONLY :  xL,                 &
                    xR,                 &
                    MaxLevel,           &
                    nCells
            
USE External_CCS_Solution_Module, &
            ONLY :  CCS_Density

IMPLICIT NONE


CONTAINS




 !+101+################################################################!
!                                                                       !
!          Poseidon_Init_Data_On_Level                                  !
!                                                                       !
 !#####################################################################!
SUBROUTINE Poseidon_Init_Data_On_Level(   Level,          &
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
INTEGER                                     ::  Here


REAL(idp)                                   ::  x
REAL(idp)                                   ::  xwidth
REAL(idp)                                   ::  xloc
REAL(idp), DIMENSION(0:1)                   ::  xlocs
REAL(idp), DIMENSION(1:Driver_NQ(1))        ::  Cur_R_Locs
REAL(idp)                                   ::  DROT

REAL(idp)                                   ::  Density

101 FORMAT (a128)



CALL TimerStart( Timer_Core_Init_Test_Problem )


Num_DOF = nComps

xlocs(0) = -1.0_idp
xlocs(1) = +1.0_idp
xwidth  = Driver_xL(2)-Driver_xL(1)

DROT = (((xR(1)-xL(1))*Centimeter)/nCells(1))/xwidth

Src = 0.0_idp
DO pe = BLo(3),BHi(3)
DO te = BLo(2),BHi(2)
DO re = BLo(1),BHi(1)
    Cur_R_Locs(:) = DROT * (Driver_RQ_xlocs(:) - Driver_xL(1) + re*xwidth)
    DO rd = 1,Driver_NQ(1)
        Density = CCS_Density(Cur_R_Locs(rd))
        
        DO pd = 1,Driver_NQ(3)
        DO td = 1,Driver_NQ(2)
            
            here = Quad_Map_Long_Array(rd,td,pd,Driver_NQ)

            Src(re,te,pe,Here) = Density
            
            
        END DO ! td
        END DO ! pd
    END DO ! rd
END DO ! re
END DO ! te
END DO ! pe



CALL TimerStop( Timer_Core_Init_Test_Problem )



END SUBROUTINE Poseidon_Init_Data_On_Level






END MODULE Driver_Init_Data_On_Level_Module

