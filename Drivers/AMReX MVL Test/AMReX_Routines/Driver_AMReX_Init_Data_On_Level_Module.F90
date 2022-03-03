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

USE Timer_Routines_Module, &
            ONLY :  TimerStart,                     &
                    TimerStop

USE Timer_VAriables_Module, &
            ONLY :  Timer_Core_Init_Test_Problem

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


CALL TimerStart( Timer_Core_Init_Test_Problem )
PRINT*,"Here"
Src = 0.0_idp
PRINT*,"There"
CALL TimerStop( Timer_Core_Init_Test_Problem )


END SUBROUTINE Poseidon_Init_Data_On_Level











END MODULE Driver_Init_Data_On_Level_Module
