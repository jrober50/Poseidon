   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Driver_InitSource_Module                                              !##!
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
USE Poseidon_Kinds_Module, &
            ONLY :  idp
            
USE Poseidon_Message_Routines_Module, &
            ONLY :  Driver_Init_Message

USE amrex_base_module

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
USE amrex_amrcore_module, ONLY: &
  amrex_amrcore_init, &
  amrex_init_virtual_functions, &
  amrex_init_from_scratch, &
  amrex_ref_ratio, &
  amrex_max_level


USE Poseidon_Units_Module, &
            ONLY :  Centimeter,         &
                    Density_Units,      &
                    Solar_Radius

USE Variables_AMReX_Core, &
            ONLY :  MF_Source,          &
                    AMReX_Num_Levels,   &
                    Source_PTR,         &
                    Mask_PTR

USE Parameters_AMReX, &
            ONLY :  iLeaf,            &
                    iTrunk

USE Variables_Driver_AMReX, &
            ONLY :  MF_Driver_Source,   &
                    MF_Src_nComps,      &
                    MF_Src_nGhost

USE Poseidon_AMReX_MakeFineMask_Module, &
            ONLY :  AMReX_MakeFineMask
            
USE Variables_External, &
            ONLY :  CCS_SurfaceRadius,              &
                    CCS_CoreRadius,                 &
                    CCS_CoreDensity

USE Driver_AMReX_Virtual_Functions_Module, &
            ONLY :  VF_Make_New_Level_From_Scratch, &
                    VF_Make_New_Level_From_Coarse,  &
                    VF_Remake_Level,                &
                    VF_Clear_Level,                 &
                    VF_Error_Estimate

USE Timer_Routines_Module, &
            ONLY :  TimerStart,     &
                    TimerSTop

USE Timer_Variables_Module, &
            ONLY :  Timer_Driver_SetSource_InitTest,        &
                    Timer_Driver_SetSource_SetSource,       &
                    Timer_Driver_SetSource_Scale

USE Driver_Variables, &
            ONLY :  Driver_NQ

USE MPI

IMPLICIT NONE


CONTAINS



 !+101+########################################################!
!                                                               !
!     Driver_InitSource                                         !
!                                                               !
 !#############################################################!
SUBROUTINE Driver_InitSource( Verbose_Flag )

LOGICAL,    INTENT(IN)              ::  Verbose_Flag
INTEGER                             ::  nVars_Source


IF ( Verbose_Flag ) CALL Driver_Init_Message('Creating AMReX source variables.')
IF ( Verbose_Flag ) CALL Driver_Init_Message('Initializing Centrally Condensed Sphere Source Multifab.')

CALL TimerStart( Timer_Driver_SetSource_InitTest )




CCS_SurfaceRadius = 1.0_idp * Solar_Radius

CCS_CoreRadius    = 0.20_idp * Solar_Radius
CCS_CoreDensity   = 150 * Density_Units


IF ( Verbose_Flag ) THEN
    WRITE(*,'(A)')'------------- Test Parameters ----------------'
    WRITE(*,'(A)')' Source Configuration : Centrally Condensed Sphere'
    WRITE(*,'(A,ES12.5,A)') ' - Surface Radius        : ', CCS_SurfaceRadius/Centimeter," cm"
    WRITE(*,'(A,ES12.5,A)') ' - Core Radius           : ', CCS_CoreRadius/Centimeter," cm"
    WRITE(*,'(A,ES12.5,A)') ' - Core Density          : ', CCS_CoreDensity/Density_Units," g/cm^3"
    WRITE(*,'(/)')
END IF




CALL amrex_init_virtual_functions &
       ( VF_Make_New_Level_From_Scratch, &
         VF_Make_New_Level_From_Coarse, &
         VF_Remake_Level, &
         VF_Clear_Level, &
         VF_Error_Estimate )



nVars_Source    = 1
MF_Src_nComps   = nVars_Source*Driver_NQ(1)*Driver_NQ(2)*Driver_NQ(3)
MF_Src_nGhost   = 0


ALLOCATE( MF_Driver_Source(0:amrex_max_level) )
CALL amrex_init_from_scratch( 0.0_idp )
CALL TimerStop( Timer_Driver_SetSource_InitTest )


END SUBROUTINE Driver_InitSource









END MODULE Driver_InitSource_Module

