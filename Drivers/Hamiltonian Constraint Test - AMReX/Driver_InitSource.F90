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

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag

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
  amrex_max_level

USE amrex_amrcore_module, &
            ONLY :  amrex_get_numlevels

USE Poseidon_Units_Module, &
            ONLY :  Grav_Constant_G,    &
                    Speed_of_Light,     &
                    C_Square,           &
                    GR_Source_Scalar,   &
                    Centimeter,         &
                    Second,             &
                    Millisecond,         &
                    Erg,                &
                    Gram,               &
                    E_Units


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
                    MF_Src_nGhost,      &
                    dt, t_old, t_new,   &
                    StepNo,             &
                    nLevels,            &
                    MaxLevel



USE Poseidon_AMReX_MakeFineMask_Module, &
            ONLY :  AMReX_MakeFineMask

USE Variables_MPI, &
            ONLY :  myID_Poseidon,      &
                    Poseidon_Comm_World

USE Poseidon_Interface_Source_Input, &
            ONLY :  Poseidon_Input_Sources


USE Poseidon_MPI_Utilities_Module, &
            ONLY :  STOP_MPI,               &
                    MPI_Master_Print,       &
                    MPI_All_Print

USE Variables_External, &
            ONLY :  HCT_Alpha,              &
                    HCT_Star_Radius

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
SUBROUTINE Driver_InitSource( iter_option )

INTEGER, intent(IN), OPTIONAL   ::  iter_option
INTEGER                         ::  nVars_Source
INTEGER                         ::  Iter

IF ( Verbose_Flag ) CALL Driver_Init_Message('Creating AMReX source variables.')
IF ( Verbose_Flag ) CALL Driver_Init_Message('Initializing Hamiltonian Constraint Test Source Multifab.')


CALL TimerStart( Timer_Driver_SetSource_InitTest )



IF ( Verbose_Flag ) THEN
    WRITE(*,'(A)')'------------- Test Parameters ----------------'
    WRITE(*,'(A)')' Source Configuration : Hamiltonian Constraint Test'
    WRITE(*,'(A,ES12.5,A)') ' - Alpha       : ', HCT_Alpha
    WRITE(*,'(A,ES12.5)')   ' - Star Radius : ', HCT_Star_Radius
    WRITE(*,'(/)')
END IF


!PRINT*,"Levels",nLevels,AMReX_Num_Levels,amrex_get_numlevels(),MaxLevel

IF (present(iter_option) ) THEN
    iter = iter_option
else
    iter = 1
end if

IF ( iter == 1 ) THEN
    CALL amrex_init_virtual_functions &
           ( VF_Make_New_Level_From_Scratch, &
             VF_Make_New_Level_From_Coarse, &
             VF_Remake_Level, &
             VF_Clear_Level, &
             VF_Error_Estimate )
end IF

nVars_Source    = 5
MF_Src_nComps   = nVars_Source*Driver_NQ(1)*Driver_NQ(2)*Driver_NQ(3)
MF_Src_nGhost   = 0

!PRINT*,"nLevels",nLevels
IF ( .NOT. Allocated(MF_Driver_Source) ) THEN
    ALLOCATE( MF_Driver_Source(0:nLevels-1) )
END IF

CALL amrex_init_from_scratch( 0.0_idp )
CALL TimerStop( Timer_Driver_SetSource_InitTest )


END SUBROUTINE Driver_InitSource






END MODULE Driver_InitSource_Module

