   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Driver_SetSource_Module                                              !##!
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
USE amrex_amrcore_module, ONLY: &
  amrex_amrcore_init, &
  amrex_init_virtual_functions, &
  amrex_init_from_scratch, &
  amrex_ref_ratio, &
  amrex_max_level


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
                    AMReX_Num_Levels


USE Variables_Driver_AMReX, &
            ONLY :  MF_Driver_Source,   &
                    MF_Src_nComps,      &
                    MF_Src_nGhost,      &
                    dt, t_old, t_new,   &
                    StepNo,             &
                    nLevels

USE Poseidon_AMReX_MakeFineMask_Module, &
            ONLY :  AMReX_MakeFineMask

USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag


USE Variables_MPI, &
            ONLY :  myID_Poseidon,      &
                    Poseidon_Comm_World

USE Source_Input_AMReX, &
            ONLY :  Poseidon_Input_Sources_AMREX


USE Poseidon_MPI_Utilities_Module, &
            ONLY :  STOP_MPI,               &
                    MPI_Master_Print,       &
                    MPI_All_Print

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

USE Functions_Quadrature, &
            ONLY :  Initialize_LG_Quadrature_Locations,         &
                    Initialize_Trapezoid_Quadrature_Locations


USE MPI

IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!     Driver_SetSource                                                			!
!                                                                               !
!###############################################################################!
SUBROUTINE Driver_SetSource( NQ )


INTEGER,    INTENT(IN), DIMENSION(3)                    ::  NQ


REAL(idp),  DIMENSION(NQ(1))                            ::  R_Quad
REAL(idp),  DIMENSION(NQ(2))                            ::  T_Quad
REAL(idp),  DIMENSION(NQ(3))                            ::  P_Quad
REAL(idp),  DIMENSION(2)                                ::  xL

INTEGER                                                 ::  Num_DOF
INTEGER                                                 ::  nVars_Source


IF ( Verbose_Flag ) THEN
    WRITE(*,'(A)')"In Driver, Creating AMReX Source Variables."
END IF


CALL TimerStart( Timer_Driver_SetSource_InitTest )


CALL amrex_init_virtual_functions &
       ( VF_Make_New_Level_From_Scratch, &
         VF_Make_New_Level_From_Coarse, &
         VF_Remake_Level, &
         VF_Clear_Level, &
         VF_Error_Estimate )

Num_DOF = NQ(1)*NQ(2)*NQ(3)

nVars_Source    = 5
MF_Src_nComps   = nVars_Source*Num_DOF
MF_Src_nGhost   = 0


ALLOCATE( MF_Driver_Source(0:nLevels-1) )

CALL amrex_init_from_scratch( 0.0_idp )



IF ( Verbose_Flag ) THEN
    WRITE(*,'(A)')"In Driver, Inputing AMReX Source Variables."
END IF

xL(1) = -1.0_idp
xL(2) = +1.0_idp
R_Quad = Initialize_LG_Quadrature_Locations(NQ(1))
T_Quad = Initialize_LG_Quadrature_Locations(NQ(2))
P_Quad = Initialize_Trapezoid_Quadrature_Locations(NQ(3))

CALL Poseidon_Input_Sources_AMREX( MF_Driver_Source,    &
                                   MF_Src_nComps,       &
                                   nLevels,             &
                                   NQ,                  &
                                   R_Quad,              &
                                   T_Quad,              &
                                   P_Quad,              &
                                   xL                   )


CALL TimerStop( Timer_Driver_SetSource_InitTest )





DEALLOCATE(MF_Driver_Source)



END SUBROUTINE Driver_SetSource






























END MODULE Driver_SetSource_Module
