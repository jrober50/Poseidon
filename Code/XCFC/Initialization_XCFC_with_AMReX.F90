   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Initialization_XCFC_with_AMReX_Module                                 !##!
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
  amrex_multifab_build, &
  amrex_mfiter, &
  amrex_mfiter_build, &
  amrex_mfiter_destroy

USE Poseidon_Message_Routines_Module, &
            ONLY :  Init_Message

USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag

USE Allocation_XCFC_Linear_Systems, &
            ONLY :  Allocate_XCFC_Linear_Systems,   &
                    Reallocate_XCFC_Linear_Systems

USE Allocation_Mesh, &
            ONLY :  Allocate_Mesh,          &
                    Reallocate_Mesh

USE Initialization_Mesh_AMReX_Module, &
            ONLY :  Determine_AMReX_Mesh

USE Initialization_Derived, &
            ONLY :  Initialize_Derived_AMReX_Part2
            
USE Initialization_Tables, &
            ONLY :  Initialize_Level_Tables

USE Matrix_Initialization_Module, &
            ONLY :  Initialize_XCFC_Matrices

USE Poseidon_Remesh_Module, &
            ONLY :  Make_Remesh_Copies,             &
                    Fill_Coeff_Vector_From_Copy,    &
                    Destroy_Remesh_Copies

USE Maps_AMReX, &
            ONLY :  Initialize_AMReX_Maps,  &
                    Reinitialize_AMReX_Maps

USE IO_Setup_Report_Module, &
            ONLY :  PRINT_AMReX_Setup

USE Timer_Routines_Module, &
            ONLY :  TimerStart,                     &
                    TimerStop

USE Timer_Variables_Module, &
            ONLY :  Timer_Initialization_XCFC, &
                    Timer_Matrix_Init

USE Flags_Initialization_Module, &
            ONLY :  lPF_Init_Flags,     &
                    iPF_Init_Method_Vars
                    
USE Flags_IO_Module, &
            ONLY :  lPF_IO_Flags,               &
                    iPF_IO_Print_Setup


IMPLICIT NONE


CONTAINS



 !+101+############################################################!
!                                                                   !
!              Initialization_XCFC_with_AMReX                       !
!                                                                   !
 !#################################################################!
SUBROUTINE  Initialization_XCFC_with_AMReX( )


! Determine Radial Base Variables from Multifab
IF ( .NOT. lPF_Init_Flags(iPF_Init_Method_Vars) ) THEN
    
    IF ( Verbose_Flag ) CALL Init_Message('Initializing XCFC system with AMReX.')
    CALL TimerStart(Timer_Initialization_XCFC)

    
    CALL Initialize_AMReX_Maps()
    CALL Initialize_Level_Tables()
    CALL Initialize_Derived_AMReX_Part2()

    CALL Allocate_Mesh()
    CALL Determine_AMReX_Mesh()

    ! Allocate Arrays
    CALL Allocate_XCFC_Linear_Systems()

    ! Construct Matrices
    CALL TimerStart( Timer_Matrix_Init )
    CALL Initialize_XCFC_Matrices()
    CALL TimerStop( Timer_Matrix_Init )

    CALL TimerStop(Timer_Initialization_XCFC)

    lPF_Init_Flags(iPF_Init_Method_Vars) = .TRUE.

    IF ( lPF_IO_Flags(iPF_IO_Print_Setup) ) THEN
        CALL PRINT_AMReX_Setup()
    END IF
END IF

END SUBROUTINE  Initialization_XCFC_with_AMReX






 !+101+############################################################!
!                                                                   !
!              Reinitialization_XCFC_with_AMReX                       !
!                                                                   !
 !#################################################################!
SUBROUTINE  Reinitialization_XCFC_with_AMReX( )


! Determine Radial Base Variables from Multifab
IF ( lPF_Init_Flags(iPF_Init_Method_Vars) ) THEN

    lPF_Init_Flags(iPF_Init_Method_Vars) = .FALSE.
    
    IF ( Verbose_Flag ) CALL Init_Message('Initializing XCFC system with AMReX.')
    CALL TimerStart(Timer_Initialization_XCFC)

    CALL Make_Remesh_Copies()
    CALL Reinitialize_AMReX_Maps()
    CALL Initialize_Derived_AMReX_Part2()
    CALL Reallocate_Mesh()
    CALL Determine_AMReX_Mesh()
    CALL Reallocate_XCFC_Linear_Systems()

    
    
    CALL Fill_Coeff_Vector_From_Copy()
    CALL Destroy_Remesh_Copies

    ! Construct Matrices
    CALL TimerStart( Timer_Matrix_Init )
    CALL Initialize_XCFC_Matrices()
    CALL TimerStop( Timer_Matrix_Init )


    CALL TimerStop(Timer_Initialization_XCFC)

    lPF_Init_Flags(iPF_Init_Method_Vars) = .TRUE.

    IF ( Verbose_Flag ) THEN
        CALL PRINT_AMReX_Setup()
    END IF

END IF

END SUBROUTINE  Reinitialization_XCFC_with_AMReX

#endif



END MODULE Initialization_XCFC_with_AMReX_Module
