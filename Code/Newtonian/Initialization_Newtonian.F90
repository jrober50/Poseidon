   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Initialization_Poisson                                                       !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains the top level subroutines needed to inialize, run, and close       !##!
!##!        Poseidon.                                                               !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!    +101+   Poseidon_Initialize                                                 !##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!


!*D*================================!
!                                   !
!           Dependencies            !
!                                   !
!===================================!
USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Poseidon_Message_Routines_Module, &
            ONLY :  Init_Message

USE Poseidon_Parameters, &
            ONLY :  Degree,                 &
                    Verbose_Flag

USE Variables_Mesh, &
            ONLY :  Num_R_Elements

USE Variables_Matrices, &
            ONLY :  Laplace_NNZ

USE Allocation_Poisson_Linear_System, &
            ONLY :  Allocate_Poisson_Linear_System,     &
                    Reallocate_Poisson_Linear_System
            
USE Allocation_Mesh, &
            ONLY :  Allocate_Mesh,          &
                    Reallocate_Mesh
                    
USE Initialization_Mesh_AMReX_Module, &
            ONLY :  Determine_AMReX_Mesh
            
USE Initialization_Derived, &
            ONLY :  Initialize_Derived_AMReX_Part2
            
USE Initialization_Tables, &
            ONLY :  Initialize_Level_Tables
            
USE Maps_AMReX, &
            ONLY :  Initialize_AMReX_Maps,  &
                    Reinitialize_AMReX_Maps

USE Matrix_Initialization_Module, &
            ONLY :  Initialize_Poisson_Matrix
                    
USE IO_Setup_Report_Module, &
            ONLY :  PRINT_AMReX_Setup

USE Timer_Routines_Module, &
            ONLY :  TimerStart,                     &
                    TimerStop

USE Timer_Variables_Module, &
            ONLY :  Timer_Newtonian_Initialization,      &
                    Timer_Newtonian_Matrix_Init

USE Flags_Initialization_Module, &
            ONLY :  lPF_Init_Flags,     &
                    iPF_Init_Method_Vars
                    
USE Flags_IO_Module, &
            ONLY :  lPF_IO_Flags,               &
                    iPF_IO_Print_Setup
              
#ifdef POSEIDON_AMREX_FLAG
USE Poseidon_Remesh_Module, &
            ONLY :  Make_Remesh_Copies,             &
                    Fill_Coeff_Vector_From_Copy,    &
                    Destroy_Remesh_Copies
#endif
IMPLICIT NONE




                    !*F&S*==========================================!
                    !                                               !
                    !           Functions & Subroutines             !
                    !                                               !
                    !===============================================!
CONTAINS










 !+101+################################################################!
!                                                                       !
!          Initialize_Poisson                                           !
!                                                                       !
!=======================================================================!
!                                                                       !
 !#####################################################################!
SUBROUTINE Initialize_Poisson( )


IF ( .NOT. lPF_Init_Flags(iPF_Init_Method_Vars) ) THEN


    IF ( Verbose_Flag ) CALL Init_Message('Initializing XCFC system variables.')

    CALL TimerStart( Timer_Newtonian_Initialization )


    Laplace_NNZ = NUM_R_ELEMENTS*(DEGREE + 1)*(DEGREE + 1) - NUM_R_ELEMENTS + 1

    CALL Allocate_Poisson_Linear_System()

    CALL TimerStart( Timer_Newtonian_Matrix_Init )
    CALL Initialize_Poisson_Matrix()
    CALL TimerStop( Timer_Newtonian_Matrix_Init )

    CALL TimerStop( Timer_Newtonian_Initialization )

    lPF_Init_Flags(iPF_Init_Method_Vars) = .TRUE.

END IF

END SUBROUTINE Initialize_Poisson




#ifdef POSEIDON_AMREX_FLAG
 !+201+############################################################!
!                                                                   !
!              Initialization_Poisson_with_AMReX                    !
!                                                                   !
 !#################################################################!
SUBROUTINE  Initialization_Poisson_with_AMReX( )

! Determine Radial Base Variables from Multifab
IF ( .NOT. lPF_Init_Flags(iPF_Init_Method_Vars) ) THEN
    
    IF ( Verbose_Flag ) CALL Init_Message('Initializing Poisson system with AMReX.')
!    CALL TimerStart(Timer_Initialization_Poisson)

    
    CALL Initialize_AMReX_Maps()
    CALL Initialize_Level_Tables()
    CALL Initialize_Derived_AMReX_Part2()

    CALL Allocate_Mesh()
    CALL Determine_AMReX_Mesh()

    ! Allocate Arrays
    CALL Allocate_Poisson_Linear_System()

    ! Construct Matrices
    CALL TimerStart( Timer_Newtonian_Matrix_Init )
    CALL Initialize_Poisson_Matrix()
    CALL TimerStop( Timer_Newtonian_Matrix_Init )

!    CALL TimerStop(Timer_Initialization_Poisson)

    lPF_Init_Flags(iPF_Init_Method_Vars) = .TRUE.

    IF ( lPF_IO_Flags(iPF_IO_Print_Setup) ) THEN
        CALL PRINT_AMReX_Setup()
    END IF
END IF

END SUBROUTINE  Initialization_Poisson_with_AMReX


 !+101+############################################################!
!                                                                   !
!           Reinitialization_Poisson_with_AMReX                     !
!                                                                   !
 !#################################################################!
SUBROUTINE  Reinitialization_Poisson_with_AMReX( )


! Determine Radial Base Variables from Multifab
IF ( lPF_Init_Flags(iPF_Init_Method_Vars) ) THEN

    lPF_Init_Flags(iPF_Init_Method_Vars) = .FALSE.
    
    IF ( Verbose_Flag ) CALL Init_Message('Initializing XCFC system with AMReX.')
!    CALL TimerStart(Timer_Initialization_XCFC)

    CALL Make_Remesh_Copies()
    CALL Reinitialize_AMReX_Maps()
    CALL Initialize_Derived_AMReX_Part2()
    CALL Reallocate_Mesh()
    CALL Determine_AMReX_Mesh()
    CALL Reallocate_Poisson_Linear_System()

    
    
    CALL Fill_Coeff_Vector_From_Copy()
    CALL Destroy_Remesh_Copies

    ! Construct Matrices
    CALL TimerStart( Timer_Newtonian_Matrix_Init )
    CALL Initialize_Poisson_Matrix()
    CALL TimerStop( Timer_Newtonian_Matrix_Init )


!    CALL TimerStop(Timer_Initialization_XCFC)

    lPF_Init_Flags(iPF_Init_Method_Vars) = .TRUE.

    IF ( Verbose_Flag ) THEN
        CALL PRINT_AMReX_Setup()
    END IF

END IF

END SUBROUTINE  Reinitialization_Poisson_with_AMReX

#endif

END MODULE Initialization_Poisson


