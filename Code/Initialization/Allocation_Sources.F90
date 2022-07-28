   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Allocation_Sources                                                           !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   Allocate_Poseidon_Source_Variables                                  !##!
!##!    +102+   Deallocate_Poseidon_Source_Variables                                !##!
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
            ONLY : idp

USE Poseidon_Message_Routines_Module, &
            ONLY :  Init_Message

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag

USE Variables_Quadrature, &
            ONLY :  Local_Quad_DOF

USE Variables_Mesh, &
            ONLY :  Num_R_Elements,     &
                    Num_T_Elements,     &
                    Num_P_Elements

USE Variables_Source, &
            ONLY :  Source_Rho,             &
                    Block_Source_E,         &
                    Block_Source_S,         &
                    Block_Source_Si

USE Flags_Initialization_Module, &
            ONLY :  lPF_Init_Flags,         &
                    iPF_Init_Alloc_Source

USE Flags_Core_Module, &
            ONLY :  iPF_Core_Flags,         &
                    iPF_Core_Method_Mode,   &
                    iPF_Core_Method_Newtonian


#ifdef POSEIDON_AMREX_FLAG
USE Variables_AMReX_Core, &
            ONLY :  MF_Source,              &
                    BA_Source,              &
                    DM_Source,              &
                    GM_Source,              &
                    iLeafElementsPerLvl,    &
                    AMReX_Num_Levels,       &
                    AMReX_Max_Level
#endif

IMPLICIT NONE


CONTAINS



!+101+###########################################################################!
!                                                                                !
!                            Allocate_Poseidon_Variables                         !
!                                                                                !
!################################################################################!
SUBROUTINE Allocate_Poseidon_Source_Variables()

IF ( Verbose_Flag ) CALL Init_Message('Allocating Source Variables.')

#ifdef POSEIDON_AMREX_FLAG

ALLOCATE( MF_Source(0:AMReX_Max_Level))
ALLOCATE( BA_Source(0:AMReX_Max_Level))
ALLOCATE( DM_Source(0:AMReX_Max_Level))
ALLOCATE( GM_Source(0:AMReX_Max_Level))

ALLOCATE( iLeafElementsPerLvl(0:AMReX_Num_Levels-1))

#else

IF ( iPF_Core_Flags(iPF_Core_Method_Mode) == iPF_Core_Method_Newtonian ) THEN

    ALLOCATE(Source_Rho(    1:Local_Quad_DOF,       &
                            0:NUM_R_ELEMENTS-1,     &
                            0:NUM_T_ELEMENTS-1,     &
                            0:NUM_P_ELEMENTS-1  )   )

ELSE

    ALLOCATE(Block_Source_E(    1:Local_Quad_DOF,       &
                                0:NUM_R_ELEMENTS-1,     &
                                0:NUM_T_ELEMENTS-1,     &
                                0:NUM_P_ELEMENTS-1  )   )

    ALLOCATE(Block_Source_S(    1:Local_Quad_DOF,       &
                                0:NUM_R_ELEMENTS-1,     &
                                0:NUM_T_ELEMENTS-1,     &
                                0:NUM_P_ELEMENTS-1  )   )

    ALLOCATE(Block_Source_Si(   1:Local_Quad_DOF,       &
                                0:NUM_R_ELEMENTS-1,     &
                                0:NUM_T_ELEMENTS-1,     &
                                0:NUM_P_ELEMENTS-1,     &
                                1:3                 )   )

END IF
#endif


lPF_Init_Flags(iPF_Init_Alloc_Source) = .TRUE.


END SUBROUTINE Allocate_Poseidon_Source_Variables











!+102+###########################################################################!
!                                                                                !
!                           Deallocate_Poseidon_Variables                        !
!                                                                                !
!################################################################################!
SUBROUTINE Deallocate_Poseidon_Source_Variables()

#ifdef POSEIDON_AMREX_FLAG


DEALLOCATE( MF_Source )
DEALLOCATE( BA_Source )
DEALLOCATE( DM_Source )
DEALLOCATE( GM_Source )

DEALLOCATE( iLeafElementsPerLvl )

#else
IF ( iPF_Core_Flags(iPF_Core_Method_Mode) == iPF_Core_Method_Newtonian ) THEN

    DEALLOCATE( Source_Rho )

ELSE

    DEALLOCATE( Block_Source_E )
    DEALLOCATE( Block_Source_S )
    DEALLOCATE( Block_Source_Si )

END IF

#endif




END SUBROUTINE Deallocate_Poseidon_Source_Variables





END MODULE Allocation_Sources
