   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Flags_Main_Module                                                     !##!
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
            ONLY : idp

USE Flags_Run_Check_Module, &
            ONLY :  lPF_RC_Flags


USE Flags_Initialization_Module, &
            ONLY :  lPF_Init_Flags,                 &
                    lPF_Init_Matrices_Flags,        &
                    lPF_Init_Quad_Flags,            &
                    lPF_Init_MTGV_Flags

USE Flags_Allocation_Module, &
            ONLY :  lPF_Alloc_Flags

USE Flags_Boundary_Conditions_Module, &
            ONLY :  lPF_BC_Flags

USE Flags_Initial_Guess_Module, &
            ONLY :  lPF_IG_Flags

USE Flags_Source_Input_Module, &
            ONLY :  lPF_SI_Flags

#ifdef POSEIDON_DEBUG_FLAG
USE Flags_Initial_Guess_Module, &
            ONLY :  lPF_Debug_Flags
#endif

USE Flag_Routines_Module, &
            ONLY :  Poseidon_Clear_Flag

IMPLICIT NONE


CONTAINS

!+101+####################################################!
!                                                           !
!          Poseidon_Clear_All_Flags                            !
!                                                           !
!#########################################################!
SUBROUTINE Poseidon_Clear_All_Flags()

CALL Poseidon_Clear_Flag( lPF_RC_Flags )

CALL Poseidon_Clear_Flag( lPF_Init_Flags )
CALL Poseidon_Clear_Flag( lPF_Init_Matrices_Flags )
CALL Poseidon_Clear_Flag( lPF_Init_Quad_Flags )
CALL Poseidon_Clear_Flag( lPF_Init_MTGV_Flags )

CALL Poseidon_Clear_Flag( lPF_Alloc_Flags )

CALL Poseidon_Clear_Flag( lPF_BC_Flags )

CALL Poseidon_Clear_Flag( lPF_IG_Flags )

CALL Poseidon_Clear_Flag( lPF_SI_Flags )


#ifdef POSEIDON_DEBUG_FLAG
CALL Poseidon_Clear_Flag( lPF_Debug_Flags )
#endif

END SUBROUTINE Poseidon_Clear_All_Flags








END MODULE Flags_Main_Module
