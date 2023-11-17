    !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_Interface_Close                                                     !##!
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
!##!    +102+   Poseidon_Run                                                        !##!
!##!    +103+   Poseidon_Close                                                      !##!
!##!    +104+   Poseidon_Set_Mesh                                                   !##!
!##!    +105+   Poseidon_CFA_Set_Boundary_Condtions                                 !##!
!##!                                                                                !##!
!##!    +201+   Poseidon_Readiness_Check                                            !##!
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
USE Poseidon_Parameters, &
            ONLY :  Method_Flag

USE Variables_Interface, &
            ONLY :  Caller_Set,                     &
                    Caller_RQ_xlocs,                &
                    Caller_TQ_xlocs,                &
                    Caller_PQ_xlocs,                &
                    Translation_Matrix

USE Allocation_Sources, &
            ONLY :  Deallocate_Poseidon_Source_Variables

USE Load_Vector_Allocation_Module, &
            ONLY :  Deallocate_Load_Vector_Construction_Variables

USE Allocation_Mesh, &
            ONLY : Deallocate_Mesh

USE Allocation_Quadrature, &
            ONLY :  Deallocate_Quadrature

USE Allocation_Tables, &
            ONLY :  Deallocate_Tables

USE Allocation_Poisson_Linear_System, &
            ONLY :  Deallocate_Poisson_Linear_System

USE Allocation_CFA_Linear_Systems, &
            ONLY :  Deallocate_CFA_Linear_Systems

USE Allocation_XCFC_Linear_Systems, &
            ONLY :  Deallocate_XCFC_Linear_Systems

USE Allocation_FEM, &
            ONLY :  Deallocate_FEM

USE Timer_Routines_Module, &
            ONLY :  Finalize_Timers
            
USE Allocation_Maps, &
            ONLY :  Deallocate_AMReX_Maps

USE Flags_Main_Module, &
            ONLY : Poseidon_Clear_All_Flags

USE Flags_Core_Module, &
            ONLY :  iPF_Core_Flags,         &
                    iPF_Core_Method_Mode,   &
                    iPF_Core_Method_Newtonian


IMPLICIT NONE




                    !*F&S*==========================================!
                    !                                               !
                    !           Functions & Subroutines             !
                    !                                               !
                    !===============================================!
CONTAINS



!+103+######################################################################################!
!                                                                                           !
!       Poseidon_Close - Deallocate all Poseidon variables.                                 !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Close()


CALL Deallocate_Mesh()

IF ( Caller_Set ) THEN
    DEALLOCATE( Caller_RQ_xlocs )
    DEALLOCATE( Caller_TQ_xlocs )
    DEALLOCATE( Caller_PQ_xlocs )
    DEALLOCATE( Translation_Matrix )
    Caller_Set = .FALSE.
END IF


IF ( iPF_Core_Flags(iPF_Core_Method_Mode) == iPF_Core_Method_Newtonian) THEN
    CALL Deallocate_Poseidon_Source_Variables
    CALL Deallocate_Poisson_Linear_System

    CALL Deallocate_Quadrature()
    CALL Deallocate_Tables()
    CALL Deallocate_FEM()

ELSE
    !!!!  Deallocate Data Space !!!!
    CALL Deallocate_Poseidon_Source_Variables
    Call Deallocate_Load_Vector_Construction_Variables

    CALL Deallocate_Quadrature()
    CALL Deallocate_Tables()
    CALL Deallocate_FEM()

    IF ( Method_Flag == 1 ) THEN
        WRITE(*,'(A)')"The Newton-Raphson method is not currently available in Poseidon. STOPING"
        STOP
        
    ELSE IF ( Method_Flag == 2 ) THEN
        CALL Deallocate_CFA_Linear_Systems
    ELSE IF ( Method_Flag == 3 ) THEN
        CALL Deallocate_XCFC_Linear_Systems
    END IF

END IF

#ifdef POSEIDON_AMREX_FLAG
CALL Deallocate_AMReX_Maps
#endif

CALL Finalize_Timers()
CALL Poseidon_Clear_All_Flags()


END SUBROUTINE Poseidon_Close








END MODULE Poseidon_Interface_Close

