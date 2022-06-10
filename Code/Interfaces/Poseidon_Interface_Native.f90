    !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_Main_Module                                                         !##!
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
            ONLY :  Poseidon_Frame,             &
                    Poisson_Mode,               &
                    Method_Flag

USE Variables_Interface, &
            ONLY :  Caller_Set,                     &
                    Caller_RQ_xlocs,                &
                    Caller_TQ_xlocs,                &
                    Caller_PQ_xlocs,                &
                    Translation_Matrix


USE Poisson_Main_Module, &
            ONLY :  Poisson_Solve

USE FP_AndersonM_Module, &
            ONLY : Fixed_Point_AndersonM

USE XCFC_Method_Module, &
            ONLY : XCFC_Method


USE Allocation_Sources, &
            ONLY :  Deallocate_Poseidon_Source_Variables

USE Allocation_Poisson, &
            ONLY :  Deallocate_Poseidon_Poisson_Variables


USE XCFC_Source_Routine_Variables_Module, &
            ONLY :  Deallocate_XCFC_Source_Routine_Variables

USE Allocation_Mesh, &
            ONLY : Deallocate_Mesh

USE Allocation_Quadrature, &
            ONLY :  Deallocate_Quadrature

USE Allocation_Tables, &
            ONLY :  Deallocate_Tables

USE Allocation_FP, &
            ONLY :  Deallocate_FP

USE Allocation_XCFC_Linear_Systems, &
            ONLY :  Deallocate_XCFC_Linear_Systems

USE IO_Print_Results, &
            ONLY :  Print_Results

USE IO_Write_Final_Results, &
            ONLY :  Write_Final_Results

USE Timer_Routines_Module, &
            ONLY :  Finalize_Timers

USE Flags_Main_Module, &
            ONLY : Poseidon_Clear_All_Flags

USE Flags_Check_Routines, &
            ONLY :  Poseidon_Run_Check

USE Flags_Boundary_Conditions_Module, &
            ONLY :  lPF_BC_Flags,           &
                    iPF_BC_Outer_Set,       &
                    iPF_BC_Inner_Set

USE Flags_IO_Module, &
            ONLY :  lPF_IO_Flags,           &
                    iPF_IO_Print_Results,   &
                    iPF_IO_Write_Results

IMPLICIT NONE




                    !*F&S*==========================================!
                    !                                               !
                    !           Functions & Subroutines             !
                    !                                               !
                    !===============================================!
CONTAINS





!+102+#######################################################################################!
!                                                                                           !
!       Poseidon_Run - Calculates Source Vector and Solves for solution coefficients        !
!                                                                                           !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Run()

LOGICAL                                             ::  Readiness_Flag


!Readiness_Flag = .TRUE.
Readiness_Flag = Poseidon_Run_Check()


IF ( Readiness_Flag ) THEN

    IF ( Poisson_Mode .eqv. .TRUE. ) THEN

        CALL Poisson_Solve()
    
    ELSE IF ( Method_Flag == 1 ) THEN
        WRITE(*,'(A)')"The Newton-Raphson method is not currently available in Poseidon. STOPING"
        STOP

!        CALL CFA_Newton_Raphson_3D()

    ELSE IF ( Method_Flag == 2 ) THEN

        Call Fixed_Point_AndersonM()

    ELSE IF ( Method_Flag == 3 ) THEN
        CALL XCFC_Method()
    ELSE

        PRINT*,"ERROR IN POSEIDON : Solver Type Flag has invalid value. "
        STOP

    END IF

ELSE

    PRINT*, "ERROR IN POSEIDON : There was an error in setting up Poseidon, therefore it did not run."
    STOP

END IF
Poseidon_Frame = Poseidon_Frame + 1



IF ( lPF_IO_Flags(iPF_IO_Print_Results) ) THEN
    Call Print_Results()
END IF
IF ( lPF_IO_Flags(iPF_IO_Write_Results) ) THEN
    Call Write_Final_Results()
END IF


END SUBROUTINE Poseidon_Run
















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


IF ( Poisson_Mode ) THEN
    Call Deallocate_Poseidon_Poisson_Variables

ELSE
    !!!!  Deallocate Data Space !!!!
    CALL Deallocate_Poseidon_Source_Variables
    Call Deallocate_XCFC_Source_Routine_Variables

    CALL Deallocate_Quadrature()
    CALL Deallocate_Tables()


    IF ( Method_Flag == 1 ) THEN
        WRITE(*,'(A)')"The Newton-Raphson method is not currently available in Poseidon. STOPING"
        STOP
!        CALL Deallocate_NR
    ELSE IF ( Method_Flag == 2 ) THEN
        CALL Deallocate_FP
    ELSE IF ( Method_Flag == 3 ) THEN
        CALL Deallocate_XCFC_Linear_Systems
    END IF

END IF


CALL Finalize_Timers()
CALL Poseidon_Clear_All_Flags()


END SUBROUTINE Poseidon_Close








END MODULE Poseidon_Main_Module
