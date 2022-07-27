   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Poseidon_Interface_Run                                          	     !##!
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
USE Poseidon_Parameters, &
            ONLY :  Poseidon_Frame,             &
                    Method_Flag

USE Poisson_Main_Module, &
            ONLY :  Poisson_Solve

USE FP_AndersonM_Module, &
            ONLY :  Fixed_Point_AndersonM

USE XCFC_Method_Module, &
            ONLY :  XCFC_Method,                &
                    XCFC_Method_Part1,          &
                    XCFC_Method_Part2

USE Variables_MPI, &
            ONLY :  myID_Poseidon,              &
                    MasterID_Poseidon,          &
                    nPROCS_Poseidon,            &
                    Poseidon_Comm_World

USE IO_Print_Results, &
            ONLY :  Print_Results

USE IO_Write_Final_Results, &
            ONLY :  Write_Final_Results

USE Poseidon_Interface_Initial_Guess,       &
            ONLY :  Poseidon_Initialize_Flat_Guess

USE Flags_Check_Routines, &
            ONLY :  Poseidon_Run_Check

USE Flags_IO_Module, &
            ONLY :  lPF_IO_Flags,           &
                    iPF_IO_Print_Results,   &
                    iPF_IO_Write_Results

USE Flags_Initial_Guess_Module, &
            ONLY :  lPF_IG_Flags,           &
                    iPF_IG_Set,             &
                    iPF_IG_Flat_Guess

USE Flags_Core_Module, &
            ONLY :  lPF_Core_Flags,         &
                    iPF_Core_Newtonian_Mode

IMPLICIT NONE



CONTAINS

!+101+#######################################################################################!
!                                                                                           !
!       Poseidon_Run - Calculates Source Vector and Solves for solution coefficients        !
!                                                                                           !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Run()

LOGICAL                                             ::  Readiness_Flag


IF (( .NOT. lPF_IG_Flags(iPF_IG_Set)) .AND. lPF_IG_Flags(iPF_IG_Flat_Guess) ) THEN
    CALL Poseidon_Initialize_Flat_Guess
END IF



!Readiness_Flag = .TRUE.
Readiness_Flag = Poseidon_Run_Check()


IF ( Readiness_Flag ) THEN

    IF ( lPF_Core_Flags(iPF_Core_Newtonian_Mode) .eqv. .TRUE. ) THEN

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



 !+201+######################################################################################!
 !                                                                                           !
 !       Poseidon_XCFC_Run_Part1 - Performs the steps to calculate the XCFC Conformal        !
 !                                 Factor.  Upon completion of this call, the routine        !
 !                                 Poseidon_Return_ConFactor() can be used to acquire        !
 !                                 values of the conformal factor at desired locations.      !
 !                                                                                           !
 !###########################################################################################!
 SUBROUTINE Poseidon_XCFC_Run_Part1()

 INTEGER         :: i, ierr

IF (( .NOT. lPF_IG_Flags(iPF_IG_Set)) .AND. lPF_IG_Flags(iPF_IG_Flat_Guess) ) THEN
    CALL Poseidon_Initialize_Flat_Guess
END IF


 CALL XCFC_Method_Part1()



 IF ( lPF_IO_Flags(iPF_IO_Print_Results) ) THEN
 DO i = 0,nPROCs_Poseidon
     IF (myID_Poseidon == MasterID_Poseidon) THEN
         PRINT*,"myID_Poseidon :",myID_Poseidon
         Call Print_Results()
     END IF
     CALL MPI_Barrier(Poseidon_Comm_World,ierr)
 END DO
 END IF




 IF ( lPF_IO_Flags(iPF_IO_Write_Results) ) THEN
 IF ( myID_Poseidon == MasterID_Poseidon ) THEN
     Call Write_Final_Results()
 END IF
 END IF


 END SUBROUTINE Poseidon_XCFC_Run_Part1





 !+202+######################################################################################!
 !                                                                                           !
 !       Poseidon_XCFC_Run_Part2 - Performs the steps to calculate the XCFC Lapse Function,  !
 !                                 and Shift Vector.  Upon completion of this call, the      !
 !                                 routine, Poseidon_Return_Lapse() and                      !
 !                                 Poseidon_Return_Shift() can be used to acquire            !
 !                                 values of the those variables at desired locations.       !
 !                                                                                           !
 !###########################################################################################!
 SUBROUTINE Poseidon_XCFC_Run_Part2()

 INTEGER         :: i, ierr


IF (( .NOT. lPF_IG_Flags(iPF_IG_Set)) .AND. lPF_IG_Flags(iPF_IG_Flat_Guess) ) THEN
    CALL Poseidon_Initialize_Flat_Guess
END IF

 CALL XCFC_Method_Part2()



 IF ( lPF_IO_Flags(iPF_IO_Print_Results) ) THEN
 DO i = 0,nPROCs_Poseidon
     IF (myID_Poseidon == MasterID_Poseidon) THEN
         PRINT*,"myID_Poseidon :",myID_Poseidon
         Call Print_Results()
     END IF
     CALL MPI_Barrier(Poseidon_Comm_World,ierr)
 END DO
 END IF



 IF ( lPF_IO_Flags(iPF_IO_Write_Results) ) THEN
 IF ( myID_Poseidon == MasterID_Poseidon ) THEN
     Call Write_Final_Results()
 END IF
 END IF


 END SUBROUTINE Poseidon_XCFC_Run_Part2


END MODULE Poseidon_Interface_Run 
