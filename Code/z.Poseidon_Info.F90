   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_Info_Module                                                         !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
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



USE Poseidon_Constants_Module, &
                    ONLY : idp

USE Poseidon_Parameters, &
                    ONLY :  CUR_ITERATION,                                  &
                            MAX_ITERATIONS

USE Poseidon_Variables_Module, &
                    ONLY :  Iteration_Histogram,                            &
                            Frame_Time_Table,                               &
                            Run_Time_Table





IMPLICIT NONE


!*F&S*==========================================!
!                                               !
!           Functions & Subroutines             !
!                                               !
!===============================================!
CONTAINS



 !+101+############################################################################!
!                                                                                   !
!                     PQ_MAX_ITERATIONS                                             !
!                                                                                   !
 !#################################################################################!
SUBROUTINE PQ_ITERATIONS_MAX( Iters_Max )

INTEGER, INTENT(OUT)                    :: Iters_Max

Iters_Max = Max_Iterations


END SUBROUTINE PQ_ITERATIONS_MAX



 !+102+############################################################################!
!                                                                                   !
!                     PQ_ITERATIONS_USED                                            !
!                                                                                   !
 !#################################################################################!
SUBROUTINE PQ_ITERATIONS_USED( Iters_Used )

INTEGER, INTENT(OUT)                    :: Iters_Used

Iters_Used = Cur_Iteration - 1

END SUBROUTINE PQ_ITERATIONS_USED






 !+102+############################################################################!
!                                                                                   !
!                     PQ_ITERATIONS_USED                                            !
!                                                                                   !
 !#################################################################################!
SUBROUTINE PQ_ITERATIONS_HIST( Iters_Hist, Iters_Max )

INTEGER, INTENT(IN)                                         :: Iters_Max
INTEGER, DIMENSION(1:Iters_Max), INTENT(OUT)                :: Iters_Hist


IF ( Iters_Max .NE. Max_Iterations) THEN

    WRITE(*,'(A)')"PQ_ITERATIONS_HIST call failed. Storage Array is Innaproprate."
    WRITE(*,'(A,I3.3,A,I3.3)')"Expected Array Size: ",Max_Iterations,           &
                              "  Recieved Array Size: ",Iters_Max

ELSE
    
    Iters_Hist = Iteration_Histogram

END IF

END SUBROUTINE PQ_ITERATIONS_HIST




 !+102+############################################################################!
!                                                                                   !
!                     PQ_ITERATIONS_USED                                            !
!                                                                                   !
 !#################################################################################!
SUBROUTINE PQ_TIMETABLE_FRAME( Timetable )

REAL(KIND = idp), DIMENSION(1:25), INTENT(OUT)                    :: Timetable

Timetable = Frame_Time_table

END SUBROUTINE PQ_TIMETABLE_FRAME






 !+102+############################################################################!
!                                                                                   !
!                     PQ_ITERATIONS_USED                                            !
!                                                                                   !
 !#################################################################################!
SUBROUTINE PQ_TIMETABLE_RUN( Timetable )

REAL(KIND = idp), DIMENSION(1:25), INTENT(OUT)                    :: Timetable

Timetable = Run_Time_table

END SUBROUTINE PQ_TIMETABLE_RUN





END MODULE Poseidon_Info_Module
