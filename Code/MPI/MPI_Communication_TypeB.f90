   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE MPI_Communication_TypeB_Module                                        !##!
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

USE Poseidon_Message_Routines_Module, &
            ONLY :  Warning_Message

USE Variables_Derived, &
            ONLY :  Num_R_Nodes,                &
                    LM_Length

USE Variables_FP,  &
            ONLY :  FP_Coeff_Vector_B,            &
                    FP_Source_Vector_B

USE Variables_MPI, &
            ONLY :  myID_Poseidon,              &
                    nPROCS_Poseidon

USE MPI

IMPLICIT NONE


CONTAINS

!+101+##########################################################################!
!                                                                               !
!          MPI_RTM_Source_TypeB                                                 !
!                                                                               !
!###############################################################################!
SUBROUTINE MPI_RTM_Source_TypeB( iVB, LLim, ULim, MasterID, COMM, ierr )

INTEGER, INTENT(IN)                     :: iVB
INTEGER, INTENT(IN)                     :: LLim
INTEGER, INTENT(IN)                     :: ULim
INTEGER, INTENT(IN)                     :: MasterID
INTEGER, INTENT(IN)                     :: COMM
INTEGER, INTENT(INOUT)                  :: ierr

INTEGER                                 :: lm_loc
INTEGER                                 :: Send_Size

Send_Size = ULim - LLim + 1


IF ( myID_Poseidon == MasterID ) THEN
    CALL MPI_Reduce(MPI_IN_PLACE,                   &
                    FP_Source_Vector_B(LLim:ULim,iVB),      &
                    Send_Size,                  &
                    MPI_Double_Complex,             &
                    MPI_SUM,                        &
                    MasterID,              &
                    Comm,            &
                    ierr )
ELSE
CALL MPI_Reduce(FP_Source_Vector_B(LLim:ULim,iVB),      &
                FP_Source_Vector_B(LLim:ULim,iVB),      &
                Send_Size,                  &
                MPI_Double_Complex,             &
                MPI_SUM,                        &
                MasterID,              &
                Comm,            &
                ierr )


END IF




END SUBROUTINE MPI_RTM_Source_TypeB






!+101+##########################################################################!
!                                                                               !
!          MPI_RTM_Source_TypeB                                                 !
!                                                                               !
!###############################################################################!
SUBROUTINE MPI_BCAST_Coeffs_TypeB( iVB, LLim, ULim, MasterID, COMM, ierr )

INTEGER, INTENT(IN)                     :: iVB
INTEGER, INTENT(IN)                     :: LLim
INTEGER, INTENT(IN)                     :: ULim
INTEGER, INTENT(IN)                     :: MasterID
INTEGER, INTENT(IN)                     :: COMM
INTEGER, INTENT(INOUT)                  :: ierr

INTEGER                                 :: lm_loc
INTEGER                                 :: Send_Size


INTEGER                                 ::  i
CHARACTER(LEN = 300)                    ::  Message

Send_Size = ULim - LLim + 1



CALL MPI_Bcast( FP_Coeff_Vector_B(LLim:ULim,iVB),   &
                Send_Size,                          &
                MPI_Double_Complex,                 &
                MasterID,                           &
                Comm,                               &
                ierr                                )
IF (ierr .NE. 0) THEN
    WRITE(Message,'(A,I1.1)')"MPI_BCAST has failed with ierr = ",ierr
    CALL Warning_Message(TRIM(Message))
END IF




END SUBROUTINE MPI_BCAST_Coeffs_TypeB








END MODULE MPI_Communication_TypeB_Module

