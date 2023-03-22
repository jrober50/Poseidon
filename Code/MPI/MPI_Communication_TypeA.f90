   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE MPI_Communication_TypeA_Module                                        !##!
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

USE Variables_Derived, &
            ONLY :  Num_R_Nodes,                &
                    LM_Length

USE Variables_Vectors,  &
            ONLY :  dVA_Coeff_Vector,            &
                    dVA_Load_Vector

USE Variables_MPI, &
            ONLY :  myID_Poseidon

USE MPI

IMPLICIT NONE


CONTAINS

!+101+##########################################################################!
!                                                                               !
!          MPI_RTM_Source_TypeA                                         		!
!                                                                               !
!###############################################################################!
SUBROUTINE MPI_RTM_Source_TypeA( iU, LLim, ULim, MasterID, COMM, ierr )

INTEGER, INTENT(IN)                     :: iU
INTEGER, INTENT(IN)                     :: LLim
INTEGER, INTENT(IN)                     :: ULim
INTEGER, INTENT(IN)                     :: MasterID
INTEGER, INTENT(IN)                     :: COMM
INTEGER, INTENT(INOUT)                  :: ierr

INTEGER                                 :: lm_loc
INTEGER                                 :: Send_Size

Send_Size = ULim - LLim + 1


DO lm_loc = 1,LM_Length
    IF ( myID_Poseidon == MasterID ) THEN


        CALL MPI_Reduce(MPI_IN_PLACE,                   &
                        dVA_Load_Vector(LLim:ULim,lm_loc,iU),&
                        Send_Size,                      &
                        MPI_Double,             &
                        MPI_SUM,                        &
                        MasterID,                       &
                        COMM,                           &
                        ierr )


    ELSE


        CALL MPI_Reduce(dVA_Load_Vector(LLim:ULim,lm_loc,iU),&
                        dVA_Load_Vector(LLim:ULim,lm_loc,iU),&
                        Send_Size,                              &
                        MPI_Double,                     &
                        MPI_SUM,                                &
                        MasterID,                               &
                        COMM,                                   &
                        ierr )


    END IF
END DO



END SUBROUTINE MPI_RTM_Source_TypeA






!+101+##########################################################################!
!                                                                               !
!          MPI_RTM_Source_TypeA                                                 !
!                                                                               !
!###############################################################################!
SUBROUTINE MPI_BCAST_Coeffs_TypeA( iU, LLim, ULim, MasterID, COMM, ierr )

INTEGER, INTENT(IN)                     :: iU
INTEGER, INTENT(IN)                     :: LLim
INTEGER, INTENT(IN)                     :: ULim
INTEGER, INTENT(IN)                     :: MasterID
INTEGER, INTENT(IN)                     :: COMM
INTEGER, INTENT(INOUT)                  :: ierr

INTEGER                                 :: lm_loc
INTEGER                                 :: Send_Size

Send_Size = ULim - LLim + 1


DO lm_loc = 1,LM_Length
    CALL MPI_Bcast( dVA_Coeff_Vector( LLim:ULim,lm_loc,iU),&
                    Send_Size,                              &
                    MPI_Double,                     &
                    MasterID,                               &
                    Comm,                                   &
                    ierr                                    )
END DO


END SUBROUTINE MPI_BCAST_Coeffs_TypeA








END MODULE MPI_Communication_TypeA_Module
