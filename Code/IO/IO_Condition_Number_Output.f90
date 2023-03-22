   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE IO_Condition_Number_Output_Module                                     !##!
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

USE Poseidon_IO_Parameters, &
            ONLY :  Poseidon_Reports_Dir

USE Poseidon_Message_Routines_Module, &
            ONLY :  Run_Message,                &
                    Warning_Message

USE Variables_Derived, &
            ONLY :  iVB_Prob_Dim


USE Variables_Matrices, &
            ONLY :  iMB_Diagonals,             &
                    iMB_IPIV,                  &
                    dMB_Matrix_Banded


USE Variables_IO, &
            ONLY :  File_Suffix

USE Variables_MPI, &
            ONLY :  myID_Poseidon,          &
                    MasterID_Poseidon


USE Poseidon_File_Routines_Module, &
            ONLY :  Open_New_File

USE Flags_IO_Module, &
            ONLY :  lPF_IO_Flags,               &
                    iPF_IO_Write_Cond,          &
                    iPF_IO_Print_Cond

IMPLICIT NONE


CONTAINS



 !+101+####################################################!
!                                                           !
!          IO_Output_Condition_Number	                    !
!                                                           !
 !#########################################################!
SUBROUTINE IO_Output_Condition_Number( )


INTEGER                                         ::  i
INTEGER                                         ::  Info
REAL(idp),      DIMENSION(2*iVB_Prob_Dim)       ::  Work
REAL(idp),      DIMENSION(iVB_Prob_Dim)         ::  RWork
REAL(idp)                                       ::  RCond_One
REAL(idp)                                       ::  Norm
CHARACTER(LEN = 300)                            ::  Message


IF ( lPF_IO_Flags(iPF_IO_Print_Cond) .OR. lPF_IO_Flags(iPF_IO_Write_Cond) ) THEN

    NORM = 0.0_idp
    DO i = 1,iVB_Prob_Dim

        NORM = MAX( NORM, ABS(SUM(dMB_Matrix_Banded(:,i) ) ) )

    END DO

    PRINT*,"Before DGBCON"
    CALL DGBCON( '1',                   &
                 iVB_Prob_Dim,          &
                 iMB_Diagonals,         &
                 iMB_Diagonals,         &
                 dMB_Matrix_Banded,     &
                 3*iMB_Diagonals+1,     &
                 iMB_IPIV,              &
                 Norm,                  &
                 RCond_One,             &
                 Work,                  &
                 RWork,                 &
                 Info                   )

    PRINT*,"After DGBCON"
    IF (Info .NE. 0) THEN
        WRITE(Message,'(A,I1.1)')"ZGBCON has failed with INFO = ",Info
        CALL Warning_Message(TRIM(Message))
    ELSE

        PRINT*,"Here"
        CALL IO_Print_Condition_Number( RCond_One )
        PRINT*,"There"
        CALL IO_Write_Condition_Number( RCond_One )
        PRINT*,"Everywhere"
    END IF

    
END IF

END SUBROUTINE IO_Output_Condition_Number





 !+101+####################################################!
!                                                           !
!          IO_Print_Condition_Number                        !
!                                                           !
 !#########################################################!
SUBROUTINE IO_Print_Condition_Number( rcond )

REAL(idp), INTENT(IN)               ::  rcond
CHARACTER(LEN = 300)                ::  Message


IF (myID_Poseidon == MasterID_Poseidon ) THEN
IF ( lPF_IO_Flags(iPF_IO_Print_Cond) ) THEN

    WRITE(Message,'(A,ES22.15,A)')'The reciprocal condition number is ',rcond,'.'
    CALL Run_Message(TRIM(Message))

END IF
END IF

END SUBROUTINE IO_Print_Condition_Number





 !+202+####################################################!
!                                                           !
!          IO_Write_Condition_Number                        !
!                                                           !
 !#########################################################!
SUBROUTINE IO_Write_Condition_Number( rcond )

REAL(idp), INTENT(IN)               ::  rcond

CHARACTER(LEN = 200)                ::  Report_Name
INTEGER                             ::  Report_ID
INTEGER                             ::  Suggested_Number

Suggested_Number = 500

IF (myID_Poseidon == MasterID_Poseidon ) THEN
IF ( lPF_IO_Flags(iPF_IO_Write_Cond) ) THEN

    WRITE(Report_Name,'(A,A,A,A)')              &
            Poseidon_Reports_Dir,"Condition_Number_",trim(File_Suffix),".out"

    CALL Open_New_File( Report_Name, Report_ID, Suggested_Number)

    WRITE(Report_ID,'(ES22.15)') 1.0_idp/rcond

    CLOSE(Report_ID)

END IF
END IF

END SUBROUTINE IO_Write_Condition_Number



END MODULE IO_Condition_Number_Output_Module
