   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_File_Routines_Module                                                !##!
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






IMPLICIT NONE

CONTAINS





 !+501+############################################################################!
!                                                                                   !
!                     OPEN_NEW_FILE                                                 !
!                                                                                   !
 !#################################################################################!
SUBROUTINE OPEN_NEW_FILE(File_Name, File_Number, Suggested_Number)



CHARACTER(LEN = *), INTENT(IN)                          ::  File_Name
INTEGER,            INTENT(INOUT)                       ::  File_Number
INTEGER, OPTIONAL,  INTENT(IN)                          ::  Suggested_Number

INTEGER                                                 ::  Temp_Number
INTEGER                                                 ::  istat = 0
LOGICAL                                                 ::  FLAG, OP
LOGICAL                                                 ::  UNIT_FLAG, NAME_FLAG


UNIT_FLAG = .FALSE.
NAME_FLAG = .FALSE.


!  Assigned an unused number, and assign it to new file
FLAG = .TRUE.
IF ( Present(Suggested_Number) ) THEN
    Temp_Number = Suggested_Number
ELSE
    Temp_Number = 3000
END IF


DO WHILE (FLAG)
    INQUIRE( UNIT = Temp_Number, OPENED = OP )

    IF ( OP ) THEN
        Temp_Number = Temp_Number + 1
    ELSE
        File_Number = Temp_Number
        FLAG = .FALSE.
        UNIT_FLAG = .TRUE.
    END IF
END DO




! Open New File
IF ( UNIT_FLAG  ) THEN

    OPEN( Unit = File_Number, File = File_Name, IOSTAT = istat )
    IF ( istat .NE. 0 ) THEN

        PRINT*,"WARNING: Could not open file at ", File_Name, istat

    END IF
END IF


END SUBROUTINE OPEN_NEW_FILE




 !+501+############################################################################!
!                                                                                   !
!                     OPEN_NEW_FILE                                                 !
!                                                                                   !
 !#################################################################################!
SUBROUTINE Open_Existing_File_Append(File_Name, File_Number, Suggested_Number)



CHARACTER(LEN = *), INTENT(IN)                          ::  File_Name
INTEGER,            INTENT(INOUT)                       ::  File_Number
INTEGER, OPTIONAL,  INTENT(IN)                          ::  Suggested_Number

INTEGER                                                 ::  Temp_Number
INTEGER                                                 ::  istat = 0
LOGICAL                                                 ::  FLAG, OP
LOGICAL                                                 ::  UNIT_FLAG, NAME_FLAG


UNIT_FLAG = .FALSE.
NAME_FLAG = .FALSE.


!  Assigned an unused number, and assign it to new file
FLAG = .TRUE.
IF ( Present(Suggested_Number) ) THEN
    Temp_Number = Suggested_Number
ELSE
    Temp_Number = 3000
END IF


DO WHILE (FLAG)
    INQUIRE( UNIT = Temp_Number, OPENED = OP )

    IF ( OP ) THEN
        Temp_Number = Temp_Number + 1
    ELSE
        File_Number = Temp_Number
        FLAG = .FALSE.
        UNIT_FLAG = .TRUE.
    END IF
END DO




! Open New File
IF ( UNIT_FLAG  ) THEN

    OPEN( Unit = File_Number, File = File_Name, IOSTAT = istat, POSITION='APPEND' )
    IF ( istat .NE. 0 ) THEN

        PRINT*,"WARNING: Could not open file at ", File_Name, istat

    END IF
END IF




END SUBROUTINE OPEN_Existing_File_Append

 !+501+############################################################################!
!                                                                                   !
!                     OPEN_EXISTING_FILE                                            !
!                                                                                   !
 !#################################################################################!
SUBROUTINE OPEN_EXISTING_FILE_Rewind(File_Name, File_Number, istat)



CHARACTER(LEN = *), INTENT(IN)                          ::  File_Name
INTEGER,            INTENT(INOUT)                       ::  File_Number
INTEGER,            INTENT(INOUT)                       ::  istat


INTEGER                                                 ::  Temp_Number

LOGICAL                                                 ::  FLAG, OP
LOGICAL                                                 ::  UNIT_FLAG, NAME_FLAG


UNIT_FLAG = .FALSE.
NAME_FLAG = .FALSE.


!  Assigned an unused number, and assign it to new file
FLAG = .TRUE.
Temp_Number = 421
DO WHILE (FLAG)
    INQUIRE( UNIT = Temp_Number, OPENED = OP )

    IF ( OP ) THEN
        Temp_Number = Temp_Number + 1
    ELSE
        File_Number = Temp_Number
        FLAG = .FALSE.
        UNIT_FLAG = .TRUE.
    END IF
END DO




! Check if file already exists !
!INQUIRE( FILE = File_Name, EXIST = EX )
!IF ( EX ) THEN
!    PRINT*,"File ",File_Name," is already opened"
!
!ELSE
!    PRINT*,File_Name," is not already opened."
!    NAME_FLAG = .TRUE.
!END IF



! Open New File
IF ( UNIT_FLAG  ) THEN

    OPEN(unit=File_Number, file=File_Name, status='old', &
         iostat=istat, action='read', position='rewind')
    IF ( istat .NE. 0 ) THEN

        PRINT*,"WARNING: Could not open file at ", File_Name

    END IF
END IF




END SUBROUTINE OPEN_EXISTING_FILE_Rewind




 !+201+############################################################################!
!                                                                                   !
!                     OPEN_FILE_INQUISITION                                         !
!                                                                                   !
 !#################################################################################!
SUBROUTINE OPEN_FILE_INQUISITION()


INTEGER                 :: i
INTEGER                 :: i_min
INTEGER                 :: i_max
LOGICAL                 :: OP
CHARACTER(len=150)        :: name_of_file

i_min = 1
i_max = 3000

DO i = i_min,i_max
    OP = .FALSE.
    INQUIRE( UNIT = i, OPENED = OP, name=name_of_file )
    IF ( OP .EQV. .TRUE. ) THEN
        PRINT*,"File ",TRIM(name_of_file)," with Unit ",i," is still open."
    END IF
END DO


END SUBROUTINE OPEN_FILE_INQUISITION




END MODULE Poseidon_File_Routines_Module
