   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Poseidon_Memory_Routines                                              !##!
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
            ONLY :  idp

USE Memory_Variables_Module


IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!        Poseidon_Memory_Usage                                    				!
!                                                                               !
!###############################################################################!
SUBROUTINE Poseidon_Mark_Memory( RSS_Memory_Marker, HWM_Memory_Marker )

INTEGER,    INTENT(OUT)         ::  RSS_Memory_Marker
INTEGER,    INTENT(OUT), OPTIONAL :: HWM_Memory_Marker

CHARACTER(LEN = 500)            ::  Filename=''
CHARACTER(LEN = 100)            ::  Line
CHARACTER(LEN = 8)              ::  Pid_Str

INTEGER                         ::  pid

LOGICAL                         ::  ifxst

INTEGER                         ::  HWM

RSS_Memory_Marker = -1    ! return negative number if not found

!--- get process ID

pid=getpid()
WRITE(Pid_Str,'(I8)') pid
filename='/proc/'//trim(adjustl(Pid_Str))//'/status'

!--- read system file

INQUIRE (file=filename,exist=ifxst)
IF (.NOT. ifxst) THEN
  WRITE(*,*) 'system file does not exist'
  RETURN
ENDIF

OPEN(unit=100, file=filename, action='read')
DO
  READ (100,'(a)',end=120) line
  IF (line(1:6).eq.'VmRSS:') then
     READ (line(7:),*) RSS_Memory_Marker
  ENDIF
  IF (line(1:6).eq.'VmHWM:') then
     READ (line(7:),*) HWM
  ENDIF
ENDDO
120 CONTINUE
CLOSE(100)

IF (PRESENT( HWM_Memory_Marker )) THEN
        HWM_Memory_Marker = HWM
ENDIF


END SUBROUTINE Poseidon_Mark_Memory







END MODULE Poseidon_Memory_Routines
