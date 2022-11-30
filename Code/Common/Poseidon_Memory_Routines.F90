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




IMPLICIT NONE


CONTAINS


!+101+##########################################################################!
!                                                                               !
!        Poseidon_Memory_Usage                                    				!
!                                                                               !
!###############################################################################!
SUBROUTINE Poseidon_Memory_Usage( RSS_Mem )

INTEGER,    INTENT(OUT)         ::  RSS_Mem


CHARACTER(LEN = 500)            ::  Filename=''
CHARACTER(LEN = 100)            ::  Line
CHARACTER(LEN = 8)              ::  Pid_Str

INTEGER                         ::  pid

LOGICAL                         ::  ifxst


RSS_Mem=-1    ! return negative number if not found

!--- get process ID

pid=getpid()
WRITE(Pid_Str,'(I8)') pid
filename='/proc/'//trim(adjustl(Pid_Str))//'/status'

!--- read system file

inquire (file=filename,exist=ifxst)
if (.not.ifxst) then
  write (*,*) 'system file does not exist'
  return
endif

open(unit=100, file=filename, action='read')
do
  read (100,'(a)',end=120) line
  if (line(1:6).eq.'VmRSS:') then
     read (line(7:),*) RSS_Mem
     exit
  endif
enddo
120 continue
close(100)

RETURN


END SUBROUTINE Poseidon_Memory_Usage




END MODULE Poseidon_Memory_Routines
