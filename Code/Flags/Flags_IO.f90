   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Flags_IO_Module                                                 	     !##!
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


IMPLICIT NONE

INTEGER,    PUBLIC, PARAMETER       ::  iPF_IO_Num_Flags            = 18

INTEGER,    PUBLIC, PARAMETER       ::  iPF_IO_Verbose              = 1
INTEGER,    PUBLIC, PARAMETER       ::  iPF_IO_Print_Setup          = 2
INTEGER,    PUBLIC, PARAMETER       ::  iPF_IO_Write_Setup          = 3
INTEGER,    PUBLIC, PARAMETER       ::  iPF_IO_Print_Results        = 4
INTEGER,    PUBLIC, PARAMETER       ::  iPF_IO_Write_Results        = 5
INTEGER,    PUBLIC, PARAMETER       ::  iPF_IO_Print_TimeTable      = 6
INTEGER,    PUBLIC, PARAMETER       ::  iPF_IO_Write_TimeTable      = 7
INTEGER,    PUBLIC, PARAMETER       ::  iPF_IO_Print_Iter_Report    = 8
INTEGER,    PUBLIC, PARAMETER       ::  iPF_IO_Write_Iter_Report    = 9
INTEGER,    PUBLIC, PARAMETER       ::  iPF_IO_Print_Frame_Report   = 10
INTEGER,    PUBLIC, PARAMETER       ::  iPF_IO_Write_Frame_Report   = 11
INTEGER,    PUBLIC, PARAMETER       ::  iPF_IO_Print_Run_Report     = 12
INTEGER,    PUBLIC, PARAMETER       ::  iPF_IO_Write_Run_Report     = 13
INTEGER,    PUBLIC, PARAMETER       ::  iPF_IO_Write_Sources        = 14
INTEGER,    PUBLIC, PARAMETER       ::  iPF_IO_Write_Convergence    = 15
INTEGER,    PUBLIC, PARAMETER       ::  iPF_IO_Write_Mesh           = 16
INTEGER,    PUBLIC, PARAMETER       ::  iPF_IO_Print_Cond           = 17
INTEGER,    PUBLIC, PARAMETER       ::  iPF_IO_Write_Cond           = 18


LOGICAL,    PUBLIC, DIMENSION(1:iPF_IO_Num_Flags)    ::  lPF_IO_Flags

CONTAINS






END MODULE Flags_IO_Module
