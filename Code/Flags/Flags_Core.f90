   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Flags_Core_Module                                                     !##!
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


INTEGER,    PUBLIC, PARAMETER       ::  iPF_Core_Num_Flags      = 3

INTEGER,    PUBLIC, PARAMETER       ::  iPF_Core_Unit_Mode      = 1
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Core_Method_Mode    = 2
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Core_AMReX_Mode     = 3

INTEGER,    PUBLIC, DIMENSION(1:iPF_Core_Num_Flags)    ::  iPF_Core_Flags





INTEGER,    PUBLIC, PARAMETER       ::  iPF_Core_Num_Unit_Modes   = 4

CHARACTER(LEN=20), PUBLIC, DIMENSION(1:iPF_Core_Num_Unit_Modes)   ::  &
                Poseidon_Unit_Mode_Names = (/'CGS            ',         &
                                             'MKS            ',         &
                                             'Geometrized    ',         &
                                             'Unitless       '          /)

INTEGER,    PUBLIC, PARAMETER       ::  iPF_Core_Units_CGS          = 1
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Core_Units_MKS          = 2
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Core_Units_Geometrized  = 3
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Core_Units_Unitless     = 4








INTEGER,    PUBLIC, PARAMETER       ::  iPF_Core_Num_Method_Modes   = 3

CHARACTER(LEN=20), PUBLIC, DIMENSION(1:iPF_Core_Num_Method_Modes)   ::      &
                Poseidon_Method_Names = (/  'Newtonian Potential ',         &
                                            'CFA Metric          ',         &
                                            'XCFC Metric         '          /)

INTEGER,    PUBLIC, PARAMETER       ::  iPF_Core_Method_Newtonian = 1
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Core_Method_CFA       = 2
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Core_Method_XCFC      = 3







INTEGER,    PUBLIC, PARAMETER       ::  iPF_Core_Num_AMReX_Modes   = 2

CHARACTER(LEN=3), PUBLIC, DIMENSION(1:iPF_Core_Num_AMReX_Modes)   ::       &
                Poseidon_AMReX_Mode = (/ 'Off', 'ON ' /)

INTEGER,    PUBLIC, PARAMETER       ::  iPF_Core_AMReX_Off       = 1
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Core_AMReX_On        = 2


CONTAINS







END MODULE Flags_Core_Module

