   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Variables_Interface                                             	     !##!
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

IMPLICIT NONE

LOGICAL                                     ::  Caller_Set = .FALSE.

INTEGER,    DIMENSION(3)                    ::  Caller_NQ
REAL(idp),  DIMENSION(2)                    ::  Caller_xL
INTEGER                                     ::  Caller_Quad_DOF

REAL(idp),  DIMENSION(:),   ALLOCATABLE     ::  Caller_RQ_xlocs
REAL(idp),  DIMENSION(:),   ALLOCATABLE     ::  Caller_TQ_xlocs
REAL(idp),  DIMENSION(:),   ALLOCATABLE     ::  Caller_PQ_xlocs

REAL(idp)                                   ::  Caller_R_Units

REAL(idp),  DIMENSION(:,:), ALLOCATABLE     ::  Translation_Matrix


END MODULE Variables_Interface
