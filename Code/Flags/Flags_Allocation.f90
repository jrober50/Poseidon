   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Flags_Allocation_Module                                                      !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains parameters used to define the running of Poseidon.                 !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!

USE Poseidon_Kinds_Module, &
                ONLY : idp

IMPLICIT NONE


! Allocation Flags
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Alloc_Num_Flags         = 3

INTEGER,    PUBLIC, PARAMETER       ::  iPF_Alloc_Variables         = 1
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Alloc_XCFC_Variables    = 2
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Alloc_XCFC_Source_Vars  = 3

LOGICAL,    PUBLIC, DIMENSION(1:iPF_Alloc_Num_Flags)    ::  lPF_Alloc_Flags

LOGICAL, PUBLIC         ::  Poseidon_Variables_Allocated_Flag = .FALSE.





CONTAINS








END MODULE Flags_Allocation_Module



