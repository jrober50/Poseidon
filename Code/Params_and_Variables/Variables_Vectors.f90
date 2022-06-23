   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Variables_Vectors                                               	     !##!
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


COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:,:)        ::  cVA_Source_Vector
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:)          ::  cVB_Source_Vector
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:)          ::  cVP_Source_Vector

COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:,:)        ::  cVA_Coeff_Vector
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:)          ::  cVB_Coeff_Vector
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:)          ::  cVP_Coeff_Vector


CONTAINS

END MODULE Variables_Vectors
