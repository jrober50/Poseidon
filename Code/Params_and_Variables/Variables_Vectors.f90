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


COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:,:)        ::  cVA_Load_Vector
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:)          ::  cVB_Load_Vector

COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:,:)        ::  cVA_Coeff_Vector
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:)          ::  cVB_Coeff_Vector


REAL(idp),   ALLOCATABLE,    DIMENSION(:,:,:)           ::  dVA_Load_Vector
REAL(idp),   ALLOCATABLE,    DIMENSION(:,:)             ::  dVB_Load_Vector

REAL(idp),   ALLOCATABLE,    DIMENSION(:,:,:)           ::  dVA_Coeff_Vector
REAL(idp),   ALLOCATABLE,    DIMENSION(:,:)             ::  dVB_Coeff_Vector

CONTAINS

END MODULE Variables_Vectors
