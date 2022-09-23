   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Variables_Matrices                                             	     !##!
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


CHARACTER(LEN = 4)                                      ::  Matrix_Format = 'CCS'
CHARACTER(LEN = 4)                                      ::  Linear_Solver = 'CHOL'

INTEGER                                                 ::  Laplace_NNZ
INTEGER                                                 ::  Factored_NNZ
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:,:)        ::  Laplace_Matrix_Full ! Should these be reals?
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:)          ::  Laplace_Matrix_Beta

COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:)          ::  Laplace_Matrix_VAL
INTEGER,        ALLOCATABLE,    DIMENSION(:,:)          ::  Laplace_Matrix_ROW
INTEGER,        ALLOCATABLE,    DIMENSION(:,:)          ::  Laplace_Matrix_COL

COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:)          ::  Laplace_Factored_VAL
INTEGER,        ALLOCATABLE,    DIMENSION(:,:)          ::  Laplace_Factored_ROW
INTEGER,        ALLOCATABLE,    DIMENSION(:,:)          ::  Laplace_Factored_COL


COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:)          ::  zMA_First_Col_Storage
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:)          ::  zMA_Last_Col_Storage

!
! Type B Matrix Variables
!

INTEGER                                                 ::  iMB_Diagonals
INTEGER                                                 ::  iMB_Bandwidth
INTEGER,        ALLOCATABLE,    DIMENSION(:)            ::  iMB_IPIV
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:)          ::  zMB_Matrix_Banded
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:)            ::  zMB_Matrix_Diagonal

COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:,:)        ::  zMB_First_Col_Storage
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:,:)        ::  zMB_Last_Col_Storage


CONTAINS




END MODULE Variables_Matrices
