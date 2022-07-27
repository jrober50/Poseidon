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


CHARACTER(LEN = 4)                                      :: Matrix_Format = 'CCS'
CHARACTER(LEN = 4)                                      :: Linear_Solver = 'CHOL'

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


COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:)          ::  First_Column_Storage
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:)          ::  Last_Column_Storage

INTEGER                                                 ::  Beta_Diagonals
INTEGER                                                 ::  Beta_Bandwidth
INTEGER,        ALLOCATABLE,    DIMENSION(:)            ::  Beta_IPIV
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:)          ::  Beta_MVL_Banded
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:)            ::  Beta_MVL_Diagonal
LOGICAL                                                 ::  Beta_Factorized_Flag = .FALSE.

COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:,:)        ::  First_Column_Beta_Storage
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:,:)        ::  Last_Column_Beta_Storage

CONTAINS




END MODULE Variables_Matrices
