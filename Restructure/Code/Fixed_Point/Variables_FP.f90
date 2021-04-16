   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Variables_FP                                                                 !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!


USE Poseidon_Kinds_Module, &
            ONLY : idp



IMPLICIT NONE




COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:,:)        ::  FP_Source_Vector
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:,:)        ::  FP_Coeff_Vector
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:,:)        ::  FP_Update_Vector
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:,:)        ::  FP_Laplace_Vector
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:,:)        ::  FP_Residual_Vector

COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:,:)        ::  Laplace_Matrix_Full ! Should these be reals?
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:)          ::  Laplace_Matrix_Beta
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:)            ::  FP_Source_Vector_Beta
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:)            ::  FP_Coeff_Vector_Beta
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:)            ::  FP_Laplace_Vector_Beta
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:)            ::  FP_Residual_Vector_Beta

! CCS Variables
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:,:)        ::  Laplace_Matrix_VAL
INTEGER,        ALLOCATABLE,    DIMENSION(:,:)          ::  Laplace_Matrix_ROW
INTEGER,        ALLOCATABLE,    DIMENSION(:,:)          ::  Laplace_Matrix_COL

COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:,:)        ::  Laplace_Factored_VAL
INTEGER,        ALLOCATABLE,    DIMENSION(:,:)          ::  Laplace_Factored_ROW
INTEGER,        ALLOCATABLE,    DIMENSION(:,:)          ::  Laplace_Factored_COL


COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:,:)        ::  First_Column_Storage
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:,:)        ::  Last_Column_Storage

INTEGER                                                 ::  Beta_Diagonals
INTEGER                                                 ::  Beta_Bandwidth
INTEGER,        ALLOCATABLE,    DIMENSION(:)            ::  Beta_IPIV
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:)          ::  Beta_MVL_Banded
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:)            ::  Beta_MVL_Diagonal
LOGICAL                                                 ::  Beta_Factorized_Flag = .FALSE.

COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:,:)        ::  First_Column_Beta_Storage
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:,:)        ::  Last_Column_Beta_Storage


! Flags & Maps
INTEGER, DIMENSION(1:5)                                 ::  CFA_EQ_Flags
INTEGER, DIMENSION(1:5)                                 ::  CFA_EQ_Map
INTEGER, DIMENSION(1:5)                                 ::  CFA_MAT_Map


INTEGER                                                 ::  Laplace_NNZ
INTEGER                                                 ::  Factored_NNZ


CHARACTER(LEN = 4)                                      :: Matrix_Format = 'Full'
CHARACTER(LEN = 4)                                      :: Linear_Solver = 'Full'

!CHARACTER(LEN = 4)                                      :: Matrix_Format = 'CCS'
!CHARACTER(LEN = 4)                                      :: Linear_Solver = 'CHOL'


! Cholesky Factorizatio Flag
INTEGER                                                 ::  MCF_Flag

INTEGER                                                 ::  FP_Anderson_M = 3
INTEGER                                                 ::  Num_Matrices

END MODULE Variables_FP 

