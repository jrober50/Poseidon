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




COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)              ::  FP_Source_Vector
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)              ::  FP_Coeff_Vector
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)              ::  FP_Update_Vector
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)              ::  FP_Laplace_Vector
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)              ::  FP_Residual_Vector

COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)              ::  Laplace_Matrix_Full ! Should these be reals?
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:)                ::  Laplace_Matrix_Beta
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)                  ::  FP_Source_Vector_Beta
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)                  ::  FP_Coeff_Vector_Beta
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)                  ::  FP_Laplace_Vector_Beta
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)                  ::  FP_Residual_Vector_Beta

! CCS Variables
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)              ::  Laplace_Matrix_VAL
INTEGER, ALLOCATABLE, DIMENSION(:,:)                            ::  Laplace_Matrix_ROW
INTEGER, ALLOCATABLE, DIMENSION(:,:)                            ::  Laplace_Matrix_COL

COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)              ::  Laplace_Factored_VAL
INTEGER, ALLOCATABLE, DIMENSION(:,:)                            ::  Laplace_Factored_ROW
INTEGER, ALLOCATABLE, DIMENSION(:,:)                            ::  Laplace_Factored_COL

COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)              ::  First_Column_Storage
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)              ::  Last_Column_Storage



! Flags & Maps
INTEGER, DIMENSION(1:5)                                         ::  CFA_EQ_Flags
INTEGER, DIMENSION(1:5)                                         ::  CFA_EQ_Map
INTEGER, DIMENSION(1:5)                                         ::  CFA_MAT_Map


INTEGER                                                         ::  Laplace_NNZ
INTEGER                                                         ::  Factored_NNZ


CHARACTER(LEN = 4)                                              :: Matrix_Format = 'Full'
CHARACTER(LEN = 4)                                              :: Linear_Solver = 'Full'

INTEGER                                                         :: Num_Matrices

! Cholesky Factorizatio Flag
INTEGER                                                         ::  MCF_Flag

INTEGER                                                         ::  FP_Anderson_M = 1


END MODULE Variables_FP 

