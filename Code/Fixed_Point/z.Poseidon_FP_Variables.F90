   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_FP_Variables_Module                                                 !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!


USE Poseidon_Constants_Module, &
            ONLY : idp, pi



IMPLICIT NONE




COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)              ::  FP_Source_Vector
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)              ::  FP_Coeff_Vector
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)              ::  FP_Update_Vector
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)              ::  FP_Laplace_Vector
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)              ::  FP_Residual_Vector

COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:,:)            ::  Laplace_Matrix_Full

COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)              ::  Laplace_Matrix_VAL
INTEGER, ALLOCATABLE, DIMENSION(:,:)                            ::  Laplace_Matrix_ROW
INTEGER, ALLOCATABLE, DIMENSION(:,:)                            ::  Laplace_Matrix_COL

COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)              ::  Laplace_Factored_VAL
INTEGER, ALLOCATABLE, DIMENSION(:,:)                            ::  Laplace_Factored_ROW
INTEGER, ALLOCATABLE, DIMENSION(:,:)                            ::  Laplace_Factored_COL

COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)              ::  First_Column_Storage
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)              ::  Last_Column_Storage

INTEGER, DIMENSION(1:5)                                         ::  CFA_EQ_Flags
INTEGER, DIMENSION(1:5)                                         ::  CFA_EQ_Map
INTEGER, DIMENSION(1:5)                                         ::  CFA_MAT_Map
INTEGER                                                         ::  Laplace_NNZ
INTEGER                                                         ::  Factored_NNZ


CHARACTER(LEN = 4)                                              :: Matrix_Format = 'Full'
CHARACTER(LEN = 4)                                              :: Linear_Solver = 'Full'

INTEGER                                                         :: Num_Matrices


INTEGER                                                         ::  MCF_Flag




END MODULE Poseidon_FP_Variables_Module

