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


INTEGER, DIMENSION(1:5)                                         ::  CFA_EQ_Flags
INTEGER, DIMENSION(1:5)                                         ::  CFA_EQ_Map
INTEGER                                                         ::  Laplace_NNZ

CHARACTER(LEN = 4)                                              :: Matrix_Format = 'FULL'
CHARACTER(LEN = 4)                                              :: Linear_Solver = 'FULL'

INTEGER                                                         :: Num_Matrices







END MODULE Poseidon_FP_Variables_Module

