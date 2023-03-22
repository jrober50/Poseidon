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
REAL(idp),      ALLOCATABLE,    DIMENSION(:,:,:)        ::  Laplace_Matrix_Full
REAL(idp),      ALLOCATABLE,    DIMENSION(:,:)          ::  Laplace_Matrix_Beta

REAL(idp),      ALLOCATABLE,    DIMENSION(:,:)          ::  Laplace_Matrix_VAL
INTEGER,        ALLOCATABLE,    DIMENSION(:,:)          ::  Laplace_Matrix_ROW
INTEGER,        ALLOCATABLE,    DIMENSION(:,:)          ::  Laplace_Matrix_COL

REAL(idp),      ALLOCATABLE,    DIMENSION(:,:)          ::  Laplace_Factored_VAL
INTEGER,        ALLOCATABLE,    DIMENSION(:,:)          ::  Laplace_Factored_ROW
INTEGER,        ALLOCATABLE,    DIMENSION(:,:)          ::  Laplace_Factored_COL


REAL(idp),      ALLOCATABLE,    DIMENSION(:,:)          ::  dMA_First_Col_Storage
REAL(idp),      ALLOCATABLE,    DIMENSION(:,:)          ::  dMA_Last_Col_Storage


!
! Type B Matrix Variables
!

INTEGER                                                 ::  iMB_Diagonals
INTEGER                                                 ::  iMB_Bandwidth
INTEGER,        ALLOCATABLE,    DIMENSION(:)            ::  iMB_IPIV

REAL(idp),      ALLOCATABLE,    DIMENSION(:,:)          ::  dMB_Matrix_Banded
REAL(idp),      ALLOCATABLE,    DIMENSION(:)            ::  dMB_Matrix_Diagonal

REAL(idp),      ALLOCATABLE,    DIMENSION(:,:,:)        ::  dMB_First_Col_Storage
REAL(idp),      ALLOCATABLE,    DIMENSION(:,:,:)        ::  dMB_Last_Col_Storage


CONTAINS




END MODULE Variables_Matrices
