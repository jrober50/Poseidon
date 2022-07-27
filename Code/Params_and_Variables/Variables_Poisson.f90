   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Variables_Poisson                                                     !##!
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



COMPLEX(idp), ALLOCATABLE, DIMENSION(:,:,:)     ::  Coefficient_Vector
COMPLEX(idp), ALLOCATABLE, DIMENSION(:,:)       ::  Source_Vector
REAL(idp),    ALLOCATABLE, DIMENSION(:,:,:,:)   ::  Source_Terms


REAL(idp),    ALLOCATABLE, DIMENSION(:,:)       ::  First_Column_Storage
REAL(idp),    ALLOCATABLE, DIMENSION(:,:)       ::  Last_Column_Storage

LOGICAL                         ::  Stiffness_Matrix_Initialized_Flag = .FALSE.
LOGICAL                         ::  Test_Space_Allocated_Flag = .FALSE.
LOGICAL                         ::  Test_Run_Flag = .FALSE.
LOGICAL                         ::  Matrix_Cholesky_Factorized_Flag = .FALSE.



REAL(idp), ALLOCATABLE, DIMENSION(:,:,:)            :: STF_MAT_Integrals


!!! STF_MAT in full matrix form !!!
REAL(idp), ALLOCATABLE, DIMENSION(:,:,:)            :: STF_MAT



!!! STF_MAT in CCS form !!!
INTEGER                                             :: STF_NNZ
REAL(idp), ALLOCATABLE, DIMENSION(:,:)              :: STF_ELEM_VAL
INTEGER, ALLOCATABLE, DIMENSION(:)                  :: STF_COL_PTR, STF_ROW_IND



LOGICAL                                         :: INNER_BC_SET_FLAG = .FALSE.
LOGICAL                                         :: OUTER_BC_SET_FLAG = .FALSE.




END MODULE Variables_Poisson
