   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Variables_NR                                                                 !##!
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



COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)      ::  Coefficient_Vector
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)      ::  Update_Vector

COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)      ::  RHS_Vector
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)      ::  Block_RHS_Vector

COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:)    ::  Block_STF_Mat
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:)    ::  BLOCK_ELEM_STF_MATVEC


END MODULE Variables_NR


