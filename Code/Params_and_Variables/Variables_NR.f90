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



COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)      ::  NR_Coeff_Vector
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)      ::  NR_Update_Vector
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)      ::  NR_Source_Vector
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)      ::  Block_RHS_Vector

COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:)    ::  Block_STF_Mat
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:)    ::  BLOCK_ELEM_STF_MATVEC


END MODULE Variables_NR


