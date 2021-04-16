   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Variables_Tables                                                             !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains parameters used to define the running of Poseidon.                 !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!

USE Poseidon_Kinds_Module, &
                ONLY : idp

!===================================================================!
!                                                                   !
!   Ylm Table                                                       !
!                                                                   !
!===================================================================!
COMPLEX(KIND = idp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:,:,:,:)    :: Ylm_Table_Block

INTEGER,             PUBLIC, ALLOCATABLE, DIMENSION(:)              :: M_VALUES

COMPLEX(KIND = idp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:,:)        :: Ylm_Values
COMPLEX(KIND = idp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:,:)        :: Ylm_dt_Values
COMPLEX(KIND = idp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:,:)        :: Ylm_dp_Values

COMPLEX(KIND = idp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:,:)        :: Ylm_CC_Values
COMPLEX(KIND = idp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:,:)        :: Ylm_CC_DT_Values
COMPLEX(KIND = idp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:,:)        :: Ylm_CC_DP_Values





!===================================================================!
!                                                                   !
!   Lagrange_Poly Table                                             !
!                                                                   !
!===================================================================!
REAL(KIND = idp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:)             :: Lagrange_Poly_Table
REAL(KIND = idp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:,:,:)         :: LPT_LPT

END MODULE Variables_Tables
