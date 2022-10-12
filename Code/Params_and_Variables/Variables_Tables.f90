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
REAL(idp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:)            :: Lagrange_Poly_Table
REAL(idp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:,:,:)        :: LPT_LPT



!===================================================================!
!                                                                   !
!   AMReX Tables                                                    !
!                                                                   !
!===================================================================!
COMPLEX(idp),   PUBLIC, ALLOCATABLE, DIMENSION(:,:)         ::  Ylm_Sqrt_Table
REAL(idp),      PUBLIC, ALLOCATABLE, DIMENSION(:,:)         ::  Ylm_Norm_Table
REAL(idp),      PUBLIC, ALLOCATABLE, DIMENSION(:,:,:,:)     ::  rBT_NormedLegendre
!REAL(idp),      PUBLIC, ALLOCATABLE, DIMENSION(:,:,:,:)     ::  rBT_NormedLegendre_dt
!REAL(idp),      PUBLIC, ALLOCATABLE, DIMENSION(:,:,:,:)     ::  rBT_NormedLegendre_CC

COMPLEX(idp),   PUBLIC, ALLOCATABLE, DIMENSION(:,:)         ::  Ylm_Elem_Values
COMPLEX(idp),   PUBLIC, ALLOCATABLE, DIMENSION(:,:)         ::  Ylm_Elem_dt_Values
COMPLEX(idp),   PUBLIC, ALLOCATABLE, DIMENSION(:,:)         ::  Ylm_Elem_dp_Values
COMPLEX(idp),   PUBLIC, ALLOCATABLE, DIMENSION(:,:)         ::  Ylm_Elem_CC_Values

REAL(idp),      PUBLIC, ALLOCATABLE, DIMENSION(:,:)         ::  Level_dx
INTEGER,        PUBLIC, ALLOCATABLE, DIMENSION(:)           ::  Level_Ratios


END MODULE Variables_Tables
