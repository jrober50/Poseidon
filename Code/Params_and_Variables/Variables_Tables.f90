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
INTEGER,             PUBLIC, ALLOCATABLE, DIMENSION(:)              ::  M_VALUES

REAL(idp),  PUBLIC, ALLOCATABLE, DIMENSION(:,:,:)                   ::  Plm_Values
REAL(idp),  PUBLIC, ALLOCATABLE, DIMENSION(:,:,:)                   ::  Plm_dt_Values

REAL(idp),  PUBLIC, ALLOCATABLE, DIMENSION(:)                       ::  Nlm_Values

REAL(idp),  PUBLIC, ALLOCATABLE, DIMENSION(:,:,:)                   ::  Am_Values
REAL(idp),  PUBLIC, ALLOCATABLE, DIMENSION(:,:,:)                   ::  Am_dp_Values

REAL(idp),   PUBLIC, ALLOCATABLE, DIMENSION(:,:)                    ::  Slm_Elem_Values
REAL(idp),   PUBLIC, ALLOCATABLE, DIMENSION(:,:)                    ::  Slm_Elem_dt_Values
REAL(idp),   PUBLIC, ALLOCATABLE, DIMENSION(:,:)                    ::  Slm_Elem_dp_Values


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
REAL(idp),      PUBLIC, ALLOCATABLE, DIMENSION(:,:)         ::  Level_dx
INTEGER,        PUBLIC, ALLOCATABLE, DIMENSION(:)           ::  Level_Ratios


END MODULE Variables_Tables
