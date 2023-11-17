   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Load_Vector_Variables_Module                                          !##!
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
            ONLY :  idp

IMPLICIT NONE



REAL( idp ), ALLOCATABLE, DIMENSION(:)             :: Cur_R_Locs
REAL( idp ), ALLOCATABLE, DIMENSION(:)             :: Cur_T_Locs
REAL( idp ), ALLOCATABLE, DIMENSION(:)             :: Cur_P_Locs

REAL( idp ), ALLOCATABLE, DIMENSION(:)             :: R_Square

REAL( idp ), ALLOCATABLE, DIMENSION(:)             :: TP_Sin_Val
REAL( idp ), ALLOCATABLE, DIMENSION(:)             :: TP_Sin_Square
REAL( idp ), ALLOCATABLE, DIMENSION(:)             :: TP_Cotan_Val
REAL( idp ), ALLOCATABLE, DIMENSION(:,:)           :: TP_RSin_Square

REAL( idp ), ALLOCATABLE, DIMENSION(:)             :: R_Int_Weights
REAL( idp ), ALLOCATABLE, DIMENSION(:)             :: TP_Int_Weights


REAL( idp ), ALLOCATABLE, DIMENSION(:,:)           :: Cur_Val_Psi
REAL( idp ), ALLOCATABLE, DIMENSION(:,:)           :: Cur_Val_AlphaPsi
REAL( idp ), ALLOCATABLE, DIMENSION(:,:,:)         :: Cur_Val_X
REAL( idp ), ALLOCATABLE, DIMENSION(:,:,:)         :: Cur_Val_Beta

REAL( idp ), ALLOCATABLE, DIMENSION(:,:,:)         :: Cur_Drv_Psi
REAL( idp ), ALLOCATABLE, DIMENSION(:,:,:)         :: Cur_Drv_AlphaPsi
REAL( idp ), ALLOCATABLE, DIMENSION(:,:,:,:)       :: Cur_Drv_X
REAL( idp ), ALLOCATABLE, DIMENSION(:,:,:,:)       :: Cur_Drv_Beta


REAL( idp ), ALLOCATABLE, DIMENSION(:,:,:)         :: SourceTerm



CONTAINS

END MODULE Load_Vector_Variables_Module
