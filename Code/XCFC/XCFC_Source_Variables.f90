   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE XCFC_Source_Variables_Module                                          !##!
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

USE Variables_Quadrature, &
            ONLY :  Num_R_Quad_Points,          &
                    Num_T_Quad_Points,          &
                    Num_P_Quad_Points,          &
                    Num_TP_Quad_Points


IMPLICIT NONE


REAL( idp )                                        :: Time_S

REAL( idp ), ALLOCATABLE, DIMENSION(:)             :: CUR_R_LOCS
REAL( idp ), ALLOCATABLE, DIMENSION(:)             :: CUR_T_LOCS
REAL( idp ), ALLOCATABLE, DIMENSION(:)             :: CUR_P_LOCS

REAL( idp ), ALLOCATABLE, DIMENSION(:)             :: R_SQUARE
REAL( idp ), ALLOCATABLE, DIMENSION(:)             :: R_CUBED

REAL( idp ), ALLOCATABLE, DIMENSION(:)             :: TP_Sin_Val
REAL( idp ), ALLOCATABLE, DIMENSION(:)             :: SIN_VAL
REAL( idp ), ALLOCATABLE, DIMENSION(:)             :: SIN_SQUARE
REAL( idp ), ALLOCATABLE, DIMENSION(:)             :: COS_VAL
REAL( idp ), ALLOCATABLE, DIMENSION(:)             :: COS_SQUARE
REAL( idp ), ALLOCATABLE, DIMENSION(:)             :: CSC_VAL
REAL( idp ), ALLOCATABLE, DIMENSION(:)             :: CSC_SQUARE
REAL( idp ), ALLOCATABLE, DIMENSION(:)             :: COTAN_VAL

REAL( idp ), ALLOCATABLE, DIMENSION(:)             :: TP_SIN_SQUARE
REAL( idp ), ALLOCATABLE, DIMENSION(:)             :: TP_Cotan_Val

REAL( idp ), ALLOCATABLE, DIMENSION(:,:)           :: RSIN_SQUARE
REAL( idp ), ALLOCATABLE, DIMENSION(:,:)           :: TP_RSIN_SQUARE

REAL( idp ), ALLOCATABLE, DIMENSION(:)             :: R_Int_Weights
REAL( idp ), ALLOCATABLE, DIMENSION(:)             :: TP_Int_Weights


REAL( idp ), ALLOCATABLE, DIMENSION(:,:)           :: Orig_VAL_PSI
REAL( idp ), ALLOCATABLE, DIMENSION(:,:)           :: CUR_VAL_PSI
REAL( idp ), ALLOCATABLE, DIMENSION(:,:)           :: CUR_VAL_ALPHAPSI
REAL( idp ), ALLOCATABLE, DIMENSION(:,:,:)         :: CUR_VAL_X
REAL( idp ), ALLOCATABLE, DIMENSION(:,:,:)         :: CUR_VAL_BETA

REAL( idp ), ALLOCATABLE, DIMENSION(:,:,:)         :: CUR_DRV_PSI
REAL( idp ), ALLOCATABLE, DIMENSION(:,:,:)         :: CUR_DRV_ALPHAPSI
REAL( idp ), ALLOCATABLE, DIMENSION(:,:,:,:)       :: CUR_DRV_X
REAL( idp ), ALLOCATABLE, DIMENSION(:,:,:,:)       :: CUR_DRV_BETA


REAL( idp ), ALLOCATABLE, DIMENSION(:,:,:)         :: SourceTerm



CONTAINS


!+701+###########################################################################!
!                                                                                !
!           Allocate_Master_Build_Variables                                      !
!                                                                                !
!################################################################################!
SUBROUTINE Allocate_XCFC_Source_Variables()

ALLOCATE( CUR_R_LOCS(1:NUM_R_QUAD_POINTS) )
ALLOCATE( CUR_T_LOCS(1:NUM_T_QUAD_POINTS) )
ALLOCATE( CUR_P_LOCS(1:NUM_P_QUAD_POINTS) )


ALLOCATE( R_SQUARE(1:NUM_R_QUAD_POINTS) )
ALLOCATE( R_CUBED(1:NUM_R_QUAD_POINTS) )


ALLOCATE( SIN_VAL( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( SIN_SQUARE( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( COS_VAL( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( COS_SQUARE( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( CSC_VAL( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( CSC_SQUARE( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( COTAN_VAL( 1:NUM_T_QUAD_POINTS ) )

ALLOCATE( RSIN_SQUARE( 1:NUM_T_QUAD_POINTS, 1:NUM_R_QUAD_POINTS ) )

ALLOCATE( TP_SIN_VAL( 1:NUM_TP_QUAD_POINTS ) )
ALLOCATE( TP_SIN_SQUARE( 1:NUM_TP_QUAD_POINTS ) )
ALLOCATE( TP_Cotan_VAL( 1:NUM_TP_QUAD_POINTS ) )
ALLOCATE( TP_RSIN_SQUARE( 1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS ) )


ALLOCATE( R_Int_Weights( 1:NUM_R_QUAD_POINTS ) )
ALLOCATE( TP_Int_Weights( 1:NUM_TP_QUAD_POINTS) )


ALLOCATE( Orig_VAL_PSI(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS) )

ALLOCATE( CUR_VAL_PSI(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS) )
ALLOCATE( CUR_VAL_ALPHAPSI(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS) )
ALLOCATE( CUR_VAL_BETA(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:3) )
ALLOCATE( CUR_VAL_X(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:3) )

ALLOCATE( CUR_DRV_PSI(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:3) )
ALLOCATE( CUR_DRV_ALPHAPSI(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:3) )
ALLOCATE( CUR_DRV_BETA(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:3, 1:3) )
ALLOCATE( CUR_DRV_X(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:3, 1:3) )


ALLOCATE( SourceTerm( 1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:8 ) )

END SUBROUTINE Allocate_XCFC_Source_Variables



!+702+###########################################################################!
!                                                                                !
!           Deallocate_Master_Build_Variables                                    !
!                                                                                !
!################################################################################!
SUBROUTINE Deallocate_XCFC_Source_Variables()

DEALLOCATE( CUR_R_LOCS )
DEALLOCATE( CUR_T_LOCS )
DEALLOCATE( CUR_P_LOCS )

DEALLOCATE( R_SQUARE )
DEALLOCATE( R_CUBED )

DEALLOCATE( TP_Sin_Val )
DEALLOCATE( SIN_VAL )
DEALLOCATE( SIN_SQUARE )
DEALLOCATE( COS_VAL )
DEALLOCATE( COS_SQUARE )
DEALLOCATE( CSC_VAL )
DEALLOCATE( CSC_SQUARE )
DEALLOCATE( COTAN_VAL )

DEALLOCATE( RSIN_SQUARE )

DEALLOCATE( TP_SIN_SQUARE )
DEALLOCATE( TP_Cotan_Val )
DEALLOCATE( TP_RSIN_SQUARE )

DEALLOCATE( R_Int_Weights )
DEALLOCATE( TP_Int_Weights )

DEALLOCATE( Orig_VAL_PSI )

DEALLOCATE( CUR_VAL_PSI )
DEALLOCATE( CUR_VAL_ALPHAPSI )
DEALLOCATE( CUR_VAL_BETA )
DEALLOCATE( CUR_VAL_X )

DEALLOCATE( CUR_DRV_PSI )
DEALLOCATE( CUR_DRV_ALPHAPSI )
DEALLOCATE( CUR_DRV_BETA )
DEALLOCATE( CUR_DRV_X )


DEALLOCATE( SourceTerm )

END SUBROUTINE Deallocate_XCFC_Source_Variables

END MODULE XCFC_Source_Variables_Module
