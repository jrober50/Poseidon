   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE XCFC_AMReX_Source_Vector_TypeB_Module                                 !##!
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

#ifdef POSEIDON_AMREX_FLAG
use amrex_base_module

USE amrex_box_module,       ONLY: &
  amrex_box
USE amrex_boxarray_module,  ONLY: &
  amrex_boxarray,         &
  amrex_boxarray_build,   &
  amrex_boxarray_destroy
USE amrex_distromap_module, ONLY: &
  amrex_distromap,       &
  amrex_distromap_build, &
  amrex_distromap_destroy
USE amrex_multifab_module,  ONLY: &
  amrex_multifab, &
  amrex_multifab_build
#endif

USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Poseidon_Numbers_Module, &
            ONLY :  pi,                         &
                    TwoPi

USE Units_Module, &
            ONLY :  GR_Source_Scalar

USE Poseidon_Parameters, &
            ONLY :  DEGREE,                     &
                    L_LIMIT,                    &
                    NUM_CFA_EQs,                &
                    NUM_CFA_VARs

USE Parameters_Variable_Indices, &
            ONLY :  iU_CF,                      &
                    iU_LF,                      &
                    iU_S1,                      &
                    iU_S2,                      &
                    iU_S3,                      &
                    iU_X1,                      &
                    iU_X2,                      &
                    iU_X3,                      &
                    iS_E,                       &
                    iS_S,                       &
                    iS_S1,                      &
                    iS_S2,                      &
                    iS_S3


USE Variables_Quadrature, &
            ONLY :  NUM_R_QUAD_POINTS,          &
                    NUM_T_QUAD_POINTS,          &
                    NUM_P_QUAD_POINTS,          &
                    Num_Quad_DOF,               &
                    NUM_TP_QUAD_POINTS,         &
                    INT_R_LOCATIONS,            &
                    INT_T_LOCATIONS,            &
                    INT_P_LOCATIONS,            &
                    INT_R_WEIGHTS,              &
                    INT_T_WEIGHTS,              &
                    INT_P_WEIGHTS,              &
                    INT_TP_WEIGHTS

USE Variables_Mesh, &
            ONLY :  NUM_R_ELEMENTS,             &
                    NUM_T_ELEMENTS,             &
                    NUM_P_ELEMENTS,             &
                    rlocs,                      &
                    tlocs,                      &
                    plocs
                  
USE Variables_Source, &
            ONLY :  Block_Source_E,             &
                    Block_Source_S,             &
                    Block_Source_Si

USE Variables_Tables, &
            ONLY :  Ylm_Values,                 &
                    Ylm_dt_Values,              &
                    Ylm_dp_Values,              &
                    Ylm_CC_Values,              &
                    Lagrange_Poly_Table

USE Variables_Derived, &
            ONLY :  NUM_R_NODES,                &
                    LM_LENGTH,                  &
                    Beta_Prob_Dim

USE Variables_FP, &
            ONLY :  FP_Coeff_Vector,            &
                    FP_Coeff_Vector_Beta,       &
                    FP_Coeff_Vector_Orig,       &
                    FP_Source_Vector_A,         &
                    FP_Source_Vector_B,         &
                    CFA_EQ_Map,                 &
                    CFA_EQ_Flags

USE Functions_Jacobian, &
            ONLY :  JCBN_kappa_FUNCTION_3D_ALL,     &
                    JCBN_BIGK_FUNCTION,         &
                    Calc_Ahat

USE Poseidon_IO_Module, &
            ONLY :  Clock_In

USE Functions_Mapping, &
            ONLY :  Map_To_X_Space

USE FP_Functions_Mapping, &
            ONLY :  FP_FEM_Node_Map,            &
                    FP_Beta_Array_Map,          &
                    FP_X_Array_Map,             &
                    FP_Array_Map,               &
                    FP_LM_Map,                  &
                    FP_tpd_Map

USE Variables_AMReX_Core,   &
            ONLY :  AMReX_Levels

#ifdef DPOSEIDON_AMREX_FLAG
USE Variables_AMReX_Multifabs, &
            ONLY :  MF_Source,  &
                    BA_Source,  &
                    DM_Source
#endif


IMPLICIT NONE

REAL(KIND = idp)                                        :: Time_S

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CUR_R_LOCS
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CUR_T_LOCS
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CUR_P_LOCS

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: R_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: R_CUBED

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: TP_Sin_Val
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: SIN_VAL
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: SIN_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: COS_VAL
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: COS_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CSC_VAL
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CSC_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: COTAN_VAL

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: TP_SIN_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: TP_Cotan_Val

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)           :: RSIN_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)           :: TP_RSIN_SQUARE

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: R_Int_Weights
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: TP_Int_Weights


REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)           :: Orig_VAL_PSI
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)           :: CUR_VAL_PSI
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)           :: CUR_VAL_ALPHAPSI
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: CUR_VAL_X
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: CUR_VAL_BETA

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: CUR_DRV_PSI
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: CUR_DRV_ALPHAPSI
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:,:)       :: CUR_DRV_X
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:,:)       :: CUR_DRV_BETA


REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: StaredSource
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: SourceTerm

TYPE(amrex_mfiter)                                      :: mfi
TYPE(amrex_box)                                         :: Box
REAL(idp), CONTIGUOUS, POINTER                          :: Source_PTR(:,:,:,:)

INTEGER, DIMENSION(3)                                   :: ELo, EHi


CONTAINS


#ifdef DPOSEIDON_AMREX_FLAG
!+101+###########################################################################!
!                                                                                !
!           XCFC_Calc_X_Source                                                   !
!                                                                                !
!################################################################################!
SUBROUTINE XCFC_AMReX_Calc_Source_Vector_TypeB( iU, iVB )

INTEGER, INTENT(IN), DIMENSION(3)               ::  iU
INTEGER, INTENT(IN)                             ::  iVB


INTEGER                                         ::  re, te, pe,     &
                                                    rd, tpd, td, pd

REAL(KIND = idp)                                ::  DROT,     &
                                                    DTOT,     &
                                                    DPOT

INTEGER                                         :: nComps

INTEGER                                         :: Here, There, lvl




nComps = 5*Num_R_Quad_Points*Num_T_Quad_Points*Num_P_Quad_Points


DO lvl = 0,AMReX_Levels-1
    CALL amrex_mfiter_build(mfi, MF_Source(lvl), tiling = .false. )
    DO WHILE(mfi%next())

        Source_PTR => MF_Source(lvl)%dataPtr(mfi)
        Box = mfi%tilebox()

        ELo = Box%lo
        EHi = Box%hi

        

        DO PE = ELo(3), EHi(3)
        DO TE = ELo(3), EHi(3)
        DO RE = ELo(3), EHi(3)


            DROT = 0.5_idp *(rlocs(re + 1) - rlocs(re))
            DPOT = 0.5_idp * (plocs(pe + 1) - plocs(pe))
            DTOT = 0.5_idp * (tlocs(te + 1) - tlocs(te))

            CUR_R_LOCS(:) = deltar_overtwo * (INT_R_LOCATIONS(:)+1.0_idp) + rlocs(re)
            CUR_T_LOCS(:) = deltat_overtwo * (INT_T_LOCATIONS(:)+1.0_idp) + tlocs(te)

            R_SQUARE(:) = CUR_R_LOCS(:)*CUR_R_LOCS(:)
            DO td = 1,NUM_T_QUAD_POINTS
            DO pd = 1,NUM_P_QUAD_POINTS
                tpd = FP_tpd_Map(td,pd)
                TP_Sin_Val(tpd)    = DSIN(CUR_T_LOCS(td))
                TP_Cotan_Val(tpd)  = 1.0_idp/DTAN(CUR_T_LOCS(td))
            END DO
            END DO
            TP_Sin_Square(:) = TP_Sin_Val(:)*TP_Sin_Val


            DO rd = 1,NUM_R_QUAD_POINTS
                TP_RSIN_SQUARE(:,rd) = R_SQUARE(rd)*TP_SIN_SQUARE(:)
            END DO


            CALL Calc_Int_Weights( DROT, DTOT,                  &
                                   R_Square, TP_Sin_Val,        &
                                   R_Int_Weights, TP_Int_Weights )

            CALL Calc_XCFC_CurVals_TypeB( re, te, pe,       &
                                          iU, iVB,          &
                                          DROT,DTOT, DPOT   )

            CALL Create_XCFC_Vector_TypeB( re, te, pe, iU, iVB )



        END DO ! RE
        END DO ! TE
        END DO ! PE

    END DO
    CALL amrex_mfiter_destroy(mfi)
END DO ! lvl


END SUBROUTINE XCFC_AMReX_Calc_X_Source










!+701+###########################################################################!
!                                                                                !
!           Get_Physical_Source                                      !
!                                                                                !
!################################################################################!
SUBROUTINE Get_Physical_Source( Source, Var, RE, TE, PE )

REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points), INTENT(OUT)    :: Source
INTEGER,                                                     INTENT(IN)     :: Var
INTEGER,                                                     INTENT(IN)     :: RE, TE, PE

INTEGER                                                         :: rd, td, pd, tpd
INTEGER                                                         :: Here, There

IF ( Var == 1 ) THEN

    DO rd = 1,NUM_R_QUAD_POINTS
    DO td = 1,NUM_T_QUAD_POINTS
    DO pd = 1,NUM_P_QUAD_POINTS

        tpd = FP_tpd_Map(td,pd)
        Here = (Var-1)*NUM_Quad_DOF                       &
             + (pd-1)*Num_R_Quad_Points*Num_T_Quad_Points   &
             + (td-1)*Num_R_Quad_Points                     &
             + (rd-1)

        Source(tpd,rd) = Source_PTR(re,te,pe,Here)


    END DO ! pd
    END DO ! td
    END DO ! rd

ELSEIF ( Var == 2 ) THEN

    DO rd = 1,NUM_R_QUAD_POINTS
    DO td = 1,NUM_T_QUAD_POINTS
    DO pd = 1,NUM_P_QUAD_POINTS

        tpd = FP_tpd_Map(td,pd)
        Here = (iS_E-1)*NUM_Quad_DOF                       &
             + (pd-1)*Num_R_Quad_Points*Num_T_Quad_Points   &
             + (td-1)*Num_R_Quad_Points                     &
             + (rd-1)
        There = (iS_S-1)*NUM_Quad_DOF                       &
                + (pd-1)*Num_R_Quad_Points*Num_T_Quad_Points   &
                + (td-1)*Num_R_Quad_Points                     &
                + (rd-1)
        Source(tpd,rd) = Source_PTR(re,te,pe,Here)  &
                    + 2.0_idp*Source_PTR(re,te,pe,There)


    END DO ! pd
    END DO ! td
    END DO ! rd

ELSE IF ( Var == 3) THEN

    DO rd = 1,NUM_R_QUAD_POINTS
    DO td = 1,NUM_T_QUAD_POINTS
    DO pd = 1,NUM_P_QUAD_POINTS

        tpd = FP_tpd_Map(td,pd)
        Here = (iS_S1-1)*NUM_Quad_DOF                       &
             + (pd-1)*Num_R_Quad_Points*Num_T_Quad_Points   &
             + (td-1)*Num_R_Quad_Points                     &
             + (rd-1)

        Source(tpd,rd) = Source_PTR(re,te,pe,Here)


    END DO ! pd
    END DO ! td
    END DO ! rd

ELSE IF ( Var == 4) THEN

    DO rd = 1,NUM_R_QUAD_POINTS
    DO td = 1,NUM_T_QUAD_POINTS
    DO pd = 1,NUM_P_QUAD_POINTS

        tpd = FP_tpd_Map(td,pd)
        Here = (iS_S1-1)*NUM_Quad_DOF                       &
             + (pd-1)*Num_R_Quad_Points*Num_T_Quad_Points   &
             + (td-1)*Num_R_Quad_Points                     &
             + (rd-1)

        Source(tpd,rd) = Source_PTR(re,te,pe,Here)


    END DO ! pd
    END DO ! td
    END DO ! rd

ELSE IF ( Var == 5 ) THEN

    DO rd = 1,NUM_R_QUAD_POINTS
    DO td = 1,NUM_T_QUAD_POINTS
    DO pd = 1,NUM_P_QUAD_POINTS

        tpd = FP_tpd_Map(td,pd)
        Here = (iS_S1-1)*NUM_Quad_DOF                       &
             + (pd-1)*Num_R_Quad_Points*Num_T_Quad_Points   &
             + (td-1)*Num_R_Quad_Points                     &
             + (rd-1)

        Source(tpd,rd) = Source_PTR(re,te,pe,Here)
        

    END DO ! pd
    END DO ! td
    END DO ! rd

END IF

END SUBROUTINE Get_Physical_Source





!+701+###########################################################################!
!                                                                                !
!           Allocate_Master_Build_Variables                                      !
!                                                                                !
!################################################################################!
SUBROUTINE Allocate_XCFC_AMReX_Source_Variables()



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


ALLOCATE( StaredSource( 1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:5) )
ALLOCATE( SourceTerm( 1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:8 ) )

END SUBROUTINE Allocate_XCFC_AMReX_Source_Variables



!+702+###########################################################################!
!                                                                                !
!           Deallocate_Master_Build_Variables                                    !
!                                                                                !
!################################################################################!
SUBROUTINE Deallocate_XCFC_AMReX_Source_Variables()

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


DEALLOCATE( StaredSource )
DEALLOCATE( SourceTerm )

END SUBROUTINE Deallocate_XCFC_AMReX_Source_Variables



#endif
END MODULE XCFC_AMReX_Source_Vector_TypeB_Module
