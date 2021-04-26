   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE FP_Source_Beta                                                               !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!

!*D*================================!
!                                   !
!           Dependencies            !
!                                   !
!===================================!
USE Poseidon_Kinds_Module, &
            ONLY :  idp


USE Poseidon_Parameters, &
            ONLY :  DEGREE,                     &
                    NUM_CFA_EQs,                &
                    NUM_CFA_VARs

USE Variables_Quadrature, &
            ONLY :  NUM_R_QUAD_POINTS,          &
                    NUM_T_QUAD_POINTS,          &
                    NUM_P_QUAD_POINTS,          &
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
                    Ylm_CC_Values,              &
                    Lagrange_Poly_Table

USE Variables_Derived, &
            ONLY :  NUM_R_NODES,                &
                    LM_LENGTH

USE Variables_FP, &
            ONLY :  FP_Coeff_Vector,            &
                    FP_Coeff_Vector_Beta,       &
                    FP_Source_Vector,           &
                    FP_Source_Vector_Beta,      &
                    CFA_EQ_Map

USE Functions_Jacobian, &
            ONLY :  JCBN_kappa_FUNCTION_3D_ALL,     &
                    JCBN_BIGK_FUNCTION

USE Poseidon_IO_Module, &
            ONLY :  Clock_In

USE Functions_Mapping, &
            ONLY :  Map_To_X_Space

USE FP_Functions_Mapping, &
            ONLY :  FP_Vector_Map

USE SubJacobian_Functions_Module_3D, &
            ONLY :  Calc_RHS_Terms

USE SubJacobian_Functions_Module_1D, &
            ONLY :  Calc_RHS_Terms_1D





IMPLICIT NONE

REAL(KIND = idp)                                        :: Time_S

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CUR_R_LOCS
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CUR_T_LOCS
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CUR_P_LOCS

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: R_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: R_CUBED

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: SIN_VAL_B
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: SIN_VAL
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: SIN_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: TP_SIN_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: COS_VAL
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: COS_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CSC_VAL
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CSC_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: COTAN_VAL

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)           :: RSIN_SQUARE

COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)          :: PHI_EXP
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)          :: PHI_TWOEXP

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: R_Int_Weights
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: TP_Int_Weights

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)           :: CUR_VAL_PSI
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)           :: CUR_VAL_ALPHAPSI
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: CUR_VAL_BETA

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: CUR_DRV_PSI
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: CUR_DRV_ALPHAPSI
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:,:)       :: CUR_DRV_BETA

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)           :: Beta_DRV_Trace

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: SOURCE_TERMS

CONTAINS



!+101+###########################################################################!
!                                                                                !
!           Calc_Source_Vector                                                   !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_FP_Source_Vector()


REAL(KIND = idp),DIMENSION(1:3)                             ::  Timer

INTEGER                                                     ::  re, te, pe,     &
                                                                d, lm,          &
                                                                rd, tpd, td, pd

REAL(KIND = idp)                                                ::  TWOOVER_DELTAR,    &
                                                                    deltar_overtwo,     &
                                                                    deltat_overtwo,     &
                                                                    deltap_overtwo

REAL(KIND = idp), DIMENSION(1:3,1:3)                            ::  JCBN_kappa_Array
REAL(KIND = idp), DIMENSION(1:3)                                ::  JCBN_n_ARRAY

REAL(KIND = idp)                                                ::  JCBN_BIGK_VALUE

Timer = 0.0_idp
FP_Source_Vector = 0.0_idp
FP_Source_Vector_Beta = 0.0_idp

DO pe = 0,NUM_P_ELEMENTS-1
    deltap_overtwo = 0.5_idp * (plocs(pe + 1)-plocs(pe))
    CUR_P_LOCS(:) = deltap_overtwo * (INT_P_LOCATIONS+1.0_idp) + plocs(pe)

    PHI_EXP(:) = EXP( CMPLX(0, -CUR_P_LOCS(:), KIND = idp) )
    PHI_TWOEXP(:) = EXP( CMPLX(0, -2.0_idp*CUR_P_LOCS(:), KIND = idp) )

    DO te = 0,NUM_T_ELEMENTS-1

        deltat_overtwo = 0.5_idp*(tlocs(te + 1) - tlocs(te))
        CUR_T_LOCS(:) = deltat_overtwo * (INT_T_LOCATIONS(:)+1.0_idp) + tlocs(te)
        

        COTAN_VAL(:) = 1.0_idp/DTAN(CUR_T_LOCS(:))
        CSC_VAL(:) = 1.0_idp/DSIN(CUR_T_LOCS(:))
        SIN_VAL(:) = DSIN(CUR_T_LOCS(:))
        COS_VAL(:) = DCOS(CUR_T_LOCS(:))

        COS_SQUARE(:) = COS_VAL(:)*COS_VAL(:)
        SIN_SQUARE(:) = SIN_VAL(:)*SIN_VAL(:)
        CSC_SQUARE(:) = CSC_VAL(:)*CSC_VAL(:)

        DO td = 1,NUM_T_QUAD_POINTS
        DO pd = 1,NUM_P_QUAD_POINTS
            tpd = (td-1)*NUM_P_QUAD_POINTS + pd
            SIN_VAL_B(tpd) =   SIN_VAL(td)
            TP_SIN_SQUARE(tpd) = SIN_SQUARE(td)
        END DO
        END DO

        DO re = 0,NUM_R_ELEMENTS-1

            deltar_overtwo = 0.5_idp *(rlocs(re + 1) - rlocs(re))
            TWOOVER_DELTAR = 1.0_idp/deltar_overtwo
            CUR_R_LOCS(:) = deltar_overtwo * (INT_R_LOCATIONS(:)+1.0_idp) + rlocs(re)


            R_SQUARE(:) = CUR_R_LOCS(:)*CUR_R_LOCS(:)
            R_CUBED(:) = R_SQUARE(:)*CUR_R_LOCS(:)
            DO rd = 1,NUM_R_QUAD_POINTS

                RSIN_SQUARE(:,rd) = R_SQUARE(rd)*SIN_SQUARE(:)

            END DO

!            PRINT*,"Before Calc_FP_Current_Values"
            CALL Calc_FP_Current_Values(re, te, pe,               &
                                        DELTAR_OVERTWO, DELTAT_OVERTWO, DELTAP_OVERTWO  )

!            PRINT*,"Before Calc_FP_Source_Terms"
            CALL Calc_FP_Source_Terms( re, te, pe )

!            PRINT*,"Before Create_FP_Source_Vector"
            CALL Create_FP_Source_Vector( re, te, pe, DELTAR_OVERTWO, SIN_VAL_B )


        END DO ! re Loop
    END DO ! te Loop
END DO ! pe Loop



END SUBROUTINE Calc_FP_Source_Vector







!+202+###########################################################################!
!                                                                                !
!                  Calc_Current_Values          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_FP_Current_Values( re, te, pe,                                  &
                                    DELTAR_OVERTWO,                             &
                                    DELTAT_OVERTWO,                             &
                                    DELTAP_OVERTWO                              )

INTEGER, INTENT(IN)                                             ::  re, te, pe


REAL(KIND = idp), INTENT(IN)                                    ::  DELTAR_OVERTWO,     &
                                                                    DELTAT_OVERTWO,     &
                                                                    DELTAP_OVERTWO



COMPLEX(KIND = idp), DIMENSION(1:5)                             ::  Tmp_U_Value,        &
                                                                    Tmp_U_R_DRV_Value,  &
                                                                    Tmp_U_T_DRV_Value,  &
                                                                    Tmp_U_P_DRV_Value



INTEGER                                                         ::  tpd, td, pd, rd,    &
                                                                    lm, d, Here, ui



                          !                                                 !
                         !!                                                 !!
                        !!!          Initialize Local Quadratures           !!!
                         !!                                                 !!
                          !                                                 !
R_Int_Weights(:) = DELTAR_OVERTWO * R_SQUARE(:) * INT_R_WEIGHTS(:)

DO td = 1,NUM_T_QUAD_POINTS
DO pd = 1,NUM_P_QUAD_POINTS
!    TP_Int_Weights( (td-1)*NUM_P_QUAD_POINTS + pd ) = SIN_VAL(td)                           &
!                                                    * DELTAT_OVERTWO * INT_T_WEIGHTS(td)    &
!                                                    * DELTAP_OVERTWO * INT_P_WEIGHTS(pd)

    TP_Int_Weights( (td-1)*NUM_P_QUAD_POINTS + pd ) = SIN_VAL(td)                           &
                                                    * DELTAT_OVERTWO * INT_T_WEIGHTS(td)    &
                                                    * INT_P_WEIGHTS(pd)
END DO
END DO



DO rd = 1,NUM_R_QUAD_POINTS
DO tpd = 1,NUM_TP_QUAD_POINTS

    Tmp_U_Value = 0.0_idp
    Tmp_U_R_DRV_Value = 0.0_idp
    Tmp_U_T_DRV_Value = 0.0_idp
    Tmp_U_P_DRV_Value = 0.0_idp


    DO d = 0,DEGREE
    DO ui = 1,NUM_CFA_VARS
        Here = FP_Vector_Map(re,d)

        TMP_U_Value(ui)         = TMP_U_Value(ui)                           &
                                + SUM( FP_Coeff_Vector( Here, :, ui )      &
                                * Ylm_Values( :, tpd, te, pe )       )      &
                                * Lagrange_Poly_Table( d, rd, 0 )

        TMP_U_R_DRV_Value(ui)   = TMP_U_R_DRV_Value(ui)                     &
                                + SUM( FP_Coeff_Vector( Here, :, ui )     &
                                * Ylm_Values( :, tpd, te, pe )       )      &
                                * Lagrange_Poly_Table( d, rd, 1 )           &
                                / DELTAR_OVERTWO



    END DO ! ui
    END DO ! d


    CUR_VAL_PSI( tpd, rd )         = REAL(Tmp_U_Value(1), KIND = idp)
    CUR_DRV_PSI( tpd, rd, 1 )      = REAL(Tmp_U_R_DRV_Value(1), KIND = idp)
    CUR_DRV_PSI( tpd, rd, 2 )      = 0.0_idp
    CUR_DRV_PSI( tpd, rd, 3 )      = 0.0_idp


    CUR_VAL_ALPHAPSI( tpd, rd )    = REAL(Tmp_U_Value(2), KIND = idp)
    CUR_DRV_ALPHAPSI( tpd, rd, 1 ) = REAL(Tmp_U_R_DRV_Value(2), KIND = idp)
    CUR_DRV_ALPHAPSI( tpd, rd, 2 ) = 0.0_idp
    CUR_DRV_ALPHAPSI( tpd, rd, 3 ) = 0.0_idp


    CUR_VAL_BETA( tpd, rd, 1 )     = REAL(Tmp_U_Value(3), KIND = idp)
    CUR_VAL_BETA( tpd, rd, 2 )     = 0.0_idp
    CUR_VAL_BETA( tpd, rd, 3 )     = 0.0_idp


    CUR_DRV_BETA( tpd, rd, 1, 1 )  = REAL(Tmp_U_R_DRV_Value(3), KIND = idp)
    CUR_DRV_BETA( tpd, rd, 2, 1 )  = 0.0_idp
    CUR_DRV_BETA( tpd, rd, 3, 1 )  = 0.0_idp
    CUR_DRV_BETA( tpd, rd, 1, 2 )  = 0.0_idp
    CUR_DRV_BETA( tpd, rd, 2, 2 )  = 0.0_idp
    CUR_DRV_BETA( tpd, rd, 3, 2 )  = 0.0_idp
    CUR_DRV_BETA( tpd, rd, 1, 3 )  = 0.0_idp
    CUR_DRV_BETA( tpd, rd, 2, 3 )  = 0.0_idp
    CUR_DRV_BETA( tpd, rd, 3, 3 )  = 0.0_idp

    Beta_DRV_Trace( tpd, rd )      = 0.0_idp

END DO ! tpd
END DO ! rd


END SUBROUTINE Calc_FP_Current_Values





!+202+###########################################################################!
!                                                                                !
!                  Calc_FP_Source_Terms          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_FP_Source_Terms( re, te, pe )

INTEGER, INTENT(IN)                                                     ::  re, te, pe

REAL(KIND = idp)                                                        ::  REUSED_VALUE

INTEGER                                                                 ::  pd, td, rd,     &
                                                                            i, tpd


REAL(KIND = idp), DIMENSION(1:11)                                       ::  PSI_POWER
REAL(KIND = idp), DIMENSION(1:4)                                        ::  ALPHAPSI_POWER


REAL(KIND = idp), DIMENSION(1:3,1:3)                                    ::  JCBN_kappa_Array
REAL(KIND = idp), DIMENSION(1:3)                                        ::  JCBN_n_ARRAY

REAL(KIND = idp)                                                        ::  JCBN_BIGK_VALUE


JCBN_Kappa_Array = 0.0_idp
JCBN_n_Array     = 0.0_idp
JCBN_BIGK_Value  = 0.0_idp


DO rd = 1,NUM_R_QUAD_POINTS
DO td = 1,NUM_T_QUAD_POINTS
DO pd = 1,NUM_P_QUAD_POINTS

    tpd = (td-1)*NUM_P_QUAD_POINTS + pd


    PSI_POWER(1) = CUR_VAL_PSI( tpd, rd)
    DO i = 2,11
        PSI_POWER(i) = PSI_POWER(i-1)*PSI_POWER(1)
    END DO


    ALPHAPSI_POWER(1) = CUR_VAL_ALPHAPSI( tpd, rd)
    DO i = 2,4
        ALPHAPSI_POWER(i) = ALPHAPSI_POWER(i-1)*ALPHAPSI_POWER(1)
    END DO

    ! K_{ij}K^{ij} = Psi^{14}/AlphaPsi^{2} * BIGK
    !        PRINT*,"Before JCBN_BIGK_VALUE"
    JCBN_BIGK_VALUE = JCBN_BIGK_FUNCTION( rd, tpd,                                                        &
                                          NUM_R_QUAD_POINTS, NUM_TP_QUAD_POINTS,                          &
                                          CUR_VAL_BETA, CUR_DRV_BETA,                                     &
                                          CUR_R_LOCS(rd), R_SQUARE(rd), SIN_SQUARE(td), CSC_SQUARE(td),   &
                                          RSIN_SQUARE(td, rd), COTAN_VAL(td)                              )



    !        PRINT*,"Before  JCBN_kappa_Array"
    JCBN_kappa_Array = JCBN_kappa_FUNCTION_3D_ALL(  rd, tpd,                                        &
                                                    NUM_R_QUAD_POINTS, NUM_TP_QUAD_POINTS,          &
                                                    CUR_R_LOCS(rd), R_SQUARE(rd), R_CUBED(rd),      &
                                                    RSIN_SQUARE(td, rd),                            &
                                                    SIN_VAL(td), SIN_SQUARE(td), CSC_SQUARE(td),    &
                                                    COS_VAL(td), COTAN_VAL(td),                     &
                                                    CUR_VAL_BETA, CUR_DRV_BETA                      )


    !        PRINT*,"Before JCBN_n_ARRAY"
    JCBN_n_ARRAY(:) = CUR_DRV_ALPHAPSI( tpd, rd, : ) / ALPHAPSI_POWER(1)   &
                        - 7.0_idp * CUR_DRV_PSI( tpd, rd, : )/ PSI_POWER(1)

    JCBN_n_ARRAY(:) = CUR_DRV_ALPHAPSI( tpd, rd, : ) / ALPHAPSI_POWER(1)   &
                        - CUR_DRV_PSI( tpd, rd, : )/ PSI_POWER(1)


    !        PRINT*,"Before Calc_RHS_Terms"
    IF ( .FALSE. ) THEN
        CALL Calc_RHS_Terms( Source_Terms,                                      &
                             re, te, pe,                                        &
                             td, pd, tpd, rd,                                   &
                             CUR_R_LOCS, R_SQUARE, RSIN_SQUARE, COTAN_VAL,      &
                             SIN_SQUARE, CSC_SQUARE,                            &
                             PSI_POWER, ALPHAPSI_POWER,                         &
                             CUR_VAL_BETA, CUR_DRV_BETA,                        &
                             JCBN_BIGK_VALUE, JCBN_n_Array, JCBN_Kappa_Array    )

    ELSE
        CALL Calc_RHS_Terms_1D( Source_Terms,                                   &
                                re, te, pe,                                     &
                                td, pd, tpd, rd,                                &
                                CUR_R_LOCS, R_SQUARE, RSIN_SQUARE, COTAN_VAL,   &
                                SIN_SQUARE, CSC_SQUARE,                         &
                                PSI_POWER, ALPHAPSI_POWER,                      &
                                CUR_VAL_BETA, CUR_DRV_BETA,                     &
                                JCBN_BIGK_VALUE, JCBN_n_Array, JCBN_Kappa_Array )

    END IF




END DO ! pd loop
END DO  ! td loop
END DO  ! rd loop





END SUBROUTINE Calc_FP_Source_Terms




!+204+###########################################################################!
!                                                                                !
!                  CREATE_3D_RHS_VECTOR                                          !
!                                                                                !
!################################################################################!
SUBROUTINE Create_FP_Source_Vector( re, te, pe, DELTAR_OVERTWO, SIN_VAL )



INTEGER, INTENT(IN)                                                     ::  re, te, pe

REAL(KIND = idp), INTENT(IN)                                            ::  DELTAR_OVERTWO
REAL(KIND = idp), DIMENSION(1:NUM_TP_QUAD_POINTS), INTENT(IN)           ::  SIN_VAL

INTEGER                                                                 ::  pd, td, rd, tpd,     &
                                                                            l, m, d,        &
                                                                            lm_loc, u,ui

INTEGER                                                                 ::  Current_i_Location

COMPLEX(KIND = idp), DIMENSION(1:5)                                     ::  RHS_TMP
COMPLEX(KIND = idp)                                                     :: Test
COMPLEX(KIND = idp)                                                     ::  Common_Basis
REAL(KIND = idp)                                                        ::  Combined_Weights

COMPLEX(KIND = idp)                                                     ::  Inner, Middle


DO ui = 1,NUM_CFA_EQs
DO d = 0,DEGREE
DO lm_loc = 1,LM_LENGTH

    RHS_TMP = 0.0_idp


    DO rd = 1,NUM_R_QUAD_POINTS


        RHS_TMP(ui) =  RHS_TMP(ui)                                          &
                        + SUM( Source_Terms( :, rd, CFA_EQ_Map(ui) )        &
                                * Ylm_CC_Values( :, lm_loc, te, pe)         &
                                * TP_Int_Weights(:)                     )   &
                        * Lagrange_Poly_Table(d, rd, 0)                     &
                        * R_Int_Weights(rd)


    END DO  ! rd Loop

    Current_i_Location = FP_Vector_Map(re,d)
    FP_Source_Vector(Current_i_Location,lm_loc,ui)                &
        = FP_Source_Vector(Current_i_Location,lm_loc,ui)          &
        + RHS_TMP(ui)

END DO  ! lm_loc Loop
END DO  ! d Loop
END DO ! ui




END SUBROUTINE Create_FP_Source_Vector






!+701+###########################################################################!
!                                                                                !
!           Allocate_Master_Build_Variables                                      !
!                                                                                !
!################################################################################!
SUBROUTINE Allocate_FP_Source_Variables()



ALLOCATE( CUR_R_LOCS(1:NUM_R_QUAD_POINTS) )
ALLOCATE( CUR_T_LOCS(1:NUM_T_QUAD_POINTS) )
ALLOCATE( CUR_P_LOCS(1:NUM_P_QUAD_POINTS) )


ALLOCATE( R_SQUARE(1:NUM_R_QUAD_POINTS) )
ALLOCATE( R_CUBED(1:NUM_R_QUAD_POINTS) )

ALLOCATE( SIN_VAL_B( 1:NUM_TP_QUAD_POINTS ) )
ALLOCATE( SIN_VAL( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( SIN_SQUARE( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( TP_SIN_SQUARE( 1:NUM_TP_QUAD_POINTS ) )
ALLOCATE( COS_VAL( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( COS_SQUARE( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( CSC_VAL( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( CSC_SQUARE( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( COTAN_VAL( 1:NUM_T_QUAD_POINTS ) )

ALLOCATE( RSIN_SQUARE( 1:NUM_T_QUAD_POINTS, 1:NUM_R_QUAD_POINTS ) )

ALLOCATE( PHI_EXP( 1:NUM_P_QUAD_POINTS ) )
ALLOCATE( PHI_TWOEXP( 1:NUM_P_QUAD_POINTS ) )

ALLOCATE( R_Int_Weights( 1:NUM_R_QUAD_POINTS ) )
ALLOCATE( TP_Int_Weights( 1:NUM_TP_QUAD_POINTS) )

ALLOCATE( CUR_VAL_PSI(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS) )
ALLOCATE( CUR_VAL_ALPHAPSI(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS) )
ALLOCATE( CUR_VAL_BETA(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:3) )

ALLOCATE( CUR_DRV_PSI(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:3) )
ALLOCATE( CUR_DRV_ALPHAPSI(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:3) )
ALLOCATE( CUR_DRV_BETA(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:3, 1:3) )

ALLOCATE( Beta_DRV_Trace(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS) )

ALLOCATE( Source_Terms( 1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:5) )


END SUBROUTINE Allocate_FP_Source_Variables



!+702+###########################################################################!
!                                                                                !
!           Deallocate_Master_Build_Variables                                    !
!                                                                                !
!################################################################################!
SUBROUTINE Deallocate_FP_Source_Variables()

DEALLOCATE( CUR_R_LOCS )
DEALLOCATE( CUR_T_LOCS )
DEALLOCATE( CUR_P_LOCS )

DEALLOCATE( R_SQUARE )
DEALLOCATE( R_CUBED )

DEALLOCATE( SIN_VAL )
DEALLOCATE( SIN_SQUARE )
DEALLOCATE( TP_SIN_SQUARE )
DEALLOCATE( COS_VAL )
DEALLOCATE( COS_SQUARE )
DEALLOCATE( CSC_VAL )
DEALLOCATE( CSC_SQUARE )
DEALLOCATE( COTAN_VAL )

DEALLOCATE( RSIN_SQUARE )

DEALLOCATE( PHI_EXP )
DEALLOCATE( PHI_TWOEXP )

DEALLOCATE( R_Int_Weights )
DEALLOCATE( TP_Int_Weights )

DEALLOCATE( CUR_VAL_PSI )
DEALLOCATE( CUR_VAL_ALPHAPSI )
DEALLOCATE( CUR_VAL_BETA )

DEALLOCATE( CUR_DRV_PSI )
DEALLOCATE( CUR_DRV_ALPHAPSI )
DEALLOCATE( CUR_DRV_BETA )

DEALLOCATE( Beta_DRV_Trace )

DEALLOCATE( Source_Terms )


END SUBROUTINE Deallocate_FP_Source_Variables

END MODULE FP_Source_Beta
