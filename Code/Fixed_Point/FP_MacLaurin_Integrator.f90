   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE FP_MacLaurin_Integrator_Module                                               !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   FP_MacLaurin_Integrator                                             !##!
!##!                                                                                !##!
!##!    +201+   Calc_FPML_Current_Values                                            !##!
!##!    +202+   Calc_FPML_Source_Terms                                              !##!
!##!    +203+   Create_FPML_Source_Vector                                           !##!
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

USE Units_Module, &
            ONLY :  Grav_Constant_G,        &
                    Gram,                  &
                    Centimeter,             &
                    C_Square

USE Poseidon_Numbers_Module, &
            ONLY :  pi

USE Poseidon_Parameters, &
            ONLY :  Degree,                     &
                    L_Limit,                    &
                    Verbose_Flag


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
                    Ylm_dt_Values,              &
                    Ylm_dp_Values,              &
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
                    CFA_EQ_Map,                 &
                    CFA_EQ_Flags


USE Functions_Jacobian, &
            ONLY :  JCBN_kappa_FUNCTION_3D_ALL,     &
                    JCBN_BIGK_FUNCTION,         &
                    JCBN_Kappa_Array_3D

USE Poseidon_IO_Module, &
            ONLY :  Clock_In


USE FP_Functions_Mapping, &
            ONLY :  FP_FEM_Node_Map,            &
                    FP_Beta_Array_Map,          &
                    FP_Array_Map,               &
                    FP_LM_Map

USE Functions_Mapping, &
            ONLY :  Map_To_X_Space,         &
                    Map_From_X_Space


IMPLICIT NONE

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)             :: CUR_R_LOCS
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CUR_T_LOCS
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CUR_P_LOCS

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)             :: R_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)             :: R_CUBED

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


REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)             :: R_Int_Weights
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
!                  FP_MacLaurin_Integrator           !
!                                                                                !
!################################################################################!
SUBROUTINE FP_MacLaurin_Integrator( re  )
                    
INTEGER, INTENT(IN)                                                 :: re

REAL(idp),DIMENSION(1:NUM_T_QUAD_POINTS)                            ::  DELTAR_OVERTWO
REAL(idp),DIMENSION(1:NUM_T_QUAD_POINTS)                            ::  TWOOVER_DELTAR

INTEGER, DIMENSION(1:NUM_T_QUAD_POINTS)                             ::  Disc_rd


INTEGER                                                             ::  te, pe
INTEGER                                                             ::  tpd, td, pd, rd

REAL(idp)                                                           ::  Deltap_overtwo
REAL(idp)                                                           ::  DeltaT_OverTwo

CALL Allocate_FPML_Source_Variables()





DO pe = 0,NUM_P_ELEMENTS-1
DO te = 0,NUM_T_ELEMENTS-1


    ! Need to id location of discontinuity.
    Disc_rd(:) = FINDLOC(Block_Source_E(:, :, 1, re, te, pe),0.0_idp,1)
    PRINT*,Disc_rd
    PRINT*,Block_Source_E(:, :, 1, re, te, pe)
    


!    DELTAR_OVERTWO = 0.5_idp *(rlocs(re + 1) - rlocs(re))
!    TWOOVER_DELTAR = 1.0_idp/deltar_overtwo
!    CUR_R_LOCS(:) = deltar_overtwo * (INT_R_LOCATIONS(:)+1.0_idp) + rlocs(re)
!
!
!    R_SQUARE(:) = CUR_R_LOCS(:)*CUR_R_LOCS(:)
!    R_CUBED(:) = R_SQUARE(:)*CUR_R_LOCS(:)



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

!    DO rd = 1,NUM_R_QUAD_POINTS
!
!        RSIN_SQUARE(:,rd) = R_SQUARE(rd)*SIN_SQUARE(:)
!
!    END DO




    
        deltap_overtwo = 0.5_idp * (plocs(pe + 1)-plocs(pe))
        CUR_P_LOCS(:) = deltap_overtwo * (INT_P_LOCATIONS+1.0_idp) + plocs(pe)



!        CALL Calc_FPML_Current_Values(  re, te, pe,     &
!                                        DELTAR_OVERTWO, &
!                                        DELTAT_OVERTWO, &
!                                        DELTAP_OVERTWO  )
!        CALL Calc_FPML_Source_Terms(    re, te, pe )
!        CALL Create_FPML_Source_Vector( re, te, pe,       &
!                                        DELTAR_OVERTWO,   &
!                                        SIN_VAL_B         )



END DO ! te Loop
END DO ! pe Loop






CALL DEAllocate_FPML_Source_Variables()

PRINT*,"STOPing in FP_MacLaurin_Integrator"
STOP

END SUBROUTINE FP_MacLaurin_Integrator














!+201+###########################################################################!
!                                                                                !
!                  Calc_FPML_Current_Values          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_FPML_Current_Values( re, te, pe,                                 &
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
!R_Int_Weights(:) = DELTAR_OVERTWO * R_SQUARE(:) * INT_R_WEIGHTS(:)

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

    
    DO ui = 1,5
    DO d = 0,DEGREE
        Here = FP_FEM_Node_Map(re,d)

    


        TMP_U_Value(ui)         = TMP_U_Value(ui)                           &
                                + SUM( FP_Coeff_Vector( Here, :, ui )       &
                                * Ylm_Values( :, tpd, te, pe )       )      &
                                * Lagrange_Poly_Table( d, rd, 0 )


        TMP_U_R_DRV_Value(ui)   = TMP_U_R_DRV_Value(ui)                     &
                                + SUM( FP_Coeff_Vector( Here, :, ui )       &
                                * Ylm_Values( :, tpd, te, pe )       )      &
                                * Lagrange_Poly_Table( d, rd, 1 )           &
                                / DELTAR_OVERTWO


        TMP_U_T_DRV_Value(ui)   = TMP_U_T_DRV_Value(ui)                     &
                                + SUM( FP_Coeff_Vector( Here, :, ui )       &
                                * Ylm_dt_Values( :, tpd, te, pe)     )      &
                                * Lagrange_Poly_Table( d, rd, 0)

        TMP_U_P_DRV_Value(ui)   = TMP_U_P_DRV_Value(ui)                     &
                                + SUM( FP_Coeff_Vector( Here, :, ui )       &
                                * Ylm_dp_Values( :, tpd, te, pe)     )      &
                                * Lagrange_Poly_Table( d, rd, 0)




    END DO  ! d
    END DO  ! ui

    CUR_VAL_PSI( tpd, rd )         = REAL(Tmp_U_Value(1), KIND = idp)
    CUR_DRV_PSI( tpd, rd, 1 )      = REAL(Tmp_U_R_DRV_Value(1), KIND = idp)
    CUR_DRV_PSI( tpd, rd, 2 )      = REAL(Tmp_U_T_DRV_Value(1), KIND = idp)
    CUR_DRV_PSI( tpd, rd, 3 )      = REAL(Tmp_U_P_DRV_Value(1), KIND = idp)


    CUR_VAL_ALPHAPSI( tpd, rd )    = REAL(Tmp_U_Value(2), KIND = idp)
    CUR_DRV_ALPHAPSI( tpd, rd, 1 ) = REAL(Tmp_U_R_DRV_Value(2), KIND = idp)
    CUR_DRV_ALPHAPSI( tpd, rd, 2 ) = REAL(Tmp_U_T_DRV_Value(2), KIND = idp)
    CUR_DRV_ALPHAPSI( tpd, rd, 3 ) = REAL(Tmp_U_P_DRV_Value(2), KIND = idp)


    CUR_VAL_BETA( tpd, rd, 1 )     = REAL(Tmp_U_Value(3), KIND = idp)
    CUR_VAL_BETA( tpd, rd, 2 )     = REAL(Tmp_U_Value(4), KIND = idp)
    CUR_VAL_BETA( tpd, rd, 3 )     = REAL(Tmp_U_Value(5), KIND = idp)


    CUR_DRV_BETA( tpd, rd, 1, 1 )  = REAL(Tmp_U_R_DRV_Value(3), KIND = idp)
    CUR_DRV_BETA( tpd, rd, 2, 1 )  = REAL(Tmp_U_R_DRV_Value(4), KIND = idp)
    CUR_DRV_BETA( tpd, rd, 3, 1 )  = REAL(Tmp_U_R_DRV_Value(5), KIND = idp)

    CUR_DRV_BETA( tpd, rd, 1, 2 )  = REAL(Tmp_U_T_DRV_Value(3), KIND = idp)
    CUR_DRV_BETA( tpd, rd, 2, 2 )  = REAL(Tmp_U_T_DRV_Value(4), KIND = idp)
    CUR_DRV_BETA( tpd, rd, 3, 2 )  = REAL(Tmp_U_T_DRV_Value(5), KIND = idp)

    CUR_DRV_BETA( tpd, rd, 1, 3 )  = REAL(Tmp_U_P_DRV_Value(3), KIND = idp)
    CUR_DRV_BETA( tpd, rd, 2, 3 )  = REAL(Tmp_U_P_DRV_Value(4), KIND = idp)
    CUR_DRV_BETA( tpd, rd, 3, 3 )  = REAL(Tmp_U_P_DRV_Value(5), KIND = idp)

    Beta_DRV_Trace( tpd, rd )      = CUR_DRV_BETA( tpd, rd, 1, 1 )              &
                                   + CUR_DRV_BETA( tpd, rd, 2, 2 )              &
                                   + CUR_DRV_BETA( tpd, rd, 3, 3 )


END DO ! tpd
END DO ! rd

END SUBROUTINE Calc_FPML_Current_Values


















!+701+###########################################################################!
!                                                                                !
!           Allocate_Master_Build_Variables                                      !
!                                                                                !
!################################################################################!
SUBROUTINE Allocate_FPML_Source_Variables()



ALLOCATE( CUR_R_LOCS(1:NUM_T_QUAD_POINTS,1:NUM_R_QUAD_POINTS) )
ALLOCATE( CUR_T_LOCS(1:NUM_T_QUAD_POINTS) )
ALLOCATE( CUR_P_LOCS(1:NUM_P_QUAD_POINTS) )


ALLOCATE( R_SQUARE(1:NUM_T_QUAD_POINTS,1:NUM_R_QUAD_POINTS) )
ALLOCATE( R_CUBED(1:NUM_T_QUAD_POINTS,1:NUM_R_QUAD_POINTS) )

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


ALLOCATE( R_Int_Weights( 1:NUM_T_QUAD_POINTS,1:NUM_R_QUAD_POINTS ) )
ALLOCATE( TP_Int_Weights( 1:NUM_TP_QUAD_POINTS) )

ALLOCATE( CUR_VAL_PSI(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS) )
ALLOCATE( CUR_VAL_ALPHAPSI(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS) )
ALLOCATE( CUR_VAL_BETA(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:3) )

ALLOCATE( CUR_DRV_PSI(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:3) )
ALLOCATE( CUR_DRV_ALPHAPSI(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:3) )
ALLOCATE( CUR_DRV_BETA(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:3, 1:3) )

ALLOCATE( Beta_DRV_Trace(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS) )

ALLOCATE( Source_Terms( 1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:5) )


END SUBROUTINE Allocate_FPML_Source_Variables



!+702+###########################################################################!
!                                                                                !
!           Deallocate_Master_Build_Variables                                    !
!                                                                                !
!################################################################################!
SUBROUTINE Deallocate_FPML_Source_Variables()

DEALLOCATE( CUR_R_LOCS )
DEALLOCATE( CUR_T_LOCS )
DEALLOCATE( CUR_P_LOCS )

DEALLOCATE( R_SQUARE )
DEALLOCATE( R_CUBED )

DEALLOCATE( SIN_VAL_B )
DEALLOCATE( SIN_VAL )
DEALLOCATE( SIN_SQUARE )
DEALLOCATE( TP_SIN_SQUARE )
DEALLOCATE( COS_VAL )
DEALLOCATE( COS_SQUARE )
DEALLOCATE( CSC_VAL )
DEALLOCATE( CSC_SQUARE )
DEALLOCATE( COTAN_VAL )

DEALLOCATE( RSIN_SQUARE )


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


END SUBROUTINE Deallocate_FPML_Source_Variables



END MODULE FP_MacLaurin_Integrator_Module
