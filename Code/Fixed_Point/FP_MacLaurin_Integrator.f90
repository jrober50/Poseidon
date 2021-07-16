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
            ONLY :  pi, twopi

USE Units_Module, &
            ONLY :  GR_Source_Scalar

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
                    Ylm_CC_Values

USE Variables_Derived, &
            ONLY :  NUM_R_NODES,                &
                    LM_LENGTH

USE Variables_FP, &
            ONLY :  FP_Coeff_Vector,            &
                    FP_Source_Vector,           &
                    FP_Source_Vector_Beta,      &
                    CFA_EQ_Map,                 &
                    CFA_EQ_Flags


USE Functions_Jacobian, &
            ONLY :  JCBN_kappa_FUNCTION_3D_ALL,     &
                    JCBN_BIGK_FUNCTION,         &
                    JCBN_Kappa_Array_3D

USE Functions_Quadrature, &
            ONLY :  Initialize_LGL_Quadrature_Locations

USE Functions_Math, &
            ONLY :  Lagrange_Poly,          &
                    Lagrange_Poly_Deriv

USE Poseidon_IO_Module, &
            ONLY :  Clock_In


USE FP_Functions_Mapping, &
            ONLY :  FP_FEM_Node_Map,            &
                    FP_Beta_Array_Map,          &
                    FP_Array_Map,               &
                    FP_LM_Map,                  &
                    FP_tpd_Map

USE Functions_Mapping, &
            ONLY :  Map_To_X_Space,         &
                    Map_From_X_Space



IMPLICIT NONE

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)             :: CUR_R_LOCS
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CUR_T_LOCS
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CUR_P_LOCS
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: xloc

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

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:,:)       :: LagPoly_Table

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

REAL(idp),DIMENSION(0:NUM_R_QUAD_POINTS)                            ::  Old_R_Locs

CALL Allocate_FPML_Source_Variables()





DO pe = 0,NUM_P_ELEMENTS-1
DO te = 0,NUM_T_ELEMENTS-1


    ! Need to id location of discontinuity.
    Disc_rd(:) = FINDLOC(Block_Source_E(:, :, 1, re, te, pe),0.0_idp,1)

    

    IF (  ANY( Disc_rd .NE. 0 ) .AND. ANY( Disc_rd .NE. 1) ) THEN
        Old_R_Locs(0) = rlocs(re+1)
        Old_R_Locs(1:Num_R_Quad_Points) = 0.5_idp *(rlocs(re + 1) - rlocs(re))       &
                                        * (INT_R_LOCATIONS(:)+1.0_idp) + rlocs(re)

        DO td = 1,Num_T_Quad_Points
            DeltaR_OverTwo(td) = 0.5_idp*(Old_R_Locs(Disc_rd(td)) - rlocs(re))
            

            TwoOver_DeltaR(td) = 1.0_idp/DeltaR_OverTwo(td)
            CUR_R_LOCS(:,td) = deltar_overtwo(td)* (INT_R_LOCATIONS(:)+1.0_idp) + rlocs(re)
            xloc(td) = Map_To_X_Space(rlocs(re),rlocs(re+1),Cur_R_Locs(Num_R_Quad_Points,td))
            R_SQUARE(:,td) = CUR_R_LOCS(:,td)*CUR_R_LOCS(:,td)
            R_CUBED(:,td) = R_SQUARE(:,td)*CUR_R_LOCS(:,td)


        END DO ! td

    ELSE

        DELTAR_OVERTWO = 0.5_idp *(rlocs(re + 1) - rlocs(re))
        TWOOVER_DELTAR = 1.0_idp/deltar_overtwo
        DO td = 1,Num_T_Quad_Points
            xloc(td) = 1.0_idp
            CUR_R_LOCS(:,td) = deltar_overtwo * (INT_R_LOCATIONS(:)+1.0_idp) + rlocs(re)
            R_SQUARE(:,td) = CUR_R_LOCS(:,td)*CUR_R_LOCS(:,td)
            R_CUBED(:,td) = R_SQUARE(:,td)*CUR_R_LOCS(:,td)
        END DO ! td

    END IF
    
    




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


    DO rd = 1,NUM_R_QUAD_POINTS

        RSIN_SQUARE(:,rd) = R_SQUARE(rd,:)*SIN_SQUARE(:)

    END DO


    deltap_overtwo = 0.5_idp * (plocs(pe + 1)-plocs(pe))
    CUR_P_LOCS(:) = deltap_overtwo * (INT_P_LOCATIONS+1.0_idp) + plocs(pe)


    CALL Calc_FPML_Current_Values(  re, te, pe,     &
                                    DELTAR_OVERTWO, &
                                    DELTAT_OVERTWO, &
                                    DELTAP_OVERTWO  )
    
    CALL Calc_FPML_Source_Terms(    re, te, pe )
    
    CALL Create_FPML_Source_Vector( re, te, pe,       &
                                    DELTAR_OVERTWO,   &
                                    SIN_VAL_B         )
    


END DO ! te Loop
END DO ! pe Loop






CALL Deallocate_FPML_Source_Variables()

!PRINT*,"STOPing in FP_MacLaurin_Integrator"
!STOP

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

REAL(KIND = idp), INTENT(IN), DIMENSION(1:Num_T_Quad_Points)    ::  DeltaR_OverTwo
REAL(KIND = idp), INTENT(IN)                                    ::  DELTAT_OVERTWO,     &
                                                                    DELTAP_OVERTWO



COMPLEX(KIND = idp), DIMENSION(1:5)                             ::  Tmp_U_Value,        &
                                                                    Tmp_U_R_DRV_Value,  &
                                                                    Tmp_U_T_DRV_Value,  &
                                                                    Tmp_U_P_DRV_Value



INTEGER                                                         ::  tpd, td, pd, rd,    &
                                                                    lm, d, Here, ui

REAL(KIND = idp), DIMENSION(1:Num_R_Quad_Points)    ::  Local_Locs
REAL(KIND = idp), DIMENSION(0:Degree)               ::  FEM_Locs
REAL(KIND = idp), DIMENSION(0:Degree)               ::  Lagrange_Poly_Values
REAL(KIND = idp), DIMENSION(0:Degree)               ::  Lagrange_DRV_Values

                          !                                                 !
                         !!                                                 !!
                        !!!          Initialize Local Quadratures           !!!
                         !!                                                 !!
                          !                                                 !



FEM_Locs = Initialize_LGL_Quadrature_Locations(Degree)


DO td = 1,NUM_T_QUAD_POINTS
    Local_Locs = Map_From_X_Space(-1.0_idp, xloc(td), Int_R_Locations)


    R_Int_Weights(:,td) = DELTAR_OVERTWO(td) * R_SQUARE(:,td) * INT_R_WEIGHTS(:)

    DO rd = 1,NUM_R_QUAD_POINTS

        Lagrange_Poly_Values = Lagrange_Poly(Local_Locs(rd), Degree, FEM_Locs)
        Lagrange_DRV_Values  = Lagrange_Poly_Deriv(Local_Locs(rd), Degree, FEM_Locs)
    

        LagPoly_Table(:, td, rd, 0) = Lagrange_Poly_Values
        LagPoly_Table(:, td, rd, 1) = Lagrange_DRV_Values

    END DO
END DO






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
DO td = 1,NUM_T_QUAD_POINTS
DO pd = 1,NUM_P_QUAD_POINTS

    tpd = FP_tpd_Map(td,pd)

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
                                * LagPoly_Table( d, td, rd, 0 )


        TMP_U_R_DRV_Value(ui)   = TMP_U_R_DRV_Value(ui)                     &
                                + SUM( FP_Coeff_Vector( Here, :, ui )       &
                                * Ylm_Values( :, tpd, te, pe )       )      &
                                * LagPoly_Table( d, td, rd, 1 )           &
                                / DELTAR_OVERTWO(td)


        TMP_U_T_DRV_Value(ui)   = TMP_U_T_DRV_Value(ui)                     &
                                + SUM( FP_Coeff_Vector( Here, :, ui )       &
                                * Ylm_dt_Values( :, tpd, te, pe)     )      &
                                * LagPoly_Table( d, td, rd, 0)

        TMP_U_P_DRV_Value(ui)   = TMP_U_P_DRV_Value(ui)                     &
                                + SUM( FP_Coeff_Vector( Here, :, ui )       &
                                * Ylm_dp_Values( :, tpd, te, pe)     )      &
                                * LagPoly_Table( d, td, rd, 0)




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

END DO ! pd
END DO ! td
END DO ! rd


END SUBROUTINE Calc_FPML_Current_Values









!+202+###########################################################################!
!                                                                                !
!                  Calc_FP_Source_Terms          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_FPML_Source_Terms( re, te, pe )

INTEGER, INTENT(IN)                                                     ::  re, te, pe

REAL(KIND = idp)                                                        ::  REUSED_VALUE

INTEGER                                                                 ::  pd, td, rd,     &
                                                                            i, tpd


REAL(KIND = idp), DIMENSION(1:11)                                       ::  PSI_POWER
REAL(KIND = idp), DIMENSION(1:4)                                        ::  ALPHAPSI_POWER


REAL(KIND = idp), DIMENSION(1:3,1:3)                                    ::  Kappa_Array
REAL(KIND = idp), DIMENSION(1:3)                                        ::  n_Array

REAL(KIND = idp)                                                        ::  BigK_Value


Kappa_Array = 0.0_idp
n_Array     = 0.0_idp
BigK_Value  = 0.0_idp


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
!    PRINT*,"Before BigK_Value"
    BigK_Value = JCBN_BIGK_FUNCTION( rd, tpd,                                                        &
                                     NUM_R_QUAD_POINTS, NUM_TP_QUAD_POINTS,                          &
                                     CUR_VAL_BETA, CUR_DRV_BETA,                                     &
                                     CUR_R_LOCS(rd,td), R_SQUARE(rd,td), SIN_SQUARE(td), CSC_SQUARE(td),   &
                                     RSIN_SQUARE(td, rd), COTAN_VAL(td)                              )



!    PRINT*,"Before  JCBN_kappa_Array"
    Kappa_Array = JCBN_kappa_FUNCTION_3D_ALL(  rd, tpd,                                        &
                                               NUM_R_QUAD_POINTS, NUM_TP_QUAD_POINTS,          &
                                               CUR_R_LOCS(rd,td), R_SQUARE(rd,td), R_CUBED(rd,td),      &
                                               RSIN_SQUARE(td, rd),                            &
                                               COTAN_VAL(td),                                   &
                                               CUR_VAL_BETA, CUR_DRV_BETA                      )

!    PRINT*,"Before n_Array"
    n_Array(:) = CUR_DRV_ALPHAPSI( tpd, rd, : ) / ALPHAPSI_POWER(1)   &
                        - 7.0_idp * CUR_DRV_PSI( tpd, rd, : )/ PSI_POWER(1)



!    PRINT*,"Before Calc_Source_Terms"
    CALL Calc_Source_Terms( Source_Terms,                                      &
                            re, te, pe,                                        &
                            td, pd, tpd, rd,                                   &
                            CUR_R_LOCS(:,td), R_SQUARE(:,td), RSIN_SQUARE, COTAN_VAL,      &
                            SIN_SQUARE, CSC_SQUARE,                            &
                            PSI_POWER, ALPHAPSI_POWER,                         &
                            CUR_VAL_BETA, CUR_DRV_BETA,                        &
                            BigK_Value, n_Array, Kappa_Array    )


    

END DO  ! pd loop
END DO  ! td loop
END DO  ! rd loop




END SUBROUTINE Calc_FPML_Source_Terms








!+301+###########################################################################!
!                                                                                !
!           Calc_Source_Terms                                                       !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Source_Terms( Source_Terms,                                      &
                              re, te, pe,                                        &
                              td, pd, tpd, rd,                                   &
                              CUR_R_LOCS, R_SQUARE, RSIN_SQUARE, COTAN_VAL,      &
                              SIN_SQUARE, CSC_SQUARE,                            &
                              PSI_POWER, ALPHAPSI_POWER,                         &
                              CUR_VAL_BETA, CUR_DRV_BETA,                        &
                              BigK_Value, n_Array, Kappa_Array    )

REAL(KIND = idp), INTENT(INOUT), DIMENSION( 1:NUM_TP_QUAD_POINTS,   &
                                            1:NUM_R_QUAD_POINTS,    &
                                            1:5                     )   ::  Source_Terms

INTEGER, INTENT(IN)                                                     ::  re, te, pe
INTEGER, INTENT(IN)                                                     ::  td, pd, tpd, rd

REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_R_QUAD_POINTS)            ::  CUR_R_LOCS

REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_R_QUAD_POINTS)            ::  R_SQUARE

REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_T_QUAD_POINTS,    &
                                        1:NUM_R_QUAD_POINTS     )       ::  RSIN_SQUARE

REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_T_QUAD_POINTS)            ::  SIN_SQUARE,         &
                                                                            CSC_SQUARE,         &
                                                                            COTAN_VAL

REAL(KIND = idp), INTENT(IN), DIMENSION(1:11)                           ::  PSI_POWER
REAL(KIND = idp), INTENT(IN), DIMENSION(1:4)                            ::  ALPHAPSI_POWER

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:NUM_TP_QUAD_POINTS,  &
                                         1:NUM_R_QUAD_POINTS,   &
                                         1:3                    )       ::  CUR_VAL_BETA

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:NUM_TP_QUAD_POINTS,  &
                                         1:NUM_R_QUAD_POINTS,   &
                                         1:3, 1:3               )       ::  CUR_DRV_BETA



REAL(KIND = idp), INTENT(IN)                                            ::  BigK_VALUE
REAL(KIND = idp), INTENT(IN), DIMENSION(1:3)                            ::  n_ARRAY
REAL(KIND = idp), INTENT(IN), DIMENSION(1:3,1:3)                        ::  Kappa_Array

REAL(KIND = idp)                                                        ::  Beta_Source_Prefix



!PRINT*,"Calc_Source_Terms has been altered"
Source_Terms(tpd, rd, :) = 0.0_idp
Source_Terms(tpd, rd, 1) = - TwoPi                                      &
                           * GR_Source_Scalar                          &
                           * Block_Source_E(1, td, pd, re, te, pe)


!Source_Terms(tpd, rd, 1) = - TwoPi                                      &
!                            * GR_Source_Scalar                          &
!                            * Block_Source_E(rd, td, pd, re, te, pe)    &
!                            * PSI_POWER(5)                              &
!                         - PSI_POWER(7)                                 &
!                            / ( 16.0_idp * ALPHAPSI_POWER(2) )          &
!                            * BigK_Value
     


!PRINT*,Block_Source_E(rd, td, pd, re, te, pe) ,PSI_POWER(1),ALPHAPSI_POWER(1)



!Source_Terms(tpd, rd, 2) = TwoPi                                                        &
!                            * ALPHAPSI_POWER(1)                                         &
!                            * PSI_POWER(4)                                              &
!                            * GR_Source_Scalar                                          &
!                                * ( Block_Source_E(rd, td, pd, re, te, pe)              &
!                                    + 2.0_idp                                           &
!                                        * Block_Source_S(rd, td, pd, re, te, pe)  )     &
!                         + 7.0_idp                                                      &
!                            * PSI_POWER(6)                                              &
!                            / ( 16.0_idp * ALPHAPSI_POWER(1) )                          &
!                            * BigK_Value
!
!
!
!
!
!Beta_Source_Prefix = 16.0_idp * pi * ALPHAPSI_POWER(1) * PSI_POWER(3) * GR_Source_Scalar
!
!Source_Terms(tpd, rd, 3) = Beta_Source_Prefix * Block_Source_Si(rd, td, pd, re, te, pe, 1)      &
!                         + Kappa_Array(1,1) * n_Array(1)                                        &
!                         + Kappa_Array(2,1) * n_Array(2)                                        &
!                         + Kappa_Array(3,1) * n_Array(3)
!
!
!Source_Terms(tpd, rd, 4) = Beta_Source_Prefix * Block_Source_Si(rd, td, pd, re, te, pe, 2)      &
!                         + Kappa_Array(1,2) * n_Array(1)                                        &
!                         + Kappa_Array(2,2) * n_Array(2)                                        &
!                         + Kappa_Array(3,2) * n_Array(3)
!
!
!Source_Terms(tpd, rd, 5) = Beta_Source_Prefix * Block_Source_Si(rd, td, pd, re, te, pe, 3)      &
!                         + Kappa_Array(1,3) * n_Array(1)                                        &
!                         + Kappa_Array(2,3) * n_Array(2)                                        &
!                         + Kappa_Array(3,3) * n_Array(3)




END SUBROUTINE Calc_Source_Terms
















!+204+###########################################################################!
!                                                                                !
!                  CREATE_3D_RHS_VECTOR                                          !
!                                                                                !
!################################################################################!
SUBROUTINE Create_FPML_Source_Vector( re, te, pe, DELTAR_OVERTWO, SIN_VAL )



INTEGER, INTENT(IN)                                                     ::  re, te, pe

REAL(KIND = idp), DIMENSION(1:NUM_T_QUAD_POINTS), INTENT(IN)            ::  DELTAR_OVERTWO
REAL(KIND = idp), DIMENSION(1:NUM_TP_QUAD_POINTS), INTENT(IN)           ::  SIN_VAL

INTEGER                                                                 ::  pd, td, rd, tpd,     &
                                                                            l, m, d,        &
                                                                            lm_loc, u,ui

INTEGER                                                                 ::  Current_i_Location

COMPLEX(KIND = idp), DIMENSION(1:5)                                     ::  RHS_TMP
COMPLEX(KIND = idp)                                                     ::  Test
COMPLEX(KIND = idp)                                                     ::  Common_Basis
REAL(KIND = idp)                                                        ::  Combined_Weights

COMPLEX(KIND = idp)                                                     ::  Inner, Middle




DO ui = 1,2
    IF ( CFA_EQ_Flags(ui) == 1 ) THEN

        
        DO lm_loc = 1,LM_LENGTH
!
!        DO l = 0,L_Limit
!        lm_loc = FP_LM_Map(l,0)
        DO d = 0,DEGREE

            RHS_TMP = 0.0_idp

            DO rd = 1,NUM_R_QUAD_POINTS
            DO td = 1,Num_T_Quad_Points
            DO pd = 1,Num_P_Quad_Points

                tpd = FP_tpd_Map(td,pd)

                RHS_TMP(ui) =  RHS_TMP(ui)                                      &
                            + Source_Terms( tpd, rd, ui )                       &
                                * Ylm_CC_Values( tpd, lm_loc, te, pe)           &
                                * TP_Int_Weights(tpd)                           &
                                * LagPoly_Table(d, td,rd, 0)              &
                                * R_Int_Weights(rd,td)


            END DO
            END DO
            END DO  ! rd Loop
            

            Current_i_Location = FP_FEM_Node_Map(re,d)
            FP_Source_Vector(Current_i_Location,lm_loc,ui)                &
                = FP_Source_Vector(Current_i_Location,lm_loc,ui)          &
                + RHS_TMP(ui)

        END DO  ! d Loop
        END DO  ! lm_loc Loop
    END IF


END DO





DO ui = 3,5
    IF( CFA_EQ_Flags(ui) == 1 ) THEN

        DO lm_loc = 1,LM_LENGTH
        DO d = 0,DEGREE
        
            RHS_TMP = 0.0_idp
            DO rd = 1,NUM_R_QUAD_POINTS
            DO td = 1,Num_T_Quad_Points
            DO pd = 1,Num_P_Quad_Points

                tpd = FP_tpd_Map(td,pd)

                RHS_TMP(ui) =  RHS_TMP(ui)                                  &
                            + Source_Terms( tpd, rd, ui )                   &
                                * Ylm_CC_Values( tpd, lm_loc, te, pe)       &
                                * TP_Int_Weights(tpd)                       &
                                * LagPoly_Table(d, td,rd, 0)          &
                                * R_Int_Weights(rd,td)


            END DO
            END DO
            END DO  ! rd Loop

            Current_i_Location = FP_Beta_Array_Map(re,d,ui-2,lm_loc)

            FP_Source_Vector_Beta(Current_i_Location)                &
                = FP_Source_Vector_Beta(Current_i_Location)          &
                + RHS_TMP(ui)

        END DO  ! d Loop
        END DO  ! lm_loc Loop

    END IF

END DO ! ui





END SUBROUTINE Create_FPML_Source_Vector




















!+701+###########################################################################!
!                                                                                !
!           Allocate_Master_Build_Variables                                      !
!                                                                                !
!################################################################################!
SUBROUTINE Allocate_FPML_Source_Variables()



ALLOCATE( CUR_R_LOCS(1:NUM_T_QUAD_POINTS,1:NUM_R_QUAD_POINTS) )
ALLOCATE( CUR_T_LOCS(1:NUM_T_QUAD_POINTS) )
ALLOCATE( CUR_P_LOCS(1:NUM_P_QUAD_POINTS) )
ALLOCATE( xloc(1:NUM_T_QUAD_POINTS) )

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


ALLOCATE( LagPoly_Table( 0:Degree, 1:Num_T_Quad_Points, 1:Num_R_Quad_Points, 0:2) )

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
DEALLOCATE( xloc )

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

DEALLOCATE( LagPoly_Table )

END SUBROUTINE Deallocate_FPML_Source_Variables



END MODULE FP_MacLaurin_Integrator_Module
