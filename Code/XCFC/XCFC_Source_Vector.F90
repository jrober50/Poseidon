   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE XCFC_Source_Vector_Module                                             !##!
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
                    FP_Coeff_Vector_Orig,       &
                    FP_Source_Vector,           &
                    FP_Source_Vector_Beta,      &
                    FP_Source_Vector_X,         &
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
                    FP_Array_Map,               &
                    FP_LM_Map,                  &
                    FP_tpd_Map


USE MPI




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
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: CUR_VAL_BETA

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: CUR_DRV_PSI
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: CUR_DRV_ALPHAPSI
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:,:)       :: CUR_DRV_BETA

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)           :: Beta_DRV_Trace

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: StaredSource
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: SourceTerm


CONTAINS




!+101+###########################################################################!
!                                                                                !
!           XCFC_Calc_X_Source                                                   !
!                                                                                !
!################################################################################!
SUBROUTINE XCFC_Calc_X_Source()



INTEGER                                                     ::  re, te, pe,     &
                                                                d, lm,          &
                                                                rd, tpd, td, pd

REAL(KIND = idp)                                            ::  deltar_overtwo,     &
                                                                deltat_overtwo,     &
                                                                deltap_overtwo


FP_Source_Vector_X = 0.0_idp

DO re = 0,NUM_R_ELEMENTS-1

    DELTAR_OVERTWO = 0.5_idp *(rlocs(re + 1) - rlocs(re))
    CUR_R_LOCS(:) = deltar_overtwo * (INT_R_LOCATIONS(:)+1.0_idp) + rlocs(re)


    R_SQUARE(:) = CUR_R_LOCS(:)*CUR_R_LOCS(:)



    DO pe = 0,NUM_P_ELEMENTS-1
        deltap_overtwo = 0.5_idp * (plocs(pe + 1)-plocs(pe))


        DO te = 0,NUM_T_ELEMENTS-1

            deltat_overtwo = 0.5_idp*(tlocs(te + 1) - tlocs(te))
            CUR_T_LOCS(:) = deltat_overtwo * (INT_T_LOCATIONS(:)+1.0_idp) + tlocs(te)

            SIN_VAL(:) = DSIN(CUR_T_LOCS(:))
            SIN_SQUARE(:) = SIN_VAL(:)*SIN_VAL(:)
            DO td = 1,NUM_T_QUAD_POINTS
            DO pd = 1,NUM_P_QUAD_POINTS
                tpd = (td-1)*NUM_P_QUAD_POINTS + pd
                TP_Sin_Val(tpd) =   SIN_VAL(td)
                TP_SIN_SQUARE(tpd) = SIN_SQUARE(td)
            END DO
            END DO


            CALL Calc_XCFC_X_CurVals( re, te, pe,       &
                                      DELTAR_OVERTWO,   &
                                      DELTAT_OVERTWO,   &
                                      DELTAP_OVERTWO    )

            CALL Create_XCFC_X_Vector( re, te, pe       )



        END DO ! te Loop
    END DO ! pe Loop
END DO ! re Loop



END SUBROUTINE XCFC_Calc_X_Source



!+202+###########################################################################!
!                                                                                !
!                  Calc_Current_Values          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_XCFC_X_CurVals( re, te, pe,                                  &
                                DELTAR_OVERTWO,                             &
                                DELTAT_OVERTWO,                             &
                                DELTAP_OVERTWO                              )

INTEGER, INTENT(IN)                                             ::  re, te, pe


REAL(KIND = idp), INTENT(IN)                                    ::  DELTAR_OVERTWO,     &
                                                                    DELTAT_OVERTWO,     &
                                                                    DELTAP_OVERTWO


COMPLEX(KIND = idp), DIMENSION(1:5)                             ::  Tmp_U_Value

INTEGER                                                         ::  tpd, td, pd, rd
INTEGER                                                         ::  d, Here
INTEGER                                                         ::  ui, s



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

    tpd = FP_tpd_Map(td,pd)
    TP_Int_Weights( tpd ) = SIN_VAL(td)                           &
                          * DELTAT_OVERTWO * INT_T_WEIGHTS(td)    &
                          * INT_P_WEIGHTS(pd)

END DO
END DO



DO rd = 1,NUM_R_QUAD_POINTS
DO td = 1,NUM_T_QUAD_POINTS
DO pd = 1,NUM_P_QUAD_POINTS


    tpd = FP_tpd_Map(td,pd)
!    ui = 1
!
!    Tmp_U_Value = 0.0_idp
!    DO d = 0,DEGREE
!        Here = FP_FEM_Node_Map(re,d)
!
!        TMP_U_Value(ui) = TMP_U_Value(ui)                           &
!                        + SUM( FP_Coeff_Vector_Orig( Here, :, ui )  &
!                                * Ylm_Values( :, tpd, te, pe )       )      &
!                        * Lagrange_Poly_Table( d, rd, 0 )
!
!    END DO  ! d


    DO s = 3,5
!        StaredSource(tpd, rd, s ) = REAL(Tmp_U_Value(ui), KIND = idp)**6         &
!                                    * Block_Source_Si(rd,td,pd,re,te,pe,s-2)
        StaredSource(tpd, rd, s ) = Block_Source_Si(rd,td,pd,re,te,pe,s-2)

    END DO

END DO ! pd
END DO ! td
END DO ! rd







END SUBROUTINE Calc_XCFC_X_CurVals







!+204+###########################################################################!
!                                                                                !
!                  Create_XCFC_X_Vector                                         !
!                                                                                !
!################################################################################!
SUBROUTINE Create_XCFC_X_Vector( re, te, pe )



INTEGER, INTENT(IN)                                                     ::  re, te, pe

INTEGER                                                                 ::  pd, td, rd, tpd,     &
                                                                            l, m, d,        &
                                                                            lm_loc, u,ui

INTEGER                                                                 ::  Current_i_Location

COMPLEX(KIND = idp), DIMENSION(1:5)                                     ::  RHS_TMP
COMPLEX(KIND = idp)                                                     ::  Test
COMPLEX(KIND = idp)                                                     ::  Common_Basis
REAL(KIND = idp)                                                        ::  Combined_Weights

COMPLEX(KIND = idp)                                                     ::  Inner, Middle

DO ui = 3,5
DO d = 0,DEGREE
DO lm_loc = 1,LM_LENGTH

    RHS_TMP = 0.0_idp
    DO rd = 1,NUM_R_QUAD_POINTS
        RHS_TMP(ui) =  RHS_TMP(ui)                                      &
                        + 8.0_idp * pi * GR_Source_Scalar               &
                        * SUM(StaredSource(:, rd, ui)                   &
                                * Ylm_CC_Values( :, lm_loc, te, pe)     &
                                * TP_Int_Weights(:)                 )   &
                        * Lagrange_Poly_Table(d, rd, 0)                 &
                        * R_Int_Weights(rd)


    END DO  ! rd Loop

    Current_i_Location = FP_Beta_Array_Map(re,d,ui-2,lm_loc)

    FP_Source_Vector_X(Current_i_Location)                &
        = FP_Source_Vector_X(Current_i_Location)          &
        + RHS_TMP(ui)

    
    
END DO  ! lm_loc Loop
END DO  ! d Loop
END DO  ! ui




END SUBROUTINE Create_XCFC_X_Vector






















!+102+###########################################################################!
!                                                                                !
!           Calc_XCFC_ConFacto_Source                                                   !
!                                                                                !
!################################################################################!
SUBROUTINE XCFC_Calc_ConFactor_Source()

INTEGER                                                     ::  re, te, pe,     &
                                                                d, lm,          &
                                                                rd, tpd, td, pd

REAL(KIND = idp)                                            ::  deltar_overtwo,     &
                                                                deltat_overtwo,     &
                                                                deltap_overtwo


FP_Source_Vector = 0.0_idp

DO re = 0,NUM_R_ELEMENTS-1

    DELTAR_OVERTWO = 0.5_idp *(rlocs(re + 1) - rlocs(re))
    CUR_R_LOCS(:) = deltar_overtwo * (INT_R_LOCATIONS(:)+1.0_idp) + rlocs(re)
    R_SQUARE(:) = CUR_R_LOCS(:)*CUR_R_LOCS(:)


    DO pe = 0,NUM_P_ELEMENTS-1

        deltap_overtwo = 0.5_idp * (plocs(pe + 1)-plocs(pe))

        DO te = 0,NUM_T_ELEMENTS-1

            deltat_overtwo = 0.5_idp*(tlocs(te + 1) - tlocs(te))
            CUR_T_LOCS(:) = deltat_overtwo * (INT_T_LOCATIONS(:)+1.0_idp) + tlocs(te)

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

            CALL Calc_XCFC_ConFact_CurVals( re, te, pe,       &
                                              DELTAR_OVERTWO,   &
                                              DELTAT_OVERTWO,   &
                                              DELTAP_OVERTWO    )

            CALL Create_XCFC_ConFact_Vector( re, te, pe       )



        END DO ! te Loop
    END DO ! pe Loop
END DO ! re Loop


END SUBROUTINE XCFC_Calc_ConFactor_Source






!+202+###########################################################################!
!                                                                                !
!                  Calc_Current_Values          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_XCFC_ConFact_CurVals( re, te, pe,                                  &
                                DELTAR_OVERTWO,                             &
                                DELTAT_OVERTWO,                             &
                                DELTAP_OVERTWO                              )

INTEGER, INTENT(IN)                                             ::  re, te, pe


REAL(KIND = idp), INTENT(IN)                                    ::  DELTAR_OVERTWO,     &
                                                                    DELTAT_OVERTWO,     &
                                                                    DELTAP_OVERTWO


COMPLEX(KIND = idp), DIMENSION(1:8)                             ::  Tmp_U_Value
COMPLEX(KIND = idp), DIMENSION(1:8,1:3)                         ::  Tmp_U_DRV_Value

REAL(KIND = idp), DIMENSION( 1:NUM_TP_QUAD_POINTS,  &
                             1:NUM_R_QUAD_POINTS,   &
                             1:3                    )          ::  CUR_VAL_X


REAL(KIND = idp), DIMENSION( 1:NUM_TP_QUAD_POINTS,  &
                             1:NUM_R_QUAD_POINTS,   &
                             1:3, 1:3               )          ::  CUR_DRV_X


REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3,3) ::  Ahat_Array
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3 )  ::  f
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points )    ::  Tmp

INTEGER                                                         ::  tpd, td, pd, rd
INTEGER                                                         ::  d, Here
INTEGER                                                         ::  ui, s
INTEGER                                                         ::  i, j

                          !                                                 !
                         !!                                                 !!
                        !!!          Initialize Local Quadratures           !!!
                         !!                                                 !!
                          !                                                 !
R_Int_Weights(:) = DELTAR_OVERTWO * R_SQUARE(:) * INT_R_WEIGHTS(:)

DO td = 1,NUM_T_QUAD_POINTS
DO pd = 1,NUM_P_QUAD_POINTS


    tpd = FP_tpd_Map(td,pd)
    TP_Int_Weights( tpd ) = TP_SIN_VAL(tpd)                       &
                          * DELTAT_OVERTWO * INT_T_WEIGHTS(td)    &
                          * INT_P_WEIGHTS(pd)
END DO
END DO



DO rd = 1,NUM_R_QUAD_POINTS
DO td = 1,NUM_T_QUAD_POINTS
DO pd = 1,NUM_P_QUAD_POINTS


    tpd = FP_tpd_Map(td,pd)
    

    ! Calc Flat Metric Terms
    f(tpd,rd,1) = 1.0_idp
    f(tpd,rd,2) = R_Square(rd)
    f(tpd,rd,3) = R_Square(rd) * TP_SIN_SQUARE(tpd)
    

    ! Calc Origional Conformal Factor For Source Scaling
    
    Tmp_U_Value = 0.0_idp
    Tmp_U_DRV_Value = 0.0_idp
    DO d = 0,DEGREE
        Here = FP_FEM_Node_Map(re,d)
        ui = 1

!        TMP_U_Value(1) = TMP_U_Value(1)                           &
!                        + SUM( FP_Coeff_Vector_Orig( Here, :, ui )  &
!                                * Ylm_Values( :, tpd, te, pe )       )      &
!                        * Lagrange_Poly_Table( d, rd, 0 )

        TMP_U_Value(2) = TMP_U_Value(2)                           &
                        + SUM( FP_Coeff_Vector( Here, :, ui )  &
                                * Ylm_Values( :, tpd, te, pe )       )      &
                        * Lagrange_Poly_Table( d, rd, 0 )


    END DO  ! d

    ! Calculate Current Psi Value and Scaled Source Term
    s = 1
!    Orig_Val_Psi(tpd, rd )    = REAL(Tmp_U_Value(1), KIND = idp)
    Cur_Val_Psi(tpd,rd)       = REAL(Tmp_U_Value(2), KIND = idp)
    StaredSource(tpd, rd, s ) = Orig_Val_Psi(tpd, rd )**6         &
                              * Block_Source_E(rd,td,pd,re,te,pe)

    StaredSource(tpd, rd, s ) = Block_Source_E(rd,td,pd,re,te,pe)


!    PRINT*,"Here"
    ! Calculate Current X Values
    DO ui = 6,8
    DO d = 0,DEGREE
        Here = FP_FEM_Node_Map(re,d)

        TMP_U_Value(ui) = TMP_U_Value(ui)                           &
                        + SUM( FP_Coeff_Vector( Here, :, ui )         &
                                * Ylm_Values( :, tpd, te, pe )       )  &
                        * Lagrange_Poly_Table( d, rd, 0 )

        TMP_U_DRV_Value(ui,1)   = TMP_U_DRV_Value(ui,1)                     &
                                + SUM( FP_Coeff_Vector( Here, :, ui )       &
                                * Ylm_Values( :, tpd, te, pe )       )      &
                                * Lagrange_Poly_Table( d, rd, 1 )           &
                                / DELTAR_OVERTWO


        TMP_U_DRV_Value(ui,2)   = TMP_U_DRV_Value(ui,2)                     &
                                + SUM( FP_Coeff_Vector( Here, :, ui )       &
                                * Ylm_dt_Values( :, tpd, te, pe)     )      &
                                * Lagrange_Poly_Table( d, rd, 0)

        TMP_U_DRV_Value(ui,3)   = TMP_U_DRV_Value(ui,3)                     &
                                + SUM( FP_Coeff_Vector( Here, :, ui )       &
                                * Ylm_dp_Values( :, tpd, te, pe)     )      &
                                * Lagrange_Poly_Table( d, rd, 0)


    END DO  ! d
    END DO  ! ui
!    PRINT*,"There"

    CUR_Val_X(tpd,rd,1:3)   = REAL(TMP_U_Value(6:8), kind = idp)

    CUR_DRV_X(tpd,rd,1:3,1) = REAL(TMP_U_DRV_Value(6:8,1), kind = idp)
    CUR_DRV_X(tpd,rd,1:3,2) = REAL(TMP_U_DRV_Value(6:8,2), kind = idp)
    CUR_DRV_X(tpd,rd,1:3,3) = REAL(TMP_U_DRV_Value(6:8,3), kind = idp)

!    PRINT*,"Cur_Val_X",re,te,pe
!    PRINT*,Cur_Val_X(tpd,rd,1:3)
!    PRINT*,"Cur_DRV_X"
!    PRINT*,CUR_DRV_X(tpd,rd,1:3,1:3)

END DO ! pd
END DO ! td
END DO ! rd




!CALL Calc_Ahat( Ahat )

CALL Calc_Ahat( Ahat_Array,                             &
                NUM_R_QUAD_POINTS, NUM_TP_QUAD_POINTS,  &
                Cur_R_Locs, R_SQUARE,                   &
                TP_RSIN_SQUARE, TP_COTAN_VAL,           &
                CUR_VAL_X, CUR_DRV_X                  )




TMP = 0.0_idp
DO i = 1,3
DO j = 1,3
    TMP(:,:) = TMP(:,:)           &
        + f(:,:,i)*f(:,:,j) * (Ahat_Array(:,:,i,j))**2

END DO ! i
END DO ! j









ui = 1
SourceTerm(:,:,ui) = -2.0_idp * pi * GR_Source_Scalar                          &
                    / Cur_Val_Psi(:,:) * StaredSource(:,:,1)                &   ! Physical Source
                  - 1.0_idp / ( 8.0_idp * Cur_Val_Psi(:,:)**7) * TMP(:,:)       ! Geometry Source

!ui = 1
!SourceTerm(:,:,ui) = -2.0_idp * pi * GR_Source_Scalar                          &
!                    * Cur_Val_Psi(:,:)**5 * StaredSource(:,:,1)/(Orig_Val_Psi(:,: )**6)    &   ! Physical Source
!                  - 1.0_idp / ( 8.0_idp * Cur_Val_Psi(:,:)**7) * TMP(:,:)       ! Geometry Source




END SUBROUTINE Calc_XCFC_ConFact_CurVals







!+204+###########################################################################!
!                                                                                !
!                  CREATE_3D_RHS_VECTOR                                          !
!                                                                                !
!################################################################################!
SUBROUTINE Create_XCFC_ConFact_Vector( re, te, pe )



INTEGER, INTENT(IN)                                                     ::  re, te, pe

INTEGER                                                                 ::  pd, td, rd, tpd,     &
                                                                            l, m, d,        &
                                                                            lm_loc, u,ui

INTEGER                                                                 ::  Current_i_Location

COMPLEX(KIND = idp), DIMENSION(1:5)                                     ::  RHS_TMP
COMPLEX(KIND = idp)                                                     ::  Test
COMPLEX(KIND = idp)                                                     ::  Common_Basis
REAL(KIND = idp)                                                        ::  Combined_Weights

COMPLEX(KIND = idp)                                                     ::  Inner, Middle

ui = 1
DO lm_loc = 1,LM_LENGTH
DO d = 0,DEGREE


    RHS_TMP = 0.0_idp
    DO rd = 1,NUM_R_QUAD_POINTS

        RHS_TMP(ui) =  RHS_TMP(ui)                                          &
                         + SUM( SourceTerm( :, rd,1 )                         &
                               * Ylm_CC_Values( :, lm_loc, te, pe)          &
                               * TP_Int_Weights(:)                     )    &
                       * Lagrange_Poly_Table(d, rd, 0)                      &
                       * R_Int_Weights(rd)

    END DO  ! rd Loop
    

    Current_i_Location = FP_FEM_Node_Map(re,d)
    FP_Source_Vector(Current_i_Location,lm_loc,ui)                &
        = FP_Source_Vector(Current_i_Location,lm_loc,ui)          &
        + RHS_TMP(ui)

END DO  ! d Loop
END DO  ! lm_loc Loop


END SUBROUTINE Create_XCFC_ConFact_Vector















!+102+###########################################################################!
!                                                                                !
!           XCFC_Calc_Lapse_Source                                                   !
!                                                                                !
!################################################################################!
SUBROUTINE XCFC_Calc_Lapse_Source()

INTEGER                                                     ::  re, te, pe,     &
                                                                d, lm,          &
                                                                rd, tpd, td, pd

REAL(KIND = idp)                                            ::  deltar_overtwo,     &
                                                                deltat_overtwo,     &
                                                                deltap_overtwo


FP_Source_Vector = 0.0_idp

DO re = 0,NUM_R_ELEMENTS-1

    DELTAR_OVERTWO = 0.5_idp *(rlocs(re + 1) - rlocs(re))
    CUR_R_LOCS(:) = deltar_overtwo * (INT_R_LOCATIONS(:)+1.0_idp) + rlocs(re)
    R_SQUARE(:) = CUR_R_LOCS(:)*CUR_R_LOCS(:)


    DO pe = 0,NUM_P_ELEMENTS-1

        deltap_overtwo = 0.5_idp * (plocs(pe + 1)-plocs(pe))

        DO te = 0,NUM_T_ELEMENTS-1

            deltat_overtwo = 0.5_idp*(tlocs(te + 1) - tlocs(te))
            CUR_T_LOCS(:) = deltat_overtwo * (INT_T_LOCATIONS(:)+1.0_idp) + tlocs(te)

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

            CALL Calc_XCFC_Lapse_CurVals( re, te, pe,       &
                                              DELTAR_OVERTWO,   &
                                              DELTAT_OVERTWO,   &
                                              DELTAP_OVERTWO    )

            CALL Create_XCFC_Lapse_Vector( re, te, pe       )



        END DO ! te Loop
    END DO ! pe Loop
END DO ! re Loop


END SUBROUTINE XCFC_Calc_Lapse_Source






!+202+###########################################################################!
!                                                                                !
!                  Calc_Current_Values          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_XCFC_Lapse_CurVals( re, te, pe,                                  &
                                DELTAR_OVERTWO,                             &
                                DELTAT_OVERTWO,                             &
                                DELTAP_OVERTWO                              )

INTEGER, INTENT(IN)                                             ::  re, te, pe


REAL(KIND = idp), INTENT(IN)                                    ::  DELTAR_OVERTWO,     &
                                                                    DELTAT_OVERTWO,     &
                                                                    DELTAP_OVERTWO


COMPLEX(KIND = idp), DIMENSION(1:8)                             ::  Tmp_U_Value
COMPLEX(KIND = idp), DIMENSION(1:8,1:3)                         ::  Tmp_U_DRV_Value

REAL(KIND = idp), DIMENSION( 1:NUM_TP_QUAD_POINTS,  &
                             1:NUM_R_QUAD_POINTS,   &
                             1:3                    )          ::  CUR_VAL_X


REAL(KIND = idp), DIMENSION( 1:NUM_TP_QUAD_POINTS,  &
                             1:NUM_R_QUAD_POINTS,   &
                             1:3, 1:3               )          ::  CUR_DRV_X


REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3,3) ::  Ahat_Array
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3 )  ::  f
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points )    ::  Tmp

INTEGER                                                         ::  tpd, td, pd, rd
INTEGER                                                         ::  d, Here
INTEGER                                                         ::  ui, s
INTEGER                                                         ::  i, j

                          !                                                 !
                         !!                                                 !!
                        !!!          Initialize Local Quadratures           !!!
                         !!                                                 !!
                          !                                                 !
R_Int_Weights(:) = DELTAR_OVERTWO * R_SQUARE(:) * INT_R_WEIGHTS(:)

DO td = 1,NUM_T_QUAD_POINTS
DO pd = 1,NUM_P_QUAD_POINTS


    tpd = FP_tpd_Map(td,pd)
    TP_Int_Weights( tpd ) = TP_SIN_VAL(tpd)                       &
                          * DELTAT_OVERTWO * INT_T_WEIGHTS(td)    &
                          * INT_P_WEIGHTS(pd)
END DO
END DO



DO rd = 1,NUM_R_QUAD_POINTS
DO td = 1,NUM_T_QUAD_POINTS
DO pd = 1,NUM_P_QUAD_POINTS


    tpd = FP_tpd_Map(td,pd)
    ui = 2

    ! Calc Flat Metric Terms
    f(tpd,rd,1) = 1.0_idp
    f(tpd,rd,2) = R_Square(rd)
    f(tpd,rd,3) = R_Square(rd) * TP_SIN_SQUARE(tpd)
    

    ! Calc Origional Conformal Factor For Source Scaling
    Tmp_U_Value = 0.0_idp
    Tmp_U_DRV_Value = 0.0_idp
    DO d = 0,DEGREE
        Here = FP_FEM_Node_Map(re,d)

!        TMP_U_Value(1) = TMP_U_Value(1)                           &
!                        + SUM( FP_Coeff_Vector_Orig( Here, :, 1 )  &
!                                * Ylm_Values( :, tpd, te, pe )       )      &
!                        * Lagrange_Poly_Table( d, rd, 0 )

        TMP_U_Value(2) = TMP_U_Value(2)                           &
                        + SUM( FP_Coeff_Vector( Here, :, 1 )  &
                                * Ylm_Values( :, tpd, te, pe )       )      &
                        * Lagrange_Poly_Table( d, rd, 0 )


        TMP_U_Value(3) = TMP_U_Value(3)                           &
                        + SUM( FP_Coeff_Vector( Here, :, 2 )  &
                                * Ylm_Values( :, tpd, te, pe )       )      &
                        * Lagrange_Poly_Table( d, rd, 0 )

    END DO  ! d

    ! Calculate Current Psi Value and Scaled Source Term
    s = 2
!    Orig_Val_Psi(tpd, rd )    = REAL(Tmp_U_Value(1), KIND = idp)
    Cur_Val_Psi(tpd,rd)       = REAL(Tmp_U_Value(2), KIND = idp)
    Cur_Val_AlphaPsi(tpd,rd)       = REAL(Tmp_U_Value(3), KIND = idp)

!    StaredSource(tpd, rd, s ) = Orig_Val_Psi(tpd, rd )**6                           &
!                                * (Block_Source_E(rd,td,pd,re,te,pe)                &
!                                    + 2.0_idp*Block_Source_S(rd,td,pd,re,te,pe)     )
    
    StaredSource(tpd, rd, s ) = Block_Source_E(rd,td,pd,re,te,pe)           &
                                + 2.0_idp*Block_Source_S(rd,td,pd,re,te,pe)


    ! Calculate Current X Values
    DO ui = 6,8
    DO d = 0,DEGREE
        Here = FP_FEM_Node_Map(re,d)

        TMP_U_Value(ui) = TMP_U_Value(ui)                           &
                        + SUM( FP_Coeff_Vector( Here, :, ui )         &
                                * Ylm_Values( :, tpd, te, pe )       )  &
                        * Lagrange_Poly_Table( d, rd, 0 )

        TMP_U_DRV_Value(ui,1)   = TMP_U_DRV_Value(ui,1)                     &
                                + SUM( FP_Coeff_Vector( Here, :, ui )       &
                                * Ylm_Values( :, tpd, te, pe )       )      &
                                * Lagrange_Poly_Table( d, rd, 1 )           &
                                / DELTAR_OVERTWO


        TMP_U_DRV_Value(ui,2)   = TMP_U_DRV_Value(ui,2)                     &
                                + SUM( FP_Coeff_Vector( Here, :, ui )       &
                                * Ylm_dt_Values( :, tpd, te, pe)     )      &
                                * Lagrange_Poly_Table( d, rd, 0)

        TMP_U_DRV_Value(ui,3)   = TMP_U_DRV_Value(ui,3)                     &
                                + SUM( FP_Coeff_Vector( Here, :, ui )       &
                                * Ylm_dp_Values( :, tpd, te, pe)     )      &
                                * Lagrange_Poly_Table( d, rd, 0)


    END DO  ! d
    END DO  ! ui


    CUR_Val_X(tpd,rd,1:3)   = REAL(TMP_U_Value(6:8), kind = idp)

    CUR_DRV_X(tpd,rd,1:3,1) = REAL(TMP_U_DRV_Value(6:8,1), kind = idp)
    CUR_DRV_X(tpd,rd,1:3,2) = REAL(TMP_U_DRV_Value(6:8,2), kind = idp)
    CUR_DRV_X(tpd,rd,1:3,3) = REAL(TMP_U_DRV_Value(6:8,3), kind = idp)

END DO ! pd
END DO ! td
END DO ! rd




!CALL Calc_Ahat( Ahat )

CALL Calc_Ahat( Ahat_Array,                             &
                NUM_R_QUAD_POINTS, NUM_TP_QUAD_POINTS,  &
                Cur_R_Locs, R_SQUARE,                   &
                TP_RSIN_SQUARE, TP_COTAN_VAL,           &
                CUR_VAL_X, CUR_DRV_X                  )






TMP = 0.0_idp
DO i = 1,3
DO j = 1,3
    TMP(:,:) = TMP(:,:)           &
        + f(:,:,i)*f(:,:,j) * (Ahat_Array(:,:,i,j))**2

END DO ! i
END DO ! j



SourceTerm(:,:,2) = 2.0_idp * pi * GR_Source_Scalar                           &
                    * Cur_Val_AlphaPsi(:,:)/(Cur_Val_Psi(:,:)**2)           &
                    * StaredSource(:,:,2)                                   &   ! Physical Source
                  + (7.0_idp*Cur_Val_AlphaPsi(:,:))                         &
                    / ( 8.0_idp * Cur_Val_Psi(:,:)**8) * TMP(:,:)               ! Geometry Source




END SUBROUTINE Calc_XCFC_Lapse_CurVals







!+204+###########################################################################!
!                                                                                !
!                  CREATE_3D_RHS_VECTOR                                          !
!                                                                                !
!################################################################################!
SUBROUTINE Create_XCFC_Lapse_Vector( re, te, pe )



INTEGER, INTENT(IN)                                                     ::  re, te, pe

INTEGER                                                                 ::  pd, td, rd, tpd,     &
                                                                            l, m, d,        &
                                                                            lm_loc, u,ui

INTEGER                                                                 ::  Current_i_Location

COMPLEX(KIND = idp), DIMENSION(1:5)                                     ::  RHS_TMP
COMPLEX(KIND = idp)                                                     ::  Test
COMPLEX(KIND = idp)                                                     ::  Common_Basis
REAL(KIND = idp)                                                        ::  Combined_Weights

COMPLEX(KIND = idp)                                                     ::  Inner, Middle

ui = 2
DO lm_loc = 1,LM_LENGTH
DO d = 0,DEGREE


    RHS_TMP = 0.0_idp
    DO rd = 1,NUM_R_QUAD_POINTS

        RHS_TMP(ui) =  RHS_TMP(ui)                                          &
                         + SUM( SourceTerm( :, rd, 2 )                         &
                               * Ylm_CC_Values( :, lm_loc, te, pe)          &
                               * TP_Int_Weights(:)                     )    &
                       * Lagrange_Poly_Table(d, rd, 0)                      &
                       * R_Int_Weights(rd)

    END DO  ! rd Loop
    

    Current_i_Location = FP_FEM_Node_Map(re,d)
    FP_Source_Vector(Current_i_Location,lm_loc,ui)                &
        = FP_Source_Vector(Current_i_Location,lm_loc,ui)          &
        + RHS_TMP(ui)

END DO  ! d Loop
END DO  ! lm_loc Loop


END SUBROUTINE Create_XCFC_Lapse_Vector

















!+102+###########################################################################!
!                                                                                !
!           XCFC_Calc_Shift_Source                                                   !
!                                                                                !
!################################################################################!
SUBROUTINE XCFC_Calc_Shift_Source()

INTEGER                                                     ::  re, te, pe,     &
                                                                d, lm,          &
                                                                rd, tpd, td, pd

REAL(KIND = idp)                                            ::  deltar_overtwo,     &
                                                                deltat_overtwo,     &
                                                                deltap_overtwo


FP_Source_Vector = 0.0_idp

DO re = 0,NUM_R_ELEMENTS-1

    DELTAR_OVERTWO = 0.5_idp *(rlocs(re + 1) - rlocs(re))
    CUR_R_LOCS(:) = deltar_overtwo * (INT_R_LOCATIONS(:)+1.0_idp) + rlocs(re)
    R_SQUARE(:) = CUR_R_LOCS(:)*CUR_R_LOCS(:)


    DO pe = 0,NUM_P_ELEMENTS-1

        deltap_overtwo = 0.5_idp * (plocs(pe + 1)-plocs(pe))

        DO te = 0,NUM_T_ELEMENTS-1

            deltat_overtwo = 0.5_idp*(tlocs(te + 1) - tlocs(te))
            CUR_T_LOCS(:) = deltat_overtwo * (INT_T_LOCATIONS(:)+1.0_idp) + tlocs(te)

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

            CALL Calc_XCFC_Shift_CurVals( re, te, pe,       &
                                              DELTAR_OVERTWO,   &
                                              DELTAT_OVERTWO,   &
                                              DELTAP_OVERTWO    )

            CALL Calc_XCFC_Shift_Vector( re, te, pe       )



        END DO ! te Loop
    END DO ! pe Loop
END DO ! re Loop


END SUBROUTINE XCFC_Calc_Shift_Source








!+202+###########################################################################!
!                                                                                !
!                  Calc_Current_Values          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_XCFC_Shift_CurVals( re, te, pe,                                  &
                                DELTAR_OVERTWO,                             &
                                DELTAT_OVERTWO,                             &
                                DELTAP_OVERTWO                              )

INTEGER, INTENT(IN)                                             ::  re, te, pe


REAL(KIND = idp), INTENT(IN)                                    ::  DELTAR_OVERTWO,     &
                                                                    DELTAT_OVERTWO,     &
                                                                    DELTAP_OVERTWO


COMPLEX(KIND = idp), DIMENSION(1:8)                             ::  Tmp_U_Value
COMPLEX(KIND = idp), DIMENSION(1:8,1:3)                         ::  Tmp_U_DRV_Value

REAL(KIND = idp), DIMENSION( 1:NUM_TP_QUAD_POINTS,  &
                             1:NUM_R_QUAD_POINTS,   &
                             1:3                    )          ::  CUR_VAL_X


REAL(KIND = idp), DIMENSION( 1:NUM_TP_QUAD_POINTS,  &
                             1:NUM_R_QUAD_POINTS,   &
                             1:3, 1:3               )          ::  CUR_DRV_X


REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3,3) ::  Ahat_Array
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3 )  ::  f
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points )    ::  Tmp

REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3 )  ::  n_Array

INTEGER                                                         ::  tpd, td, pd, rd
INTEGER                                                         ::  d, Here
INTEGER                                                         ::  ui, s
INTEGER                                                         ::  i, j

                          !                                                 !
                         !!                                                 !!
                        !!!          Initialize Local Quadratures           !!!
                         !!                                                 !!
                          !                                                 !
R_Int_Weights(:) = DELTAR_OVERTWO * R_SQUARE(:) * INT_R_WEIGHTS(:)

DO td = 1,NUM_T_QUAD_POINTS
DO pd = 1,NUM_P_QUAD_POINTS


    tpd = FP_tpd_Map(td,pd)
    TP_Int_Weights( tpd ) = TP_SIN_VAL(tpd)                       &
                          * DELTAT_OVERTWO * INT_T_WEIGHTS(td)    &
                          * INT_P_WEIGHTS(pd)
END DO
END DO



DO rd = 1,NUM_R_QUAD_POINTS
DO td = 1,NUM_T_QUAD_POINTS
DO pd = 1,NUM_P_QUAD_POINTS


    tpd = FP_tpd_Map(td,pd)


    ! Calc Flat Metric Terms
    f(tpd,rd,1) = 1.0_idp
    f(tpd,rd,2) = R_Square(rd)
    f(tpd,rd,3) = R_Square(rd) * TP_SIN_SQUARE(tpd)
    

    ! Calc Origional Conformal Factor For Source Scaling
    Tmp_U_Value = 0.0_idp
    TMP_U_DRV_Value = 0.0_idp
    DO ui = 1,2
    DO d  = 0,Degree
        Here = FP_FEM_Node_Map(re,d)
        TMP_U_Value(ui) = TMP_U_Value(ui)                                   &
                        + SUM( FP_Coeff_Vector( Here, :, ui )                &
                                * Ylm_Values( :, tpd, te, pe )       )      &
                        * Lagrange_Poly_Table( d, rd, 0 )


        TMP_U_DRV_Value(ui,1) = TMP_U_DRV_Value(ui,1)                        &
                            + SUM( FP_Coeff_Vector( Here, :, ui )            &
                                    * Ylm_Values( :, tpd, te, pe )       )  &
                            * Lagrange_Poly_Table( d, rd, 1 )


        TMP_U_DRV_Value(ui,2)   = TMP_U_DRV_Value(ui,2)                     &
                                + SUM( FP_Coeff_Vector( Here, :, ui )       &
                                * Ylm_dt_Values( :, tpd, te, pe)     )      &
                                * Lagrange_Poly_Table( d, rd, 0)

        TMP_U_DRV_Value(ui,3)   = TMP_U_DRV_Value(ui,3)                     &
                                + SUM( FP_Coeff_Vector( Here, :, ui )       &
                                * Ylm_dp_Values( :, tpd, te, pe)     )      &
                                * Lagrange_Poly_Table( d, rd, 0)

    END DO  ! d
    END DO  ! UI

    Cur_Val_Psi(tpd,rd)         = REAL(Tmp_U_Value(1), KIND = idp)
    Cur_Val_AlphaPsi(tpd,rd)    = REAL(Tmp_U_Value(2), KIND = idp)

    Cur_DRV_Psi(tpd,rd,:)       =  REAL(Tmp_U_DRV_Value(1,:), KIND = idp )
    Cur_DRV_AlphaPsi(tpd,rd,:)  =  REAL(Tmp_U_DRV_Value(2,:), KIND = idp )



!    TMP_U_Value(1) = 0.0_idp
!    DO d = 0,DEGREE
!        Here = FP_FEM_Node_Map(re,d)
!
!        TMP_U_Value(1) = TMP_U_Value(1)                           &
!                        + SUM( FP_Coeff_Vector_Orig( Here, :, 1 )  &
!                                * Ylm_Values( :, tpd, te, pe )       )      &
!                        * Lagrange_Poly_Table( d, rd, 0 )
!
!    END DO
!    Orig_Val_Psi(tpd, rd )      = REAL(Tmp_U_Value(1), KIND = idp)






    ! Calculate Current X Values
    DO ui = 6,8
    DO d = 0,DEGREE
        Here = FP_FEM_Node_Map(re,d)

        TMP_U_Value(ui) = TMP_U_Value(ui)                           &
                        + SUM( FP_Coeff_Vector( Here, :, ui )         &
                                * Ylm_Values( :, tpd, te, pe )       )  &
                        * Lagrange_Poly_Table( d, rd, 0 )

        TMP_U_DRV_Value(ui,1)   = TMP_U_DRV_Value(ui,1)                     &
                                + SUM( FP_Coeff_Vector( Here, :, ui )       &
                                * Ylm_Values( :, tpd, te, pe )       )      &
                                * Lagrange_Poly_Table( d, rd, 1 )           &
                                / DELTAR_OVERTWO


        TMP_U_DRV_Value(ui,2)   = TMP_U_DRV_Value(ui,2)                     &
                                + SUM( FP_Coeff_Vector( Here, :, ui )       &
                                * Ylm_dt_Values( :, tpd, te, pe)     )      &
                                * Lagrange_Poly_Table( d, rd, 0)

        TMP_U_DRV_Value(ui,3)   = TMP_U_DRV_Value(ui,3)                     &
                                + SUM( FP_Coeff_Vector( Here, :, ui )       &
                                * Ylm_dp_Values( :, tpd, te, pe)     )      &
                                * Lagrange_Poly_Table( d, rd, 0)


    END DO  ! d
    END DO  ! ui


    CUR_Val_X(tpd,rd,1:3)   = REAL(TMP_U_Value(6:8), kind = idp)

    CUR_DRV_X(tpd,rd,1:3,1) = REAL(TMP_U_DRV_Value(6:8,1), kind = idp)
    CUR_DRV_X(tpd,rd,1:3,2) = REAL(TMP_U_DRV_Value(6:8,2), kind = idp)
    CUR_DRV_X(tpd,rd,1:3,3) = REAL(TMP_U_DRV_Value(6:8,3), kind = idp)

    DO s = 3,5
!        StaredSource(tpd, rd, s ) = Orig_Val_Psi(tpd, rd )**6         &
!                                    * Block_Source_Si(rd,td,pd,re,te,pe,s-2)

        StaredSource(tpd, rd, s ) = Block_Source_Si(rd,td,pd,re,te,pe,s-2)
    END DO


END DO ! pd
END DO ! td
END DO ! rd




CALL Calc_Ahat( Ahat_Array,                             &
                NUM_R_QUAD_POINTS, NUM_TP_QUAD_POINTS,  &
                Cur_R_Locs, R_SQUARE,                   &
                TP_RSIN_SQUARE, TP_COTAN_VAL,           &
                CUR_VAL_X, CUR_DRV_X                  )





DO i = 1,3
    n_Array(:,:,i) = Cur_DRV_AlphaPsi(:,:,i)/(Cur_VAL_Psi(:,:)**7)       &
                   - 7.0_idp * Cur_DRV_Psi(:,:,i)                      &
                     * Cur_Val_AlphaPsi(:,:)/(Cur_Val_psi(:,:)**8)
END DO




DO ui = 3,5

SourceTerm(:,:,ui) = 16.0_idp * pi * GR_Source_Scalar               &   ! Physical Source
                    * Cur_Val_AlphaPsi(:,:)/(Cur_Val_Psi(:,:)**7)   &
                    * StaredSource(:,:,ui)                          &
                  + 2.0_idp                                         &   ! Geometry Source
                    * ( N_Array(:,:,1)*Ahat_Array(:,:,ui-2,1)       &
                        + N_Array(:,:,2)*Ahat_Array(:,:,ui-2,2)     &
                        + N_Array(:,:,3)*Ahat_Array(:,:,ui-2,3)     )

END DO



END SUBROUTINE Calc_XCFC_Shift_CurVals




!+104+###########################################################################!
!                                                                                !
!           XCFC_Calc_Shift_Source                                                   !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_XCFC_Shift_Vector(re, te, pe)


INTEGER, INTENT(IN)                                                     ::  re, te, pe

INTEGER                                                                 ::  pd, td, rd, tpd,     &
                                                                            l, m, d,        &
                                                                            lm_loc, u,ui

INTEGER                                                                 ::  Current_i_Location

COMPLEX(KIND = idp), DIMENSION(1:5)                                     ::  RHS_TMP
COMPLEX(KIND = idp)                                                     ::  Test
COMPLEX(KIND = idp)                                                     ::  Common_Basis
REAL(KIND = idp)                                                        ::  Combined_Weights

COMPLEX(KIND = idp)                                                     ::  Inner, Middle

DO ui = 3,5
DO d = 0,DEGREE
DO lm_loc = 1,LM_LENGTH

    RHS_TMP = 0.0_idp
    DO rd = 1,NUM_R_QUAD_POINTS
        RHS_TMP(ui) =  RHS_TMP(ui)                              &
                        + SUM( SourceTerm( :, rd, ui )              &
                        * Ylm_CC_Values( :, lm_loc, te, pe)     &
                        * TP_Int_Weights(:)                 )   &
                        * Lagrange_Poly_Table(d, rd, 0)         &
                        * R_Int_Weights(rd)

!        PRINT*,TP_Int_Weights(:)
    END DO  ! rd Loop

    Current_i_Location = FP_Beta_Array_Map(re,d,ui-2,lm_loc)

    FP_Source_Vector_Beta(Current_i_Location)                &
        = FP_Source_Vector_Beta(Current_i_Location)          &
        + RHS_TMP(ui)

    
END DO  ! lm_loc Loop
END DO  ! d Loop
END DO  ! ui




END SUBROUTINE Calc_XCFC_Shift_Vector





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

ALLOCATE( CUR_DRV_PSI(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:3) )
ALLOCATE( CUR_DRV_ALPHAPSI(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:3) )
ALLOCATE( CUR_DRV_BETA(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:3, 1:3) )

ALLOCATE( Beta_DRV_Trace(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS) )

ALLOCATE( StaredSource( 1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:5) )
ALLOCATE( SourceTerm( 1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:5 ) )

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

DEALLOCATE( CUR_DRV_PSI )
DEALLOCATE( CUR_DRV_ALPHAPSI )
DEALLOCATE( CUR_DRV_BETA )

DEALLOCATE( Beta_DRV_Trace )

DEALLOCATE( StaredSource )
DEALLOCATE( SourceTerm )

END SUBROUTINE Deallocate_XCFC_Source_Variables

END MODULE XCFC_Source_Vector_Module
