   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Komar_Mass_Utilities_Module                                           !##!
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
            ONLY : idp

USE Poseidon_Numbers_Module, &
            ONLY :  Pi,         &
                    TwoThirds,  &
                    FourThirds

USE Poseidon_Parameters, &
            ONLY :  Degree

USE Variables_Derived,  &
            ONLY :  LM_Length

USE Parameters_Variable_Indices, &
            ONLY :  iU_CF,                      &
                    iU_X1,                      &
                    iU_X2,                      &
                    iU_X3,                      &
                    iVB_X


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
                    plocs,                      &
                    R_Outer


USE Maps_Fixed_Point, &
            ONLY :  FP_Array_Map_TypeB

USE Maps_Domain, &
            ONLY :  Map_TO_FEM_Node

USE Functions_Jacobian, &
            ONLY :  Calc_Ahat

USE Variables_Vectors, &
            ONLY :  cVA_Coeff_Vector,          &
                    cVB_Coeff_Vector

USE Variables_Tables, &
            ONLY :  Ylm_Values,                 &
                    Ylm_dt_Values,              &
                    Ylm_dp_Values

USE Functions_Math, &
            ONLY :  Lagrange_Poly_Deriv

USE Functions_Quadrature, &
            ONLY :  Initialize_LGL_Quadrature_Locations



IMPLICIT NONE


CONTAINS


!+202+###########################################################################!
!                                                                                !
!                  Calc_ConFact_Values         !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Val_On_Edge_TypeA( TE, PE, Val, iU )

INTEGER,   INTENT(IN)                       :: TE, PE

REAL(idp), DIMENSION(Num_TP_Quad_Points),   INTENT(OUT)  :: Val
INTEGER,   INTENT(IN)                       :: iU

REAL(idp)                                   :: TMP_Val
INTEGER                                     :: tpd, Here


DO tpd = 1,NUM_TP_QUAD_POINTS

   
    Here = Map_To_FEM_Node(Num_R_Elements-1,Degree)
    
    TMP_Val = SUM( cVA_Coeff_Vector( Here, :, iU )     &
                  * Ylm_Values( :, tpd, te, pe )        )


    Val(tpd)       = REAL(TMP_Val, KIND = idp)

END DO ! td


END SUBROUTINE Calc_Val_On_Edge_TypeA



!+202+###########################################################################!
!                                                                                !
!                  Calc_ConFact_Values         !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Val_And_Drv_On_Edge_TypeA( DROT, TE, PE, Val, Drv, iU )

INTEGER,   INTENT(IN)                       :: TE, PE
REAL(idp),                                   INTENT(IN)     :: DROT
REAL(idp), DIMENSION(Num_TP_Quad_Points),    INTENT(OUT)    :: Val
REAL(idp), DIMENSION(Num_TP_Quad_Points,3),  INTENT(OUT)    :: Drv
INTEGER,   INTENT(IN)                       :: iU

REAL(idp)                                   :: TMP_Val
REAL(idp), DIMENSION(3)                     :: TMP_Drv

INTEGER                                     :: tpd, Here

REAL(idp), DIMENSION(0:DEGREE)              ::  Local_Locs
REAL(idp), DIMENSION(0:DEGREE)              ::  Lagrange_DRV_Values

Local_Locs = Initialize_LGL_Quadrature_Locations(Degree)

DO tpd = 1,NUM_TP_QUAD_POINTS

    TMP_Val = 0.0_idp
    TMP_Drv = 0.0_idp

    Here = Map_To_FEM_Node(Num_R_Elements-1,Degree)

    Lagrange_DRV_Values = Lagrange_Poly_Deriv(1.0_idp, Degree, Local_Locs)


    TMP_Val = SUM( cVA_Coeff_Vector( Here, :, iU )     &
                    * Ylm_Values( :, tpd, te, pe )        )


    TMP_Drv(1) = SUM( cVA_Coeff_Vector( Here, :, iU )      &
                    * Ylm_Values( :, tpd, te, pe )     )    &
                * Lagrange_DRV_Values ( Degree )            &
                / DROT


    TMP_Drv(2) = SUM( cVA_Coeff_Vector( Here, :, iU )      &
                    * Ylm_dt_Values( :, tpd, te, pe)   )

    TMP_Drv(3) = SUM( cVA_Coeff_Vector( Here, :, iU )      &
                    * Ylm_dp_Values( :, tpd, te, pe)   )


    Val(tpd)         = REAL(TMP_Val,    KIND = idp)
    Drv(tpd,1)       = REAL(TMP_Drv(1), KIND = idp)
    Drv(tpd,2)       = REAL(TMP_Drv(2), KIND = idp)
    Drv(tpd,3)       = REAL(TMP_Drv(2), KIND = idp)

END DO ! td



END SUBROUTINE Calc_Val_And_Drv_On_Edge_TypeA



!+202+###########################################################################!
!                                                                                !
!                  Calc_ConFact_Values         !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Val_On_Edge_TypeB( TE, PE, Val, iU, iVB )

INTEGER,                                    INTENT(IN)  :: TE, PE
REAL(idp),  DIMENSION(Num_TP_Quad_Points),  INTENT(OUT) :: Val
INTEGER,                                    INTENT(IN)  :: iU, iVB


INTEGER                                     ::  tpd, Here, There
REAL(idp)                                   ::  TMP_Val
REAL(idp), DIMENSION(0:DEGREE)              ::  Local_Locs
REAL(idp), DIMENSION(0:DEGREE)              ::  Lagrange_DRV_Values

Local_Locs = Initialize_LGL_Quadrature_Locations(Degree)

DO tpd = 1,NUM_TP_QUAD_POINTS

    TMP_Val = 0.0_idp


    Lagrange_DRV_Values = Lagrange_Poly_Deriv(1.0_idp, Degree, Local_Locs)

    Here  = FP_Array_Map_TypeB(iU,iVB,Num_R_Elements-1,Degree,1)
    There = FP_Array_Map_TypeB(iU,iVB,Num_R_Elements-1,Degree,LM_Length)

    TMP_Val = SUM( cVB_Coeff_Vector( Here:There, iVB )     &
                   * Ylm_Values( :, tpd, te, pe )           )



    Val(tpd)         = REAL(TMP_Val,    KIND = idp)

END DO ! td



END SUBROUTINE Calc_Val_On_Edge_TypeB



!+202+###########################################################################!
!                                                                                !
!                  Calc_ConFact_Values         !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Val_And_Drv_On_Edge_TypeB( DROT, TE, PE, Val, Drv, iU, iVB )

INTEGER,                                        INTENT(IN)      ::  TE, PE
REAL(idp),                                      INTENT(IN)      ::  DROT
REAL(idp),  DIMENSION(Num_TP_Quad_Points),      INTENT(OUT)     ::  Val
REAL(idp),  DIMENSION(Num_TP_Quad_Points,3),    INTENT(OUT)     ::  Drv
INTEGER,                                        INTENT(IN)      ::  iU, iVB

REAL(idp)                                       ::  TMP_Val
REAL(idp),  DIMENSION(3)                        ::  TMP_Drv

INTEGER                                         ::  tpd, Here, There

REAL(idp),  DIMENSION(0:DEGREE)                 ::  Local_Locs
REAL(idp),  DIMENSION(0:DEGREE)                 ::  Lagrange_DRV_Values

Local_Locs = Initialize_LGL_Quadrature_Locations(Degree)

DO tpd = 1,NUM_TP_QUAD_POINTS

    TMP_Val = 0.0_idp
    TMP_Drv = 0.0_idp


    Lagrange_DRV_Values = Lagrange_Poly_Deriv(1.0_idp, Degree, Local_Locs)

    Here  = FP_Array_Map_TypeB(iU,iVB,Num_R_Elements-1,Degree,1)
    There = FP_Array_Map_TypeB(iU,iVB,Num_R_Elements-1,Degree,LM_Length)

    TMP_Val = SUM( cVB_Coeff_Vector( Here:There, iVB )      &
                   * Ylm_Values( :, tpd, te, pe )   )


    TMP_Drv(1) = SUM( cVB_Coeff_Vector( Here:There, iVB )  &
                    * Ylm_Values( :, tpd, te, pe )     )    &
                * Lagrange_DRV_Values ( Degree )            &
                / DROT


    TMP_Drv(2) = SUM( cVB_Coeff_Vector( Here:There, iVB )   &
                    * Ylm_dt_Values( :, tpd, te, pe)   )

    TMP_Drv(3) = SUM( cVB_Coeff_Vector( Here:There, iVB )   &
                    * Ylm_dp_Values( :, tpd, te, pe)   )


    Val(tpd)         = REAL(TMP_Val,    KIND = idp)
    Drv(tpd,1)       = REAL(TMP_Drv(1), KIND = idp)
    Drv(tpd,2)       = REAL(TMP_Drv(2), KIND = idp)
    Drv(tpd,3)       = REAL(TMP_Drv(2), KIND = idp)

END DO ! td



END SUBROUTINE Calc_Val_And_Drv_On_Edge_TypeB







!+301+###########################################################################!
!                                                                                !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Ahat_Edge( Ahat_Array,                              &
                            NUM_TP_QUAD_POINTS,                     &
                            rVal, rSqr,                             &
                            RSIN_SQUARE, COTAN_VAL,                 &
                            CUR_VAL_X, CUR_DRV_X                    )

REAL(idp),  DIMENSION(Num_TP_Quad_Points, 3, 3 ),   INTENT(OUT) ::  Ahat_Array

INTEGER,                                            INTENT(IN)  ::  NUM_TP_QUAD_POINTS

REAL(idp),                                          INTENT(IN)  ::  rVal
REAL(idp),                                          INTENT(IN)  ::  rSqr
REAL(idp),  DIMENSION(Num_TP_Quad_Points),          INTENT(IN)  ::  COTAN_VAL
REAL(idp),  DIMENSION(1:Num_TP_Quad_Points),        INTENT(IN)  ::  RSIN_SQUARE

REAL(idp),  DIMENSION(1:NUM_TP_QUAD_POINTS, 1:3),   INTENT(IN)  ::  CUR_VAL_X
REAL(idp),  DIMENSION(1:NUM_TP_QUAD_POINTS, 1:3, 1:3),INTENT(IN)::  CUR_DRV_X


INTEGER                                                         ::  tpd



DO tpd = 1,Num_TP_Quad_Points
    
    Ahat_Array(tpd, 1, 1) =  FourThirds                 * CUR_DRV_X( tpd, 1, 1 )    &
                            - FourThirds/ rVal          * CUR_VAL_X( tpd, 1 )       &
                            - TwoThirds * COTAN_VAL(tpd)* CUR_VAL_X( tpd, 2 )       &
                            - TwoThirds                 * CUR_DRV_X( tpd, 2, 2 )    &
                            - TwoThirds                 * CUR_DRV_X( tpd, 3, 3 )

END DO



DO tpd = 1,Num_TP_Quad_Points

    Ahat_Array(tpd, 2, 1) =  CUR_DRV_X( tpd, 1, 2 )                                 &
                            +  CUR_DRV_X( tpd,  2, 1 )/rSqr

END DO


DO tpd = 1,Num_TP_Quad_Points

    Ahat_Array(tpd, 3, 1) =  CUR_DRV_X( tpd, 1, 3 )                                 &
                            +  CUR_DRV_X( tpd, 3, 1 )/RSIN_SQUARE(tpd)
END DO




Ahat_Array(:,1,2) = Ahat_Array(:,2,1)

DO tpd = 1,Num_TP_Quad_Points
    Ahat_Array(tpd, 2, 2) = TwoThirds   /(rSqr*rVal)            * CUR_VAL_X( tpd, 1 )       &
                            - (TwoThirds/ rSqr)*COTAN_VAL(tpd)  * CUR_VAL_X( tpd, 2 )       &
                            + FourThirds/ rSqr                  * CUR_DRV_X( tpd, 2, 2 )    &
                            - TwoThirds / rSqr                  * CUR_DRV_X( tpd, 1, 1 )    &
                            - TwoThirds / rSqr                  * CUR_DRV_X( tpd, 3, 3 )
END DO


DO tpd = 1,Num_TP_Quad_Points
    Ahat_Array(tpd, 3, 2) = CUR_DRV_X( tpd, 2, 3 )/rSqr                             &
                            + CUR_DRV_X( tpd, 3, 2 )/RSIN_SQUARE(tpd)
END DO


Ahat_Array(:,1,3) = Ahat_Array(:,3,1)
Ahat_Array(:,2,3) = Ahat_Array(:,3,2)



DO tpd = 1,Num_TP_Quad_Points
    Ahat_Array(tpd, 3, 3) = TwoThirds/( rVal* RSIN_SQUARE(tpd))                         &
                                                           * CUR_VAL_X( tpd, 1 )        &
                            + FourThirds * COTAN_VAL(tpd)/RSIN_SQUARE(tpd)              &
                                                           * CUR_VAL_X( tpd, 2 )        &
                            + FourThirds/ RSIN_SQUARE(tpd) * CUR_DRV_X( tpd, 3, 3 )     &
                            - TwoThirds / RSIN_SQUARE(tpd) * CUR_DRV_X( tpd, 1, 1 )     &
                            - TwoThirds / RSIN_SQUARE(tpd) * CUR_DRV_X( tpd, 2, 2 )
END DO


END SUBROUTINE Calc_Ahat_Edge






END MODULE Komar_Mass_Utilities_Module
