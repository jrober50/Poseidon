   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE XCFC_Functions_Calc_Values_Module                                      !##!
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

USE Poseidon_Parameters, &
            ONLY :  DEGREE

USE Variables_Derived, &
            ONLY :  LM_Length

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

USE Variables_Tables, &
            ONLY :  Ylm_Values,                 &
                    Ylm_dt_Values,              &
                    Ylm_dp_Values,              &
                    Lagrange_Poly_Table

USE Variables_FP, &
            ONLY :  FP_Coeff_Vector_A,          &
                    FP_Coeff_Vector_B


USE FP_Functions_Mapping, &
            ONLY :  FP_tpd_Map,                 &
                    FP_FEM_Node_Map,            &
                    FP_Array_Map_TypeB


IMPLICIT NONE


CONTAINS


!+202+###########################################################################!
!                                                                                !
!                  Calc_Int_Weights          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Int_Weights( DROT, DTOT,            &
                             R_Square, Sin_Val,     &
                             rWeights, tpWeights    )


REAL(idp), INTENT(IN)                                       ::  DROT, DTOT
REAL(idp), INTENT(IN),    DIMENSION(1:Num_R_Quad_Points)    ::  R_Square
REAL(idp), INTENT(IN),    DIMENSION(1:Num_TP_Quad_Points)   ::  Sin_Val
REAL(idp), INTENT(INOUT), DIMENSION(1:Num_R_Quad_Points)    ::  rWeights
REAL(idp), INTENT(INOUT), DIMENSION(1:Num_TP_Quad_Points)   ::  tpWeights

INTEGER                                                     ::  tpd, td, pd



                         !                                                 !
                        !!                                                 !!
                       !!!          Initialize Local Quadratures           !!!
                        !!                                                 !!
                         !                                                 !
rWeights(:) = DROT * R_SQUARE(:) * INT_R_WEIGHTS(:)


DO td = 1,NUM_T_QUAD_POINTS
DO pd = 1,NUM_P_QUAD_POINTS

   tpd = FP_tpd_Map(td,pd)
   tpWeights( tpd ) = SIN_VAL(tpd)                   &
                    * DTOT * INT_T_WEIGHTS(td)      &
                    * INT_P_WEIGHTS(pd)

END DO
END DO

END SUBROUTINE Calc_Int_Weights













!+202+###########################################################################!
!                                                                                !
!                  Calc_ConFact_Values         !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Val_On_Elem_TypeA( RE, TE, PE, Val, iU )

INTEGER,   INTENT(IN)                       :: RE, TE, PE

REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points),   INTENT(OUT)  :: Val
INTEGER,   INTENT(IN)                       :: iU

REAL(idp)                                   :: TMP_Val
INTEGER                                     :: rd, tpd, d, Here


DO rd = 1,NUM_R_QUAD_POINTS
DO tpd = 1,NUM_TP_QUAD_POINTS

   
    TMP_Val = 0.0_idp
    DO d = 0,DEGREE
        Here = FP_FEM_Node_Map(re,d)
!        PRINT*,iU,tpd,te,pe, Ylm_Values( :, tpd, te, pe ), Lagrange_Poly_Table( d, rd, 0 )
    
        TMP_Val = TMP_Val                               &
            + SUM( FP_Coeff_Vector_A( Here, :, iU )     &
                   * Ylm_Values( :, tpd, te, pe )   )   &
            * Lagrange_Poly_Table( d, rd, 0 )


    END DO  ! d

    Val(tpd,rd)       = REAL(TMP_Val, KIND = idp)

END DO ! td
END DO ! rd

END SUBROUTINE Calc_Val_On_Elem_TypeA




!+202+###########################################################################!
!                                                                                !
!                  Calc_ConFact_Values         !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Drv_On_Elem_TypeA( RE, TE, PE, DROT, Drv, iU )

INTEGER,   INTENT(IN)                       :: RE, TE, PE
REAL(idp), INTENT(IN)                       :: DROT

REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3),  INTENT(OUT)  :: Drv
INTEGER,   INTENT(IN)                       :: iU

REAL(idp), DIMENSION(3)                     :: TMP_Drv
INTEGER                                     :: rd, tpd, d, Here


DO rd = 1,NUM_R_QUAD_POINTS
DO tpd = 1,NUM_TP_QUAD_POINTS

   
   TMP_Drv = 0.0_idp
   DO d = 0,DEGREE
       Here = FP_FEM_Node_Map(re,d)


       TMP_Drv(1) = TMP_Drv(1)                                 &
                  + SUM( FP_Coeff_Vector_A( Here, :, iU )      &
                        * Ylm_Values( :, tpd, te, pe )     )   &
                  * Lagrange_Poly_Table( d, rd, 1 )            &
                  / DROT


       TMP_Drv(2) = TMP_Drv(2)                                 &
                  + SUM( FP_Coeff_Vector_A( Here, :, iU )      &
                        * Ylm_dt_Values( :, tpd, te, pe)   )   &
                  * Lagrange_Poly_Table( d, rd, 0)

       TMP_Drv(3) = TMP_Drv(3)                                 &
                  + SUM( FP_Coeff_Vector_A( Here, :, iU )      &
                        * Ylm_dp_Values( :, tpd, te, pe)   )   &
                  * Lagrange_Poly_Table( d, rd, 0)

   END DO  ! d

   Drv(tpd,rd,1)       = REAL(TMP_Drv(1), KIND = idp)
   Drv(tpd,rd,2)       = REAL(TMP_Drv(2), KIND = idp)
   Drv(tpd,rd,3)       = REAL(TMP_Drv(2), KIND = idp)

END DO ! td
END DO ! rd

END SUBROUTINE Calc_Drv_On_Elem_TypeA



!+202+###########################################################################!
!                                                                                !
!                  Calc_ConFact_Values         !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Val_And_Drv_On_Elem_TypeA( RE, TE, PE, DROT, Val, Drv, iU )

INTEGER,   INTENT(IN)                       :: RE, TE, PE
REAL(idp), INTENT(IN)                       :: DROT

REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points),    INTENT(OUT)  :: Val
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3),  INTENT(OUT)  :: Drv
INTEGER,   INTENT(IN)                       :: iU

REAL(idp)                                   :: TMP_Val
REAL(idp), DIMENSION(3)                     :: TMP_Drv

INTEGER                                     :: rd, tpd, d, Here




DO rd = 1,NUM_R_QUAD_POINTS
DO tpd = 1,NUM_TP_QUAD_POINTS

   TMP_Val = 0.0_idp
   TMP_Drv = 0.0_idp
   DO d = 0,DEGREE
       Here = FP_FEM_Node_Map(re,d)

       TMP_Val = TMP_Val                                        &
               + SUM( FP_Coeff_Vector_A( Here, :, iU )          &
                       * Ylm_Values( :, tpd, te, pe )   )       &
               * Lagrange_Poly_Table( d, rd, 0 )


       TMP_Drv(1) = TMP_Drv(1)                                  &
                  + SUM( FP_Coeff_Vector_A( Here, :, iU )       &
                        * Ylm_Values( :, tpd, te, pe )     )    &
                  * Lagrange_Poly_Table( d, rd, 1 )             &
                  / DROT


       TMP_Drv(2) = TMP_Drv(2)                                  &
                  + SUM( FP_Coeff_Vector_A( Here, :, iU )       &
                        * Ylm_dt_Values( :, tpd, te, pe)   )    &
                  * Lagrange_Poly_Table( d, rd, 0)

       TMP_Drv(3) = TMP_Drv(3)                                  &
                  + SUM( FP_Coeff_Vector_A( Here, :, iU )       &
                        * Ylm_dp_Values( :, tpd, te, pe)   )    &
                  * Lagrange_Poly_Table( d, rd, 0)

   END DO  ! d

   Val(tpd,rd)         = REAL(TMP_Val,    KIND = idp)
   Drv(tpd,rd,1)       = REAL(TMP_Drv(1), KIND = idp)
   Drv(tpd,rd,2)       = REAL(TMP_Drv(2), KIND = idp)
   Drv(tpd,rd,3)       = REAL(TMP_Drv(2), KIND = idp)

END DO ! td
END DO ! rd

END SUBROUTINE Calc_Val_And_Drv_On_Elem_TypeA















!+202+###########################################################################!
!                                                                                !
!                  Calc_ConFact_Values         !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Val_On_Elem_TypeB( RE, TE, PE, Val, iU, iVB )

INTEGER,   INTENT(IN)                       :: RE, TE, PE

REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points),   INTENT(OUT)  :: Val
INTEGER,   INTENT(IN)                       :: iU, iVB

REAL(idp)                                   :: TMP_Val
INTEGER                                     :: rd, tpd, d, Here, There


DO rd = 1,NUM_R_QUAD_POINTS
DO tpd = 1,NUM_TP_QUAD_POINTS

   
    TMP_Val = 0.0_idp
    DO d = 0,DEGREE

        Here  = FP_Array_Map_TypeB(iU,iVB,re,d,1)
        There = FP_Array_Map_TypeB(iU,iVB,re,d,LM_Length)

        TMP_Val = TMP_Val                                   &
               + SUM( FP_Coeff_Vector_B( Here:There, iVB )  &
                       * Ylm_Values( :, tpd, te, pe )   )   &
               * Lagrange_Poly_Table( d, rd, 0 )


    END DO  ! d

    Val(tpd,rd)       = REAL(TMP_Val, KIND = idp)

END DO ! td
END DO ! rd

END SUBROUTINE Calc_Val_On_Elem_TypeB




!+202+###########################################################################!
!                                                                                !
!                  Calc_ConFact_Values         !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Drv_On_Elem_TypeB( RE, TE, PE, DROT, Drv, iU, iVB )

INTEGER,   INTENT(IN)                       :: RE, TE, PE
REAL(idp), INTENT(IN)                       :: DROT

REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3),  INTENT(OUT)  :: Drv
INTEGER,   INTENT(IN)                       :: iU, iVB

REAL(idp), DIMENSION(3)                     :: TMP_Drv
INTEGER                                     :: rd, tpd, d, Here, There


DO rd = 1,NUM_R_QUAD_POINTS
DO tpd = 1,NUM_TP_QUAD_POINTS

   
    TMP_Drv = 0.0_idp
    DO d = 0,DEGREE
        Here  = FP_Array_Map_TypeB(iU,iVB,re,d,1)
        There = FP_Array_Map_TypeB(iU,iVB,re,d,LM_Length)


        TMP_Drv(1) = TMP_Drv(1)                                 &
                  + SUM( FP_Coeff_Vector_B( Here:There, iVB )   &
                        * Ylm_Values( :, tpd, te, pe )     )    &
                  * Lagrange_Poly_Table( d, rd, 1 )             &
                  / DROT


        TMP_Drv(2) = TMP_Drv(2)                                 &
                  + SUM( FP_Coeff_Vector_B( Here:There, iVB )   &
                        * Ylm_dt_Values( :, tpd, te, pe)   )    &
                  * Lagrange_Poly_Table( d, rd, 0)

        TMP_Drv(3) = TMP_Drv(3)                                 &
                  + SUM( FP_Coeff_Vector_B( Here:There, iVB )   &
                        * Ylm_dp_Values( :, tpd, te, pe)   )    &
                  * Lagrange_Poly_Table( d, rd, 0)

    END DO  ! d

    Drv(tpd,rd,1)       = REAL(TMP_Drv(1), KIND = idp)
    Drv(tpd,rd,2)       = REAL(TMP_Drv(2), KIND = idp)
    Drv(tpd,rd,3)       = REAL(TMP_Drv(2), KIND = idp)

END DO ! td
END DO ! rd

END SUBROUTINE Calc_Drv_On_Elem_TypeB



!+202+###########################################################################!
!                                                                                !
!                  Calc_ConFact_Values         !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Val_And_Drv_On_Elem_TypeB( RE, TE, PE, DROT, Val, Drv, iU, iVB )

INTEGER,   INTENT(IN)                       :: RE, TE, PE
REAL(idp), INTENT(IN)                       :: DROT

REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points),    INTENT(OUT)  :: Val
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3),  INTENT(OUT)  :: Drv
INTEGER,   INTENT(IN)                       :: iU, iVB

REAL(idp)                                   :: TMP_Val
REAL(idp), DIMENSION(3)                     :: TMP_Drv

INTEGER                                     :: rd, tpd, d, Here, There




DO rd = 1,NUM_R_QUAD_POINTS
DO tpd = 1,NUM_TP_QUAD_POINTS

    TMP_Val = 0.0_idp
    TMP_Drv = 0.0_idp
    DO d = 0,DEGREE
        Here  = FP_Array_Map_TypeB(iU,iVB,re,d,1)
        There = FP_Array_Map_TypeB(iU,iVB,re,d,LM_Length)

        TMP_Val = TMP_Val                                       &
               + SUM( FP_Coeff_Vector_B( Here:There, iVB )      &
                       * Ylm_Values( :, tpd, te, pe )   )       &
               * Lagrange_Poly_Table( d, rd, 0 )


        TMP_Drv(1) = TMP_Drv(1)                                 &
                  + SUM( FP_Coeff_Vector_B( Here:There, iVB )   &
                        * Ylm_Values( :, tpd, te, pe )     )    &
                  * Lagrange_Poly_Table( d, rd, 1 )             &
                  / DROT


        TMP_Drv(2) = TMP_Drv(2)                                 &
                  + SUM( FP_Coeff_Vector_B( Here:There, iVB )   &
                        * Ylm_dt_Values( :, tpd, te, pe)   )    &
                  * Lagrange_Poly_Table( d, rd, 0)

        TMP_Drv(3) = TMP_Drv(3)                                 &
                  + SUM( FP_Coeff_Vector_B( Here:There, iVB )   &
                        * Ylm_dp_Values( :, tpd, te, pe)   )    &
                  * Lagrange_Poly_Table( d, rd, 0)

   END DO  ! d

   Val(tpd,rd)         = REAL(TMP_Val,    KIND = idp)
   Drv(tpd,rd,1)       = REAL(TMP_Drv(1), KIND = idp)
   Drv(tpd,rd,2)       = REAL(TMP_Drv(2), KIND = idp)
   Drv(tpd,rd,3)       = REAL(TMP_Drv(2), KIND = idp)

END DO ! td
END DO ! rd

END SUBROUTINE Calc_Val_And_Drv_On_Elem_TypeB


END MODULE XCFC_Functions_Calc_Values_Module
