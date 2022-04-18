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

USE Poseidon_Numbers_Module, &
            ONLY :  pi

USE Poseidon_Parameters, &
            ONLY :  DEGREE

USE Variables_Derived, &
            ONLY :  LM_Length

USE Variables_Mesh, &
            ONLY :  Num_P_Elements
            

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
                    Ylm_Elem_Values,            &
                    Ylm_Elem_dt_Values,         &
                    Ylm_Elem_dp_Values,         &
                    Lagrange_Poly_Table,        &
                    LagPoly_MultiLayer_Table

USE Variables_FP, &
            ONLY :  FP_Coeff_Vector_A,          &
                    FP_Coeff_Vector_B

USE Parameters_Variable_Indices, &
            ONLY :  iU_CF,                      &
                    iU_LF,                      &
                    iU_S1,                      &
                    iU_S2,                      &
                    iU_S3,                      &
                    iU_X1,                      &
                    iU_X2,                      &
                    iU_X2

USE Maps_Fixed_Point, &
            ONLY :  FP_Array_Map_TypeB

USE Maps_Domain, &
            ONLY :  Map_To_FEM_Node,            &
                    FEM_Elem_Map

USE Maps_Quadrature, &
            ONLY :  Map_To_tpd

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

    
   tpd = Map_To_tpd(td,pd)
   tpWeights( tpd ) = SIN_VAL(tpd)                   &
                    * DTOT * INT_T_WEIGHTS(td)      &
                    * INT_P_WEIGHTS(pd)



END DO
END DO

END SUBROUTINE Calc_Int_Weights





!+202+###########################################################################!
!                                                                                !
!                  Calc_Int_Weights          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Int_Weights_AMReX( DROT, DTOT,              &
                                    R_Square, Sin_Val,      &
                                    rWeights, tpWeights,    &
                                    Level                   )


REAL(idp), INTENT(IN)                                       ::  DROT, DTOT
REAL(idp), INTENT(IN),    DIMENSION(1:Num_R_Quad_Points)    ::  R_Square
REAL(idp), INTENT(IN),    DIMENSION(1:Num_TP_Quad_Points)   ::  Sin_Val
REAL(idp), INTENT(INOUT), DIMENSION(1:Num_R_Quad_Points)    ::  rWeights
REAL(idp), INTENT(INOUT), DIMENSION(1:Num_TP_Quad_Points)   ::  tpWeights
INTEGER,   INTENT(IN)                                       ::  Level


INTEGER                                                     ::  tpd, td, pd
!REAL(idp)                                                   ::  Int_P_Weight

!Int_P_Weight = (2*pi/(Num_P_Elements*2**Level))/Num_P_Quad_Points


                         !                                                 !
                        !!                                                 !!
                       !!!          Initialize Local Quadratures           !!!
                        !!                                                 !!
                         !                                                 !
rWeights(:) = DROT * R_SQUARE(:) * INT_R_WEIGHTS(:)


DO td = 1,NUM_T_QUAD_POINTS
DO pd = 1,NUM_P_QUAD_POINTS

    
    
   tpd = Map_To_tpd(td,pd)
   tpWeights( tpd ) = SIN_VAL(tpd)                   &
                    * DTOT * INT_T_WEIGHTS(td)      &
                    * Int_P_Weights(pd)

!    PRINT*,tpWeights(tpd),SIN_VAL(tpd)
END DO
END DO

END SUBROUTINE Calc_Int_Weights_AMReX








!+202+###########################################################################!
!                                                                                !
!                  Calc_ConFact_Values         !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Val_On_Elem_TypeA( iE, Val, iU, Level )

INTEGER,    INTENT(IN), DIMENSION(3)        ::  iE

REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points),   INTENT(OUT)  :: Val
INTEGER,    INTENT(IN)                      :: iU
INTEGER,    INTENT(IN)                      :: Level

REAL(idp)                                   :: TMP_Val
INTEGER                                     :: rd, tpd, d
INTEGER,        DIMENSION(0:DEGREE)         :: Here
INTEGER                                     :: iRE, iCT


!PRINT*,"1A"
#ifdef POSEIDON_AMREX_FLAG
    iRE = FEM_Elem_Map(iE(1),Level)
#else
    iRE = iE(1)
#endif




!PRINT*,"2A"

DO d = 0,DEGREE
    Here(d) = Map_To_FEM_Node(iRE,d)
END DO ! d




!PRINT*,"3A"
iCT = 0

DO rd = 1,NUM_R_QUAD_POINTS
DO tpd = 1,NUM_TP_QUAD_POINTS
    
    TMP_Val = 0.0_idp
    DO d = 0,DEGREE


#if POSEIDON_AMREX_FLAG

        TMP_Val = TMP_Val                               &
            + SUM( FP_Coeff_Vector_A( Here(d), :, iU )  &
                   * Ylm_Elem_Values( :, tpd )   )      &
            * LagPoly_MultiLayer_Table( d, rd, 0, iCT)

#else


        TMP_Val = TMP_Val                                   &
            + SUM( FP_Coeff_Vector_A( Here(d), :, iU )      &
                   * Ylm_Values( :, tpd, iE(2), iE(3) ) )   &
            * Lagrange_Poly_Table( d, rd, 0 )

#endif


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
SUBROUTINE Calc_Drv_On_Elem_TypeA( iE, DROT, Drv, iU, Level )

INTEGER,    INTENT(IN), DIMENSION(3)        ::  iE
REAL(idp), INTENT(IN)                       ::  DROT

REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3),  INTENT(OUT)  :: Drv
INTEGER,   INTENT(IN)                       ::  iU
INTEGER,    INTENT(IN)                      ::  Level

REAL(idp), DIMENSION(3)                     ::  TMP_Drv
INTEGER                                     ::  rd, tpd, d
INTEGER                                     ::  Here
INTEGER                                     ::  iRE, iCT



#ifdef POSEIDON_AMREX_FLAG
    iRE = FEM_Elem_Map(iE(1),Level)
#else
    iRE = iE(1)
#endif





!iCT = 2**(level+1) - mod(iE(1),2**level) - 2
iCT = 0


DO rd = 1,NUM_R_QUAD_POINTS
DO tpd = 1,NUM_TP_QUAD_POINTS

   
   TMP_Drv = 0.0_idp
   DO d = 0,DEGREE
       Here = Map_To_FEM_Node(iRE,d)

#if POSEIDON_AMREX_FLAG

        TMP_Drv(1) = TMP_Drv(1)                                 &
                   + SUM( FP_Coeff_Vector_A( Here, :, iU )      &
                         * Ylm_Elem_Values( :, tpd ) )          &
                   * LagPoly_MultiLayer_Table( d, rd, 1, iCT)   &
                   / DROT


        TMP_Drv(2) = TMP_Drv(2)                                 &
                   + SUM( FP_Coeff_Vector_A( Here, :, iU )      &
                         * Ylm_Elem_dt_Values( :, tpd )     )   &
                   * LagPoly_MultiLayer_Table( d, rd, 0, iCT)

        TMP_Drv(3) = TMP_Drv(3)                                 &
                   + SUM( FP_Coeff_Vector_A( Here, :, iU )      &
                         * Ylm_Elem_dp_Values( :, tpd )     )   &
                   * LagPoly_MultiLayer_Table( d, rd, 0, iCT)

#else

       TMP_Drv(1) = TMP_Drv(1)                                      &
                  + SUM( FP_Coeff_Vector_A( Here, :, iU )           &
                        * Ylm_Values( :, tpd, iE(2), iE(3)))        &
                  * Lagrange_Poly_Table( d, rd, 1 )                 &
                  / DROT


       TMP_Drv(2) = TMP_Drv(2)                                      &
                  + SUM( FP_Coeff_Vector_A( Here, :, iU )           &
                        * Ylm_dt_Values( :, tpd, iE(2), iE(3))   )  &
                  * Lagrange_Poly_Table( d, rd, 0)

       TMP_Drv(3) = TMP_Drv(3)                                      &
                  + SUM( FP_Coeff_Vector_A( Here, :, iU )           &
                        * Ylm_dp_Values( :, tpd, iE(2), iE(3))   )  &
                  * Lagrange_Poly_Table( d, rd, 0)

#endif


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
SUBROUTINE Calc_Val_And_Drv_On_Elem_TypeA( iE, DROT, Val, Drv, iU, Level)

INTEGER,   INTENT(IN), DIMENSION(3)         ::  iE
REAL(idp), INTENT(IN)                       ::  DROT

REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points),    INTENT(OUT)  :: Val
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3),  INTENT(OUT)  :: Drv
INTEGER,   INTENT(IN)                       ::  iU
INTEGER,   INTENT(IN)                       ::  Level

REAL(idp)                                   ::  TMP_Val
REAL(idp), DIMENSION(3)                     ::  TMP_Drv

INTEGER                                     ::  rd, tpd, d
INTEGER                                     ::  Here
INTEGER                                     ::  iRE, iCT


#ifdef POSEIDON_AMREX_FLAG
    iRE = FEM_Elem_Map(iE(1),Level)
#else
    iRE = iE(1)
#endif




!iCT = 2**(level+1) - mod(iE(1),2**level) - 2
iCT = 0

DO rd = 1,NUM_R_QUAD_POINTS
DO tpd = 1,NUM_TP_QUAD_POINTS

   TMP_Val = 0.0_idp
   TMP_Drv = 0.0_idp
   DO d = 0,DEGREE
       Here = Map_To_FEM_Node(iRE,d)


#ifdef POSEIDON_AMREX_FLAG
        TMP_Val = TMP_Val                                       &
            + SUM( FP_Coeff_Vector_A( Here, :, iU )             &
                   * Ylm_Elem_Values( :, tpd )   )              &
            * LagPoly_MultiLayer_Table( d, rd, 0, iCT)


        TMP_Drv(1) = TMP_Drv(1)                                 &
                   + SUM( FP_Coeff_Vector_A( Here, :, iU )      &
                         * Ylm_Elem_Values( :, tpd ) )          &
                   * LagPoly_MultiLayer_Table( d, rd, 1, iCT)   &
                   / DROT


        TMP_Drv(2) = TMP_Drv(2)                                 &
                   + SUM( FP_Coeff_Vector_A( Here, :, iU )      &
                         * Ylm_Elem_dt_Values( :, tpd )     )   &
                   * LagPoly_MultiLayer_Table( d, rd, 0, iCT)

        TMP_Drv(3) = TMP_Drv(3)                                 &
                   + SUM( FP_Coeff_Vector_A( Here, :, iU )      &
                         * Ylm_Elem_dp_Values( :, tpd )     )   &
                   * LagPoly_MultiLayer_Table( d, rd, 0, iCT)

#else

       TMP_Val = TMP_Val                                        &
               + SUM( FP_Coeff_Vector_A( Here, :, iU )          &
                       * Ylm_Values( :, tpd, iE(2), iE(3) )   ) &
               * Lagrange_Poly_Table( d, rd, 0 )


       TMP_Drv(1) = TMP_Drv(1)                                  &
                  + SUM( FP_Coeff_Vector_A( Here, :, iU )       &
                        * Ylm_Values( :, tpd, iE(2), iE(3) )  ) &
                  * Lagrange_Poly_Table( d, rd, 1 )             &
                  / DROT


       TMP_Drv(2) = TMP_Drv(2)                                  &
                  + SUM( FP_Coeff_Vector_A( Here, :, iU )       &
                        * Ylm_dt_Values( :, tpd, iE(2), iE(3))   )    &
                  * Lagrange_Poly_Table( d, rd, 0)

       TMP_Drv(3) = TMP_Drv(3)                                  &
                  + SUM( FP_Coeff_Vector_A( Here, :, iU )       &
                        * Ylm_dp_Values( :, tpd, iE(2), iE(3))   )    &
                  * Lagrange_Poly_Table( d, rd, 0)

#endif



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
SUBROUTINE Calc_Val_On_Elem_TypeB( iE, Val, iU, iVB, Level )

INTEGER,    INTENT(IN), DIMENSION(3)        ::  iE

REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points),   INTENT(OUT)  :: Val
INTEGER,   INTENT(IN)                       ::  iU, iVB
INTEGER,    INTENT(IN)                      ::  Level

REAL(idp)                                   ::  TMP_Val
INTEGER                                     ::  rd, tpd, d
INTEGER                                     ::  Here, There
INTEGER                                     ::  iRE, iCT



#ifdef POSEIDON_AMREX_FLAG
    iRE = FEM_Elem_Map(iE(1),Level)
#else
    iRE = iE(1)
#endif




!iCT = 2**(level+1) - mod(iE(1),2**level) - 2
iCT = 0

DO rd = 1,NUM_R_QUAD_POINTS
DO tpd = 1,NUM_TP_QUAD_POINTS

   
    TMP_Val = 0.0_idp
    DO d = 0,DEGREE

        Here  = FP_Array_Map_TypeB(iU,iVB,iRE,d,1)
        There = FP_Array_Map_TypeB(iU,iVB,iRE,d,LM_Length)


#ifdef POSEIDON_AMREX_FLAG

        TMP_Val = TMP_Val                                   &
               + SUM( FP_Coeff_Vector_B( Here:There, iVB )  &
                       * Ylm_Elem_Values( :, tpd )      )   &
               * LagPoly_MultiLayer_Table( d, rd, 0, iCT)

#else

        TMP_Val = TMP_Val                                   &
               + SUM( FP_Coeff_Vector_B( Here:There, iVB )  &
                       * Ylm_Values( :, tpd, iE(2), iE(3) )   )   &
               * Lagrange_Poly_Table( d, rd, 0 )
#endif


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
SUBROUTINE Calc_Drv_On_Elem_TypeB( iE, DROT, Drv, iU, iVB, Level )

INTEGER,    INTENT(IN), DIMENSION(3)        ::  iE
REAL(idp), INTENT(IN)                       :: DROT

REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3),  INTENT(OUT)  :: Drv
INTEGER,   INTENT(IN)                       :: iU, iVB
INTEGER,    INTENT(IN)                      :: Level


REAL(idp), DIMENSION(3)                     ::  TMP_Drv
INTEGER                                     ::  rd, tpd, d
INTEGER                                     ::  Here, There
INTEGER                                     ::  iRE, iCT



#ifdef POSEIDON_AMREX_FLAG
    iRE = FEM_Elem_Map(iE(1),Level)
#else
    iRE = iE(1)
#endif



!iCT = 2**(level+1) - mod(iE(1),2**level) - 2
iCT = 0


DO rd = 1,NUM_R_QUAD_POINTS
DO tpd = 1,NUM_TP_QUAD_POINTS

   
    TMP_Drv = 0.0_idp
    DO d = 0,DEGREE
        Here  = FP_Array_Map_TypeB(iU,iVB,iRE,d,1)
        There = FP_Array_Map_TypeB(iU,iVB,iRE,d,LM_Length)

#ifdef POSEIDON_AMREX_FLAG

        TMP_Drv(1) = TMP_Drv(1)                                 &
                  + SUM( FP_Coeff_Vector_B( Here:There, iVB )   &
                        * Ylm_Elem_Values( :, tpd )         )   &
                  * LagPoly_MultiLayer_Table( d, rd, 1, iCT)             &
                  / DROT


        TMP_Drv(2) = TMP_Drv(2)                                 &
                  + SUM( FP_Coeff_Vector_B( Here:There, iVB )   &
                        * Ylm_Elem_dt_Values( :, tpd )      )   &
                  * LagPoly_MultiLayer_Table( d, rd, 0, iCT)

        TMP_Drv(3) = TMP_Drv(3)                                 &
                  + SUM( FP_Coeff_Vector_B( Here:There, iVB )   &
                        * Ylm_Elem_dp_Values( :, tpd )      )   &
                  * LagPoly_MultiLayer_Table( d, rd, 0, iCT)

#else

        TMP_Drv(1) = TMP_Drv(1)                                 &
                  + SUM( FP_Coeff_Vector_B( Here:There, iVB )   &
                        * Ylm_Values( :, tpd, iE(2), iE(3) )     )    &
                  * Lagrange_Poly_Table( d, rd, 1 )             &
                  / DROT


        TMP_Drv(2) = TMP_Drv(2)                                 &
                  + SUM( FP_Coeff_Vector_B( Here:There, iVB )   &
                        * Ylm_dt_Values( :, tpd, iE(2), iE(3))   )    &
                  * Lagrange_Poly_Table( d, rd, 0)

        TMP_Drv(3) = TMP_Drv(3)                                 &
                  + SUM( FP_Coeff_Vector_B( Here:There, iVB )   &
                        * Ylm_dp_Values( :, tpd, iE(2), iE(3))   )    &
                  * Lagrange_Poly_Table( d, rd, 0)

#endif


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
SUBROUTINE Calc_Val_And_Drv_On_Elem_TypeB( iE, DROT, Val, Drv, iU, iVB, Level )

INTEGER,    INTENT(IN), DIMENSION(3)        ::  iE
REAL(idp), INTENT(IN)                       :: DROT

REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points),    INTENT(OUT)  :: Val
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3),  INTENT(OUT)  :: Drv
INTEGER,   INTENT(IN)                       :: iU, iVB
INTEGER,    INTENT(IN)                      :: Level



REAL(idp)                                   :: TMP_Val
REAL(idp), DIMENSION(3)                     :: TMP_Drv

INTEGER                                     ::  rd, tpd, d
INTEGER                                     ::  Here, There
INTEGER                                     ::  iRE, iCT



#ifdef POSEIDON_AMREX_FLAG
    iRE = FEM_Elem_Map(iE(1),Level)
#else
    iRE = iE(1)
#endif


!iCT = 2**(level+1) - mod(iE(1),2**level) - 2
iCT = 0



DO rd = 1,NUM_R_QUAD_POINTS
DO tpd = 1,NUM_TP_QUAD_POINTS


    TMP_Val = 0.0_idp
    TMP_Drv = 0.0_idp
    DO d = 0,DEGREE
        Here  = FP_Array_Map_TypeB(iU,iVB,iRE,d,1)
        There = FP_Array_Map_TypeB(iU,iVB,iRE,d,LM_Length)
        
#ifdef POSEIDON_AMREX_FLAG
        TMP_Val = TMP_Val                                       &
               + SUM( FP_Coeff_Vector_B( Here:There, iVB )      &
                       * Ylm_Elem_Values( :, tpd )          )   &
               * LagPoly_MultiLayer_Table( d, rd, 0, iCT)


        TMP_Drv(1) = TMP_Drv(1)                                 &
                  + SUM( FP_Coeff_Vector_B( Here:There, iVB )   &
                        * Ylm_Elem_Values( :, tpd )         )   &
                  * LagPoly_MultiLayer_Table( d, rd, 1, iCT)             &
                  / DROT


        TMP_Drv(2) = TMP_Drv(2)                                 &
                  + SUM( FP_Coeff_Vector_B( Here:There, iVB )   &
                        * Ylm_Elem_dt_Values( :, tpd )      )   &
                  * LagPoly_MultiLayer_Table( d, rd, 0, iCT)


        TMP_Drv(3) = TMP_Drv(3)                                 &
                  + SUM( FP_Coeff_Vector_B( Here:There, iVB )   &
                        * Ylm_Elem_dp_Values( :, tpd )      )   &
                  * LagPoly_MultiLayer_Table( d, rd, 0, iCT)

#else

        TMP_Val = TMP_Val                                       &
               + SUM( FP_Coeff_Vector_B( Here:There, iVB )      &
                       * Ylm_Values( :, tpd, iE(2), iE(3) )   )       &
               * Lagrange_Poly_Table( d, rd, 0 )



        TMP_Drv(1) = TMP_Drv(1)                                 &
                  + SUM( FP_Coeff_Vector_B( Here:There, iVB )   &
                        * Ylm_Values( :, tpd, iE(2), iE(3) )     )    &
                  * Lagrange_Poly_Table( d, rd, 1 )             &
                  / DROT

        TMP_Drv(2) = TMP_Drv(2)                                 &
                  + SUM( FP_Coeff_Vector_B( Here:There, iVB )   &
                        * Ylm_dt_Values( :, tpd, iE(2), iE(3))   )    &
                  * Lagrange_Poly_Table( d, rd, 0)



        TMP_Drv(3) = TMP_Drv(3)                                 &
                  + SUM( FP_Coeff_Vector_B( Here:There, iVB )   &
                        * Ylm_dp_Values( :, tpd, iE(2), iE(3))   )    &
                  * Lagrange_Poly_Table( d, rd, 0)



#endif


   END DO  ! d

   Val(tpd,rd)         = REAL(TMP_Val,    KIND = idp)
   Drv(tpd,rd,1)       = REAL(TMP_Drv(1), KIND = idp)
   Drv(tpd,rd,2)       = REAL(TMP_Drv(2), KIND = idp)
   Drv(tpd,rd,3)       = REAL(TMP_Drv(2), KIND = idp)

END DO ! td
END DO ! rd







END SUBROUTINE Calc_Val_And_Drv_On_Elem_TypeB
























END MODULE XCFC_Functions_Calc_Values_Module
