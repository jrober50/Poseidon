   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Komar_Mass_Module                                                     !##!
!##!                                                                         !##!
!##!_________________________________________________________________________!##!
!##!                                                                         !##!
!##!                                                                         !##!
!##!=========================================================================!##!
!##!                                                                         !##!
!##!    Contains:                                                            !##!
!##!                                                                         !##!
!##!        +101+                   Calc_Komar_Mass                          !##!
!##!                                                                         !##!
!##!        +201+                   Calc_Cur_Values                          !##!
!##!        +202+                   Calc_Int_Weights                         !##!
!##!        +203+                   Calc_Int_Source                          !##!
!##!        +204+                   Calc_Element_ADM_Integral                !##!
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
                    iU_LF,                      &
                    iU_S1,                      &
                    iU_S2,                      &
                    iU_S3,                      &
                    iU_X1,                      &
                    iU_X2,                      &
                    iU_X3,                      &
                    iVB_S,                      &
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


USE Maps_Quadrature, &
            ONLY :  Map_To_tpd

USE Komar_Mass_Utilities_Module, &
            ONLY :  Calc_Val_and_Drv_On_Edge_TypeA, &
                    Calc_Val_On_Edge_TypeB,         &
                    Calc_Val_And_Drv_On_Edge_TypeB, &
                    Calc_Ahat_Edge


IMPLICIT NONE




CONTAINS
!+101+##################################################################!
!                                                                       !
!          Calc_Komar_Mass                                                  !
!                                                                       !
!#######################################################################!
SUBROUTINE Calc_Komar_Mass( Komar_Mass )

REAL(idp), INTENT(OUT)                              ::  Komar_Mass




INTEGER                                             ::  te, pe

REAL(idp)                                           ::  DROT, DTOT, DPOT
REAL(idp)                                           ::  crlocs
REAL(idp)                                           ::  rSquare


REAL(idp), DIMENSION(1:Num_T_Quad_Points)           ::  ctlocs

REAL(idp), DIMENSION(1:Num_TP_Quad_Points)          ::  TP_Sin_Val
REAL(idp), DIMENSION(1:Num_TP_Quad_Points)          ::  TP_Cotan_Val

REAL(idp), DIMENSION(1:Num_TP_Quad_Points)          ::  TP_Rsin_Square


REAL(idp), DIMENSION(1:Num_TP_Quad_Points)          ::  Int_Weights
REAL(idp), DIMENSION(1:Num_TP_Quad_Points)          ::  Int_Source

REAL(idp)                                           ::  Int_Val


PRINT*,"In Komar_Mass"
Komar_Mass = 0.0_idp


DO te = 0,Num_T_Elements-1
DO pe = 0,Num_P_Elements-1

    PRINT*,"Before Calc_Cur_Values",te,pe
    CALL Calc_Cur_Values( te, pe,           &
                          DROT, DTOT, DPOT,     &
                          crlocs, rSquare,       &
                          ctlocs,              &
                          TP_Sin_Val,           &
                          TP_Cotan_Val,         &
                          TP_rSin_Square        )


    PRINT*,"Before Calc_Int_Weights",te,pe
    CALL Calc_Int_Weights( DTOT,                    &
                           rSquare, TP_Sin_Val,     &
                           Int_Weights              )


    PRINT*,"Before Calc_Int_Source",te,pe
    CALL Calc_Int_Source( te, pe,           &
                          DROT, DTOT, DPOT, &
                          crlocs, rSquare,  &
                          TP_Sin_Val,       &
                          TP_Cotan_Val,     &
                          TP_rSin_Square,   &
                          Int_Source        )

    
    Int_Val = SUM( Int_Source( : ) * Int_Weights(:)  )

    PRINT*,"Int_Source"
    PRint*,Int_Source

    PRINT*,"Int_Weights"
    pRINT*,Int_Weights

    Komar_Mass = Komar_Mass + Int_Val

END DO ! pe Loop
END DO ! te Loop


END SUBROUTINE Calc_Komar_Mass



!+201+###########################################################################!
!                                                                                !
!                  Calc_Cur_Values          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Cur_Values( te, pe,                 &
                            DROT, DTOT, DPOT,       &
                            rVal, rSqr,             &
                            ctlocs,                 &
                            TP_Sin_Val,             &
                            TP_Cotan_Val,           &
                            TP_RSin_Square          )

INTEGER,   INTENT(IN)                                       ::  te, pe

REAL(idp), INTENT(OUT)                                      ::  DROT, DTOT, DPOT
REAL(idp), INTENT(OUT), DIMENSION(1:Num_T_Quad_Points)      ::  ctlocs

REAL(idp), INTENT(OUT)                                      ::  rVal
REAL(idp), INTENT(OUT)                                      ::  rSqr

REAL(idp), INTENT(OUT), DIMENSION(1:Num_TP_Quad_Points)     ::  TP_Sin_Val
REAL(idp), INTENT(OUT), DIMENSION(1:Num_TP_Quad_Points)     ::  TP_Cotan_Val

REAL(idp), INTENT(OUT), DIMENSION(1:Num_TP_Quad_Points)     ::  TP_Rsin_Square


INTEGER                                                     ::  td, pd, tpd


DROT = 0.5_idp * (rlocs(Num_R_Elements) - rlocs(Num_R_Elements-1))
DTOT = 0.5_idp * (tlocs(te + 1) - tlocs(te))
DPOT = 0.5_idp * (plocs(pe + 1) - plocs(pe))


ctlocs(:) = DTOT * (INT_T_LOCATIONS(:)+1.0_idp) + tlocs(te)


rVal = R_Outer
rSqr = rVal*rVal

DO td = 1,NUM_T_QUAD_POINTS
DO pd = 1,NUM_P_QUAD_POINTS
    tpd = Map_To_tpd(td,pd)
    TP_Sin_Val(tpd)    = DSIN(ctlocs(td))
    TP_Cotan_Val(tpd)  = 1.0_idp/DTAN(ctlocs(td))
END DO
END DO


TP_RSIN_SQUARE(:) = rSqr*TP_Sin_Val(:)*TP_Sin_Val(:)



END SUBROUTINE Calc_Cur_Values






!+202+###########################################################################!
!                                                                                !
!                  Calc_Int_Weights          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Int_Weights( DTOT,                  &
                             R_Square, Sin_Val,     &
                             Int_Weights            )


REAL(idp), INTENT(IN)                                       ::  DTOT
REAL(idp), INTENT(IN)                                       ::  R_Square
REAL(idp), INTENT(IN),    DIMENSION(1:Num_TP_Quad_Points)   ::  Sin_Val
REAL(idp), INTENT(INOUT), DIMENSION(1:Num_TP_Quad_Points)   ::  Int_Weights

INTEGER                                                     ::  tpd, td, pd



DO td = 1,NUM_T_QUAD_POINTS
DO pd = 1,NUM_P_QUAD_POINTS

   tpd = Map_To_tpd(td,pd)
   Int_Weights( tpd) = R_SQUARE                                 &
                     * DTOT * SIN_VAL(tpd) * INT_T_WEIGHTS(td)  &
                     * INT_P_WEIGHTS(pd)

END DO ! pd Loop
END DO ! td Loop

END SUBROUTINE Calc_Int_Weights








!+203+###########################################################################!
!                                                                                !
!                  Calc_Current_Values          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Int_Source( te, pe,             &
                            DROT, DTOT, DPOT,   &
                            rVal, rSqr,         &
                            TP_Sin_Val,         &
                            TP_Cotan_Val,       &
                            TP_RSin_Square,     &
                            Int_Source          )

INTEGER,    INTENT(IN)                                      ::  te, pe
REAL(idp),  INTENT(IN)                                      ::  DROT, DTOT, DPOT

REAL(idp),  INTENT(IN)                                      ::  rVal
REAL(idp),  INTENT(IN)                                      ::  rSqr

REAL(idp),  INTENT(IN), DIMENSION(1:Num_TP_Quad_Points)     ::  TP_Sin_Val
REAL(idp),  INTENT(IN), DIMENSION(1:Num_TP_Quad_Points)     ::  TP_Cotan_Val

REAL(idp),  INTENT(IN), DIMENSION(1:Num_TP_Quad_Points,     &
                                  1:Num_R_Quad_Points)      ::  TP_Rsin_Square

REAL(idp),  INTENT(INOUT), DIMENSION(1:Num_TP_Quad_Points)  ::  Int_Source


REAL(idp), DIMENSION(Num_TP_Quad_Points )                   ::  Cur_Val_Psi
REAL(idp), DIMENSION(Num_TP_Quad_Points, 3 )                ::  Cur_Drv_Psi

REAL(idp), DIMENSION(Num_TP_Quad_Points )                   ::  Cur_Val_AlphaPsi
REAL(idp), DIMENSION(Num_TP_Quad_Points, 3 )                ::  Cur_Drv_AlphaPsi

REAL(idp), DIMENSION(Num_TP_Quad_Points, 3 )                ::  Cur_Val_X
REAL(idp), DIMENSION(Num_TP_Quad_Points, 3, 3 )             ::  Cur_Drv_X

REAL(idp), DIMENSION(Num_TP_Quad_Points, 3 )                ::  Cur_Val_Beta

REAL(idp), DIMENSION(Num_TP_Quad_Points )                   ::  AB_Array
REAL(idp), DIMENSION(Num_TP_Quad_Points, 3 )                ::  f
REAL(idp), DIMENSION(Num_TP_Quad_Points, 3, 3 )             ::  Ahat_Array


INTEGER                                                     ::  tpd, td, pd
INTEGER                                                     ::  i




CALL Calc_Val_and_Drv_On_Edge_TypeA( DROT, TE, PE,          &
                                     Cur_Val_Psi,           &
                                     Cur_DRV_Psi,           &
                                     iU_CF                  )


CALL Calc_Val_and_Drv_On_Edge_TypeA( DROT, TE, PE,          &
                                     Cur_Val_AlphaPsi,      &
                                     Cur_Drv_AlphaPsi,      &
                                     iU_LF                  )

CALL Calc_Val_And_Drv_On_Edge_TypeB( DROT, TE, PE,          &
                                     CUR_Val_X(:,1),        &
                                     CUR_DRV_X(:,:,1),      &
                                     iU_X1, iVB_X           )

CALL Calc_Val_And_Drv_On_Edge_TypeB( DROT, TE, PE,          &
                                     CUR_Val_X(:,2),        &
                                     CUR_DRV_X(:,:,2),      &
                                     iU_X2, iVB_X           )

CALL Calc_Val_And_Drv_On_Edge_TypeB( DROT, TE, PE,          &
                                     CUR_Val_X(:,3),        &
                                     CUR_DRV_X(:,:,3),      &
                                     iU_X3, iVB_X           )

CALL Calc_Val_On_Edge_TypeB( TE, PE,                &
                             CUR_Val_Beta(:,1),     &
                             iU_S1, iVB_S           )

CALL Calc_Val_On_Edge_TypeB( TE, PE,                &
                             CUR_Val_Beta(:,3),     &
                             iU_S2, iVB_S           )

CALL Calc_Val_On_Edge_TypeB( TE, PE,                &
                             CUR_Val_Beta(:,3),     &
                             iU_S3, iVB_S           )



CALL Calc_Ahat_Edge( Ahat_Array,                    &
                     NUM_TP_QUAD_POINTS,            &
                     rVal, rSqr,                    &
                     TP_RSIN_SQUARE, TP_COTAN_VAL,  &
                     CUR_VAL_X, CUR_DRV_X           )





DO tpd = 1,NUM_T_QUAD_POINTS

    ! Calc Flat Metric Terms
    f(tpd,1) = 1.0_idp
    f(tpd,2) = rSqr
    f(tpd,3) = rSqr * TP_Sin_Val(tpd) * TP_Sin_Val(tpd)
    
END DO ! tpd


AB_Array = 0.0_idp
DO i = 1,3
    AB_Array(:) = AB_Array(:)           &
        + f(:,i) * Ahat_Array(:,1,i) * Cur_Val_Beta(:,i)
    PRINT*,"Beta ",i
    PRINT*,Cur_Val_Beta(:,i)
END DO ! i




DO td = 1,NUM_T_QUAD_POINTS
DO pd = 1,NUM_P_QUAD_POINTS

    tpd = Map_To_tpd(td,pd)

!    PRINT*,"Cur_Vals",Cur_Val_psi(tpd),Cur_Val_Alphapsi(tpd)
!    PRINT*,"Cur_Drvs",Cur_Drv_Psi(tpd,1), Cur_DRV_Alphapsi(tpd,1)

    Int_Source(tpd) = Cur_Val_Psi(tpd)**4                               &
                    * ( Cur_Val_Psi(tpd) * Cur_Drv_AlphaPsi(tpd,1)      &
                         - Cur_Val_AlphaPsi(tpd) * Cur_Drv_Psi(tpd,1)   &
                         - AB_Array(tpd)                                )
                            
    PRINT*,"Int_Source",Int_Source(tpd)
    PRINT*,  Cur_Drv_AlphaPsi(tpd,1),   &
             Cur_Drv_Psi(tpd,1),   &
             AB_Array(tpd)
END DO ! pd
END DO ! td


END SUBROUTINE Calc_Int_Source


















END MODULE Komar_Mass_Module

