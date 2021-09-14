   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE ADM_Mass_Module                                                       !##!
!##!                                                                         !##!
!##!_________________________________________________________________________!##!
!##!                                                                         !##!
!##!                                                                         !##!
!##!=========================================================================!##!
!##!                                                                         !##!
!##!    Contains:                                                            !##!
!##!                                                                         !##!
!##!        +101+                   Calc_ADM_Mass                            !##!
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
            ONLY :  Pi

USE Units_Module, &
            ONLY :  GR_Source_Scalar

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
                    plocs

USE Variables_Source, &
            ONLY :  Block_Source_E

USE FP_Functions_Mapping, &
            ONLY :  FP_tpd_Map

USE XCFC_Source_Functions_Module, &
            ONLY :  Calc_Val_On_Elem_TypeA,         &
                    Calc_Val_And_Drv_On_Elem_TypeB

USE Functions_Jacobian, &
            ONLY :  Calc_Ahat

IMPLICIT NONE




CONTAINS
!+101+##################################################################!
!                                                                       !
!          Calc_ADM_Mass                          				        !
!                                                                       !
!#######################################################################!
SUBROUTINE Calc_ADM_Mass( ADM_Mass )

REAL(idp), INTENT(OUT)                              ::  ADM_Mass




INTEGER                                             ::  re, te, pe
INTEGER                                             ::  rd, td, pd, tpd

REAL(idp)                                           ::  DROT, DTOT, DPOT
REAL(idp), DIMENSION(1:Num_R_Quad_Points)           ::  crlocs
REAL(idp), DIMENSION(1:Num_T_Quad_Points)           ::  ctlocs
REAL(idp), DIMENSION(1:Num_R_Quad_Points)           ::  rSquare

REAL(idp), DIMENSION(1:Num_TP_Quad_Points)          ::  TP_Sin_Val
REAL(idp), DIMENSION(1:Num_TP_Quad_Points)          ::  TP_Cotan_Val

REAL(idp), DIMENSION(1:Num_TP_Quad_Points,          &
                     1:Num_R_Quad_Points)           ::  TP_Rsin_Square

REAL(idp), DIMENSION(1:Num_TP_Quad_Points,          &
                     1:Num_R_Quad_Points)           ::  Int_Weights

REAL(idp), DIMENSION(1:Num_TP_Quad_Points,          &
                     1:Num_R_Quad_Points)           ::  Int_Source

REAL(idp)                                           ::  Int_Val



ADM_Mass = 0.0_idp


DO re = 0,Num_R_Elements-1
DO te = 0,Num_T_Elements-1
DO pe = 0,Num_P_Elements-1


    CALL Calc_Cur_Values( re, te, pe,           &
                          DROT, DTOT, DPOT,     &
                          crlocs, ctlocs,       &
                          rSquare,              &
                          TP_Sin_Val,           &
                          TP_Cotan_Val,         &
                          TP_rSin_Square        )


    CALL Calc_Int_Weights( DROT, DTOT,              &
                           rSquare, TP_Sin_Val,     &
                           Int_Weights              )


    CALL Calc_Int_Source( re, te, pe,      &
                          DROT, DTOT, DPOT,&
                          crlocs, rSquare, &
                          TP_Sin_Val,      &
                          TP_Cotan_Val,    &
                          TP_rSin_Square,  &
                          Int_Source       )


    CALL Calc_Element_ADM_Integral( Int_Weights,    &
                                    Int_Source,     &
                                    Int_Val         )


    ADM_Mass = ADM_Mass + Int_Val

END DO ! pe Loop
END DO ! te Loop
END DO ! re Loop


END SUBROUTINE Calc_ADM_Mass



!+201+###########################################################################!
!                                                                                !
!                  Calc_Cur_Values          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Cur_Values( re, te, pe,             &
                            DROT, DTOT, DPOT,       &
                            crlocs, ctlocs,         &
                            rSquare,                &
                            TP_Sin_Val,             &
                            TP_Cotan_Val,           &
                            TP_RSin_Square          )

INTEGER,   INTENT(IN)                                       ::  re, te, pe

REAL(idp), INTENT(OUT)                                      ::  DROT, DTOT, DPOT
REAL(idp), INTENT(OUT), DIMENSION(1:Num_R_Quad_Points)      ::  crlocs
REAL(idp), INTENT(OUT), DIMENSION(1:Num_T_Quad_Points)      ::  ctlocs
REAL(idp), INTENT(OUT), DIMENSION(1:Num_R_Quad_Points)      ::  rSquare

REAL(idp), INTENT(OUT), DIMENSION(1:Num_TP_Quad_Points)     ::  TP_Sin_Val
REAL(idp), INTENT(OUT), DIMENSION(1:Num_TP_Quad_Points)     ::  TP_Cotan_Val

REAL(idp), INTENT(OUT), DIMENSION(1:Num_TP_Quad_Points,     &
                                  1:Num_R_Quad_Points)      ::  TP_Rsin_Square


INTEGER                                                     ::  rd, td, pd, tpd


DROT = 0.5_idp * (rlocs(re + 1) - rlocs(re))
DTOT = 0.5_idp * (tlocs(te + 1) - tlocs(te))
DPOT = 0.5_idp * (plocs(pe + 1) - plocs(pe))

crlocs(:) = DROT * (INT_R_LOCATIONS(:)+1.0_idp) + rlocs(re)
ctlocs(:) = DTOT * (INT_T_LOCATIONS(:)+1.0_idp) + tlocs(te)

rSquare(:) = crlocs(:)*crlocs(:)

DO td = 1,NUM_T_QUAD_POINTS
DO pd = 1,NUM_P_QUAD_POINTS
    tpd = FP_tpd_Map(td,pd)
    TP_Sin_Val(tpd)    = DSIN(ctlocs(td))
    TP_Cotan_Val(tpd)  = 1.0_idp/DTAN(ctlocs(td))
END DO
END DO


DO rd = 1,NUM_R_QUAD_POINTS
    TP_RSIN_SQUARE(:,rd) = rSquare(rd)*TP_Sin_Val(:)*TP_Sin_Val(:)
END DO



END SUBROUTINE Calc_Cur_Values

!+202+###########################################################################!
!                                                                                !
!                  Calc_Int_Weights          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Int_Weights( DROT, DTOT,            &
                             R_Square, Sin_Val,     &
                             Int_Weights            )


REAL(idp), INTENT(IN)                                       ::  DROT, DTOT
REAL(idp), INTENT(IN),    DIMENSION(1:Num_R_Quad_Points)    ::  R_Square
REAL(idp), INTENT(IN),    DIMENSION(1:Num_TP_Quad_Points)   ::  Sin_Val
REAL(idp), INTENT(INOUT), DIMENSION(1:Num_TP_Quad_Points,   &
                                    1:Num_R_Quad_Points)    ::  Int_Weights

INTEGER                                                     ::  tpd, rd, td, pd



DO rd = 1,Num_R_Quad_Points
DO td = 1,NUM_T_QUAD_POINTS
DO pd = 1,NUM_P_QUAD_POINTS

   tpd = FP_tpd_Map(td,pd)
   Int_Weights( tpd, rd ) = DROT * R_SQUARE(rd) * INT_R_WEIGHTS(rd)  &
                          * DTOT * SIN_VAL(tpd) * INT_T_WEIGHTS(td)  &
                          * INT_P_WEIGHTS(pd)

END DO ! pd Loop
END DO ! td Loop
END DO ! rd Loop

END SUBROUTINE Calc_Int_Weights



!+203+###########################################################################!
!                                                                                !
!                  Calc_Current_Values          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Int_Source( re, te, pe,        &
                             DROT, DTOT, DPOT,  &
                             rlocs, rSquare,    &
                             TP_Sin_Val,        &
                             TP_Cotan_Val,      &
                             TP_RSin_Square,    &
                             Int_Source         )

INTEGER,    INTENT(IN)                                      ::  re, te, pe
REAL(idp),  INTENT(IN)                                      ::  DROT, DTOT, DPOT

REAL(idp),  INTENT(IN), DIMENSION(1:Num_R_Quad_Points)      ::  rlocs
REAL(idp),  INTENT(IN), DIMENSION(1:Num_R_Quad_Points)      ::  rSquare

REAL(idp),  INTENT(IN), DIMENSION(1:Num_TP_Quad_Points)     ::  TP_Sin_Val
REAL(idp),  INTENT(IN), DIMENSION(1:Num_TP_Quad_Points)     ::  TP_Cotan_Val

REAL(idp),  INTENT(IN), DIMENSION(1:Num_TP_Quad_Points,     &
                                  1:Num_R_Quad_Points)      ::  TP_Rsin_Square

REAL(idp),  INTENT(INOUT), DIMENSION(1:Num_TP_Quad_Points,  &
                                     1:Num_R_Quad_Points)   ::  Int_Source




REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3,3) ::  Ahat_Array
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3 )  ::  f
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points )    ::  AA_Array


REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points )    ::  Cur_Val_Psi
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3 )  ::  Cur_Val_X
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3,3) ::  Cur_Drv_X

INTEGER                                                         ::  tpd, rd, td, pd
INTEGER                                                         ::  i, j




CALL Calc_Val_On_Elem_TypeA( RE, TE, PE, Cur_Val_Psi, iU_CF )


CALL Calc_Val_And_Drv_On_Elem_TypeB( RE, TE, PE, DROT,      &
                                     CUR_Val_X(:,:,1),      &
                                     CUR_DRV_X(:,:,:,1),    &
                                     iU_X1, iVB_X           )

CALL Calc_Val_And_Drv_On_Elem_TypeB( RE, TE, PE, DROT,      &
                                     CUR_Val_X(:,:,2),      &
                                     CUR_DRV_X(:,:,:,2),    &
                                     iU_X2, iVB_X           )

CALL Calc_Val_And_Drv_On_Elem_TypeB( RE, TE, PE, DROT,      &
                                     CUR_Val_X(:,:,3),      &
                                     CUR_DRV_X(:,:,:,3),    &
                                     iU_X3, iVB_X           )

CALL Calc_Ahat( Ahat_Array,                             &
                NUM_R_QUAD_POINTS, NUM_TP_QUAD_POINTS,  &
                rlocs, rSquare,                         &
                TP_RSIN_SQUARE, TP_COTAN_VAL,           &
                CUR_VAL_X, CUR_DRV_X                  )





DO rd = 1,NUM_R_QUAD_POINTS
DO tpd = 1,NUM_T_QUAD_POINTS

    ! Calc Flat Metric Terms
    f(tpd,rd,1) = 1.0_idp
    f(tpd,rd,2) = rSquare(rd)
    f(tpd,rd,3) = rSquare(rd) * TP_Sin_Val(tpd) * TP_Sin_Val(tpd)
    
END DO ! tpd
END DO ! rd


AA_Array = 0.0_idp
DO i = 1,3
DO j = 1,3
    AA_Array(:,:) = AA_Array(:,:)           &
        + f(:,:,i)*f(:,:,j) * (Ahat_Array(:,:,i,j))**2

END DO ! i
END DO ! j




DO rd = 1,NUM_R_QUAD_POINTS
DO td = 1,NUM_T_QUAD_POINTS
DO pd = 1,NUM_P_QUAD_POINTS

    tpd = FP_tpd_Map(td,pd)

    Int_Source(tpd,rd) = GR_Source_Scalar                               &
                            / Cur_Val_Psi(tpd,rd)                       &
                            * Block_Source_E(rd,td,pd,re,te,pe)         &
                         + AA_Array(tpd,rd)                             &
                            / ( 16.0_idp * pi * Cur_Val_Psi(tpd,rd)**7)
                            

END DO ! pd
END DO ! td
END DO ! rd


END SUBROUTINE Calc_Int_Source







!+204+##########################################################################!
!                                                                               !
!                 Calc_Element_ADM_Integral                                     !
!                                                                               !
!###############################################################################!
SUBROUTINE Calc_Element_ADM_Integral( Int_Weights,      &
                                      Int_Source,       &
                                      Int_Val           )

REAL(idp),  INTENT(IN), DIMENSION(1:Num_TP_Quad_Points, &
                                  1:Num_R_Quad_Points)  ::  Int_Weights

REAL(idp),  INTENT(IN), DIMENSION(1:Num_TP_Quad_Points, &
                                  1:Num_R_Quad_Points)  ::  Int_Source


REAL(idp), INTENT(OUT)                                  ::  Int_Val

INTEGER                                                 ::  rd
COMPLEX(idp)                                            ::  Tmp_Val




Tmp_Val = 0.0_idp
DO rd = 1,NUM_R_QUAD_POINTS

    Tmp_Val =  Tmp_Val                      &
             + SUM( Int_Source( :, rd )     &
                  * Int_Weights(:, rd )     )

END DO  ! rd Loop
    

Int_Val = Tmp_Val


END SUBROUTINE Calc_Element_ADM_Integral









END MODULE ADM_Mass_Module
