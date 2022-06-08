   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE ADM_Mass_In_Parts_Module                                              !##!
!##!                                                                         !##!
!##!_________________________________________________________________________!##!
!##!                                                                         !##!
!##!                                                                         !##!
!##!=========================================================================!##!
!##!                                                                         !##!
!##!    Contains:                                                            !##!
!##!                                                                         !##!
!##!        +101+                   Calc_ADM_Mass_In_Parts                   !##!
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

USE Poseidon_Units_Module, &
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

USE Maps_Quadrature, &
            ONLY :  Map_To_tpd,                 &
                    Quad_Map

USE XCFC_Functions_Calc_Values_Module, &
            ONLY :  Calc_Val_On_Elem_TypeA,         &
                    Calc_Val_And_Drv_On_Elem_TypeB

USE Functions_Jacobian, &
            ONLY :  Calc_Ahat

USE ADM_Mass_Module, &
            ONLY :  Calc_Cur_Values,                &
                    Calc_Int_Weights

IMPLICIT NONE




CONTAINS
!+101+##################################################################!
!                                                                       !
!          Calc_ADM_Mass_In_Parts                          				!
!                                                                       !
!#######################################################################!
SUBROUTINE Calc_ADM_Mass_In_Parts( ADM_Mass, ADM_Phys, ADM_Curve )

REAL(idp), INTENT(OUT)                              ::  ADM_Mass
REAL(idp), INTENT(OUT)                              ::  ADM_Phys
REAL(idp), INTENT(OUT)                              ::  ADM_Curve



INTEGER                                             ::  re, te, pe

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
                     1:Num_R_Quad_Points, 3 )       ::  Int_Source

REAL(idp), DIMENSION(3)                             ::  Int_Val



ADM_Mass = 0.0_idp
ADM_Phys = 0.0_idp
ADM_Curve = 0.0_idp

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


    CALL Calc_Int_Source_In_Parts( re, te, pe,      &
                                   DROT, DTOT, DPOT,&
                                   crlocs, rSquare, &
                                   TP_Sin_Val,      &
                                   TP_Cotan_Val,    &
                                   TP_rSin_Square,  &
                                   Int_Source       )


    CALL Calc_Element_ADM_Parts_Integral( Int_Weights,    &
                                          Int_Source,     &
                                          Int_Val         )


    ADM_Mass = ADM_Mass + Int_Val(1)
    ADM_Phys = ADM_Phys + Int_Val(2)
    ADM_Curve = ADM_Curve + Int_Val(3)

END DO ! pe Loop
END DO ! te Loop
END DO ! re Loop


END SUBROUTINE Calc_ADM_Mass_In_Parts






!+203+###########################################################################!
!                                                                                !
!                  Calc_Current_Values          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Int_Source_In_Parts( re, te, pe,        &
                                     DROT, DTOT, DPOT,  &
                                     rlocs, rSquare,    &
                                     TP_Sin_Val,        &
                                     TP_Cotan_Val,      &
                                     TP_RSin_Square,    &
                                     Int_Source         )

INTEGER,    INTENT(IN)                                          ::  re, te, pe
REAL(idp),  INTENT(IN)                                          ::  DROT, DTOT, DPOT

REAL(idp),  INTENT(IN), DIMENSION(1:Num_R_Quad_Points)          ::  rlocs
REAL(idp),  INTENT(IN), DIMENSION(1:Num_R_Quad_Points)          ::  rSquare

REAL(idp),  INTENT(IN), DIMENSION(1:Num_TP_Quad_Points)         ::  TP_Sin_Val
REAL(idp),  INTENT(IN), DIMENSION(1:Num_TP_Quad_Points)         ::  TP_Cotan_Val

REAL(idp),  INTENT(IN), DIMENSION(1:Num_TP_Quad_Points,         &
                                  1:Num_R_Quad_Points)          ::  TP_Rsin_Square

REAL(idp),  INTENT(INOUT), DIMENSION(1:Num_TP_Quad_Points,      &
                                     1:Num_R_Quad_Points, 3 )   ::  Int_Source




REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3,3) ::  Ahat_Array
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3 )  ::  f
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points )    ::  AA_Array


REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points )    ::  Cur_Val_Psi
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3 )  ::  Cur_Val_X
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3,3) ::  Cur_Drv_X

INTEGER                                                         ::  tpd, rd, td, pd, Here
INTEGER                                                         ::  i, j

INTEGER, DIMENSION(3)                                           ::  iE
INTEGER                                                         ::  Level

Level = 0
iE = [re, te, pe]


CALL Calc_Val_On_Elem_TypeA( iE, Cur_Val_Psi, iU_CF, Level )


CALL Calc_Val_And_Drv_On_Elem_TypeB( iE, DROT,      &
                                     CUR_Val_X(:,:,1),      &
                                     CUR_DRV_X(:,:,:,1),    &
                                     iU_X1, iVB_X,          &
                                     Level                  )

CALL Calc_Val_And_Drv_On_Elem_TypeB( iE, DROT,      &
                                     CUR_Val_X(:,:,2),      &
                                     CUR_DRV_X(:,:,:,2),    &
                                     iU_X2, iVB_X,          &
                                     Level                  )

CALL Calc_Val_And_Drv_On_Elem_TypeB( iE, DROT,      &
                                     CUR_Val_X(:,:,3),      &
                                     CUR_DRV_X(:,:,:,3),    &
                                     iU_X3, iVB_X,          &
                                     Level                  )

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

    tpd = Map_To_tpd(td,pd)
    Here = Quad_Map(rd,td,pd)

    Int_Source(tpd,rd,1) = GR_Source_Scalar                                 &
                            / Cur_Val_Psi(tpd,rd)                           &
                            * Block_Source_E(Here,re,te,pe)             &
                          + AA_Array(tpd,rd)                                &
                            / ( 16.0_idp * pi * Cur_Val_Psi(tpd,rd)**7)


    Int_Source(tpd,rd,2) = GR_Source_Scalar                                 &
                            / Cur_Val_Psi(tpd,rd)                           &
                            * Block_Source_E(Here,re,te,pe)

    Int_Source(tpd,rd,3) = AA_Array(tpd,rd)                                 &
                            / ( 16.0_idp * pi * Cur_Val_Psi(tpd,rd)**7)

END DO ! pd
END DO ! td
END DO ! rd


END SUBROUTINE Calc_Int_Source_In_Parts







!+204+##########################################################################!
!                                                                               !
!                 Calc_Element_ADM_Integral                                     !
!                                                                               !
!###############################################################################!
SUBROUTINE Calc_Element_ADM_Parts_Integral( Int_Weights,      &
                                            Int_Source,       &
                                            Int_Val           )

REAL(idp),  INTENT(IN), DIMENSION(1:Num_TP_Quad_Points, &
                                  1:Num_R_Quad_Points)  ::  Int_Weights

REAL(idp),  INTENT(IN), DIMENSION(1:Num_TP_Quad_Points, &
                                  1:Num_R_Quad_Points,3)::  Int_Source


REAL(idp),    DIMENSION(3), INTENT(OUT)                 ::  Int_Val

INTEGER                                                 ::  rd
COMPLEX(idp), DIMENSION(3)                              ::  Tmp_Val




Tmp_Val = 0.0_idp
DO rd = 1,NUM_R_QUAD_POINTS

    Tmp_Val(1) =  Tmp_Val(1)                    &
               + SUM( Int_Source( :, rd, 1 )    &
               * Int_Weights(:, rd )            )

    Tmp_Val(2) =  Tmp_Val(2)                    &
               + SUM( Int_Source( :, rd, 2 )    &
               * Int_Weights(:, rd )            )

    Tmp_Val(3) =  Tmp_Val(3)                    &
               + SUM( Int_Source( :, rd, 3 )    &
               * Int_Weights(:, rd )            )

END DO  ! rd Loop
    

Int_Val = REAL(Tmp_Val, Kind = idp)


END SUBROUTINE Calc_Element_ADM_Parts_Integral









END MODULE ADM_Mass_In_Parts_Module
