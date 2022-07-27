   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Functions_Jacobian                                                           !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!        Contains the functions used to calculate the components of the          !##!
!##!    extended Jacobian matrix as well as the derivative coefficients. These      !##!
!##!    are used to construct the stiffness matrix.                                 !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +301+   JCBN_kappa_FUNCTION_3D_ALL                                          !##!
!##!    +302+   JCBN_BIGK_FUNCTION                                                  !##!
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
                        ONLY : idp

USE Poseidon_Numbers_Module, &
                        ONLY : pi,                      &
                               TwoThirds,               &
                               FourThirds




IMPLICIT NONE


CONTAINS











!+301+###########################################################################!
!                                                                                !
!                                                                                !
!################################################################################!
PURE FUNCTION JCBN_kappa_FUNCTION_3D_ALL(rd, tpd,                               &
                                        NUM_R_QUAD_POINTS, NUM_TP_QUAD_POINTS,  &
                                        r, R_SQUARE, R_CUBED, RSIN_SQUARE,      &
                                        COTAN_VAL,                     &
                                        CUR_VAL_BETA, CUR_DRV_BETA              )

REAL(KIND = idp), DIMENSION(1:3,1:3)      ::  JCBN_kappa_FUNCTION_3D_ALL

INTEGER, INTENT(IN)                     ::  rd, tpd
INTEGER, INTENT(IN)                     ::  NUM_R_QUAD_POINTS, NUM_TP_QUAD_POINTS

REAL(KIND = idp), INTENT(IN)            ::  r, R_SQUARE, RSIN_SQUARE, R_CUBED,    &
                                            COTAN_VAL





REAL(KIND = idp), INTENT(IN), DIMENSION( 1:NUM_TP_QUAD_POINTS,  &
                                         1:NUM_R_QUAD_POINTS,   &
                                         1:3                    )          ::  CUR_VAL_BETA


REAL(KIND = idp), INTENT(IN), DIMENSION( 1:NUM_TP_QUAD_POINTS,  &
                                         1:NUM_R_QUAD_POINTS,   &
                                         1:3, 1:3               )          ::  CUR_DRV_BETA






JCBN_kappa_FUNCTION_3D_ALL(1,1) = FourThirds                * CUR_DRV_BETA( tpd, rd, 1, 1 )       &
                                - FourThirds/r              * CUR_VAL_BETA( tpd, rd, 1 )          &
                                - TwoThirds * COTAN_VAL     * CUR_VAL_BETA( tpd, rd, 2 )          &
                                - TwoThirds                 * CUR_DRV_BETA( tpd, rd, 2, 2 )       &
                                - TwoThirds                 * CUR_DRV_BETA( tpd, rd, 3, 3 )

JCBN_kappa_FUNCTION_3D_ALL(2,1) =  CUR_DRV_BETA( tpd, rd, 2, 1 )                                     &
                                +  CUR_DRV_BETA( tpd, rd, 1, 2 )/R_SQUARE

JCBN_kappa_FUNCTION_3D_ALL(3,1) =  CUR_DRV_BETA( tpd, rd, 3, 1 )                                     &
                                +  CUR_DRV_BETA( tpd, rd, 1, 3 )/RSIN_SQUARE





JCBN_kappa_FUNCTION_3D_ALL(1,2) = JCBN_kappa_FUNCTION_3D_ALL(2,1)

JCBN_kappa_FUNCTION_3D_ALL(2,2) = TwoThirds / R_CUBED              * CUR_VAL_BETA( tpd, rd, 1 )      &
                                -(TwoThirds / R_SQUARE)*COTAN_VAL  * CUR_VAL_BETA( tpd, rd, 2 )      &
                                + FourThirds/ R_SQUARE             * CUR_DRV_BETA( tpd, rd, 2, 2 )   &
                                - TwoThirds / R_SQUARE             * CUR_DRV_BETA( tpd, rd, 1, 1 )   &
                                - TwoThirds / R_SQUARE             * CUR_DRV_BETA( tpd, rd, 3, 3 )

JCBN_kappa_FUNCTION_3D_ALL(3,2) = CUR_DRV_BETA( tpd, rd, 3, 2 )/R_SQUARE                             &
                                + CUR_DRV_BETA( tpd, rd, 2, 3 )/RSIN_SQUARE



JCBN_kappa_FUNCTION_3D_ALL(1,3) = JCBN_kappa_FUNCTION_3D_ALL(3,1)

JCBN_kappa_FUNCTION_3D_ALL(2,3) = JCBN_kappa_FUNCTION_3D_ALL(3,2)

JCBN_kappa_FUNCTION_3D_ALL(3,3) = TwoThirds/( r * RSIN_SQUARE)          * CUR_VAL_BETA( tpd, rd, 1 )      &
                                + FourThirds * COTAN_VAL/RSIN_SQUARE    * CUR_VAL_BETA( tpd, rd, 2 )      &
                                + FourThirds / RSIN_SQUARE              * CUR_DRV_BETA( tpd, rd, 3, 3 )   &
                                - TwoThirds / RSIN_SQUARE               * CUR_DRV_BETA( tpd, rd, 1, 1 )   &
                                - TwoThirds / RSIN_SQUARE               * CUR_DRV_BETA( tpd, rd, 2, 2 )





!JCBN_kappa_Function_3D_All(1,2) = 0.0_idp
!JCBN_kappa_Function_3D_All(2,1) = 0.0_idp
!
!JCBN_kappa_Function_3D_All(3,1) = 0.0_idp
!JCBN_kappa_Function_3D_All(1,3) = 0.0_idp


!JCBN_kappa_Function_3D_All(2,3) = 0.0_idp
!JCBN_kappa_Function_3D_All(3,2) = 0.0_idp
!
!JCBN_kappa_Function_3D_All(2,2) = 0.0_idp
!JCBN_kappa_Function_3D_All(3,3) = 0.0_idp



END FUNCTION JCBN_kappa_FUNCTION_3D_ALL






!+301+###########################################################################!
!                                                                                !
!                                                                                !
!################################################################################!
SUBROUTINE JCBN_Kappa_Array_3D( Kappa_Array,                                &
                                NUM_R_QUAD_POINTS, NUM_TP_QUAD_POINTS,      &
                                Cur_R_Locs, R_SQUARE, R_CUBED, RSIN_SQUARE, &
                                COTAN_VAL,                         &
                                CUR_VAL_BETA, CUR_DRV_BETA                  )

REAL(KIND = idp), DIMENSION(1:Num_TP_Quad_Points, 1:Num_R_Quad_Points, 1:3,1:3), INTENT(OUT) ::  Kappa_Array

INTEGER, INTENT(IN)                     ::  NUM_R_QUAD_POINTS, NUM_TP_QUAD_POINTS

REAL(KIND = idp), DIMENSION(1:Num_R_Quad_Points), INTENT(IN)        ::  Cur_R_Locs
REAL(KIND = idp), DIMENSION(1:Num_R_Quad_Points), INTENT(IN)        ::  R_SQUARE
REAL(KIND = idp), DIMENSION(1:Num_R_Quad_Points), INTENT(IN)        ::  R_CUBED
REAL(KIND = idp), DIMENSION(1:Num_TP_Quad_Points), INTENT(IN)       ::  COTAN_VAL

REAL(KIND = idp), DIMENSION(1:Num_TP_Quad_Points, 1:Num_R_Quad_Points), INTENT(IN)  ::  RSIN_SQUARE


REAL(KIND = idp), INTENT(IN), DIMENSION( 1:NUM_TP_QUAD_POINTS,  &
                                         1:NUM_R_QUAD_POINTS,   &
                                         1:3                    )          ::  CUR_VAL_BETA


REAL(KIND = idp), INTENT(IN), DIMENSION( 1:NUM_TP_QUAD_POINTS,  &
                                         1:NUM_R_QUAD_POINTS,   &
                                         1:3, 1:3               )          ::  CUR_DRV_BETA


INTEGER                                                                     ::  rd, tpd


DO rd = 1,Num_R_Quad_Points
DO tpd = 1,Num_TP_Quad_Points
    
    Kappa_Array(tpd, rd, 1,1) = FourThirds                  * CUR_DRV_BETA( tpd, rd, 1, 1 )       &
                              - FourThirds/ Cur_R_Locs(rd)  * CUR_VAL_BETA( tpd, rd, 1 )          &
                              - TwoThirds * COTAN_VAL(tpd)  * CUR_VAL_BETA( tpd, rd, 2 )          &
                              - TwoThirds                   * CUR_DRV_BETA( tpd, rd, 2, 2 )       &
                              - TwoThirds                   * CUR_DRV_BETA( tpd, rd, 3, 3 )

END DO
END DO



DO rd = 1,Num_R_Quad_Points
DO tpd = 1,Num_TP_Quad_Points

    Kappa_Array(tpd, rd, 2, 1) =  CUR_DRV_BETA( tpd, rd, 1, 2 )                                     &
                                +  CUR_DRV_BETA( tpd, rd, 2, 1 )/R_SQUARE(rd)

END DO
END DO


DO rd = 1,Num_R_Quad_Points
DO tpd = 1,Num_TP_Quad_Points

    Kappa_Array(tpd, rd, 3, 1) =  CUR_DRV_BETA( tpd, rd, 1, 3 )                                     &
                                +  CUR_DRV_BETA( tpd, rd, 3, 1 )/RSIN_SQUARE(tpd,rd)
END DO
END DO




Kappa_Array(:,:,1,2) = Kappa_Array(:,:,2,1)

DO rd = 1,Num_R_Quad_Points
DO tpd = 1,Num_TP_Quad_Points
    Kappa_Array(tpd, rd, 2, 2) = TwoThirds / R_CUBED(rd)                   * CUR_VAL_BETA( tpd, rd, 1 )      &
                                - (TwoThirds/ R_SQUARE(rd))*COTAN_VAL(tpd)  * CUR_VAL_BETA( tpd, rd, 2 )      &
                                + FourThirds/ R_SQUARE(rd)                  * CUR_DRV_BETA( tpd, rd, 2, 2 )   &
                                - TwoThirds / R_SQUARE(rd)                  * CUR_DRV_BETA( tpd, rd, 1, 1 )   &
                                - TwoThirds / R_SQUARE(rd)                  * CUR_DRV_BETA( tpd, rd, 3, 3 )
END DO
END DO


DO rd = 1,Num_R_Quad_Points
DO tpd = 1,Num_TP_Quad_Points
    Kappa_Array(tpd, rd, 3, 2) = CUR_DRV_BETA( tpd, rd, 2, 3 )/R_SQUARE(rd)                             &
                                + CUR_DRV_BETA( tpd, rd, 3, 2 )/RSIN_SQUARE(tpd,rd)
END DO
END DO


Kappa_Array(:,:,1,3) = Kappa_Array(:,:,3,1)
Kappa_Array(:,:,2,3) = Kappa_Array(:,:,3,2)



DO rd = 1,Num_R_Quad_Points
DO tpd = 1,Num_TP_Quad_Points
    Kappa_Array(tpd, rd, 3, 3) = TwoThirds/( Cur_R_Locs(rd) * RSIN_SQUARE(tpd,rd))                 &
                                                                    * CUR_VAL_BETA( tpd, rd, 1 )    &
                                + FourThirds * COTAN_VAL(tpd)/RSIN_SQUARE(tpd,rd)                   &
                                                                    * CUR_VAL_BETA( tpd, rd, 2 )    &
                                + FourThirds / RSIN_SQUARE(tpd,rd)  * CUR_DRV_BETA( tpd, rd, 3, 3 ) &
                                - TwoThirds / RSIN_SQUARE(tpd,rd)   * CUR_DRV_BETA( tpd, rd, 1, 1 ) &
                                - TwoThirds / RSIN_SQUARE(tpd,rd)   * CUR_DRV_BETA( tpd, rd, 2, 2 )
END DO
END DO

END SUBROUTINE JCBN_Kappa_Array_3D




!+302+###########################################################################!
!                                                                                !
!                                                                                !
!################################################################################!
PURE FUNCTION JCBN_kappa_FUNCTION_1D_ALL(rd, r, CUR_VAL_BETA, CUR_DRV_BETA, NUM_R_QUAD_POINTS )

REAL(KIND = idp), DIMENSION(1:3,1:3)                                ::  JCBN_kappa_FUNCTION_1D_ALL

INTEGER, INTENT(IN)                                                 ::  rd
REAL(KIND = idp), INTENT(IN)                                        ::  r

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:NUM_R_QUAD_POINTS )      ::  CUR_VAL_BETA
REAL(KIND = idp), INTENT(IN), DIMENSION( 1:NUM_R_QUAD_POINTS )      ::  CUR_DRV_BETA

INTEGER, INTENT(IN)                     ::  NUM_R_QUAD_POINTS


JCBN_kappa_FUNCTION_1D_ALL = FourThirds     * CUR_DRV_BETA( rd )       &
                           - FourThirds/r   * CUR_VAL_BETA( rd )


END FUNCTION JCBN_kappa_FUNCTION_1D_ALL











!+401+###########################################################################!
!                                                                                !
!                                                                                !
!################################################################################!
PURE FUNCTION JCBN_BIGK_FUNCTION(   rd, tpd,                                    &
                                    NUM_R_QUAD_POINTS, NUM_TP_QUAD_POINTS,      &
                                    CUR_VAL_BETA, CUR_DRV_BETA,                 &
                                    R_VAL, R_SQUARE, SIN_SQUARE,  CSC_SQUARE,   &
                                    RSIN_SQUARE, COTAN_VAL                      )

REAL(KIND = idp)                    ::  JCBN_BIGK_FUNCTION

INTEGER, INTENT(IN)                 ::  rd, tpd
INTEGER, INTENT(IN)                     ::  NUM_R_QUAD_POINTS, NUM_TP_QUAD_POINTS

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:NUM_TP_QUAD_POINTS,  &
                                         1:NUM_R_QUAD_POINTS,   &
                                         1:3                    )   ::  CUR_VAL_BETA



REAL(KIND = idp), INTENT(IN), DIMENSION( 1:NUM_TP_QUAD_POINTS,  &
                                         1:NUM_R_QUAD_POINTS,   &
                                         1:3,                   &
                                         1:3                    )   ::  CUR_DRV_BETA



REAL(KIND = idp), INTENT(IN)                                    ::  R_VAL,                  &
                                                                    R_SQUARE,               &
                                                                    SIN_SQUARE,             &
                                                                    CSC_SQUARE,             &
                                                                    RSIN_SQUARE,            &
                                                                    COTAN_VAL




REAL(KIND = idp)                    ::  EightThirds





EightThirds = 2.0_idp*FourThirds

JCBN_BIGK_FUNCTION =                                                                                    &
        FourThirds                      * CUR_DRV_BETA(tpd, rd, 1,1)    * CUR_DRV_BETA(tpd, rd, 1,1)    &
      + R_SQUARE                        * CUR_DRV_BETA(tpd, rd, 2,1)    * CUR_DRV_BETA(tpd, rd, 2,1)    &
      + 1.0_idp/R_SQUARE                * CUR_DRV_BETA(tpd, rd, 1,2)    * CUR_DRV_BETA(tpd, rd, 1,2)    &
      + 2.0_idp                         * CUR_DRV_BETA(tpd, rd, 2,1)    * CUR_DRV_BETA(tpd, rd, 1,2)    &
      - FourThirds                      * CUR_DRV_BETA(tpd, rd, 1,1)    * CUR_DRV_BETA(tpd, rd, 2,2)    &
      + FourThirds                      * CUR_DRV_BETA(tpd, rd, 2,2)    * CUR_DRV_BETA(tpd, rd, 2,2)    &
      + SIN_SQUARE                      * CUR_DRV_BETA(tpd, rd, 3,2)    * CUR_DRV_BETA(tpd, rd, 3,2)    &
      + CSC_SQUARE                      * CUR_DRV_BETA(tpd, rd, 2,3)    * CUR_DRV_BETA(tpd, rd, 2,3)    &
      + 2.0_idp                         * CUR_DRV_BETA(tpd, rd, 3,2)    * CUR_DRV_BETA(tpd, rd, 2,3)    &
      - FourThirds                      * CUR_DRV_BETA(tpd, rd, 2,2)    * CUR_DRV_BETA(tpd, rd, 3,3)    &
      + FourThirds                      * CUR_DRV_BETA(tpd, rd, 3,3)    * CUR_DRV_BETA(tpd, rd, 3,3)    &
      + RSIN_SQUARE                     * CUR_DRV_BETA(tpd, rd, 3,1)    * CUR_DRV_BETA(tpd, rd, 3,1)    &
      + 1.0_idp/RSIN_SQUARE             * CUR_DRV_BETA(tpd, rd, 1,3)    * CUR_DRV_BETA(tpd, rd, 1,3)    &
      + 2.0_idp                         * CUR_DRV_BETA(tpd, rd, 3,1)    * CUR_DRV_BETA(tpd, rd, 1,3)    &
      - FourThirds                      * CUR_DRV_BETA(tpd, rd, 1,1)    * CUR_DRV_BETA(tpd, rd, 3,3)    &
      +(- EightThirds/R_VAL             * CUR_DRV_BETA(tpd, rd, 1,1)                                    &
      + FourThirds/R_VAL                * CUR_DRV_BETA(tpd, rd, 2,2)                                    &
      + FourThirds/R_VAL                * CUR_DRV_BETA(tpd, rd, 3,3) )  * CUR_VAL_BETA(tpd, rd, 1)      &
      + (- FourThirds*COTAN_VAL         * CUR_DRV_BETA(tpd, rd, 1,1)                                    &
      - FourThirds*COTAN_VAL            * CUR_DRV_BETA(tpd, rd, 2,2)                                    &
      + EightThirds*COTAN_VAL           * CUR_DRV_BETA(tpd, rd, 3,3) )  * CUR_VAL_BETA(tpd, rd, 2)      &
      + FourThirds/R_SQUARE             * CUR_VAL_BETA(tpd, rd, 1)      * CUR_VAL_BETA(tpd, rd, 1)      &
      + FourThirds*COTAN_VAL/R_VAL      * CUR_VAL_BETA(tpd, rd, 1)      * CUR_VAL_BETA(tpd, rd, 2)      &
      + FourThirds*COTAN_VAL*COTAN_VAL  * CUR_VAL_BETA(tpd, rd, 2)      * CUR_VAL_BETA(tpd, rd, 2)


END FUNCTION JCBN_BIGK_FUNCTION




!+401+###########################################################################!
!                                                                                !
!                                                                                !
!################################################################################!
PURE FUNCTION JCBN_BIGK_FUNCTION_1D( rd, CUR_VAL_BETA, CUR_DRV_BETA,        &
                                     R_VAL, R_SQUARE,                         &
                                     NUM_R_QUAD_POINTS )

REAL(KIND = idp)                    ::  JCBN_BIGK_FUNCTION_1D

INTEGER, INTENT(IN)                 ::  rd

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:NUM_R_QUAD_POINTS )      ::  CUR_VAL_BETA
REAL(KIND = idp), INTENT(IN), DIMENSION( 1:NUM_R_QUAD_POINTS )      ::  CUR_DRV_BETA

REAL(KIND = idp), INTENT(IN)                                        ::  R_VAL, R_SQUARE

INTEGER, INTENT(IN)                     ::  NUM_R_QUAD_POINTS


REAL(KIND = idp)                    ::  EightThirds

EightThirds = 2.0_idp*FourThirds

JCBN_BIGK_FUNCTION_1D =                                                     &
        FourThirds          * CUR_DRV_BETA( rd )    * CUR_DRV_BETA( rd )    &
      - EightThirds/R_VAL   * CUR_DRV_BETA( rd )    * CUR_VAL_BETA( rd )    &
      + FourThirds/R_SQUARE * CUR_VAL_BETA( rd )    * CUR_VAL_BETA( rd )




END FUNCTION JCBN_BIGK_FUNCTION_1D












!+301+###########################################################################!
!                                                                                !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Ahat( Ahat_Array,                               &
                        NUM_R_QUAD_POINTS, NUM_TP_QUAD_POINTS,  &
                        Cur_R_Locs, R_SQUARE,                   &
                        RSIN_SQUARE, COTAN_VAL,                 &
                        CUR_VAL_X, CUR_DRV_X                    )

REAL(KIND = idp), DIMENSION(Num_TP_Quad_Points,Num_R_Quad_Points,3,3), INTENT(OUT) ::  Ahat_Array

INTEGER, INTENT(IN)                                                 ::  NUM_R_QUAD_POINTS,      &
                                                                        NUM_TP_QUAD_POINTS

REAL(KIND = idp), DIMENSION(Num_R_Quad_Points), INTENT(IN)          ::  Cur_R_Locs
REAL(KIND = idp), DIMENSION(Num_R_Quad_Points), INTENT(IN)          ::  R_SQUARE

REAL(KIND = idp), DIMENSION(Num_TP_Quad_Points), INTENT(IN)         ::  COTAN_VAL

REAL(KIND = idp), DIMENSION(1:Num_TP_Quad_Points, 1:Num_R_Quad_Points), INTENT(IN)  ::  RSIN_SQUARE


REAL(KIND = idp), INTENT(IN), DIMENSION( 1:NUM_TP_QUAD_POINTS,  &
                                         1:NUM_R_QUAD_POINTS,   &
                                         1:3                    )          ::  CUR_VAL_X


REAL(KIND = idp), INTENT(IN), DIMENSION( 1:NUM_TP_QUAD_POINTS,  &
                                         1:NUM_R_QUAD_POINTS,   &
                                         1:3, 1:3               )          ::  CUR_DRV_X


INTEGER                                                                     ::  rd, tpd





DO rd = 1,Num_R_Quad_Points
DO tpd = 1,Num_TP_Quad_Points
    
        Ahat_Array(tpd, rd, 1,1) =  FourThirds                  * CUR_DRV_X( tpd, rd, 1, 1 )       &
                                  - FourThirds/ Cur_R_Locs(rd)  * CUR_VAL_X( tpd, rd, 1 )          &
                                  - TwoThirds * COTAN_VAL(tpd)  * CUR_VAL_X( tpd, rd, 2 )          &
                                  - TwoThirds                   * CUR_DRV_X( tpd, rd, 2, 2 )       &
                                  - TwoThirds                   * CUR_DRV_X( tpd, rd, 3, 3 )

END DO
END DO



DO rd = 1,Num_R_Quad_Points
DO tpd = 1,Num_TP_Quad_Points

!    PRINT*,CUR_DRV_X( tpd, rd, 1, 2 ),CUR_DRV_X( tpd, rd, 2, 1 ),R_SQUARE(rd)

    Ahat_Array(tpd, rd, 2, 1) =  CUR_DRV_X( tpd, rd, 1, 2 )                                     &
                                +  CUR_DRV_X( tpd, rd, 2, 1 )/R_SQUARE(rd)

END DO
END DO


DO rd = 1,Num_R_Quad_Points
DO tpd = 1,Num_TP_Quad_Points

    Ahat_Array(tpd, rd, 3, 1) =  CUR_DRV_X( tpd, rd, 1, 3 )                                     &
                                +  CUR_DRV_X( tpd, rd, 3, 1 )/RSIN_SQUARE(tpd,rd)
END DO
END DO




Ahat_Array(:,:,1,2) = Ahat_Array(:,:,2,1)

DO rd = 1,Num_R_Quad_Points
DO tpd = 1,Num_TP_Quad_Points
    Ahat_Array(tpd, rd, 2, 2) = TwoThirds    /(R_SQUARE(rd)*Cur_R_Locs(rd)) * CUR_VAL_X( tpd, rd, 1 )      &
                                - (TwoThirds/ R_SQUARE(rd))*COTAN_VAL(tpd)  * CUR_VAL_X( tpd, rd, 2 )      &
                                + FourThirds/ R_SQUARE(rd)                  * CUR_DRV_X( tpd, rd, 2, 2 )   &
                                - TwoThirds / R_SQUARE(rd)                  * CUR_DRV_X( tpd, rd, 1, 1 )   &
                                - TwoThirds / R_SQUARE(rd)                  * CUR_DRV_X( tpd, rd, 3, 3 )
END DO
END DO


DO rd = 1,Num_R_Quad_Points
DO tpd = 1,Num_TP_Quad_Points
    Ahat_Array(tpd, rd, 3, 2) = CUR_DRV_X( tpd, rd, 2, 3 )/R_SQUARE(rd)                             &
                                + CUR_DRV_X( tpd, rd, 3, 2 )/RSIN_SQUARE(tpd,rd)
END DO
END DO


Ahat_Array(:,:,1,3) = Ahat_Array(:,:,3,1)
Ahat_Array(:,:,2,3) = Ahat_Array(:,:,3,2)



DO rd = 1,Num_R_Quad_Points
DO tpd = 1,Num_TP_Quad_Points
    Ahat_Array(tpd, rd, 3, 3) = TwoThirds/( Cur_R_Locs(rd) * RSIN_SQUARE(tpd,rd))                   &
                                                                    * CUR_VAL_X( tpd, rd, 1 )    &
                                + FourThirds * COTAN_VAL(tpd)/RSIN_SQUARE(tpd,rd)                   &
                                                                    * CUR_VAL_X( tpd, rd, 2 )    &
                                + FourThirds / RSIN_SQUARE(tpd,rd)  * CUR_DRV_X( tpd, rd, 3, 3 ) &
                                - TwoThirds / RSIN_SQUARE(tpd,rd)   * CUR_DRV_X( tpd, rd, 1, 1 ) &
                                - TwoThirds / RSIN_SQUARE(tpd,rd)   * CUR_DRV_X( tpd, rd, 2, 2 )
END DO
END DO

END SUBROUTINE Calc_Ahat





END MODULE Functions_Jacobian
