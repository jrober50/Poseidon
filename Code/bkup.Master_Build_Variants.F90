

!+203b+###########################################################################!
!                                                                                !
!                  Calc_3D_SubJcbn_Terms                                         !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_3D_SubJcbn_Terms_B(   re, te, pe,                                     &
                                    CUR_VAL_PSI, CUR_DRV_PSI,                       &
                                    CUR_VAL_ALPHAPSI, CUR_DRV_ALPHAPSI,             &
                                    CUR_VAL_BETA, CUR_DRV_BETA, CUR_DDRV_BETA,      &
                                    RHS_TERMS,                                      &
                                    SUBJCBN_PSI_TERMS_B, SUBJCBN_ALPHAPSI_TERMS_B,      &
                                    SUBJCBN_BETA1_TERMS_B, SUBJCBN_BETA2_TERMS_B,       &
                                    SUBJCBN_BETA3_TERMS_B,                            &
                                    CUR_R_LOCS, CUR_T_LOCS, CUR_P_LOCS,             &
                                    R_SQUARE, R_CUBE, R_INVERSE,                    &
                                    SIN_VAL, COS_VAL, COTAN_VAL,                    &
                                    SIN_SQUARE, RSIN_SQUARE                         )




INTEGER, INTENT(IN)                                                     ::  re, te, pe

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:NUM_R_QUAD_POINTS,   &
                                         1:NUM_T_QUAD_POINTS,   &
                                         1:NUM_P_QUAD_POINTS    )       ::  CUR_VAL_PSI

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:NUM_R_QUAD_POINTS,   &
                                         1:NUM_T_QUAD_POINTS,   &
                                         1:NUM_P_QUAD_POINTS    )       ::  CUR_VAL_ALPHAPSI


REAL(KIND = idp), INTENT(IN), DIMENSION( 1:3,                   &
                                         1:NUM_R_QUAD_POINTS,   &
                                         1:NUM_T_QUAD_POINTS,   &
                                         1:NUM_P_QUAD_POINTS    )       ::  CUR_VAL_BETA


REAL(KIND = idp), INTENT(IN), DIMENSION( 1:3,                   &
                                         1:NUM_R_QUAD_POINTS,   &
                                         1:NUM_T_QUAD_POINTS,   &
                                         1:NUM_P_QUAD_POINTS    )       ::  CUR_DRV_PSI

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:3,                   &
                                         1:NUM_R_QUAD_POINTS,   &
                                         1:NUM_T_QUAD_POINTS,   &
                                         1:NUM_P_QUAD_POINTS    )       ::  CUR_DRV_ALPHAPSI

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:3, 1:3,              &
                                         1:NUM_R_QUAD_POINTS,   &
                                         1:NUM_T_QUAD_POINTS,   &
                                         1:NUM_P_QUAD_POINTS    )       ::  CUR_DRV_BETA

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:3, 1:3, 1:3,         &
                                         1:NUM_R_QUAD_POINTS,   &
                                         1:NUM_T_QUAD_POINTS,   &
                                         1:NUM_P_QUAD_POINTS    )       ::  CUR_DDRV_BETA





REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_R_QUAD_POINTS)            ::  CUR_R_LOCS
REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_T_QUAD_POINTS)            ::  CUR_T_LOCS
REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_P_QUAD_POINTS)            ::  CUR_P_LOCS


REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_R_QUAD_POINTS,    &
                                        1:NUM_T_QUAD_POINTS     )       ::  RSIN_SQUARE

REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_R_QUAD_POINTS)            ::  R_SQUARE,           &
                                                                            R_CUBE,             &
                                                                            R_INVERSE


REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_T_QUAD_POINTS)            ::  SIN_VAL,            &
                                                                            SIN_SQUARE,         &
                                                                            COTAN_VAL,          &
                                                                            COS_VAL



REAL(KIND = idp), INTENT(INOUT), DIMENSION( 1:5,                        &
                                            1:NUM_R_QUAD_POINTS,        &
                                            1:NUM_T_QUAD_POINTS,        &
                                            1:NUM_P_QUAD_POINTS         )   ::  RHS_TERMS


REAL(KIND = idp), INTENT(INOUT), DIMENSION( NUM_QUAD_DOF, 1:14 )            ::  SUBJCBN_PSI_TERMS_B

REAL(KIND = idp), INTENT(INOUT), DIMENSION( NUM_QUAD_DOF, 1:14 )            ::  SUBJCBN_ALPHAPSI_TERMS_B

REAL(KIND = idp), INTENT(INOUT), DIMENSION( NUM_QUAD_DOF, 1:20 )            ::  SUBJCBN_BETA1_TERMS_B

REAL(KIND = idp), INTENT(INOUT), DIMENSION( NUM_QUAD_DOF, 1:20 )            ::  SUBJCBN_BETA2_TERMS_B

REAL(KIND = idp), INTENT(INOUT), DIMENSION( NUM_QUAD_DOF, 1:20 )            ::  SUBJCBN_BETA3_TERMS_B


REAL(KIND = idp)                                                            ::  REUSED_VALUE

INTEGER                                                                     ::  pd, td, rd,     &
                                                                                i,qd


REAL(KIND = idp), DIMENSION(1:8)                                            ::  PSI_POWER
REAL(KIND = idp), DIMENSION(1:4)                                            ::  ALPHAPSI_POWER


REAL(KIND = idp), DIMENSION(1:3,1:3)                                        ::  JCBN_mu_Array
REAL(KIND = idp), DIMENSION(1:3)                                            ::  JCBN_n_ARRAY

REAL(KIND = idp)                                                            ::  JCBN_BIGK_VALUE












!$OMP PARALLEL DEFAULT(none)                                            &
!$OMP PRIVATE( pd, td, rd, qd,                                          &
!$OMP           PSI_POWER, ALPHAPSI_POWER,                              &
!$OMP           JCBN_BIGK_VALUE, JCBN_mu_ARRAY,JCBN_n_ARRAY,            &
!$OMP           REUSED_VALUE                                        )   &
!$OMP SHARED( re, te, pe,                                               &
!$OMP           CUR_VAL_PSI, CUR_VAL_ALPHAPSI, CUR_VAL_BETA,            &
!$OMP           CUR_DRV_PSI, CUR_DRV_ALPHAPSI, CUR_DRV_BETA,            &
!$OMP           CUR_DDRV_BETA,                                          &
!$OMP           CUR_R_LOCS, CUR_T_LOCS, CUR_P_LOCS,                     &
!$OMP           RSIN_SQUARE, R_SQUARE, R_CUBE, R_INVERSE,               &
!$OMP           SIN_VAL, SIN_SQUARE, COTAN_VAL, COS_VAL,                &
!$OMP           NUM_OFF_DIAGONALS,                                      &
!$OMP           RHS_TERMS,                                              &
!$OMP           SUBJCBN_PSI_TERMS_B, SUBJCBN_ALPHAPSI_TERMS_B,              &
!$OMP           SUBJCBN_BETA1_TERMS_B, SUBJCBN_BETA2_TERMS_B,               &
!$OMP           SUBJCBN_BETA3_TERMS_B,                                    &
!$OMP           Block_SOURCE_E, Block_SOURCE_S, Block_SOURCE_Si,        &
!$OMP           OneThird, OneThirtySecond, FourThirds, TwoThirds,       &
!$OMP           LM_Length                           )



!$OMP DO SCHEDULE(dynamic), COLLAPSE(3)
DO pd = 1,NUM_P_QUAD_POINTS
 DO td = 1,NUM_T_QUAD_POINTS
  DO rd = 1,NUM_R_QUAD_POINTS

        qd = ((pd-1) * NUM_T_QUAD_POINTS + td-1 ) * NUM_R_QUAD_POINTS + rd


        PSI_POWER(1) = CUR_VAL_PSI(rd, td, pd)
        DO i = 2,8

            PSI_POWER(i) = PSI_POWER(i-1)*PSI_POWER(1)

        END DO


        ALPHAPSI_POWER(1) = CUR_VAL_ALPHAPSI(rd, td, pd)
        DO i = 2,4

            ALPHAPSI_POWER(i) = ALPHAPSI_POWER(i-1)*ALPHAPSI_POWER(1)

        END DO




        ! K_{ij}K^{ij} = Psi^{14}/AlphaPsi^{2} * BIGK
        JCBN_BIGK_VALUE = JCBN_BIGK_FUNCTION(rd, td, pd,                                        &
                                                CUR_VAL_BETA, CUR_DRV_BETA,                     &
                                                CUR_R_LOCS(rd), R_SQUARE(rd), SIN_SQUARE(td),   &
                                                RSIN_SQUARE(rd, td), COTAN_VAL(td)              )




        JCBN_mu_Array = JCBN_mu_FUNCTION_3D_ALL( rd, td, pd,                                &
                                                CUR_R_LOCS(rd), R_SQUARE(rd), R_CUBE(rd),   &
                                                R_INVERSE(rd), RSIN_SQUARE(rd, td),         &
                                                SIN_VAL(td), SIN_SQUARE(td),                &
                                                COS_VAL(td), COTAN_VAL(td),                 &
                                                CUR_VAL_BETA, CUR_DRV_BETA                          )





        JCBN_n_ARRAY(:) = CUR_DRV_ALPHAPSI(:, rd, td, pd) / ALPHAPSI_POWER(1)   &
                            - 7 * CUR_DRV_PSI(:, rd, td, pd )/ PSI_POWER(1)











        RHS_Terms(1, rd, td, pd) = - TwoPi * Block_Source_E(rd, td, pd, re, te, pe) * PSI_POWER(5)            &
                                 - PSI_POWER(7)/ (16.0_idp * ALPHAPSI_POWER(2)) * JCBN_BIGK_VALUE



        RHS_Terms(2, rd, td, pd) = TwoPi * ALPHAPSI_POWER(1) * PSI_POWER(4)                                     &
                        * ( Block_Source_E(rd, td, pd, re, te, pe) + 2.0_idp * Block_Source_S(rd, td, pd, re, te, pe)  )    &
                        + 7.0_idp*PSI_POWER(6)/ (16.0_idp * ALPHAPSI_POWER(1)) * JCBN_BIGK_VALUE






        RHS_Terms(3, rd, td, pd) = 16.0_idp * pi                                                &
                          * ALPHAPSI_POWER(1)                                                   &
                          * PSI_POWER(3)                                                        &
                          * Block_Source_Si(rd, td, pd, re, te, pe, 1)                          &
                     - OneThird                                                                 &
                          * ( CUR_DDRV_BETA(1, 1, 1, rd, td, pd )                               &
                            + CUR_DDRV_BETA(1, 2, 2, rd, td, pd )                               &
                            + CUR_DDRV_BETA(1, 3, 3, rd, td, pd )                               &
                            + 2.0_idp/CUR_R_LOCS(rd) * CUR_DRV_BETA(1, 1, rd, td, pd)           &
                            + COTAN_VAL(td) * CUR_DRV_BETA(1, 2, rd, td, pd)                    &
                            + COTAN_VAL(td) * CUR_DRV_BETA(1, 3, rd, td, pd)                    &
                            - 2.0_idp /R_SQUARE(rd) * CUR_VAL_BETA(1, rd, td, pd )              &
                            )                                                                   &
                      + JCBN_n_ARRAY(1)* JCBN_mu_Array(1,1)                                     &
                      + JCBN_n_ARRAY(2)* JCBN_mu_Array(2,1)                                     &
                      + JCBN_n_ARRAY(3)* JCBN_mu_Array(3,1)





        RHS_Terms(4, rd, td, pd) =  16.0_idp * pi                                               &
                           * ALPHAPSI_POWER(1)                                                  &
                           * PSI_POWER(3)                                                       &
                           * Block_Source_Si(rd, td, pd, re, te, pe, 2)                         &
                        - OneThird                                                              &
                            * ( R_SQUARE(rd)                                                    &
                             * ( CUR_DDRV_BETA(2, 1, 1, rd, td, pd )                            &
                               + CUR_DDRV_BETA(2, 2, 2, rd, td, pd )                            &
                               + CUR_DDRV_BETA(2, 3, 3, rd, td, pd )                            &
                               + 2.0_idp/CUR_R_LOCS(rd) * CUR_DRV_BETA(2, 1, rd, td, pd)        &
                               + COTAN_VAL(td)                                                  &
                                * ( CUR_DRV_BETA(2, 2, rd, td, pd)                              &
                                    + CUR_DRV_BETA(2, 3, rd, td, pd)  )                         &
                               - ( 1.0_idp / SIN_SQUARE(td) )                                   &
                                * ( CUR_VAL_BETA(2, rd, td, pd)                                 &
                                    - CUR_VAL_BETA(3, rd, td, pd)     )                         &
                               )                                                                &
                             )                                                                  &
                        + JCBN_n_ARRAY(1) * JCBN_mu_Array(1,2)                                  &
                        + JCBN_n_ARRAY(2) * JCBN_mu_Array(2,2)                                  &
                        + JCBN_n_ARRAY(3) * JCBN_mu_Array(3,2)





        RHS_Terms(5, rd, td, pd) = 16.0_idp * pi                                            &
                          * ALPHAPSI_POWER(1)                                               &
                          * PSI_POWER(3)                                                    &
                          * Block_Source_Si(rd, td, pd, re, te, pe, 3)                      &
                        - OneThird                                                          &
                            * ( RSIN_SQUARE(rd, td)                                         &
                             * ( CUR_DDRV_BETA(3, 1, 1, rd, td, pd )                        &
                               + CUR_DDRV_BETA(3, 2, 2, rd, td, pd )                        &
                               + CUR_DDRV_BETA(3, 3, 3, rd, td, pd )                        &
                               + 2.0_idp/CUR_R_LOCS(rd) * CUR_DRV_BETA(3, 1, rd, td, pd)    &
                               + COTAN_VAL(td) * CUR_DRV_BETA(3, 2, rd, td, pd)             &
                               + COTAN_VAL(td) * CUR_DRV_BETA(3, 3, rd, td, pd)             &
                               )                                                            &
                             )                                                              &
                        + JCBN_n_ARRAY(1) * JCBN_mu_Array(1,3)                              &
                        + JCBN_n_ARRAY(2) * JCBN_mu_Array(2,3)                              &
                        + JCBN_n_ARRAY(3) * JCBN_mu_Array(3,3)







        !
        !   With Respect to Psi Terms Jacobian Terms
        !

        ! J_{1,1lmn}
        SUBJCBN_PSI_TERMS_B(qd, 1 ) = 10.0_idp*pi                                                  &
                                                * Block_Source_E(rd, td, pd, re, te, pe)* PSI_POWER(4)  &
                                         + (7.0_idp/16.0_idp)                                           &
                                                * PSI_POWER(6)/ALPHAPSI_POWER(2)                        &
                                                * JCBN_BIGK_VALUE

        ! J_{2,1lmn}
        SUBJCBN_PSI_TERMS_B(qd, 2 ) = -8.0_idp * pi * ALPHAPSI_POWER(1) * PSI_POWER(3)             &
                                            * ( Block_Source_E(rd, td, pd, re, te, pe)                  &
                                                + 2 * Block_Source_S(rd, td, pd, re, te, pe) )          &
                                            - (42.0_idp / 16.0_idp )                                    &
                                                * PSI_POWER(5)/ALPHAPSI_POWER(1)                        &
                                                * JCBN_BIGK_VALUE

        ! J_{3,1lmn} Non-Derivative Term
        SUBJCBN_PSI_TERMS_B(qd, 3 ) = ( 7.0_idp / PSI_POWER(2) )                                   &
                                                * (CUR_DRV_PSI(1, rd, td, pd )*JCBN_mu_ARRAY(1,1)       &
                                                  + CUR_DRV_PSI(2, rd, td, pd )*JCBN_mu_ARRAY(2,1)      &
                                                  + CUR_DRV_PSI(3, rd, td, pd )*JCBN_mu_ARRAY(3,1) )    &
                                            - 48.0_idp * pi                                             &
                                                * Block_Source_Si(rd, td, pd, re, te, pe, 1)            &
                                                * ALPHAPSI_POWER(1)                                     &
                                                * PSI_POWER(2)


        ! J_{3,1lmn} Derivative Terms
        SUBJCBN_PSI_TERMS_B(qd, 4:6 ) = (7.0_idp/ PSI_POWER(1) )*JCBN_mu_ARRAY(1:3,1)



        ! J_{4,1lmn} Non-Derivative Term
        SUBJCBN_PSI_TERMS_B(qd, 7 ) = (-7.0_idp / PSI_POWER(2) )                                   &
                                                * (CUR_DRV_PSI(1, rd, td, pd)*JCBN_mu_ARRAY(1,2)        &
                                                 + CUR_DRV_PSI(2, rd, td, pd)*JCBN_mu_ARRAY(2,2)        &
                                                 + CUR_DRV_PSI(3, rd, td, pd)*JCBN_mu_ARRAY(3,2) )      &
                                            - 48.0_idp * pi                                             &
                                                * Block_Source_Si(rd, td, pd, re, te, pe, 2)            &
                                                * ALPHAPSI_POWER(1)                                     &
                                                * PSI_POWER(2)

        ! J_{4,1lmn} Derivative Terms
        SUBJCBN_PSI_TERMS_B(qd, 8:10 ) = (7.0_idp/ PSI_POWER(1) )*JCBN_mu_ARRAY(1:3,2)






        ! J_{5,1lmn} Non-Derivative Term
        SUBJCBN_PSI_TERMS_B(qd, 11 ) = (-7.0_idp / PSI_POWER(2) )                  &
                                * (CUR_DRV_PSI(1, rd, td, pd)*JCBN_mu_ARRAY(1,3)        &
                                 + CUR_DRV_PSI(2, rd, td, pd)*JCBN_mu_ARRAY(2,3)        &
                                 + CUR_DRV_PSI(3, rd, td, pd)*JCBN_mu_ARRAY(3,3) )      &
                            - 48.0_idp * pi                                             &
                                * Block_Source_Si(rd, td, pd, re, te, pe, 3)            &
                                * ALPHAPSI_POWER(1)                                     &
                                * PSI_POWER(2)

        ! J_{5,1lmn} Derivative Terms
        SUBJCBN_PSI_TERMS_B(qd, 12:14 ) = (7.0_idp/ PSI_POWER(1) )*JCBN_mu_ARRAY(1:3,3)






        !
        ! Alpha Psi Terms
        !


        ! J_{1,2lmn}
        SUBJCBN_ALPHAPSI_TERMS_B(qd, 1 ) = -PSI_POWER(7)/(8.0_idp*ALPHAPSI_POWER(3))               &
                                                * JCBN_BIGK_VALUE

        ! J_{2,2lmn}
        SUBJCBN_ALPHAPSI_TERMS_B(qd, 2 ) =  -TwoPi * PSI_POWER(4)                                  &
                                                *( Block_Source_E(rd, td, pd, re, te, pe)               &
                                                    + 2.0_idp * Block_Source_S(rd, td, pd, re, te, pe) )      &
                                                - (7.0_idp * OneThirtySecond)                           &
                                                    * PSI_POWER(3)/ALPHAPSI_POWER(3)                    &
                                                    * JCBN_BIGK_VALUE

        ! J_{3,2lmn} Non-Derivative Term
        SUBJCBN_ALPHAPSI_TERMS_B(qd, 3 ) = -(CUR_DRV_ALPHAPSI(1, rd, td, pd )*JCBN_mu_Array(1,1)       &
                                                + CUR_DRV_ALPHAPSI(2, rd, td, pd )*JCBN_mu_Array(2,1)       &
                                                + CUR_DRV_ALPHAPSI(3, rd, td, pd )*JCBN_mu_Array(3,1)   )   &
                                                /ALPHAPSI_POWER(1)                                          &
                                                -16.0_idp * pi * Block_Source_Si(rd,td,pd,re,te,pe,1)       &
                                                    * PSI_POWER(3)


        ! J_{3,2lmn} Derivative Terms
        SUBJCBN_ALPHAPSI_TERMS_B(qd, 4:6 ) = (-1.0_idp/ALPHAPSI_POWER(1) )*JCBN_mu_ARRAY(1:3,1)




        ! J_{4,2lmn} Non-Derivative Term
        SUBJCBN_ALPHAPSI_TERMS_B(qd, 7 ) = -( CUR_DRV_ALPHAPSI(1, rd, td, pd )*JCBN_mu_Array(1,2)      &
                                                + CUR_DRV_ALPHAPSI(2, rd, td, pd )*JCBN_mu_Array(2,2)       &
                                                + CUR_DRV_ALPHAPSI(3, rd, td, pd )*JCBN_mu_Array(3,2)   )   &
                                                /ALPHAPSI_POWER(1)                                          &
                                                -16.0_idp * pi * Block_Source_Si(rd, td, pd, re, te, pe, 2) &
                                                    * PSI_POWER(3)


        ! J_{4,2lmn} Derivative Terms
        SUBJCBN_ALPHAPSI_TERMS_B(qd, 8:10 ) = (-1.0_idp/ ALPHAPSI_POWER(1) )*JCBN_mu_ARRAY(1:3,2)



        ! J_{5,2lmn} Non-Derivative Term
        SUBJCBN_ALPHAPSI_TERMS_B(qd, 11 ) = -(CUR_DRV_ALPHAPSI(1, rd, td, pd ) * JCBN_mu_Array(1,3)    &
                                                 + CUR_DRV_ALPHAPSI(2, rd, td, pd ) * JCBN_mu_Array(2,3)    &
                                                 + CUR_DRV_ALPHAPSI(3, rd, td, pd ) * JCBN_mu_Array(3,3)  ) &
                                                /ALPHAPSI_POWER(1)                                          &
                                                -16.0_idp * pi * Block_Source_Si(rd,td,pd,re,te,pe,3)       &
                                                    * PSI_POWER(3)

        ! J_{5,2lmn} Derivative Terms
        SUBJCBN_ALPHAPSI_TERMS_B(qd, 12:14 ) = (-1.0_idp/ ALPHAPSI_POWER(1) )*JCBN_mu_ARRAY(1:3,3)






        !
        !   Beta 1 Terms
        !




        REUSED_VALUE = (1.0_idp/16_idp) * ( PSI_POWER(7) / ALPHAPSI_POWER(2) )


        ! J_{1,3lmn}
        ! Reused_Value * kappa_{11}
        SUBJCBN_BETA1_TERMS_B(qd, 1 ) = REUSED_VALUE * FourThirds                               &
                                * ( 2.0_idp / R_SQUARE(rd)           * CUR_VAL_BETA(1,rd,td,pd)      &
                                + COTAN_VAL(td) /CUR_R_LOCS(rd)      * CUR_VAL_BETA(2,rd,td,pd)      &
                                - 2.0_idp/CUR_R_LOCS(rd)             * CUR_DRV_BETA(1,1,rd,td,pd)    &
                                + 1.0_idp/CUR_R_LOCS(rd)             * CUR_DRV_BETA(2,2,rd,td,pd)    &
                                + 1.0_idp/CUR_R_LOCS(rd)             * CUR_DRV_BETA(3,3,rd,td,pd)    )
        ! J_{1,3lmn}
        ! Reused_Value * kappa_{12}
        SUBJCBN_BETA1_TERMS_B(qd, 2 ) = REUSED_VALUE * FourThirds                              &
                                * ( - 2.0_idp /CUR_R_LOCS(rd)       * CUR_VAL_BETA(1,rd,td,pd)      &
                                - COTAN_VAL(td)                     * CUR_VAL_BETA(2,rd,td,pd)      &
                                + 2.0_idp                           * CUR_DRV_BETA(1,1,rd,td,pd)    &
                                - 1.0_idp                           * CUR_DRV_BETA(2,2,rd,td,pd)    &
                                - 1.0_idp                           * CUR_DRV_BETA(3,3,rd,td,pd)    )
        ! J_{1,3lmn}
        ! Reused_Value * kappa_{13}
        SUBJCBN_BETA1_TERMS_B(qd, 3 ) = REUSED_VALUE                                           &
                                * ( 2.0_idp / R_SQUARE(rd)      * CUR_DRV_BETA(2,1,rd,td,pd)        &
                                + 2.0_idp                       * CUR_DRV_BETA(1,2,rd,td,pd)    )
        ! J_{1,3lmn}
        ! Reused_Value * kappa_{14}
        SUBJCBN_BETA1_TERMS_B(qd, 4 ) = REUSED_VALUE                                           &
                                * ( 2.0_idp / RSIN_SQUARE(rd, td)   * CUR_DRV_BETA(3,1,rd,td,pd)    &
                                + 2.0_idp                           * CUR_DRV_BETA(1,3,rd,td,pd)    )


        ! J_{2,4lmn}
        SUBJCBN_BETA1_TERMS_B(qd, 5:8 ) = -7.0_idp * PSI_POWER(1)                          &
                                                    * SUBJCBN_BETA1_TERMS_B(qd, 1:4)



        ! J_{3,3lmn} Non-Derivative Term
        SUBJCBN_BETA1_TERMS_B(qd, 9 ) = (4.0_idp/(3.0_idp *CUR_R_LOCS(rd))) * JCBN_n_ARRAY(1)      &
                                            - ( 2.0_idp/(3.0_idp * R_SQUARE(rd) ) )

        ! J_{3,3lmn} Derivative Terms
        SUBJCBN_BETA1_TERMS_B(qd, 10 ) = OneThird * ((2.0_idp/CUR_R_LOCS(rd)) - 4.0_idp * JCBN_n_ARRAY(1) )
        SUBJCBN_BETA1_TERMS_B(qd, 11 ) = -(1.0_idp/R_SQUARE(rd)) * JCBN_n_ARRAY(2)
        SUBJCBN_BETA1_TERMS_B(qd, 12 ) = -(1.0_idp/RSIN_SQUARE(rd, td)) * JCBN_n_ARRAY(3)

        ! J_{3,3lmn} Double Derivative Term



        ! J_{4,3lmn} Non-Derivative Term
        SUBJCBN_BETA1_TERMS_B(qd, 13 ) = -TwoThirds*CUR_R_LOCS(rd)* JCBN_n_ARRAY(2)

        ! J_{4,3lmn} Derivative Terms
        SUBJCBN_BETA1_TERMS_B(qd, 14 ) = -TwoThirds * R_SQUARE(rd) * JCBN_n_ARRAY(2)
        SUBJCBN_BETA1_TERMS_B(qd, 15 ) = (2.0_idp / (3.0_idp * R_CUBE(rd))) - JCBN_n_ARRAY(1)*(1.0_idp/R_SQUARE(rd))
        SUBJCBN_BETA1_TERMS_B(qd, 16 ) = 0.0_idp


        ! J_{5,3lmn} Non-Derivative Term
        SUBJCBN_BETA1_TERMS_B(qd, 17 ) = -TwoThirds*CUR_R_LOCS(rd)*SIN_SQUARE(td)* JCBN_n_ARRAY(3)

        ! J_{5,3lmn} Derivative Terms
        SUBJCBN_BETA1_TERMS_B(qd, 18 ) = -FourThirds * RSIN_SQUARE(rd, td) * JCBN_n_ARRAY(3)
        SUBJCBN_BETA1_TERMS_B(qd, 19 ) = 0.0_idp
        SUBJCBN_BETA1_TERMS_B(qd, 20 ) = OneThird * RSIN_SQUARE(rd, td) * (( 2.0_idp /CUR_R_LOCS(rd) )             &
                                                + 3.0_idp * JCBN_n_ARRAY(1) )


        !
        !   Beta 2 Terms
        !




        ! J_{1,4lmn}
        ! Reused_Value * kappa_{20}
        SUBJCBN_BETA2_TERMS_B(qd, 1 ) = REUSED_VALUE * FourThirds                                      &
                                    * ( COTAN_VAL(td) /CUR_R_LOCS(rd)       * CUR_VAL_BETA(1,rd,td,pd)      &
                                            + 2.0_idp * COTAN_VAL(td)       * CUR_VAL_BETA(2,rd,td,pd)      &
                                            - COTAN_VAL(td)                 * CUR_DRV_BETA(2,2,rd,td,pd)    &
                                            + 2.0_idp * COTAN_VAL(td)       * CUR_DRV_BETA(3,3,rd,td,pd)    )
        ! J_{1,4lmn}
        ! Reused_Value * kappa_{21}
        SUBJCBN_BETA2_TERMS_B(qd, 2 ) = REUSED_VALUE                                                   &
                                            * ( 2.0_idp * R_SQUARE(rd)      * CUR_DRV_BETA(1,2,rd,td,pd)    &
                                                + 2.0_idp                   * CUR_DRV_BETA(2,1,rd,td,pd)    )
        ! J_{1,4lmn}
        ! Reused_Value * kappa_{22}
        SUBJCBN_BETA2_TERMS_B(qd, 3 ) = REUSED_VALUE * FourThirds                                      &
                                        * ( 1.0_idp /CUR_R_LOCS(rd)         * CUR_VAL_BETA(1,rd,td,pd)      &
                                            - COTAN_VAL(td)                 * CUR_VAL_BETA(2,rd,td,pd)      &
                                            - 1.0_idp                       * CUR_DRV_BETA(1,1,rd,td,pd)    &
                                            + 2.0_idp                       * CUR_DRV_BETA(2,2,rd,td,pd)    &
                                            - 1.0_idp                       * CUR_DRV_BETA(3,3,rd,td,pd)    )
            ! J_{1,4lmn}
        ! Reused_Value * kappa_{23}
        SUBJCBN_BETA2_TERMS_B(qd, 4 ) = REUSED_VALUE                                                   &
                                            * ( 2.0_idp / SIN_SQUARE(td)    * CUR_DRV_BETA(3,2,rd,td,pd)    &
                                                + 2.0_idp                   * CUR_DRV_BETA(2,3,rd,td,pd)    )


        ! J_{2,4lmn}
        SUBJCBN_BETA2_TERMS_B(qd, 5:8 ) = -7.0_idp * PSI_POWER(1)                      &
                                                * SUBJCBN_BETA2_TERMS_B(qd, 1:4)



        ! J_{3,4lmn} Non-Derivative Term
        SUBJCBN_BETA2_TERMS_B(qd, 9 ) = TwoThirds*COTAN_VAL(td) * JCBN_n_ARRAY(1)

        ! J_{3,4lmn} Derivative Terms
        SUBJCBN_BETA2_TERMS_B(qd, 10 ) = OneThird * COTAN_VAL(td) - JCBN_n_ARRAY(2)
        SUBJCBN_BETA2_TERMS_B(qd, 11 ) = TwoThirds * JCBN_n_ARRAY(1)
        SUBJCBN_BETA2_TERMS_B(qd, 12 ) = 0.0_idp


        ! J_{4,4lmn} Non-Derivative Term
        SUBJCBN_BETA2_TERMS_B(qd, 13 ) = TwoThirds * R_SQUARE(rd) * COTAN_VAL(td) * JCBN_n_ARRAY(2)               &
                                                - 1.0_idp /( 3.0_idp * RSIN_SQUARE(rd, td) )

        ! J_{4,4lmn} Derivative Terms
        SUBJCBN_BETA2_TERMS_B(qd, 14 ) = -JCBN_n_ARRAY(1)
        SUBJCBN_BETA2_TERMS_B(qd, 15 ) = OneThird * R_SQUARE(rd) * ( COTAN_VAL(td) - 4.0_idp * JCBN_n_ARRAY(2) )
        SUBJCBN_BETA2_TERMS_B(qd, 16 ) = - (1.0_idp/RSIN_SQUARE(rd, td)) * JCBN_n_ARRAY(3)


        ! J_{5,4lmn} Non-Derivative Term
        SUBJCBN_BETA2_TERMS_B(qd, 17 ) = -FourThirds*R_SQUARE(rd)*SIN_VAL(td)*COS_VAL(td)* JCBN_n_ARRAY(3)

        ! J_{5,4lmn} Derivative Terms
        SUBJCBN_BETA2_TERMS_B(qd, 18 ) = 0.0_idp
        SUBJCBN_BETA2_TERMS_B(qd, 19 ) = -FourThirds * RSIN_SQUARE(rd, td) * JCBN_n_ARRAY(3)
        SUBJCBN_BETA2_TERMS_B(qd, 20 ) = OneThird * RSIN_SQUARE(rd, td) * ( COTAN_VAL(td) + 3.0_idp * JCBN_n_ARRAY(2)  )






        !
        !   Beta 3 Terms
        !




        ! Reused_Value * kappa_{30}
        SUBJCBN_BETA3_TERMS_B(qd, 1 ) = 0.0_idp

        ! Reused_Value * kappa_{31}
        SUBJCBN_BETA3_TERMS_B(qd, 2 ) = REUSED_VALUE                                           &
                                * ( 2.0_idp * RSIN_SQUARE(rd, td)   * CUR_DRV_BETA(1,3,rd,td,pd)    &
                                + 2.0_idp                   * CUR_DRV_BETA(3,1,rd,td,pd)    )

        ! Reused_Value * kappa_{32}
        SUBJCBN_BETA3_TERMS_B(qd, 3 ) = REUSED_VALUE                                           &
                                * ( 2.0_idp * SIN_SQUARE(td)    * CUR_DRV_BETA(2,3,rd,td,pd)        &
                                + 2.0_idp                   * CUR_DRV_BETA(3,2,rd,td,pd)    )

        ! Reused_Value * kappa_{33}
        SUBJCBN_BETA3_TERMS_B(qd, 4 ) = REUSED_VALUE * FourThirds                              &
                                * ( 1.0_idp /CUR_R_LOCS(rd)         * CUR_VAL_BETA(1,rd,td,pd)      &
                                + 2.0_idp * COTAN_VAL(td)       * CUR_VAL_BETA(2,rd,td,pd)          &
                                - 1.0_idp                   * CUR_DRV_BETA(1,1,rd,td,pd)            &
                                - 1.0_idp                   * CUR_DRV_BETA(2,2,rd,td,pd)            &
                                + 2.0_idp                   * CUR_DRV_BETA(3,3,rd,td,pd)    )


        ! J_{2,5lmn}
        SUBJCBN_BETA3_TERMS_B(qd, 5:8 ) = -7.0_idp * PSI_POWER(1)                              &
                                                * SUBJCBN_BETA3_TERMS_B(qd, 1:4)










        ! J_{3,5lmn} Non-Derivative Term
        SUBJCBN_BETA3_TERMS_B(qd, 9 ) =  (2.0_idp/3.0_idp)*COTAN_VAL(td)* JCBN_n_ARRAY(1)

        ! J_{3,5lmn} Derivative Terms
        SUBJCBN_BETA3_TERMS_B(qd, 10 ) = -JCBN_n_ARRAY(3) + COTAN_VAL(td)/3.0_idp
        SUBJCBN_BETA3_TERMS_B(qd, 11 ) = 0.0_idp
        SUBJCBN_BETA3_TERMS_B(qd, 12 ) = -TwoThirds * JCBN_n_ARRAY(1)


        ! J_{4,5lmn} Non-Derivative Term
        SUBJCBN_BETA3_TERMS_B(qd, 13 ) = ( RSIN_SQUARE(rd,td)*SIN_VAL(td)*COS_VAL(td) - R_SQUARE(rd)*COTAN_VAL(td) )    &
                                    * JCBN_n_ARRAY(3)

        ! J_{4,5lmn} Derivative Terms
        SUBJCBN_BETA3_TERMS_B(qd, 14 ) = 0.0_idp
        SUBJCBN_BETA3_TERMS_B(qd, 15 ) = FourThirds * R_SQUARE(rd) * JCBN_n_ARRAY(3)
        SUBJCBN_BETA3_TERMS_B(qd, 16 ) = -FourThirds * R_SQUARE(rd) * JCBN_n_ARRAY(2)


        ! J_{5,5lmn} Non-Derivative Term
        SUBJCBN_BETA3_TERMS_B(qd, 17 ) = (R_CUBE(rd)*SIN_SQUARE(td)*SIN_SQUARE(td) - R_INVERSE(rd) )  &
                                    * JCBN_n_ARRAY(1)                                           &
                              + ( RSIN_SQUARE(rd, td)*SIN_VAL(td)*COS_VAL(td) - R_SQUARE(rd)*COTAN_VAL(td) )    &
                                    * JCBN_n_ARRAY(2)
        ! J_{5,5lmn} Derivative Terms
        SUBJCBN_BETA3_TERMS_B(qd, 18 ) = JCBN_n_ARRAY(1)
        SUBJCBN_BETA3_TERMS_B(qd, 19 ) = R_SQUARE(rd) * JCBN_n_ARRAY(2)
        SUBJCBN_BETA3_TERMS_B(qd, 20 ) = FourThirds * R_SQUARE(rd) * JCBN_n_ARRAY(3)







        END DO ! pd loop
    END DO  ! td loop
END DO  ! rd loop


!$OMP END DO

!$OMP END PARALLEL

END SUBROUTINE Calc_3D_SubJcbn_Terms_B







!+203+###########################################################################!
!                                                                                !
!                  Calc_3D_SubJcbn_Terms                                         !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_3D_SubJcbn_Terms_C(   re, te, pe,                                     &
                                    CUR_VAL_PSI, CUR_DRV_PSI,                       &
                                    CUR_VAL_ALPHAPSI, CUR_DRV_ALPHAPSI,             &
                                    CUR_VAL_BETA, CUR_DRV_BETA, CUR_DDRV_BETA,      &
                                    RHS_TERMS,                                      &
                                    SUBJCBN_PSI_TERMS_B, SUBJCBN_ALPHAPSI_TERMS_B,      &
                                    SUBJCBN_BETA1_TERMS_B, SUBJCBN_BETA2_TERMS_B,       &
                                    SUBJCBN_BETA3_TERMS_B,                            &
                                    CUR_R_LOCS, CUR_T_LOCS, CUR_P_LOCS,             &
                                    R_SQUARE, R_CUBE, R_INVERSE,                    &
                                    SIN_VAL, COS_VAL, COTAN_VAL,                    &
                                    SIN_SQUARE, RSIN_SQUARE                         )




INTEGER, INTENT(IN)                                                     ::  re, te, pe

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:NUM_R_QUAD_POINTS,   &
                                         1:NUM_T_QUAD_POINTS,   &
                                         1:NUM_P_QUAD_POINTS    )       ::  CUR_VAL_PSI

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:NUM_R_QUAD_POINTS,   &
                                         1:NUM_T_QUAD_POINTS,   &
                                         1:NUM_P_QUAD_POINTS    )       ::  CUR_VAL_ALPHAPSI


REAL(KIND = idp), INTENT(IN), DIMENSION( 1:3,                   &
                                         1:NUM_R_QUAD_POINTS,   &
                                         1:NUM_T_QUAD_POINTS,   &
                                         1:NUM_P_QUAD_POINTS    )       ::  CUR_VAL_BETA


REAL(KIND = idp), INTENT(IN), DIMENSION( 1:3,                   &
                                         1:NUM_R_QUAD_POINTS,   &
                                         1:NUM_T_QUAD_POINTS,   &
                                         1:NUM_P_QUAD_POINTS    )       ::  CUR_DRV_PSI

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:3,                   &
                                         1:NUM_R_QUAD_POINTS,   &
                                         1:NUM_T_QUAD_POINTS,   &
                                         1:NUM_P_QUAD_POINTS    )       ::  CUR_DRV_ALPHAPSI

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:3, 1:3,              &
                                         1:NUM_R_QUAD_POINTS,   &
                                         1:NUM_T_QUAD_POINTS,   &
                                         1:NUM_P_QUAD_POINTS    )       ::  CUR_DRV_BETA

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:3, 1:3, 1:3,         &
                                         1:NUM_R_QUAD_POINTS,   &
                                         1:NUM_T_QUAD_POINTS,   &
                                         1:NUM_P_QUAD_POINTS    )       ::  CUR_DDRV_BETA





REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_R_QUAD_POINTS)            ::  CUR_R_LOCS
REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_T_QUAD_POINTS)            ::  CUR_T_LOCS
REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_P_QUAD_POINTS)            ::  CUR_P_LOCS


REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_R_QUAD_POINTS,    &
                                        1:NUM_T_QUAD_POINTS     )       ::  RSIN_SQUARE

REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_R_QUAD_POINTS)            ::  R_SQUARE,           &
                                                                            R_CUBE,             &
                                                                            R_INVERSE


REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_T_QUAD_POINTS)            ::  SIN_VAL,            &
                                                                            SIN_SQUARE,         &
                                                                            COTAN_VAL,          &
                                                                            COS_VAL



REAL(KIND = idp), INTENT(INOUT), DIMENSION( 1:5,                        &
                                            1:NUM_R_QUAD_POINTS,        &
                                            1:NUM_T_QUAD_POINTS,        &
                                            1:NUM_P_QUAD_POINTS         )   ::  RHS_TERMS


REAL(KIND = idp), INTENT(INOUT), DIMENSION( 1:14, 1:NUM_QUAD_DOF )          ::  SUBJCBN_PSI_TERMS_B

REAL(KIND = idp), INTENT(INOUT), DIMENSION( 1:14, 1:NUM_QUAD_DOF )          ::  SUBJCBN_ALPHAPSI_TERMS_B

REAL(KIND = idp), INTENT(INOUT), DIMENSION( 1:20, 1:NUM_QUAD_DOF )          ::  SUBJCBN_BETA1_TERMS_B

REAL(KIND = idp), INTENT(INOUT), DIMENSION( 1:20, 1:NUM_QUAD_DOF )          ::  SUBJCBN_BETA2_TERMS_B

REAL(KIND = idp), INTENT(INOUT), DIMENSION( 1:20, 1:NUM_QUAD_DOF )          ::  SUBJCBN_BETA3_TERMS_B


REAL(KIND = idp)                                                            ::  REUSED_VALUE

INTEGER                                                                     ::  pd, td, rd,     &
                                                                                i,qd


REAL(KIND = idp), DIMENSION(1:8)                                            ::  PSI_POWER
REAL(KIND = idp), DIMENSION(1:4)                                            ::  ALPHAPSI_POWER


REAL(KIND = idp), DIMENSION(1:3,1:3)                                        ::  JCBN_mu_Array
REAL(KIND = idp), DIMENSION(1:3)                                            ::  JCBN_n_ARRAY

REAL(KIND = idp)                                                            ::  JCBN_BIGK_VALUE












!$OMP PARALLEL DEFAULT(none)                                            &
!$OMP PRIVATE( pd, td, rd, qd,                                          &
!$OMP           PSI_POWER, ALPHAPSI_POWER,                              &
!$OMP           JCBN_BIGK_VALUE, JCBN_mu_ARRAY,JCBN_n_ARRAY,            &
!$OMP           REUSED_VALUE                                        )   &
!$OMP SHARED( re, te, pe,                                               &
!$OMP           CUR_VAL_PSI, CUR_VAL_ALPHAPSI, CUR_VAL_BETA,            &
!$OMP           CUR_DRV_PSI, CUR_DRV_ALPHAPSI, CUR_DRV_BETA,            &
!$OMP           CUR_DDRV_BETA,                                          &
!$OMP           CUR_R_LOCS, CUR_T_LOCS, CUR_P_LOCS,                     &
!$OMP           RSIN_SQUARE, R_SQUARE, R_CUBE, R_INVERSE,               &
!$OMP           SIN_VAL, SIN_SQUARE, COTAN_VAL, COS_VAL,                &
!$OMP           NUM_OFF_DIAGONALS,                                      &
!$OMP           RHS_TERMS,                                              &
!$OMP           SUBJCBN_PSI_TERMS_B, SUBJCBN_ALPHAPSI_TERMS_B,              &
!$OMP           SUBJCBN_BETA1_TERMS_B, SUBJCBN_BETA2_TERMS_B,               &
!$OMP           SUBJCBN_BETA3_TERMS_B,                                    &
!$OMP           Block_SOURCE_E, Block_SOURCE_S, Block_SOURCE_Si,        &
!$OMP           OneThird, OneThirtySecond, FourThirds, TwoThirds,       &
!$OMP           LM_Length                           )



!$OMP DO SCHEDULE(dynamic), COLLAPSE(3)
DO pd = 1,NUM_P_QUAD_POINTS
 DO td = 1,NUM_T_QUAD_POINTS
  DO rd = 1,NUM_R_QUAD_POINTS

        qd = ((pd-1) * NUM_T_QUAD_POINTS + td-1 ) * NUM_R_QUAD_POINTS + rd


        PSI_POWER(1) = CUR_VAL_PSI(rd, td, pd)
        DO i = 2,8

            PSI_POWER(i) = PSI_POWER(i-1)*PSI_POWER(1)

        END DO


        ALPHAPSI_POWER(1) = CUR_VAL_ALPHAPSI(rd, td, pd)
        DO i = 2,4

            ALPHAPSI_POWER(i) = ALPHAPSI_POWER(i-1)*ALPHAPSI_POWER(1)

        END DO




        ! K_{ij}K^{ij} = Psi^{14}/AlphaPsi^{2} * BIGK
        JCBN_BIGK_VALUE = JCBN_BIGK_FUNCTION(rd, td, pd,                                        &
                                                CUR_VAL_BETA, CUR_DRV_BETA,                     &
                                                CUR_R_LOCS(rd), R_SQUARE(rd), SIN_SQUARE(td),   &
                                                RSIN_SQUARE(rd, td), COTAN_VAL(td)              )




        JCBN_mu_Array = JCBN_mu_FUNCTION_3D_ALL( rd, td, pd,                                &
                                                CUR_R_LOCS(rd), R_SQUARE(rd), R_CUBE(rd),   &
                                                R_INVERSE(rd), RSIN_SQUARE(rd, td),         &
                                                SIN_VAL(td), SIN_SQUARE(td),                &
                                                COS_VAL(td), COTAN_VAL(td),                 &
                                                CUR_VAL_BETA, CUR_DRV_BETA                          )





        JCBN_n_ARRAY(:) = CUR_DRV_ALPHAPSI(:, rd, td, pd) / ALPHAPSI_POWER(1)   &
                            - 7 * CUR_DRV_PSI(:, rd, td, pd )/ PSI_POWER(1)










        RHS_Terms(1, rd, td, pd) = - TwoPi * Block_Source_E(rd, td, pd, re, te, pe) * PSI_POWER(5)            &
                                 - PSI_POWER(7)/ (16.0_idp * ALPHAPSI_POWER(2)) * JCBN_BIGK_VALUE



        RHS_Terms(2, rd, td, pd) = TwoPi * ALPHAPSI_POWER(1) * PSI_POWER(4)                                     &
                        * ( Block_Source_E(rd, td, pd, re, te, pe) + 2.0_idp * Block_Source_S(rd, td, pd, re, te, pe)  )    &
                        + 7.0_idp*PSI_POWER(6)/ (16.0_idp * ALPHAPSI_POWER(1)) * JCBN_BIGK_VALUE






        RHS_Terms(3, rd, td, pd) = 16.0_idp * pi                                                &
                          * ALPHAPSI_POWER(1)                                                   &
                          * PSI_POWER(3)                                                        &
                          * Block_Source_Si(rd, td, pd, re, te, pe, 1)                          &
                     - OneThird                                                                 &
                          * ( CUR_DDRV_BETA(1, 1, 1, rd, td, pd )                               &
                            + CUR_DDRV_BETA(1, 2, 2, rd, td, pd )                               &
                            + CUR_DDRV_BETA(1, 3, 3, rd, td, pd )                               &
                            + 2.0_idp/CUR_R_LOCS(rd) * CUR_DRV_BETA(1, 1, rd, td, pd)           &
                            + COTAN_VAL(td) * CUR_DRV_BETA(1, 2, rd, td, pd)                    &
                            + COTAN_VAL(td) * CUR_DRV_BETA(1, 3, rd, td, pd)                    &
                            - 2.0_idp /R_SQUARE(rd) * CUR_VAL_BETA(1, rd, td, pd )              &
                            )                                                                   &
                      + JCBN_n_ARRAY(1)* JCBN_mu_Array(1,1)                                     &
                      + JCBN_n_ARRAY(2)* JCBN_mu_Array(2,1)                                     &
                      + JCBN_n_ARRAY(3)* JCBN_mu_Array(3,1)





        RHS_Terms(4, rd, td, pd) =  16.0_idp * pi                                               &
                           * ALPHAPSI_POWER(1)                                                  &
                           * PSI_POWER(3)                                                       &
                           * Block_Source_Si(rd, td, pd, re, te, pe, 2)                         &
                        - OneThird                                                              &
                            * ( R_SQUARE(rd)                                                    &
                             * ( CUR_DDRV_BETA(2, 1, 1, rd, td, pd )                            &
                               + CUR_DDRV_BETA(2, 2, 2, rd, td, pd )                            &
                               + CUR_DDRV_BETA(2, 3, 3, rd, td, pd )                            &
                               + 2.0_idp/CUR_R_LOCS(rd) * CUR_DRV_BETA(2, 1, rd, td, pd)        &
                               + COTAN_VAL(td)                                                  &
                                * ( CUR_DRV_BETA(2, 2, rd, td, pd)                              &
                                    + CUR_DRV_BETA(2, 3, rd, td, pd)  )                         &
                               - ( 1.0_idp / SIN_SQUARE(td) )                                   &
                                * ( CUR_VAL_BETA(2, rd, td, pd)                                 &
                                    - CUR_VAL_BETA(3, rd, td, pd)     )                         &
                               )                                                                &
                             )                                                                  &
                        + JCBN_n_ARRAY(1) * JCBN_mu_Array(1,2)                                  &
                        + JCBN_n_ARRAY(2) * JCBN_mu_Array(2,2)                                  &
                        + JCBN_n_ARRAY(3) * JCBN_mu_Array(3,2)





        RHS_Terms(5, rd, td, pd) = 16.0_idp * pi                                            &
                          * ALPHAPSI_POWER(1)                                               &
                          * PSI_POWER(3)                                                    &
                          * Block_Source_Si(rd, td, pd, re, te, pe, 3)                      &
                        - OneThird                                                          &
                            * ( RSIN_SQUARE(rd, td)                                         &
                             * ( CUR_DDRV_BETA(3, 1, 1, rd, td, pd )                        &
                               + CUR_DDRV_BETA(3, 2, 2, rd, td, pd )                        &
                               + CUR_DDRV_BETA(3, 3, 3, rd, td, pd )                        &
                               + 2.0_idp/CUR_R_LOCS(rd) * CUR_DRV_BETA(3, 1, rd, td, pd)    &
                               + COTAN_VAL(td) * CUR_DRV_BETA(3, 2, rd, td, pd)             &
                               + COTAN_VAL(td) * CUR_DRV_BETA(3, 3, rd, td, pd)             &
                               )                                                            &
                             )                                                              &
                        + JCBN_n_ARRAY(1) * JCBN_mu_Array(1,3)                              &
                        + JCBN_n_ARRAY(2) * JCBN_mu_Array(2,3)                              &
                        + JCBN_n_ARRAY(3) * JCBN_mu_Array(3,3)







        !
        !   With Respect to Psi Terms Jacobian Terms
        !

        ! J_{1,1lmn}
        SUBJCBN_PSI_TERMS_B(1, qd ) = 10.0_idp*pi                                                  &
                                                * Block_Source_E(rd, td, pd, re, te, pe)* PSI_POWER(4)  &
                                         + (7.0_idp/16.0_idp)                                           &
                                                * PSI_POWER(6)/ALPHAPSI_POWER(2)                        &
                                                * JCBN_BIGK_VALUE

        ! J_{2,1lmn}
        SUBJCBN_PSI_TERMS_B(2, qd ) = -8.0_idp * pi * ALPHAPSI_POWER(1) * PSI_POWER(3)             &
                                            * ( Block_Source_E(rd, td, pd, re, te, pe)                  &
                                                + 2 * Block_Source_S(rd, td, pd, re, te, pe) )          &
                                            - (42.0_idp / 16.0_idp )                                    &
                                                * PSI_POWER(5)/ALPHAPSI_POWER(1)                        &
                                                * JCBN_BIGK_VALUE

        ! J_{3,1lmn} Non-Derivative Term
        SUBJCBN_PSI_TERMS_B(3, qd ) = ( 7.0_idp / PSI_POWER(2) )                                   &
                                                * (CUR_DRV_PSI(1, rd, td, pd )*JCBN_mu_ARRAY(1,1)       &
                                                  + CUR_DRV_PSI(2, rd, td, pd )*JCBN_mu_ARRAY(2,1)      &
                                                  + CUR_DRV_PSI(3, rd, td, pd )*JCBN_mu_ARRAY(3,1) )    &
                                            - 48.0_idp * pi                                             &
                                                * Block_Source_Si(rd, td, pd, re, te, pe, 1)            &
                                                * ALPHAPSI_POWER(1)                                     &
                                                * PSI_POWER(2)


        ! J_{3,1lmn} Derivative Terms
        SUBJCBN_PSI_TERMS_B(4:6, qd ) = (7.0_idp/ PSI_POWER(1) )*JCBN_mu_ARRAY(1:3,1)



        ! J_{4,1lmn} Non-Derivative Term
        SUBJCBN_PSI_TERMS_B(7, qd ) = (-7.0_idp / PSI_POWER(2) )                                   &
                                                * (CUR_DRV_PSI(1, rd, td, pd)*JCBN_mu_ARRAY(1,2)        &
                                                 + CUR_DRV_PSI(2, rd, td, pd)*JCBN_mu_ARRAY(2,2)        &
                                                 + CUR_DRV_PSI(3, rd, td, pd)*JCBN_mu_ARRAY(3,2) )      &
                                            - 48.0_idp * pi                                             &
                                                * Block_Source_Si(rd, td, pd, re, te, pe, 2)            &
                                                * ALPHAPSI_POWER(1)                                     &
                                                * PSI_POWER(2)

        ! J_{4,1lmn} Derivative Terms
        SUBJCBN_PSI_TERMS_B(8:10, qd ) = (7.0_idp/ PSI_POWER(1) )*JCBN_mu_ARRAY(1:3,2)






        ! J_{5,1lmn} Non-Derivative Term
        SUBJCBN_PSI_TERMS_B(11, qd ) = (-7.0_idp / PSI_POWER(2) )                  &
                                * (CUR_DRV_PSI(1, rd, td, pd)*JCBN_mu_ARRAY(1,3)        &
                                 + CUR_DRV_PSI(2, rd, td, pd)*JCBN_mu_ARRAY(2,3)        &
                                 + CUR_DRV_PSI(3, rd, td, pd)*JCBN_mu_ARRAY(3,3) )      &
                            - 48.0_idp * pi                                             &
                                * Block_Source_Si(rd, td, pd, re, te, pe, 3)            &
                                * ALPHAPSI_POWER(1)                                     &
                                * PSI_POWER(2)

        ! J_{5,1lmn} Derivative Terms
        SUBJCBN_PSI_TERMS_B( 12:14, qd ) = (7.0_idp/ PSI_POWER(1) )*JCBN_mu_ARRAY(1:3,3)






        !
        ! Alpha Psi Terms
        !


        ! J_{1,2lmn}
        SUBJCBN_ALPHAPSI_TERMS_B(1, qd ) = -PSI_POWER(7)/(8.0_idp*ALPHAPSI_POWER(3))               &
                                                * JCBN_BIGK_VALUE

        ! J_{2,2lmn}
        SUBJCBN_ALPHAPSI_TERMS_B(2, qd ) =  -TwoPi * PSI_POWER(4)                                  &
                                                *( Block_Source_E(rd, td, pd, re, te, pe)               &
                                                    + 2.0_idp * Block_Source_S(rd, td, pd, re, te, pe) )      &
                                                - (7.0_idp * OneThirtySecond)                           &
                                                    * PSI_POWER(3)/ALPHAPSI_POWER(3)                    &
                                                    * JCBN_BIGK_VALUE

        ! J_{3,2lmn} Non-Derivative Term
        SUBJCBN_ALPHAPSI_TERMS_B(3, qd ) = -(CUR_DRV_ALPHAPSI(1, rd, td, pd )*JCBN_mu_Array(1,1)       &
                                                + CUR_DRV_ALPHAPSI(2, rd, td, pd )*JCBN_mu_Array(2,1)       &
                                                + CUR_DRV_ALPHAPSI(3, rd, td, pd )*JCBN_mu_Array(3,1)   )   &
                                                /ALPHAPSI_POWER(1)                                          &
                                                -16.0_idp * pi * Block_Source_Si(rd,td,pd,re,te,pe,1)       &
                                                    * PSI_POWER(3)


        ! J_{3,2lmn} Derivative Terms
        SUBJCBN_ALPHAPSI_TERMS_B(4:6, qd ) = (-1.0_idp/ALPHAPSI_POWER(1) )*JCBN_mu_ARRAY(1:3,1)




        ! J_{4,2lmn} Non-Derivative Term
        SUBJCBN_ALPHAPSI_TERMS_B(7, qd ) = -( CUR_DRV_ALPHAPSI(1, rd, td, pd )*JCBN_mu_Array(1,2)      &
                                                + CUR_DRV_ALPHAPSI(2, rd, td, pd )*JCBN_mu_Array(2,2)       &
                                                + CUR_DRV_ALPHAPSI(3, rd, td, pd )*JCBN_mu_Array(3,2)   )   &
                                                /ALPHAPSI_POWER(1)                                          &
                                                -16.0_idp * pi * Block_Source_Si(rd, td, pd, re, te, pe, 2) &
                                                    * PSI_POWER(3)


        ! J_{4,2lmn} Derivative Terms
        SUBJCBN_ALPHAPSI_TERMS_B(8:10, qd ) = (-1.0_idp/ ALPHAPSI_POWER(1) )*JCBN_mu_ARRAY(1:3,2)



        ! J_{5,2lmn} Non-Derivative Term
        SUBJCBN_ALPHAPSI_TERMS_B(11, qd ) = -(CUR_DRV_ALPHAPSI(1, rd, td, pd ) * JCBN_mu_Array(1,3)    &
                                                 + CUR_DRV_ALPHAPSI(2, rd, td, pd ) * JCBN_mu_Array(2,3)    &
                                                 + CUR_DRV_ALPHAPSI(3, rd, td, pd ) * JCBN_mu_Array(3,3)  ) &
                                                /ALPHAPSI_POWER(1)                                          &
                                                -16.0_idp * pi * Block_Source_Si(rd,td,pd,re,te,pe,3)       &
                                                    * PSI_POWER(3)

        ! J_{5,2lmn} Derivative Terms
        SUBJCBN_ALPHAPSI_TERMS_B(12:14, qd ) = (-1.0_idp/ ALPHAPSI_POWER(1) )*JCBN_mu_ARRAY(1:3,3)






        !
        !   Beta 1 Terms
        !




        REUSED_VALUE = (1.0_idp/16_idp) * ( PSI_POWER(7) / ALPHAPSI_POWER(2) )


        ! J_{1,3lmn}
        ! Reused_Value * kappa_{11}
        SUBJCBN_BETA1_TERMS_B(1, qd ) = REUSED_VALUE * FourThirds                               &
                                * ( 2.0_idp / R_SQUARE(rd)           * CUR_VAL_BETA(1,rd,td,pd)      &
                                + COTAN_VAL(td) /CUR_R_LOCS(rd)      * CUR_VAL_BETA(2,rd,td,pd)      &
                                - 2.0_idp/CUR_R_LOCS(rd)             * CUR_DRV_BETA(1,1,rd,td,pd)    &
                                + 1.0_idp/CUR_R_LOCS(rd)             * CUR_DRV_BETA(2,2,rd,td,pd)    &
                                + 1.0_idp/CUR_R_LOCS(rd)             * CUR_DRV_BETA(3,3,rd,td,pd)    )
        ! J_{1,3lmn}
        ! Reused_Value * kappa_{12}
        SUBJCBN_BETA1_TERMS_B(2, qd ) = REUSED_VALUE * FourThirds                              &
                                * ( - 2.0_idp /CUR_R_LOCS(rd)       * CUR_VAL_BETA(1,rd,td,pd)      &
                                - COTAN_VAL(td)                     * CUR_VAL_BETA(2,rd,td,pd)      &
                                + 2.0_idp                           * CUR_DRV_BETA(1,1,rd,td,pd)    &
                                - 1.0_idp                           * CUR_DRV_BETA(2,2,rd,td,pd)    &
                                - 1.0_idp                           * CUR_DRV_BETA(3,3,rd,td,pd)    )
        ! J_{1,3lmn}
        ! Reused_Value * kappa_{13}
        SUBJCBN_BETA1_TERMS_B(3, qd ) = REUSED_VALUE                                           &
                                * ( 2.0_idp / R_SQUARE(rd)      * CUR_DRV_BETA(2,1,rd,td,pd)        &
                                + 2.0_idp                       * CUR_DRV_BETA(1,2,rd,td,pd)    )
        ! J_{1,3lmn}
        ! Reused_Value * kappa_{14}
        SUBJCBN_BETA1_TERMS_B(4, qd ) = REUSED_VALUE                                           &
                                * ( 2.0_idp / RSIN_SQUARE(rd, td)   * CUR_DRV_BETA(3,1,rd,td,pd)    &
                                + 2.0_idp                           * CUR_DRV_BETA(1,3,rd,td,pd)    )


        ! J_{2,4lmn}
        SUBJCBN_BETA1_TERMS_B(5:8, qd ) = -7.0_idp * PSI_POWER(1)                          &
                                                    * SUBJCBN_BETA1_TERMS_B(1:4, qd)



        ! J_{3,3lmn} Non-Derivative Term
        SUBJCBN_BETA1_TERMS_B(9, qd ) = (4.0_idp/(3.0_idp *CUR_R_LOCS(rd))) * JCBN_n_ARRAY(1)      &
                                            - ( 2.0_idp/(3.0_idp * R_SQUARE(rd) ) )

        ! J_{3,3lmn} Derivative Terms
        SUBJCBN_BETA1_TERMS_B(10, qd ) = OneThird * ((2.0_idp/CUR_R_LOCS(rd)) - 4.0_idp * JCBN_n_ARRAY(1) )
        SUBJCBN_BETA1_TERMS_B(11, qd ) = -(1.0_idp/R_SQUARE(rd)) * JCBN_n_ARRAY(2)
        SUBJCBN_BETA1_TERMS_B(12, qd ) = -(1.0_idp/RSIN_SQUARE(rd, td)) * JCBN_n_ARRAY(3)

        ! J_{3,3lmn} Double Derivative Term



        ! J_{4,3lmn} Non-Derivative Term
        SUBJCBN_BETA1_TERMS_B(13, qd ) = -TwoThirds*CUR_R_LOCS(rd)* JCBN_n_ARRAY(2)

        ! J_{4,3lmn} Derivative Terms
        SUBJCBN_BETA1_TERMS_B(14, qd ) = -TwoThirds * R_SQUARE(rd) * JCBN_n_ARRAY(2)
        SUBJCBN_BETA1_TERMS_B(15, qd ) = (2.0_idp / (3.0_idp * R_CUBE(rd))) - JCBN_n_ARRAY(1)*(1.0_idp/R_SQUARE(rd))
        SUBJCBN_BETA1_TERMS_B(16, qd ) = 0.0_idp


        ! J_{5,3lmn} Non-Derivative Term
        SUBJCBN_BETA1_TERMS_B(17, qd ) = -TwoThirds*CUR_R_LOCS(rd)*SIN_SQUARE(td)* JCBN_n_ARRAY(3)

        ! J_{5,3lmn} Derivative Terms
        SUBJCBN_BETA1_TERMS_B(18, qd ) = -FourThirds * RSIN_SQUARE(rd, td) * JCBN_n_ARRAY(3)
        SUBJCBN_BETA1_TERMS_B(19, qd ) = 0.0_idp
        SUBJCBN_BETA1_TERMS_B(20, qd ) = OneThird * RSIN_SQUARE(rd, td) * (( 2.0_idp /CUR_R_LOCS(rd) )             &
                                                + 3.0_idp * JCBN_n_ARRAY(1) )


        !
        !   Beta 2 Terms
        !




        ! J_{1,4lmn}
        ! Reused_Value * kappa_{20}
        SUBJCBN_BETA2_TERMS_B(1, qd ) = REUSED_VALUE * FourThirds                                      &
                                    * ( COTAN_VAL(td) /CUR_R_LOCS(rd)       * CUR_VAL_BETA(1,rd,td,pd)      &
                                            + 2.0_idp * COTAN_VAL(td)       * CUR_VAL_BETA(2,rd,td,pd)      &
                                            - COTAN_VAL(td)                 * CUR_DRV_BETA(2,2,rd,td,pd)    &
                                            + 2.0_idp * COTAN_VAL(td)       * CUR_DRV_BETA(3,3,rd,td,pd)    )
        ! J_{1,4lmn}
        ! Reused_Value * kappa_{21}
        SUBJCBN_BETA2_TERMS_B(2, qd ) = REUSED_VALUE                                                   &
                                            * ( 2.0_idp * R_SQUARE(rd)      * CUR_DRV_BETA(1,2,rd,td,pd)    &
                                                + 2.0_idp                   * CUR_DRV_BETA(2,1,rd,td,pd)    )
        ! J_{1,4lmn}
        ! Reused_Value * kappa_{22}
        SUBJCBN_BETA2_TERMS_B(3, qd ) = REUSED_VALUE * FourThirds                                      &
                                        * ( 1.0_idp /CUR_R_LOCS(rd)         * CUR_VAL_BETA(1,rd,td,pd)      &
                                            - COTAN_VAL(td)                 * CUR_VAL_BETA(2,rd,td,pd)      &
                                            - 1.0_idp                       * CUR_DRV_BETA(1,1,rd,td,pd)    &
                                            + 2.0_idp                       * CUR_DRV_BETA(2,2,rd,td,pd)    &
                                            - 1.0_idp                       * CUR_DRV_BETA(3,3,rd,td,pd)    )
            ! J_{1,4lmn}
        ! Reused_Value * kappa_{23}
        SUBJCBN_BETA2_TERMS_B(4, qd ) = REUSED_VALUE                                                   &
                                            * ( 2.0_idp / SIN_SQUARE(td)    * CUR_DRV_BETA(3,2,rd,td,pd)    &
                                                + 2.0_idp                   * CUR_DRV_BETA(2,3,rd,td,pd)    )


        ! J_{2,4lmn}
        SUBJCBN_BETA2_TERMS_B(5:8, qd ) = -7.0_idp * PSI_POWER(1)                      &
                                                * SUBJCBN_BETA2_TERMS_B(1:4, qd)



        ! J_{3,4lmn} Non-Derivative Term
        SUBJCBN_BETA2_TERMS_B(9, qd ) = TwoThirds*COTAN_VAL(td) * JCBN_n_ARRAY(1)

        ! J_{3,4lmn} Derivative Terms
        SUBJCBN_BETA2_TERMS_B(10, qd ) = OneThird * COTAN_VAL(td) - JCBN_n_ARRAY(2)
        SUBJCBN_BETA2_TERMS_B(11, qd ) = TwoThirds * JCBN_n_ARRAY(1)
        SUBJCBN_BETA2_TERMS_B(12, qd ) = 0.0_idp


        ! J_{4,4lmn} Non-Derivative Term
        SUBJCBN_BETA2_TERMS_B(13, qd ) = TwoThirds * R_SQUARE(rd) * COTAN_VAL(td) * JCBN_n_ARRAY(2)               &
                                                - 1.0_idp /( 3.0_idp * RSIN_SQUARE(rd, td) )

        ! J_{4,4lmn} Derivative Terms
        SUBJCBN_BETA2_TERMS_B(14, qd ) = -JCBN_n_ARRAY(1)
        SUBJCBN_BETA2_TERMS_B(15, qd ) = OneThird * R_SQUARE(rd) * ( COTAN_VAL(td) - 4.0_idp * JCBN_n_ARRAY(2) )
        SUBJCBN_BETA2_TERMS_B(16, qd ) = - (1.0_idp/RSIN_SQUARE(rd, td)) * JCBN_n_ARRAY(3)


        ! J_{5,4lmn} Non-Derivative Term
        SUBJCBN_BETA2_TERMS_B(17, qd ) = -FourThirds*R_SQUARE(rd)*SIN_VAL(td)*COS_VAL(td)* JCBN_n_ARRAY(3)

        ! J_{5,4lmn} Derivative Terms
        SUBJCBN_BETA2_TERMS_B(18, qd ) = 0.0_idp
        SUBJCBN_BETA2_TERMS_B(19, qd ) = -FourThirds * RSIN_SQUARE(rd, td) * JCBN_n_ARRAY(3)
        SUBJCBN_BETA2_TERMS_B(20, qd ) = OneThird * RSIN_SQUARE(rd, td) * ( COTAN_VAL(td) + 3.0_idp * JCBN_n_ARRAY(2)  )






        !
        !   Beta 3 Terms
        !




        ! Reused_Value * kappa_{30}
        SUBJCBN_BETA3_TERMS_B(1, qd ) = 0.0_idp

        ! Reused_Value * kappa_{31}
        SUBJCBN_BETA3_TERMS_B(2, qd ) = REUSED_VALUE                                           &
                                * ( 2.0_idp * RSIN_SQUARE(rd, td)   * CUR_DRV_BETA(1,3,rd,td,pd)    &
                                + 2.0_idp                   * CUR_DRV_BETA(3,1,rd,td,pd)    )

        ! Reused_Value * kappa_{32}
        SUBJCBN_BETA3_TERMS_B(3, qd ) = REUSED_VALUE                                           &
                                * ( 2.0_idp * SIN_SQUARE(td)    * CUR_DRV_BETA(2,3,rd,td,pd)        &
                                + 2.0_idp                   * CUR_DRV_BETA(3,2,rd,td,pd)    )

        ! Reused_Value * kappa_{33}
        SUBJCBN_BETA3_TERMS_B(4, qd ) = REUSED_VALUE * FourThirds                              &
                                * ( 1.0_idp /CUR_R_LOCS(rd)         * CUR_VAL_BETA(1,rd,td,pd)      &
                                + 2.0_idp * COTAN_VAL(td)       * CUR_VAL_BETA(2,rd,td,pd)          &
                                - 1.0_idp                   * CUR_DRV_BETA(1,1,rd,td,pd)            &
                                - 1.0_idp                   * CUR_DRV_BETA(2,2,rd,td,pd)            &
                                + 2.0_idp                   * CUR_DRV_BETA(3,3,rd,td,pd)    )


        ! J_{2,5lmn}
        SUBJCBN_BETA3_TERMS_B(5:8, qd ) = -7.0_idp * PSI_POWER(1)                              &
                                                * SUBJCBN_BETA3_TERMS_B(1:4, qd)










        ! J_{3,5lmn} Non-Derivative Term
        SUBJCBN_BETA3_TERMS_B(9, qd ) =  (2.0_idp/3.0_idp)*COTAN_VAL(td)* JCBN_n_ARRAY(1)

        ! J_{3,5lmn} Derivative Terms
        SUBJCBN_BETA3_TERMS_B(10, qd ) = -JCBN_n_ARRAY(3) + COTAN_VAL(td)/3.0_idp
        SUBJCBN_BETA3_TERMS_B(11, qd ) = 0.0_idp
        SUBJCBN_BETA3_TERMS_B(12, qd ) = -TwoThirds * JCBN_n_ARRAY(1)


        ! J_{4,5lmn} Non-Derivative Term
        SUBJCBN_BETA3_TERMS_B(13, qd ) = ( RSIN_SQUARE(rd,td)*SIN_VAL(td)*COS_VAL(td) - R_SQUARE(rd)*COTAN_VAL(td) )    &
                                    * JCBN_n_ARRAY(3)

        ! J_{4,5lmn} Derivative Terms
        SUBJCBN_BETA3_TERMS_B(14, qd ) = 0.0_idp
        SUBJCBN_BETA3_TERMS_B(15, qd ) = FourThirds * R_SQUARE(rd) * JCBN_n_ARRAY(3)
        SUBJCBN_BETA3_TERMS_B(16, qd ) = -FourThirds * R_SQUARE(rd) * JCBN_n_ARRAY(2)


        ! J_{5,5lmn} Non-Derivative Term
        SUBJCBN_BETA3_TERMS_B(17, qd ) = (R_CUBE(rd)*SIN_SQUARE(td)*SIN_SQUARE(td) - R_INVERSE(rd) )  &
                                    * JCBN_n_ARRAY(1)                                           &
                              + ( RSIN_SQUARE(rd, td)*SIN_VAL(td)*COS_VAL(td) - R_SQUARE(rd)*COTAN_VAL(td) )    &
                                    * JCBN_n_ARRAY(2)
        ! J_{5,5lmn} Derivative Terms
        SUBJCBN_BETA3_TERMS_B(18, qd ) = JCBN_n_ARRAY(1)
        SUBJCBN_BETA3_TERMS_B(19, qd ) = R_SQUARE(rd) * JCBN_n_ARRAY(2)
        SUBJCBN_BETA3_TERMS_B(20, qd ) = FourThirds * R_SQUARE(rd) * JCBN_n_ARRAY(3)







        END DO ! pd loop
    END DO  ! td loop
END DO  ! rd loop


!$OMP END DO

!$OMP END PARALLEL

END SUBROUTINE Calc_3D_SubJcbn_Terms_C










!+205+###########################################################################!
!                                                                                !
!                  CREATE_3D_JCBN_MATRIX          !
!                                                                                !
!################################################################################!
SUBROUTINE CREATE_3D_JCBN_MATRIX(   Local_re, Local_te, Local_pe,                   &
                                    Global_re, Global_te, Global_pe,                &
                                    TWOOVER_DELTAR,                                 &
                                    A_FACTORS,                                      &
                                    CC_Ylm_Term,                                    &
                                    SUBJCBN_PSI_TERMS, SUBJCBN_ALPHAPSI_TERMS,      &
                                    SUBJCBN_BETA1_TERMS, SUBJCBN_BETA2_TERMS,       &
                                    SUBJCBN_BETA3_TERMS,                            &
                                    Int_Factor                                      )



INTEGER, INTENT(IN)                                                     ::  Local_re, Local_te, Local_pe,       &
                                                                            Global_re, Global_te, Global_pe



REAL(KIND = idp), INTENT(IN)                                            ::  TWOOVER_DELTAR



COMPLEX(KIND = idp), INTENT(IN), DIMENSION( 1:6,                    &
                                            1:NUM_T_QUAD_POINTS,    &
                                            1:NUM_P_QUAD_POINTS,    &
                                            0:LM_LENGTH             )   :: A_FACTORS

COMPLEX(KIND = idp), INTENT(IN), DIMENSION( 1:NUM_T_QUAD_POINTS,    &
                                            1:NUM_P_QUAD_POINTS,    &
                                            0:LM_LENGTH             )   :: CC_Ylm_Term

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:14,                       &
                                         1:NUM_R_QUAD_POINTS,        &
                                         1:NUM_T_QUAD_POINTS,        &
                                         1:NUM_P_QUAD_POINTS         )  ::  SUBJCBN_PSI_TERMS

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:14,                       &
                                         1:NUM_R_QUAD_POINTS,        &
                                         1:NUM_T_QUAD_POINTS,        &
                                         1:NUM_P_QUAD_POINTS         )  ::  SUBJCBN_ALPHAPSI_TERMS

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:20,                       &
                                         1:NUM_R_QUAD_POINTS,        &
                                         1:NUM_T_QUAD_POINTS,        &
                                         1:NUM_P_QUAD_POINTS         )  ::  SUBJCBN_BETA1_TERMS

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:20,                       &
                                         1:NUM_R_QUAD_POINTS,        &
                                         1:NUM_T_QUAD_POINTS,        &
                                         1:NUM_P_QUAD_POINTS         )  ::  SUBJCBN_BETA2_TERMS

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:20,                       &
                                         1:NUM_R_QUAD_POINTS,        &
                                         1:NUM_T_QUAD_POINTS,        &
                                         1:NUM_P_QUAD_POINTS         )  ::  SUBJCBN_BETA3_TERMS




REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_R_QUAD_POINTS,        &
                                        1:NUM_T_QUAD_POINTS,        &
                                        1:NUM_P_QUAD_POINTS         )   ::  Int_Factor









INTEGER                                                                 ::  pd, td, rd,     &
                                                                            l, m, d,        &
                                                                            lp, mp, dp



INTEGER                                                                 ::  Current_i_Location,         &
                                                                            Current_j_Location




INTEGER                                                                 ::  SPARSE_SHIFT, L_SHIFT
INTEGER                                                                 ::  lm_loc, lpmp_loc

INTEGER                                                                 ::  MATVEC_LOC
INTEGER                                                                 ::  iloc, jloc


REAL(KIND = idp)                                                        ::  Common_Term_A,  &
                                                                            Common_Term_B,  &
                                                                            Common_Term_C

COMPLEX(KIND = idp)                                                     ::  REUSED_VALUE_A, &
                                                                            REUSED_VALUE_B, &
                                                                            REUSED_VALUE_C


COMPLEX(KIND = idp), DIMENSION(1:25)                                    ::  Uncommon_Terms





REAL(KIND = idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:2)                   ::  Current_Lag_Polys




INTEGER                                                                 ::  i, ierr





L_SHIFT = NUM_CFA_VARS*LM_LENGTH



!$OMP PARALLEL DEFAULT(none)                                            &
!$OMP PRIVATE( pd, td, rd, l, m, d, lp, mp, dp, i, ierr,          &
!$OMP           Current_i_Location, Current_j_Location,                 &
!$OMP           MATVEC_LOC, iloc, jloc,                                 &
!$OMP           SPARSE_SHIFT, lm_loc, lpmp_loc,                         &
!$OMP           Common_Term_A, Common_Term_B, Common_Term_C,            &
!$OMP           Current_Lag_Polys,                                      &
!$OMP           Uncommon_Terms                                      )   &
!$OMP SHARED(   Local_re, Local_te, Local_pe,                           &
!$OMP           Global_re, Global_te, Global_pe,                        &
!$OMP           TWOOVER_DELTAR,                                         &
!$OMP           A_FACTORS, CC_Ylm_Term,                                 &
!$OMP           SUBJCBN_PSI_TERMS, SUBJCBN_ALPHAPSI_TERMS,              &
!$OMP           SUBJCBN_BETA1_TERMS, SUBJCBN_BETA2_TERMS,               &
!$OMP           SUBJCBN_BETA3_TERMS,                                    &
!$OMP           LPT_LPT, Int_Factor,                                    &
!$OMP           BLOCK_ELEM_STF_MATVEC,                                  &
!$OMP           NUM_OFF_DIAGONALS,                                      &
!$OMP           ELEM_PROB_DIM,                                          &
!$OMP           Ylm_Table_Block,                                        &
!$OMP           OneThird,                                               &
!$OMP           M_VALUES,                                               &
!$OMP           POSEIDON_COMM_WORLD, myID,                              &
!$OMP           LM_Length, L_SHIFT                                  )



!$OMP DO SCHEDULE(dynamic), COLLAPSE(3)         !, REDUCTION(+:BLOCK_ELEM_STF_MATVEC)
DO dp = 0,DEGREE

    DO d = 0,DEGREE

        DO lp = 0,L_LIMIT


            Current_Lag_Polys(:,0) = LPT_LPT(:, d, dp, 0, 0 )
            Current_Lag_Polys(:,1) = LPT_LPT(:, d, dp, 0, 1 )               &
                                    * TWOOVER_DELTAR
            Current_Lag_Polys(:,2) = LPT_LPT(:, d, dp, 0, 2 )               &
                                    * TWOOVER_DELTAR                     &
                                    * TWOOVER_DELTAR



            DO mp = -lp,lp


                lpmp_loc = lp*(lp+1)+mp

!                Current_i_Location = (Local_re*DEGREE + dp)*L_SHIFT  &
!                                    + NUM_CFA_VARS *lpmp_loc

!                SPARSE_SHIFT = NUM_OFF_DIAGONALS - Current_i_Location
!                SPARSE_SHIFT = - Current_i_Location

                iloc = dp*L_SHIFT + NUM_CFA_VARS * lpmp_loc


                DO l = 0,L_LIMIT

                    DO m = -l,l

                        lm_loc = l*(l+1) + m
!                        Current_j_Location = (Local_re*DEGREE + d) * L_SHIFT      &
!                                           + (l*(l+1) + m) * 5              &
!                                           + SPARSE_SHIFT


!                        Current_j_Location = (Local_re*DEGREE + d) * L_SHIFT      &
!                                           + NUM_CFA_VARS * (l*(l+1) + m)

                        jloc = d*L_SHIFT + NUM_CFA_VARS *(l*(l+1) + m)










                        Uncommon_Terms = 0.0_idp


                        DO pd = 1,NUM_P_QUAD_POINTS
                            DO td = 1,NUM_T_QUAD_POINTS
                                DO rd = 1,NUM_R_QUAD_POINTS



                                    Common_Term_A = CC_Ylm_Term(td, pd, lm_loc)          &
                                                * Current_Lag_Polys(rd,0)                   &
                                                * Int_Factor(rd, td, pd)

                                    Common_Term_B = CC_Ylm_Term(td, pd, lm_loc)          &
                                                * Current_Lag_Polys(rd,1)                   &
                                                * Int_Factor(rd, td, pd)

                                    Common_Term_C = CC_Ylm_Term(td, pd, lm_loc)          &
                                                * Current_Lag_Polys(rd,2)                   &
                                                * Int_Factor(rd, td, pd)







                                    Uncommon_Terms(1) = Uncommon_Terms(1)                                       &
                                                      + SUBJCBN_PSI_TERMS(1, rd, td, pd)                        &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_A



                                    Uncommon_Terms(2) = Uncommon_Terms(2)                                       &
                                                      + SUBJCBN_ALPHAPSI_TERMS(1, rd, td, pd)                   &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_A

                                    Uncommon_Terms(3) = Uncommon_Terms(3)                                     &
                                                       + SUBJCBN_BETA1_TERMS(1, rd, td, pd)                     &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + SUBJCBN_BETA1_TERMS(2, rd, td, pd)                     &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_B      &
                                                       + SUBJCBN_BETA1_TERMS(3, rd, td, pd)                     &
                                                            * A_Factors(2, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + SUBJCBN_BETA1_TERMS(4, rd, td, pd)                     &
                                                            * A_Factors(3, td, pd, lpmp_loc) * Common_Term_A

                                    Uncommon_Terms(4) = Uncommon_Terms(4)                                     &
                                                       + SUBJCBN_BETA2_TERMS(1, rd, td, pd)                     &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + SUBJCBN_BETA2_TERMS(2, rd, td, pd)                     &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_B      &
                                                       + SUBJCBN_BETA2_TERMS(3, rd, td, pd)                     &
                                                            * A_Factors(2, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + SUBJCBN_BETA2_TERMS(4, rd, td, pd)                     &
                                                            * A_Factors(3, td, pd, lpmp_loc) * Common_Term_A

                                    Uncommon_Terms(5) = Uncommon_Terms(5)                                     &
                                                       + SUBJCBN_BETA3_TERMS(1, rd, td, pd)                     &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + SUBJCBN_BETA3_TERMS(2, rd, td, pd)                     &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_B      &
                                                       + SUBJCBN_BETA3_TERMS(3, rd, td, pd)                     &
                                                            * A_Factors(2, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + SUBJCBN_BETA3_TERMS(4, rd, td, pd)                     &
                                                            * A_Factors(3, td, pd, lpmp_loc) * Common_Term_A













                                    Uncommon_Terms(6) = Uncommon_Terms(6)                                       &
                                                      + SUBJCBN_PSI_TERMS(2, rd, td, pd)                        &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_A

                                    Uncommon_Terms(7) = Uncommon_Terms(7)                                       &
                                                      + SUBJCBN_ALPHAPSI_TERMS(2, rd, td, pd)                   &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_A

                                    Uncommon_Terms(8) = Uncommon_Terms(8)                                     &
                                                       + SUBJCBN_BETA1_TERMS(5, rd, td, pd)                     &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + SUBJCBN_BETA1_TERMS(6, rd, td, pd)                     &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_B      &
                                                       + SUBJCBN_BETA1_TERMS(7, rd, td, pd)                     &
                                                            * A_Factors(2, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + SUBJCBN_BETA1_TERMS(8, rd, td, pd)                     &
                                                            * A_Factors(3, td, pd, lpmp_loc) * Common_Term_A

                                    Uncommon_Terms(9) = Uncommon_Terms(9)                                     &
                                                       + SUBJCBN_BETA2_TERMS(5, rd, td, pd)                     &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + SUBJCBN_BETA2_TERMS(6, rd, td, pd)                     &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_B      &
                                                       + SUBJCBN_BETA2_TERMS(7, rd, td, pd)                     &
                                                            * A_Factors(2, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + SUBJCBN_BETA2_TERMS(8, rd, td, pd)                     &
                                                            * A_Factors(3, td, pd, lpmp_loc) * Common_Term_A

                                    Uncommon_Terms(10) = Uncommon_Terms(10)                                     &
                                                       + SUBJCBN_BETA3_TERMS(5, rd, td, pd)                     &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + SUBJCBN_BETA3_TERMS(6, rd, td, pd)                     &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_B      &
                                                       + SUBJCBN_BETA3_TERMS(7, rd, td, pd)                     &
                                                            * A_Factors(2, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + SUBJCBN_BETA3_TERMS(8, rd, td, pd)                     &
                                                            * A_Factors(3, td, pd, lpmp_loc) * Common_Term_A















                                    Uncommon_Terms(11) = Uncommon_Terms(11)                                       &
                                                      + SUBJCBN_PSI_TERMS(3, rd, td, pd)                        &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_A      &
                                                      + SUBJCBN_PSI_TERMS(4, rd, td, pd)                        &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_B      &
                                                      + SUBJCBN_PSI_TERMS(5, rd, td, pd)                        &
                                                            * A_Factors(2, td, pd, lpmp_loc) * Common_Term_A      &
                                                      + SUBJCBN_PSI_TERMS(6, rd, td, pd)                        &
                                                            * A_Factors(3, td, pd, lpmp_loc) * Common_Term_A

                                    Uncommon_Terms(12) = Uncommon_Terms(12)                                       &
                                                      + SUBJCBN_ALPHAPSI_TERMS(3, rd, td, pd)                   &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_A      &
                                                      + SUBJCBN_ALPHAPSI_TERMS(4, rd, td, pd)                   &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_B      &
                                                      + SUBJCBN_ALPHAPSI_TERMS(5, rd, td, pd)                   &
                                                            * A_Factors(2, td, pd, lpmp_loc) * Common_Term_A      &
                                                      + SUBJCBN_ALPHAPSI_TERMS(6, rd, td, pd)                   &
                                                            * A_Factors(3, td, pd, lpmp_loc) * Common_Term_A

                                    Uncommon_Terms(13) = Uncommon_Terms(13)                                     &
                                                       + SUBJCBN_BETA1_TERMS(9, rd, td, pd)                     &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + SUBJCBN_BETA1_TERMS(10, rd, td, pd)                    &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_B      &
                                                       + SUBJCBN_BETA1_TERMS(11, rd, td, pd)                    &
                                                            * A_Factors(2, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + SUBJCBN_BETA1_TERMS(12, rd, td, pd)                    &
                                                            * A_Factors(3, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + OneThird * A_FACTORS(1,td, pd, lpmp_loc) * Common_Term_C

                                    Uncommon_Terms(14) = Uncommon_Terms(14)                                     &
                                                       + SUBJCBN_BETA2_TERMS(9, rd, td, pd)                     &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + SUBJCBN_BETA2_TERMS(10, rd, td, pd)                    &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_B      &
                                                       + SUBJCBN_BETA2_TERMS(11, rd, td, pd)                    &
                                                            * A_Factors(2, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + OneThird * A_FACTORS(2,td, pd, lpmp_loc) * Common_Term_B


                                    Uncommon_Terms(15) = Uncommon_Terms(15)                                     &
                                                       + SUBJCBN_BETA3_TERMS(9, rd, td, pd)                     &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + SUBJCBN_BETA3_TERMS(11, rd, td, pd)                    &
                                                            * A_Factors(2, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + SUBJCBN_BETA3_TERMS(12, rd, td, pd)                    &
                                                            * A_Factors(3, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + OneThird * A_FACTORS(3,td, pd, lpmp_loc) * Common_Term_B








                                    Uncommon_Terms(16) = Uncommon_Terms(16)                       &
                                                      + SUBJCBN_PSI_TERMS(7, rd, td, pd)        &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_A      &
                                                      + SUBJCBN_PSI_TERMS(8, rd, td, pd)                        &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_B      &
                                                      + SUBJCBN_PSI_TERMS(9, rd, td, pd)                        &
                                                            * A_Factors(2, td, pd, lpmp_loc) * Common_Term_A      &
                                                      + SUBJCBN_PSI_TERMS(10, rd, td, pd)                       &
                                                            * A_Factors(3, td, pd, lpmp_loc) * Common_Term_A


                                    Uncommon_Terms(17) = Uncommon_Terms(17)                                       &
                                                      + SUBJCBN_ALPHAPSI_TERMS(7, rd, td, pd)                   &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_A      &
                                                      + SUBJCBN_ALPHAPSI_TERMS(8, rd, td, pd)                   &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_B      &
                                                      + SUBJCBN_ALPHAPSI_TERMS(9, rd, td, pd)                   &
                                                            * A_Factors(2, td, pd, lpmp_loc) * Common_Term_A      &
                                                      + SUBJCBN_ALPHAPSI_TERMS(10, rd, td, pd)                  &
                                                            * A_Factors(3, td, pd, lpmp_loc) * Common_Term_A

                                    Uncommon_Terms(18) = Uncommon_Terms(18)                                     &
                                                       + SUBJCBN_BETA1_TERMS(13, rd, td, pd)                    &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + SUBJCBN_BETA1_TERMS(14, rd, td, pd)                    &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_B      &
                                                       + SUBJCBN_BETA1_TERMS(16, rd, td, pd)                    &
                                                            * A_Factors(3, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + OneThird * A_FACTORS(2,td, pd, lpmp_loc) * Common_Term_B

                                    Uncommon_Terms(19) = Uncommon_Terms(19)                                     &
                                                       + SUBJCBN_BETA2_TERMS(13, rd, td, pd)                    &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + SUBJCBN_BETA2_TERMS(14, rd, td, pd)                    &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_B      &
                                                       + SUBJCBN_BETA2_TERMS(15, rd, td, pd)                    &
                                                            * A_Factors(2, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + SUBJCBN_BETA2_TERMS(16, rd, td, pd)                    &
                                                            * A_Factors(3, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + OneThird * A_FACTORS(4,td, pd, lpmp_loc) * Common_Term_A

                                    Uncommon_Terms(20) = Uncommon_Terms(20)                                     &
                                                       + SUBJCBN_BETA3_TERMS(13, rd, td, pd)                    &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + SUBJCBN_BETA3_TERMS(14, rd, td, pd)                    &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_B      &
                                                       + SUBJCBN_BETA3_TERMS(15, rd, td, pd)                    &
                                                            * A_Factors(2, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + OneThird * A_FACTORS(5,td, pd, lpmp_loc) * Common_Term_A






                                    Uncommon_Terms(21) = Uncommon_Terms(21)                                       &
                                                      + SUBJCBN_PSI_TERMS(11, rd, td, pd)                       &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_A      &
                                                      + SUBJCBN_PSI_TERMS(12, rd, td, pd)                       &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_B      &
                                                      + SUBJCBN_PSI_TERMS(13, rd, td, pd)                       &
                                                            * A_Factors(2, td, pd, lpmp_loc) * Common_Term_A      &
                                                      + SUBJCBN_PSI_TERMS(14, rd, td, pd)                       &
                                                            * A_Factors(3, td, pd, lpmp_loc) * Common_Term_A

                                    Uncommon_Terms(22) = Uncommon_Terms(22)                                     &
                                                       + SUBJCBN_ALPHAPSI_TERMS(11, rd, td, pd)                 &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + SUBJCBN_ALPHAPSI_TERMS(12, rd, td, pd)                 &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_B      &
                                                       + SUBJCBN_ALPHAPSI_TERMS(13, rd, td, pd)                 &
                                                            * A_Factors(2, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + SUBJCBN_ALPHAPSI_TERMS(14, rd, td, pd)                 &
                                                            * A_Factors(3, td, pd, lpmp_loc) * Common_Term_A

                                    Uncommon_Terms(23) = Uncommon_Terms(23)                                     &
                                                       + SUBJCBN_BETA1_TERMS(17, rd, td, pd)                    &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + SUBJCBN_BETA1_TERMS(19, rd, td, pd)                    &
                                                            * A_Factors(2, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + SUBJCBN_BETA1_TERMS(20, rd, td, pd)                    &
                                                            * A_Factors(3, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + OneThird * A_FACTORS(3,td, pd, lpmp_loc) * Common_Term_B

                                    Uncommon_Terms(24) = Uncommon_Terms(24)                                     &
                                                       + SUBJCBN_BETA2_TERMS(17, rd, td, pd)                    &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + SUBJCBN_BETA2_TERMS(18, rd, td, pd)                    &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_B      &
                                                       + SUBJCBN_BETA2_TERMS(19, rd, td, pd)                    &
                                                            * A_Factors(2, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + OneThird * A_FACTORS(5,td, pd, lpmp_loc) * Common_Term_A

                                    Uncommon_Terms(25) = Uncommon_Terms(25)                                     &
                                                       + SUBJCBN_BETA3_TERMS(17, rd, td, pd)                    &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + SUBJCBN_BETA3_TERMS(18, rd, td, pd)                    &
                                                            * A_Factors(1, td, pd, lpmp_loc) * Common_Term_B      &
                                                       + SUBJCBN_BETA3_TERMS(19, rd, td, pd)                    &
                                                            * A_Factors(2, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + SUBJCBN_BETA3_TERMS(20, rd, td, pd)                    &
                                                            * A_Factors(3, td, pd, lpmp_loc) * Common_Term_A      &
                                                       + OneThird * A_FACTORS(6,td, pd, lpmp_loc) * Common_Term_A





                                END DO  ! rd Loop
                            END DO  ! td Loop
                        END DO  ! pd Loop


!                        DO i = 0,NUM_SHELLS
!
!                            IF (myID_PETSC == i ) THEN
!
!                                PRINT*,"myID",myID,"Global Elements",Global_re, Global_te, Global_pe
!                                PRINT*,"l,m,d,lp,mp,dp",l,m,d,lp,mp,dp
!                                PRINT*,Uncommon_Terms
!
!                            END IF
!                            CALL MPI_BARRIER(POSEIDON_COMM_PETSC,ierr)
!                       END DO


!                        MATVEC_LOC = jloc*ELEM_PROB_DIM + iloc
!                        PRINT*,"Before",myID,MATVEC_LOC,MATVEC_LOC+4
!                        PRINT*,BLOCK_ELEM_STF_MATVEC(MATVEC_LOC:MATVEC_LOC+4,Local_re)
!                        PRINT*,"   "

                        MATVEC_LOC = jloc*ELEM_PROB_DIM + iloc
                        BLOCK_ELEM_STF_MATVEC(MATVEC_LOC:MATVEC_LOC+4,Local_re)                             &
                                = BLOCK_ELEM_STF_MATVEC(MATVEC_LOC:MATVEC_LOC+4,Local_re)                   &
                                + Uncommon_Terms(1:5)

                        MATVEC_LOC = (jloc+1)*ELEM_PROB_DIM + iloc
                        BLOCK_ELEM_STF_MATVEC(MATVEC_LOC:MATVEC_LOC+4,Local_re)                             &
                                = BLOCK_ELEM_STF_MATVEC(MATVEC_LOC:MATVEC_LOC+4,Local_re)                   &
                                + Uncommon_Terms(6:10)

                        MATVEC_LOC = (jloc+2)*ELEM_PROB_DIM + iloc
                        BLOCK_ELEM_STF_MATVEC(MATVEC_LOC:MATVEC_LOC+4,Local_re)                             &
                                = BLOCK_ELEM_STF_MATVEC(MATVEC_LOC:MATVEC_LOC+4,Local_re)                   &
                                + Uncommon_Terms(11:15)

                        MATVEC_LOC = (jloc+3)*ELEM_PROB_DIM + iloc
                        BLOCK_ELEM_STF_MATVEC(MATVEC_LOC:MATVEC_LOC+4,Local_re)                             &
                                = BLOCK_ELEM_STF_MATVEC(MATVEC_LOC:MATVEC_LOC+4,Local_re)                   &
                                + Uncommon_Terms(16:20)

                        MATVEC_LOC = (jloc+4)*ELEM_PROB_DIM + iloc
                        BLOCK_ELEM_STF_MATVEC(MATVEC_LOC:MATVEC_LOC+4,Local_re)                             &
                                = BLOCK_ELEM_STF_MATVEC(MATVEC_LOC:MATVEC_LOC+4,Local_re)                   &
                                + Uncommon_Terms(21:25)



                    END DO  ! m Loop
                END DO  ! l Loop
            END DO  ! mp Loop
        END DO  ! lp Loop
    END DO  ! dp Loop
END DO  ! d Loop
!$OMP END DO

!$OMP END PARALLEL




END SUBROUTINE CREATE_3D_JCBN_MATRIX








!+205e+###########################################################################!
!                                                                                !
!                  CREATE_3D_JCBN_MATRIX          !
!                                                                                !
!################################################################################!
SUBROUTINE CREATE_3D_JCBN_MATRIX_E( Local_re, Local_te, Local_pe,                   &
                                    Global_re, Global_te, Global_pe,                &
                                    TWOOVER_DELTAR,                                 &
                                    A_FACTORS,                                      &
                                    CC_Ylm_Term,                                    &
                                    SUBJCBN_PSI_TERMS_B, SUBJCBN_ALPHAPSI_TERMS_B,      &
                                    SUBJCBN_BETA1_TERMS_B, SUBJCBN_BETA2_TERMS_B,       &
                                    SUBJCBN_BETA3_TERMS_B,                            &
                                    Int_Factor                                      )



INTEGER, INTENT(IN)                                                     ::  Local_re, Local_te, Local_pe,       &
                                                                            Global_re, Global_te, Global_pe



REAL(KIND = idp), INTENT(IN)                                            ::  TWOOVER_DELTAR




COMPLEX(KIND = idp), INTENT(IN), DIMENSION( 1:6,                    &
                                            1:NUM_T_QUAD_POINTS,    &
                                            1:NUM_P_QUAD_POINTS,    &
                                            0:LM_LENGTH             )   :: A_FACTORS

COMPLEX(KIND = idp), INTENT(IN), DIMENSION( 1:NUM_T_QUAD_POINTS,    &
                                            1:NUM_P_QUAD_POINTS,    &
                                            0:LM_LENGTH             )   :: CC_Ylm_Term

REAL(KIND = idp), INTENT(IN), DIMENSION( NUM_QUAD_DOF, 1:14 )           ::  SUBJCBN_PSI_TERMS_B

REAL(KIND = idp), INTENT(IN), DIMENSION( NUM_QUAD_DOF, 1:14 )           ::  SUBJCBN_ALPHAPSI_TERMS_B

REAL(KIND = idp), INTENT(IN), DIMENSION( NUM_QUAD_DOF, 1:20 )           ::  SUBJCBN_BETA1_TERMS_B

REAL(KIND = idp), INTENT(IN), DIMENSION( NUM_QUAD_DOF, 1:20 )           ::  SUBJCBN_BETA2_TERMS_B

REAL(KIND = idp), INTENT(IN), DIMENSION( NUM_QUAD_DOF, 1:20 )           ::  SUBJCBN_BETA3_TERMS_B



REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_R_QUAD_POINTS,        &
                                        1:NUM_T_QUAD_POINTS,        &
                                        1:NUM_P_QUAD_POINTS         )   ::  Int_Factor









INTEGER                                                                 ::  pd, td, rd,     &
                                                                            l, m, d,        &
                                                                            lp, mp, dp,     &
                                                                            qd



INTEGER                                                                 ::  Current_i_Location,         &
                                                                            Current_j_Location




INTEGER                                                                 ::  SPARSE_SHIFT, L_SHIFT
INTEGER                                                                 ::  lm_loc, lpmp_loc

INTEGER                                                                 ::  MATVEC_LOC
INTEGER                                                                 ::  iloc, jloc


REAL(KIND = idp), DIMENSION(1:NUM_R_QUAD_POINTS)                        ::  Common_Term_A,  &
                                                                            Common_Term_B,  &
                                                                            Common_Term_C



COMPLEX(KIND = idp), DIMENSION(1:25)                                    ::  Uncommon_Terms

COMPLEX(KIND = idp), DIMENSION(1:NUM_QUAD_DOF, 1:10)                    ::  Integral_Term


REAL(KIND = idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:2)                   ::  Current_Lag_Polys

COMPLEX(KIND = idp), DIMENSION(1:6)                                     ::  Ylm_Cross_Terms



INTEGER                                                                 ::  i, ierr





L_SHIFT = NUM_CFA_VARS*LM_LENGTH



!$OMP PARALLEL DEFAULT(none)                                            &
!$OMP PRIVATE( pd, td, rd, l, m, d, lp, mp, dp, i, qd,  ierr,          &
!$OMP           Current_i_Location, Current_j_Location,                 &
!$OMP           MATVEC_LOC, iloc, jloc,                                 &
!$OMP           SPARSE_SHIFT, lm_loc, lpmp_loc,                         &
!$OMP           Common_Term_A, Common_Term_B, Common_Term_C,            &
!$OMP           Current_Lag_Polys,                                      &
!$OMP           Ylm_Cross_Terms,                                        &
!$OMP           Integral_Term,                                          &
!$OMP           Uncommon_Terms                                      )   &
!$OMP SHARED(   Local_re, Local_te, Local_pe,                           &
!$OMP           Global_re, Global_te, Global_pe,                        &
!$OMP           TWOOVER_DELTAR,                                         &
!$OMP           CC_Ylm_Term,                                            &
!$OMP           A_FACTORS,                                              &
!$OMP           SUBJCBN_PSI_TERMS_B, SUBJCBN_ALPHAPSI_TERMS_B,              &
!$OMP           SUBJCBN_BETA1_TERMS_B, SUBJCBN_BETA2_TERMS_B,               &
!$OMP           SUBJCBN_BETA3_TERMS_B,                                    &
!$OMP           LPT_LPT, Int_Factor,                                    &
!$OMP           BLOCK_ELEM_STF_MATVEC,                                  &
!$OMP           NUM_OFF_DIAGONALS,                                      &
!$OMP           ELEM_PROB_DIM,                                          &
!$OMP           Ylm_Table_Block,                                        &
!$OMP           OneThird,                                               &
!$OMP           M_VALUES,                                               &
!$OMP           POSEIDON_COMM_WORLD, myID,                              &
!$OMP           LM_Length, L_SHIFT                                  )



!$OMP DO SCHEDULE(dynamic), COLLAPSE(3)         !, REDUCTION(+:BLOCK_ELEM_STF_MATVEC)
DO dp = 0,DEGREE

    DO d = 0,DEGREE

        DO lpmp_loc = 0,LM_Length - 1

            Current_Lag_Polys(:,0) = LPT_LPT(:, d, dp, 0, 0 )
            Current_Lag_Polys(:,1) = LPT_LPT(:, d, dp, 0, 1 )            &
                                    * TWOOVER_DELTAR
            Current_Lag_Polys(:,2) = LPT_LPT(:, d, dp, 0, 2 )            &
                                    * TWOOVER_DELTAR                     &
                                    * TWOOVER_DELTAR

            iloc = dp*L_SHIFT + NUM_CFA_VARS * lpmp_loc


            DO lm_loc = 0,LM_Length - 1


                jloc = d*L_SHIFT + NUM_CFA_VARS * lm_loc


                DO pd = 1,NUM_P_QUAD_POINTS
                    DO td = 1,NUM_T_QUAD_POINTS

                        qd = ((pd-1) * NUM_T_QUAD_POINTS + td-1 ) * NUM_R_QUAD_POINTS

                        Integral_Term(qd+1:qd+NUM_R_QUAD_POINTS,1)                                          &
                                            = Current_Lag_Polys(:,0) * Int_Factor(:, td, pd)                &
                                            * A_Factors(1, td, pd, lpmp_loc)*CC_Ylm_Term(td, pd, lm_loc)

                        Integral_Term(qd+1:qd+NUM_R_QUAD_POINTS,2)                                          &
                                            = Current_Lag_Polys(:,1) * Int_Factor(:, td, pd)                &
                                            * A_Factors(1, td, pd, lpmp_loc)*CC_Ylm_Term(td, pd, lm_loc)

                        Integral_Term(qd+1:qd+NUM_R_QUAD_POINTS,3)                                          &
                                            = Current_Lag_Polys(:,0) * Int_Factor(:, td, pd)                &
                                            * A_Factors(2, td, pd, lpmp_loc)*CC_Ylm_Term(td, pd, lm_loc)

                        Integral_Term(qd+1:qd+NUM_R_QUAD_POINTS,4)                                          &
                                            = Current_Lag_Polys(:,0) * Int_Factor(:, td, pd)                &
                                            * A_Factors(3, td, pd, lpmp_loc)*CC_Ylm_Term(td, pd, lm_loc)

                        Integral_Term(qd+1:qd+NUM_R_QUAD_POINTS,5)                                          &
                                            = Current_Lag_Polys(:,2) * Int_Factor(:, td, pd)                &
                                            * A_Factors(1, td, pd, lpmp_loc)*CC_Ylm_Term(td, pd, lm_loc)

                        Integral_Term(qd+1:qd+NUM_R_QUAD_POINTS,6)                                          &
                                            = Current_Lag_Polys(:,1) * Int_Factor(:, td, pd)                &
                                            * A_Factors(2, td, pd, lpmp_loc)*CC_Ylm_Term(td, pd, lm_loc)

                        Integral_Term(qd+1:qd+NUM_R_QUAD_POINTS,7)                                          &
                                            = Current_Lag_Polys(:,1) * Int_Factor(:, td, pd)                &
                                            * A_Factors(3, td, pd, lpmp_loc)*CC_Ylm_Term(td, pd, lm_loc)

                        Integral_Term(qd+1:qd+NUM_R_QUAD_POINTS,8)                                          &
                                            = Current_Lag_Polys(:,0) * Int_Factor(:, td, pd)                &
                                            * A_Factors(4, td, pd, lpmp_loc)*CC_Ylm_Term(td, pd, lm_loc)

                        Integral_Term(qd+1:qd+NUM_R_QUAD_POINTS,9)                                          &
                                            = Current_Lag_Polys(:,0) * Int_Factor(:, td, pd)                &
                                            * A_Factors(5, td, pd, lpmp_loc)*CC_Ylm_Term(td, pd, lm_loc)

                        Integral_Term(qd+1:qd+NUM_R_QUAD_POINTS,10)                                         &
                                            = Current_Lag_Polys(:,0) * Int_Factor(:, td, pd)                &
                                            * A_Factors(6, td, pd, lpmp_loc)*CC_Ylm_Term(td, pd, lm_loc)

                    END DO
                END DO






                Uncommon_Terms(1) = DOT_PRODUCT(SUBJCBN_PSI_TERMS_B( :, 1 ) , Integral_Term(:,1) )

                Uncommon_Terms(2) = DOT_PRODUCT( SUBJCBN_ALPHAPSI_TERMS_B( :, 1 ) , Integral_Term(:,1) )


                Uncommon_Terms(3) = DOT_PRODUCT(SUBJCBN_BETA1_TERMS_B( :, 1 ) , Integral_Term(:,1) )    &
                                  + DOT_PRODUCT(SUBJCBN_BETA1_TERMS_B( :, 2 ) , Integral_Term(:,2) )    &
                                  + DOT_PRODUCT(SUBJCBN_BETA1_TERMS_B( :, 3 ) , Integral_Term(:,3) )    &
                                  + DOT_PRODUCT(SUBJCBN_BETA1_TERMS_B( :, 4 ) , Integral_Term(:,4) )


                Uncommon_Terms(4) = DOT_PRODUCT(SUBJCBN_BETA2_TERMS_B( :, 1 ) , Integral_Term(:,1) )    &
                                  + DOT_PRODUCT(SUBJCBN_BETA2_TERMS_B( :, 2 ) , Integral_Term(:,2) )    &
                                  + DOT_PRODUCT(SUBJCBN_BETA2_TERMS_B( :, 3 ) , Integral_Term(:,3) )    &
                                  + DOT_PRODUCT(SUBJCBN_BETA2_TERMS_B( :, 4 ) , Integral_Term(:,4) )

                Uncommon_Terms(5) = DOT_PRODUCT(SUBJCBN_BETA3_TERMS_B( :, 1 ) , Integral_Term(:,1) )    &
                                  + DOT_PRODUCT(SUBJCBN_BETA3_TERMS_B( :, 2 ) , Integral_Term(:,2) )    &
                                  + DOT_PRODUCT(SUBJCBN_BETA3_TERMS_B( :, 3 ) , Integral_Term(:,3) )    &
                                  + DOT_PRODUCT(SUBJCBN_BETA3_TERMS_B( :, 4 ) , Integral_Term(:,4) )









                Uncommon_Terms(6) = DOT_PRODUCT(SUBJCBN_PSI_TERMS_B( :, 2 ) , Integral_Term(:,1) )

                Uncommon_Terms(7) = DOT_PRODUCT(SUBJCBN_ALPHAPSI_TERMS_B( :, 2 ) , Integral_Term(:,1) )

                Uncommon_Terms(8) = DOT_PRODUCT(SUBJCBN_BETA1_TERMS_B( :, 5 ), Integral_Term(:,1) )     &
                                  + DOT_PRODUCT(SUBJCBN_BETA1_TERMS_B( :, 6 ), Integral_Term(:,2) )     &
                                  + DOT_PRODUCT(SUBJCBN_BETA1_TERMS_B( :, 7 ), Integral_Term(:,3) )     &
                                  + DOT_PRODUCT(SUBJCBN_BETA1_TERMS_B( :, 8 ), Integral_Term(:,4) )


                Uncommon_Terms(9) = DOT_PRODUCT(SUBJCBN_BETA2_TERMS_B( :, 5 ), Integral_Term(:,1) )     &
                                  + DOT_PRODUCT(SUBJCBN_BETA2_TERMS_B( :, 6 ), Integral_Term(:,2) )     &
                                  + DOT_PRODUCT(SUBJCBN_BETA2_TERMS_B( :, 7 ), Integral_Term(:,3) )     &
                                  + DOT_PRODUCT(SUBJCBN_BETA2_TERMS_B( :, 8 ), Integral_Term(:,4) )

                Uncommon_Terms(10) = DOT_PRODUCT(SUBJCBN_BETA3_TERMS_B( :, 5 ), Integral_Term(:,1) )    &
                                   + DOT_PRODUCT(SUBJCBN_BETA3_TERMS_B( :, 6 ), Integral_Term(:,2) )    &
                                   + DOT_PRODUCT(SUBJCBN_BETA3_TERMS_B( :, 7 ), Integral_Term(:,3) )    &
                                   + DOT_PRODUCT(SUBJCBN_BETA3_TERMS_B( :, 8 ), Integral_Term(:,4) )







                Uncommon_Terms(11) = DOT_PRODUCT(SUBJCBN_PSI_TERMS_B( :, 3 ), Integral_Term(:,1) )      &
                                   + DOT_PRODUCT(SUBJCBN_PSI_TERMS_B( :, 4 ), Integral_Term(:,2) )      &
                                   + DOT_PRODUCT(SUBJCBN_PSI_TERMS_B( :, 5 ), Integral_Term(:,3) )      &
                                   + DOT_PRODUCT(SUBJCBN_PSI_TERMS_B( :, 6 ), Integral_Term(:,4) )

                Uncommon_Terms(12) = DOT_PRODUCT(SUBJCBN_ALPHAPSI_TERMS_B( :, 3 ), Integral_Term(:,1) ) &
                                   + DOT_PRODUCT(SUBJCBN_ALPHAPSI_TERMS_B( :, 4 ), Integral_Term(:,2) ) &
                                   + DOT_PRODUCT(SUBJCBN_ALPHAPSI_TERMS_B( :, 5 ), Integral_Term(:,3) ) &
                                   + DOT_PRODUCT(SUBJCBN_ALPHAPSI_TERMS_B( :, 6 ), Integral_Term(:,4) )


                Uncommon_Terms(13) = DOT_PRODUCT(SUBJCBN_BETA1_TERMS_B( :, 9 ), Integral_Term(:,1) )    &
                                   + DOT_PRODUCT(SUBJCBN_BETA1_TERMS_B( :, 10), Integral_Term(:,2) )    &
                                   + DOT_PRODUCT(SUBJCBN_BETA1_TERMS_B( :, 11), Integral_Term(:,3) )    &
                                   + DOT_PRODUCT(SUBJCBN_BETA1_TERMS_B( :, 12), Integral_Term(:,4) )    &
                                   + OneThird * SUM( Integral_Term(:,5) )


                Uncommon_Terms(14) = DOT_PRODUCT(SUBJCBN_BETA2_TERMS_B( :, 9 ), Integral_Term(:,1) )    &
                                   + DOT_PRODUCT(SUBJCBN_BETA2_TERMS_B( :, 10), Integral_Term(:,2) )    &
                                   + DOT_PRODUCT(SUBJCBN_BETA2_TERMS_B( :, 11), Integral_Term(:,3) )    &
                                   + DOT_PRODUCT(SUBJCBN_BETA2_TERMS_B( :, 12), Integral_Term(:,4) )    &
                                   + OneThird * SUM( Integral_Term(:,6)  )

                Uncommon_Terms(15) = DOT_PRODUCT(SUBJCBN_BETA3_TERMS_B( :, 9 ), Integral_Term(:,1) )    &
                                   + DOT_PRODUCT(SUBJCBN_BETA3_TERMS_B( :, 11), Integral_Term(:,3) )    &
                                   + DOT_PRODUCT(SUBJCBN_BETA3_TERMS_B( :, 12), Integral_Term(:,4) )    &
                                   + OneThird * SUM( Integral_Term(:,7) )






                Uncommon_Terms(16) = DOT_PRODUCT(SUBJCBN_PSI_TERMS_B( :, 7 ), Integral_Term(:,1) )      &
                                   + DOT_PRODUCT(SUBJCBN_PSI_TERMS_B( :, 8 ), Integral_Term(:,2) )      &
                                   + DOT_PRODUCT(SUBJCBN_PSI_TERMS_B( :, 9 ), Integral_Term(:,3) )      &
                                   + DOT_PRODUCT(SUBJCBN_PSI_TERMS_B( :, 10), Integral_Term(:,4) )

                Uncommon_Terms(17) = DOT_PRODUCT(SUBJCBN_ALPHAPSI_TERMS_B( :, 7 ), Integral_Term(:,1) ) &
                                   + DOT_PRODUCT(SUBJCBN_ALPHAPSI_TERMS_B( :, 8 ), Integral_Term(:,2) ) &
                                   + DOT_PRODUCT(SUBJCBN_ALPHAPSI_TERMS_B( :, 9 ), Integral_Term(:,3) ) &
                                   + DOT_PRODUCT(SUBJCBN_ALPHAPSI_TERMS_B( :, 10), Integral_Term(:,4) )


                Uncommon_Terms(18) = DOT_PRODUCT(SUBJCBN_BETA1_TERMS_B( :, 13), Integral_Term(:,1) )    &
                                   + DOT_PRODUCT(SUBJCBN_BETA1_TERMS_B( :, 14), Integral_Term(:,2) )    &
                                   + DOT_PRODUCT(SUBJCBN_BETA1_TERMS_B( :, 16), Integral_Term(:,4) )    &
                                   + OneThird * SUM( Integral_Term(:,6) )


                Uncommon_Terms(19) = DOT_PRODUCT(SUBJCBN_BETA2_TERMS_B( :, 13), Integral_Term(:,1) )    &
                                   + DOT_PRODUCT(SUBJCBN_BETA2_TERMS_B( :, 14), Integral_Term(:,2) )    &
                                   + DOT_PRODUCT(SUBJCBN_BETA2_TERMS_B( :, 15), Integral_Term(:,3) )    &
                                   + DOT_PRODUCT(SUBJCBN_BETA2_TERMS_B( :, 16), Integral_Term(:,4) )    &
                                   + OneThird * SUM( Integral_Term(:,8)  )

                Uncommon_Terms(20) = DOT_PRODUCT(SUBJCBN_BETA3_TERMS_B( :, 13), Integral_Term(:,1) )    &
                                   + DOT_PRODUCT(SUBJCBN_BETA3_TERMS_B( :, 14), Integral_Term(:,2) )    &
                                   + DOT_PRODUCT(SUBJCBN_BETA3_TERMS_B( :, 15), Integral_Term(:,3) )    &
                                   + OneThird * SUM( Integral_Term(:,9)  )








                Uncommon_Terms(21) = DOT_PRODUCT(SUBJCBN_PSI_TERMS_B( :, 11 ), Integral_Term(:,1) )     &
                                   + DOT_PRODUCT(SUBJCBN_PSI_TERMS_B( :, 12 ), Integral_Term(:,2) )     &
                                   + DOT_PRODUCT(SUBJCBN_PSI_TERMS_B( :, 13 ), Integral_Term(:,3) )     &
                                   + DOT_PRODUCT(SUBJCBN_PSI_TERMS_B( :, 14 ), Integral_Term(:,4) )

                Uncommon_Terms(22) = DOT_PRODUCT(SUBJCBN_ALPHAPSI_TERMS_B( :, 11), Integral_Term(:,1) )&
                                   + DOT_PRODUCT(SUBJCBN_ALPHAPSI_TERMS_B( :, 12), Integral_Term(:,2) )&
                                   + DOT_PRODUCT(SUBJCBN_ALPHAPSI_TERMS_B( :, 13), Integral_Term(:,3) )&
                                   + DOT_PRODUCT(SUBJCBN_ALPHAPSI_TERMS_B( :, 14), Integral_Term(:,4) )


                Uncommon_Terms(23) = DOT_PRODUCT(SUBJCBN_BETA1_TERMS_B( :, 17), Integral_Term(:,1) )   &
                                   + DOT_PRODUCT(SUBJCBN_BETA1_TERMS_B( :, 19), Integral_Term(:,3) )   &
                                   + DOT_PRODUCT(SUBJCBN_BETA1_TERMS_B( :, 20), Integral_Term(:,4) )   &
                                   + OneThird * SUM( Integral_Term(:,7)  )


                Uncommon_Terms(24) = DOT_PRODUCT(SUBJCBN_BETA2_TERMS_B( :, 17 ), Integral_Term(:,1) )   &
                                   + DOT_PRODUCT(SUBJCBN_BETA2_TERMS_B( :, 18 ), Integral_Term(:,2) )   &
                                   + DOT_PRODUCT(SUBJCBN_BETA2_TERMS_B( :, 19 ), Integral_Term(:,3) )   &
                                   + OneThird * SUM( Integral_Term(:,9)  )

                Uncommon_Terms(25) = DOT_PRODUCT(SUBJCBN_BETA3_TERMS_B( :, 17 ), Integral_Term(:,1) )   &
                                   + DOT_PRODUCT(SUBJCBN_BETA3_TERMS_B( :, 18 ), Integral_Term(:,2) )   &
                                   + DOT_PRODUCT(SUBJCBN_BETA3_TERMS_B( :, 19 ), Integral_Term(:,3) )   &
                                   + DOT_PRODUCT(SUBJCBN_BETA3_TERMS_B( :, 20 ), Integral_Term(:,4) )   &
                                   + OneThird * SUM( Integral_Term(:,10)  )




                MATVEC_LOC = jloc*ELEM_PROB_DIM + iloc
                BLOCK_ELEM_STF_MATVEC(MATVEC_LOC:MATVEC_LOC+4,Local_re)                             &
                        = BLOCK_ELEM_STF_MATVEC(MATVEC_LOC:MATVEC_LOC+4,Local_re)                   &
                        + Uncommon_Terms(1:5)

                MATVEC_LOC = (jloc+1)*ELEM_PROB_DIM + iloc
                BLOCK_ELEM_STF_MATVEC(MATVEC_LOC:MATVEC_LOC+4,Local_re)                             &
                        = BLOCK_ELEM_STF_MATVEC(MATVEC_LOC:MATVEC_LOC+4,Local_re)                   &
                        + Uncommon_Terms(6:10)

                MATVEC_LOC = (jloc+2)*ELEM_PROB_DIM + iloc
                BLOCK_ELEM_STF_MATVEC(MATVEC_LOC:MATVEC_LOC+4,Local_re)                             &
                        = BLOCK_ELEM_STF_MATVEC(MATVEC_LOC:MATVEC_LOC+4,Local_re)                   &
                        + Uncommon_Terms(11:15)

                MATVEC_LOC = (jloc+3)*ELEM_PROB_DIM + iloc
                BLOCK_ELEM_STF_MATVEC(MATVEC_LOC:MATVEC_LOC+4,Local_re)                             &
                        = BLOCK_ELEM_STF_MATVEC(MATVEC_LOC:MATVEC_LOC+4,Local_re)                   &
                        + Uncommon_Terms(16:20)

                MATVEC_LOC = (jloc+4)*ELEM_PROB_DIM + iloc
                BLOCK_ELEM_STF_MATVEC(MATVEC_LOC:MATVEC_LOC+4,Local_re)                             &
                        = BLOCK_ELEM_STF_MATVEC(MATVEC_LOC:MATVEC_LOC+4,Local_re)                   &
                        + Uncommon_Terms(21:25)



            END DO  ! lm_loc Loop
        END DO  ! lpmp_loc Loop
    END DO  ! dp Loop
END DO  ! d Loop
!$OMP END DO

!$OMP END PARALLEL




END SUBROUTINE CREATE_3D_JCBN_MATRIX_E

















!+205G+###########################################################################!
!                                                                                !
!                  CREATE_3D_JCBN_MATRIX          !
!                                                                                !
!################################################################################!
SUBROUTINE CREATE_3D_JCBN_MATRIX_G( Local_re, Local_te, Local_pe,                   &
                                    Global_re, Global_te, Global_pe,                &
                                    TWOOVER_DELTAR,                                 &
                                    A_FACTORS,                                      &
                                    CC_Ylm_Term,                                    &
                                    SUBJCBN_PSI_TERMS_B, SUBJCBN_ALPHAPSI_TERMS_B,      &
                                    SUBJCBN_BETA1_TERMS_B, SUBJCBN_BETA2_TERMS_B,       &
                                    SUBJCBN_BETA3_TERMS_B,                            &
                                    Int_Factor                                      )



INTEGER, INTENT(IN)                                                     ::  Local_re, Local_te, Local_pe,       &
                                                                            Global_re, Global_te, Global_pe



REAL(KIND = idp), INTENT(IN)                                            ::  TWOOVER_DELTAR




COMPLEX(KIND = idp), INTENT(IN), DIMENSION( 1:6,                    &
                                            1:NUM_T_QUAD_POINTS,    &
                                            1:NUM_P_QUAD_POINTS,    &
                                            0:LM_LENGTH             )   :: A_FACTORS

COMPLEX(KIND = idp), INTENT(IN), DIMENSION( 1:NUM_T_QUAD_POINTS,    &
                                            1:NUM_P_QUAD_POINTS,    &
                                            0:LM_LENGTH             )   :: CC_Ylm_Term

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:14, 1:NUM_QUAD_DOF )         ::  SUBJCBN_PSI_TERMS_B

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:14, 1:NUM_QUAD_DOF )         ::  SUBJCBN_ALPHAPSI_TERMS_B

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:20, 1:NUM_QUAD_DOF )         ::  SUBJCBN_BETA1_TERMS_B

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:20, 1:NUM_QUAD_DOF )         ::  SUBJCBN_BETA2_TERMS_B

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:20, 1:NUM_QUAD_DOF )         ::  SUBJCBN_BETA3_TERMS_B



REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_R_QUAD_POINTS,        &
                                        1:NUM_T_QUAD_POINTS,        &
                                        1:NUM_P_QUAD_POINTS         )   ::  Int_Factor









INTEGER                                                                 ::  pd, td, rd,     &
                                                                            l, m, d,        &
                                                                            lp, mp, dp,     &
                                                                            qd



INTEGER                                                                 ::  Current_i_Location,         &
                                                                            Current_j_Location




INTEGER                                                                 ::  SPARSE_SHIFT, L_SHIFT
INTEGER                                                                 ::  lm_loc, lpmp_loc

INTEGER                                                                 ::  MATVEC_LOC
INTEGER                                                                 ::  iloc, jloc






COMPLEX(KIND = idp), DIMENSION(1:25)                                    ::  Uncommon_Terms

COMPLEX(KIND = idp), DIMENSION(1:10, 1:NUM_QUAD_DOF)                    ::  Integral_Term



REAL(KIND = idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:2)                   ::  Current_Lag_Polys

COMPLEX(KIND = idp), DIMENSION(1:6)                                     ::  Ylm_Cross_Terms



INTEGER                                                                 ::  i, ierr





L_SHIFT = NUM_CFA_VARS*LM_LENGTH



!$OMP PARALLEL DEFAULT(none)                                            &
!$OMP PRIVATE( pd, td, rd, lm_loc, d, lpmp_loc, dp, i, ierr, qd,        &
!$OMP           Current_i_Location, Current_j_Location,                 &
!$OMP           MATVEC_LOC, iloc, jloc,                                 &
!$OMP           Current_Lag_Polys,                                      &
!$OMP           Ylm_Cross_Terms,                                        &
!$OMP           Integral_Term,                                          &
!$OMP           Uncommon_Terms                                      )   &
!$OMP SHARED(   Local_re, Local_te, Local_pe,                           &
!$OMP           TWOOVER_DELTAR,                                         &
!$OMP           CC_Ylm_Term,                                            &
!$OMP           A_FACTORS,                                              &
!$OMP           SUBJCBN_PSI_TERMS_B, SUBJCBN_ALPHAPSI_TERMS_B,          &
!$OMP           SUBJCBN_BETA1_TERMS_B, SUBJCBN_BETA2_TERMS_B,           &
!$OMP           SUBJCBN_BETA3_TERMS_B,                                  &
!$OMP           LPT_LPT, Int_Factor,                                    &
!$OMP           BLOCK_ELEM_STF_MATVEC,                                  &
!$OMP           ELEM_PROB_DIM,                                          &
!$OMP           Ylm_Table_Block,                                        &
!$OMP           OneThird,                                               &
!$OMP           LM_Length, L_SHIFT                                  )




!$OMP DO SCHEDULE(dynamic), COLLAPSE(3)
DO dp = 0,DEGREE

    DO d = 0,DEGREE

        DO lpmp_loc = 0,LM_Length - 1

            Current_Lag_Polys(:,0) = LPT_LPT(:, d, dp, 0, 0 )
            Current_Lag_Polys(:,1) = LPT_LPT(:, d, dp, 0, 1 )            &
                                    * TWOOVER_DELTAR
            Current_Lag_Polys(:,2) = LPT_LPT(:, d, dp, 0, 2 )            &
                                    * TWOOVER_DELTAR                     &
                                    * TWOOVER_DELTAR

            iloc = dp*L_SHIFT + NUM_CFA_VARS * lpmp_loc


            DO lm_loc = 0,LM_Length - 1


                jloc = d*L_SHIFT + NUM_CFA_VARS * lm_loc


                DO pd = 1,NUM_P_QUAD_POINTS
                    DO td = 1,NUM_T_QUAD_POINTS

                        qd = ((pd-1) * NUM_T_QUAD_POINTS + td-1 ) * NUM_R_QUAD_POINTS

                        Integral_Term(1, qd+1:qd+NUM_R_QUAD_POINTS)                                          &
                                            = Current_Lag_Polys(:,0) * Int_Factor(:, td, pd)                &
                                            * A_Factors(1, td, pd, lpmp_loc)*CC_Ylm_Term(td, pd, lm_loc)

                        Integral_Term(2, qd+1:qd+NUM_R_QUAD_POINTS)                                          &
                                            = Current_Lag_Polys(:,1) * Int_Factor(:, td, pd)                &
                                            * A_Factors(1, td, pd, lpmp_loc)*CC_Ylm_Term(td, pd, lm_loc)

                        Integral_Term(3, qd+1:qd+NUM_R_QUAD_POINTS)                                          &
                                            = Current_Lag_Polys(:,0) * Int_Factor(:, td, pd)                &
                                            * A_Factors(2, td, pd, lpmp_loc)*CC_Ylm_Term(td, pd, lm_loc)

                        Integral_Term(4, qd+1:qd+NUM_R_QUAD_POINTS)                                          &
                                            = Current_Lag_Polys(:,0) * Int_Factor(:, td, pd)                &
                                            * A_Factors(3, td, pd, lpmp_loc)*CC_Ylm_Term(td, pd, lm_loc)

                        Integral_Term(5, qd+1:qd+NUM_R_QUAD_POINTS)                                          &
                                            = Current_Lag_Polys(:,2) * Int_Factor(:, td, pd)                &
                                            * A_Factors(1, td, pd, lpmp_loc)*CC_Ylm_Term(td, pd, lm_loc)

                        Integral_Term(6, qd+1:qd+NUM_R_QUAD_POINTS)                                          &
                                            = Current_Lag_Polys(:,1) * Int_Factor(:, td, pd)                &
                                            * A_Factors(2, td, pd, lpmp_loc)*CC_Ylm_Term(td, pd, lm_loc)

                        Integral_Term(7, qd+1:qd+NUM_R_QUAD_POINTS)                                          &
                                            = Current_Lag_Polys(:,1) * Int_Factor(:, td, pd)                &
                                            * A_Factors(3, td, pd, lpmp_loc)*CC_Ylm_Term(td, pd, lm_loc)

                        Integral_Term(8, qd+1:qd+NUM_R_QUAD_POINTS)                                          &
                                            = Current_Lag_Polys(:,0) * Int_Factor(:, td, pd)                &
                                            * A_Factors(4, td, pd, lpmp_loc)*CC_Ylm_Term(td, pd, lm_loc)

                        Integral_Term(9, qd+1:qd+NUM_R_QUAD_POINTS)                                          &
                                            = Current_Lag_Polys(:,0) * Int_Factor(:, td, pd)                &
                                            * A_Factors(5, td, pd, lpmp_loc)*CC_Ylm_Term(td, pd, lm_loc)

                        Integral_Term(10, qd+1:qd+NUM_R_QUAD_POINTS)                                         &
                                            = Current_Lag_Polys(:,0) * Int_Factor(:, td, pd)                &
                                            * A_Factors(6, td, pd, lpmp_loc)*CC_Ylm_Term(td, pd, lm_loc)

                    END DO
                END DO


                Uncommon_Terms = 0.0_idp

                DO i = 1,NUM_QUAD_DOF

                    Uncommon_Terms(1) = Uncommon_Terms(1)                                                   &
                                      + SUBJCBN_PSI_TERMS_B( 1 , i ) * Integral_Term( 1 , i )

                    Uncommon_Terms(2) = Uncommon_Terms(2)                                                   &
                                      + SUBJCBN_ALPHAPSI_TERMS_B( 1 , i ) * Integral_Term( 1 , i )


                    Uncommon_Terms(3) = Uncommon_Terms(3)                                                   &
                                      + SUBJCBN_BETA1_TERMS_B( 1 , i ) * Integral_Term( 1 , i )     &
                                      + SUBJCBN_BETA1_TERMS_B( 2 , i ) * Integral_Term( 2 , i )     &
                                      + SUBJCBN_BETA1_TERMS_B( 3 , i ) * Integral_Term( 3 , i )     &
                                      + SUBJCBN_BETA1_TERMS_B( 4 , i ) * Integral_Term( 4 , i )


                    Uncommon_Terms(4) = Uncommon_Terms(4)                                                   &
                                      + SUBJCBN_BETA2_TERMS_B( 1 , i ) * Integral_Term( 1 , i )     &
                                      + SUBJCBN_BETA2_TERMS_B( 2 , i ) * Integral_Term( 2 , i )     &
                                      + SUBJCBN_BETA2_TERMS_B( 3 , i ) * Integral_Term( 3 , i )     &
                                      + SUBJCBN_BETA2_TERMS_B( 4 , i ) * Integral_Term( 4 , i )

                    Uncommon_Terms(5) = Uncommon_Terms(5)                                                   &
                                      + SUBJCBN_BETA3_TERMS_B( 1 , i ) * Integral_Term( 1 , i )     &
                                      + SUBJCBN_BETA3_TERMS_B( 2 , i ) * Integral_Term( 2 , i )     &
                                      + SUBJCBN_BETA3_TERMS_B( 3 , i ) * Integral_Term( 3 , i )     &
                                      + SUBJCBN_BETA3_TERMS_B( 4 , i ) * Integral_Term( 4 , i )









                    Uncommon_Terms(6) = Uncommon_Terms(6)                                                   &
                                      + SUBJCBN_PSI_TERMS_B( 2 , i ) * Integral_Term( 1 , i )

                    Uncommon_Terms(7) = Uncommon_Terms(7)                                                   &
                                      + SUBJCBN_ALPHAPSI_TERMS_B( 2 , i ) * Integral_Term( 1 , i )

                    Uncommon_Terms(8) = Uncommon_Terms(8)                                                   &
                                      + SUBJCBN_BETA1_TERMS_B( 5 , i ) * Integral_Term( 1 , i )      &
                                      + SUBJCBN_BETA1_TERMS_B( 6 , i ) * Integral_Term( 2 , i )      &
                                      + SUBJCBN_BETA1_TERMS_B( 7 , i ) * Integral_Term( 3 , i )      &
                                      + SUBJCBN_BETA1_TERMS_B( 8 , i ) * Integral_Term( 4 , i )


                    Uncommon_Terms(9) = Uncommon_Terms(9)                                                   &
                                      + SUBJCBN_BETA2_TERMS_B( 5 , i ) * Integral_Term( 1 , i )      &
                                      + SUBJCBN_BETA2_TERMS_B( 6 , i ) * Integral_Term( 2 , i )      &
                                      + SUBJCBN_BETA2_TERMS_B( 7 , i ) * Integral_Term( 3 , i )      &
                                      + SUBJCBN_BETA2_TERMS_B( 8 , i ) * Integral_Term( 4 , i )

                    Uncommon_Terms(10) = Uncommon_Terms(10)                                                  &
                                       + SUBJCBN_BETA3_TERMS_B( 5 , i ) * Integral_Term( 1 , i )     &
                                       + SUBJCBN_BETA3_TERMS_B( 6 , i ) * Integral_Term( 2 , i )     &
                                       + SUBJCBN_BETA3_TERMS_B( 7 , i ) * Integral_Term( 3 , i )     &
                                       + SUBJCBN_BETA3_TERMS_B( 8 , i ) * Integral_Term( 4 , i )







                    Uncommon_Terms(11) = Uncommon_Terms(11)                                                  &
                                       + SUBJCBN_PSI_TERMS_B( 3 , i ) * Integral_Term( 1 , i )       &
                                       + SUBJCBN_PSI_TERMS_B( 4 , i ) * Integral_Term( 2 , i )       &
                                       + SUBJCBN_PSI_TERMS_B( 5 , i ) * Integral_Term( 3 , i )       &
                                       + SUBJCBN_PSI_TERMS_B( 6 , i ) * Integral_Term( 4 , i )

                    Uncommon_Terms(12) = Uncommon_Terms(12)                                                  &
                                       + SUBJCBN_ALPHAPSI_TERMS_B( 3 , i ) * Integral_Term( 1 , i )  &
                                       + SUBJCBN_ALPHAPSI_TERMS_B( 4 , i ) * Integral_Term( 2 , i )  &
                                       + SUBJCBN_ALPHAPSI_TERMS_B( 5 , i ) * Integral_Term( 3 , i )  &
                                       + SUBJCBN_ALPHAPSI_TERMS_B( 6 , i ) * Integral_Term( 4 , i )


                    Uncommon_Terms(13) = Uncommon_Terms(13)                                                  &
                                       + SUBJCBN_BETA1_TERMS_B( 9 , i ) * Integral_Term( 1 , i )     &
                                       + SUBJCBN_BETA1_TERMS_B( 10 , i ) * Integral_Term( 2 , i )     &
                                       + SUBJCBN_BETA1_TERMS_B( 11 , i ) * Integral_Term( 3 , i )     &
                                       + SUBJCBN_BETA1_TERMS_B( 12 , i ) * Integral_Term( 4 , i )     &
                                       + OneThird * Integral_Term( 5 , i )


                    Uncommon_Terms(14) = Uncommon_Terms(14)                                                  &
                                       + SUBJCBN_BETA2_TERMS_B( 9 , i ) * Integral_Term( 1 , i )     &
                                       + SUBJCBN_BETA2_TERMS_B( 10 , i ) * Integral_Term( 2 , i )     &
                                       + SUBJCBN_BETA2_TERMS_B( 11 , i ) * Integral_Term( 3 , i )     &
                                       + SUBJCBN_BETA2_TERMS_B( 12 , i ) * Integral_Term( 4 , i )     &
                                       + OneThird * Integral_Term( 6 , i )

                    Uncommon_Terms(15) = Uncommon_Terms(15)                                                  &
                                       + SUBJCBN_BETA3_TERMS_B( 9 , i ) * Integral_Term( 1 , i )     &
                                       + SUBJCBN_BETA3_TERMS_B( 11 , i ) * Integral_Term( 3 , i )     &
                                       + SUBJCBN_BETA3_TERMS_B( 12 , i ) * Integral_Term( 4 , i )     &
                                       + OneThird * Integral_Term( 7 , i )






                    Uncommon_Terms(16) = Uncommon_Terms(16)                                                  &
                                       + SUBJCBN_PSI_TERMS_B( 7 , i ) * Integral_Term( 1 , i )       &
                                       + SUBJCBN_PSI_TERMS_B( 8 , i ) * Integral_Term( 2 , i )       &
                                       + SUBJCBN_PSI_TERMS_B( 9 , i ) * Integral_Term( 3 , i )       &
                                       + SUBJCBN_PSI_TERMS_B( 10 , i ) * Integral_Term( 4 , i )

                    Uncommon_Terms(17) = Uncommon_Terms(17)                                                  &
                                       + SUBJCBN_ALPHAPSI_TERMS_B( 7 , i ) * Integral_Term( 1 , i )  &
                                       + SUBJCBN_ALPHAPSI_TERMS_B( 8 , i ) * Integral_Term( 2 , i )  &
                                       + SUBJCBN_ALPHAPSI_TERMS_B( 9 , i ) * Integral_Term( 3 , i )  &
                                       + SUBJCBN_ALPHAPSI_TERMS_B( 10 , i ) * Integral_Term( 4 , i )


                    Uncommon_Terms(18) = Uncommon_Terms(18)                                                  &
                                       + SUBJCBN_BETA1_TERMS_B( 13 , i ) * Integral_Term( 1 , i )     &
                                       + SUBJCBN_BETA1_TERMS_B( 14 , i ) * Integral_Term( 2 , i )     &
                                       + SUBJCBN_BETA1_TERMS_B( 16 , i ) * Integral_Term( 4 , i )     &
                                       + OneThird * Integral_Term( 6 , i )


                    Uncommon_Terms(19) = Uncommon_Terms(19)                                                  &
                                       + SUBJCBN_BETA2_TERMS_B( 13 , i ) * Integral_Term( 1 , i )     &
                                       + SUBJCBN_BETA2_TERMS_B( 14 , i ) * Integral_Term( 2 , i )     &
                                       + SUBJCBN_BETA2_TERMS_B( 15 , i ) * Integral_Term( 3 , i )     &
                                       + SUBJCBN_BETA2_TERMS_B( 16 , i ) * Integral_Term( 4 , i )     &
                                       + OneThird * Integral_Term( 8 , i )

                    Uncommon_Terms(20) = Uncommon_Terms(20)                                                  &
                                       + SUBJCBN_BETA3_TERMS_B( 13 , i ) * Integral_Term( 1 , i )     &
                                       + SUBJCBN_BETA3_TERMS_B( 14 , i ) * Integral_Term( 2 , i )     &
                                       + SUBJCBN_BETA3_TERMS_B( 15 , i ) * Integral_Term( 3 , i )     &
                                       + OneThird * Integral_Term( 9 , i )








                    Uncommon_Terms(21) = Uncommon_Terms(21)                                                  &
                                       + SUBJCBN_PSI_TERMS_B( 11, i ) * Integral_Term( 1 , i )      &
                                       + SUBJCBN_PSI_TERMS_B( 12, i ) * Integral_Term( 2 , i )      &
                                       + SUBJCBN_PSI_TERMS_B( 13, i ) * Integral_Term( 3 , i )      &
                                       + SUBJCBN_PSI_TERMS_B( 14, i ) * Integral_Term( 4 , i )

                    Uncommon_Terms(22) = Uncommon_Terms(22)                                                  &
                                       + SUBJCBN_ALPHAPSI_TERMS_B( 11 , i ) * Integral_Term( 1 , i ) &
                                       + SUBJCBN_ALPHAPSI_TERMS_B( 12 , i ) * Integral_Term( 2 , i ) &
                                       + SUBJCBN_ALPHAPSI_TERMS_B( 13 , i ) * Integral_Term( 3 , i ) &
                                       + SUBJCBN_ALPHAPSI_TERMS_B( 14 , i ) * Integral_Term( 4 , i )


                    Uncommon_Terms(23) = Uncommon_Terms(23)                                                  &
                                       + SUBJCBN_BETA1_TERMS_B( 17 , i ) * Integral_Term( 1 , i )    &
                                       + SUBJCBN_BETA1_TERMS_B( 19 , i ) * Integral_Term( 3 , i )    &
                                       + SUBJCBN_BETA1_TERMS_B( 20 , i ) * Integral_Term( 4 , i )    &
                                       + OneThird * Integral_Term( 7 , i )


                    Uncommon_Terms(24) = Uncommon_Terms(24)                                                  &
                                       + SUBJCBN_BETA2_TERMS_B( 17 , i ) * Integral_Term( 1 , i )    &
                                       + SUBJCBN_BETA2_TERMS_B( 18 , i ) * Integral_Term( 2 , i )    &
                                       + SUBJCBN_BETA2_TERMS_B( 19 , i ) * Integral_Term( 3 , i )    &
                                       + OneThird * Integral_Term( 9 , i )

                    Uncommon_Terms(25) = Uncommon_Terms(25)                                                  &
                                       + SUBJCBN_BETA3_TERMS_B( 17 , i ) * Integral_Term( 1 , i )    &
                                       + SUBJCBN_BETA3_TERMS_B( 18 , i ) * Integral_Term( 2 , i )    &
                                       + SUBJCBN_BETA3_TERMS_B( 19 , i ) * Integral_Term( 3 , i )    &
                                       + SUBJCBN_BETA3_TERMS_B( 20 , i ) * Integral_Term( 4 , i )    &
                                       + OneThird * Integral_Term( 10 , i )

                END DO


                MATVEC_LOC = jloc*ELEM_PROB_DIM + iloc
                BLOCK_ELEM_STF_MATVEC(MATVEC_LOC:MATVEC_LOC+4,Local_re)                             &
                        = BLOCK_ELEM_STF_MATVEC(MATVEC_LOC:MATVEC_LOC+4,Local_re)                   &
                        + Uncommon_Terms(1:5)

                MATVEC_LOC = (jloc+1)*ELEM_PROB_DIM + iloc
                BLOCK_ELEM_STF_MATVEC(MATVEC_LOC:MATVEC_LOC+4,Local_re)                             &
                        = BLOCK_ELEM_STF_MATVEC(MATVEC_LOC:MATVEC_LOC+4,Local_re)                   &
                        + Uncommon_Terms(6:10)

                MATVEC_LOC = (jloc+2)*ELEM_PROB_DIM + iloc
                BLOCK_ELEM_STF_MATVEC(MATVEC_LOC:MATVEC_LOC+4,Local_re)                             &
                        = BLOCK_ELEM_STF_MATVEC(MATVEC_LOC:MATVEC_LOC+4,Local_re)                   &
                        + Uncommon_Terms(11:15)

                MATVEC_LOC = (jloc+3)*ELEM_PROB_DIM + iloc
                BLOCK_ELEM_STF_MATVEC(MATVEC_LOC:MATVEC_LOC+4,Local_re)                             &
                        = BLOCK_ELEM_STF_MATVEC(MATVEC_LOC:MATVEC_LOC+4,Local_re)                   &
                        + Uncommon_Terms(16:20)

                MATVEC_LOC = (jloc+4)*ELEM_PROB_DIM + iloc
                BLOCK_ELEM_STF_MATVEC(MATVEC_LOC:MATVEC_LOC+4,Local_re)                             &
                        = BLOCK_ELEM_STF_MATVEC(MATVEC_LOC:MATVEC_LOC+4,Local_re)                   &
                        + Uncommon_Terms(21:25)



            END DO  ! lm_loc Loop
        END DO  ! lpmp_loc Loop
    END DO  ! dp Loop
END DO  ! d Loop
!$OMP END DO

!$OMP END PARALLEL




END SUBROUTINE CREATE_3D_JCBN_MATRIX_G


