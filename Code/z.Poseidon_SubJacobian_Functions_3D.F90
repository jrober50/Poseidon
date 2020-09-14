   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE SubJacobian_Functions_Module_3D                                              !##!
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

USE Units_Module, &
                                ONLY :  GR_Source_Scalar

USE Poseidon_Constants_Module, &
                                ONLY :  idp,                &
                                        pi,                 &
                                        TwoPi,              &
                                        OneThird,           &
                                        TwoThirds,          &
                                        FourThirds,         &
                                        OneThirtySecond

USE Poseidon_Parameters, &
                                ONLY :  DOMAIN_DIM,                 &
                                        DEGREE,                     &
                                        L_LIMIT,                    &
                                        NUM_SHELLS,                 &
                                        NUM_SUBSHELLS,              &
                                        NUM_SUBSHELLS_PER_SHELL,    &
                                        NUM_CFA_VARS,               &
                                        nPROCS_POSEIDON,            &
                                        NUM_R_QUAD_POINTS,          &
                                        NUM_T_QUAD_POINTS,          &
                                        NUM_P_QUAD_POINTS,          &
                                        NUM_QUAD_DOF,               &
                                        NUM_R_ELEMS_PER_BLOCK,      &
                                        NUM_T_ELEMS_PER_BLOCK,      &
                                        NUM_P_ELEMS_PER_BLOCK,      &
                                        NUM_R_ELEMS_PER_SUBSHELL,   &
                                        NUM_BLOCK_THETA_ROWS,       &
                                        NUM_BLOCK_PHI_COLUMNS,      &
                                        NUM_R_ELEMS_PER_SHELL

USE Poseidon_Variables_Module, &
                                ONLY :  NUM_R_ELEMENTS,             &
                                        NUM_T_ELEMENTS,             &
                                        NUM_P_ELEMENTS,             &
                                        NUM_TP_QUAD_POINTS,         &
                                        rlocs, tlocs, plocs,        &
                                        NUM_R_NODES,                &
                                        INT_R_LOCATIONS,            &
                                        INT_T_LOCATIONS,            &
                                        INT_P_LOCATIONS,            &
                                        INT_R_WEIGHTS,              &
                                        INT_T_WEIGHTS,              &
                                        INT_P_WEIGHTS,              &
                                        VAR_DIM,                    &
                                        Block_Source_E,             &
                                        Block_Source_S,             &
                                        Block_Source_Si,            &
                                        Coefficient_Vector,         &
                                        nPROCS_SHELL,               &
                                        myID_Poseidon,              &
                                        myID_Shell,                 &
                                        myID_SubShell,              &
                                        myID_PETSc,                 &
                                        INNER_CFA_BC_VALUES,        &
                                        OUTER_CFA_BC_VALUES,        &
                                        INNER_CFA_BC_TYPE,          &
                                        OUTER_CFA_BC_TYPE,          &
                                        PROB_DIM,                   &
                                        Block_PROB_DIM,             &
                                        Coefficient_Vector,         &
                                        NUM_OFF_DIAGONALS,          &
                                        LM_LENGTH,                  &
                                        M_VALUES,                   &
                                        ULM_LENGTH,                 &
                                        Local_Length,               &
                                        POSEIDON_COMM_WORLD,        &
                                        POSEIDON_COMM_SHELL,        &
                                        POSEIDON_COMM_PETSC,        &
                                        Block_STF_MAT,              &
                                        ELEM_PROB_DIM,              &
                                        ELEM_PROB_DIM_SQR,          &
                                        SUBSHELL_PROB_DIM,          &
                                        BLOCK_ELEM_STF_MATVEC,      &
                                        myShell,                    &
                                        Block_RHS_Vector,           &
                                        Matrix_Location,            &
                                        LM_Location

CONTAINS





!+101+###########################################################################!
!                                                                                !
!           Calc_EQ1_SubJacobian                                                 !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_EQ1_SubJacobian( SubJacobian_EQ1_Term,                              &
                                 re, te, pe,                                        &
                                 td, pd, tpd, rd,                                   &
                                 CUR_R_LOCS, R_SQUARE, RSIN_SQUARE, COTAN_VAL,      &
                                 SIN_SQUARE, CSC_SQUARE,                            &
                                 PSI_POWER, ALPHAPSI_POWER,                         &
                                 CUR_VAL_BETA, CUR_DRV_BETA,                        &
                                 JCBN_BIGK_VALUE                                    )

REAL(KIND = idp), INTENT(INOUT), DIMENSION( 1:NUM_TP_QUAD_POINTS,   &
                                            1:NUM_R_QUAD_POINTS,    &
                                            1:14                    )   ::  SubJacobian_EQ1_Term

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


REAL(KIND = idp), INTENT(IN)                                            ::  JCBN_BIGK_VALUE





REAL(KIND = idp)                                                        ::  REUSED_VALUE


!PRINT*,"Calc_EQ1_SubJacobian Has been altered"

REUSED_VALUE = PSI_POWER(7)/(16.0_idp* ALPHAPSI_POWER(2) )



! d F_1 / d u_1

SubJacobian_EQ1_Term( tpd, rd, 1)= 10.0_idp*pi                                                  &
                                        * GR_Source_Scalar                                      &
                                        * PSI_POWER(4)                                          &
                                        * Block_Source_E(rd, td, pd, re, te, pe)                &
                                 + (7.0_idp/16.0_idp)                                           &
                                        * PSI_POWER(6)/ALPHAPSI_POWER(2)                        &
                                        * JCBN_BIGK_VALUE

! d F_1 / d u_2
SubJacobian_EQ1_Term( tpd, rd, 2) = -PSI_POWER(7)/(8.0_idp*ALPHAPSI_POWER(3))                   &
                                        * JCBN_BIGK_VALUE





! d F_1 / d u_3
! Reused_Value * kappa_{10}
SubJacobian_EQ1_Term( tpd, rd, 3) = REUSED_VALUE * FourThirds                               &
                        * ( 2.0_idp / R_SQUARE(rd)           * CUR_VAL_BETA(tpd, rd, 1 )      &
                        + COTAN_VAL(td) /CUR_R_LOCS(rd)      * CUR_VAL_BETA(tpd, rd, 2 )      &
                        - 2.0_idp/CUR_R_LOCS(rd)             * CUR_DRV_BETA(tpd, rd, 1,1 )    &
                        + 1.0_idp/CUR_R_LOCS(rd)             * CUR_DRV_BETA(tpd, rd, 2,2 )    &
                        + 1.0_idp/CUR_R_LOCS(rd)             * CUR_DRV_BETA(tpd, rd, 3,3 )    )
!  F_1 / d u_3
! Reused_Value * kappa_{11}
SubJacobian_EQ1_Term( tpd, rd, 4) = REUSED_VALUE * FourThirds                              &
                        * ( - 2.0_idp /CUR_R_LOCS(rd)       * CUR_VAL_BETA(tpd, rd, 1 )      &
                        - COTAN_VAL(td)                     * CUR_VAL_BETA(tpd, rd, 2 )      &
                        + 2.0_idp                           * CUR_DRV_BETA(tpd, rd, 1,1 )    &
                        - 1.0_idp                           * CUR_DRV_BETA(tpd, rd, 2,2 )    &
                        - 1.0_idp                           * CUR_DRV_BETA(tpd, rd, 3,3 )    )
!  F_1 / d u_3
! Reused_Value * kappa_{12}
SubJacobian_EQ1_Term( tpd, rd, 5) = REUSED_VALUE                                           &
                        * ( 2.0_idp / R_SQUARE(rd)      * CUR_DRV_BETA(tpd, rd, 2,1 )        &
                        + 2.0_idp                       * CUR_DRV_BETA(tpd, rd, 1,2 )    )
!  F_1 / d u_3
! Reused_Value * kappa_{13}
SubJacobian_EQ1_Term( tpd, rd, 6) = REUSED_VALUE                                           &
                        * ( 2.0_idp/RSIN_SQUARE(td, rd)      * CUR_DRV_BETA(tpd, rd, 3,1 )    &
                        + 2.0_idp                           * CUR_DRV_BETA(tpd, rd, 1,3 )    )





!  F_1 / d u_4
! Reused_Value * kappa_{20}
SubJacobian_EQ1_Term( tpd, rd, 7) = REUSED_VALUE * FourThirds * COTAN_VAL(td)                      &
                                    * ( 1.0_idp / CUR_R_LOCS(rd)    * CUR_VAL_BETA(tpd, rd, 1 )      &
                                      + 2.0_idp                     * CUR_VAL_BETA(tpd, rd, 2 )      &
                                      - 1.0_idp                     * CUR_DRV_BETA(tpd, rd, 2,2 )    &
                                      + 2.0_idp                     * CUR_DRV_BETA(tpd, rd, 3,3 )    )
! F_1 / d u_4
! Reused_Value * kappa_{21}
SubJacobian_EQ1_Term( tpd, rd, 8) = REUSED_VALUE                                                   &
                                    * ( 2.0_idp * R_SQUARE(rd)      * CUR_DRV_BETA(tpd, rd, 1,2 )    &
                                        + 2.0_idp                   * CUR_DRV_BETA(tpd, rd, 2,1 )    )
! F_1 / d u_4
! Reused_Value * kappa_{22}
SubJacobian_EQ1_Term( tpd, rd, 9) = REUSED_VALUE * FourThirds                                      &
                                * ( 1.0_idp /CUR_R_LOCS(rd)         * CUR_VAL_BETA(tpd, rd, 1 )      &
                                    - COTAN_VAL(td)                 * CUR_VAL_BETA(tpd, rd, 2 )      &
                                    - 1.0_idp                       * CUR_DRV_BETA(tpd, rd, 1,1 )    &
                                    + 2.0_idp                       * CUR_DRV_BETA(tpd, rd, 2,2 )    &
                                    - 1.0_idp                       * CUR_DRV_BETA(tpd, rd, 3,3 )    )
! F_1 / d u_4
! Reused_Value * kappa_{23}
SubJacobian_EQ1_Term( tpd, rd, 10) = REUSED_VALUE                                                   &
                                    * ( 2.0_idp * CSC_SQUARE(td)    * CUR_DRV_BETA(tpd, rd, 3,2 )    &
                                        + 2.0_idp                   * CUR_DRV_BETA(tpd, rd, 2,3 )    )



! F_1 / d u_5
! Reused_Value * kappa_{30}
SubJacobian_EQ1_Term( tpd, rd, 11) = 0.0_idp

! F_1 / d u_5
! Reused_Value * kappa_{31}
SubJacobian_EQ1_Term( tpd, rd, 12) = REUSED_VALUE                                           &
                        * ( 2.0_idp * RSIN_SQUARE(td, rd)   * CUR_DRV_BETA(tpd, rd, 1,3 )    &
                            + 2.0_idp                       * CUR_DRV_BETA(tpd, rd, 3,1 )    )

! F_1 / d u_5
! Reused_Value * kappa_{32}
SubJacobian_EQ1_Term( tpd, rd, 13) = REUSED_VALUE                                           &
                        * ( 2.0_idp * SIN_SQUARE(td)    * CUR_DRV_BETA(tpd, rd, 2,3 )        &
                            + 2.0_idp                   * CUR_DRV_BETA(tpd, rd, 3,2 )        )

! F_1 / d u_5
! Reused_Value * kappa_{33}
SubJacobian_EQ1_Term( tpd, rd, 14) = REUSED_VALUE * FourThirds                              &
                        * ( 1.0_idp /CUR_R_LOCS(rd)         * CUR_VAL_BETA(tpd, rd, 1 )      &
                            + 2.0_idp * COTAN_VAL(td)       * CUR_VAL_BETA(tpd, rd, 2 )      &
                            - 1.0_idp                       * CUR_DRV_BETA(tpd, rd, 1,1 )    &
                            - 1.0_idp                       * CUR_DRV_BETA(tpd, rd, 2,2 )    &
                            + 2.0_idp                       * CUR_DRV_BETA(tpd, rd, 3,3 )    )





END SUBROUTINE Calc_EQ1_SubJacobian









!+102+###########################################################################!
!                                                                                !
!           Calc_EQ2_SubJacobian                                                 !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_EQ2_SubJacobian( SubJacobian_EQ2_Term,                              &
                                 re, te, pe,                                        &
                                 td, pd, tpd, rd,                                   &
                                 CUR_R_LOCS, R_SQUARE, RSIN_SQUARE, COTAN_VAL,      &
                                 SIN_SQUARE, CSC_SQUARE,                            &
                                 PSI_POWER, ALPHAPSI_POWER,                         &
                                 CUR_VAL_BETA, CUR_DRV_BETA,                        &
                                 JCBN_BIGK_VALUE                                    )


REAL(KIND = idp), INTENT(INOUT), DIMENSION( 1:NUM_TP_QUAD_POINTS,   &
                                            1:NUM_R_QUAD_POINTS,    &
                                            1:14                    )   ::  SubJacobian_EQ2_Term

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


REAL(KIND = idp), INTENT(IN)                                            ::  JCBN_BIGK_VALUE



REAL(KIND = idp)                                                        ::  REUSED_VALUE

REUSED_VALUE = (-7.0_idp* PSI_POWER(6))/(16.0_idp* ALPHAPSI_POWER(1) )



! d F_2 / d u_1
SubJacobian_EQ2_Term(tpd, rd, 1) = -8.0_idp * pi * ALPHAPSI_POWER(1) * PSI_POWER(3)             &
                                    * GR_Source_Scalar                                          &
                                    * ( Block_Source_E(rd, td, pd, re, te, pe)                  &
                                        + 2.0_idp * Block_Source_S(rd, td, pd, re, te, pe) )    &
                                    - (42.0_idp / 16.0_idp )                                    &
                                        * PSI_POWER(5)/ALPHAPSI_POWER(1)                        &
                                        * JCBN_BIGK_VALUE


! d F_2 / d u_2
SubJacobian_EQ2_Term(tpd, rd, 2) = -TwoPi * PSI_POWER(4)                                        &
                                    * GR_Source_Scalar                                          &
                                    *( Block_Source_E(rd, td, pd, re, te, pe)                   &
                                        + 2.0_idp * Block_Source_S(rd, td, pd, re, te, pe) )    &
                                    + (7.0_idp / 16.0_idp)                                      &
                                        * PSI_POWER(6)/ALPHAPSI_POWER(2)                        &
                                        * JCBN_BIGK_VALUE




! d F_2 / d u_3
! Reused_Value * kappa_{10}
SubJacobian_EQ2_Term( tpd, rd, 3) = REUSED_VALUE * FourThirds                               &
                        * ( 2.0_idp / R_SQUARE(rd)           * CUR_VAL_BETA(tpd, rd, 1 )      &
                        + COTAN_VAL(td) /CUR_R_LOCS(rd)      * CUR_VAL_BETA(tpd, rd, 2 )      &
                        - 2.0_idp/CUR_R_LOCS(rd)             * CUR_DRV_BETA(tpd, rd, 1,1 )    &
                        + 1.0_idp/CUR_R_LOCS(rd)             * CUR_DRV_BETA(tpd, rd, 2,2 )    &
                        + 1.0_idp/CUR_R_LOCS(rd)             * CUR_DRV_BETA(tpd, rd, 3,3 )    )
!  F_2 / d u_3
! Reused_Value * kappa_{11}
SubJacobian_EQ2_Term( tpd, rd, 4) = REUSED_VALUE * FourThirds                              &
                        * ( - 2.0_idp /CUR_R_LOCS(rd)       * CUR_VAL_BETA(tpd, rd, 1 )      &
                        - COTAN_VAL(td)                     * CUR_VAL_BETA(tpd, rd, 2 )      &
                        + 2.0_idp                           * CUR_DRV_BETA(tpd, rd, 1,1 )    &
                        - 1.0_idp                           * CUR_DRV_BETA(tpd, rd, 2,2 )    &
                        - 1.0_idp                           * CUR_DRV_BETA(tpd, rd, 3,3 )    )
!  F_2 / d u_3
! Reused_Value * kappa_{12}
SubJacobian_EQ2_Term( tpd, rd, 5) = REUSED_VALUE                                           &
                        * ( 2.0_idp / R_SQUARE(rd)      * CUR_DRV_BETA(tpd, rd, 2,1 )        &
                        + 2.0_idp                       * CUR_DRV_BETA(tpd, rd, 1,2 )    )
!  F_2 / d u_3
! Reused_Value * kappa_{13}
SubJacobian_EQ2_Term( tpd, rd, 6) = REUSED_VALUE                                           &
                        * ( 2.0_idp/RSIN_SQUARE(td, rd)      * CUR_DRV_BETA(tpd, rd, 3,1 )    &
                        + 2.0_idp                           * CUR_DRV_BETA(tpd, rd, 1,3 )    )





!  F_2 / d u_4
! Reused_Value * kappa_{20}
SubJacobian_EQ2_Term( tpd, rd, 7) = REUSED_VALUE * FourThirds * COTAN_VAL(td)                      &
                                    * ( 1.0_idp / CUR_R_LOCS(rd)    * CUR_VAL_BETA(tpd, rd, 1 )      &
                                      + 2.0_idp                     * CUR_VAL_BETA(tpd, rd, 2 )      &
                                      - 1.0_idp                     * CUR_DRV_BETA(tpd, rd, 2,2 )    &
                                      + 2.0_idp                     * CUR_DRV_BETA(tpd, rd, 3,3 )    )
! F_2 / d u_4
! Reused_Value * kappa_{21}
SubJacobian_EQ2_Term( tpd, rd, 8) = REUSED_VALUE                                                   &
                                    * ( 2.0_idp * R_SQUARE(rd)      * CUR_DRV_BETA(tpd, rd, 1,2 )    &
                                        + 2.0_idp                   * CUR_DRV_BETA(tpd, rd, 2,1 )    )
! F_2 / d u_4
! Reused_Value * kappa_{22}
SubJacobian_EQ2_Term( tpd, rd, 9) = REUSED_VALUE * FourThirds                                      &
                                * ( 1.0_idp /CUR_R_LOCS(rd)         * CUR_VAL_BETA(tpd, rd, 1 )      &
                                    - COTAN_VAL(td)                 * CUR_VAL_BETA(tpd, rd, 2 )      &
                                    - 1.0_idp                       * CUR_DRV_BETA(tpd, rd, 1,1 )    &
                                    + 2.0_idp                       * CUR_DRV_BETA(tpd, rd, 2,2 )    &
                                    - 1.0_idp                       * CUR_DRV_BETA(tpd, rd, 3,3 )    )
! F_2 / d u_4
! Reused_Value * kappa_{23}
SubJacobian_EQ2_Term( tpd, rd, 10) = REUSED_VALUE                                                   &
                                    * ( 2.0_idp * CSC_SQUARE(td)    * CUR_DRV_BETA(tpd, rd, 3,2 )    &
                                        + 2.0_idp                   * CUR_DRV_BETA(tpd, rd, 2,3 )    )



! F_2 / d u_5
! Reused_Value * kappa_{30}
SubJacobian_EQ2_Term( tpd, rd, 11) = 0.0_idp

! F_2 / d u_5
! Reused_Value * kappa_{31}
SubJacobian_EQ2_Term( tpd, rd, 12) = REUSED_VALUE                                           &
                        * ( 2.0_idp * RSIN_SQUARE(td, rd)   * CUR_DRV_BETA(tpd, rd, 1,3 )    &
                            + 2.0_idp                       * CUR_DRV_BETA(tpd, rd, 3,1 )    )

! F_2 / d u_5
! Reused_Value * kappa_{32}
SubJacobian_EQ2_Term( tpd, rd, 13) = REUSED_VALUE                                           &
                        * ( 2.0_idp * SIN_SQUARE(td)    * CUR_DRV_BETA(tpd, rd, 2,3 )        &
                            + 2.0_idp                   * CUR_DRV_BETA(tpd, rd, 3,2 )        )

! F_2 / d u_5
! Reused_Value * kappa_{33}
SubJacobian_EQ2_Term( tpd, rd, 14) = REUSED_VALUE * FourThirds                              &
                        * ( 1.0_idp /CUR_R_LOCS(rd)         * CUR_VAL_BETA(tpd, rd, 1 )      &
                            + 2.0_idp * COTAN_VAL(td)       * CUR_VAL_BETA(tpd, rd, 2 )      &
                            - 1.0_idp                       * CUR_DRV_BETA(tpd, rd, 1,1 )    &
                            - 1.0_idp                       * CUR_DRV_BETA(tpd, rd, 2,2 )    &
                            + 2.0_idp                       * CUR_DRV_BETA(tpd, rd, 3,3 )    )




END SUBROUTINE Calc_EQ2_SubJacobian








!+103+###########################################################################!
!                                                                                !
!           Calc_EQ3_SubJacobian                                                 !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_EQ3_SubJacobian( SubJacobian_EQ3_Term,                              &
                                 re, te, pe,                                        &
                                 td, pd, tpd, rd,                                   &
                                 CUR_R_LOCS, R_SQUARE, RSIN_SQUARE, COTAN_VAL,      &
                                 SIN_SQUARE, CSC_SQUARE,                            &
                                 PSI_POWER, ALPHAPSI_POWER,                         &
                                 CUR_DRV_PSI, CUR_DRV_ALPHAPSI,                     &
                                 JCBN_BIGK_VALUE, JCBN_n_Array, JCBN_Kappa_Array    )

REAL(KIND = idp), INTENT(INOUT), DIMENSION( 1:NUM_TP_QUAD_POINTS,   &
                                            1:NUM_R_QUAD_POINTS,    &
                                            1:20                    )   ::  SubJacobian_EQ3_Term

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
                                         1:3                    )       ::  CUR_DRV_PSI

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:NUM_TP_QUAD_POINTS,  &
                                         1:NUM_R_QUAD_POINTS,   &
                                         1:3                    )       ::  CUR_DRV_ALPHAPSI



REAL(KIND = idp), INTENT(IN)                                            ::  JCBN_BIGK_VALUE

REAL(KIND = idp), INTENT(IN), DIMENSION(1:3)                            ::  JCBN_n_ARRAY
REAL(KIND = idp), INTENT(IN), DIMENSION(1:3,1:3)                        ::  JCBN_kappa_Array


INTEGER                                                                 ::  testa, testb



!! d F_3 / d u_1 Non-Derivative Term
!SubJacobian_EQ3_Term( tpd, rd, 1) = (-7.0_idp / PSI_POWER(2) )                                       &
!                                        * ( CUR_DRV_PSI(tpd, rd, 1 )*JCBN_kappa_ARRAY(1,1)       &
!                                          + CUR_DRV_PSI(tpd, rd, 2 )*JCBN_kappa_ARRAY(2,1)       &
!                                          + CUR_DRV_PSI(tpd, rd, 3 )*JCBN_kappa_ARRAY(3,1) )     &
!                                    - 48.0_idp * pi                                                 &
!                                        * GR_Source_Scalar                                          &
!                                        * Block_Source_Si(rd, td, pd, re, te, pe, 1)                &
!                                        * ALPHAPSI_POWER(1)                                         &
!                                        * PSI_POWER(2)
!
!
!! d F_3 / d u_1 Derivative Terms
!SubJacobian_EQ3_Term( tpd, rd, 2:4)= (7.0_idp/ PSI_POWER(1) )*JCBN_kappa_ARRAY(1:3,1)

! d F_3 / d u_1 Non-Derivative Term
SubJacobian_EQ3_Term( tpd, rd, 1) = (- 1.0_idp / PSI_POWER(2) )                                     &
                                        * ( CUR_DRV_PSI(tpd, rd, 1 )*JCBN_kappa_ARRAY(1,1)       &
                                          + CUR_DRV_PSI(tpd, rd, 2 )*JCBN_kappa_ARRAY(2,1)       &
                                          + CUR_DRV_PSI(tpd, rd, 3 )*JCBN_kappa_ARRAY(3,1) )     &
                                    - 48.0_idp * pi                                                 &
                                        * GR_Source_Scalar                                          &
                                        * Block_Source_Si(rd, td, pd, re, te, pe, 1)                &
                                        * ALPHAPSI_POWER(1)                                         &
                                        * PSI_POWER(2)


! d F_3 / d u_1 Derivative Terms
SubJacobian_EQ3_Term( tpd, rd, 2:4)= JCBN_kappa_ARRAY(1:3,1)/ PSI_POWER(1)




!! d F_3 / d u_2 Non-Derivative Term
!SubJacobian_EQ3_Term( tpd, rd, 5) = ( CUR_DRV_ALPHAPSI(tpd, rd, 1 )*JCBN_kappa_Array(1,1)       &
!                                        + CUR_DRV_ALPHAPSI(tpd, rd, 2 )*JCBN_kappa_Array(2,1)       &
!                                        + CUR_DRV_ALPHAPSI(tpd, rd, 3 )*JCBN_kappa_Array(3,1)   )   &
!                                        /ALPHAPSI_POWER(2)                                          &
!                                        -16.0_idp * pi                                              &
!                                            * GR_Source_Scalar                                      &
!                                            * Block_Source_Si(rd,td,pd,re,te,pe,1)                  &
!                                            * PSI_POWER(3)
!
!
!! d F_3 / d u_2 Derivative Terms
!SubJacobian_EQ3_Term( tpd, rd, 6:8 ) = (-1.0_idp/ALPHAPSI_POWER(1) )*JCBN_kappa_ARRAY(1:3,1)


! d F_3 / d u_2 Non-Derivative Term
SubJacobian_EQ3_Term( tpd, rd, 5) = (- 2.0_idp /ALPHAPSI_POWER(2))                                  &
                                      * ( CUR_DRV_ALPHAPSI(tpd, rd, 1 )*JCBN_kappa_Array(1,1)       &
                                        + CUR_DRV_ALPHAPSI(tpd, rd, 2 )*JCBN_kappa_Array(2,1)       &
                                        + CUR_DRV_ALPHAPSI(tpd, rd, 3 )*JCBN_kappa_Array(3,1)   )   &
                                        -16.0_idp * pi                                              &
                                            * GR_Source_Scalar                                      &
                                            * Block_Source_Si(rd,td,pd,re,te,pe,1)                  &
                                            * PSI_POWER(3)


! d F_3 / d u_2 Derivative Terms
SubJacobian_EQ3_Term( tpd, rd, 6:8 ) = (1.0_idp/ALPHAPSI_POWER(1) )*JCBN_kappa_ARRAY(1:3,1)




!  F_3 / d u_3 Non-Derivative Term
SubJacobian_EQ3_Term( tpd, rd, 9) = (4.0_idp/(3.0_idp *CUR_R_LOCS(rd))) * JCBN_n_ARRAY(1)      &
                                    - ( 8.0_idp/(3.0_idp * R_SQUARE(rd) ) )

!  F_3 / d u_3 Derivative Terms
SubJacobian_EQ3_Term( tpd, rd, 10) = OneThird * ((2.0_idp/CUR_R_LOCS(rd)) - 4.0_idp * JCBN_n_ARRAY(1) )
SubJacobian_EQ3_Term( tpd, rd, 11) = -(1.0_idp/R_SQUARE(rd)) * JCBN_n_ARRAY(2)
SubJacobian_EQ3_Term( tpd, rd, 12) = -(1.0_idp/RSIN_SQUARE(td, rd)) * JCBN_n_ARRAY(3)

! F_3 / d u_4 Non-Derivative Term
SubJacobian_EQ3_Term( tpd, rd, 13) = -2.0_idp * COTAN_VAL(td) / CUR_R_LOCS(rd)                      &
                                     + TwoThirds*COTAN_VAL(td) * JCBN_n_ARRAY(1)

! F_3 / d u_4 Derivative Terms
SubJacobian_EQ3_Term( tpd, rd, 14) = OneThird * COTAN_VAL(td) - JCBN_n_ARRAY(2)
SubJacobian_EQ3_Term( tpd, rd, 15) = -2.0_idp/CUR_R_LOCS(rd) + TwoThirds * JCBN_n_ARRAY(1)
SubJacobian_EQ3_Term( tpd, rd, 16) = 0.0_idp


! F_3 / d u_5 Non-Derivative Term
SubJacobian_EQ3_Term( tpd, rd, 17) =  0.0_idp

! F_3 / d u_5 Derivative Terms
SubJacobian_EQ3_Term( tpd, rd, 18) = -JCBN_n_ARRAY(3)
SubJacobian_EQ3_Term( tpd, rd, 19) = 0.0_idp
SubJacobian_EQ3_Term( tpd, rd, 20) = -2.0_idp/CUR_R_LOCS(rd) - TwoThirds * JCBN_n_ARRAY(1)


!testa = 1
!testb = 1
!SubJacobian_EQ3_Term( tpd, rd, 1:testa)  = 0.0_idp
!SubJacobian_EQ3_Term( tpd, rd, testb:20) = 0.0_idp


END SUBROUTINE Calc_EQ3_SubJacobian








!+104+###########################################################################!
!                                                                                !
!           Calc_EQ4_SubJacobian                                                 !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_EQ4_SubJacobian( SubJacobian_EQ4_Term,                              &
                                 re, te, pe,                                        &
                                 td, pd, tpd, rd,                                   &
                                 CUR_R_LOCS, R_SQUARE, R_CUBED, RSIN_SQUARE,        &
                                 COTAN_VAL, SIN_SQUARE, CSC_SQUARE,                 &
                                 PSI_POWER, ALPHAPSI_POWER,                         &
                                 CUR_DRV_PSI, CUR_DRV_ALPHAPSI,                     &
                                 JCBN_BIGK_VALUE, JCBN_n_Array, JCBN_Kappa_Array    )

REAL(KIND = idp), INTENT(INOUT), DIMENSION( 1:NUM_TP_QUAD_POINTS,   &
                                            1:NUM_R_QUAD_POINTS,    &
                                            1:20                    )   ::  SubJacobian_EQ4_Term

INTEGER, INTENT(IN)                                                     ::  re, te, pe
INTEGER, INTENT(IN)                                                     ::  td, pd, tpd, rd

REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_R_QUAD_POINTS)            ::  CUR_R_LOCS

REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_R_QUAD_POINTS)            ::  R_SQUARE

REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_R_QUAD_POINTS)            ::  R_CUBED

REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_T_QUAD_POINTS,    &
                                        1:NUM_R_QUAD_POINTS     )       ::  RSIN_SQUARE

REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_T_QUAD_POINTS)            ::  SIN_SQUARE,         &
                                                                            CSC_SQUARE,         &
                                                                            COTAN_VAL

REAL(KIND = idp), INTENT(IN), DIMENSION(1:11)                           ::  PSI_POWER
REAL(KIND = idp), INTENT(IN), DIMENSION(1:4)                            ::  ALPHAPSI_POWER

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:NUM_TP_QUAD_POINTS,  &
                                         1:NUM_R_QUAD_POINTS,   &
                                         1:3                    )       ::  CUR_DRV_PSI

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:NUM_TP_QUAD_POINTS,  &
                                         1:NUM_R_QUAD_POINTS,   &
                                         1:3                    )       ::  CUR_DRV_ALPHAPSI



REAL(KIND = idp), INTENT(IN)                                            ::  JCBN_BIGK_VALUE

REAL(KIND = idp), INTENT(IN), DIMENSION(1:3)                            ::  JCBN_n_ARRAY
REAL(KIND = idp), INTENT(IN), DIMENSION(1:3,1:3)                        ::  JCBN_kappa_Array




! d F_4 / d u_1 Non-Derivative Term
SubJacobian_EQ4_Term( tpd, rd, 1) = (-7.0_idp / PSI_POWER(2) )                                       &
                                       * ( CUR_DRV_PSI(tpd, rd, 1 )*JCBN_kappa_ARRAY(1,2)         &
                                         + CUR_DRV_PSI(tpd, rd, 2 )*JCBN_kappa_ARRAY(2,2)         &
                                         + CUR_DRV_PSI(tpd, rd, 3 )*JCBN_kappa_ARRAY(3,2) )       &
                                    - 48.0_idp * pi                                                 &
                                        * GR_Source_Scalar                                          &
                                        * Block_Source_Si(rd, td, pd, re, te, pe, 2)                &
                                        * ALPHAPSI_POWER(1)                                         &
                                        * PSI_POWER(2)

! d F_4 / d u_1 Derivative Terms
SubJacobian_EQ4_Term( tpd, rd, 2:4) = (7.0_idp/ PSI_POWER(1) )*JCBN_kappa_ARRAY(1:3,2)

! d F_4 / d u_2 Non-Derivative Term
SubJacobian_EQ4_Term( tpd, rd, 5) = ( CUR_DRV_ALPHAPSI(tpd, rd, 1 )*JCBN_kappa_Array(1,2)       &
                                        + CUR_DRV_ALPHAPSI(tpd, rd, 2 )*JCBN_kappa_Array(2,2)       &
                                        + CUR_DRV_ALPHAPSI(tpd, rd, 3 )*JCBN_kappa_Array(3,2)   )   &
                                        /ALPHAPSI_POWER(2)                                          &
                                        -16.0_idp * pi                                              &
                                            * GR_Source_Scalar                                      &
                                            * Block_Source_Si(rd, td, pd, re, te, pe, 2)            &
                                            * PSI_POWER(3)


! d F_4 / d u_2 Derivative Terms
SubJacobian_EQ4_Term( tpd, rd, 6:8) = (-1.0_idp/ ALPHAPSI_POWER(1) )*JCBN_kappa_ARRAY(1:3,2)


!  F_4 / d u_3 Non-Derivative Term
SubJacobian_EQ4_Term( tpd, rd, 9) = -(TwoThirds/R_CUBED(rd)) * JCBN_n_ARRAY(2)

!  F_4 / d u_3 Derivative Terms
SubJacobian_EQ4_Term( tpd, rd, 10) = (TwoThirds / R_SQUARE(rd)) * JCBN_n_ARRAY(2) - COTAN_VAL(td)/(3.0_idp * R_SQUARE(rd))
SubJacobian_EQ4_Term( tpd, rd, 11) = (8.0_idp / (3.0_idp * R_CUBED(rd))) - JCBN_n_ARRAY(1)/R_SQUARE(rd)
SubJacobian_EQ4_Term( tpd, rd, 12) = 0.0_idp


! F_4 / d u_4 Non-Derivative Term
SubJacobian_EQ4_Term( tpd, rd, 13) = (TwoThirds/R_SQUARE(rd)) * COTAN_VAL(td) * JCBN_n_ARRAY(2)               &
                                    + 1.0_idp/R_SQUARE(rd) - TwoThirds/RSIN_SQUARE(td,rd)

! F_4 / d u_4 Derivative Terms
SubJacobian_EQ4_Term( tpd, rd, 14) = -JCBN_n_ARRAY(1) + 2.0_idp/CUR_R_LOCS(rd)
SubJacobian_EQ4_Term( tpd, rd, 15) = (OneThird / R_SQUARE(rd) )* ( 4.0_idp * JCBN_n_ARRAY(2) )
SubJacobian_EQ4_Term( tpd, rd, 16) = - JCBN_n_ARRAY(3)/RSIN_SQUARE(td, rd)


! F_4 / d u_5 Non-Derivative Term
SubJacobian_EQ4_Term( tpd, rd, 17) = 0.0_idp

! F_4 / d u_5 Derivative Terms
SubJacobian_EQ4_Term( tpd, rd, 18) = 0.0_idp
SubJacobian_EQ4_Term( tpd, rd, 19) = (-1.0_idp/ R_SQUARE(rd)) * JCBN_n_ARRAY(3)
SubJacobian_EQ4_Term( tpd, rd, 20) = (TwoThirds / R_SQUARE(rd)) * JCBN_n_ARRAY(2)                           &
                                    -  3.5_idp*COTAN_VAL(td)/R_SQUARE(rd)



END SUBROUTINE Calc_EQ4_SubJacobian







!+105+###########################################################################!
!                                                                                !
!           Calc_EQ5_SubJacobian                                                 !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_EQ5_SubJacobian( SubJacobian_EQ5_Term,                              &
                                 re, te, pe,                                        &
                                 td, pd, tpd, rd,                                   &
                                 CUR_R_LOCS, R_SQUARE, R_CUBED, RSIN_SQUARE,        &
                                 COTAN_VAL, SIN_SQUARE, CSC_SQUARE,                 &
                                 PSI_POWER, ALPHAPSI_POWER,                         &
                                 CUR_DRV_PSI, CUR_DRV_ALPHAPSI,                     &
                                 JCBN_BIGK_VALUE, JCBN_n_Array, JCBN_Kappa_Array    )

REAL(KIND = idp), INTENT(INOUT), DIMENSION( 1:NUM_TP_QUAD_POINTS,   &
                                            1:NUM_R_QUAD_POINTS,    &
                                            1:20                    )   ::  SubJacobian_EQ5_Term

INTEGER, INTENT(IN)                                                     ::  re, te, pe
INTEGER, INTENT(IN)                                                     ::  td, pd, tpd, rd

REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_R_QUAD_POINTS)            ::  CUR_R_LOCS

REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_R_QUAD_POINTS)            ::  R_SQUARE

REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_R_QUAD_POINTS)            ::  R_CUBED

REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_T_QUAD_POINTS,    &
                                        1:NUM_R_QUAD_POINTS     )       ::  RSIN_SQUARE

REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_T_QUAD_POINTS)            ::  SIN_SQUARE,         &
                                                                            CSC_SQUARE,         &
                                                                            COTAN_VAL

REAL(KIND = idp), INTENT(IN), DIMENSION(1:11)                           ::  PSI_POWER
REAL(KIND = idp), INTENT(IN), DIMENSION(1:4)                            ::  ALPHAPSI_POWER

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:NUM_TP_QUAD_POINTS,  &
                                         1:NUM_R_QUAD_POINTS,   &
                                         1:3                    )       ::  CUR_DRV_PSI

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:NUM_TP_QUAD_POINTS,  &
                                         1:NUM_R_QUAD_POINTS,   &
                                         1:3                    )       ::  CUR_DRV_ALPHAPSI







REAL(KIND = idp), INTENT(IN)                                            ::  JCBN_BIGK_VALUE

REAL(KIND = idp), INTENT(IN), DIMENSION(1:3)                            ::  JCBN_n_ARRAY
REAL(KIND = idp), INTENT(IN), DIMENSION(1:3,1:3)                        ::  JCBN_kappa_Array



! d F_5 / d u_1 Non-Derivative Term
SubJacobian_EQ5_Term( tpd, rd, 1) = (-7.0_idp / PSI_POWER(2) )                                      &
                                    * ( CUR_DRV_PSI(tpd, rd, 1 )*JCBN_kappa_ARRAY(1,3)            &
                                      + CUR_DRV_PSI(tpd, rd, 2 )*JCBN_kappa_ARRAY(2,3)            &
                                      + CUR_DRV_PSI(tpd, rd, 3 )*JCBN_kappa_ARRAY(3,3) )          &
                                    - 48.0_idp * pi                                                 &
                                        * GR_Source_Scalar                                          &
                                        * Block_Source_Si(rd, td, pd, re, te, pe, 3)                &
                                        * ALPHAPSI_POWER(1)                                         &
                                        * PSI_POWER(2)

! d F_5 / d u_1 Derivative Terms
SubJacobian_EQ5_Term( tpd, rd, 2:4) = (7.0_idp/ PSI_POWER(1) )*JCBN_kappa_ARRAY(1:3,3)





! d F_5 / d u_2 Non-Derivative Term
SubJacobian_EQ5_Term( tpd, rd, 5) = ( CUR_DRV_ALPHAPSI(tpd, rd, 1 ) * JCBN_kappa_Array(1,3)    &
                                         + CUR_DRV_ALPHAPSI(tpd, rd, 2 ) * JCBN_kappa_Array(2,3)    &
                                         + CUR_DRV_ALPHAPSI(tpd, rd, 3 ) * JCBN_kappa_Array(3,3)  ) &
                                        /ALPHAPSI_POWER(2)                                          &
                                        -16.0_idp * pi                                              &
                                            * GR_Source_Scalar                                      &
                                            * Block_Source_Si(rd,td,pd,re,te,pe,3)                  &
                                            * PSI_POWER(3)

! d F_5 / d u_2 Derivative Terms
SubJacobian_EQ5_Term( tpd, rd, 6:8)  = (-1.0_idp/ ALPHAPSI_POWER(1) )*JCBN_kappa_ARRAY(1:3,3)


!  F_5 / d u_3 Non-Derivative Term
SubJacobian_EQ5_Term( tpd, rd, 9) = -(TwoThirds*CSC_SQUARE(td)/R_CUBED(rd))* JCBN_n_ARRAY(3)

! J F_5 / d u_3 Derivative Terms
SubJacobian_EQ5_Term( tpd, rd, 10) = (TwoThirds/RSIN_SQUARE(td, rd)) * JCBN_n_ARRAY(3)
SubJacobian_EQ5_Term( tpd, rd, 11) = 0.0_idp
SubJacobian_EQ5_Term( tpd, rd, 12) = 8.0_idp * OneThird /(RSIN_SQUARE(td, rd)*CUR_R_LOCS(rd) )             &
                                        - JCBN_n_ARRAY(1)/RSIN_SQUARE(td, rd)


! F_5 / d u_4 Non-Derivative Term
SubJacobian_EQ5_Term( tpd, rd, 13) = -FourThirds*(COTAN_VAL(td)/RSIN_SQUARE(td, rd))* JCBN_n_ARRAY(3)

! F_5 / d u_4 Derivative Terms
SubJacobian_EQ5_Term( tpd, rd, 14) = 0.0_idp
SubJacobian_EQ5_Term( tpd, rd, 15) = (TwoThirds / RSIN_SQUARE(td, rd)) * JCBN_n_ARRAY(3)
SubJacobian_EQ5_Term( tpd, rd, 16) = (OneThird / RSIN_SQUARE(td, rd))                                      &
                                    * ( 4.0_idp*COTAN_VAL(td) - 3.0_idp * JCBN_n_ARRAY(2)  )


! F_5 / d u_5 Non-Derivative Term
SubJacobian_EQ5_Term( tpd, rd, 17) = 0.0_idp

! F_5 / d u_5 Derivative Terms
SubJacobian_EQ5_Term( tpd, rd, 18) = -JCBN_n_ARRAY(1) + 2.0_idp/CUR_R_LOCS(rd)
SubJacobian_EQ5_Term( tpd, rd, 19) = 2.0_idp*COTAN_VAL(td)/R_SQUARE(rd) - JCBN_n_ARRAY(2)/R_SQUARE(rd)
SubJacobian_EQ5_Term( tpd, rd, 20) = -(FourThirds / RSIN_SQUARE(td, rd)) * JCBN_n_ARRAY(3)


END SUBROUTINE Calc_EQ5_SubJacobian










!+106+###########################################################################!
!                                                                                !
!           Calc_RHS_Terms                                                       !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_RHS_Terms( RHS_Terms,                                         &
                           re, te, pe,                                        &
                           td, pd, tpd, rd,                                   &
                           CUR_R_LOCS, R_SQUARE, RSIN_SQUARE, COTAN_VAL,      &
                           SIN_SQUARE, CSC_SQUARE,                            &
                           PSI_POWER, ALPHAPSI_POWER,                         &
                           CUR_VAL_BETA, CUR_DRV_BETA,                        &
                           JCBN_BIGK_VALUE, JCBN_n_Array, JCBN_Kappa_Array    )

REAL(KIND = idp), INTENT(INOUT), DIMENSION( 1:NUM_TP_QUAD_POINTS,   &
                                            1:NUM_R_QUAD_POINTS,    &
                                            1:5                     )   ::  RHS_Terms

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



REAL(KIND = idp), INTENT(IN)                                            ::  JCBN_BIGK_VALUE
REAL(KIND = idp), INTENT(IN), DIMENSION(1:3)                            ::  JCBN_n_ARRAY
REAL(KIND = idp), INTENT(IN), DIMENSION(1:3,1:3)                        ::  JCBN_kappa_Array

REAL(KIND = idp)                                                        ::  Beta_Source_Prefix



!PRINT*,"Calc_RHS_Terms has been altered"

Beta_Source_Prefix = 16.0_idp * pi * ALPHAPSI_POWER(1) * PSI_POWER(3) * GR_Source_Scalar


RHS_Terms(tpd, rd, 1) = - TwoPi * GR_Source_Scalar * Block_Source_E(rd, td, pd, re, te, pe) * PSI_POWER(5)   &
                         - PSI_POWER(7)/ (16.0_idp * ALPHAPSI_POWER(2)) * JCBN_BIGK_VALUE


RHS_Terms(tpd, rd, 2) = TwoPi * ALPHAPSI_POWER(1) * PSI_POWER(4)                             &
                * GR_Source_Scalar * ( Block_Source_E(rd, td, pd, re, te, pe)                   &
                                       + 2.0_idp * Block_Source_S(rd, td, pd, re, te, pe)  )    &
                + 7.0_idp*PSI_POWER(6)/ (16.0_idp * ALPHAPSI_POWER(1)) * JCBN_BIGK_VALUE



!RHS_Terms(tpd, rd, 3) = Beta_Source_Prefix * Block_Source_Si(rd, td, pd, re, te, pe, 1)             &
!                      + ( 8.0_idp/(3.0_idp * R_SQUARE(rd) )                                         &
!                            - 4.0_idp* JCBN_n_ARRAY(1)/(3.0_idp * CUR_R_LOCS(rd)) )                 &
!                        * CUR_VAL_BETA(tpd, rd, 1)                                                  &
!                      + ( 2.0_idp * COTAN_VAL(td)/CUR_R_LOCS(rd)                                    &
!                          - TwoThirds * JCBN_n_ARRAY(1)*COTAN_VAL(td)  )                            &
!                        * CUR_VAL_BETA(tpd, rd, 2 )                                                 &
!                      + ( FourThirds * JCBN_n_ARRAY(1) )                                            &
!                        * CUR_DRV_BETA(tpd, rd, 1, 1 )                                              &
!                      + ( JCBN_n_ARRAY(2)/R_SQUARE(rd) )                                            &
!                        * CUR_DRV_BETA(tpd, rd, 2, 1 )                                              &
!                      + ( JCBN_n_ARRAY(3)/RSIN_SQUARE(td,rd) )                                      &
!                        * CUR_DRV_BETA(tpd, rd, 3, 1 )                                              &
!                      + ( JCBN_n_ARRAY(2) - COTAN_VAL(td)/3.0_idp )                                 &
!                        * CUR_DRV_BETA(tpd, rd, 1, 2 )                                              &
!                      + ( 2.0_idp * FourThirds / CUR_R_LOCS(rd) - TwoThirds * JCBN_n_ARRAY(1) )    &
!                        * CUR_DRV_BETA(tpd, rd, 2, 2 )                                              &
!                      + ( JCBN_n_ARRAY(3)  )                                                        &
!                        * CUR_DRV_BETA(tpd, rd, 1, 3 )                                              &
!                      + ( 2.0_idp * FourThirds / CUR_R_LOCS(rd) - TwoThirds * JCBN_n_Array(1) )     &
!                        * CUR_DRV_BETA(tpd, rd, 3, 3 )

RHS_Terms(tpd, rd, 3) = Beta_Source_Prefix * Block_Source_Si(rd, td, pd, re, te, pe, 1)             &
                      + ( 8.0_idp/(3.0_idp * R_SQUARE(rd) )                                         &
                            - 4.0_idp* JCBN_n_ARRAY(1)/(3.0_idp * CUR_R_LOCS(rd)) )                 &
                        * CUR_VAL_BETA(tpd, rd, 1)                                                  &
                      + ( 2.0_idp * COTAN_VAL(td)/CUR_R_LOCS(rd)                                    &
                          - TwoThirds * JCBN_n_ARRAY(1)*COTAN_VAL(td)  )                            &
                        * CUR_VAL_BETA(tpd, rd, 2 )                                                 &
                      + ( FourThirds * JCBN_n_ARRAY(1) )                                            &
                        * CUR_DRV_BETA(tpd, rd, 1, 1 )                                              &
                      + ( JCBN_n_ARRAY(2)/R_SQUARE(rd) )                                            &
                        * CUR_DRV_BETA(tpd, rd, 2, 1 )                                              &
                      + ( JCBN_n_ARRAY(3)/RSIN_SQUARE(td,rd) )                                      &
                        * CUR_DRV_BETA(tpd, rd, 3, 1 )                                              &
                      + ( JCBN_n_ARRAY(2) - COTAN_VAL(td)/3.0_idp )                                 &
                        * CUR_DRV_BETA(tpd, rd, 1, 2 )                                              &
                      + ( 2.0_idp * FourThirds / CUR_R_LOCS(rd) - TwoThirds * JCBN_n_ARRAY(1) )    &
                        * CUR_DRV_BETA(tpd, rd, 2, 2 )                                              &
                      + ( JCBN_n_ARRAY(3)  )                                                        &
                        * CUR_DRV_BETA(tpd, rd, 1, 3 )                                              &
                      + ( 2.0_idp * FourThirds / CUR_R_LOCS(rd) - TwoThirds * JCBN_n_Array(1) )     &
                        * CUR_DRV_BETA(tpd, rd, 3, 3 )



RHS_Terms(tpd, rd, 4) = Beta_Source_Prefix * Block_Source_Si(rd, td, pd, re, te, pe, 2)             &
                    + ( TwoThirds * JCBN_n_ARRAY(2)/R_SQUARE(rd) )                                  &
                        * CUR_VAL_BETA(tpd, rd, 1)                                                  &
                    + ( TwoThirds/R_SQUARE(rd)*( 2.0_idp*COTAN_VAL(td) - JCBN_n_ARRAY(2) ) )        &
                        * CUR_DRV_BETA(tpd, rd, 1, 1)                                               &
                    + ( (JCBN_n_ARRAY(1) - 2.0_idp * TwoThirds/CUR_R_LOCS(rd))/R_SQUARE(rd) )       &
                        * CUR_DRV_BETA(tpd, rd, 2, 1)                                               &
                    + ( TwoThirds/RSIN_SQUARE(td,rd) - 1.0_idp/R_SQUARE(rd)                         &
                       - TwoThirds*COTAN_VAL(td)*JCBN_n_ARRAY(2)/R_SQUARE(rd) )                     &
                        * CUR_VAL_BETA(tpd, rd, 2)                                                  &
                    + ( JCBN_n_ARRAY(1) - 2.0_idp/CUR_R_LOCS(rd) )                                  &
                        * CUR_DRV_BETA(tpd, rd, 1, 2 )                                              &
                    + ( OneThird/R_SQUARE(rd) * ( - 4.0_idp*JCBN_n_ARRAY(2)) )                      &
                        * CUR_DRV_BETA(tpd, rd, 2, 2 )                                              &
                    + ( JCBN_n_ARRAY(3)/RSIN_SQUARE(td, rd) )                                       &
                        * CUR_DRV_BETA(tpd, rd, 3, 2 )                                              &
                    + ( JCBN_n_ARRAY(3)/R_SQUARE(rd) )                                              &
                        * CUR_DRV_BETA(tpd, rd, 2, 3 )                                              &
                    + ( TwoThirds/R_SQUARE(rd)*(3.5_idp*COTAN_VAL(td) - JCBN_n_ARRAY(2) )  )        &
                        * CUR_DRV_BETA(tpd, rd, 3, 3 )





RHS_Terms(tpd, rd, 5) = Beta_Source_Prefix * Block_Source_Si(rd, td, pd, re, te, pe, 3)             &
                    + ( TwoThirds * JCBN_n_ARRAY(3)/RSIN_SQUARE(td, rd) )                           &
                        * CUR_VAL_BETA(tpd, rd, 1 )                                                 &
                    - ( TwoThirds * JCBN_n_ARRAY(3)/RSIN_SQUARE(td, rd) )                           &
                        * CUR_DRV_BETA(tpd, rd, 1, 1 )                                              &
                    + ( (JCBN_n_ARRAY(1) + 2.0_idp*FourThirds/CUR_R_LOCS(rd))/RSIN_SQUARE(td, rd))  &
                        * CUR_DRV_BETA(tpd, rd, 3, 1)                                               &
                    + ( FourThirds*JCBN_n_ARRAY(3)*COTAN_VAL(td)/RSIN_SQUARE(td, rd) )              &
                        * CUR_VAL_BETA(tpd, rd, 2)                                                  &
                    - ( TwoThirds * JCBN_n_ARRAY(3)/RSIN_SQUARE(td, rd) )                           &
                        * CUR_DRV_BETA(tpd, rd, 2, 2 )                                              &
                    + ( (JCBN_n_ARRAY(2) - FourThirds*COTAN_VAL(td)/CUR_R_LOCS(rd))/RSIN_SQUARE(td, rd)) &
                        * CUR_DRV_BETA(tpd, rd, 3, 2 )                                              &
                    + ( JCBN_n_ARRAY(1) - 2.0_idp/CUR_R_LOCS(rd))                                   &
                        * CUR_DRV_BETA(tpd, rd, 1, 3 )                                              &
                    + ( ( JCBN_n_ARRAY(2) - 2.0_idp*COTAN_VAL(td))/R_SQUARE(rd))                    &
                        * CUR_DRV_BETA(tpd, rd, 2, 3 )                                              &
                    + ( FourThirds * JCBN_n_ARRAY(3)/RSIN_SQUARE(td, rd) )                          &
                        * CUR_DRV_BETA(tpd, rd, 3, 3)




END SUBROUTINE Calc_RHS_Terms










END MODULE SubJacobian_Functions_Module_3D
