   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE New_PETSc_Build_Module                                                       !##!
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
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!

!*D*================================!
!                                   !
!           Dependencies            !
!                                   !
!===================================!

USE MPI

USE OMP_LIB

USE Units_Module, &
                ONLY :  GR_Source_Scalar



USE SelfSimilar_Module, &
                ONLY :  SELFSIM_SHIFT_SOL



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
                        Ylm_Table_Block,            &
                        Ylm_Values,                 &
                        Ylm_CC_Values,              &
                        Ylm_dt_Values,              &
                        Ylm_dp_Values,              &
                        Ylm_dtt_Values,             &
                        Ylm_dtp_Values,             &
                        Ylm_dpp_Values,             &
                        Lagrange_Poly_Table,        &
                        LPT_LPT,                    &
                        Coefficient_Vector,         &
                        RHS_Vector,                 &
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
                        BLOCK_NONZEROS,             &
                        Matrix_Location,            &
                        LM_Location






USE Additional_Functions_Module, &
                ONLY :  Lagrange_Poly,              &
                        Spherical_Harmonic,         &
                        Initialize_LGL_Quadrature,  &
                        Map_To_X_Space,             &
                        CFA_Matrix_Map,             &
                        CFA_ALL_Matrix_Map


USE Jacobian_Internal_Functions_Module, &
                ONLY :  JCBN_kappa_FUNCTION_3D_ALL,    &
                        JCBN_BIGK_SUBROUTINE,       &
                        JCBN_BIGK_FUNCTION


USE IO_Functions_Module, &
                ONLY :  Clock_In



IMPLICIT NONE





CONTAINS














!+203+###########################################################################!
!                                                                                !
!                  Calc_3D_SubJcbn_Terms                                         !
!                                                                                !
!################################################################################!
SUBROUTINE Create_Resdiual_Vector(   re, te, pe,                                    &
                                    CUR_VAL_PSI, CUR_DRV_PSI,                       &
                                    CUR_VAL_ALPHAPSI, CUR_DRV_ALPHAPSI,             &
                                    CUR_VAL_BETA, CUR_DRV_BETA, CUR_DDRV_BETA,      &
                                    INT_FACTOR,                                     &
                                    RHS_TERMS,                                      &
                                    CUR_R_LOCS, CUR_T_LOCS, CUR_P_LOCS,             &
                                    R_SQUARE, R_CUBED, R_INVERSE,                   &
                                    SIN_VAL, COS_VAL, CSC_VAL, COTAN_VAL,           &
                                    SIN_SQUARE, CSC_SQUARE, RSIN_SQUARE             )




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


REAL(KIND = idp), INTENT(IN), DIMENSION( 1:NUM_R_QUAD_POINTS,        &
                                         1:NUM_T_QUAD_POINTS,        &
                                         1:NUM_P_QUAD_POINTS         )      ::  Int_Factor


REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_R_QUAD_POINTS)            ::  CUR_R_LOCS
REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_T_QUAD_POINTS)            ::  CUR_T_LOCS
REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_P_QUAD_POINTS)            ::  CUR_P_LOCS


REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_R_QUAD_POINTS,    &
                                        1:NUM_T_QUAD_POINTS     )       ::  RSIN_SQUARE

REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_R_QUAD_POINTS)            ::  R_SQUARE,           &
                                                                            R_CUBED,             &
                                                                            R_INVERSE


REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_T_QUAD_POINTS)            ::  SIN_VAL,            &
                                                                            CSC_VAL,            &
                                                                            SIN_SQUARE,         &
                                                                            CSC_SQUARE,         &
                                                                            COTAN_VAL,          &
                                                                            COS_VAL



REAL(KIND = idp), INTENT(INOUT), DIMENSION( 1:5,                        &
                                            1:NUM_R_QUAD_POINTS,        &
                                            1:NUM_T_QUAD_POINTS,        &
                                            1:NUM_P_QUAD_POINTS         )   ::  RHS_TERMS



REAL(KIND = idp)                                                            ::  REUSED_VALUE

INTEGER                                                                     ::  pd, td, rd,     &
                                                                                i, lm_loc, d,   &
                                                                                current_i_location


REAL(KIND = idp), DIMENSION(1:8)                                            ::  PSI_POWER
REAL(KIND = idp), DIMENSION(1:4)                                            ::  ALPHAPSI_POWER


REAL(KIND = idp), DIMENSION(1:3,1:3)                                        ::  JCBN_kappa_Array
REAL(KIND = idp), DIMENSION(1:3)                                            ::  JCBN_n_ARRAY

REAL(KIND = idp)                                                            ::  JCBN_BIGK_VALUE

COMPLEX(KIND = idp), DIMENSION(1:5)                                         ::  RHS_TMP











!$OMP PARALLEL DEFAULT(none)                                            &
!$OMP PRIVATE( pd, td, rd,                                              &
!$OMP           PSI_POWER, ALPHAPSI_POWER,                              &
!$OMP           JCBN_BIGK_VALUE, JCBN_kappa_ARRAY,JCBN_n_ARRAY,         &
!$OMP           RHS_TMP,                                                &
!$OMP           REUSED_VALUE                                        )   &
!$OMP SHARED( re, te, pe,                                               &
!$OMP           NUM_P_QUAD_POINTS, NUM_T_QUAD_POINTS, NUM_R_QUAD_POINTS,&
!$OMP           CUR_VAL_PSI, CUR_VAL_ALPHAPSI, CUR_VAL_BETA,            &
!$OMP           CUR_DRV_PSI, CUR_DRV_ALPHAPSI, CUR_DRV_BETA,            &
!$OMP           CUR_DDRV_BETA,                                          &
!$OMP           CUR_R_LOCS, CUR_T_LOCS, CUR_P_LOCS,                     &
!$OMP           RSIN_SQUARE, R_SQUARE, R_CUBED, R_INVERSE,              &
!$OMP           SIN_VAL, SIN_SQUARE, COTAN_VAL, COS_VAL,                &
!$OMP           CSC_VAL, CSC_SQUARE,                                    &
!$OMP           GR_Source_Scalar,                                       &
!$OMP           NUM_OFF_DIAGONALS,                                      &
!$OMP           RHS_TERMS,                                              &
!$OMP           Block_SOURCE_E, Block_SOURCE_S, Block_SOURCE_Si,        &
!$OMP           OneThird, OneThirtySecond, FourThirds, TwoThirds,       &
!$OMP           BLOCK_RHS_Vector,                                       &
!$OMP           LM_Length                           )



!$OMP DO SCHEDULE(dynamic), COLLAPSE(3)
DO pd = 1,NUM_P_QUAD_POINTS
 DO td = 1,NUM_T_QUAD_POINTS
  DO rd = 1,NUM_R_QUAD_POINTS




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
                                                CUR_R_LOCS(rd), R_SQUARE(rd), SIN_SQUARE(td), CSC_SQUARE(td),   &
                                                RSIN_SQUARE(rd, td), COTAN_VAL(td)              )




        JCBN_kappa_Array = JCBN_kappa_FUNCTION_3D_ALL( rd, td, pd,                                    &
                                                CUR_R_LOCS(rd), R_SQUARE(rd), R_CUBED(rd),      &
                                                R_INVERSE(rd), RSIN_SQUARE(rd, td),             &
                                                SIN_VAL(td), SIN_SQUARE(td), CSC_SQUARE(td),    &
                                                COS_VAL(td), COTAN_VAL(td),                     &
                                                CUR_VAL_BETA, CUR_DRV_BETA                      )




        JCBN_n_ARRAY(:) = CUR_DRV_ALPHAPSI(:, rd, td, pd) / ALPHAPSI_POWER(1)   &
                            - 7 * CUR_DRV_PSI(:, rd, td, pd )/ PSI_POWER(1)









        RHS_Terms(1, rd, td, pd) = - TwoPi * GR_Source_Scalar                                           &
                                           * Block_Source_E(rd, td, pd, re, te, pe) * PSI_POWER(5)      &
                                   - PSI_POWER(7)/ (16.0_idp * ALPHAPSI_POWER(2)) * JCBN_BIGK_VALUE





        RHS_Terms(2, rd, td, pd) = TwoPi * ALPHAPSI_POWER(1) * PSI_POWER(4)                             &
                        * GR_Source_Scalar * ( Block_Source_E(rd, td, pd, re, te, pe)                   &
                                               + 2.0_idp * Block_Source_S(rd, td, pd, re, te, pe)  )    &
                        + 7.0_idp*PSI_POWER(6)/ (16.0_idp * ALPHAPSI_POWER(1)) * JCBN_BIGK_VALUE






        RHS_Terms(3, rd, td, pd) = 16.0_idp * pi                                                &
                          * ALPHAPSI_POWER(1)                                                   &
                          * PSI_POWER(3)                                                        &
                          * GR_Source_Scalar                                                    &
                          * Block_Source_Si(rd, td, pd, re, te, pe, 1)                          &
                     + 2.0_idp/R_SQUARE(rd) * CUR_VAL_BETA(1,rd,td,pd)                          &  ! New
                     + 2.0_idp*COTAN_VAL(td)/CUR_R_LOCS(rd) * CUR_VAL_BETA(2,rd,td,pd)          &  ! New
                     + 2.0_idp/CUR_R_LOCS(rd) * CUR_DRV_BETA(2,2,rd,td,pd)                      &  ! New
                     + 2.0_idp/CUR_R_LOCS(rd) * CUR_DRV_BETA(3,3,rd,td,pd)                      &  ! New
                     - OneThird                                                                 &
                          * ( CUR_DDRV_BETA(1, 1, 1, rd, td, pd )                               &
                            + CUR_DDRV_BETA(1, 2, 2, rd, td, pd )                               &
                            + CUR_DDRV_BETA(1, 3, 3, rd, td, pd )                               &
                            + 2.0_idp/CUR_R_LOCS(rd) * CUR_DRV_BETA(1, 1, rd, td, pd)           &
                            + COTAN_VAL(td)          * CUR_DRV_BETA(1, 2, rd, td, pd)           &
                            - 2.0_idp /R_SQUARE(rd)  * CUR_VAL_BETA(1, rd, td, pd )             &
                            )                                                                   &
                      + JCBN_n_ARRAY(1)* JCBN_kappa_Array(1,1)                                     &
                      + JCBN_n_ARRAY(2)* JCBN_kappa_Array(2,1)                                     &
                      + JCBN_n_ARRAY(3)* JCBN_kappa_Array(3,1)

!        PRINT*,"RHS_TERM 3 ",rd ,td ,pd
!        PRINT*,RHS_TERMS(1:3, rd,td,pd)
!        PRINT*,16.0_idp * pi                                                &
!                          * ALPHAPSI_POWER(1)                                                   &
!                          * PSI_POWER(3)                                                        &
!                          * GR_Source_Scalar                                                    &
!                          * Block_Source_Si(rd, td, pd, re, te, pe, 1)
!        PRINT*,2.0_idp/R_SQUARE(rd) * CUR_VAL_BETA(1,rd,td,pd)                          &  ! New
!                     + 2.0_idp*COTAN_VAL(td)/CUR_R_LOCS(rd) * CUR_VAL_BETA(2,rd,td,pd)          &  ! New
!                     + 2.0_idp/CUR_R_LOCS(rd) * CUR_DRV_BETA(2,2,rd,td,pd)                      &  ! New
!                     + 2.0_idp/CUR_R_LOCS(rd) * CUR_DRV_BETA(3,3,rd,td,pd)
!        PRINT*, - OneThird                                                                 &
!                          * ( CUR_DDRV_BETA(1, 1, 1, rd, td, pd )                               &
!                            + CUR_DDRV_BETA(1, 2, 2, rd, td, pd )                               &
!                            + CUR_DDRV_BETA(1, 3, 3, rd, td, pd )                               &
!                            + 2.0_idp/CUR_R_LOCS(rd) * CUR_DRV_BETA(1, 1, rd, td, pd)           &
!                            + COTAN_VAL(td)          * CUR_DRV_BETA(1, 2, rd, td, pd)           &
!                            - 2.0_idp /R_SQUARE(rd)  * CUR_VAL_BETA(1, rd, td, pd )             &
!                            )
!        PRINT*,JCBN_n_ARRAY(1)* JCBN_kappa_Array(1,1)                                     &
!                      + JCBN_n_ARRAY(2)* JCBN_kappa_Array(2,1)                                     &
!                      + JCBN_n_ARRAY(3)* JCBN_kappa_Array(3,1)
!        PRINT*,"----"
!        PRINT*, CUR_DDRV_BETA(1, 1, 1, rd, td, pd )         &
!              , 2.0_idp/CUR_R_LOCS(rd) * CUR_DRV_BETA(1, 1, rd, td, pd)     &
!              , - 2.0_idp /R_SQUARE(rd)  * CUR_VAL_BETA(1, rd, td, pd )     &
!              , CUR_R_LOCS(rd)


        RHS_Terms(4, rd, td, pd) =  16.0_idp * pi                                               &
                           * ALPHAPSI_POWER(1)                                                  &
                           * PSI_POWER(3)                                                       &
                           * GR_Source_Scalar                                                   &
                           * Block_Source_Si(rd, td, pd, re, te, pe, 2)                         &
                        - 2.0_idp/CUR_R_LOCS(rd) * CUR_DRV_BETA(1,2, rd, td, pd)                &            !NEW
                        - (1.0_idp - COTAN_VAL(td)*COTAN_VAL(td))/R_SQUARE(rd) * CUR_VAL_BETA(2,rd,td,pd) &  !NEW
                        - 2.0_idp/R_CUBED(rd) * CUR_DRV_BETA(2,1, rd, td, pd)                   &            !NEW
                        + (2.0_idp*COTAN_VAL(td))/R_SQUARE(rd) * CUR_DRV_BETA(3,3,rd,td,pd)     &            !NEW
                        - ( OneThird /R_SQUARE(rd) )                                            &
                             * ( CUR_DDRV_BETA(2, 1, 1, rd, td, pd )                            &
                               + CUR_DDRV_BETA(2, 2, 2, rd, td, pd )                            &
                               + CUR_DDRV_BETA(2, 3, 3, rd, td, pd )                            &
                               + 2.0_idp/CUR_R_LOCS(rd) * CUR_DRV_BETA(2, 1, rd, td, pd)        &
                               + COTAN_VAL(td) * CUR_DRV_BETA(2, 2, rd, td, pd)                 &
                               - CSC_SQUARE(td) * CUR_VAL_BETA(2, rd, td, pd)                   &
                               )                                                                &
                        + JCBN_n_ARRAY(1) * JCBN_kappa_Array(1,2)                                  &
                        + JCBN_n_ARRAY(2) * JCBN_kappa_Array(2,2)                                  &
                        + JCBN_n_ARRAY(3) * JCBN_kappa_Array(3,2)





        RHS_Terms(5, rd, td, pd) = 16.0_idp * pi                                            &
                          * ALPHAPSI_POWER(1)                                               &
                          * PSI_POWER(3)                                                    &
                          * GR_Source_Scalar                                                &
                          * Block_Source_Si(rd, td, pd, re, te, pe, 3)                      &
                        - 2.0_idp/(CUR_R_LOCS(rd)*RSIN_SQUARE(rd,td)) * CUR_DRV_BETA(3,1,rd,td,pd) & ! NEW
                        - COTAN_VAL(td)/RSIN_SQUARE(rd,td) * CUR_DRV_BETA(3,2,rd,td,pd)     &  ! NEW
                        - 2.0_idp/CUR_R_LOCS(rd) * CUR_DRV_BETA(1,3,rd,td,pd)               &  ! NEW
                        - 2.0_idp*COTAN_VAL(td) / R_SQUARE(rd) * CUR_DRV_BETA(2,3,rd,td,pd) &  ! NEW
                        - ( OneThird /RSIN_SQUARE(rd, td ))                                 &
                             * ( CUR_DDRV_BETA(3, 1, 1, rd, td, pd )                        &
                               + CUR_DDRV_BETA(3, 2, 2, rd, td, pd )                        &
                               + CUR_DDRV_BETA(3, 3, 3, rd, td, pd )                        &
                               + 2.0_idp/CUR_R_LOCS(rd) * CUR_DRV_BETA(3, 1, rd, td, pd)    &
                               + COTAN_VAL(td) * CUR_DRV_BETA(3, 2, rd, td, pd)             &
                             )                                                              &
                        + JCBN_n_ARRAY(1) * JCBN_kappa_Array(1,3)                              &
                        + JCBN_n_ARRAY(2) * JCBN_kappa_Array(2,3)                              &
                        + JCBN_n_ARRAY(3) * JCBN_kappa_Array(3,3)



        END DO ! rd loop
    END DO  ! td loop
END DO  ! pd loop
!$OMP END DO

!$OMP DO SCHEDULE(dynamic), COLLAPSE(2)
DO d = 0,DEGREE
    DO lm_loc = 0,LM_LENGTH-1

        Current_i_Location = CFA_ALL_Matrix_Map(1, lm_loc, re, d)

        RHS_TMP = 0.0_idp
        DO pd = 1,NUM_P_QUAD_POINTS
            DO td = 1,NUM_T_QUAD_POINTS
                DO rd = 1,NUM_R_QUAD_POINTS


                    RHS_TMP(1:5) =  RHS_TMP(1:5)                                        &
                                    + RHS_TERMS(1:5, rd, td, pd)                        &
                                        * Ylm_CC_Values(lm_loc, td, pd, te, pe)         &
                                        * Lagrange_Poly_Table(d, rd, 0)                 &
                                        * Int_Factor(rd, td, pd)





                END DO  ! rd Loop
            END DO  ! td Loop
        END DO  ! pd Loop


        Block_RHS_Vector(Current_i_Location:Current_i_Location+4)                             &
            = Block_RHS_Vector(Current_i_Location:Current_i_Location+4)                       &
            + RHS_TMP(1:5)


    END DO  ! l Loop
END DO  ! d Loop
!$OMP END DO




!$OMP END PARALLEL





END SUBROUTINE Create_Resdiual_Vector





END MODULE New_PETSc_Build_Module
