   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE GMRES_Master_Build_Module                                                    !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
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


USE MPI

USE Poseidon_Constants_Module, &
                ONLY :  idp,                &
                        pi,                 &
                        TwoPi,              &
                        OneThird,           &
                        TwoThirds,          &
                        FourThirds,         &
                        OneThirtySecond

USE Poseidon_BC_Module, &
                ONLY :  CFA_3D_Apply_BCs_Part1,         &
                        CFA_3D_Apply_BCs_Part2


USE Jacobian_Internal_Functions_Module, &
                ONLY :  JCBN_kappa_FUNCTION_3D_ALL,     &
                        JCBN_BIGK_FUNCTION

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
                        NUM_R_ELEMS_PER_SHELL,      &
                        OUTPUT_RHS_VECTOR_FLAG

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
                        Ylm_Table_Block,            &
                        Ylm_Values,                 &
                        Ylm_dt_Values,              &
                        Ylm_dp_Values,              &
                        Ylm_CC_Values,              &
                        Ylm_CC_dt_Values,           &
                        Ylm_CC_dp_Values,           &
                        Lagrange_Poly_Table,        &
                        LPT_LPT,                    &
                        Coefficient_Vector,         &
                        myID_Shell,                 &
                        PROB_DIM,                   &
                        LM_LENGTH,                  &
                        myShell,                    &
                        Block_STF_Mat,              &
                        Subshell_Prob_Dim,          &
                        Num_off_Diagonals,          &
                        M_Values,                   &
                        ULM_Length,                 &
                        Matrix_Location,            &
                        myID_Subshell,              &
                        Poseidon_Comm_Petsc
        
USE Poseidon_Mapping_Functions_Module, &
                ONLY :  Map_To_X_Space,                 &
                        CFA_ALL_Matrix_Map

USE Poseidon_IO_Module, &
                ONLY :  Clock_In

USE Poseidon_LinSys_IO_Module, &
                ONLY :  OUTPUT_RHS_VECTOR_Parts

USE SubJacobian_Functions_Module_3D, &
                ONLY :  Calc_RHS_Terms

IMPLICIT NONE



REAL(KIND = idp), ALLOCATABLE, DIMENSION( :,:,: )       :: RR_Factor
REAL(KIND = idp), ALLOCATABLE, DIMENSION( :,:,: )       :: DRR_Factor
REAL(KIND = idp), ALLOCATABLE, DIMENSION( :,:,: )       :: RDR_Factor
REAL(KIND = idp), ALLOCATABLE, DIMENSION( :,:,: )       :: DRDR_Factor

COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION( :,:,: )    :: TP_TP_Factor
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION( :,:,: )    :: dTP_TP_Factor
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION( :,:,: )    :: TdP_TP_Factor
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION( :,:,: )    :: TP_dTP_Factor
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION( :,:,: )    :: TP_TdP_Factor
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION( :,:,: )    :: dTP_dTP_Factor
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION( :,:,: )    :: dTP_TdP_Factor
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION( :,:,: )    :: TdP_dTP_Factor
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION( :,:,: )    :: TdP_TdP_Factor

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: R_Int_Weights
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: TP_Int_Weights

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CUR_R_LOCS
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CUR_T_LOCS
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CUR_P_LOCS


REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: R_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: R_CUBED

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: SIN_VAL
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: SIN_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: TP_SIN_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: COS_VAL
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: COS_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CSC_VAL
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CSC_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: COTAN_VAL

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)           :: RSIN_SQUARE

COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)          :: PHI_EXP
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)          :: PHI_TWOEXP


REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)           :: CUR_VAL_PSI
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)           :: CUR_VAL_ALPHAPSI
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: CUR_VAL_BETA

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: CUR_DRV_PSI
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: CUR_DRV_ALPHAPSI
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:,:)       :: CUR_DRV_BETA

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)           :: Beta_DRV_Trace

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: RHS_TERMS


CONTAINS




!+202+###########################################################################!
!                                                                                !
!                  Calc_Equation          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Equation( Eq, Coeff_Vec, Coeff_Length )

COMPLEX(KIND = idp), DIMENSION(0:Coeff_Length-1), INTENT(OUT)   :: Eq
COMPLEX(KIND = idp), DIMENSION(0:Coeff_Length-1), INTENT(IN)    :: Coeff_Vec
INTEGER, INTENT(IN)                                             :: Coeff_Length







END SUBROUTINE Calc_Equation









!+101+###########################################################################!
!                                                                                !
!           CFA_3D_Master_Build                                                  !
!                                                                                !
!################################################################################!
SUBROUTINE CFA_GMRES_Eq_Build(Eq, Coeff_Vec, Coeff_Length)

COMPLEX(KIND = idp), DIMENSION(0:Coeff_Length-1), INTENT(OUT)   :: Eq
COMPLEX(KIND = idp), DIMENSION(0:Coeff_Length-1), INTENT(IN)    :: Coeff_Vec
INTEGER, INTENT(IN)                                             :: Coeff_Length

REAL(KIND = idp)        :: timea, timeb, timec
INTEGER                 :: i,j,k
INTEGER                 :: here, re, d, dp, lm_loc, lpmp_loc, F, u, l, m
INTEGER                 :: i_loc, j_loc, m_loc

timea = 0.0_idp
timeb = 0.0_idp
timec = 0.0_idp


!*!
!*! Allocate module variables
!*!
CALL Allocate_GMRES_Eq_Variables()



!*!
!*! Alter the Coefficient Vector to Reflect Boundary Conditions
!*!

CALL CFA_3D_Apply_BCs_Part1()
timec = MPI_Wtime()
CALL Clock_In(timec-timeb, 4)




!*!
!*! Create the Non-Laplacian Contributions to the Equation Vector
!*!
timeb = MPI_Wtime()
CALL CREATE_GMRES_Equations( Eq, Coeff_Vec, Coeff_Length )
timea = MPI_Wtime()
CALL Clock_In(timea-timeb, 9)




!*!
!*! Add the Laplacian contribution to the Equation Vector
!*!
timeb = MPI_Wtime()
CALL FINISH_GMRES_Equations( Eq, Coeff_Length)
timea = MPI_Wtime()
CALL Clock_In(timea-timeb, 12)



!CALL Jacobi_Type_B_PC()



!*!
!*! Alter the Jacobian Matrix to Reflect Boundary Conditions
!*!
timeb = MPI_Wtime()
CALL CFA_3D_Apply_BCs_Part2()
timec = MPI_Wtime()
CALL Clock_In(timec-timeb, 13)




!*!
!*! Deallocate module variables
!*!
CALL Deallocate_GMRES_Eq_Variables()


END SUBROUTINE CFA_GMRES_Eq_Build












!+201+###############################################################################!
!                                                                                    !
!                                  CREATE_3D_NONLAPLACIAN_SOE                        !
!                                                                                    !
!####################################################################################!
SUBROUTINE CREATE_GMRES_Equations( Eq, Coeff_Vec, Coeff_Length )

COMPLEX(KIND = idp), DIMENSION(0:Coeff_Length-1), INTENT(OUT)   :: Eq
COMPLEX(KIND = idp), DIMENSION(0:Coeff_Length-1), INTENT(IN)    :: Coeff_Vec
INTEGER, INTENT(IN)                                             :: Coeff_Length

INTEGER                                                         ::  Local_re,       &
                                                                    Local_te,       &
                                                                    Local_pe,       &
                                                                    Global_re,      &
                                                                    Global_te,      &
                                                                    Global_pe,      &
                                                                    d, l, m,        &
                                                                    dp, lp, mp,     &
                                                                    rd, td, pd, tpd




COMPLEX(KIND = idp)                                             ::  REUSED_VALUE

REAL(KIND = idp), DIMENSION(1:3,1:3)                            ::  JCBN_kappa_Array
REAL(KIND = idp), DIMENSION(1:3)                                ::  JCBN_n_ARRAY

REAL(KIND = idp)                                                ::  JCBN_BIGK_VALUE



INTEGER, DIMENSION(1:NUM_T_QUAD_POINTS)                         ::  here


REAL(KIND = idp)                                                ::  TWOOVER_DELTAR,    &
                                                                    deltar_overtwo,     &
                                                                    deltat_overtwo,     &
                                                                    deltap_overtwo




INTEGER                                                         ::  Block_T_Begin, Block_P_Begin




REAL(KIND = idp)                                    ::  timea, timeb, timec, timed, timee,  &
                                                        time_CURVALS, time_SJT, time_JCBNM, &
                                                        time_RHS


INTEGER                                             ::  ierr, i
INTEGER                                             ::  CUR_LOC

time_CURVALS = 0.0_idp
time_SJT = 0.0_idp
time_JCBNM = 0.0_idp
time_RHS = 0.0_idp





!
!   Identifiy process's first theta and phi elements
!
Block_T_Begin = MOD(myID_Shell,NUM_BLOCK_THETA_ROWS)*NUM_T_ELEMS_PER_BLOCK
Block_P_Begin = (myID_Shell/NUM_BLOCK_THETA_ROWS)*NUM_P_ELEMS_PER_BLOCK






!
!   Move through phi
!
DO Local_pe = 0, NUM_P_ELEMS_PER_BLOCK-1


    Global_pe = Block_P_Begin + Local_pe

    deltap_overtwo = (plocs(Global_pe + 1) - plocs(Global_pe))/2.0_idp
    CUR_P_LOCS(:) = deltap_overtwo * (INT_P_LOCATIONS(:)+1.0_idp) + plocs(Global_pe)

    PHI_EXP(:) = EXP( CMPLX(0, -CUR_P_LOCS(:), KIND = idp) )
    PHI_TWOEXP(:) = EXP( CMPLX(0, -2.0_idp*CUR_P_LOCS(:), KIND = idp) )


!    PRINT*,"Here",CUR_P_LOCS

    DO Local_te = 0, NUM_T_ELEMS_PER_BLOCK-1


        Global_te = Block_T_Begin + Local_te

        deltat_overtwo = (tlocs(Global_te + 1) - tlocs(Global_te))/2.0_idp
        CUR_T_LOCS(:) = deltat_overtwo * (INT_T_LOCATIONS(:)+1.0_idp) + tlocs(Global_te)
        

        COTAN_VAL(:) = 1.0_idp/DTAN(CUR_T_LOCS(:))
        CSC_VAL(:) = 1.0_idp/DSIN(CUR_T_LOCS(:))
        SIN_VAL(:) = DSIN(CUR_T_LOCS(:))
        COS_VAL(:) = DCOS(CUR_T_LOCS(:))

        COS_SQUARE(:) = COS_VAL(:)*COS_VAL(:)
        SIN_SQUARE(:) = SIN_VAL(:)*SIN_VAL(:)
        CSC_SQUARE(:) = CSC_VAL(:)*CSC_VAL(:)

        DO td = 1,NUM_T_QUAD_POINTS
            DO pd = 1,NUM_P_QUAD_POINTS
                tpd = (td-1)*NUM_P_QUAD_POINTS + pd
                TP_SIN_SQUARE(tpd) = SIN_SQUARE(td)
            END DO
        END DO

        !
        !    Move through radius
        !
        DO Local_re = 0,NUM_R_ELEMS_PER_BLOCK-1

            Global_re = myShell*NUM_R_ELEMS_PER_SHELL + Local_re


            deltar_overtwo = (rlocs(Global_re + 1) - rlocs(Global_re))/2.0_idp
            TWOOVER_DELTAR = 1.0_idp/deltar_overtwo
            CUR_R_LOCS(:) = deltar_overtwo * (INT_R_LOCATIONS(:)+1.0_idp) + rlocs(Global_re)


            R_SQUARE(:) = CUR_R_LOCS(:)*CUR_R_LOCS(:)
            R_CUBED(:) = R_SQUARE(:)*CUR_R_LOCS(:)


            DO rd = 1,NUM_R_QUAD_POINTS

                RSIN_SQUARE(:,rd) = R_SQUARE(rd)*SIN_SQUARE(:)

            END DO

!            IF ( Local_re == 0) THEN
!                WRITE(*,*)"R_INNER^2 = ",R_SQUARE(1)
!            ELSE IF ( Local_RE == NUM_R_ELEMS_PER_BLOCK-1) THEN
!                WRITE(*,*)"R_OUTER^2 = ",R_SQUARE(NUM_R_QUAD_POINTS)
!            END IF


            !*!
            !*! Calculate Current Values of CFA Varaiables and their Deriviatives
            !*!
            timea = MPI_Wtime()
            CALL Calc_3D_Current_Values(Global_re , Local_te  , Local_pe,               &
                                        DELTAR_OVERTWO, DELTAT_OVERTWO, DELTAP_OVERTWO, &
                                        Coeff_Vec, Coeff_Length )


            timeb = MPI_Wtime()



            CALL  Calc_GMRES_Terms( Local_re, Local_te, Local_pe )



            !*!
            !*! Create the Residual Vector ( Sans Laplacian Contribution )
            !*!
            CALL Calc_GMRES_Equations( Eq, Coeff_Length, Local_re, Local_te, Local_pe, DELTAR_OVERTWO )

            timed = MPI_Wtime()






            !*!
            !*! Update Timer Values
            !*!
            time_CurVals = time_CurVals + timeb - timea
            time_RHS = time_RHS + timed - timeb



        END DO  ! pe Loop
    END DO  ! te Loop
END DO  ! re Loop


CALL Clock_In(time_CurVals, 5)
CALL Clock_In(time_RHS, 7)





END SUBROUTINE CREATE_GMRES_Equations








!+202+###########################################################################!
!                                                                                !
!                  Calc_Current_Values          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_3D_Current_Values( re, te, pe,                                  &
                                    DELTAR_OVERTWO,                             &
                                    DELTAT_OVERTWO,                             &
                                    DELTAP_OVERTWO,                             &
                                    Coeff_Vec,                                  &
                                    Coeff_Length                                )

INTEGER, INTENT(IN)                                             ::  re, te, pe


REAL(KIND = idp), INTENT(IN)                                    ::  DELTAR_OVERTWO,     &
                                                                    DELTAT_OVERTWO,     &
                                                                    DELTAP_OVERTWO

COMPLEX(KIND = idp), DIMENSION(0:Coeff_Length-1), INTENT(IN)    :: Coeff_Vec
INTEGER, INTENT(IN)                                             :: Coeff_Length

COMPLEX(KIND = idp), DIMENSION(1:5)                             ::  Tmp_U_Value,        &
                                                                    Tmp_U_R_DRV_Value,  &
                                                                    Tmp_U_T_DRV_Value,  &
                                                                    Tmp_U_P_DRV_Value



INTEGER                                                         ::  l, m, d, dp,        &
                                                                    rd, td, pd, tpd,    &
                                                                    ui


INTEGER                                                         ::  Here, There
INTEGER                                                         ::  lm_loc, lpmp_loc

COMPLEX(KIND = idp), DIMENSION(1:5)                             ::  Local_Coefficients



                          !                                                 !
                         !!                                                 !!
                        !!!          Initialize Local Quadratures           !!!
                         !!                                                 !!
                          !                                                 !




R_Int_Weights(:) = DELTAR_OVERTWO * R_SQUARE(:) * INT_R_WEIGHTS(:)

DO td = 1,NUM_T_QUAD_POINTS
    DO pd = 1,NUM_P_QUAD_POINTS
        TP_Int_Weights( (td-1)*NUM_P_QUAD_POINTS + pd ) = SIN_VAL(td)                           &
                                                        * DELTAT_OVERTWO * INT_T_WEIGHTS(td)    &
                                                        * DELTAP_OVERTWO * INT_P_WEIGHTS(pd)
    END DO
END DO

DO d = 0,DEGREE
    DO dp = 0,DEGREE

        RR_Factor(:, d, dp)      = R_Int_Weights(:) * LPT_LPT(:, d, dp, 0, 0)
        DRR_Factor(:, d, dp)     = R_Int_Weights(:) * LPT_LPT(:, d, dp, 1, 0) / DELTAR_OVERTWO
        RDR_Factor(:, d, dp)     = R_Int_Weights(:) * LPT_LPT(:, d, dp, 0, 1) / DELTAR_OVERTWO
        DRDR_Factor(:, d, dp)    = R_Int_Weights(:) * LPT_LPT(:, d, dp, 1, 1) / (DELTAR_OVERTWO * DELTAR_OVERTWO)


    END DO
END DO



DO lm_loc = 0,LM_LENGTH-1
    DO lpmp_loc = 0,LM_LENGTH-1
    
        TP_TP_Factor( :, lm_loc, lpmp_loc )  = Ylm_Values( lm_loc, :, te, pe )              &
                                                * Ylm_CC_Values( :, lpmp_loc, te, pe)       &
                                                * TP_Int_Weights(:)

        dTP_TP_Factor( :, lm_loc, lpmp_loc )  = Ylm_dt_Values( lm_loc, :, te, pe )           &
                                                * Ylm_CC_Values( :, lpmp_loc, te, pe)       &
                                                * TP_Int_Weights(:)

        TdP_TP_Factor( :, lm_loc, lpmp_loc )  = Ylm_dp_Values( lm_loc, :, te, pe )           &
                                                * Ylm_CC_Values( :, lpmp_loc, te, pe)       &
                                                * TP_Int_Weights(:)

        TP_dTP_Factor(:, lm_loc, lpmp_loc ) = Ylm_Values( lm_loc, :, te, pe )                &
                                                * Ylm_CC_dt_Values( :, lpmp_loc, te, pe)    &
                                                * TP_Int_Weights(:)

        dTP_dTp_Factor(:, lm_loc, lpmp_loc ) = Ylm_dt_Values( lm_loc, :, te, pe )             &
                                                * Ylm_CC_dt_Values( :, lpmp_loc, te, pe)    &
                                                * TP_Int_Weights(:)

        TP_TdP_Factor(:, lm_loc, lpmp_loc ) = Ylm_Values( lm_loc, :, te, pe )                &
                                                * Ylm_CC_dp_Values( :, lpmp_loc, te, pe)    &
                                                * TP_Int_Weights(:)

        TdP_TdP_Factor(:, lm_loc, lpmp_loc ) = Ylm_dp_Values( lm_loc, :, te, pe )             &
                                                * Ylm_CC_dp_Values( :, lpmp_loc, te, pe)    &
                                                * TP_Int_Weights(:)

        dTP_TdP_Factor(:, lm_loc, lpmp_loc ) = Ylm_dt_Values( lm_loc, :, te, pe )             &
                                                * Ylm_CC_dp_Values( :, lpmp_loc, te, pe)    &
                                                * TP_Int_Weights(:)

        TdP_dTP_Factor(:, lm_loc, lpmp_loc ) = Ylm_dp_Values( lm_loc, :, te, pe )             &
                                                * Ylm_CC_dt_Values( :, lpmp_loc, te, pe)    &
                                                * TP_Int_Weights(:)

    END DO
END DO





DO rd = 1,NUM_R_QUAD_POINTS

    DO tpd = 1,NUM_TP_QUAD_POINTS

         !                                                  !
        !!   Set/Reset temporary value holder to zero.      !!
         !                                                  !
        Tmp_U_Value = 0.0_idp
        Tmp_U_R_DRV_Value = 0.0_idp
        Tmp_U_T_DRV_Value = 0.0_idp
        Tmp_U_P_DRV_Value = 0.0_idp


        DO d = 0,DEGREE

            DO ui = 1,NUM_CFA_VARS

                Here = CFA_All_Matrix_Map(ui, 0, re, d)
                There = Here + LM_LENGTH - 1

!                PRINT*,Coefficient_Vector(Here:There)

                TMP_U_Value(ui)         = TMP_U_Value(ui)                           &
                                        + SUM( Coeff_Vec( Here:There )              &
                                        * Ylm_Values( :, tpd, te, pe )       )      &
                                        * Lagrange_Poly_Table( d, rd, 0 )

                TMP_U_R_DRV_Value(ui)   = TMP_U_R_DRV_Value(ui)                     &
                                        + SUM( Coeff_Vec( Here:There )              &
                                        * Ylm_Values( :, tpd, te, pe )       )      &
                                        * Lagrange_Poly_Table( d, rd, 1 )           &
                                        / DELTAR_OVERTWO


!                IF (ui == 3 ) THEN
!
!                    PRINT*,SUM( Coefficient_Vector( Here:There )     &
!                    * Ylm_Values( :, tpd, te, pe )       )      &
!                    * Lagrange_Poly_Table( d, rd, 1 )           &
!                    / DELTAR_OVERTWO,                           &
!                    Coefficient_Vector( Here:There ),            &
!                    Lagrange_Poly_Table( d, rd, 1 )/ DELTAR_OVERTWO
!
!                END IF



                TMP_U_T_DRV_Value(ui)   = TMP_U_T_DRV_Value(ui)                     &
                                        + SUM( Coeff_Vec( Here:There )              &
                                        * Ylm_dt_Values( :, tpd, te, pe)     )      &
                                        * Lagrange_Poly_Table( d, rd, 0)

                TMP_U_P_DRV_Value(ui)   = TMP_U_P_DRV_Value(ui)                     &
                                        + SUM( Coeff_Vec( Here:There )              &
                                        * Ylm_dp_Values( :, tpd, te, pe)     )      &
                                        * Lagrange_Poly_Table( d, rd, 0)

            END DO ! ui Loop

        END DO  !   d Loop


        CUR_VAL_PSI( tpd, rd )         = REAL(Tmp_U_Value(1), KIND = idp)
        CUR_DRV_PSI( tpd, rd, 1 )      = REAL(Tmp_U_R_DRV_Value(1), KIND = idp)
        CUR_DRV_PSI( tpd, rd, 2 )      = REAL(Tmp_U_T_DRV_Value(1), KIND = idp)
        CUR_DRV_PSI( tpd, rd, 3 )      = REAL(Tmp_U_P_DRV_Value(1), KIND = idp)


        CUR_VAL_ALPHAPSI( tpd, rd )    = REAL(Tmp_U_Value(2), KIND = idp)
        CUR_DRV_ALPHAPSI( tpd, rd, 1 ) = REAL(Tmp_U_R_DRV_Value(2), KIND = idp)
        CUR_DRV_ALPHAPSI( tpd, rd, 2 ) = REAL(Tmp_U_T_DRV_Value(2), KIND = idp)
        CUR_DRV_ALPHAPSI( tpd, rd, 3 ) = REAL(Tmp_U_P_DRV_Value(2), KIND = idp)


        CUR_VAL_BETA( tpd, rd, 1 )     = REAL(Tmp_U_Value(3), KIND = idp)
        CUR_VAL_BETA( tpd, rd, 2 )     = REAL(Tmp_U_Value(4), KIND = idp)
        CUR_VAL_BETA( tpd, rd, 3 )     = REAL(Tmp_U_Value(5), KIND = idp)


        CUR_DRV_BETA( tpd, rd, 1, 1 )  = REAL(Tmp_U_R_DRV_Value(3), KIND = idp)
        CUR_DRV_BETA( tpd, rd, 2, 1 )  = REAL(Tmp_U_T_DRV_Value(3), KIND = idp)
        CUR_DRV_BETA( tpd, rd, 3, 1 )  = REAL(Tmp_U_P_DRV_Value(3), KIND = idp)
        CUR_DRV_BETA( tpd, rd, 1, 2 )  = REAL(Tmp_U_R_DRV_Value(4), KIND = idp)
        CUR_DRV_BETA( tpd, rd, 2, 2 )  = REAL(Tmp_U_T_DRV_Value(4), KIND = idp)
        CUR_DRV_BETA( tpd, rd, 3, 2 )  = REAL(Tmp_U_P_DRV_Value(4), KIND = idp)
        CUR_DRV_BETA( tpd, rd, 1, 3 )  = REAL(Tmp_U_R_DRV_Value(5), KIND = idp)
        CUR_DRV_BETA( tpd, rd, 2, 3 )  = REAL(Tmp_U_T_DRV_Value(5), KIND = idp)
        CUR_DRV_BETA( tpd, rd, 3, 3 )  = REAL(Tmp_U_P_DRV_Value(5), KIND = idp)

        Beta_DRV_Trace( tpd, rd )      = REAL(Tmp_U_R_DRV_Value(3), KIND = idp)     &
                                       + REAL(Tmp_U_T_DRV_Value(4), KIND = idp)     &
                                       + REAL(Tmp_U_P_DRV_Value(5), KIND = idp)
                

    END DO  !   tpd Loop

END DO  !   rd Loop



END SUBROUTINE Calc_3D_Current_Values






!+203+###########################################################################!
!                                                                                !
!                  Calc_3D_SubJcbn_Terms                                         !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_GMRES_Terms( re, te, pe )




INTEGER, INTENT(IN)                                                     ::  re, te, pe


REAL(KIND = idp)                                                        ::  REUSED_VALUE

INTEGER                                                                 ::  pd, td, rd,     &
                                                                            i, tpd


REAL(KIND = idp), DIMENSION(1:11)                                       ::  PSI_POWER
REAL(KIND = idp), DIMENSION(1:4)                                        ::  ALPHAPSI_POWER


REAL(KIND = idp), DIMENSION(1:3,1:3)                                    ::  JCBN_kappa_Array
REAL(KIND = idp), DIMENSION(1:3)                                        ::  JCBN_n_ARRAY

REAL(KIND = idp)                                                        ::  JCBN_BIGK_VALUE





!$OMP PARALLEL DEFAULT(none)                                            &
!$OMP PRIVATE( pd, td, rd, tpd, i,                                      &
!$OMP           PSI_POWER, ALPHAPSI_POWER,                              &
!$OMP           JCBN_BIGK_VALUE, JCBN_kappa_ARRAY,JCBN_n_ARRAY,         &
!$OMP           REUSED_VALUE                                        )   &
!$OMP SHARED( re, te, pe,                                               &
!$OMP           NUM_P_QUAD_POINTS, NUM_T_QUAD_POINTS, NUM_R_QUAD_POINTS,&
!$OMP           CUR_VAL_PSI, CUR_VAL_ALPHAPSI, CUR_VAL_BETA,            &
!$OMP           CUR_DRV_PSI, CUR_DRV_ALPHAPSI, CUR_DRV_BETA,            &
!$OMP           CUR_R_LOCS, CUR_T_LOCS, CUR_P_LOCS,                     &
!$OMP           SubJacobian_EQ1_Term, SubJacobian_EQ2_Term,             &
!$OMP           SubJacobian_EQ3_Term, SubJacobian_EQ4_Term,             &
!$OMP           SubJacobian_EQ5_Term,                                   &
!$OMP           RSIN_SQUARE, R_SQUARE, R_CUBED,                         &
!$OMP           SIN_VAL, SIN_SQUARE, COTAN_VAL, COS_VAL,                &
!$OMP           CSC_VAL, CSC_SQUARE,                                    &
!$OMP           GR_Source_Scalar,                                       &
!$OMP           RHS_TERMS,                                              &
!$OMP           Block_SOURCE_E, Block_SOURCE_S, Block_SOURCE_Si,        &
!$OMP           OneThird, OneThirtySecond, FourThirds, TwoThirds,       &
!$OMP           LM_Length                           )



!$OMP DO SCHEDULE(dynamic), COLLAPSE(3)
DO rd = 1,NUM_R_QUAD_POINTS
 DO td = 1,NUM_T_QUAD_POINTS
  DO pd = 1,NUM_P_QUAD_POINTS

        tpd = (td-1)*NUM_P_QUAD_POINTS + pd


        PSI_POWER(1) = CUR_VAL_PSI( tpd, rd)
        DO i = 2,11
            PSI_POWER(i) = PSI_POWER(i-1)*PSI_POWER(1)
        END DO


        ALPHAPSI_POWER(1) = CUR_VAL_ALPHAPSI( tpd, rd)
        DO i = 2,4
            ALPHAPSI_POWER(i) = ALPHAPSI_POWER(i-1)*ALPHAPSI_POWER(1)
        END DO


        ! K_{ij}K^{ij} = Psi^{14}/AlphaPsi^{2} * BIGK
        JCBN_BIGK_VALUE = JCBN_BIGK_FUNCTION( rd, tpd,                                                        &
                                              CUR_VAL_BETA, CUR_DRV_BETA,                                     &
                                              CUR_R_LOCS(rd), R_SQUARE(rd), SIN_SQUARE(td), CSC_SQUARE(td),   &
                                              RSIN_SQUARE(td, rd), COTAN_VAL(td)                              )


        JCBN_kappa_Array = JCBN_kappa_FUNCTION_3D_ALL(  rd, tpd,                                        &
                                                        CUR_R_LOCS(rd), R_SQUARE(rd), R_CUBED(rd),      &
                                                        RSIN_SQUARE(td, rd),                            &
                                                        SIN_VAL(td), SIN_SQUARE(td), CSC_SQUARE(td),    &
                                                        COS_VAL(td), COTAN_VAL(td),                     &
                                                        CUR_VAL_BETA, CUR_DRV_BETA                      )



        JCBN_n_ARRAY(:) = CUR_DRV_ALPHAPSI( tpd, rd, : ) / ALPHAPSI_POWER(1)   &
                            - 7 * CUR_DRV_PSI( tpd, rd, : )/ PSI_POWER(1)






        CALL Calc_RHS_Terms( RHS_Terms,                                         &
                             re, te, pe,                                        &
                             td, pd, tpd, rd,                                   &
                             CUR_R_LOCS, R_SQUARE, RSIN_SQUARE, COTAN_VAL,      &
                             SIN_SQUARE, CSC_SQUARE,                            &
                             PSI_POWER, ALPHAPSI_POWER,                         &
                             CUR_VAL_BETA, CUR_DRV_BETA,                        &
                             JCBN_BIGK_VALUE, JCBN_n_Array, JCBN_Kappa_Array    )

        
        END DO ! rd loop
    END DO  ! td loop
END DO  ! pd loop


!$OMP END DO

!$OMP END PARALLEL


END SUBROUTINE Calc_GMRES_Terms




!+204+###########################################################################!
!                                                                                !
!                  CREATE_3D_RHS_VECTOR                                          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_GMRES_Equations( Eq, Coeff_Length, re, te, pe, DELTAR_OVERTWO )

COMPLEX(KIND = idp), DIMENSION(0:Coeff_Length-1), INTENT(OUT)   :: Eq
INTEGER, INTENT(IN)                                             :: Coeff_Length

INTEGER, INTENT(IN)                                                     ::  re, te, pe

REAL(KIND = idp), INTENT(IN)                                            ::  DELTAR_OVERTWO

INTEGER                                                                 ::  pd, td, rd, tpd,     &
                                                                            l, m, d,        &
                                                                            lm_loc, u

INTEGER                                                                 ::  Current_i_Location

COMPLEX(KIND = idp), DIMENSION(1:5)                                     ::  RHS_TMP
COMPLEX(KIND = idp)                                                     :: Test
COMPLEX(KIND = idp)                                                     ::  Common_Basis
REAL(KIND = idp)                                                        ::  Combined_Weights

COMPLEX(KIND = idp)                                                     ::  Inner, Middle




!$OMP PARALLEL DEFAULT(none)                                            &
!$OMP PRIVATE( pd, td, rd, d, l, m, tpd, u,                             &
!$OMP           lm_loc,                                                 &
!$OMP           Current_i_Location, RHS_TMP, TEST,                      &
!$OMP           Common_Basis, Combined_Weights                          )   &
!$OMP SHARED( re, te, pe,                                               &
!$OMP           DEGREE,                                                 &
!$OMP           R_SQUARE, Deltar_Overtwo, Int_R_Weights,                &
!$OMP           NUM_P_QUAD_POINTS, NUM_T_QUAD_POINTS, NUM_R_QUAD_POINTS,&
!$OMP           R_Int_Weights, TP_Int_Weights,                          &
!$OMP           ULM_LENGTH, LM_LENGTH,                                  &
!$OMP           Beta_DRV_Trace,                                         &
!$OMP           OneThird,                                               &
!$OMP           Ylm_CC_Values, Lagrange_Poly_Table,                     &
!$OMP           Block_RHS_Vector, RHS_Terms                             )

!$OMP DO SCHEDULE(dynamic), COLLAPSE(2)
DO d = 0,DEGREE


    DO lm_loc = 0,LM_LENGTH-1

        RHS_TMP = 0.0_idp


        DO rd = 1,NUM_R_QUAD_POINTS


            RHS_TMP(1) =  RHS_TMP(1)                                            &
                            + SUM( RHS_TERMS( :, rd, 1 )                        &
                                    * Ylm_CC_Values( :, lm_loc, te, pe)         &
                                    * TP_Int_Weights(:)                     )   &
                            * Lagrange_Poly_Table(d, rd, 0)                     &
                            * R_Int_Weights(rd)


            RHS_TMP(2) =  RHS_TMP(2)                                            &
                            + SUM( RHS_TERMS( :, rd, 2 )                        &
                                    * Ylm_CC_Values( :, lm_loc, te, pe)         &
                                    * TP_Int_Weights(:)                     )   &
                            * Lagrange_Poly_Table(d, rd, 0)                     &
                            * R_Int_Weights(rd)


            RHS_TMP(3) =  RHS_TMP(3)                                            &
                            + SUM( RHS_TERMS( :, rd, 3 )                        &
                                    * Ylm_CC_Values( :, lm_loc, te, pe)         &
                                    * TP_Int_Weights(:)                     )   &
                            * Lagrange_Poly_Table(d, rd, 0)                     &
                            * R_Int_Weights(rd)                                 &
                            + OneThird * SUM( Beta_DRV_Trace(:,rd)              &
                                * Ylm_CC_Values( :, lm_loc, te, pe)             &
                                * TP_Int_Weights(:)                         )   &
                            * Lagrange_Poly_Table(d, rd, 1 )/ DELTAR_OVERTWO    &
                            * R_Int_Weights(rd)



!            RHS_TMP(4) =  RHS_TMP(4)                                            &
!                            + SUM( RHS_TERMS( :, rd, 4 )                        &
!                                    * Ylm_CC_Values( :, lm_loc, te, pe)         &
!                                    * TP_Int_Weights(:)                     )   &
!                            * Lagrange_Poly_Table(d, rd, 0)                     &
!                            * R_Int_Weights(rd)                                 &
!                            + OneThird * SUM( Beta_DRV_Trace(:,rd)              &
!                                * Ylm_CC_dt_Values( :, lm_loc, te, pe)          &
!                                * TP_Int_Weights(:)                         )   &
!                            * Lagrange_Poly_Table(d, rd, 0 )                    &
!                            * R_Int_Weights(rd)/ R_SQUARE(rd)
!
!
!
!            RHS_TMP(5) =  RHS_TMP(5)                                                    &
!                            + SUM( RHS_TERMS( :, rd, 5 )                                &
!                                    * Ylm_CC_Values( :, lm_loc, te, pe)                 &
!                                    * TP_Int_Weights(:)                     )           &
!                            * Lagrange_Poly_Table(d, rd, 0)                             &
!                            * R_Int_Weights(rd)                                         &
!                            + OneThird * SUM( Beta_DRV_Trace(:,rd)                      &
!                                * Ylm_CC_dp_Values( :, lm_loc, te, pe)                  &
!                                * TP_Int_Weights(:)/(R_SQUARE(rd)*TP_SIN_SQUARE(:) ) )  &
!                            * Lagrange_Poly_Table(d, rd, 0 )                            &
!                            * R_Int_Weights(rd)


        END DO  ! rd Loop
        Current_i_Location = CFA_ALL_Matrix_Map(1, lm_loc, re, d)
        Eq(Current_i_Location) = Eq(Current_i_Location) + RHS_TMP(1)

        Current_i_Location = CFA_ALL_Matrix_Map(2, lm_loc, re, d)
        Eq(Current_i_Location) = Eq(Current_i_Location) + RHS_TMP(2)

        Current_i_Location = CFA_ALL_Matrix_Map(3, lm_loc, re, d)
        Eq(Current_i_Location) = Eq(Current_i_Location) + RHS_TMP(3)


        Current_i_Location = CFA_ALL_Matrix_Map(4, lm_loc, re, d)
        Eq(Current_i_Location) = Eq(Current_i_Location) + RHS_TMP(4)

        Current_i_Location = CFA_ALL_Matrix_Map(5, lm_loc, re, d)
        Eq(Current_i_Location) = Eq(Current_i_Location) + RHS_TMP(5)


    END DO  ! l Loop
END DO  ! d Loop
!$OMP END DO

!$OMP END PARALLEL



END SUBROUTINE Calc_GMRES_Equations








!+501+##############################################################################!
!                                                                                   !
!                  FINISH_3D_RHS_VECTOR                             !
!                                                                                   !
!###################################################################################!
SUBROUTINE FINISH_GMRES_Equations( Eq, Coeff_Length)

COMPLEX(KIND = idp), DIMENSION(0:Coeff_Length-1), INTENT(OUT)   :: Eq
INTEGER, INTENT(IN)                                             :: Coeff_Length

INTEGER                                                         ::  ui, l, m, re, d, dp, rd,k

INTEGER                                                         ::  Current_i_Location, &
                                                                    Current_j_Location

REAL(KIND = idp)                                                ::  TWOOVER_DELTAR, &
                                                                    L_Lp1

REAL(KIND = idp), DIMENSION(1:NUM_R_QUAD_POINTS)                ::  CUR_R_LOCS,     &
                                                                    R_SQUARE


REAL(KIND = idp), DIMENSION(0:DEGREE)                           :: Reusable_Values


INTEGER                                                         ::  Global_re
INTEGER                                                         ::  Block_re
INTEGER                                                         ::  iloc, jloc
INTEGER                                                         ::  MATVEC_LOC
REAL(KIND = idp)                                                ::  TMP_VALUE


INTEGER                                                         ::  Start_Here,     &
                                                                    End_Here,       &
                                                                    here


COMPLEX(KIND = idp), DIMENSION(0:SUBSHELL_PROB_DIM-1)           ::  TMP_VECTOR


INTEGER                                                         :: Start_Here_b, End_Here_b

INTEGER                                                         ::  ierr



INTEGER :: i,j


IF ( POSEIDON_COMM_PETSC .NE. MPI_COMM_NULL ) THEN

    Block_STF_MAT = 0.0_idp

    !$OMP PARALLEL DEFAULT(none)                                            &
    !$OMP PRIVATE( ui, l, m, re, d, dp, rd,                                 &
    !$OMP           CUR_R_LOCS, TWOOVER_DELTAR,                             &
    !$OMP           R_SQUARE, Reusable_Values, L_Lp1,                       &
    !$OMP           Global_re,                                              &
    !$OMP           Block_re,                                               &
    !$OMP           Start_Here, End_Here,                                   &
    !$OMP           Current_i_Location, Current_j_Location )                &
    !$OMP SHARED(   NUM_OFF_DIAGONALS,                                      &
    !$OMP           NUM_R_ELEMS_PER_SUBSHELL,                               &
    !$OMP           NUM_R_ELEMS_PER_SHELL,                                  &
    !$OMP           DEGREE, L_LIMIT,                                        &
    !$OMP           rlocs,                                                  &
    !$OMP           INT_R_WEIGHTS, INT_R_LOCATIONS, Lagrange_Poly_Table,    &
    !$OMP           myShell, myID_Poseidon, myID_SUBSHELL, ierr,            &
    !$OMP           myID_SHELL,                                             &
    !$OMP           Local_Length,                                           &
    !$OMP           ULM_LENGTH,                                             &
    !$OMP           M_VALUES,                                               &
    !$OMP           Matrix_Location,                                        &
    !$OMP           Block_STF_MAT, LPT_LPT         )





    DO re = 0,NUM_R_ELEMS_PER_SUBSHELL - 1


        Global_re = myShell*NUM_R_ELEMS_PER_SHELL                 &
                  + myID_Shell*NUM_R_ELEMS_PER_SUBSHELL      &
                  + re


        TWOOVER_DELTAR = 2.0_idp/(rlocs(Global_re  + 1) - rlocs(Global_re ))
        CUR_R_LOCS = (INT_R_LOCATIONS+1.0_idp)/TWOOVER_DELTAR + rlocs(Global_re )
        R_SQUARE = CUR_R_LOCS*CUR_R_LOCS

        !$OMP DO SCHEDULE(dynamic), COLLAPSE(2)
        DO d = 0,DEGREE

            DO ui = 1,NUM_CFA_VARS

                DO l = 0,L_LIMIT

                    L_Lp1 = REAL( l*(l+1), idp )


                    DO dp = 0,DEGREE

                        Reusable_Values(dp) = SUM ( ( -R_SQUARE(:) * LPT_LPT(:,d,dp,1,1)* TWOOVER_DELTAR* TWOOVER_DELTAR      &
                                                + L_Lp1 * LPT_LPT(:,d,dp,0,0)        )    &
                                                * INT_R_WEIGHTS(:)/TWOOVER_DELTAR     )
                    END DO
                    

                    DO m = -M_VALUES(l),M_VALUES(l)



                        Current_j_Location = Matrix_Location(ui, l, m, re, d)


                        DO dp = 0,DEGREE

                            Current_i_Location = NUM_OFF_DIAGONALS + (dp - d) * ULM_LENGTH


                            Block_STF_MAT(Current_i_Location, Current_j_Location )          &
                                    = Block_STF_MAT(Current_i_Location, Current_j_Location )      &
                                        + Reusable_Values(dp)



                        END DO ! dp Loop

                    END DO ! m Loop

                END DO ! l Loop

            END DO ! ui Loop

        END DO ! d Loop
        !$OMP END DO

    END DO ! re Loop

    !$OMP END PARALLEL

!    PRINT*,BLOCK_STF_MAT

    IF ( .FALSE. ) THEN
            CALL OUTPUT_LAPLACE_MATRIX
    END IF

    Start_Here = mySHELL*NUM_R_ELEMS_PER_BLOCK*DEGREE*ULM_LENGTH                &
                 + myID_SUBSHELL*NUM_R_ELEMS_PER_SUBSHELL*DEGREE*ULM_LENGTH
    End_Here = Start_Here + SUBSHELL_Prob_Dim - 1

    CALL ZGBMV('N',                                         &   ! TRANS
                SUBSHELL_PROB_DIM,                          &   ! M
                SUBSHELL_PROB_DIM,                          &   ! N
                NUM_OFF_DIAGONALS,                          &   ! KL
                NUM_OFF_DIAGONALS,                          &   ! KU
                CMPLX(1.0_idp,0.0_idp,KIND = idp),          &   ! ALPHA
                Block_STF_MAT(:,:),                         &   ! A
                2*NUM_OFF_DIAGONALS+1,                      &   ! LDA
                -Coefficient_Vector(Start_Here:End_Here),   &   ! X
                1,                                          &   ! INCX
                CMPLX(0.0_idp,0.0_idp,KIND = idp),          &   ! BETA
                TMP_VECTOR(:),                              &   ! Y
                1                                           )   ! INCY





    Start_Here = myID_SubShell*NUM_R_ELEMS_PER_SUBSHELL*DEGREE*ULM_LENGTH
    END_HERE = Start_HERE + SUBSHELL_PROB_DIM -1

!    PRINT*,TMP_VECTOR

!    CALL ANALYZE_3D_RHS_VECTOR( TMP_VECTOR, Block_RHS_Vector)
    IF ( OUTPUT_RHS_VECTOR_FLAG == 1 ) THEN

        CALL OUTPUT_RHS_VECTOR_Parts(TMP_VECTOR, Eq)

    END IF


    Eq(Start_Here:End_Here) = Eq(Start_Here:End_Here)   &
                            + TMP_VECTOR(0:SUBSHELL_PROB_DIM-1)

END IF

END SUBROUTINE FINISH_GMRES_Equations








!+601+###########################################################################!
!                                                                                !
!           Allocate_Master_Build_Variables                                      !
!                                                                                !
!################################################################################!
SUBROUTINE Allocate_GMRES_Eq_Variables()


ALLOCATE( RR_Factor( 1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE )    )
ALLOCATE( RDR_Factor( 1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE )    )
ALLOCATE( DRR_Factor( 1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE )    )
ALLOCATE( DRDR_Factor( 1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE )    )

ALLOCATE( TP_TP_Factor( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1) )
ALLOCATE( dTP_TP_Factor( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1) )
ALLOCATE( TdP_TP_Factor( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1) )
ALLOCATE( TP_dTP_Factor( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1) )
ALLOCATE( TP_TdP_Factor( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1) )
ALLOCATE( dTP_dTP_Factor( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1) )
ALLOCATE( dTP_TdP_Factor( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1) )
ALLOCATE( TdP_dTP_Factor( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1) )
ALLOCATE( TdP_TdP_Factor( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1) )

ALLOCATE( R_Int_Weights( 1:NUM_R_QUAD_POINTS ) )
ALLOCATE( TP_Int_Weights( 1:NUM_TP_QUAD_POINTS) )

ALLOCATE( CUR_R_LOCS(1:NUM_R_QUAD_POINTS) )
ALLOCATE( CUR_T_LOCS(1:NUM_T_QUAD_POINTS) )
ALLOCATE( CUR_P_LOCS(1:NUM_P_QUAD_POINTS) )


ALLOCATE( R_SQUARE(1:NUM_R_QUAD_POINTS) )
ALLOCATE( R_CUBED(1:NUM_R_QUAD_POINTS) )

ALLOCATE( SIN_VAL( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( SIN_SQUARE( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( TP_SIN_SQUARE( 1:NUM_TP_QUAD_POINTS ) )
ALLOCATE( COS_VAL( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( COS_SQUARE( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( CSC_VAL( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( CSC_SQUARE( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( COTAN_VAL( 1:NUM_T_QUAD_POINTS ) )

ALLOCATE( RSIN_SQUARE( 1:NUM_T_QUAD_POINTS, 1:NUM_R_QUAD_POINTS ) )

ALLOCATE( PHI_EXP( 1:NUM_P_QUAD_POINTS ) )
ALLOCATE( PHI_TWOEXP( 1:NUM_P_QUAD_POINTS ) )

ALLOCATE( CUR_VAL_PSI(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS) )
ALLOCATE( CUR_VAL_ALPHAPSI(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS) )
ALLOCATE( CUR_VAL_BETA(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:3) )

ALLOCATE( CUR_DRV_PSI(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:3) )
ALLOCATE( CUR_DRV_ALPHAPSI(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:3) )
ALLOCATE( CUR_DRV_BETA(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:3, 1:3) )

ALLOCATE( Beta_DRV_Trace(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS) )

ALLOCATE( RHS_Terms( 1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:5) )

END SUBROUTINE Allocate_GMRES_Eq_Variables



!+602+###########################################################################!
!                                                                                !
!           Deallocate_Master_Build_Variables                                    !
!                                                                                !
!################################################################################!
SUBROUTINE Deallocate_GMRES_Eq_Variables()

DEALLOCATE( RR_Factor )
DEALLOCATE( DRR_Factor )
DEALLOCATE( RDR_Factor )
DEALLOCATE( DRDR_Factor )

DEALLOCATE( TP_TP_Factor )
DEALLOCATE( dTP_TP_Factor )
DEALLOCATE( TdP_TP_Factor )
DEALLOCATE( TP_dTP_Factor )
DEALLOCATE( TP_TdP_Factor )
DEALLOCATE( dTP_dTP_Factor )
DEALLOCATE( dTP_TdP_Factor )
DEALLOCATE( TdP_dTP_Factor )
DEALLOCATE( TdP_TdP_Factor )

DEALLOCATE( R_Int_Weights )
DEALLOCATE( TP_Int_Weights )

DEALLOCATE( CUR_R_LOCS )
DEALLOCATE( CUR_T_LOCS )
DEALLOCATE( CUR_P_LOCS )

DEALLOCATE( R_SQUARE )
DEALLOCATE( R_CUBED )

DEALLOCATE( SIN_VAL )
DEALLOCATE( SIN_SQUARE )
DEALLOCATE( TP_SIN_SQUARE )
DEALLOCATE( COS_VAL )
DEALLOCATE( COS_SQUARE )
DEALLOCATE( CSC_VAL )
DEALLOCATE( CSC_SQUARE )
DEALLOCATE( COTAN_VAL )

DEALLOCATE( RSIN_SQUARE )

DEALLOCATE( PHI_EXP )
DEALLOCATE( PHI_TWOEXP )

DEALLOCATE( CUR_VAL_PSI )
DEALLOCATE( CUR_VAL_ALPHAPSI )
DEALLOCATE( CUR_VAL_BETA )

DEALLOCATE( CUR_DRV_PSI )
DEALLOCATE( CUR_DRV_ALPHAPSI )
DEALLOCATE( CUR_DRV_BETA )

DEALLOCATE( Beta_DRV_Trace )

DEALLOCATE( RHS_Terms )



END SUBROUTINE Deallocate_GMRES_Eq_Variables






END MODULE GMRES_Master_Build_Module
