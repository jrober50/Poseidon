   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE CFA_3D_Master_Build_Module                                                   !##!
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
!##!    +101+   CFA_3D_Master_Build                                                 !##!
!##!                                                                                !##!
!##!    +201+   CREATE_3D_NONLAPLACIAN_SOE                                          !##!
!##!    +202+   Calc_3D_Current_Values                                              !##!
!##!    +203+   Calc_3D_SubJcbn_Terms                                               !##!
!##!    +204+   CREATE_3D_RHS_VECTOR                                                !##!
!##!    +205+   CREATE_3D_JCBN_MATRIX                                               !##!
!##!                                                                                !##!
!##!    +301+   REDUCE_3D_NONLAPLACIAN_SOE                                          !##!
!##!                                                                                !##!
!##!    +401+   FINISH_3D_JACOBIAN_MATRIX                                           !##!
!##!                                                                                !##!
!##!    +501+   FINISH_3D_RHS_VECTOR                                                !##!
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

USE MPI

USE OMP_LIB

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
                                        nPROCS_SHELL,               &
                                        myID_Poseidon,              &
                                        myID_Shell,                 &
                                        myID_SubShell,              &
                                        myID_PETSc,                 &
                                        PROB_DIM,                   &
                                        Block_PROB_DIM,             &
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






USE Poseidon_Additional_Functions_Module, &
                                ONLY :  Lagrange_Poly,              &
                                        Spherical_Harmonic,         &
                                        Initialize_LGL_Quadrature


USE Jacobian_Internal_Functions_Module, &
                                ONLY :  JCBN_kappa_FUNCTION_3D_ALL,     &
                                        JCBN_BIGK_FUNCTION


USE IO_Functions_Module, &
                                ONLY :  Clock_In,                       &
                                        OUTPUT_RHS_VECTOR_Parts

USE Poseidon_BC_Module, &
                                ONLY :  CFA_3D_Apply_BCs_Part1,         &
                                        CFA_3D_Apply_BCs_Part2

USE Poseidon_Mapping_Functions_Module, &
                                ONLY :  Map_To_X_Space,                 &
                                        CFA_ALL_Matrix_Map

USE SubJacobian_Functions_Module_3D, &
                                ONLY :  Calc_EQ1_SubJacobian,           &
                                        Calc_EQ2_SubJacobian,           &
                                        Calc_EQ3_SubJacobian,           &
                                        Calc_EQ4_SubJacobian,           &
                                        Calc_EQ5_SubJacobian,           &
                                        Calc_RHS_Terms

USE Poseidon_Preconditioner_Module, &
                                ONLY :  Jacobi_Type_A_PC,               &
                                        Jacobi_Type_B_PC

IMPLICIT NONE


!
!   Module Specific Variables
!

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

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: SubJacobian_EQ1_Term
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: SubJacobian_EQ2_Term
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: SubJacobian_EQ3_Term
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: SubJacobian_EQ4_Term
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: SubJacobian_EQ5_Term


CONTAINS

!+101+###########################################################################!
!                                                                                !
!           CFA_3D_Master_Build                                                  !
!                                                                                !
!################################################################################!
SUBROUTINE CFA_3D_Master_Build()



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
CALL Allocate_Master_Build_Variables()



!*!
!*! Alter the Coefficient Vector to Reflect Boundary Conditions
!*!
timeb = MPI_Wtime()
CALL CFA_3D_Apply_BCs_Part1()
timec = MPI_Wtime()
CALL Clock_In(timec-timeb, 4)


!*!
!*! Create the Non-Laplacian Contributions to the Linear System
!*!
timeb = MPI_Wtime()
CALL CREATE_3D_NONLAPLACIAN_SOE()
timea = MPI_Wtime()
CALL Clock_In(timea-timeb, 9)






!*!
!*! Reduce the Non-Laplacian Contributions to Processes involved in the solve
!*!
timeb = MPI_Wtime()
CALL REDUCE_3D_NONLAPLACIAN_SOE()
timea = MPI_Wtime()
CALL Clock_In(timea-timeb, 10)



!*!
!*! Add the Laplacian contribution to the Jacobian Matrix
!*!
timeb = MPI_Wtime()
CALL FINISH_3D_JACOBIAN_MATRIX()
timea = MPI_Wtime()
CALL Clock_In(timea-timeb, 11)


!*!
!*! Add the Laplacian contribution to the Residual Vector
!*!
timeb = MPI_Wtime()
CALL FINISH_3D_RHS_VECTOR()
timea = MPI_Wtime()
CALL Clock_In(timea-timeb, 12)



!CALL Jacobi_Type_A_PC()



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
CALL Deallocate_Master_Build_Variables()


END SUBROUTINE CFA_3D_Master_Build
















!+201+###############################################################################!
!                                                                                    !
!                                  CREATE_3D_NONLAPLACIAN_SOE                        !
!                                                                                    !
!####################################################################################!
SUBROUTINE CREATE_3D_NONLAPLACIAN_SOE()

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
                                        DELTAR_OVERTWO, DELTAT_OVERTWO, DELTAP_OVERTWO  )


            timeb = MPI_Wtime()




            !*!
            !*!  Calculate the Sub-Jacobian and RHS Terms
            !*!
            CALL Calc_3D_SubJcbn_Terms( Local_re, Local_te, Local_pe )

            timec = MPI_Wtime()



            !*!
            !*! Create the Residual Vector ( Sans Laplacian Contribution )
            !*!
            CALL CREATE_3D_RHS_VECTOR(  Local_re, Local_te, Local_pe, DELTAR_OVERTWO )

            timed = MPI_Wtime()


            

            !*!
            !*! Create Jacobian Matrix ( Sans Laplacian Contribution )
            !*!
            CALL CREATE_3D_JCBN_MATRIX( Local_re, Local_te, Local_pe,               &
                                        Global_re , Global_te  , Global_pe,         &
                                        TWOOVER_DELTAR                              )



            timee = MPI_Wtime()



            !*!
            !*! Update Timer Values
            !*!
            time_CurVals = time_CurVals + timeb - timea
            time_SJT = time_SJT + timec - timeb
            time_RHS = time_RHS + timed - timec
            time_JCBNM = time_JCBNM + timee - timed



        END DO  ! pe Loop
    END DO  ! te Loop
END DO  ! re Loop


CALL Clock_In(time_CurVals, 5)
CALL Clock_In(time_SJT, 6)
CALL Clock_In(time_RHS, 7)
CALL Clock_In(time_JCBNM, 8)




END SUBROUTINE CREATE_3D_NONLAPLACIAN_SOE






!+202+###########################################################################!
!                                                                                !
!                  Calc_Current_Values          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_3D_Current_Values( re, te, pe,                                  &
                                    DELTAR_OVERTWO,                             &
                                    DELTAT_OVERTWO,                             &
                                    DELTAP_OVERTWO                              )

INTEGER, INTENT(IN)                                             ::  re, te, pe


REAL(KIND = idp), INTENT(IN)                                    ::  DELTAR_OVERTWO,     &
                                                                    DELTAT_OVERTWO,     &
                                                                    DELTAP_OVERTWO



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





!$OMP PARALLEL DEFAULT(none)                                            &
!$OMP PRIVATE( tpd, rd, l, m, d, ui,                                    &
!$OMP           TMP_U_Value, Tmp_U_R_DRV_Value, Tmp_U_T_DRV_Value,      &
!$OMP           Tmp_U_P_DRV_Value, Tmp_U_RR_DDRV_Value,                 &
!$OMP           local_coefficients,                                     &
!$OMP           lm_loc, Current_Location                            )   &
!$OMP SHARED( re, te, pe,                                               &
!$OMP           DEGREE,                                                 &
!$OMP           NUM_P_QUAD_POINTS, NUM_T_QUAD_POINTS, NUM_R_QUAD_POINTS,&
!$OMP           CUR_VAL_PSI, CUR_VAL_ALPHAPSI, CUR_VAL_BETA,            &
!$OMP           CUR_DRV_PSI, CUR_DRV_ALPHAPSI, CUR_DRV_BETA,            &
!$OMP           CUR_DDRV_BETA,                                          &
!$OMP           CUR_R_LOCS, CUR_T_LOCS, CUR_P_LOCS,                     &
!$OMP           R_SQUARE,                                               &
!$OMP           SIN_VAL,                                                &
!$OMP           myID_Poseidon,                                          &
!$OMP           DELTAR_OVERTWO, DELTAT_OVERTWO, DELTAP_OVERTWO,         &
!$OMP           Lagrange_Poly_Table,                                    &
!$OMP           Ylm_Values, Ylm_dt_Values, Ylm_dp_Values,               &
!$OMP           Coefficient_Vector,                                     &
!$OMP           Matrix_Location,                                        &
!$OMP           LM_Location,                                            &
!$OMP           LM_Length, ULM_LENGTH                           )




!$OMP DO SCHEDULE(dynamic), COLLAPSE(2)
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
                                        + SUM( Coefficient_Vector( Here:There )     &
                                        * Ylm_Values( :, tpd, te, pe )       )      &
                                        * Lagrange_Poly_Table( d, rd, 0 )

                TMP_U_R_DRV_Value(ui)   = TMP_U_R_DRV_Value(ui)                     &
                                        + SUM( Coefficient_Vector( Here:There )     &
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
                                        + SUM( Coefficient_Vector( Here:There )     &
                                        * Ylm_dt_Values( :, tpd, te, pe)     )       &
                                        * Lagrange_Poly_Table( d, rd, 0)

                TMP_U_P_DRV_Value(ui)   = TMP_U_P_DRV_Value(ui)                     &
                                        + SUM( Coefficient_Vector( Here:There )     &
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
        
!        PRINT*,RE, rd, REAL(Tmp_U_R_DRV_Value(3), KIND = idp)
        

    END DO  !   tpd Loop

END DO  !   rd Loop
!$OMP END DO

!$OMP END PARALLEL








END SUBROUTINE Calc_3D_Current_Values











!+203+###########################################################################!
!                                                                                !
!                  Calc_3D_SubJcbn_Terms                                         !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_3D_SubJcbn_Terms( re, te, pe )




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
!$OMP           CUR_DDRV_BETA,                                          &
!$OMP           CUR_R_LOCS, CUR_T_LOCS, CUR_P_LOCS,                     &
!$OMP           RSIN_SQUARE, R_SQUARE, R_CUBED, R_INVERSE,              &
!$OMP           SIN_VAL, SIN_SQUARE, COTAN_VAL, COS_VAL,                &
!$OMP           CSC_VAL, CSC_SQUARE,                                    &
!$OMP           GR_Source_Scalar,                                       &
!$OMP           NUM_OFF_DIAGONALS,                                      &
!$OMP           RHS_TERMS, SUBJCBN_PSI_TERMS, SUBJCBN_ALPHAPSI_TERMS,   &
!$OMP           SUBJCBN_BETA1_TERMS, SUBJCBN_BETA2_TERMS,               &
!$OMP           SUBJCBN_BETA3_TERMS,                                    &
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

        
        !
        !   Equation 1 ( Conformal Factor ) Subjacobian Terms
        !

        CALL Calc_EQ1_SubJacobian( SubJacobian_EQ1_Term,                              &
                                   re, te, pe,                                        &
                                   td, pd, tpd, rd,                                   &
                                   CUR_R_LOCS, R_SQUARE, RSIN_SQUARE, COTAN_VAL,      &
                                   SIN_SQUARE, CSC_SQUARE,                            &
                                   PSI_POWER, ALPHAPSI_POWER,                         &
                                   CUR_VAL_BETA, CUR_DRV_BETA,                        &
                                   JCBN_BIGK_VALUE                                    )


        CALL Calc_EQ2_SubJacobian( SubJacobian_EQ2_Term,                              &
                                   SubJacobian_EQ1_Term,                              &
                                   re, te, pe,                                        &
                                   td, pd, tpd, rd,                                   &
                                   CUR_R_LOCS, R_SQUARE, RSIN_SQUARE, COTAN_VAL,      &
                                   SIN_SQUARE, CSC_SQUARE,                            &
                                   PSI_POWER, ALPHAPSI_POWER,                         &
                                   JCBN_BIGK_VALUE                                    )

!        PRINT*,rd,td,pd
!        PRINT*,SubJacobian_EQ1_Term(tpd,rd,1:14)
!        PRINT*,SubJacobian_EQ2_Term(tpd,rd,1:14)


        CALL Calc_EQ3_SubJacobian( SubJacobian_EQ3_Term,                              &
                                   re, te, pe,                                        &
                                   td, pd, tpd, rd,                                   &
                                   CUR_R_LOCS, R_SQUARE, RSIN_SQUARE, COTAN_VAL,      &
                                   SIN_SQUARE, CSC_SQUARE,                            &
                                   PSI_POWER, ALPHAPSI_POWER,                         &
                                   CUR_DRV_PSI, CUR_DRV_ALPHAPSI,                     &
                                   JCBN_BIGK_VALUE, JCBN_n_Array, JCBN_Kappa_Array    )


        CALL Calc_EQ4_SubJacobian( SubJacobian_EQ4_Term,                              &
                                   re, te, pe,                                        &
                                   td, pd, tpd, rd,                                   &
                                   CUR_R_LOCS, R_SQUARE, R_CUBED, RSIN_SQUARE,        &
                                   COTAN_VAL, SIN_SQUARE, CSC_SQUARE,                 &
                                   PSI_POWER, ALPHAPSI_POWER,                         &
                                   CUR_DRV_PSI, CUR_DRV_ALPHAPSI,                     &
                                   JCBN_BIGK_VALUE, JCBN_n_Array, JCBN_Kappa_Array    )


        CALL Calc_EQ5_SubJacobian( SubJacobian_EQ5_Term,                              &
                                   re, te, pe,                                        &
                                   td, pd, tpd, rd,                                   &
                                   CUR_R_LOCS, R_SQUARE, R_CUBED, RSIN_SQUARE,        &
                                   COTAN_VAL, SIN_SQUARE, CSC_SQUARE,                 &
                                   PSI_POWER, ALPHAPSI_POWER,                         &
                                   CUR_DRV_PSI, CUR_DRV_ALPHAPSI,                     &
                                   JCBN_BIGK_VALUE, JCBN_n_Array, JCBN_Kappa_Array    )

!        PRINT*,rd,td,pd
!        PRINT*,SubJacobian_EQ3_Term(tpd,rd,1:20)
!        PRINT*,""
!        PRINT*," "
!        PRINT*,SubJacobian_EQ4_Term(tpd,rd,1:20)
!        PRINT*,""
!        PRINT*," "
!        PRINT*,SubJacobian_EQ5_Term(tpd,rd,1:20)
!        PRINT*,""
!        PRINT*," "

        END DO ! rd loop
    END DO  ! td loop
END DO  ! pd loop


!$OMP END DO

!$OMP END PARALLEL





END SUBROUTINE Calc_3D_SubJcbn_Terms

















!+204+###########################################################################!
!                                                                                !
!                  CREATE_3D_RHS_VECTOR                                          !
!                                                                                !
!################################################################################!
SUBROUTINE CREATE_3D_RHS_VECTOR( re, te, pe, DELTAR_OVERTWO )



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
!$OMP           NUM_CFA_VARS,                                           &
!$OMP           R_Int_Weights, TP_Int_Weights,                          &
!$OMP           ULM_LENGTH, LM_LENGTH,                                  &
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



            RHS_TMP(4) =  RHS_TMP(4)                                            &
                            + SUM( RHS_TERMS( :, rd, 4 )                        &
                                    * Ylm_CC_Values( :, lm_loc, te, pe)         &
                                    * TP_Int_Weights(:)                     )   &
                            * Lagrange_Poly_Table(d, rd, 0)                     &
                            * R_Int_Weights(rd)                                 &
                            + OneThird * SUM( Beta_DRV_Trace(:,rd)              &
                                * Ylm_CC_dt_Values( :, lm_loc, te, pe)          &
                                * TP_Int_Weights(:)                         )   &
                            * Lagrange_Poly_Table(d, rd, 0 )                    &
                            * R_Int_Weights(rd)/ R_SQUARE(rd)



            RHS_TMP(5) =  RHS_TMP(5)                                                    &
                            + SUM( RHS_TERMS( :, rd, 5 )                                &
                                    * Ylm_CC_Values( :, lm_loc, te, pe)                 &
                                    * TP_Int_Weights(:)                     )           &
                            * Lagrange_Poly_Table(d, rd, 0)                             &
                            * R_Int_Weights(rd)                                         &
                            + OneThird * SUM( Beta_DRV_Trace(:,rd)                      &
                                * Ylm_CC_dp_Values( :, lm_loc, te, pe)                  &
                                * TP_Int_Weights(:)/(R_SQUARE(rd)*TP_SIN_SQUARE(:) ) )  &
                            * Lagrange_Poly_Table(d, rd, 0 )                            &
                            * R_Int_Weights(rd)


        END DO  ! rd Loop
        Current_i_Location = CFA_ALL_Matrix_Map(1, lm_loc, re, d)
        Block_RHS_Vector(Current_i_Location)                             &
            = Block_RHS_Vector(Current_i_Location)                       &
            + RHS_TMP(1)

        Current_i_Location = CFA_ALL_Matrix_Map(2, lm_loc, re, d)
        Block_RHS_Vector(Current_i_Location)                             &
            = Block_RHS_Vector(Current_i_Location)                       &
            + RHS_TMP(2)

        Current_i_Location = CFA_ALL_Matrix_Map(3, lm_loc, re, d)
        Block_RHS_Vector(Current_i_Location)                             &
            = Block_RHS_Vector(Current_i_Location)                       &
            + RHS_TMP(3)


        Current_i_Location = CFA_ALL_Matrix_Map(4, lm_loc, re, d)
        Block_RHS_Vector(Current_i_Location)                             &
            = Block_RHS_Vector(Current_i_Location)                       &
            + RHS_TMP(4)

        Current_i_Location = CFA_ALL_Matrix_Map(5, lm_loc, re, d)
        Block_RHS_Vector(Current_i_Location)                             &
            = Block_RHS_Vector(Current_i_Location)                       &
            + RHS_TMP(5)


    END DO  ! l Loop
END DO  ! d Loop
!$OMP END DO

!$OMP END PARALLEL



END SUBROUTINE CREATE_3D_RHS_VECTOR






















!+205+###########################################################################!
!                                                                                !
!                  CREATE_3D_JCBN_MATRIX          !
!                                                                                !
!################################################################################!
SUBROUTINE CREATE_3D_JCBN_MATRIX(   Local_re, Local_te, Local_pe,               &
                                    Global_re , Global_te  , Global_pe,         &
                                    TWOOVER_DELTAR                              )



INTEGER, INTENT(IN)                                                     ::  Local_re, Local_te, Local_pe,       &
                                                                            Global_re, Global_te, Global_pe


REAL(KIND = idp), INTENT(IN)                                            ::  TWOOVER_DELTAR





INTEGER                                                                 ::  pd, td, rd,     &
                                                                            l, m, d,        &
                                                                            lp, mp, dp



INTEGER                                                                 ::  Current_i_Location,         &
                                                                            Current_j_Location



INTEGER                                                                 ::  F, u
INTEGER                                                                 ::  lm_loc, lpmp_loc
INTEGER                                                                 ::  MATVEC_LOC
INTEGER                                                                 ::  i_loc, j_loc

INTEGER                                                                 ::  Start, Finish

REAL(KIND = idp), DIMENSION(1:NUM_R_QUAD_POINTS)                        ::  Common_Term_A,  &
                                                                            Common_Term_B,  &
                                                                            Common_Term_C


COMPLEX(KIND = idp), DIMENSION(0:ELEM_PROB_DIM-1)                       ::  Jacobian_Buffer
COMPLEX(KIND = idp), DIMENSION(1:25)                                    ::  Jacobian_Terms


REAL(KIND = idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:2)                   ::  Current_Lag_Polys




INTEGER                                                                 ::  i, ierr


! dp, F, and lpmp_loc choose the row
DO dp = 0,DEGREE
    DO F = 1,NUM_CFA_VARS
        DO lpmp_loc = 0,LM_LENGTH-1

            CALL CALC_JACOBIAN_A( Jacobian_Buffer, dp, lpmp_loc, F )


            ! Put Buffer into Matrix
            ! row location
            i_loc = dp*ULM_LENGTH + (F-1)*LM_LENGTH + lpmp_loc

            Start = i_loc*ELEM_PROB_DIM
            Finish = Start + ELEM_PROB_DIM-1
            BLOCK_ELEM_STF_MATVEC(Start:Finish,Local_re) = BLOCK_ELEM_STF_MATVEC(Start:Finish,Local_re)     &
                                                         + Jacobian_Buffer(:)

!            STOP
        END DO ! lpmp_loc Loop
    END DO ! F Loop
END DO ! dp Loop




END SUBROUTINE CREATE_3D_JCBN_MATRIX





















!+301+##############################################################################!
                                                                                   !
!                  REDUCE_3D_NONLAPLACIAN_SOE                             !
!                                                                                   !
!###################################################################################!
SUBROUTINE REDUCE_3D_NONLAPLACIAN_SOE()

INTEGER                                     ::  SUBSHELL
INTEGER                                     ::  SEND_LENGTH
INTEGER                                     ::  ierr

INTEGER                                     ::  Start_Here_RE,     &
                                                End_Here_RE

INTEGER                                     ::  Start_Here_ND,     &
                                                MidA_Here_ND,      &
                                                End_Here_ND

COMPLEX(KIND=idp),DIMENSION(0:SUBSHELL_PROB_DIM-1)    :: RHS_Receive_Buffer


INTEGER    :: i


IF ( nPROCS_SHELL .NE. 0 ) THEN

    IF ( NUM_SUBSHELLS_PER_SHELL == 1 ) THEN

        SEND_LENGTH = ELEM_PROB_DIM_SQR*NUM_R_ELEMS_PER_BLOCK

        IF ( myID_Shell == 0 ) THEN

            CALL MPI_REDUCE( MPI_IN_PLACE, BLOCK_ELEM_STF_MATVEC, SEND_LENGTH,          &
                             MPI_DOUBLE_COMPLEX, MPI_SUM, 0, POSEIDON_COMM_SHELL, ierr )


            CALL MPI_REDUCE( MPI_IN_PLACE, Block_RHS_VECTOR, Block_PROB_DIM,            &
                             MPI_DOUBLE_COMPLEX, MPI_SUM, 0, POSEIDON_COMM_SHELL, ierr )


        ELSE

            CALL MPI_REDUCE( BLOCK_ELEM_STF_MATVEC, BLOCK_ELEM_STF_MATVEC, SEND_LENGTH,    &
                             MPI_DOUBLE_COMPLEX, MPI_SUM, 0, POSEIDON_COMM_SHELL, ierr )


            CALL MPI_REDUCE( Block_RHS_VECTOR, Block_RHS_VECTOR, Block_PROB_DIM,           &
                             MPI_DOUBLE_COMPLEX, MPI_SUM, 0, POSEIDON_COMM_SHELL, ierr  )


        END IF


    ELSE ! NUM_SUBSHELLS_PER_SHELL > 1

        SEND_LENGTH = ELEM_PROB_DIM_SQR*NUM_R_ELEMS_PER_SUBSHELL
        DO SUBSHELL = 0,NUM_SUBSHELLS_PER_SHELL-1

            Start_Here_RE = SUBSHELL*NUM_R_ELEMS_PER_SUBSHELL 
            End_Here_RE = Start_Here_RE + NUM_R_ELEMS_PER_SUBSHELL



            Start_Here_ND = Start_Here_RE * DEGREE*ULM_LENGTH
            End_Here_ND = End_Here_RE * SUBSHELL_PROB_DIM

            



            IF ( myID_SUBSHELL == SUBSHELL ) THEN

                CALL MPI_REDUCE( MPI_IN_PLACE,                                          &
                                 BLOCK_ELEM_STF_MATVEC(0,Start_Here_RE),    &
                                 SEND_LENGTH,                                           &
                                 MPI_DOUBLE_COMPLEX,                                    &
                                 MPI_SUM,                                               &
                                 SUBSHELL,                                              &
                                 POSEIDON_COMM_SHELL,                                   &
                                 ierr                                                   )


                CALL MPI_REDUCE( Block_RHS_VECTOR(Start_Here_ND),                       &
                                 RHS_Receive_Buffer,                                    &
                                 SUBSHELL_PROB_DIM,                                     &
                                 MPI_DOUBLE_COMPLEX,                                    &
                                 MPI_SUM,                                               &
                                 SUBSHELL,                                              &
                                 POSEIDON_COMM_SHELL,                                   &
                                 ierr                                                   )


            ELSE

                CALL MPI_REDUCE( BLOCK_ELEM_STF_MATVEC(0,Start_Here_RE),    &
                                 BLOCK_ELEM_STF_MATVEC(0,Start_Here_RE),    &
                                 SEND_LENGTH,                                           &
                                 MPI_DOUBLE_COMPLEX,                                    &
                                 MPI_SUM,                                               &
                                 SUBSHELL,                                              &
                                 POSEIDON_COMM_SHELL,                                   &
                                 ierr                                                   )


                CALL MPI_REDUCE( Block_RHS_VECTOR(Start_Here_ND),           &
                                 Block_RHS_VECTOR(Start_Here_ND),           &
                                 SUBSHELL_PROB_DIM,                                     &
                                 MPI_DOUBLE_COMPLEX,                                    &
                                 MPI_SUM,                                               &
                                 SUBSHELL,                                              &
                                 POSEIDON_COMM_SHELL,                                   &
                                 ierr                                                   )


            END IF



        END DO ! SUBSHELL Loop



        IF (myID_SUBSHELL == NUM_SUBSHELLS_PER_SHELL - 1) THEN

            Start_Here_ND = myID_SUBSHELL*NUM_R_ELEMS_PER_SUBSHELL * DEGREE*ULM_LENGTH
            End_Here_ND = Start_Here_ND + SUBSHELL_PROB_DIM  - 1

            BLOCK_RHS_VECTOR(Start_Here_ND:End_Here_ND) = RHS_Receive_Buffer

        ELSE IF ( myID_SUBSHELL .NE. -1 ) THEN

            START_HERE_ND = myID_SUBSHELL*NUM_R_ELEMS_PER_SUBSHELL * DEGREE*ULM_LENGTH
            MIDA_HERE_ND = Start_Here_ND + SUBSHELL_PROB_DIM - ULM_LENGTH - 1
            END_HERE_ND = Start_Here_ND + SUBSHELL_PROB_DIM - 1

            BLOCK_RHS_VECTOR(Start_Here_ND:MidA_Here_ND) = RHS_Receive_Buffer(0:SUBSHELL_PROB_DIM - ULM_LENGTH - 1)
            BLOCK_RHS_VECTOR(MidA_HERE_ND+1:END_HERE_ND) = 0.0_idp


        END IF

    END IF 
END IF






END SUBROUTINE REDUCE_3D_NONLAPLACIAN_SOE







!+401+##############################################################################!
!                                                                                   !
!                  FINISH_3D_JACOBIAN_MATRIX                             !
!                                                                                   !
!###################################################################################!
SUBROUTINE FINISH_3D_JACOBIAN_MATRIX()



INTEGER                                             ::  ui, l, m, re, d, dp, rd, u

INTEGER                                             ::  Current_i_Location, &
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


REAL(KIND = idp)                                                :: TMP_VALUE

INTEGER :: i, ierr









IF ( myID_SubShell .NE. -1 ) THEN

    !$OMP PARALLEL DEFAULT(none)                                            &
    !$OMP PRIVATE( ui, l, m, re, d, dp, rd, u,                              &
    !$OMP           CUR_R_LOCS, TWOOVER_DELTAR,                             &
    !$OMP           R_SQUARE, Reusable_Values, L_Lp1,                       &
    !$OMP           Global_re, Block_re,                                    &
    !$OMP           iloc, jloc, MATVEC_LOC,                                 &
    !$OMP           Current_i_Location, Current_j_Location )                &
    !$OMP SHARED(   NUM_OFF_DIAGONALS,                                      &
    !$OMP           NUM_R_ELEMS_PER_SUBSHELL,                               &
    !$OMP           NUM_R_ELEMS_PER_SHELL,                                  &
    !$OMP           DEGREE, L_LIMIT,                                        &
    !$OMP           ELEM_PROB_DIM,                                          &
    !$OMP           BLOCK_ELEM_STF_MATVEC,                                  &
    !$OMP           rlocs,                                                  &
    !$OMP           INT_R_WEIGHTS, INT_R_LOCATIONS, Lagrange_Poly_Table,    &
    !$OMP           myShell, myID_Poseidon,                                 &
    !$OMP           myID_SubShell,                                          &
    !$OMP           M_VALUES,                                               &
    !$OMP           LPT_LPT, ULM_Length         )







    !$OMP DO SCHEDULE(dynamic),COLLAPSE(3)
    DO re = 0,NUM_R_ELEMS_PER_SUBSHELL - 1

        DO d = 0,DEGREE

            DO u = 1,NUM_CFA_VARS

                DO l = 0,L_LIMIT

                    Global_re = myShell*NUM_R_ELEMS_PER_SHELL                 &
                              + myID_SubShell*NUM_R_ELEMS_PER_SUBSHELL      &
                              + re

                    Block_re = myID_SubShell*NUM_R_ELEMS_PER_SUBSHELL       &
                             + re

                    TWOOVER_DELTAR = 2.0_idp/(rlocs(Global_re  + 1) - rlocs(Global_re ))
                    CUR_R_LOCS = (INT_R_LOCATIONS+1.0_idp)/TWOOVER_DELTAR + rlocs(Global_re )
                    R_SQUARE = CUR_R_LOCS*CUR_R_LOCS


                    L_Lp1 = REAL( l*(l+1), idp )

                    DO dp = 0,DEGREE

                        Reusable_Values(dp) = SUM ( ( -R_SQUARE(:) * LPT_LPT(:,d,dp,1,1)*TWOOVER_DELTAR*TWOOVER_DELTAR      &
                                                    + L_Lp1 * LPT_LPT(:,d,dp,0,0)        )    &
                                                    * INT_R_WEIGHTS(:)/TWOOVER_DELTAR     )
!                        Reusable_Values(dp) = SUM( -DRDR_Factor(:,d,dp) + L_Lp1*RR_Factor(:,d,dp) )
!                        PRINT*,Reusable_Values(dp)-SUM( -DRDR_Factor(:,d,dp) + L_Lp1*RR_Factor(:,d,dp) )
                    END DO




                    DO m = -M_VALUES(l),M_VALUES(l)

                        ! Fix This !
                       jloc = d*ULM_LENGTH + (u-1)*LM_LENGTH + (l*(l+1) + m)


                        DO dp = 0,DEGREE

                            ! Fix This !
                            iloc = dp*ULM_LENGTH + (u-1)*LM_LENGTH + (l*(l+1) + m)

                            MATVEC_LOC = jloc*ELEM_PROB_DIM + iloc

                            BLOCK_ELEM_STF_MATVEC(MATVEC_LOC,Block_re)              &
                                    = BLOCK_ELEM_STF_MATVEC(MATVEC_LOC,Block_re)    &
                                    + Reusable_Values(dp)




                        END DO ! dp Loop

                    END DO ! m Loop

                END DO ! l Loop

            END DO ! u Loop

        END DO ! d Loop

    END DO ! re Loop


    !$OMP END DO

    !$OMP END PARALLEL

END IF






END SUBROUTINE FINISH_3D_JACOBIAN_MATRIX














!+501+##############################################################################!
!                                                                                   !
!                  FINISH_3D_RHS_VECTOR                             !
!                                                                                   !
!###################################################################################!
SUBROUTINE FINISH_3D_RHS_VECTOR()


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
    !$OMP PRIVATE( ui, l, m, re, d, dp, rd, ui,                             &
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



                        Current_i_Location = Matrix_Location(ui, l, m, re, d)


                        DO dp = 0,DEGREE

                            Current_j_Location = NUM_OFF_DIAGONALS + (dp - d) * ULM_LENGTH


                            Block_STF_MAT(Current_j_Location, Current_i_Location )          &
                                    = Block_STF_MAT(Current_j_Location, Current_i_Location )      &
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

        CALL OUTPUT_RHS_VECTOR_Parts(TMP_VECTOR, BLOCK_RHS_VECTOR)

    END IF


    Block_RHS_Vector(Start_Here:End_Here) = Block_RHS_Vector(Start_Here:End_Here)   &
                                          + TMP_VECTOR(0:SUBSHELL_PROB_DIM-1)

END IF

END SUBROUTINE FINISH_3D_RHS_VECTOR









!+501+##############################################################################!
!                                                                                   !
!                  ANALYZE_3D_RHS_VECTOR                             !
!                                                                                   !
!###################################################################################!
SUBROUTINE ANALYZE_3D_RHS_VECTOR( Laplace, Source)


COMPLEX(KIND = idp), DIMENSION(0:SUBSHELL_PROB_DIM-1), INTENT(IN)   ::  Laplace
COMPLEX(KIND = idp), DIMENSION(0:SUBSHELL_PROB_DIM-1), INTENT(IN)   ::  Source

COMPLEX(KIND = idp), DIMENSION(1:NUM_CFA_VARS, 0:VAR_DIM-1)         ::  Laplace_Holder
COMPLEX(KIND = idp), DIMENSION(1:NUM_CFA_VARS, 0:VAR_DIM-1)         ::  Source_Holder

INTEGER                                                             ::  re, d, ui, lm, here


DO re = 0,NUM_R_ELEMENTS - 1
    DO D = 0,DEGREE
        DO ui = 1,3
            DO lm = 0,LM_LENGTH - 1

                Here = (re*DEGREE + D)*ULM_LENGTH       &
                     + (ui -1)*LM_LENGTH                &
                     + lm
        
                !Laplace_Holder(ui, )
            

            END DO ! lm Loop
        END DO ! ui Loop
    END DO ! D Loop
END DO ! re Loop



END SUBROUTINE ANALYZE_3D_RHS_VECTOR









!+???+###########################################################################!
!                                                                                !
!                  Jacobian_Term          !
!                                                                                !
!################################################################################!
FUNCTION Calc_Jacobian_Term( F, u, l, lp, d, dp,            &
                            SubJacobian_EQ1_Term,           &
                            SubJacobian_EQ2_Term,           &
                            SubJacobian_EQ3_Term,           &
                            SubJacobian_EQ4_Term,           &
                            SubJacobian_EQ5_Term            )


INTEGER, INTENT(IN)                         :: F, u, l, lp, d, dp


REAL(KIND = idp), INTENT(IN), DIMENSION( 1:NUM_TP_QUAD_POINTS,   &
                                         1:NUM_R_QUAD_POINTS,    &
                                         1:14                    )           ::  SubJacobian_EQ1_Term

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:NUM_TP_QUAD_POINTS,   &
                                         1:NUM_R_QUAD_POINTS,    &
                                         1:14                    )           ::  SubJacobian_EQ2_Term

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:NUM_TP_QUAD_POINTS,   &
                                         1:NUM_R_QUAD_POINTS,    &
                                         1:20                    )           ::  SubJacobian_EQ3_Term

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:NUM_TP_QUAD_POINTS,   &
                                         1:NUM_R_QUAD_POINTS,    &
                                         1:20                    )           ::  SubJacobian_EQ4_Term

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:NUM_TP_QUAD_POINTS,   &
                                         1:NUM_R_QUAD_POINTS,    &
                                         1:20                    )           ::  SubJacobian_EQ5_Term


COMPLEX                                     :: Calc_Jacobian_Term, Jacobian_Term

INTEGER                                     :: Case_Number
INTEGER                                     :: rd, td, pd, tpd

Case_Number = (F-1) * 5 + u-1
Jacobian_Term = 0.0_idp
SELECT CASE ( Case_Number )
    CASE( 0 )
        ! Jacobian Term corresponding to d Eq. 1/ d u_1
        DO rd = 1,NUM_R_QUAD_POINTS
            Jacobian_Term = Jacobian_Term                                                       &
                          + SUM( SubJacobian_EQ1_Term( :, rd, 1 ) * TP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp )
        END DO
    CASE( 1 )
        ! Jacobian Term corresponding to d Eq. 1/ d u_2
        DO rd = 1,NUM_R_QUAD_POINTS
            Jacobian_Term = Jacobian_Term                                                       &
                          + SUM( SubJacobian_EQ1_Term( :, rd, 2 ) * TP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp )
        END DO
    CASE( 2 )
        ! Jacobian Term corresponding to d Eq. 1/ d u_3
        DO rd = 1,NUM_R_QUAD_POINTS
            Jacobian_Term = Jacobian_Term                                                       &
                          + SUM( SubJacobian_EQ1_Term(:, rd, 3) * TP_TP_Factor( :, l, lp ) )    &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ1_Term(:, rd, 4) * TP_TP_Factor( :, l, lp ) )    &
                                * DRR_Factor( rd, d, dp)                                        &
                          + SUM( SubJacobian_EQ1_Term(:, rd, 5) * dTP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ1_Term(:, rd, 6) * TdP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)
        END DO
    CASE( 3 )
        ! Jacobian Term corresponding to d Eq. 1/ d u_4
        DO rd = 1,NUM_R_QUAD_POINTS
            Jacobian_Term = Jacobian_Term                                                       &
                          + SUM( SubJacobian_EQ1_Term(:, rd, 7) * TP_TP_Factor( :, l, lp ) )    &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ1_Term(:, rd, 8) * TP_TP_Factor( :, l, lp ) )    &
                                * DRR_Factor( rd, d, dp)                                        &
                          + SUM( SubJacobian_EQ1_Term(:, rd, 9) * dTP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ1_Term(:, rd, 10) * TdP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp)
        END DO
    CASE( 4 )
        ! Jacobian Term corresponding to d Eq. 1/ d u_5
        DO rd = 1,NUM_R_QUAD_POINTS
            Jacobian_Term = Jacobian_Term                                                       &
                          + SUM( SubJacobian_EQ1_Term(:, rd, 11) * TP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ1_Term(:, rd, 12) * TP_TP_Factor( :, l, lp ) )   &
                                * DRR_Factor( rd, d, dp)                                        &
                          + SUM( SubJacobian_EQ1_Term(:, rd, 13) * dTP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ1_Term(:, rd, 14) * TdP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp)
        END DO
    CASE( 5 )
        ! Jacobian Term corresponding to d Eq. 2/ d u_1
        DO rd = 1,NUM_R_QUAD_POINTS
            Jacobian_Term = Jacobian_Term                                                       &
                          + SUM( SubJacobian_EQ2_Term( :, rd, 1 ) * TP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp )
        END DO
    CASE( 6 )
        ! Jacobian Term corresponding to d Eq. 2/ d u_2
        DO rd = 1,NUM_R_QUAD_POINTS
            Jacobian_Term = Jacobian_Term                                                       &
                          + SUM( SubJacobian_EQ2_Term( :, rd, 2 ) * TP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp )
        END DO
    CASE( 7 )
        ! Jacobian Term corresponding to d Eq. 2/ d u_3
        DO rd = 1,NUM_R_QUAD_POINTS
            Jacobian_Term = Jacobian_Term                                                       &
                          + SUM( SubJacobian_EQ2_Term(:, rd, 3) * TP_TP_Factor( :, l, lp ) )    &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ2_Term(:, rd, 4) * TP_TP_Factor( :, l, lp ) )    &
                                * DRR_Factor( rd, d, dp)                                        &
                          + SUM( SubJacobian_EQ2_Term(:, rd, 5) * dTP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ2_Term(:, rd, 6) * TdP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)
        END DO
    CASE( 8 )
        ! Jacobian Term corresponding to d Eq. 2/ d u_4
        DO rd = 1,NUM_R_QUAD_POINTS
            Jacobian_Term = Jacobian_Term                                                       &
                          + SUM( SubJacobian_EQ2_Term(:, rd, 7) * TP_TP_Factor( :, l, lp ) )    &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ2_Term(:, rd, 8) * TP_TP_Factor( :, l, lp ) )    &
                                * DRR_Factor( rd, d, dp)                                        &
                          + SUM( SubJacobian_EQ2_Term(:, rd, 9) * dTP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ2_Term(:, rd, 10) * TdP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp)
        END DO
    CASE( 9 )
        ! Jacobian Term corresponding to d Eq. 2/ d u_5
        DO rd = 1,NUM_R_QUAD_POINTS
            Jacobian_Term = Jacobian_Term                                                       &
                          + SUM( SubJacobian_EQ2_Term(:, rd, 11) * TP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ2_Term(:, rd, 12) * TP_TP_Factor( :, l, lp ) )   &
                                * DRR_Factor( rd, d, dp)                                        &
                          + SUM( SubJacobian_EQ2_Term(:, rd, 13) * dTP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ2_Term(:, rd, 14) * TdP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp)
        END DO
    CASE( 10 )
        ! Jacobian Term corresponding to d Eq. 3/ d u_1
        DO rd = 1,NUM_R_QUAD_POINTS
            Jacobian_Term = Jacobian_Term                                                       &
                          + SUM( SubJacobian_EQ3_Term(:, rd, 1) * TP_TP_Factor( :, l, lp ) )    &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ3_Term(:, rd, 2) * TP_TP_Factor( :, l, lp ) )    &
                                * DRR_Factor( rd, d, dp)                                        &
                          + SUM( SubJacobian_EQ3_Term(:, rd, 3) * dTP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ3_Term(:, rd, 4) * TdP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)
        END DO
    CASE( 11 )
        ! Jacobian Term corresponding to d Eq. 3/ d u_2
       DO rd = 1,NUM_R_QUAD_POINTS
            Jacobian_Term = Jacobian_Term                                                       &
                          + SUM( SubJacobian_EQ3_Term(:, rd, 5) * TP_TP_Factor( :, l, lp ) )    &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ3_Term(:, rd, 6) * TP_TP_Factor( :, l, lp ) )    &
                                * DRR_Factor( rd, d, dp)                                        &
                          + SUM( SubJacobian_EQ3_Term(:, rd, 7) * dTP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ3_Term(:, rd, 8) * TdP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)
        END DO
    CASE( 12 )
        ! Jacobian Term corresponding to d Eq. 3/ d u_3
       DO rd = 1,NUM_R_QUAD_POINTS
            Jacobian_Term = Jacobian_Term                                                       &
                          + SUM( SubJacobian_EQ3_Term(:, rd, 9) * TP_TP_Factor( :, l, lp ) )    &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ3_Term(:, rd, 10) * TP_TP_Factor( :, l, lp ) )   &
                                * DRR_Factor( rd, d, dp)                                        &
                          + SUM( SubJacobian_EQ3_Term(:, rd, 11) * dTP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ3_Term(:, rd, 12) * TdP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp)                                         &
                          + OneThird * SUM( TP_TP_Factor(:,l,lp) * DRDR_Factor( rd, d, dp) )
        END DO
    CASE( 13 )
        ! Jacobian Term corresponding to d Eq. 3/ d u_4
       DO rd = 1,NUM_R_QUAD_POINTS
            Jacobian_Term = Jacobian_Term                                                       &
                          + SUM( SubJacobian_EQ3_Term(:, rd, 13) * TP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ3_Term(:, rd, 14) * TP_TP_Factor( :, l, lp ) )   &
                                * DRR_Factor( rd, d, dp)                                        &
                          + SUM( SubJacobian_EQ3_Term(:, rd, 15) * dTP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ3_Term(:, rd, 16) * TdP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp)                                         &
                          + OneThird * SUM( dTP_TP_Factor(:,l,lp) * RDR_Factor( rd, d, dp) )
        END DO
    CASE( 14 )
        ! Jacobian Term corresponding to d Eq. 3/ d u_5
       DO rd = 1,NUM_R_QUAD_POINTS
            Jacobian_Term = Jacobian_Term                                                       &
                          + SUM( SubJacobian_EQ3_Term(:, rd, 17) * TP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ3_Term(:, rd, 18) * TP_TP_Factor( :, l, lp ) )   &
                                * DRR_Factor( rd, d, dp)                                        &
                          + SUM( SubJacobian_EQ3_Term(:, rd, 19) * dTP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ3_Term(:, rd, 20) * TdP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp)                                         &
                          + OneThird * SUM( TdP_TP_Factor(:,l,lp) * RDR_Factor( rd, d, dp) )
        END DO
    CASE( 15 )
        ! Jacobian Term corresponding to d Eq. 4/ d u_1
        DO rd = 1,NUM_R_QUAD_POINTS
            Jacobian_Term = Jacobian_Term                                                       &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 1) * TP_TP_Factor( :, l, lp ) )    &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 2) * TP_TP_Factor( :, l, lp ) )    &
                                * DRR_Factor( rd, d, dp)                                        &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 3) * dTP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 4) * TdP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)
        END DO
    CASE( 16 )
        ! Jacobian Term corresponding to d Eq. 4/ d u_2
        DO rd = 1,NUM_R_QUAD_POINTS
            Jacobian_Term = Jacobian_Term                                                       &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 5) * TP_TP_Factor( :, l, lp ) )    &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 6) * TP_TP_Factor( :, l, lp ) )    &
                                * DRR_Factor( rd, d, dp)                                        &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 7) * dTP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 8) * TdP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)
        END DO
    CASE( 17 )
        ! Jacobian Term corresponding to d Eq. 4/ d u_3
        DO rd = 1,NUM_R_QUAD_POINTS
            Jacobian_Term = Jacobian_Term                                                       &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 9) * TP_TP_Factor( :, l, lp ) )    &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 10) * TP_TP_Factor( :, l, lp ) )   &
                                * DRR_Factor( rd, d, dp)                                        &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 11) * dTP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 12) * TdP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp)                                         &
                          + OneThird * SUM( TP_dTP_Factor(:,l,lp) )                             &
                                * DRR_Factor( rd, d, dp)/R_SQUARE(rd)
        END DO
    CASE( 18 )
        ! Jacobian Term corresponding to d Eq. 4/ d u_4
        DO rd = 1,NUM_R_QUAD_POINTS
            Jacobian_Term = Jacobian_Term                                                       &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 13) * TP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 14) * TP_TP_Factor( :, l, lp ) )   &
                                * DRR_Factor( rd, d, dp)                                        &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 15) * dTP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 16) * TdP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp)                                         &
                          + OneThird * SUM( dTP_dTP_Factor(:,l,lp) )                            &
                                * RR_Factor( rd, d, dp)/R_SQUARE(rd)
        END DO
    CASE( 19 )
        ! Jacobian Term corresponding to d Eq. 4/ d u_5
        DO rd = 1,NUM_R_QUAD_POINTS
            Jacobian_Term = Jacobian_Term                                                       &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 17) * TP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 18) * TP_TP_Factor( :, l, lp ) )   &
                                * DRR_Factor( rd, d, dp)                                        &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 19) * dTP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 20) * TdP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp)                                         &
                          + OneThird * SUM( TdP_dTP_Factor(:,l,lp) )                            &
                                * RR_Factor( rd, d, dp)/R_SQUARE(rd)
        END DO
    CASE( 20 )
        ! Jacobian Term corresponding to d Eq. 5/ d u_1
        DO rd = 1,NUM_R_QUAD_POINTS
            Jacobian_Term = Jacobian_Term                                                       &
                          + SUM( SubJacobian_EQ5_Term(:, rd, 1) * TP_TP_Factor( :, l, lp ) )    &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ5_Term(:, rd, 2) * TP_TP_Factor( :, l, lp ) )    &
                                * DRR_Factor( rd, d, dp)                                        &
                          + SUM( SubJacobian_EQ5_Term(:, rd, 3) * dTP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ5_Term(:, rd, 4) * TdP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)
        END DO
    CASE( 21 )
        ! Jacobian Term corresponding to d Eq. 5/ d u_2
        DO rd = 1,NUM_R_QUAD_POINTS
            Jacobian_Term = Jacobian_Term                                                       &
                          + SUM( SubJacobian_EQ5_Term(:, rd, 5) * TP_TP_Factor( :, l, lp ) )    &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ5_Term(:, rd, 6) * TP_TP_Factor( :, l, lp ) )    &
                                * DRR_Factor( rd, d, dp)                                        &
                          + SUM( SubJacobian_EQ5_Term(:, rd, 7) * dTP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ5_Term(:, rd, 8) * TdP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)
        END DO
    CASE( 22 )
        ! Jacobian Term corresponding to d Eq. 5/ d u_3
        DO rd = 1,NUM_R_QUAD_POINTS
            Jacobian_Term = Jacobian_Term                                                       &
                          + SUM( SubJacobian_EQ5_Term(:, rd, 9) * TP_TP_Factor( :, l, lp ) )    &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ5_Term(:, rd, 10) * TP_TP_Factor( :, l, lp ) )   &
                                * DRR_Factor( rd, d, dp)                                        &
                          + SUM( SubJacobian_EQ5_Term(:, rd, 11) * dTP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ5_Term(:, rd, 12) * TdP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp)                                         &
                          + OneThird * SUM( TP_TdP_Factor(:,l,lp)/(R_SQUARE(rd)*TP_SIN_SQUARE(:) ) )  &
                                * DRR_Factor( rd, d, dp)
        END DO
    CASE( 23 )
        ! Jacobian Term corresponding to d Eq. 5/ d u_4
        DO rd = 1,NUM_R_QUAD_POINTS
            Jacobian_Term = Jacobian_Term                                                       &
                          + SUM( SubJacobian_EQ5_Term(:, rd, 13) * TP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ5_Term(:, rd, 14) * TP_TP_Factor( :, l, lp ) )   &
                                * DRR_Factor( rd, d, dp)                                        &
                          + SUM( SubJacobian_EQ5_Term(:, rd, 15) * dTP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ5_Term(:, rd, 16) * TdP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp)                                         &
                          + OneThird * SUM( dTP_TdP_Factor(:,l,lp)/(R_SQUARE(rd)*TP_SIN_SQUARE(:) ))           &
                                * RR_Factor( rd, d, dp)
        END DO
    CASE( 24 )
        ! Jacobian Term corresponding to d Eq. 5/ d u_5
        DO rd = 1,NUM_R_QUAD_POINTS
            Jacobian_Term = Jacobian_Term                                                       &
                          + SUM( SubJacobian_EQ5_Term(:, rd, 17) * TP_TP_Factor( :, l, lp ) )   &
                                    * RR_Factor( rd, d, dp)                                     &
                          + SUM( SubJacobian_EQ5_Term(:, rd, 18) * TP_TP_Factor( :, l, lp ) )   &
                                * DRR_Factor( rd, d, dp)                                        &
                          + SUM( SubJacobian_EQ5_Term(:, rd, 19) * dTP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ5_Term(:, rd, 20) * TdP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp)                                         &
                          + OneThird * SUM( TdP_TdP_Factor(:,l,lp)/(R_SQUARE(rd)*TP_SIN_SQUARE(:) ) )          &
                                * RR_Factor( rd, d, dp)
        END DO
END SELECT


Calc_Jacobian_Term = Jacobian_Term

END FUNCTION Calc_Jacobian_Term




!+101+###########################################################################!
!                                                                                !
!                  CALC_JACOBIAN_A                                               !
!                                                                                !
!################################################################################!
SUBROUTINE CALC_JACOBIAN_A( Jacobian_Buffer, dp, lpmp_loc, F )

COMPLEX(KIND = idp), DIMENSION(0:ELEM_PROB_DIM-1), INTENT(OUT)      :: Jacobian_Buffer
INTEGER, INTENT(IN)                                                 :: dp, lpmp_loc, F

INTEGER                                                             :: j_loc, d, u, lm_loc

! d, u, and lm_loc pick the column
DO d = 0,DEGREE
    DO u = 1,NUM_CFA_VARS
        DO lm_loc = 0,LM_LENGTH-1

            

            ! column location
            j_loc = d*ULM_LENGTH + (u-1)*LM_LENGTH + lm_loc



            Jacobian_Buffer(j_loc) = Calc_Jacobian_Term( F, u, lm_loc, lpmp_loc, d, dp, &
                                                        SubJacobian_EQ1_Term,           &
                                                        SubJacobian_EQ2_Term,           &
                                                        SubJacobian_EQ3_Term,           &
                                                        SubJacobian_EQ4_Term,           &
                                                        SubJacobian_EQ5_Term            )

        END DO ! lm_loc Loop
    END DO ! u Loop
END DO ! d Loop

END SUBROUTINE CALC_JACOBIAN_A




!+101+###########################################################################!
!                                                                                !
!                  CALC_JACOBIAN_B                                               !
!                                                                                !
!################################################################################!
SUBROUTINE CALC_JACOBIAN_B( Jacobian_Buffer, dp, lp, F )

COMPLEX(KIND = idp), DIMENSION(0:ELEM_PROB_DIM-1), INTENT(OUT)    :: Jacobian_Buffer
INTEGER, INTENT(IN)                                                 :: dp, lp, F

INTEGER                                                             :: j_loc, d, u, l, rd

Jacobian_Buffer = 0.0_idp
! d, u, and l pick the column
IF ( F == 1 ) THEN
    DO d = 0,DEGREE

        ! Jacobian Term corresponding to d Eq. 1/ d u_1
        u = 1
        DO l = 0,LM_LENGTH-1
            j_loc = d*ULM_LENGTH + (u-1)*LM_LENGTH + l
            DO rd = 1,NUM_R_QUAD_POINTS
                Jacobian_Buffer(j_loc) = Jacobian_Buffer(j_loc)                                     &
                               + SUM( SubJacobian_EQ1_Term( :, rd, 1 ) * TP_TP_Factor( :, l, lp ) )  &
                                     * RR_Factor( rd, d, dp )
            END DO
        END DO

        ! Jacobian Term corresponding to d Eq. 1/ d u_2
        u = 2
        DO l = 0,LM_LENGTH-1
            j_loc = d*ULM_LENGTH + (u-1)*LM_LENGTH + l
            DO rd = 1,NUM_R_QUAD_POINTS
                 Jacobian_Buffer(j_loc) = Jacobian_Buffer(j_loc)                                        &
                                + SUM( SubJacobian_EQ1_Term( :, rd, 2 ) * TP_TP_Factor( :, l, lp ) )  &
                                    * RR_Factor( rd, d, dp )
            END DO
        END DO

        ! Jacobian Term corresponding to d Eq. 1/ d u_3
        u = 3
        DO l = 0,LM_LENGTH-1
            j_loc = d*ULM_LENGTH + (u-1)*LM_LENGTH + l
            DO rd = 1,NUM_R_QUAD_POINTS
                Jacobian_Buffer(j_loc) = Jacobian_Buffer(j_loc)                                     &
                               + SUM( SubJacobian_EQ1_Term(:, rd, 3) * TP_TP_Factor( :, l, lp ) )    &
                                     * RR_Factor( rd, d, dp)                                         &
                               + SUM( SubJacobian_EQ1_Term(:, rd, 4) * TP_TP_Factor( :, l, lp ) )    &
                                     * DRR_Factor( rd, d, dp)                                        &
                               + SUM( SubJacobian_EQ1_Term(:, rd, 5) * dTP_TP_Factor( :, l, lp ) )   &
                                     * RR_Factor( rd, d, dp)                                         &
                               + SUM( SubJacobian_EQ1_Term(:, rd, 6) * TdP_TP_Factor( :, l, lp ) )   &
                                     * RR_Factor( rd, d, dp)
            END DO
        END DO



        ! Jacobian Term corresponding to d Eq. 1/ d u_4
        u = 4
        DO l = 0,LM_LENGTH-1
            j_loc = d*ULM_LENGTH + (u-1)*LM_LENGTH + l
            DO rd = 1,NUM_R_QUAD_POINTS
                Jacobian_Buffer(j_loc) = Jacobian_Buffer(j_loc)                                      &
                               + SUM( SubJacobian_EQ1_Term(:, rd, 7) * TP_TP_Factor( :, l, lp ) )    &
                                     * RR_Factor( rd, d, dp)                                         &
                               + SUM( SubJacobian_EQ1_Term(:, rd, 8) * TP_TP_Factor( :, l, lp ) )    &
                                     * DRR_Factor( rd, d, dp)                                        &
                               + SUM( SubJacobian_EQ1_Term(:, rd, 9) * dTP_TP_Factor( :, l, lp ) )   &
                                     * RR_Factor( rd, d, dp)                                         &
                               + SUM( SubJacobian_EQ1_Term(:, rd, 10) * TdP_TP_Factor( :, l, lp ) )  &
                                     * RR_Factor( rd, d, dp)
            END DO
        END DO



        ! Jacobian Term corresponding to d Eq. 1/ d u_5
        u = 5
        DO l = 0,LM_LENGTH-1
            j_loc = d*ULM_LENGTH + (u-1)*LM_LENGTH + l
            DO rd = 1,NUM_R_QUAD_POINTS
                Jacobian_Buffer(j_loc) = Jacobian_Buffer(j_loc)                                      &
                               + SUM( SubJacobian_EQ1_Term(:, rd, 11) * TP_TP_Factor( :, l, lp ) )   &
                                     * RR_Factor( rd, d, dp)                                         &
                               + SUM( SubJacobian_EQ1_Term(:, rd, 12) * TP_TP_Factor( :, l, lp ) )   &
                                     * DRR_Factor( rd, d, dp)                                        &
                               + SUM( SubJacobian_EQ1_Term(:, rd, 13) * dTP_TP_Factor( :, l, lp ) )  &
                                     * RR_Factor( rd, d, dp)                                         &
                               + SUM( SubJacobian_EQ1_Term(:, rd, 14) * TdP_TP_Factor( :, l, lp ) )  &
                                     * RR_Factor( rd, d, dp)
            END DO
        END DO

    END DO ! d Loop
ELSE IF ( F == 2 ) THEN
    DO d = 0,DEGREE

        ! Jacobian Term corresponding to d Eq. 2/ d u_1
        u = 1
        DO l = 0,LM_LENGTH-1
            j_loc = d*ULM_LENGTH + (u-1)*LM_LENGTH + l
            DO rd = 1,NUM_R_QUAD_POINTS
                Jacobian_Buffer(j_loc) = Jacobian_Buffer(j_loc)                                     &
                              + SUM( SubJacobian_EQ2_Term( :, rd, 1 ) * TP_TP_Factor( :, l, lp ) )  &
                                    * RR_Factor( rd, d, dp )
            END DO  ! rd Loop
        END DO ! l Loop


    ! Jacobian Term corresponding to d Eq. 2/ d u_2
    u = 2
    DO l = 0,LM_LENGTH-1
        j_loc = d*ULM_LENGTH + (u-1)*LM_LENGTH + l
        DO rd = 1,NUM_R_QUAD_POINTS
            Jacobian_Buffer(j_loc) = Jacobian_Buffer(j_loc)                                     &
                          + SUM( SubJacobian_EQ2_Term( :, rd, 2 ) * TP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp )
        END DO ! rd Loop
    END DO ! l Loop

    ! Jacobian Term corresponding to d Eq. 2/ d u_3
    u = 3
    DO l = 0,LM_LENGTH-1
        j_loc = d*ULM_LENGTH + (u-1)*LM_LENGTH + l
        DO rd = 1,NUM_R_QUAD_POINTS
            Jacobian_Buffer(j_loc) = Jacobian_Buffer(j_loc)                                     &
                          + SUM( SubJacobian_EQ2_Term(:, rd, 3) * TP_TP_Factor( :, l, lp ) )    &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ2_Term(:, rd, 4) * TP_TP_Factor( :, l, lp ) )    &
                                * DRR_Factor( rd, d, dp)                                        &
                          + SUM( SubJacobian_EQ2_Term(:, rd, 5) * dTP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ2_Term(:, rd, 6) * TdP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)
        END DO ! rd Loop
    END DO ! l
    
    ! Jacobian Term corresponding to d Eq. 2/ d u_4
    u = 4
    DO l = 0,LM_LENGTH-1
        j_loc = d*ULM_LENGTH + (u-1)*LM_LENGTH + l
        DO rd = 1,NUM_R_QUAD_POINTS
            Jacobian_Buffer(j_loc) = Jacobian_Buffer(j_loc)                                     &
                          + SUM( SubJacobian_EQ2_Term(:, rd, 7) * TP_TP_Factor( :, l, lp ) )    &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ2_Term(:, rd, 8) * TP_TP_Factor( :, l, lp ) )    &
                                * DRR_Factor( rd, d, dp)                                        &
                          + SUM( SubJacobian_EQ2_Term(:, rd, 9) * dTP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ2_Term(:, rd, 10) * TdP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp)
        END DO ! rd Loop
    END DO ! l Loop


    ! Jacobian Term corresponding to d Eq. 2/ d u_5
    u = 5
    DO l = 0,LM_LENGTH-1
        j_loc = d*ULM_LENGTH + (u-1)*LM_LENGTH + l
        DO rd = 1,NUM_R_QUAD_POINTS
            Jacobian_Buffer(j_loc) = Jacobian_Buffer(j_loc)                                     &
                          + SUM( SubJacobian_EQ2_Term(:, rd, 11) * TP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ2_Term(:, rd, 12) * TP_TP_Factor( :, l, lp ) )   &
                                * DRR_Factor( rd, d, dp)                                        &
                          + SUM( SubJacobian_EQ2_Term(:, rd, 13) * dTP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ2_Term(:, rd, 14) * TdP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp)
        END DO ! rd Loop
    END DO ! l Loop

    END DO ! d Loop


ELSE IF ( F == 3 ) THEN
    
    DO d = 0,DEGREE

        ! Jacobian Term corresponding to d Eq. 3/ d u_1
        u = 1
        DO l = 0,LM_LENGTH-1
            j_loc = d*ULM_LENGTH + (u-1)*LM_LENGTH + l
            DO rd = 1,NUM_R_QUAD_POINTS
                Jacobian_Buffer(j_loc) = Jacobian_Buffer(j_loc)                                                       &
                              + SUM( SubJacobian_EQ3_Term(:, rd, 1) * TP_TP_Factor( :, l, lp ) )    &
                                    * RR_Factor( rd, d, dp)                                         &
                              + SUM( SubJacobian_EQ3_Term(:, rd, 2) * TP_TP_Factor( :, l, lp ) )    &
                                    * DRR_Factor( rd, d, dp)                                        &
                              + SUM( SubJacobian_EQ3_Term(:, rd, 3) * dTP_TP_Factor( :, l, lp ) )   &
                                    * RR_Factor( rd, d, dp)                                         &
                              + SUM( SubJacobian_EQ3_Term(:, rd, 4) * TdP_TP_Factor( :, l, lp ) )   &
                                    * RR_Factor( rd, d, dp)
            END DO ! rd Loop
        END DO ! l Loop



        ! Jacobian Term corresponding to d Eq. 3/ d u_2
        u = 2
        DO l = 0,LM_LENGTH-1
            j_loc = d*ULM_LENGTH + (u-1)*LM_LENGTH + l
            DO rd = 1,NUM_R_QUAD_POINTS
                Jacobian_Buffer(j_loc) = Jacobian_Buffer(j_loc)                                                       &
                              + SUM( SubJacobian_EQ3_Term(:, rd, 5) * TP_TP_Factor( :, l, lp ) )    &
                                    * RR_Factor( rd, d, dp)                                         &
                              + SUM( SubJacobian_EQ3_Term(:, rd, 6) * TP_TP_Factor( :, l, lp ) )    &
                                    * DRR_Factor( rd, d, dp)                                        &
                              + SUM( SubJacobian_EQ3_Term(:, rd, 7) * dTP_TP_Factor( :, l, lp ) )   &
                                    * RR_Factor( rd, d, dp)                                         &
                              + SUM( SubJacobian_EQ3_Term(:, rd, 8) * TdP_TP_Factor( :, l, lp ) )   &
                                    * RR_Factor( rd, d, dp)
            END DO ! rd Loop
        END DO ! l Loop

        ! Jacobian Term corresponding to d Eq. 3/ d u_3
        u = 3
        DO l = 0,LM_LENGTH-1
            j_loc = d*ULM_LENGTH + (u-1)*LM_LENGTH + l
            DO rd = 1,NUM_R_QUAD_POINTS
                Jacobian_Buffer(j_loc) = Jacobian_Buffer(j_loc)                                                       &
                              + SUM( SubJacobian_EQ3_Term(:, rd, 9) * TP_TP_Factor( :, l, lp ) )    &
                                    * RR_Factor( rd, d, dp)                                         &
                              + SUM( SubJacobian_EQ3_Term(:, rd, 10) * TP_TP_Factor( :, l, lp ) )   &
                                    * DRR_Factor( rd, d, dp)                                        &
                              + SUM( SubJacobian_EQ3_Term(:, rd, 11) * dTP_TP_Factor( :, l, lp ) )  &
                                    * RR_Factor( rd, d, dp)                                         &
                              + SUM( SubJacobian_EQ3_Term(:, rd, 12) * TdP_TP_Factor( :, l, lp ) )  &
                                    * RR_Factor( rd, d, dp)                                         &
                              + OneThird * SUM( TP_TP_Factor(:,l,lp) * DRDR_Factor( rd, d, dp) )
            END DO ! rd Loop
        END DO ! l Loop

        ! Jacobian Term corresponding to d Eq. 3/ d u_4
        u = 4
        DO l = 0,LM_LENGTH-1
            j_loc = d*ULM_LENGTH + (u-1)*LM_LENGTH + l
            DO rd = 1,NUM_R_QUAD_POINTS
                Jacobian_Buffer(j_loc) = Jacobian_Buffer(j_loc)                                                       &
                              + SUM( SubJacobian_EQ3_Term(:, rd, 13) * TP_TP_Factor( :, l, lp ) )   &
                                    * RR_Factor( rd, d, dp)                                         &
                              + SUM( SubJacobian_EQ3_Term(:, rd, 14) * TP_TP_Factor( :, l, lp ) )   &
                                    * DRR_Factor( rd, d, dp)                                        &
                              + SUM( SubJacobian_EQ3_Term(:, rd, 15) * dTP_TP_Factor( :, l, lp ) )  &
                                    * RR_Factor( rd, d, dp)                                         &
                              + SUM( SubJacobian_EQ3_Term(:, rd, 16) * TdP_TP_Factor( :, l, lp ) )  &
                                    * RR_Factor( rd, d, dp)                                         &
                              + OneThird * SUM( dTP_TP_Factor(:,l,lp) * RDR_Factor( rd, d, dp) )
            END DO ! rd Loop
        END DO ! l Loop

        ! Jacobian Term corresponding to d Eq. 3/ d u_5
        u = 5
        DO l = 0,LM_LENGTH-1
            j_loc = d*ULM_LENGTH + (u-1)*LM_LENGTH + l
            DO rd = 1,NUM_R_QUAD_POINTS
                Jacobian_Buffer(j_loc) = Jacobian_Buffer(j_loc)                                                       &
                              + SUM( SubJacobian_EQ3_Term(:, rd, 17) * TP_TP_Factor( :, l, lp ) )   &
                                    * RR_Factor( rd, d, dp)                                         &
                              + SUM( SubJacobian_EQ3_Term(:, rd, 18) * TP_TP_Factor( :, l, lp ) )   &
                                    * DRR_Factor( rd, d, dp)                                        &
                              + SUM( SubJacobian_EQ3_Term(:, rd, 19) * dTP_TP_Factor( :, l, lp ) )  &
                                    * RR_Factor( rd, d, dp)                                         &
                              + SUM( SubJacobian_EQ3_Term(:, rd, 20) * TdP_TP_Factor( :, l, lp ) )  &
                                    * RR_Factor( rd, d, dp)                                         &
                              + OneThird * SUM( TdP_TP_Factor(:,l,lp) * RDR_Factor( rd, d, dp) )
            END DO ! rd Loop
        END DO ! l Loop

    END DO ! d Loop

ELSE IF ( F == 4 ) THEN

    DO d = 0,DEGREE


    ! Jacobian Term corresponding to d Eq. 4/ d u_1
    u = 1
    DO l = 0,LM_LENGTH-1
        j_loc = d*ULM_LENGTH + (u-1)*LM_LENGTH + l
        DO rd = 1,NUM_R_QUAD_POINTS
            Jacobian_Buffer(j_loc) = Jacobian_Buffer(j_loc)                                                       &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 1) * TP_TP_Factor( :, l, lp ) )    &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 2) * TP_TP_Factor( :, l, lp ) )    &
                                * DRR_Factor( rd, d, dp)                                        &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 3) * dTP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 4) * TdP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)
        END DO ! rd Loop
    END DO ! l Loop

    ! Jacobian Term corresponding to d Eq. 4/ d u_2
    u = 2
    DO l = 0,LM_LENGTH-1
        j_loc = d*ULM_LENGTH + (u-1)*LM_LENGTH + l
        DO rd = 1,NUM_R_QUAD_POINTS
            Jacobian_Buffer(j_loc) = Jacobian_Buffer(j_loc)                                                       &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 5) * TP_TP_Factor( :, l, lp ) )    &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 6) * TP_TP_Factor( :, l, lp ) )    &
                                * DRR_Factor( rd, d, dp)                                        &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 7) * dTP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 8) * TdP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)
        END DO ! rd Loop
    END DO ! l Loop

    ! Jacobian Term corresponding to d Eq. 4/ d u_3
    u = 3
    DO l = 0,LM_LENGTH-1
        j_loc = d*ULM_LENGTH + (u-1)*LM_LENGTH + l
        DO rd = 1,NUM_R_QUAD_POINTS
            Jacobian_Buffer(j_loc) = Jacobian_Buffer(j_loc)                                                       &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 9) * TP_TP_Factor( :, l, lp ) )    &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 10) * TP_TP_Factor( :, l, lp ) )   &
                                * DRR_Factor( rd, d, dp)                                        &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 11) * dTP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 12) * TdP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp)                                         &
                          + OneThird * SUM( TP_dTP_Factor(:,l,lp) )                             &
                                * DRR_Factor( rd, d, dp)/R_SQUARE(rd)
        END DO ! rd Loop
    END DO ! l Loop
    
    ! Jacobian Term corresponding to d Eq. 4/ d u_4
    u = 4
    DO l = 0,LM_LENGTH-1
        j_loc = d*ULM_LENGTH + (u-1)*LM_LENGTH + l
        DO rd = 1,NUM_R_QUAD_POINTS
            Jacobian_Buffer(j_loc) = Jacobian_Buffer(j_loc)                                                       &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 13) * TP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 14) * TP_TP_Factor( :, l, lp ) )   &
                                * DRR_Factor( rd, d, dp)                                        &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 15) * dTP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 16) * TdP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp)                                         &
                          + OneThird * SUM( dTP_dTP_Factor(:,l,lp) )                            &
                                * RR_Factor( rd, d, dp)/R_SQUARE(rd)
        END DO ! rd Loop
    END DO ! l Loop

    ! Jacobian Term corresponding to d Eq. 4/ d u_5
    u = 5
    DO l = 0,LM_LENGTH-1
        j_loc = d*ULM_LENGTH + (u-1)*LM_LENGTH + l
        DO rd = 1,NUM_R_QUAD_POINTS
            Jacobian_Buffer(j_loc) = Jacobian_Buffer(j_loc)                                                       &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 17) * TP_TP_Factor( :, l, lp ) )   &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 18) * TP_TP_Factor( :, l, lp ) )   &
                                * DRR_Factor( rd, d, dp)                                        &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 19) * dTP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp)                                         &
                          + SUM( SubJacobian_EQ4_Term(:, rd, 20) * TdP_TP_Factor( :, l, lp ) )  &
                                * RR_Factor( rd, d, dp)                                         &
                          + OneThird * SUM( TdP_dTP_Factor(:,l,lp) )                            &
                                * RR_Factor( rd, d, dp)/R_SQUARE(rd)
        END DO ! rd Loop
    END DO ! l Loop

    END DO ! d Loop

ELSE IF ( F == 5 ) THEN
    DO d = 0,DEGREE

        ! Jacobian Term corresponding to d Eq. 5/ d u_1
        u = 1
        DO l = 0,LM_LENGTH-1
            j_loc = d*ULM_LENGTH + (u-1)*LM_LENGTH + l
            DO rd = 1,NUM_R_QUAD_POINTS
                Jacobian_Buffer(j_loc) = Jacobian_Buffer(j_loc)                                                       &
                              + SUM( SubJacobian_EQ5_Term(:, rd, 1) * TP_TP_Factor( :, l, lp ) )    &
                                    * RR_Factor( rd, d, dp)                                         &
                              + SUM( SubJacobian_EQ5_Term(:, rd, 2) * TP_TP_Factor( :, l, lp ) )    &
                                    * DRR_Factor( rd, d, dp)                                        &
                              + SUM( SubJacobian_EQ5_Term(:, rd, 3) * dTP_TP_Factor( :, l, lp ) )   &
                                    * RR_Factor( rd, d, dp)                                         &
                              + SUM( SubJacobian_EQ5_Term(:, rd, 4) * TdP_TP_Factor( :, l, lp ) )   &
                                    * RR_Factor( rd, d, dp)
            END DO ! rd Loop
        END DO ! l Loop


        ! Jacobian Term corresponding to d Eq. 5/ d u_2
        u = 2
        DO l = 0,LM_LENGTH-1
            j_loc = d*ULM_LENGTH + (u-1)*LM_LENGTH + l
            DO rd = 1,NUM_R_QUAD_POINTS
                Jacobian_Buffer(j_loc) = Jacobian_Buffer(j_loc)                                                       &
                              + SUM( SubJacobian_EQ5_Term(:, rd, 5) * TP_TP_Factor( :, l, lp ) )    &
                                    * RR_Factor( rd, d, dp)                                         &
                              + SUM( SubJacobian_EQ5_Term(:, rd, 6) * TP_TP_Factor( :, l, lp ) )    &
                                    * DRR_Factor( rd, d, dp)                                        &
                              + SUM( SubJacobian_EQ5_Term(:, rd, 7) * dTP_TP_Factor( :, l, lp ) )   &
                                    * RR_Factor( rd, d, dp)                                         &
                              + SUM( SubJacobian_EQ5_Term(:, rd, 8) * TdP_TP_Factor( :, l, lp ) )   &
                                    * RR_Factor( rd, d, dp)
            END DO ! rd Loop
        END DO ! l Loop


        ! Jacobian Term corresponding to d Eq. 5/ d u_3
        u = 3
        DO l = 0,LM_LENGTH-1
            j_loc = d*ULM_LENGTH + (u-1)*LM_LENGTH + l
            DO rd = 1,NUM_R_QUAD_POINTS
                Jacobian_Buffer(j_loc) = Jacobian_Buffer(j_loc)                                              &
                              + SUM( SubJacobian_EQ5_Term(:, rd, 9) * TP_TP_Factor( :, l, lp ) )    &
                                    * RR_Factor( rd, d, dp)                                         &
                              + SUM( SubJacobian_EQ5_Term(:, rd, 10) * TP_TP_Factor( :, l, lp ) )   &
                                    * DRR_Factor( rd, d, dp)                                        &
                              + SUM( SubJacobian_EQ5_Term(:, rd, 11) * dTP_TP_Factor( :, l, lp ) )  &
                                    * RR_Factor( rd, d, dp)                                         &
                              + SUM( SubJacobian_EQ5_Term(:, rd, 12) * TdP_TP_Factor( :, l, lp ) )  &
                                    * RR_Factor( rd, d, dp)                                         &
                              + OneThird * SUM( TP_TdP_Factor(:,l,lp)/(R_SQUARE(rd)*TP_SIN_SQUARE(:) ) )  &
                                    * DRR_Factor( rd, d, dp)
            END DO ! rd Loop
        END DO ! l Loop

        ! Jacobian Term corresponding to d Eq. 5/ d u_4
        u = 4
        DO l = 0,LM_LENGTH-1
            j_loc = d*ULM_LENGTH + (u-1)*LM_LENGTH + l
            DO rd = 1,NUM_R_QUAD_POINTS
                Jacobian_Buffer(j_loc) = Jacobian_Buffer(j_loc)                                           &
                              + SUM( SubJacobian_EQ5_Term(:, rd, 13) * TP_TP_Factor( :, l, lp ) )   &
                                    * RR_Factor( rd, d, dp)                                         &
                              + SUM( SubJacobian_EQ5_Term(:, rd, 14) * TP_TP_Factor( :, l, lp ) )   &
                                    * DRR_Factor( rd, d, dp)                                        &
                              + SUM( SubJacobian_EQ5_Term(:, rd, 15) * dTP_TP_Factor( :, l, lp ) )  &
                                    * RR_Factor( rd, d, dp)                                         &
                              + SUM( SubJacobian_EQ5_Term(:, rd, 16) * TdP_TP_Factor( :, l, lp ) )  &
                                    * RR_Factor( rd, d, dp)                                         &
                              + OneThird * SUM( dTP_TdP_Factor(:,l,lp)/(R_SQUARE(rd)*TP_SIN_SQUARE(:) ))           &
                                    * RR_Factor( rd, d, dp)
            END DO ! rd Loop
        END DO ! l Loop

        ! Jacobian Term corresponding to d Eq. 5/ d u_5
        u = 5
        DO l = 0,LM_LENGTH-1
            j_loc = d*ULM_LENGTH + (u-1)*LM_LENGTH + l
            DO rd = 1,NUM_R_QUAD_POINTS
                Jacobian_Buffer(j_loc) = Jacobian_Buffer(j_loc)                                        &
                              + SUM( SubJacobian_EQ5_Term(:, rd, 17) * TP_TP_Factor( :, l, lp ) )   &
                                        * RR_Factor( rd, d, dp)                                     &
                              + SUM( SubJacobian_EQ5_Term(:, rd, 18) * TP_TP_Factor( :, l, lp ) )   &
                                    * DRR_Factor( rd, d, dp)                                        &
                              + SUM( SubJacobian_EQ5_Term(:, rd, 19) * dTP_TP_Factor( :, l, lp ) )  &
                                    * RR_Factor( rd, d, dp)                                         &
                              + SUM( SubJacobian_EQ5_Term(:, rd, 20) * TdP_TP_Factor( :, l, lp ) )  &
                                    * RR_Factor( rd, d, dp)                                         &
                              + OneThird * SUM( TdP_TdP_Factor(:,l,lp)/(R_SQUARE(rd)*TP_SIN_SQUARE(:) ) )          &
                                    * RR_Factor( rd, d, dp)
            END DO ! rd Loop
        END DO ! l Loop

    END DO ! d Loop
END IF

END SUBROUTINE CALC_JACOBIAN_B








!+601+###########################################################################!
!                                                                                !
!           Allocate_Master_Build_Variables                                      !
!                                                                                !
!################################################################################!
SUBROUTINE Allocate_Master_Build_Variables()


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

ALLOCATE( SubJacobian_EQ1_Term( 1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:14) )
ALLOCATE( SubJacobian_EQ2_Term( 1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:14) )
ALLOCATE( SubJacobian_EQ3_Term( 1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:20) )
ALLOCATE( SubJacobian_EQ4_Term( 1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:20) )
ALLOCATE( SubJacobian_EQ5_Term( 1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:20) )


END SUBROUTINE Allocate_Master_Build_Variables



!+602+###########################################################################!
!                                                                                !
!           Deallocate_Master_Build_Variables                                    !
!                                                                                !
!################################################################################!
SUBROUTINE Deallocate_Master_Build_Variables()

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

DEALLOCATE( SubJacobian_EQ1_Term )
DEALLOCATE( SubJacobian_EQ2_Term )
DEALLOCATE( SubJacobian_EQ3_Term )
DEALLOCATE( SubJacobian_EQ4_Term )
DEALLOCATE( SubJacobian_EQ5_Term )


END SUBROUTINE Deallocate_Master_Build_Variables











END MODULE CFA_3D_Master_Build_Module
