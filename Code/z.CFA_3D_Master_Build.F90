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
!##!    +601+   CFA_3D_Apply_BCs_Part1                                              !##!
!##!    +602+   CFA_3D_Apply_BCs_Part2                                              !##!
!##!    +603+   CFA_3D_Dirichlet_BCs_Part1                                          !##!
!##!    +604+   CFA_3D_Dirichlet_BCs_Part2                                          !##!
!##!    +605+   CFA_3D_Neumann_BCs                                                  !##!
!##!                                                                                !##!
!##!    +701+   Calc_3D_Values_At_Location                                          !##!
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

USE Global_Variables_And_Parameters, &
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
                                ONLY :  JCBN_mu_FUNCTION_3D_ALL,    &
                                        JCBN_BIGK_SUBROUTINE,       &
                                        JCBN_BIGK_FUNCTION


USE IO_Functions_Module, &
                                ONLY :  Clock_In



IMPLICIT NONE





CONTAINS

!+101+###########################################################################!
!                                                                                !
!           CFA_3D_Master_Build                                                  !
!                                                                                !
!################################################################################!
SUBROUTINE CFA_3D_Master_Build()



REAL(KIND = idp)        :: timea, timeb, timec
INTEGER                 :: i,j

timea = 0.0_idp
timeb = 0.0_idp
timec = 0.0_idp

PRINT*,"Before CFA_3D_Apply_BCs_Part1"
PRINT*,Coefficient_Vector(0:4)
!*!
!*! Alter the Coefficient Vector to Reflect Boundary Conditions
!*!
timeb = MPI_Wtime()
CALL CFA_3D_Apply_BCs_Part1()
timec = MPI_Wtime()
CALL Clock_In(timec-timeb, 4)

PRINT*,"AFTER CFA_3D_Apply_BCs_Part1"
PRINT*,Coefficient_Vector(0:4)




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


PRINT*,"Before CFA_3D_Apply_BCs_Part2"
PRINT*,Coefficient_Vector(0:4)

!*!
!*! Alter the Jacobian Matrix to Reflect Boundary Conditions
!*!
timeb = MPI_Wtime()
CALL CFA_3D_Apply_BCs_Part2()
timec = MPI_Wtime()
CALL Clock_In(timec-timeb, 13)



PRINT*,"AFTER CFA_3D_Apply_BCs_Part2"
PRINT*,Coefficient_Vector(0:4)


END SUBROUTINE CFA_3D_Master_Build
















!+201+##########################################################################!
!                                                                               !
!                  CREATE_3D_NONLAPLACIAN_SOE                                   !
!                                                                               !
!###############################################################################!
SUBROUTINE CREATE_3D_NONLAPLACIAN_SOE()








INTEGER                                                         ::  Local_re,       &
                                                                    Local_te,       &
                                                                    Local_pe,       &
                                                                    Global_re,      &
                                                                    Global_te,      &
                                                                    Global_pe,      &
                                                                    d, l, m,        &
                                                                    dp, lp, mp,     &
                                                                    rd, td, pd




REAL(KIND = idp), DIMENSION( 1:NUM_R_QUAD_POINTS,   &
                             1:NUM_T_QUAD_POINTS,   &
                             1:NUM_P_QUAD_POINTS)               ::  CUR_VAL_PSI

REAL(KIND = idp), DIMENSION( 1:NUM_R_QUAD_POINTS,   &
                             1:NUM_T_QUAD_POINTS,   &
                             1:NUM_P_QUAD_POINTS)               ::  CUR_VAL_ALPHAPSI

REAL(KIND = idp), DIMENSION( 1:3,                   &
                             1:NUM_R_QUAD_POINTS,   &
                             1:NUM_T_QUAD_POINTS,   &
                             1:NUM_P_QUAD_POINTS )              ::  CUR_VAL_BETA

REAL(KIND = idp), DIMENSION( 1:3,                   &
                             1:NUM_R_QUAD_POINTS,   &
                             1:NUM_T_QUAD_POINTS,   &
                             1:NUM_P_QUAD_POINTS    )           ::  CUR_DRV_PSI

REAL(KIND = idp), DIMENSION( 1:3,                   &
                             1:NUM_R_QUAD_POINTS,   &
                             1:NUM_T_QUAD_POINTS,   &
                             1:NUM_P_QUAD_POINTS    )           ::  CUR_DRV_ALPHAPSI

REAL(KIND = idp), DIMENSION( 1:3, 1:3,              &
                             1:NUM_R_QUAD_POINTS,   &
                             1:NUM_T_QUAD_POINTS,   &
                             1:NUM_P_QUAD_POINTS    )           ::  CUR_DRV_BETA

REAL(KIND = idp), DIMENSION( 1:3, 1:3, 1:3,         &
                             1:NUM_R_QUAD_POINTS,   &
                             1:NUM_T_QUAD_POINTS,   &
                             1:NUM_P_QUAD_POINTS    )           ::  CUR_DDRV_BETA


REAL(KIND = idp), DIMENSION(1:NUM_R_QUAD_POINTS)                ::  CUR_R_LOCS
REAL(KIND = idp), DIMENSION(1:NUM_T_QUAD_POINTS)                ::  CUR_T_LOCS
REAL(KIND = idp), DIMENSION(1:NUM_P_QUAD_POINTS)                ::  CUR_P_LOCS


REAL(KIND = idp), DIMENSION(1:8)                                ::  PSI_POWER
REAL(KIND = idp), DIMENSION(1:4)                                ::  ALPHAPSI_POWER







COMPLEX(KIND = idp)                                             ::  Common_Term_A
REAL(KIND = idp)                                                ::  Common_Term_B,  &
                                                                    Common_Term_C,  &
                                                                    Common_Term_D



REAL(KIND = idp),  DIMENSION(1:5,                                    &
                                         1:NUM_R_QUAD_POINTS,        &
                                         1:NUM_T_QUAD_POINTS,        &
                                         1:NUM_P_QUAD_POINTS         )  ::  RHS_TERMS


REAL(KIND = idp),  DIMENSION(1:14,                                   &
                                         1:NUM_R_QUAD_POINTS,        &
                                         1:NUM_T_QUAD_POINTS,        &
                                         1:NUM_P_QUAD_POINTS         )  ::  SUBJCBN_PSI_TERMS

REAL(KIND = idp), DIMENSION(1:14,                                    &
                                         1:NUM_R_QUAD_POINTS,        &
                                         1:NUM_T_QUAD_POINTS,        &
                                         1:NUM_P_QUAD_POINTS         )  ::  SUBJCBN_ALPHAPSI_TERMS

REAL(KIND = idp), DIMENSION(1:20,                                    &
                                         1:NUM_R_QUAD_POINTS,        &
                                         1:NUM_T_QUAD_POINTS,        &
                                         1:NUM_P_QUAD_POINTS         )  ::  SUBJCBN_BETA1_TERMS

REAL(KIND = idp), DIMENSION(1:20,                                    &
                                         1:NUM_R_QUAD_POINTS,        &
                                         1:NUM_T_QUAD_POINTS,        &
                                         1:NUM_P_QUAD_POINTS         )  ::  SUBJCBN_BETA2_TERMS

REAL(KIND = idp), DIMENSION(1:20,                                    &
                                         1:NUM_R_QUAD_POINTS,        &
                                         1:NUM_T_QUAD_POINTS,        &
                                         1:NUM_P_QUAD_POINTS         )  ::  SUBJCBN_BETA3_TERMS



COMPLEX(KIND = idp)                                                     ::  REUSED_VALUE

REAL(KIND = idp), DIMENSION(1:3,1:3)                                    ::  JCBN_mu_Array
REAL(KIND = idp), DIMENSION(1:3)                                        ::  JCBN_n_ARRAY

REAL(KIND = idp)                                                        ::  JCBN_BIGK_VALUE



INTEGER, DIMENSION(1:NUM_T_QUAD_POINTS)                                 ::  here

REAL(KIND = idp), DIMENSION(1:NUM_R_QUAD_POINTS, 1:NUM_T_QUAD_POINTS)   ::  RSIN_SQUARE

REAL(KIND = idp), DIMENSION(1:NUM_R_QUAD_POINTS)                        ::  R_VAL,              &
                                                                            R_SQUARE,           &
                                                                            R_CUBED,             &
                                                                            R_INVERSE


REAL(KIND = idp), DIMENSION(1:NUM_T_QUAD_POINTS)                        ::  SIN_VAL,            &
                                                                            SIN_SQUARE,         &
                                                                            CSC_VAL,            &
                                                                            CSC_SQUARE,         &
                                                                            COTAN_VAL,          &
                                                                            COS_VAL,            &
                                                                            COS_SQUARE


COMPLEX(KIND = idp), DIMENSION(1:NUM_P_QUAD_POINTS)                     ::  PHI_EXP,            &
                                                                            PHI_TWOEXP



REAL(KIND = idp)                                                        ::  TWOOVER_DELTAR,    &
                                                                            deltar_overtwo,     &
                                                                            deltat_overtwo,     &
                                                                            deltap_overtwo

REAL(KIND = idp)                                                        ::  REAL_L



REAL(KIND = idp), DIMENSION( 1:NUM_R_QUAD_POINTS,        &
                             1:NUM_T_QUAD_POINTS,        &
                             1:NUM_P_QUAD_POINTS         )              ::  Int_Factor


INTEGER                                                                 ::  Block_T_Begin, Block_P_Begin




REAL(KIND = idp)                                    ::  timea, timeb, timec, timed, timee,  &
                                                        time_CURVALS, time_SJT, time_JCBNM, &
                                                        time_RHS


INTEGER                                             ::  ierr, i


time_CURVALS = 0.0_idp
time_SJT = 0.0_idp
time_JCBNM = 0.0_idp
time_RHS = 0.0_idp







Block_T_Begin = MOD(myID_Shell,NUM_BLOCK_THETA_ROWS)*NUM_T_ELEMS_PER_BLOCK
Block_P_Begin = (myID_Shell/NUM_BLOCK_THETA_ROWS)*NUM_P_ELEMS_PER_BLOCK





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


    R_INVERSE(:) = 1.0_idp/CUR_R_LOCS(:)



    !
    !   Move through phi
    !
    DO Local_pe = 0, NUM_P_ELEMS_PER_BLOCK-1



        Global_pe = Block_P_Begin + Local_pe

        deltap_overtwo = (plocs(Global_pe + 1) - plocs(Global_pe))/2.0_idp
        CUR_P_LOCS(:) = deltap_overtwo * (INT_P_LOCATIONS(:)+1.0_idp) + plocs(Global_pe)

        PHI_EXP(:) = EXP( CMPLX(0, -CUR_P_LOCS(:), KIND = idp) )
        PHI_TWOEXP(:) = EXP( CMPLX(0, -2.0_idp*CUR_P_LOCS(:), KIND = idp) )


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

                RSIN_SQUARE(:,td) = R_SQUARE(:)*SIN_SQUARE(td)

            END DO




            !*!
            !*! Calculate Current Values of CFA Varaiables and their Deriviatives
            !*!
            timea = MPI_Wtime()
            CALL Calc_3D_Current_Values(Global_re , Local_te  , Local_pe,       &
                                        CUR_VAL_PSI, CUR_DRV_PSI,                       &
                                        CUR_VAL_ALPHAPSI, CUR_DRV_ALPHAPSI,             &
                                        CUR_VAL_BETA, CUR_DRV_BETA, CUR_DDRV_BETA,      &
                                        Int_Factor,                                     &
                                        CUR_R_LOCS, CUR_T_LOCS, CUR_P_LOCS,             &
                                        R_SQUARE,                                       &
                                        SIN_VAL,                                        &
                                        DELTAR_OVERTWO, DELTAT_OVERTWO, DELTAP_OVERTWO  )


            timeb = MPI_Wtime()





            !*!
            !*!  Calculate the Sub-Jacobian and RHS Terms
            !*!
            CALL Calc_3D_SubJcbn_Terms( Local_re, Local_te, Local_pe,             &
                                        CUR_VAL_PSI, CUR_DRV_PSI,                       &
                                        CUR_VAL_ALPHAPSI, CUR_DRV_ALPHAPSI,             &
                                        CUR_VAL_BETA, CUR_DRV_BETA, CUR_DDRV_BETA,      &
                                        RHS_TERMS,                                      &
                                        SUBJCBN_PSI_TERMS, SUBJCBN_ALPHAPSI_TERMS,      &
                                        SUBJCBN_BETA1_TERMS, SUBJCBN_BETA2_TERMS,       &
                                        SUBJCBN_BETA3_TERMS,                            &
                                        CUR_R_LOCS, CUR_T_LOCS, CUR_P_LOCS,             &
                                        R_SQUARE, R_CUBED, R_INVERSE,                   &
                                        SIN_VAL, COS_VAL, CSC_VAL, COTAN_VAL,           &
                                        SIN_SQUARE, CSC_SQUARE, RSIN_SQUARE             )




            timec = MPI_Wtime()
            



            !*!
            !*! Create the Residual Vector ( Sans Laplacian Contribution )
            !*!
            CALL CREATE_3D_RHS_VECTOR(  Local_re, Local_te, Local_pe,           &
                                        RHS_TERMS,                              &
                                        R_SQUARE, deltar_OVERtwo,        &
                                        Int_Factor                              )


            timed = MPI_Wtime()


            !*!
            !*! Create Jacobian Matrix ( Sans Laplacian Contribution )
            !*!
            CALL CREATE_3D_JCBN_MATRIX( Local_re, Local_te, Local_pe,               &
                                        Global_re , Global_te  , Global_pe,         &
                                        TWOOVER_DELTAR,                             &
                                        SUBJCBN_PSI_TERMS, SUBJCBN_ALPHAPSI_TERMS,  &
                                        SUBJCBN_BETA1_TERMS, SUBJCBN_BETA2_TERMS,   &
                                        SUBJCBN_BETA3_TERMS,                        &
                                        Int_Factor                                  )



            timee = MPI_Wtime()







            !*!
            !*! Update Timer Values
            !*!
            time_CurVals = time_CurVals + timeb - timea
            time_SJT = time_SJT + timec - timeb
            time_RHS = time_RHS + timed - timec
            time_JCBNM = time_JCBNM + timee - timed



        END DO  ! re Loop
    END DO  ! te Loop
END DO  ! pe Loop


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
SUBROUTINE Calc_3D_Current_Values(  re, te, pe,                                     &
                                    CUR_VAL_PSI, CUR_DRV_PSI,                       &
                                    CUR_VAL_ALPHAPSI, CUR_DRV_ALPHAPSI,             &
                                    CUR_VAL_BETA, CUR_DRV_BETA, CUR_DDRV_BETA,      &
                                    Int_Factor,                                     &
                                    CUR_R_LOCS, CUR_T_LOCS, CUR_P_LOCS,             &
                                    R_SQUARE,                                       &
                                    SIN_VAL,                                        &
                                    DELTAR_OVERTWO, DELTAT_OVERTWO, DELTAP_OVERTWO  )

INTEGER, INTENT(IN)                                             ::  re, te, pe


REAL(KIND = idp), INTENT(INOUT), DIMENSION(1:NUM_R_QUAD_POINTS,   &
                                           1:NUM_T_QUAD_POINTS,   &
                                           1:NUM_P_QUAD_POINTS    )       ::  CUR_VAL_PSI

REAL(KIND = idp), INTENT(INOUT), DIMENSION(1:NUM_R_QUAD_POINTS,   &
                                           1:NUM_T_QUAD_POINTS,   &
                                           1:NUM_P_QUAD_POINTS    )       ::  CUR_VAL_ALPHAPSI


REAL(KIND = idp), INTENT(INOUT), DIMENSION(1:3,                   &
                                           1:NUM_R_QUAD_POINTS,   &
                                           1:NUM_T_QUAD_POINTS,   &
                                           1:NUM_P_QUAD_POINTS    )       ::  CUR_VAL_BETA


REAL(KIND = idp), INTENT(INOUT), DIMENSION(1:3,                   &
                                           1:NUM_R_QUAD_POINTS,   &
                                           1:NUM_T_QUAD_POINTS,   &
                                           1:NUM_P_QUAD_POINTS    )       ::  CUR_DRV_PSI

REAL(KIND = idp), INTENT(INOUT), DIMENSION(1:3,                   &
                                           1:NUM_R_QUAD_POINTS,   &
                                           1:NUM_T_QUAD_POINTS,   &
                                           1:NUM_P_QUAD_POINTS    )       ::  CUR_DRV_ALPHAPSI

REAL(KIND = idp), INTENT(INOUT), DIMENSION( 1:3, 1:3,              &
                                            1:NUM_R_QUAD_POINTS,   &
                                            1:NUM_T_QUAD_POINTS,   &
                                            1:NUM_P_QUAD_POINTS    )        ::  CUR_DRV_BETA

REAL(KIND = idp), INTENT(INOUT), DIMENSION( 1:3, 1:3, 1:3,         &
                                            1:NUM_R_QUAD_POINTS,   &
                                            1:NUM_T_QUAD_POINTS,   &
                                            1:NUM_P_QUAD_POINTS    )        ::  CUR_DDRV_BETA



REAL(KIND = idp), INTENT(INOUT), DIMENSION( 1:NUM_R_QUAD_POINTS,        &
                                            1:NUM_T_QUAD_POINTS,        &
                                            1:NUM_P_QUAD_POINTS         )   ::  Int_Factor




REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_R_QUAD_POINTS)                ::  CUR_R_LOCS
REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_T_QUAD_POINTS)                ::  CUR_T_LOCS
REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_P_QUAD_POINTS)                ::  CUR_P_LOCS



REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_R_QUAD_POINTS)            ::  R_SQUARE


REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_T_QUAD_POINTS)            ::  SIN_VAL


REAL(KIND = idp), INTENT(IN)                                            ::  DELTAR_OVERTWO,     &
                                                                            DELTAT_OVERTWO,     &
                                                                            DELTAP_OVERTWO





COMPLEX(KIND = idp), DIMENSION(1:5)                             ::  Tmp_U_Value,        &
                                                                    Tmp_U_R_DRV_Value,  &
                                                                    Tmp_U_T_DRV_Value,  &
                                                                    Tmp_U_P_DRV_Value,  &
                                                                    Tmp_U_RR_DDRV_Value


COMPLEX(KIND = idp), DIMENSION(1:5)                             ::  Tmp_BETA_TR_DDRV_Value,    &
                                                                    TMP_BETA_TT_DDRV_Value,    &
                                                                    Tmp_BETA_PR_DDRV_Value,    &
                                                                    Tmp_BETA_PT_DDRV_Value,    &
                                                                    Tmp_BETA_PP_DDRV_Value


INTEGER                                                         ::  l, m, d,    &
                                                                    rd, td, pd, &
                                                                    ui


INTEGER                                                         ::  Current_Location
INTEGER                                                         ::  lm_loc

COMPLEX(KIND = idp), DIMENSION(1:5)                             ::  Local_Coefficients




REAL(KIND = idp)                                                :: TWO_OVER_DELTAR


REAL(KIND = idp)                                                :: SQRT_NORM_TERM, REAL_L


                          !                                                 !
                         !!                                                 !!
                        !!!          Initialize Local Quadratures           !!!
                         !!                                                 !!
                          !                                                 !





TWO_OVER_DELTAR = 2.0_idp/(rlocs(re+1) - rlocs(re))


!$OMP PARALLEL DEFAULT(none)                                            &
!$OMP PRIVATE( pd, td, rd, l, m, d,                                     &
!$OMP           TMP_U_Value, Tmp_U_R_DRV_Value, Tmp_U_T_DRV_Value,      &
!$OMP           Tmp_U_P_DRV_Value, Tmp_U_RR_DDRV_Value,                 &
!$OMP           Tmp_Beta_TR_DDRV_Value, Tmp_Beta_TT_DDRV_Value,         &
!$OMP           Tmp_Beta_PR_DDRV_Value, Tmp_Beta_PT_DDRV_Value,         &
!$OMP           Tmp_Beta_PP_DDRV_Value,                                 &
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
!$OMP           TWO_OVER_DELTAR,                                        &
!$OMP           DELTAR_OVERTWO, DELTAT_OVERTWO, DELTAP_OVERTWO,         &
!$OMP           INT_R_WEIGHTS, INT_T_WEIGHTS, INT_P_WEIGHTS,            &
!$OMP           Lagrange_Poly_Table,                                    &
!$OMP           Ylm_Values, Ylm_dt_Values, Ylm_dp_Values,               &
!$OMP           Ylm_dtt_Values, Ylm_dtp_Values, Ylm_dpp_Values,         &
!$OMP           Ylm_CC_Values,                                          &
!$OMP           Coefficient_Vector,                                     &
!$OMP           Int_Factor,                                             &
!$OMP           Matrix_Location,                                        &
!$OMP           LM_Location,                                            &
!$OMP           LM_Length, ULM_LENGTH                           )









!$OMP DO SCHEDULE(dynamic), COLLAPSE(3)
DO pd = 1,NUM_P_QUAD_POINTS

    DO td = 1,NUM_T_QUAD_POINTS

        DO rd = 1,NUM_R_QUAD_POINTS

            


            Int_Factor(rd,td,pd) = SIN_VAL(td)*R_SQUARE(rd)                 &
                                    * DELTAR_OVERTWO * INT_R_WEIGHTS(rd)    &
                                    * DELTAT_OVERTWO * INT_T_WEIGHTS(td)    &
                                    * DELTAP_OVERTWO * INT_P_WEIGHTS(pd)


             !                                                  !
            !!   Set/Reset temporary value holder to zero.      !!
             !                                                  !
            Tmp_U_Value = 0.0_idp
            Tmp_U_R_DRV_Value = 0.0_idp
            Tmp_U_T_DRV_Value = 0.0_idp
            Tmp_U_P_DRV_Value = 0.0_idp
            TMP_U_RR_DDRV_Value = 0.0_idp
            Tmp_BETA_TR_DDRV_Value = 0.0_idp
            Tmp_BETA_TT_DDRV_Value = 0.0_idp
            Tmp_BETA_PR_DDRV_Value = 0.0_idp
            Tmp_BETA_PT_DDRV_Value = 0.0_idp
            Tmp_BETA_PP_DDRV_Value = 0.0_idp


            DO d = 0,DEGREE

                DO lm_loc = 0,LM_LENGTH-1



                    Current_Location = CFA_ALL_Matrix_Map(1, lm_loc, re, d)


                    Local_Coefficients(1:5) = Coefficient_Vector(Current_Location:Current_Location+4)


                    TMP_U_Value(:) = TMP_U_Value(:)                                         &
                                        + Local_Coefficients(:)                             &
                                                * Ylm_Values(lm_loc, td, pd, te, pe)        &
                                                * Lagrange_Poly_Table(d, rd, 0)



                    Tmp_U_R_DRV_Value(:) = Tmp_U_R_DRV_Value(:)                     &
                                            + Local_Coefficients(:)                 &
                                                * Ylm_Values(lm_loc, td, pd, te, pe)&
                                                * Lagrange_Poly_Table(d, rd, 1)     &
                                                * TWO_OVER_DELTAR

 


                    Tmp_U_T_DRV_Value(:) = Tmp_U_T_DRV_Value(:)                             &
                                            + Local_Coefficients(:)                         &
                                                * Ylm_dt_Values(lm_loc, td, pd, te, pe)     &
                                                * Lagrange_Poly_Table(d,rd,0)


                    Tmp_U_P_DRV_Value(:) = Tmp_U_P_DRV_Value(:)                     &
                                            + Local_Coefficients(:)                 &
                                                * Ylm_dp_Values(lm_loc, td, pd, te, pe)     &
                                                * Lagrange_Poly_Table(d, rd, 0)




                    Tmp_U_RR_DDRV_Value(:) = Tmp_U_RR_DDRV_Value(:)                 &
                                            + Local_Coefficients(:)                 &
                                                * Ylm_Values(lm_loc, td, pd, te, pe)&
                                                * Lagrange_Poly_Table(d, rd, 2)     &
                                                * TWO_OVER_DELTAR                   &
                                                * TWO_OVER_DELTAR



                    TMP_BETA_TR_DDRV_VALUE(3:5) = TMP_BETA_TR_DDRV_VALUE(3:5)               &
                                                + Local_Coefficients(3:5)                   &
                                                    * Ylm_dt_Values(lm_loc, td, pd, te, pe) &
                                                    * Lagrange_Poly_Table(d,rd,1)           &
                                                    * TWO_OVER_DELTAR

                    TMP_BETA_TT_DDRV_VALUE(3:5) = TMP_BETA_TT_DDRV_VALUE(3:5)       &
                                                + Local_Coefficients(3:5)           &
                                                    * Ylm_dtt_Values(lm_loc, td, pd, te, pe)&
                                                    * Lagrange_Poly_Table(d,rd,0)

                    TMP_BETA_PR_DDRV_VALUE(3:5) = TMP_BETA_PR_DDRV_VALUE(3:5)       &
                                                + Local_Coefficients(3:5)           &
                                                    * Ylm_dp_Values(lm_loc, td, pd, te, pe)&
                                                    * Lagrange_Poly_Table(d, rd, 1) &
                                                    * TWO_OVER_DELTAR

                    TMP_BETA_PT_DDRV_VALUE(3:5) = TMP_BETA_PT_DDRV_VALUE(3:5)       &
                                                + Local_Coefficients(3:5)           &
                                                    * Ylm_dtp_Values(lm_loc, td, pd, te, pe) &
                                                    * Lagrange_Poly_Table(d,rd,0)


                    TMP_BETA_PP_DDRV_VALUE(3:5) = TMP_BETA_PP_DDRV_VALUE(3:5)       &
                                                + Local_Coefficients(3:5)           &
                                                    * Ylm_dpp_Values(lm_loc, td, pd, te, pe)&
                                                    * Lagrange_Poly_Table(d, rd, 0)


                END DO  !   lm_loc Loop
            END DO  !   d Loop


 


            CUR_VAL_PSI(rd, td, pd)          = REAL(Tmp_U_Value(1), KIND = idp)
            CUR_DRV_PSI(1, rd, td, pd )      = REAL(Tmp_U_R_DRV_Value(1), KIND = idp)
            CUR_DRV_PSI(2, rd, td, pd )      = REAL(Tmp_U_T_DRV_Value(1), KIND = idp)
            CUR_DRV_PSI(3, rd, td, pd )      = REAL(Tmp_U_P_DRV_Value(1), KIND = idp)



            CUR_VAL_ALPHAPSI(rd, td, pd)     = REAL(Tmp_U_Value(2), KIND = idp)
            CUR_DRV_ALPHAPSI(1, rd, td, pd ) = REAL(Tmp_U_R_DRV_Value(2), KIND = idp)
            CUR_DRV_ALPHAPSI(2, rd, td, pd ) = REAL(Tmp_U_T_DRV_Value(2), KIND = idp)
            CUR_DRV_ALPHAPSI(3, rd, td, pd ) = REAL(Tmp_U_P_DRV_Value(2), KIND = idp)




            CUR_VAL_BETA(1, rd, td, pd )     = REAL(Tmp_U_Value(3), KIND = idp)
            CUR_VAL_BETA(2, rd, td, pd )     = REAL(Tmp_U_Value(4), KIND = idp)
            CUR_VAL_BETA(3, rd, td, pd )     = REAL(Tmp_U_Value(5), KIND = idp)



            CUR_DRV_BETA(1, 1, rd, td, pd )  = REAL(Tmp_U_R_DRV_Value(3), KIND = idp)
            CUR_DRV_BETA(2, 1, rd, td, pd )  = REAL(Tmp_U_T_DRV_Value(3), KIND = idp)
            CUR_DRV_BETA(3, 1, rd, td, pd )  = REAL(Tmp_U_P_DRV_Value(3), KIND = idp)
            CUR_DRV_BETA(1, 2, rd, td, pd )  = REAL(Tmp_U_R_DRV_Value(4), KIND = idp)
            CUR_DRV_BETA(2, 2, rd, td, pd )  = REAL(Tmp_U_T_DRV_Value(4), KIND = idp)
            CUR_DRV_BETA(3, 2, rd, td, pd )  = REAL(Tmp_U_P_DRV_Value(4), KIND = idp)
            CUR_DRV_BETA(1, 3, rd, td, pd )  = REAL(Tmp_U_R_DRV_Value(5), KIND = idp)
            CUR_DRV_BETA(2, 3, rd, td, pd )  = REAL(Tmp_U_T_DRV_Value(5), KIND = idp)
            CUR_DRV_BETA(3, 3, rd, td, pd )  = REAL(Tmp_U_P_DRV_Value(5), KIND = idp)


            CUR_DDRV_BETA(1, 1, 1, rd, td, pd )   = REAL(Tmp_U_RR_DDRV_Value(3), KIND = idp)
            CUR_DDRV_BETA(2, 1, 1, rd, td, pd )   = REAL(Tmp_BETA_TR_DDRV_Value(3), KIND = idp)
            CUR_DDRV_BETA(3, 1, 1, rd, td, pd )   = REAL(Tmp_BETA_PR_DDRV_Value(3), KIND = idp)
            CUR_DDRV_BETA(1, 2, 1, rd, td, pd )   = REAL(Tmp_BETA_TR_DDRV_Value(3), KIND = idp)
            CUR_DDRV_BETA(2, 2, 1, rd, td, pd )   = REAL(Tmp_BETA_TT_DDRV_Value(3), KIND = idp)
            CUR_DDRV_BETA(3, 2, 1, rd, td, pd )   = REAL(Tmp_BETA_PT_DDRV_Value(3), KIND = idp)
            CUR_DDRV_BETA(1, 3, 1, rd, td, pd )   = REAL(Tmp_BETA_PR_DDRV_Value(3), KIND = idp)
            CUR_DDRV_BETA(2, 3, 1, rd, td, pd )   = REAL(Tmp_BETA_PT_DDRV_Value(3), KIND = idp)
            CUR_DDRV_BETA(3, 3, 1, rd, td, pd )   = REAL(Tmp_BETA_PP_DDRV_Value(3), KIND = idp)
            CUR_DDRV_BETA(1, 1, 2, rd, td, pd )   = REAL(Tmp_U_RR_DDRV_Value(4), KIND = idp)
            CUR_DDRV_BETA(2, 1, 2, rd, td, pd )   = REAL(Tmp_BETA_TR_DDRV_Value(4), KIND = idp)
            CUR_DDRV_BETA(3, 1, 2, rd, td, pd )   = REAL(Tmp_BETA_PR_DDRV_Value(4), KIND = idp)
            CUR_DDRV_BETA(1, 2, 2, rd, td, pd )   = REAL(Tmp_BETA_TR_DDRV_Value(4), KIND = idp)
            CUR_DDRV_BETA(2, 2, 2, rd, td, pd )   = REAL(Tmp_BETA_TT_DDRV_Value(4), KIND = idp)
            CUR_DDRV_BETA(3, 2, 2, rd, td, pd )   = REAL(Tmp_BETA_PT_DDRV_Value(4), KIND = idp)
            CUR_DDRV_BETA(1, 3, 2, rd, td, pd )   = REAL(Tmp_BETA_PR_DDRV_Value(4), KIND = idp)
            CUR_DDRV_BETA(2, 3, 2, rd, td, pd )   = REAL(Tmp_BETA_PT_DDRV_Value(4), KIND = idp)
            CUR_DDRV_BETA(3, 3, 2, rd, td, pd )   = REAL(Tmp_BETA_PP_DDRV_Value(4), KIND = idp)
            CUR_DDRV_BETA(1, 1, 3, rd, td, pd )   = REAL(Tmp_U_RR_DDRV_Value(5), KIND = idp)
            CUR_DDRV_BETA(2, 1, 3, rd, td, pd )   = REAL(Tmp_BETA_TR_DDRV_Value(5), KIND = idp)
            CUR_DDRV_BETA(3, 1, 3, rd, td, pd )   = REAL(Tmp_BETA_PR_DDRV_Value(5), KIND = idp)
            CUR_DDRV_BETA(1, 2, 3, rd, td, pd )   = REAL(Tmp_BETA_TR_DDRV_Value(5), KIND = idp)
            CUR_DDRV_BETA(2, 2, 3, rd, td, pd )   = REAL(Tmp_BETA_TT_DDRV_Value(5), KIND = idp)
            CUR_DDRV_BETA(3, 2, 3, rd, td, pd )   = REAL(Tmp_BETA_PT_DDRV_Value(5), KIND = idp)
            CUR_DDRV_BETA(1, 3, 3, rd, td, pd )   = REAL(Tmp_BETA_PR_DDRV_Value(5), KIND = idp)
            CUR_DDRV_BETA(2, 3, 3, rd, td, pd )   = REAL(Tmp_BETA_PT_DDRV_Value(5), KIND = idp)
            CUR_DDRV_BETA(3, 3, 3, rd, td, pd )   = REAL(Tmp_BETA_PP_DDRV_Value(5), KIND = idp)


        END DO  !   rd Loop

    END DO  !   td Loop

END DO  !   pd Loop
!$OMP END DO

!$OMP END PARALLEL





END SUBROUTINE Calc_3D_Current_Values











!+203+###########################################################################!
!                                                                                !
!                  Calc_3D_SubJcbn_Terms                                         !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_3D_SubJcbn_Terms(   re, te, pe,                                     &
                                    CUR_VAL_PSI, CUR_DRV_PSI,                       &
                                    CUR_VAL_ALPHAPSI, CUR_DRV_ALPHAPSI,             &
                                    CUR_VAL_BETA, CUR_DRV_BETA, CUR_DDRV_BETA,      &
                                    RHS_TERMS,                                      &
                                    SUBJCBN_PSI_TERMS, SUBJCBN_ALPHAPSI_TERMS,      &
                                    SUBJCBN_BETA1_TERMS, SUBJCBN_BETA2_TERMS,       &
                                    SUBJCBN_BETA3_TERMS,                            &
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


REAL(KIND = idp), INTENT(INOUT), DIMENSION( 1:14,                       &
                                            1:NUM_R_QUAD_POINTS,        &
                                            1:NUM_T_QUAD_POINTS,        &
                                            1:NUM_P_QUAD_POINTS         )   ::  SUBJCBN_PSI_TERMS

REAL(KIND = idp), INTENT(INOUT), DIMENSION( 1:14,                       &
                                            1:NUM_R_QUAD_POINTS,        &
                                            1:NUM_T_QUAD_POINTS,        &
                                            1:NUM_P_QUAD_POINTS         )   ::  SUBJCBN_ALPHAPSI_TERMS

REAL(KIND = idp), INTENT(INOUT), DIMENSION( 1:20,                       &
                                            1:NUM_R_QUAD_POINTS,        &
                                            1:NUM_T_QUAD_POINTS,        &
                                            1:NUM_P_QUAD_POINTS         )   ::  SUBJCBN_BETA1_TERMS

REAL(KIND = idp), INTENT(INOUT), DIMENSION( 1:20,                       &
                                            1:NUM_R_QUAD_POINTS,        &
                                            1:NUM_T_QUAD_POINTS,        &
                                            1:NUM_P_QUAD_POINTS         )   ::  SUBJCBN_BETA2_TERMS

REAL(KIND = idp), INTENT(INOUT), DIMENSION( 1:20,                       &
                                            1:NUM_R_QUAD_POINTS,        &
                                            1:NUM_T_QUAD_POINTS,        &
                                            1:NUM_P_QUAD_POINTS         )   ::  SUBJCBN_BETA3_TERMS


REAL(KIND = idp)                                                            ::  REUSED_VALUE

INTEGER                                                                     ::  pd, td, rd,     &
                                                                                i


REAL(KIND = idp), DIMENSION(1:8)                                            ::  PSI_POWER
REAL(KIND = idp), DIMENSION(1:4)                                            ::  ALPHAPSI_POWER


REAL(KIND = idp), DIMENSION(1:3,1:3)                                        ::  JCBN_mu_Array
REAL(KIND = idp), DIMENSION(1:3)                                            ::  JCBN_n_ARRAY

REAL(KIND = idp)                                                            ::  JCBN_BIGK_VALUE












!$OMP PARALLEL DEFAULT(none)                                            &
!$OMP PRIVATE( pd, td, rd,                                              &
!$OMP           PSI_POWER, ALPHAPSI_POWER,                              &
!$OMP           JCBN_BIGK_VALUE, JCBN_mu_ARRAY,JCBN_n_ARRAY,            &
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




        JCBN_mu_Array = JCBN_mu_FUNCTION_3D_ALL( rd, td, pd,                                    &
                                                CUR_R_LOCS(rd), R_SQUARE(rd), R_CUBED(rd),      &
                                                R_INVERSE(rd), RSIN_SQUARE(rd, td),             &
                                                SIN_VAL(td), SIN_SQUARE(td), CSC_SQUARE(td),    &
                                                COS_VAL(td), COTAN_VAL(td),                     &
                                                CUR_VAL_BETA, CUR_DRV_BETA                      )




        JCBN_n_ARRAY(:) = CUR_DRV_ALPHAPSI(:, rd, td, pd) / ALPHAPSI_POWER(1)   &
                            - 7 * CUR_DRV_PSI(:, rd, td, pd )/ PSI_POWER(1)









        RHS_Terms(1, rd, td, pd) = - TwoPi * GR_Source_Scalar * Block_Source_E(rd, td, pd, re, te, pe) * PSI_POWER(5)            &
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
                            + COTAN_VAL(td) * CUR_DRV_BETA(1, 2, rd, td, pd)                    &
                            - 2.0_idp /R_SQUARE(rd) * CUR_VAL_BETA(1, rd, td, pd )              &
                            )                                                                   &
                      + JCBN_n_ARRAY(1)* JCBN_mu_Array(1,1)                                     &
                      + JCBN_n_ARRAY(2)* JCBN_mu_Array(2,1)                                     &
                      + JCBN_n_ARRAY(3)* JCBN_mu_Array(3,1)





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
                               - CSC_SQUARE(td) * CUR_VAL_BETA(2, rd, td, pd)     &
                             )                                                                  &
                        + JCBN_n_ARRAY(1) * JCBN_mu_Array(1,2)                                  &
                        + JCBN_n_ARRAY(2) * JCBN_mu_Array(2,2)                                  &
                        + JCBN_n_ARRAY(3) * JCBN_mu_Array(3,2)





        RHS_Terms(5, rd, td, pd) = 16.0_idp * pi                                            &
                          * ALPHAPSI_POWER(1)                                               &
                          * PSI_POWER(3)                                                    &
                          * GR_Source_Scalar                                                &
                          * Block_Source_Si(rd, td, pd, re, te, pe, 3)                      &
                        - 2.0_idp/(CUR_R_LOCS(rd)*RSIN_SQUARE(rd,td)) * CUR_DRV_BETA(3,1,rd,td,pd) & ! NEW
                        - COTAN_VAL(td)/RSIN_SQUARE(rd,td) * CUR_DRV_BETA(3,2,rd,td,pd)     &  ! NEW
                        - 2.0_idp/CUR_R_LOCS(rd) * CUR_DRV_BETA(1,3,rd,td,pd)               &  ! NEW
                        - 2.0_idp*COTAN_VAL(td) / R_SQUARE(rd) * CUR_DRV_BETA(2,3,rd,td,pd) &  ! NEW
                        - ( OneThird*CSC_SQUARE(td)/R_SQUARE(rd ))                          &
                             * ( CUR_DDRV_BETA(3, 1, 1, rd, td, pd )                        &
                               + CUR_DDRV_BETA(3, 2, 2, rd, td, pd )                        &
                               + CUR_DDRV_BETA(3, 3, 3, rd, td, pd )                        &
                               + 2.0_idp/CUR_R_LOCS(rd) * CUR_DRV_BETA(3, 1, rd, td, pd)    &
                               + COTAN_VAL(td) * CUR_DRV_BETA(3, 2, rd, td, pd)             &
                             )                                                              &
                        + JCBN_n_ARRAY(1) * JCBN_mu_Array(1,3)                              &
                        + JCBN_n_ARRAY(2) * JCBN_mu_Array(2,3)                              &
                        + JCBN_n_ARRAY(3) * JCBN_mu_Array(3,3)





        !
        !   With Respect to Psi Terms Jacobian Terms
        !

        ! J_{1,1lmn}
        SUBJCBN_PSI_TERMS(1, rd, td, pd) = 10.0_idp*pi                                                  &
                                                * GR_Source_Scalar                                      &
                                                * PSI_POWER(4)                                          &
                                                * Block_Source_E(rd, td, pd, re, te, pe)                &
                                         + (7.0_idp/16.0_idp)                                           &
                                                * PSI_POWER(6)/ALPHAPSI_POWER(2)                        &
                                                * JCBN_BIGK_VALUE

        ! J_{2,1lmn}
        SUBJCBN_PSI_TERMS(2, rd, td, pd) = -8.0_idp * pi * ALPHAPSI_POWER(1) * PSI_POWER(3)             &
                                            * GR_Source_Scalar                                          &
                                            * ( Block_Source_E(rd, td, pd, re, te, pe)                  &
                                                + 2 * Block_Source_S(rd, td, pd, re, te, pe) )          &
                                            - (42.0_idp / 16.0_idp )                                    &
                                                * PSI_POWER(5)/ALPHAPSI_POWER(1)                        &
                                                * JCBN_BIGK_VALUE



        ! J_{3,1lmn} Non-Derivative Term
        SUBJCBN_PSI_TERMS(3, rd, td, pd) = (-7.0_idp / PSI_POWER(2) )                                   &
                                                * (CUR_DRV_PSI(1, rd, td, pd )*JCBN_mu_ARRAY(1,1)       &
                                                  + CUR_DRV_PSI(2, rd, td, pd )*JCBN_mu_ARRAY(2,1)      &
                                                  + CUR_DRV_PSI(3, rd, td, pd )*JCBN_mu_ARRAY(3,1) )    &
                                            - 48.0_idp * pi                                             &
                                                * GR_Source_Scalar                                      &
                                                * Block_Source_Si(rd, td, pd, re, te, pe, 1)            &
                                                * ALPHAPSI_POWER(1)                                     &
                                                * PSI_POWER(2)


        ! J_{3,1lmn} Derivative Terms
        SUBJCBN_PSI_TERMS(4:6, rd, td, pd) = (7.0_idp/ PSI_POWER(1) )*JCBN_mu_ARRAY(1:3,1)



        ! J_{4,1lmn} Non-Derivative Term
        SUBJCBN_PSI_TERMS(7, rd, td, pd) = (-7.0_idp / PSI_POWER(2) )                                   &
                                                * (CUR_DRV_PSI(1, rd, td, pd)*JCBN_mu_ARRAY(1,2)        &
                                                 + CUR_DRV_PSI(2, rd, td, pd)*JCBN_mu_ARRAY(2,2)        &
                                                 + CUR_DRV_PSI(3, rd, td, pd)*JCBN_mu_ARRAY(3,2) )      &
                                            - 48.0_idp * pi                                             &
                                                * GR_Source_Scalar                                      &
                                                * Block_Source_Si(rd, td, pd, re, te, pe, 2)            &
                                                * ALPHAPSI_POWER(1)                                     &
                                                * PSI_POWER(2)

        ! J_{4,1lmn} Derivative Terms
        SUBJCBN_PSI_TERMS(8:10, rd, td, pd) = (7.0_idp/ PSI_POWER(1) )*JCBN_mu_ARRAY(1:3,2)






        ! J_{5,1lmn} Non-Derivative Term
        SUBJCBN_PSI_TERMS(11, rd, td, pd) = (-7.0_idp / PSI_POWER(2) )                  &
                                * (CUR_DRV_PSI(1, rd, td, pd)*JCBN_mu_ARRAY(1,3)        &
                                 + CUR_DRV_PSI(2, rd, td, pd)*JCBN_mu_ARRAY(2,3)        &
                                 + CUR_DRV_PSI(3, rd, td, pd)*JCBN_mu_ARRAY(3,3) )      &
                            - 48.0_idp * pi                                             &
                                * GR_Source_Scalar                                      &
                                * Block_Source_Si(rd, td, pd, re, te, pe, 3)            &
                                * ALPHAPSI_POWER(1)                                     &
                                * PSI_POWER(2)

        ! J_{5,1lmn} Derivative Terms
        SUBJCBN_PSI_TERMS(12:14, rd, td, pd) = (7.0_idp/ PSI_POWER(1) )*JCBN_mu_ARRAY(1:3,3)






        !
        ! Alpha Psi Terms
        !


        ! J_{1,2lmn}
        SUBJCBN_ALPHAPSI_TERMS(1, rd, td, pd) = -PSI_POWER(7)/(8.0_idp*ALPHAPSI_POWER(3))               &
                                                * JCBN_BIGK_VALUE

        ! J_{2,2lmn}
        SUBJCBN_ALPHAPSI_TERMS(2, rd, td, pd) =  -TwoPi * PSI_POWER(4)                                      &
                                                * GR_Source_Scalar                                          &
                                                *( Block_Source_E(rd, td, pd, re, te, pe)                   &
                                                    + 2.0_idp * Block_Source_S(rd, td, pd, re, te, pe) )    &
                                                + (7.0_idp / 16.0_idp)                                      &
                                                    * PSI_POWER(5)/ALPHAPSI_POWER(2)                        &
                                                    * JCBN_BIGK_VALUE

        ! J_{3,2lmn} Non-Derivative Term
        SUBJCBN_ALPHAPSI_TERMS(3, rd, td, pd) = (CUR_DRV_ALPHAPSI(1, rd, td, pd )*JCBN_mu_Array(1,1)        &
                                                + CUR_DRV_ALPHAPSI(2, rd, td, pd )*JCBN_mu_Array(2,1)       &
                                                + CUR_DRV_ALPHAPSI(3, rd, td, pd )*JCBN_mu_Array(3,1)   )   &
                                                /ALPHAPSI_POWER(2)                                          &
                                                -16.0_idp * pi                                              &
                                                    * GR_Source_Scalar                                      &
                                                    * Block_Source_Si(rd,td,pd,re,te,pe,1)                  &
                                                    * PSI_POWER(3)


        ! J_{3,2lmn} Derivative Terms
        SUBJCBN_ALPHAPSI_TERMS(4:6, rd, td, pd) = (-1.0_idp/ALPHAPSI_POWER(1) )*JCBN_mu_ARRAY(1:3,1)




        ! J_{4,2lmn} Non-Derivative Term
        SUBJCBN_ALPHAPSI_TERMS(7, rd, td, pd) = ( CUR_DRV_ALPHAPSI(1, rd, td, pd )*JCBN_mu_Array(1,2)       &
                                                + CUR_DRV_ALPHAPSI(2, rd, td, pd )*JCBN_mu_Array(2,2)       &
                                                + CUR_DRV_ALPHAPSI(3, rd, td, pd )*JCBN_mu_Array(3,2)   )   &
                                                /ALPHAPSI_POWER(2)                                          &
                                                -16.0_idp * pi                                              &
                                                    * GR_Source_Scalar                                      &
                                                    * Block_Source_Si(rd, td, pd, re, te, pe, 2)            &
                                                    * PSI_POWER(3)


        ! J_{4,2lmn} Derivative Terms
        SUBJCBN_ALPHAPSI_TERMS(8:10, rd, td, pd) = (-1.0_idp/ ALPHAPSI_POWER(1) )*JCBN_mu_ARRAY(1:3,2)



        ! J_{5,2lmn} Non-Derivative Term
        SUBJCBN_ALPHAPSI_TERMS(11, rd, td, pd) = (CUR_DRV_ALPHAPSI(1, rd, td, pd ) * JCBN_mu_Array(1,3)     &
                                                 + CUR_DRV_ALPHAPSI(2, rd, td, pd ) * JCBN_mu_Array(2,3)    &
                                                 + CUR_DRV_ALPHAPSI(3, rd, td, pd ) * JCBN_mu_Array(3,3)  ) &
                                                /ALPHAPSI_POWER(2)                                          &
                                                -16.0_idp * pi                                              &
                                                    * GR_Source_Scalar                                      &
                                                    * Block_Source_Si(rd,td,pd,re,te,pe,3)                  &
                                                    * PSI_POWER(3)

        ! J_{5,2lmn} Derivative Terms
        SUBJCBN_ALPHAPSI_TERMS(12:14, rd, td, pd) = (-1.0_idp/ ALPHAPSI_POWER(1) )*JCBN_mu_ARRAY(1:3,3)






        !
        !   Beta 1 Terms
        !




        REUSED_VALUE = PSI_POWER(7)/(16.0_idp* ALPHAPSI_POWER(2) )


        ! J_{1,3lmn}
        ! Reused_Value * kappa_{10}
        SUBJCBN_BETA1_TERMS(1, rd, td, pd) = REUSED_VALUE * FourThirds                               &
                                * ( 2.0_idp / R_SQUARE(rd)           * CUR_VAL_BETA(1,rd,td,pd)      &
                                + COTAN_VAL(td) /CUR_R_LOCS(rd)      * CUR_VAL_BETA(2,rd,td,pd)      &
                                - 2.0_idp/CUR_R_LOCS(rd)             * CUR_DRV_BETA(1,1,rd,td,pd)    &
                                + 1.0_idp/CUR_R_LOCS(rd)             * CUR_DRV_BETA(2,2,rd,td,pd)    &
                                + 1.0_idp/CUR_R_LOCS(rd)             * CUR_DRV_BETA(3,3,rd,td,pd)    )
        ! J_{1,3lmn}
        ! Reused_Value * kappa_{11}
        SUBJCBN_BETA1_TERMS(2, rd, td, pd) = REUSED_VALUE * FourThirds                              &
                                * ( - 2.0_idp /CUR_R_LOCS(rd)       * CUR_VAL_BETA(1,rd,td,pd)      &
                                - COTAN_VAL(td)                     * CUR_VAL_BETA(2,rd,td,pd)      &
                                + 2.0_idp                           * CUR_DRV_BETA(1,1,rd,td,pd)    &
                                - 1.0_idp                           * CUR_DRV_BETA(2,2,rd,td,pd)    &
                                - 1.0_idp                           * CUR_DRV_BETA(3,3,rd,td,pd)    )
        ! J_{1,3lmn}
        ! Reused_Value * kappa_{12}
        SUBJCBN_BETA1_TERMS(3, rd, td, pd) = REUSED_VALUE                                           &
                                * ( 2.0_idp / R_SQUARE(rd)      * CUR_DRV_BETA(2,1,rd,td,pd)        &
                                + 2.0_idp                       * CUR_DRV_BETA(1,2,rd,td,pd)    )
        ! J_{1,3lmn}
        ! Reused_Value * kappa_{13}
        SUBJCBN_BETA1_TERMS(4, rd, td, pd) = REUSED_VALUE                                           &
                                * ( 2.0_idp*CSC_SQUARE(td)/R_SQUARE(rd)   * CUR_DRV_BETA(3,1,rd,td,pd)    &
                                + 2.0_idp                                 * CUR_DRV_BETA(1,3,rd,td,pd)    )


        ! J_{2,3lmn}
        SUBJCBN_BETA1_TERMS(5:8, rd, td, pd) = -7.0_idp * ( ALPHAPSI_POWER(1)/PSI_POWER(1) )        &
                                                    * SUBJCBN_BETA1_TERMS(1:4, rd, td, pd)



        ! J_{3,3lmn} Non-Derivative Term
        SUBJCBN_BETA1_TERMS(9, rd, td, pd) = (4.0_idp/(3.0_idp *CUR_R_LOCS(rd))) * JCBN_n_ARRAY(1)      &
                                            - ( 8.0_idp/(3.0_idp * R_SQUARE(rd) ) )

        ! J_{3,3lmn} Derivative Terms
        SUBJCBN_BETA1_TERMS(10, rd, td, pd) = OneThird * ((2.0_idp/CUR_R_LOCS(rd)) - 4.0_idp * JCBN_n_ARRAY(1) )
        SUBJCBN_BETA1_TERMS(11, rd, td, pd) = -(1.0_idp/R_SQUARE(rd)) * JCBN_n_ARRAY(2)
        SUBJCBN_BETA1_TERMS(12, rd, td, pd) = -(CSC_SQUARE(td)/R_SQUARE(rd)) * JCBN_n_ARRAY(3)

        ! J_{3,3lmn} Double Derivative Term



        ! J_{4,3lmn} Non-Derivative Term
        SUBJCBN_BETA1_TERMS(13, rd, td, pd) = -(TwoThirds/R_CUBED(rd)) * JCBN_n_ARRAY(2)

        ! J_{4,3lmn} Derivative Terms
        SUBJCBN_BETA1_TERMS(14, rd, td, pd) = - (TwoThirds / R_SQUARE(rd)) * JCBN_n_ARRAY(2)
        SUBJCBN_BETA1_TERMS(15, rd, td, pd) = (8.0_idp / (3.0_idp * R_CUBED(rd))) - JCBN_n_ARRAY(1)*(1.0_idp/R_SQUARE(rd))
        SUBJCBN_BETA1_TERMS(16, rd, td, pd) = 0.0_idp


        ! J_{5,3lmn} Non-Derivative Term
        SUBJCBN_BETA1_TERMS(17, rd, td, pd) = -(TwoThirds*CSC_SQUARE(td)/R_CUBED(rd))* JCBN_n_ARRAY(3)

        ! J_{5,3lmn} Derivative Terms
        SUBJCBN_BETA1_TERMS(18, rd, td, pd) = -(TwoThirds*CSC_SQUARE(td)/R_SQUARE(rd)) * JCBN_n_ARRAY(3)
        SUBJCBN_BETA1_TERMS(19, rd, td, pd) = 0.0_idp
        SUBJCBN_BETA1_TERMS(20, rd, td, pd) = (OneThird*CSC_SQUARE(td) / R_SQUARE(rd)) * (( 8.0_idp /CUR_R_LOCS(rd) )             &
                                                + 3.0_idp * JCBN_n_ARRAY(1) )


        !
        !   Beta 2 Terms
        !




        ! J_{1,4lmn}
        ! Reused_Value * kappa_{20}
        SUBJCBN_BETA2_TERMS(1, rd, td, pd) = REUSED_VALUE * FourThirds * COTAN_VAL(td)                      &
                                            * ( 1.0_idp /CUR_R_LOCS(rd)     * CUR_VAL_BETA(1,rd,td,pd)      &
                                              + 2.0_idp                     * CUR_VAL_BETA(2,rd,td,pd)      &
                                              - 1.0_idp                     * CUR_DRV_BETA(2,2,rd,td,pd)    &
                                              + 2.0_idp                     * CUR_DRV_BETA(3,3,rd,td,pd)    )
        ! J_{1,4lmn}
        ! Reused_Value * kappa_{21}
        SUBJCBN_BETA2_TERMS(2, rd, td, pd) = REUSED_VALUE                                                   &
                                            * ( 2.0_idp * R_SQUARE(rd)      * CUR_DRV_BETA(1,2,rd,td,pd)    &
                                                + 2.0_idp                   * CUR_DRV_BETA(2,1,rd,td,pd)    )
        ! J_{1,4lmn}
        ! Reused_Value * kappa_{22}
        SUBJCBN_BETA2_TERMS(3, rd, td, pd) = REUSED_VALUE * FourThirds                                      &
                                        * ( 1.0_idp /CUR_R_LOCS(rd)         * CUR_VAL_BETA(1,rd,td,pd)      &
                                            - COTAN_VAL(td)                 * CUR_VAL_BETA(2,rd,td,pd)      &
                                            - 1.0_idp                       * CUR_DRV_BETA(1,1,rd,td,pd)    &
                                            + 2.0_idp                       * CUR_DRV_BETA(2,2,rd,td,pd)    &
                                            - 1.0_idp                       * CUR_DRV_BETA(3,3,rd,td,pd)    )
            ! J_{1,4lmn}
        ! Reused_Value * kappa_{23}
        SUBJCBN_BETA2_TERMS(4, rd, td, pd) = REUSED_VALUE                                                   &
                                            * ( 2.0_idp * CSC_SQUARE(td)    * CUR_DRV_BETA(3,2,rd,td,pd)    &
                                                + 2.0_idp                   * CUR_DRV_BETA(2,3,rd,td,pd)    )


        ! J_{2,4lmn}
        SUBJCBN_BETA2_TERMS(5:8, rd, td, pd) = -7.0_idp * ( ALPHAPSI_POWER(1)/PSI_POWER(1) )                &
                                                * SUBJCBN_BETA2_TERMS(1:4, rd, td, pd)



        ! J_{3,4lmn} Non-Derivative Term
        SUBJCBN_BETA2_TERMS(9, rd, td, pd) = -2.0_idp * COTAN_VAL(td) / CUR_R_LOCS(rd)                      &
                                             + TwoThirds*COTAN_VAL(td) * JCBN_n_ARRAY(1)

        ! J_{3,4lmn} Derivative Terms
        SUBJCBN_BETA2_TERMS(10, rd, td, pd) = OneThird * COTAN_VAL(td) - JCBN_n_ARRAY(2)
        SUBJCBN_BETA2_TERMS(11, rd, td, pd) = -2.0_idp/CUR_R_LOCS(rd) + TwoThirds * JCBN_n_ARRAY(1)
        SUBJCBN_BETA2_TERMS(12, rd, td, pd) = 0.0_idp


        ! J_{4,4lmn} Non-Derivative Term
        SUBJCBN_BETA2_TERMS(13, rd, td, pd) = (TwoThirds / R_SQUARE(rd)) * COTAN_VAL(td) * JCBN_n_ARRAY(2)               &
                                                - OneThird * (2.0_idp - 4.0_idp*COTAN_VAL(td)*COTAN_VAL(td))/R_SQUARE(rd)

        ! J_{4,4lmn} Derivative Terms
        SUBJCBN_BETA2_TERMS(14, rd, td, pd) = -JCBN_n_ARRAY(1) + 2.0_idp/CUR_R_LOCS(rd)
        SUBJCBN_BETA2_TERMS(15, rd, td, pd) = (OneThird / R_SQUARE(rd) )* ( COTAN_VAL(td) - 4.0_idp * JCBN_n_ARRAY(2) )
        SUBJCBN_BETA2_TERMS(16, rd, td, pd) = - (1.0_idp*CSC_SQUARE(td)/R_SQUARE(rd)) * JCBN_n_ARRAY(3)


        ! J_{5,4lmn} Non-Derivative Term
        SUBJCBN_BETA2_TERMS(17, rd, td, pd) = -FourThirds*(COTAN_VAL(td)*CSC_SQUARE(td)/R_SQUARE(rd))* JCBN_n_ARRAY(3)

        ! J_{5,4lmn} Derivative Terms
        SUBJCBN_BETA2_TERMS(18, rd, td, pd) = 0.0_idp
        SUBJCBN_BETA2_TERMS(19, rd, td, pd) = (TwoThirds*CSC_SQUARE(td) / R_SQUARE(rd)) * JCBN_n_ARRAY(3)
        SUBJCBN_BETA2_TERMS(20, rd, td, pd) = (OneThird*CSC_SQUARE(td) / R_SQUARE(rd))                                      &
                                            * ( 7.0_idp*COTAN_VAL(td) + 3.0_idp * JCBN_n_ARRAY(2)  )







        !
        !   Beta 3 Terms
        !




        ! Reused_Value * kappa_{30}
        SUBJCBN_BETA3_TERMS(1, rd, td, pd) = 0.0_idp

        ! Reused_Value * kappa_{31}
        SUBJCBN_BETA3_TERMS(2, rd, td, pd) = REUSED_VALUE                                           &
                                * ( 2.0_idp * RSIN_SQUARE(rd, td)   * CUR_DRV_BETA(1,3,rd,td,pd)    &
                                + 2.0_idp                   * CUR_DRV_BETA(3,1,rd,td,pd)    )

        ! Reused_Value * kappa_{32}
        SUBJCBN_BETA3_TERMS(3, rd, td, pd) = REUSED_VALUE                                           &
                                * ( 2.0_idp * SIN_SQUARE(td)    * CUR_DRV_BETA(2,3,rd,td,pd)        &
                                + 2.0_idp                   * CUR_DRV_BETA(3,2,rd,td,pd)    )

        ! Reused_Value * kappa_{33}
        SUBJCBN_BETA3_TERMS(4, rd, td, pd) = REUSED_VALUE * FourThirds                              &
                                * ( 1.0_idp /CUR_R_LOCS(rd)         * CUR_VAL_BETA(1,rd,td,pd)      &
                                + 2.0_idp * COTAN_VAL(td)       * CUR_VAL_BETA(2,rd,td,pd)          &
                                - 1.0_idp                   * CUR_DRV_BETA(1,1,rd,td,pd)            &
                                - 1.0_idp                   * CUR_DRV_BETA(2,2,rd,td,pd)            &
                                + 2.0_idp                   * CUR_DRV_BETA(3,3,rd,td,pd)    )


        ! J_{2,5lmn}
        SUBJCBN_BETA3_TERMS(5:8, rd, td, pd) = -7.0_idp * ( ALPHAPSI_POWER(1)/PSI_POWER(1) )        &
                                                * SUBJCBN_BETA3_TERMS(1:4, rd, td, pd)










        ! J_{3,5lmn} Non-Derivative Term
        SUBJCBN_BETA3_TERMS(9, rd, td, pd) =  0.0_idp

        ! J_{3,5lmn} Derivative Terms
        SUBJCBN_BETA3_TERMS(10, rd, td, pd) = -JCBN_n_ARRAY(3)
        SUBJCBN_BETA3_TERMS(11, rd, td, pd) = 0.0_idp
        SUBJCBN_BETA3_TERMS(12, rd, td, pd) = -2.0_idp/CUR_R_LOCS(rd) - TwoThirds * JCBN_n_ARRAY(1)


        ! J_{4,5lmn} Non-Derivative Term
        SUBJCBN_BETA3_TERMS(13, rd, td, pd) = 0.0_idp

        ! J_{4,5lmn} Derivative Terms
        SUBJCBN_BETA3_TERMS(14, rd, td, pd) = 0.0_idp
        SUBJCBN_BETA3_TERMS(15, rd, td, pd) = (-1.0_idp/ R_SQUARE(rd)) * JCBN_n_ARRAY(3)
        SUBJCBN_BETA3_TERMS(16, rd, td, pd) = (TwoThirds / R_SQUARE(rd)) * JCBN_n_ARRAY(2)- 2.0_idp*COTAN_VAL(td)/R_SQUARE(rd)


        ! J_{5,5lmn} Non-Derivative Term
        SUBJCBN_BETA3_TERMS(17, rd, td, pd) = 0.0_idp
        ! J_{5,5lmn} Derivative Terms
        SUBJCBN_BETA3_TERMS(18, rd, td, pd) = -JCBN_n_ARRAY(1) + 2.0_idp/CUR_R_LOCS(rd)
        SUBJCBN_BETA3_TERMS(19, rd, td, pd) = (-1.0_idp/R_SQUARE(rd)) * JCBN_n_ARRAY(2) + 2.0_idp*COTAN_VAL(td)/R_SQUARE(rd)
        SUBJCBN_BETA3_TERMS(20, rd, td, pd) = -(FourThirds*CSC_SQUARE(td) / R_SQUARE(rd)) * JCBN_n_ARRAY(3)







        END DO ! pd loop
    END DO  ! td loop
END DO  ! rd loop


!$OMP END DO

!$OMP END PARALLEL






END SUBROUTINE Calc_3D_SubJcbn_Terms

















!+204+###########################################################################!
!                                                                                !
!                  CREATE_3D_RHS_VECTOR                                          !
!                                                                                !
!################################################################################!
SUBROUTINE CREATE_3D_RHS_VECTOR(    re, te, pe,                             &
                                    RHS_TERMS,                              &
                                    R_SQUARE, DELTAR_OVERTWO,       &
                                    Int_Factor                              )



INTEGER, INTENT(IN)                                                         ::  re, te, pe

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:5,                        &
                                         1:NUM_R_QUAD_POINTS,        &
                                         1:NUM_T_QUAD_POINTS,        &
                                         1:NUM_P_QUAD_POINTS         )      ::  RHS_TERMS





REAL(KIND = idp), INTENT(IN), DIMENSION( 1:NUM_R_QUAD_POINTS,        &
                                         1:NUM_T_QUAD_POINTS,        &
                                         1:NUM_P_QUAD_POINTS         )      ::  Int_Factor


REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_R_QUAD_POINTS)            ::  R_SQUARE
REAL(KIND = idp), INTENT(IN)                                            ::  DELTAR_OVERTWO

INTEGER                                                                 ::  pd, td, rd,     &
                                                                            l, m, d,        &
                                                                            lm_loc

INTEGER                                                                 ::  Current_i_Location

COMPLEX(KIND = idp), DIMENSION(1:5)                                     ::  RHS_TMP
COMPLEX(KIND = idp)                      :: Test





!$OMP PARALLEL DEFAULT(none)                                            &
!$OMP PRIVATE( pd, td, rd, d, l, m,                                     &
!$OMP           lm_loc,                                                 &
!$OMP           Current_i_Location, RHS_TMP, TEST                         )   &
!$OMP SHARED( re, te, pe,                                                &
!$OMP           DEGREE,                                                 &
!$OMP           R_SQUARE, Deltar_Overtwo, Int_R_Weights,                &
!$OMP           NUM_P_QUAD_POINTS, NUM_T_QUAD_POINTS, NUM_R_QUAD_POINTS,&
!$OMP           Int_Factor,                                             &
!$OMP           ULM_LENGTH, LM_LENGTH,                                  &
!$OMP           Ylm_CC_Values, Lagrange_Poly_Table,                     &
!$OMP           LM_Location,                                            &
!$OMP           Block_RHS_Vector, RHS_Terms                             )

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



END SUBROUTINE CREATE_3D_RHS_VECTOR






















!+205+###########################################################################!
!                                                                                !
!                  CREATE_3D_JCBN_MATRIX          !
!                                                                                !
!################################################################################!
SUBROUTINE CREATE_3D_JCBN_MATRIX( Local_re, Local_te, Local_pe,                   &
                                    Global_re, Global_te, Global_pe,                &
                                    TWOOVER_DELTAR,                                 &
                                    SUBJCBN_PSI_TERMS, SUBJCBN_ALPHAPSI_TERMS,      &
                                    SUBJCBN_BETA1_TERMS, SUBJCBN_BETA2_TERMS,       &
                                    SUBJCBN_BETA3_TERMS,                            &
                                    Int_Factor                                      )



INTEGER, INTENT(IN)                                                     ::  Local_re, Local_te, Local_pe,       &
                                                                            Global_re, Global_te, Global_pe



REAL(KIND = idp), INTENT(IN)                                            ::  TWOOVER_DELTAR





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




INTEGER                                                                 ::  lm_loc, lpmp_loc

INTEGER                                                                 ::  MATVEC_LOC
INTEGER                                                                 ::  iloc, jloc


REAL(KIND = idp), DIMENSION(1:NUM_R_QUAD_POINTS)                        ::  Common_Term_A,  &
                                                                            Common_Term_B,  &
                                                                            Common_Term_C



COMPLEX(KIND = idp), DIMENSION(1:25)                                    ::  Uncommon_Terms


REAL(KIND = idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:2)                   ::  Current_Lag_Polys

COMPLEX(KIND = idp), DIMENSION(1:6)                                     ::  Ylm_Cross_Terms



INTEGER                                                                 ::  i, ierr








!$OMP PARALLEL DEFAULT(none)                                            &
!$OMP PRIVATE( pd, td, rd,  d, dp, i, ierr,                             &
!$OMP           Current_i_Location, Current_j_Location,                 &
!$OMP           MATVEC_LOC, iloc, jloc,                                 &
!$OMP           lm_loc, lpmp_loc,                                       &
!$OMP           Common_Term_A, Common_Term_B, Common_Term_C,            &
!$OMP           Current_Lag_Polys,                                      &
!$OMP           Ylm_Cross_Terms,                                        &
!$OMP           Uncommon_Terms                                      )   &
!$OMP SHARED(   Local_re, Local_te, Local_pe,                           &
!$OMP           DEGREE,                                                 &
!$OMP           NUM_P_QUAD_POINTS, NUM_T_QUAD_POINTS, NUM_R_QUAD_POINTS,&
!$OMP           Global_re, Global_te, Global_pe,                        &
!$OMP           TWOOVER_DELTAR,                                         &
!$OMP           Ylm_Values, Ylm_dt_Values, Ylm_dp_Values,               &
!$OMP           Ylm_dtt_Values, Ylm_dtp_Values, Ylm_dpp_Values,         &
!$OMP           Ylm_CC_Values,                                          &
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
!$OMP           POSEIDON_COMM_WORLD, myID_Poseidon,                     &
!$OMP           LM_Length, ULM_Length                                   )



!$OMP DO SCHEDULE(dynamic), COLLAPSE(3)      
DO d = 0,DEGREE

    DO dp = 0,DEGREE

         DO lm_loc = 0,LM_Length-1

            Current_Lag_Polys(:,0) = LPT_LPT(:, d, dp, 0, 0 )
            Current_Lag_Polys(:,1) = LPT_LPT(:, d, dp, 0, 1 )            &
                                    * TWOOVER_DELTAR
            Current_Lag_Polys(:,2) = LPT_LPT(:, d, dp, 0, 2 )            &
                                    * TWOOVER_DELTAR                     &
                                    * TWOOVER_DELTAR


            jloc = d*ULM_LENGTH + NUM_CFA_VARS*lm_loc

            DO lpmp_loc = 0,LM_Length - 1


                iloc = dp*ULM_LENGTH + NUM_CFA_VARS * lpmp_loc

                Uncommon_Terms = 0.0_idp


                DO pd = 1,NUM_P_QUAD_POINTS
                    DO td = 1,NUM_T_QUAD_POINTS


                        Common_Term_A(:) = Current_Lag_Polys(:,0) * Int_Factor(:, td, pd)
                        Common_Term_B(:) = Current_Lag_Polys(:,1) * Int_Factor(:, td, pd)
                        Common_Term_C(:) = Current_Lag_Polys(:,2) * Int_Factor(:, td, pd)


                        Ylm_Cross_Terms(1) = Ylm_Values(lm_loc, td, pd, Local_te, Local_pe)       &
                                           *  Ylm_CC_Values(lm_loc, td, pd, Local_te, Local_pe)
                        Ylm_Cross_Terms(2) = Ylm_dt_Values(lm_loc, td, pd, Local_te, Local_pe)    &
                                           *  Ylm_CC_Values(lm_loc, td, pd, Local_te, Local_pe)
                        Ylm_Cross_Terms(3) = Ylm_dp_Values(lm_loc, td, pd, Local_te, Local_pe)    &
                                           *  Ylm_CC_Values(lm_loc, td, pd, Local_te, Local_pe)
                        Ylm_Cross_Terms(4) = Ylm_dtt_Values(lm_loc, td, pd, Local_te, Local_pe)   &
                                           *  Ylm_CC_Values(lm_loc, td, pd, Local_te, Local_pe)
                        Ylm_Cross_Terms(5) = Ylm_dtp_Values(lm_loc, td, pd, Local_te, Local_pe)   &
                                           *  Ylm_CC_Values(lm_loc, td, pd, Local_te, Local_pe)
                        Ylm_Cross_Terms(6) = Ylm_dpp_Values(lm_loc, td, pd, Local_te, Local_pe)   &
                                           *  Ylm_CC_Values(lm_loc, td, pd, Local_te, Local_pe)


                        Uncommon_Terms(1) = Uncommon_Terms(1)                                                   &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_PSI_TERMS(1, :, td, pd), Common_Term_A(:))



                        Uncommon_Terms(2) = Uncommon_Terms(2)                                                   &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT( SUBJCBN_ALPHAPSI_TERMS(1, :, td, pd), Common_Term_A(:) )


                        Uncommon_Terms(3) = Uncommon_Terms(3)                                                   &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA1_TERMS(1, :, td, pd), Common_Term_A(:))    &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA1_TERMS(2, :, td, pd), Common_Term_B(:))    &
                                          + Ylm_Cross_Terms(2)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA1_TERMS(3, :, td, pd), Common_Term_A(:))    &
                                          + Ylm_Cross_Terms(3)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA1_TERMS(4, :, td, pd), Common_Term_A(:))


                        Uncommon_Terms(4) = Uncommon_Terms(4)                                                   &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA2_TERMS(1, :, td, pd), Common_Term_A(:))    &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA2_TERMS(2, :, td, pd), Common_Term_B(:))    &
                                          + Ylm_Cross_Terms(2)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA2_TERMS(3, :, td, pd), Common_Term_A(:))    &
                                          + Ylm_Cross_Terms(3)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA2_TERMS(4, :, td, pd), Common_Term_A(:))

                        Uncommon_Terms(5) = Uncommon_Terms(5)                                                   &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA3_TERMS(1, :, td, pd), Common_Term_A(:))    &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA3_TERMS(2, :, td, pd), Common_Term_B(:))    &
                                          + Ylm_Cross_Terms(2)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA3_TERMS(3, :, td, pd), Common_Term_A(:))    &
                                          + Ylm_Cross_Terms(3)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA3_TERMS(4, :, td, pd), Common_Term_A(:))









                        Uncommon_Terms(6) = Uncommon_Terms(6)                                                   &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_PSI_TERMS(2, :, td, pd), Common_Term_A(:))



                        Uncommon_Terms(7) = Uncommon_Terms(7)                                                   &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_ALPHAPSI_TERMS(2, :, td, pd), Common_Term_A(:))

                        Uncommon_Terms(8) = Uncommon_Terms(8)                                                   &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA1_TERMS(5, :, td, pd), Common_Term_A(:))    &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA1_TERMS(6, :, td, pd), Common_Term_B(:))    &
                                          + Ylm_Cross_Terms(2)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA1_TERMS(7, :, td, pd), Common_Term_A(:))    &
                                          + Ylm_Cross_Terms(3)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA1_TERMS(8, :, td, pd), Common_Term_A(:))


                        Uncommon_Terms(9) = Uncommon_Terms(9)                                                   &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA2_TERMS(5, :, td, pd), Common_Term_A(:))    &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA2_TERMS(6, :, td, pd), Common_Term_B(:))    &
                                          + Ylm_Cross_Terms(2)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA2_TERMS(7, :, td, pd), Common_Term_A(:))    &
                                          + Ylm_Cross_Terms(3)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA2_TERMS(8, :, td, pd), Common_Term_A(:))

                        Uncommon_Terms(10) = Uncommon_Terms(10)                                                 &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA3_TERMS(5, :, td, pd), Common_Term_A(:))    &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA3_TERMS(6, :, td, pd), Common_Term_B(:))    &
                                          + Ylm_Cross_Terms(2)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA3_TERMS(7, :, td, pd), Common_Term_A(:))    &
                                          + Ylm_Cross_Terms(3)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA3_TERMS(8, :, td, pd), Common_Term_A(:))







                        Uncommon_Terms(11) = Uncommon_Terms(11)                                             &
                                          + Ylm_Cross_Terms(1)                                              &
                                          * DOT_PRODUCT(SUBJCBN_PSI_TERMS(3, :, td, pd), Common_Term_A(:))  &
                                          + Ylm_Cross_Terms(1)                                              &
                                          * DOT_PRODUCT(SUBJCBN_PSI_TERMS(4, :, td, pd), Common_Term_B(:))  &
                                          + Ylm_Cross_Terms(2)                                              &
                                          * DOT_PRODUCT(SUBJCBN_PSI_TERMS(5, :, td, pd), Common_Term_A(:))  &
                                          + Ylm_Cross_Terms(3)                                              &
                                          * DOT_PRODUCT(SUBJCBN_PSI_TERMS(6, :, td, pd), Common_Term_A(:))

                        Uncommon_Terms(12) = Uncommon_Terms(12)                                                 &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_ALPHAPSI_TERMS(3, :, td, pd), Common_Term_A(:)) &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_ALPHAPSI_TERMS(4, :, td, pd), Common_Term_B(:)) &
                                          + Ylm_Cross_Terms(2)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_ALPHAPSI_TERMS(5, :, td, pd), Common_Term_A(:)) &
                                          + Ylm_Cross_Terms(3)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_ALPHAPSI_TERMS(6, :, td, pd), Common_Term_A(:))


                        Uncommon_Terms(13) = Uncommon_Terms(13)                                                 &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA1_TERMS(9, :, td, pd), Common_Term_A(:))    &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA1_TERMS(10, :, td, pd), Common_Term_B(:))   &
                                          + Ylm_Cross_Terms(2)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA1_TERMS(11, :, td, pd), Common_Term_A(:))   &
                                          + Ylm_Cross_Terms(3)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA1_TERMS(12, :, td, pd), Common_Term_A(:))   &
                                          + OneThird * Ylm_Cross_Terms(1) * SUM( Common_Term_C )


                        Uncommon_Terms(14) = Uncommon_Terms(14)                                                 &
                                          + Ylm_Cross_Terms(1)                                      &
                                          * DOT_PRODUCT(SUBJCBN_BETA2_TERMS(9, :, td, pd), Common_Term_A(:))    &
                                          + Ylm_Cross_Terms(1)                                      &
                                          * DOT_PRODUCT(SUBJCBN_BETA2_TERMS(10, :, td, pd), Common_Term_B(:))   &
                                          + Ylm_Cross_Terms(2)                                      &
                                          * DOT_PRODUCT(SUBJCBN_BETA2_TERMS(11, :, td, pd), Common_Term_A(:))   &
                                          + Ylm_Cross_Terms(3)                                      &
                                          * DOT_PRODUCT(SUBJCBN_BETA2_TERMS(12, :, td, pd), Common_Term_A(:))   &
                                          + OneThird * Ylm_Cross_Terms(2) * SUM( Common_Term_B )



                        Uncommon_Terms(15) = Uncommon_Terms(15)                                                 &
                                          + Ylm_Cross_Terms(1)                                      &
                                          * DOT_PRODUCT(SUBJCBN_BETA3_TERMS(9, :, td, pd), Common_Term_A(:))    &
                                          + Ylm_Cross_Terms(2)                                     &
                                          * DOT_PRODUCT(SUBJCBN_BETA3_TERMS(11, :, td, pd), Common_Term_A(:))   &
                                          + Ylm_Cross_Terms(3)                                      &
                                          * DOT_PRODUCT(SUBJCBN_BETA3_TERMS(12, :, td, pd), Common_Term_A(:))   &
                                          + OneThird * Ylm_Cross_Terms(3) * SUM(Common_Term_B)










                        Uncommon_Terms(16) = Uncommon_Terms(16)                                                 &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_PSI_TERMS(7, :, td, pd), Common_Term_A(:))      &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_PSI_TERMS(8, :, td, pd), Common_Term_B(:))      &
                                          + Ylm_Cross_Terms(2)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_PSI_TERMS(9, :, td, pd), Common_Term_A(:))      &
                                          + Ylm_Cross_Terms(3)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_PSI_TERMS(10, :, td, pd), Common_Term_A(:))

                        Uncommon_Terms(17) = Uncommon_Terms(17)                                                 &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_ALPHAPSI_TERMS(7, :, td, pd), Common_Term_A(:)) &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_ALPHAPSI_TERMS(8, :, td, pd), Common_Term_B(:)) &
                                          + Ylm_Cross_Terms(2)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_ALPHAPSI_TERMS(9, :, td, pd), Common_Term_A(:)) &
                                          + Ylm_Cross_Terms(3)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_ALPHAPSI_TERMS(10, :, td, pd), Common_Term_A(:))


                        Uncommon_Terms(18) = Uncommon_Terms(18)                                                 &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA1_TERMS(13, :, td, pd), Common_Term_A(:))   &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA1_TERMS(14, :, td, pd), Common_Term_B(:))   &
                                          + Ylm_Cross_Terms(3)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA1_TERMS(16, :, td, pd), Common_Term_A(:))   &
                                          + OneThird * Ylm_Cross_Terms(2) * SUM( Common_Term_B )


                        Uncommon_Terms(19) = Uncommon_Terms(19)                                                 &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA2_TERMS(13, :, td, pd), Common_Term_A(:))   &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA2_TERMS(14, :, td, pd), Common_Term_B(:))   &
                                          + Ylm_Cross_Terms(2)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA2_TERMS(15, :, td, pd), Common_Term_A(:))   &
                                          + Ylm_Cross_Terms(3)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA2_TERMS(16, :, td, pd), Common_Term_A(:))   &
                                          + OneThird * Ylm_Cross_Terms(4)  * SUM( Common_Term_A )

                        Uncommon_Terms(20) = Uncommon_Terms(20)                                                 &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA3_TERMS(13, :, td, pd), Common_Term_A(:))   &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA3_TERMS(14, :, td, pd), Common_Term_B(:))   &
                                          + Ylm_Cross_Terms(2)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA3_TERMS(15, :, td, pd), Common_Term_A(:))   &
                                          + OneThird * Ylm_Cross_Terms(5) * SUM( Common_Term_A )








                        Uncommon_Terms(21) = Uncommon_Terms(21)                                                 &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_PSI_TERMS(11, :, td, pd), Common_Term_A(:))     &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_PSI_TERMS(12, :, td, pd), Common_Term_B(:))     &
                                          + Ylm_Cross_Terms(2)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_PSI_TERMS(13, :, td, pd), Common_Term_A(:))     &
                                          + Ylm_Cross_Terms(3)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_PSI_TERMS(14, :, td, pd), Common_Term_A(:))

                        Uncommon_Terms(22) = Uncommon_Terms(22)                                                 &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_ALPHAPSI_TERMS(11, :, td, pd), Common_Term_A(:))&
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_ALPHAPSI_TERMS(12, :, td, pd), Common_Term_B(:))&
                                          + Ylm_Cross_Terms(2)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_ALPHAPSI_TERMS(13, :, td, pd), Common_Term_A(:))&
                                          + Ylm_Cross_Terms(3)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_ALPHAPSI_TERMS(14, :, td, pd), Common_Term_A(:))


                        Uncommon_Terms(23) = Uncommon_Terms(23)                                                 &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA1_TERMS(17, :, td, pd), Common_Term_A(:))   &
                                          + Ylm_Cross_Terms(2)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA1_TERMS(19, :, td, pd), Common_Term_A(:))   &
                                          + Ylm_Cross_Terms(3)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA1_TERMS(20, :, td, pd), Common_Term_A(:))   &
                                          + OneThird * Ylm_Cross_Terms(3) * SUM( Common_Term_B )


                        Uncommon_Terms(24) = Uncommon_Terms(24)                                                 &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA2_TERMS(17, :, td, pd), Common_Term_A(:))   &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA2_TERMS(18, :, td, pd), Common_Term_B(:))   &
                                          + Ylm_Cross_Terms(2)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA2_TERMS(19, :, td, pd), Common_Term_A(:))   &
                                          + OneThird * Ylm_Cross_Terms(4) * SUM( Common_Term_A )

                        Uncommon_Terms(25) = Uncommon_Terms(25)                                                 &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA3_TERMS(17, :, td, pd), Common_Term_A(:))   &
                                          + Ylm_Cross_Terms(1)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA3_TERMS(18, :, td, pd), Common_Term_B(:))   &
                                          + Ylm_Cross_Terms(2)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA3_TERMS(19, :, td, pd), Common_Term_A(:))   &
                                          + Ylm_Cross_Terms(3)                                                  &
                                          * DOT_PRODUCT(SUBJCBN_BETA3_TERMS(20, :, td, pd), Common_Term_A(:))   &
                                          + OneThird * Ylm_Cross_Terms(6) * SUM( Common_Term_A )






                    END DO  ! td Loop
                END DO  ! pd Loop


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



INTEGER                                             ::  ui, l, m, re, d, dp, rd

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
    !$OMP PRIVATE( ui, l, m, re, d, dp, rd,                                 &
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

!                    Reusable_Values(dp) = SUM ( ( -R_SQUARE(:) * LPT_LPT(:,dp,d,1,1)       &
!                                                + L_Lp1 * LPT_LPT(:,dp,d,0,0)        )    &
!                                                * INT_R_WEIGHTS(:) * TWOOVER_DELTAR     )


                    Reusable_Values(dp) = SUM ( ( -R_SQUARE(:) * LPT_LPT(:,dp,d,1,1)* TWOOVER_DELTAR* TWOOVER_DELTAR      &
                                                + L_Lp1 * LPT_LPT(:,dp,d,0,0)        )    &
                                                * INT_R_WEIGHTS(:)/TWOOVER_DELTAR     )




                END DO


                DO m = -M_VALUES(l),M_VALUES(l)

                   jloc = d*ULM_LENGTH + NUM_CFA_VARS *(l*(l+1) + m)


                    DO dp = 0,DEGREE


                        iloc = dp*ULM_LENGTH + NUM_CFA_VARS *(l*(l+1) + m)


                        ! Psi Laplacian Terms
                        MATVEC_LOC = jloc*ELEM_PROB_DIM + iloc
                        BLOCK_ELEM_STF_MATVEC(MATVEC_LOC,Block_re)              &
                                = BLOCK_ELEM_STF_MATVEC(MATVEC_LOC,Block_re)    &
                                + Reusable_Values(dp)


                        ! Alpha Psi Laplacian Terms 
                        MATVEC_LOC = (jloc+1)*ELEM_PROB_DIM + iloc + 1
                        BLOCK_ELEM_STF_MATVEC(MATVEC_LOC,Block_re)              &
                                = BLOCK_ELEM_STF_MATVEC(MATVEC_LOC,Block_re)    &
                                + Reusable_Values(dp)



                        ! Beta 1 Laplacian Terms 
                        MATVEC_LOC = (jloc+2)*ELEM_PROB_DIM + iloc + 2
                        BLOCK_ELEM_STF_MATVEC(MATVEC_LOC,Block_re)              &
                                = BLOCK_ELEM_STF_MATVEC(MATVEC_LOC,Block_re)    &
                                + Reusable_Values(dp)


                        ! Beta 2 Laplacian Terms
                        MATVEC_LOC = (jloc+3)*ELEM_PROB_DIM + iloc + 3
                        BLOCK_ELEM_STF_MATVEC(MATVEC_LOC,Block_re)              &
                                = BLOCK_ELEM_STF_MATVEC(MATVEC_LOC,Block_re)    &
                                + Reusable_Values(dp)


                        !  Beta 3 Laplacian Terms
                        MATVEC_LOC = (jloc+4)*ELEM_PROB_DIM + iloc + 4
                        BLOCK_ELEM_STF_MATVEC(MATVEC_LOC,Block_re)              &
                                = BLOCK_ELEM_STF_MATVEC(MATVEC_LOC,Block_re)    &
                                + Reusable_Values(dp)




                    END DO ! dp Loop

                END DO ! m Loop

            END DO ! l Loop

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


INTEGER                                                         ::  ui, l, m, re, d, dp, rd

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
                                                                    End_Here


COMPLEX(KIND = idp), DIMENSION(0:SUBSHELL_PROB_DIM-1)                ::  TMP_VECTOR


INTEGER :: Start_Here_b, End_Here_b

INTEGER                                                         ::  ierr

REAL(KIND = idp)           ::  timea, timeb, timec



INTEGER :: i,j


IF ( POSEIDON_COMM_PETSC .NE. MPI_COMM_NULL ) THEN

    Block_STF_MAT = 0.0_idp


    timea = MPI_Wtime()

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

            DO l = 0,L_LIMIT



                L_Lp1 = REAL( l*(l+1), idp )




                DO dp = 0,DEGREE


                    Reusable_Values(dp) = SUM ( ( -R_SQUARE(:) * LPT_LPT(:,dp,d,1,1)* TWOOVER_DELTAR* TWOOVER_DELTAR      &
                                                + L_Lp1 * LPT_LPT(:,dp,d,0,0)        )    &
                                                * INT_R_WEIGHTS(:)/TWOOVER_DELTAR     )



                END DO


                DO m = -M_VALUES(l),M_VALUES(l)



                    Current_i_Location = Matrix_Location(1, l, m, re, d)


                    DO dp = 0,DEGREE

                        Current_j_Location = NUM_OFF_DIAGONALS + (dp - d) * ULM_LENGTH


                        Block_STF_MAT(Current_j_Location, Current_i_Location )          &
                                = Block_STF_MAT(Current_j_Location, Current_i_Location )      &
                                    + Reusable_Values(dp)

                        Block_STF_MAT(Current_j_Location, Current_i_Location + 1)       &
                                = Block_STF_MAT(Current_j_Location, Current_i_Location + 1)   &
                                    + Reusable_Values(dp)


                        Block_STF_MAT(Current_j_Location, Current_i_Location + 2 )      &
                                = Block_STF_MAT(Current_j_Location, Current_i_Location + 2 )  &
                                    + Reusable_Values(dp)

                        Block_STF_MAT(Current_j_Location, Current_i_Location + 3 )      &
                                = Block_STF_MAT(Current_j_Location, Current_i_Location + 3 )  &
                                    + Reusable_Values(dp)

                        Block_STF_MAT(Current_j_Location, Current_i_Location + 4 )      &
                                = Block_STF_MAT(Current_j_Location, Current_i_Location + 4 )  &
                                    + Reusable_Values(dp)



                    END DO ! dp Loop

                END DO ! m Loop

            END DO ! l Loop

        END DO ! d Loop
        !$OMP END DO

    END DO ! re Loop
!    !$OMP END DO

    !$OMP END PARALLEL


    timeb=MPI_Wtime()


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
    Block_RHS_Vector(Start_Here:End_Here) = Block_RHS_Vector(Start_Here:End_Here)   &
                                          + TMP_VECTOR(0:SUBSHELL_PROB_DIM-1)

    timec=MPI_Wtime()

 
END IF

END SUBROUTINE FINISH_3D_RHS_VECTOR









!+601+###########################################################################!
!                                                                                !
!                  CFA_3D_Apply_BCs                                                  !
!                                                                                !
!################################################################################!
SUBROUTINE CFA_3D_Apply_BCs_Part1(  )

CALL CFA_3D_Dirichlet_BCs_Part1()
CALL CFA_3D_Neumann_BCs()

END SUBROUTINE CFA_3D_Apply_BCs_Part1



!+602+###########################################################################!
!                                                                                !
!                  CFA_3D_Apply_BCs                                                  !
!                                                                                !
!################################################################################!
SUBROUTINE CFA_3D_Apply_BCs_Part2(  )

CALL CFA_3D_Dirichlet_BCs_Part2()
CALL CFA_3D_Neumann_BCs()

END SUBROUTINE CFA_3D_Apply_BCs_Part2


!+603+###########################################################################!
!                                                                                !
!                  CFA_3D_Dirichlet_BCs                                                  !
!                                                                                !
!################################################################################!
SUBROUTINE CFA_3D_Dirichlet_BCs_Part1( )

INTEGER                                             ::  i, l, m, ui, d

INTEGER                                             ::  Value_Location

!*!
!*!     Steps 1
!*!
DO ui = 1,NUM_CFA_VARS


    !*!
    !*!     Innner BCs
    !*!
    IF ( INNER_CFA_BC_TYPE(ui) == 'D' ) THEN

        DO l = 0,L_LIMIT
            DO m = -M_VALUES(l),M_VALUES(l)

                Value_Location =  MAtrix_Location( ui, l, m, 0, 0 )
                Coefficient_Vector( Value_Location ) = 0.0_idp

            END DO
        END DO

        Value_Location =  MAtrix_Location( ui, 0, 0, 0, 0 )
        Coefficient_Vector( Value_Location ) = 2.0_idp * sqrt(pi) * INNER_CFA_BC_VALUES(ui)
    END IF





    !*!
    !*!     Outer BCs
    !*!
    IF ( OUTER_CFA_BC_TYPE(ui) == 'D' ) THEN

        DO l = 1,L_LIMIT
            DO m = -M_VALUES(l),M_VALUES(l)

                Value_Location =  Matrix_Location( ui, l, m, NUM_R_ELEMENTS-1, DEGREE )
                Coefficient_Vector( Value_Location ) = 0.0_idp

            END DO
        END DO

        Value_Location =  MAtrix_Location( ui, 0, 0, NUM_R_ELEMENTS-1, DEGREE  )
        Coefficient_Vector( Value_Location ) = 2.0_idp * sqrt(pi) * OUTER_CFA_BC_VALUES(ui)

    END IF

END DO


END SUBROUTINE CFA_3D_Dirichlet_BCs_Part1






!+604+###########################################################################!
!                                                                                !
!                  CFA_Dirichlet_BCs                                                  !
!                                                                                !
!################################################################################!
SUBROUTINE CFA_3D_Dirichlet_BCs_Part2( )



INTEGER                                             ::  i, l, m, ui, d

INTEGER                                             ::  Value_Location



!*!
!*!     Step 3
!*!
DO ui = 1,NUM_CFA_VARS



    !*!
    !*!     Innner BCs
    !*!

    IF ( INNER_CFA_BC_TYPE(ui) == 'D' .AND. myID_PETSc == 0 ) THEN

        DO l = 0,L_LIMIT
            DO m = -M_VALUES(l),M_VALUES(l)

                Value_Location =  MAtrix_Location( ui, l, m, 0, 0 )
                RHS_VECTOR(Value_Location ) = 0.0_idp

             END DO
         END DO
                !*!
                !*!     Modify the Stiffness Matrix !
                !*!


         Value_Location =  MAtrix_Location( ui, 0, 0, 0, 0 )
         DO i = 0,ELEM_PROB_DIM-1

              BLOCK_ELEM_STF_MATVEC( i*ELEM_PROB_DIM+Value_Location, 0)=0.0_idp

         END DO
         BLOCK_ELEM_STF_MATVEC(Value_Location*ELEM_PROB_DIM:(Value_Location+1)*ELEM_PROB_DIM-1, 0) = 0.0_idp
         BLOCK_ELEM_STF_MATVEC(Value_Location*ELEM_PROB_DIM+ Value_Location, 0) = 1.0_idp


    END IF






    !*!
    !*!     Outer BCs
    !*!

    IF ( OUTER_CFA_BC_TYPE(ui) == 'D' .AND. myID_PETSC == NUM_SUBSHELLS-1 ) THEN

        DO l = 0,L_LIMIT
            DO m = -M_VALUES(l),M_VALUES(l)

                Value_Location =  Matrix_Location( ui, l, m, NUM_R_ELEMS_PER_BLOCK-1, DEGREE )
                Block_RHS_VECTOR(Value_Location ) = 0.0_idp

             END DO
         END DO

                !*!
                !*!     Modify the Stiffness Matrix !
                !*!
          Value_Location =  Matrix_Location( ui, 0, 0, 0, DEGREE )

          DO i = 0,ELEM_PROB_DIM-1

              ! Clear the Column !
              BLOCK_ELEM_STF_MATVEC( i*ELEM_PROB_DIM + Value_Location, NUM_R_ELEMS_PER_BLOCK-1) = 0.0_idp

          END DO
          ! Clear the Row !
          BLOCK_ELEM_STF_MATVEC( Value_Location*ELEM_PROB_DIM:(Value_Location+1)*ELEM_PROB_DIM-1,       &
                                 NUM_R_ELEMS_PER_BLOCK-1) = 0.0_idp

          BLOCK_ELEM_STF_MATVEC( Value_Location*ELEM_PROB_DIM+Value_Location, NUM_R_ELEMS_PER_BLOCK-1) = 1.0_idp

    END IF

END DO


END SUBROUTINE CFA_3D_Dirichlet_BCs_Part2






!+605+###########################################################################!
!                                                                                !
!                  CFA_Neumann_BCs                                               !
!                                                                                !
!################################################################################!
SUBROUTINE CFA_3D_Neumann_BCs( )


! Nothing needs to be done. Score!


END SUBROUTINE CFA_3D_Neumann_BCs










!+701+###########################################################################!
!                                                                                !
!                  Calc_3D_Values_At_Location          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_3D_Values_At_Location( r, theta, phi, Return_Psi, Return_AlphaPsi,  &
                                        Return_Beta1, Return_Beta2, Return_Beta3    )


REAL(KIND = idp), INTENT(IN)                                ::  r, theta, phi
REAL(KIND = idp), INTENT(INOUT)                             ::  Return_Psi,         &
                                                                Return_AlphaPsi,    &
                                                                Return_Beta1,       &
                                                                Return_Beta2,       &
                                                                Return_Beta3



COMPLEX(KIND = idp), DIMENSION(1:5)                         ::  Tmp_U_Value


REAL(KIND = idp)                                                ::  r_tmp
REAL(KIND = idp), DIMENSION(0:DEGREE)                           ::  LagP

INTEGER                                                         ::  re, l, m, d


INTEGER                                                         :: Current_Location



COMPLEX(KIND = idp)                                             ::  TMP_VALUE_A
REAL(KIND = idp), DIMENSION(0:DEGREE)                           ::  xlocP, weightP



Tmp_U_Value = 0.0_idp

IF ( r == rlocs(0) ) THEN

    DO l = 0,L_Limit
        DO m = -M_VALUES(l),M_VALUES(l)

            Current_Location =  Matrix_Location( 1, l, m, 0, 0 )

            Tmp_U_Value(1) = Tmp_U_Value(1)                                 &
                            + Coefficient_Vector( Current_Location + 0 )    &
                            * Spherical_Harmonic(l,m,theta,phi)

            Tmp_U_Value(2) = Tmp_U_Value(2)                                 &
                            + Coefficient_Vector( Current_Location + 1 )    &
                            * Spherical_Harmonic(l,m,theta,phi)

            Tmp_U_Value(3) = Tmp_U_Value(3)                                 &
                            + Coefficient_Vector( Current_Location + 2 )    &
                            * Spherical_Harmonic(l,m,theta,phi)

            Tmp_U_Value(4) = Tmp_U_Value(4)                                 &
                            + Coefficient_Vector( Current_Location + 3 )    &
                            * Spherical_Harmonic(l,m,theta,phi)

            Tmp_U_Value(5) = Tmp_U_Value(5)                                 &
                            + Coefficient_Vector( Current_Location + 4 )    &
                            * Spherical_Harmonic(l,m,theta,phi)

        END DO
    END DO

ELSE

    CALL Initialize_LGL_Quadrature(DEGREE,xlocP,weightP)

    DO re = 0,NUM_R_ELEMENTS-1

        IF ( r > rlocs(re) .AND. r <= rlocs(re+1) ) THEN

            r_tmp = Map_To_X_Space(rlocs(re),rlocs(re+1),r)
            LagP = Lagrange_Poly(r_tmp,DEGREE,xlocP)


            DO l = 0,L_Limit
                DO m = -M_VALUES(l),M_VALUES(l)
                    DO d = 0,DEGREE

                        TMP_VALUE_A = Spherical_Harmonic(l,m,theta,phi) * LagP(d)

                        Current_Location = Matrix_Location( 1, l, m, re, d )

                        Tmp_U_Value(1) = Tmp_U_Value(1) + Coefficient_Vector( Current_Location ) * TMP_VALUE_A
                        Tmp_U_Value(2) = Tmp_U_Value(2) + Coefficient_Vector( Current_Location + 1 ) * TMP_VALUE_A
                        Tmp_U_Value(3) = Tmp_U_Value(3) + Coefficient_Vector( Current_Location + 2 ) * TMP_VALUE_A
                        Tmp_U_Value(4) = Tmp_U_Value(4) + Coefficient_Vector( Current_Location + 3 ) * TMP_VALUE_A
                        Tmp_U_Value(5) = Tmp_U_Value(5) + Coefficient_Vector( Current_Location + 4 ) * TMP_VALUE_A


                    END DO  !   d Loop
                END DO  !   m Loop
            END DO  !   l Loop
            
            EXIT
        END IF

    END DO

END IF


Return_Psi = REAL(Tmp_U_Value(1), KIND = idp)
Return_AlphaPsi = REAL(Tmp_U_Value(2), KIND = idp)
Return_Beta1 = REAL(Tmp_U_Value(3), KIND = idp)
Return_Beta2 = REAL(Tmp_U_Value(4), KIND = idp)
Return_Beta3 = REAL(Tmp_U_Value(5), KIND = idp)



END SUBROUTINE Calc_3D_Values_At_Location



















END MODULE CFA_3D_Master_Build_Module
