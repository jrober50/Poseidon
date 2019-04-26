   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_Newton_Module                                                       !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!================================================================================!##!
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


USE Poseidon_Constants_Module, &
                                ONLY :  idp

USE Poseidon_Parameters, &
                                ONLY :  DOMAIN_DIM,                 &
                                        DEGREE,                     &
                                        L_LIMIT,                    &
                                        nPROCS_POSEIDON,            &
                                        NUM_R_ELEMS_PER_BLOCK,      &
                                        NUM_T_ELEMS_PER_BLOCK,      &
                                        NUM_P_ELEMS_PER_BLOCK,      &
                                        NUM_R_ELEMS_PER_SUBSHELL,   &
                                        NUM_BLOCK_THETA_ROWS,       &
                                        NUM_BLOCK_PHI_COLUMNS





!+101+###########################################################################!
!                                                                                !
!           Newtonian_3D_Master_Build                                                  !
!                                                                                !
!################################################################################!
SUBROUTINE Newtonian_3D_Master_Build()



REAL(KIND = idp)        :: timea, timeb, timec

INTEGER                 :: re

timea = 0.0_idp
timeb = 0.0_idp
timec = 0.0_idp


!Block_STF_MATVEC = 0.0_idp
Block_RHS_Vector = 0.0_idp
!BLOCK_ELEM_STF_MATVEC = 0.0_idp




!timeb = MPI_Wtime()
CALL CREATE_3D_NEWTONAIN_RHS()
!timea = MPI_Wtime()
!CALL Clock_In(timea-timeb, 9)


CALL CREATE_3D_NEWTONIAN_LAPLACIAN_MATRIX()




END SUBROUTINE Newtonian_3D_Master_Build





!+201+##########################################################################!
!                                                                               !
!                  CREATE_3D_NONLAPLACIAN_SOE                                   !
!                                                                               !
!###############################################################################!
SUBROUTINE CREATE_3D_NEWTONAIN_RHS()




INTEGER                                                         ::  Local_re,       &
                                                                    Local_te,       &
                                                                    Local_pe,       &
                                                                    Global_re,      &
                                                                    Global_te,      &
                                                                    Global_pe,      &
                                                                    d, l, m,        &
                                                                    dp, lp, mp,     &
                                                                    rd, td, pd




REAL(KIND = idp), DIMENSION(1:NUM_R_QUAD_POINTS)                ::  CUR_R_LOCS
REAL(KIND = idp), DIMENSION(1:NUM_T_QUAD_POINTS)                ::  CUR_T_LOCS
REAL(KIND = idp), DIMENSION(1:NUM_P_QUAD_POINTS)                ::  CUR_P_LOCS






COMPLEX(KIND = idp)                                                     ::  REUSED_VALUE



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



    !
    !   Move through phi
    !
    DO Local_pe = 0, NUM_P_ELEMS_PER_BLOCK-1



        Global_pe = Block_P_Begin + Local_pe


        deltap_overtwo = (plocs(Global_pe + 1) - plocs(Global_pe))/2.0_idp
        CUR_P_LOCS(:) = deltap_overtwo * (INT_P_LOCATIONS(:)+1.0_idp) + plocs(Global_pe)



        DO Local_te = 0, NUM_T_ELEMS_PER_BLOCK-1


            Global_te = Block_T_Begin + Local_te

            deltat_overtwo = (tlocs(Global_te + 1) - tlocs(Global_te))/2.0_idp
            CUR_T_LOCS(:) = deltat_overtwo * (INT_T_LOCATIONS(:)+1.0_idp) + tlocs(Global_te)




            DO pd = 1,NUM_P_QUAD_POINTS
                DO td = 1,NUM_T_QUAD_POINTS
                    DO rd = 1,NUM_R_QUAD_POINTS

                        Int_Factor(rd,td,pd) = SIN_VAL(td)*R_SQUARE(rd)                 &
                                                * DELTAR_OVERTWO * INT_R_WEIGHTS(rd)    &
                                                * DELTAT_OVERTWO * INT_T_WEIGHTS(td)    &
                                                * DELTAP_OVERTWO * INT_P_WEIGHTS(pd)

                    END DO  !   rd Loop
                END DO  !   td Loop
            END DO  !   pd Loop



            DO d = 0,DEGREE
                DO lm_loc = 0,LM_LENGTH-1

                    Current_i_Location = CFA_ALL_Matrix_Map(1, lm_loc, re, d)

                    RHS_TMP = 0.0_idp
                    DO pd = 1,NUM_P_QUAD_POINTS
                        DO td = 1,NUM_T_QUAD_POINTS
                            DO rd = 1,NUM_R_QUAD_POINTS


                                RHS_TMP =  RHS_TMP                                          &
                                        + Block_Source_E( rd, td, pd, re, te, pe )          &
                                            * Ylm_CC_Values(lm_loc, td, pd, te, pe)         &
                                            * Lagrange_Poly_Table(d, rd, 0)                 &
                                            * Int_Factor(rd, td, pd)



                            END DO  ! rd Loop
                        END DO  ! td Loop
                    END DO  ! pd Loop


                    Block_RHS_Vector(Current_i_Location)                             &
                        = Block_RHS_Vector(Current_i_Location)                       &
                        + RHS_TMP

                END DO  ! l Loop
            END DO  ! d Loop

        END DO  ! re Loop
    END DO  ! te Loop
END DO  ! pe Loop







END SUBROUTINE CREATE_3D_NEWTONAIN_RHS







!+204+###########################################################################!
!                                                                                !
!                  CREATE_3D_NEWTON_RHS_VECTOR                                   !
!                                                                                !
!################################################################################!
SUBROUTINE CREATE_3D_NEWTON_RHS_VECTOR(    re, te, pe,                          &
                                           R_SQUARE, DELTAR_OVERTWO,            &
                                           Int_Factor                           )



INTEGER, INTENT(IN)                                                         ::  re, te, pe



REAL(KIND = idp), INTENT(IN), DIMENSION( 1:NUM_R_QUAD_POINTS,        &
                                         1:NUM_T_QUAD_POINTS,        &
                                         1:NUM_P_QUAD_POINTS         )      ::  Int_Factor


REAL(KIND = idp), INTENT(IN), DIMENSION(1:NUM_R_QUAD_POINTS)            ::  R_SQUARE
REAL(KIND = idp), INTENT(IN)                                            ::  DELTAR_OVERTWO

INTEGER                                                                 ::  pd, td, rd,     &
                                                                            l, m, d,        &
                                                                            lm_loc

INTEGER                                                                 ::  Current_i_Location

COMPLEX(KIND = idp)                                                     ::  RHS_TMP





!$OMP PARALLEL DEFAULT(none)                                            &
!$OMP PRIVATE( pd, td, rd, d, l, m,                                     &
!$OMP           lm_loc,                                                 &
!$OMP           Current_i_Location, RHS_TMP                         )   &
!$OMP SHARED( re, te, pe,                                               &
!$OMP           DEGREE,                                                 &
!$OMP           R_SQUARE, Deltar_Overtwo, Int_R_Weights,                &
!$OMP           NUM_P_QUAD_POINTS, NUM_T_QUAD_POINTS, NUM_R_QUAD_POINTS,&
!$OMP           Int_Factor,                                             &
!$OMP           ULM_LENGTH, LM_LENGTH,                                  &
!$OMP           Ylm_CC_Values, Lagrange_Poly_Table,                     &
!$OMP           LM_Location,                                            &
!$OMP           Block_Source_E,                                         &
!$OMP           Block_RHS_Vector, RHS_Terms                             )

!$OMP DO SCHEDULE(dynamic), COLLAPSE(2)
DO d = 0,DEGREE
    DO lm_loc = 0,LM_LENGTH-1

        Current_i_Location = CFA_ALL_Matrix_Map(1, lm_loc, re, d)

        RHS_TMP = 0.0_idp
        DO pd = 1,NUM_P_QUAD_POINTS
            DO td = 1,NUM_T_QUAD_POINTS
                DO rd = 1,NUM_R_QUAD_POINTS


                    RHS_TMP =  RHS_TMP                                          &
                            + Block_Source_E( rd, td, pd, re, te, pe )          &
                                * Ylm_CC_Values(lm_loc, td, pd, te, pe)         &
                                * Lagrange_Poly_Table(d, rd, 0)                 &
                                * Int_Factor(rd, td, pd)



                END DO  ! rd Loop
            END DO  ! td Loop
        END DO  ! pd Loop


        Block_RHS_Vector(Current_i_Location)                             &
            = Block_RHS_Vector(Current_i_Location)                       &
            + RHS_TMP

    END DO  ! l Loop
END DO  ! d Loop
!$OMP END DO

!$OMP END PARALLEL



END SUBROUTINE CREATE_3D_RHS_VECTOR















!+401+##############################################################################!
!                                                                                   !
!                  FINISH_3D_NONLAPLACIAN_MATRIX                             !
!                                                                                   !
!###################################################################################!
SUBROUTINE CREATE_LAPLACIAN_MATRIX()



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


                    Reusable_Values(dp) = SUM ( ( -R_SQUARE(:) * LPT_LPT(:,dp,d,1,1)*TWOOVER_DELTAR*TWOOVER_DELTAR      &
                                                + L_Lp1 * LPT_LPT(:,dp,d,0,0)        )    &
                                                * INT_R_WEIGHTS(:)/TWOOVER_DELTAR     )


                END DO


                DO m = -M_VALUES(l),M_VALUES(l)

                   jloc = d*ULM_LENGTH + (l*(l+1) + m)


                    DO dp = 0,DEGREE


                        iloc = dp*ULM_LENGTH + (l*(l+1) + m)


                        ! Psi Laplacian Terms
                        MATVEC_LOC = jloc*ELEM_PROB_DIM + iloc
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


END SUBROUTINE CREATE_LAPLACIAN_MATRIX









END MODULE Poseidon_Newton_Module
