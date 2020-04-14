   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_Jacobian_Functions                                                  !##!
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


USE Poseidon_Constants_Module, &
                                ONLY : idp, pi

USE Poseidon_Parameters, &
                                ONLY :  DOMAIN_DIM,                 &
                                        DEGREE,                     &
                                        L_LIMIT,                    &
                                        NUM_CFA_VARS

USE Poseidon_Variables_Module, &
                                ONLY :  ELEM_PROB_DIM


IMPLICIT NONE



CONTAINS



!+101+###########################################################################!
!                                                                                !
!                  CALC_JACOBIAN_A                                               !
!                                                                                !
!################################################################################!
SUBROUTINE CALC_JACOBIAN_A( Jacobian_Buffer, dp, lpmp_loc, F            &
                            SUBJACOBIAN_EQ1_TERM,                       &
                            SUBJACOBIAN_EQ2_TERM,                       &
                            SUBJACOBIAN_EQ3_TERM,                       &
                            SUBJACOBIAN_EQ4_TERM,                       &
                            SUBJACOBIAN_EQ5_TERM,                       )

COMPLEX, DIMENSION(0:ELEM_PROB_DIM-1), INTENT(INOUT)                                     :: Jacobian_Buffer
INTEGER, INTENT(IN)                                                                      :: dp, lpmp_loc, F

REAL(KIND = idp), DIMENSION(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:14), INTENT(IN) :: SUBJACOBIAN_EQ1_TERM
REAL(KIND = idp), DIMENSION(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:14), INTENT(IN) :: SUBJACOBIAN_EQ2_TERM
REAL(KIND = idp), DIMENSION(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:20), INTENT(IN) :: SUBJACOBIAN_EQ3_TERM
REAL(KIND = idp), DIMENSION(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:20), INTENT(IN) :: SUBJACOBIAN_EQ4_TERM
REAL(KIND = idp), DIMENSION(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:20), INTENT(IN) :: SUBJACOBIAN_EQ5_TERM


! d, u, and lm_loc pick the column
DO d = 0,DEGREE
    DO u = 1,NUM_CFA_VARS
        DO lm_loc = 0,LM_LENGTH-1

            

            ! column location
            j_loc = d*ULM_LENGTH + (u-1)*LM_LENGTH + lm_loc



            Jacobian_Buffer(j_loc) = Calc_Jacobian_Term( F, u, lm_loc, lpmp_loc, d, dp,             &
                                                              SubJacobian_EQ1_Term,           &
                                                              SubJacobian_EQ2_Term,           &
                                                              SubJacobian_EQ3_Term,           &
                                                              SubJacobian_EQ4_Term,           &
                                                              SubJacobian_EQ5_Term           )

        END DO ! lm_loc Loop
    END DO ! u Loop
END DO ! d Loop

END SUBROUTINE CALC_JACOBIAN_A







END MODULE Poseidon_Jacobian_Functions
