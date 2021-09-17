   !#################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE FP_Functions_Mapping                                                         !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !#################################################################################!

!*D*================================!
!                                   !
!           Dependencies            !
!                                   !
!===================================!
USE Poseidon_Kinds_Module, &
            ONLY :  idp


USE Poseidon_Parameters, &
            ONLY :  DEGREE

USE Variables_Derived, &
            ONLY :  LM_Length

USE Variables_Quadrature, &
            ONLY :  Num_R_Quad_Points,          &
                    Num_T_Quad_Points,          &
                    Num_P_Quad_Points

IMPLICIT NONE


CONTAINS



!+201+###########################################################################!
!                                                                                !
!                  FP_Array_Map                                                  !
!                                                                                !
!################################################################################!
PURE FUNCTION FP_Array_Map(re,d,ui,lm_loc, m_opt)

INTEGER                                     ::  FP_Array_Map

INTEGER, INTENT(IN)                         ::  re
INTEGER, INTENT(IN)                         ::  d
INTEGER, INTENT(IN)                         ::  ui
INTEGER, INTENT(IN)                         ::  lm_loc
INTEGER, INTENT(IN), OPTIONAL               ::  m_opt


IF ( PRESENT(m_opt) ) THEN

    FP_Array_Map = (re*Degree + d) * 5 * LM_Length   &
                 + (ui - 1) * LM_Length              &
                 + lm_loc*(lm_loc+1) + m_opt + 1


ELSE
    FP_Array_Map = (re*Degree + d) * 5 * LM_Length   &
                 + (ui - 1) * LM_Length              &
                 + lm_loc

END IF


END FUNCTION FP_Array_Map


!+201+###########################################################################!
!                                                                                !
!                  FP_Beta_Array_Map                                             !
!                                                                                !
!################################################################################!
PURE FUNCTION FP_Beta_Array_Map(re,d,ui,lm_loc, m_opt)

INTEGER                                     ::  FP_Beta_Array_Map

INTEGER, INTENT(IN)                         ::  re
INTEGER, INTENT(IN)                         ::  d
INTEGER, INTENT(IN)                         ::  ui
INTEGER, INTENT(IN)                         ::  lm_loc
INTEGER, INTENT(IN), OPTIONAL               ::  m_opt


IF ( PRESENT(m_opt) ) THEN

    FP_Beta_Array_Map = (re*Degree + d) * 3 * LM_Length   &
                      + (ui - 1) * LM_Length              &
                      + lm_loc*(lm_loc+1) + m_opt + 1


ELSE
    FP_Beta_Array_Map = (re*Degree + d) * 3 * LM_Length   &
                      + (ui - 1) * LM_Length              &
                      + lm_loc

END IF


END FUNCTION FP_Beta_Array_Map



!+201+###########################################################################!
!                                                                                !
!                  FP_Beta_Array_Map                                                  !
!                                                                                !
!################################################################################!
PURE FUNCTION FP_X_Array_Map(re,d,ui,lm_loc, m_opt)

INTEGER                                     ::  FP_X_Array_Map

INTEGER, INTENT(IN)                         ::  re
INTEGER, INTENT(IN)                         ::  d
INTEGER, INTENT(IN)                         ::  ui
INTEGER, INTENT(IN)                         ::  lm_loc
INTEGER, INTENT(IN), OPTIONAL               ::  m_opt


IF ( PRESENT(m_opt) ) THEN

    FP_X_Array_Map = (re*Degree + d) * 3 * LM_Length   &
                      + (ui - 6) * LM_Length              &
                      + lm_loc*(lm_loc+1) + m_opt + 1


ELSE
    FP_X_Array_Map = (re*Degree + d) * 3 * LM_Length   &
                      + (ui - 6) * LM_Length              &
                      + lm_loc

END IF


END FUNCTION FP_X_Array_Map







!+201+###########################################################################!
!                                                                                !
!                  FP_Beta_Array_Map                                                  !
!                                                                                !
!################################################################################!
FUNCTION FP_Array_Map_TypeB(ui,iVB,re,d,lm_loc, m_opt)

INTEGER                                     ::  FP_Array_Map_TypeB

INTEGER, INTENT(IN)                         ::  ui
INTEGER, INTENT(IN)                         ::  iVB
INTEGER, INTENT(IN)                         ::  re
INTEGER, INTENT(IN)                         ::  d
INTEGER, INTENT(IN)                         ::  lm_loc
INTEGER, INTENT(IN), OPTIONAL               ::  m_opt

INTEGER                                     ::  Offset

Offset = 3*iVB


IF ( PRESENT(m_opt) ) THEN

    FP_Array_Map_TypeB = (re*Degree + d) * 3 * LM_Length   &
                       + (ui - Offset) * LM_Length         &
                       + lm_loc*(lm_loc+1) + m_opt + 1


ELSE

    FP_Array_Map_TypeB = (re*Degree + d) * 3 * LM_Length   &
                       + (ui - Offset) * LM_Length         &
                       + lm_loc


END IF


END FUNCTION FP_Array_Map_TypeB














!+201+###########################################################################!
!                                                                                !
!                  FP_Vector_Map                                                  !
!                                                                                !
!################################################################################!
PURE FUNCTION FP_FEM_Node_Map( re, d )

INTEGER                                     :: FP_FEM_Node_Map

INTEGER, INTENT(IN)                         :: re, d


FP_FEM_Node_Map = re*DEGREE + d + 1


END FUNCTION FP_FEM_Node_Map






!+201+###########################################################################!
!                                                                                !
!                  FP_Vector_Map                                                  !
!                                                                                !
!################################################################################!
PURE FUNCTION FP_LM_Map( l, m )

INTEGER                                     :: FP_LM_Map

INTEGER, INTENT(IN)                         :: l, m


FP_LM_Map = l*(l+1) + m + 1


END FUNCTION FP_LM_Map




!+201+###########################################################################!
!                                                                                !
!                  FP_Vector_Map                                                  !
!                                                                                !
!################################################################################!
PURE FUNCTION FP_tpd_Map( td, pd )

INTEGER                                     :: FP_tpd_Map

INTEGER, INTENT(IN)                         :: td, pd


FP_tpd_Map = (td-1)*Num_P_Quad_Points + pd


END FUNCTION FP_tpd_Map




END MODULE FP_Functions_Mapping


