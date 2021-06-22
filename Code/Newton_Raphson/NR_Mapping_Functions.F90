   !#################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE NR_Mapping_Functions                                                         !##!
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
!                  NR_Beta_Array_Map                                                  !
!                                                                                !
!################################################################################!
PURE FUNCTION NR_Array_Map(re,d,ui,lm_loc, m_opt)

INTEGER                                     ::  NR_Array_Map

INTEGER, INTENT(IN)                         ::  re
INTEGER, INTENT(IN)                         ::  d
INTEGER, INTENT(IN)                         ::  ui
INTEGER, INTENT(IN)                         ::  lm_loc
INTEGER, INTENT(IN), OPTIONAL               ::  m_opt


IF ( PRESENT(m_opt) ) THEN

    NR_Array_Map = (re*Degree + d) * 5 * LM_Length   &
                 + (ui - 1) * LM_Length              &
                 + lm_loc*(lm_loc+1) + m_opt


ELSE
    NR_Array_Map = (re*Degree + d) * 5 * LM_Length   &
                 + (ui - 1) * LM_Length              &
                 + lm_loc

END IF


END FUNCTION NR_Array_Map





!+201+###########################################################################!
!                                                                                !
!                  NR_Vector_Map                                                  !
!                                                                                !
!################################################################################!
PURE FUNCTION NR_FEM_Node_Map( re, d )

INTEGER                                     :: NR_FEM_Node_Map

INTEGER, INTENT(IN)                         :: re, d


NR_FEM_Node_Map = re*DEGREE + d + 1


END FUNCTION NR_FEM_Node_Map






!+201+###########################################################################!
!                                                                                !
!                  NR_Vector_Map                                                  !
!                                                                                !
!################################################################################!
PURE FUNCTION NR_LM_Map( l, m )

INTEGER                                     :: NR_LM_Map

INTEGER, INTENT(IN)                         :: l, m


NR_LM_Map = l*(l+1) + m + 1


END FUNCTION NR_LM_Map




!+201+###########################################################################!
!                                                                                !
!                  NR_Vector_Map                                                  !
!                                                                                !
!################################################################################!
PURE FUNCTION NR_tpd_Map( td, pd )

INTEGER                                     :: NR_tpd_Map

INTEGER, INTENT(IN)                         :: td, pd


NR_tpd_Map = (td-1)*Num_P_Quad_Points + pd


END FUNCTION NR_tpd_Map




END MODULE NR_Mapping_Functions


