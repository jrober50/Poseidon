   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Source_Input_Module                                                          !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!


USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Poseidon_Internal_Communication_Module, &
            ONLY :  Poseidon_CFA_Block_Share

USE Poseidon_Parameters, &
            ONLY :  DOMAIN_DIM

USE Variables_Mesh, &
            ONLY :  NUM_R_ELEMENTS,     &
                    NUM_T_ELEMENTS,     &
                    NUM_P_ELEMENTS,     &
                    drlocs,             &
                    dtlocs,             &
                    dplocs

USE Variables_Source, &
            ONLY :  Block_Source_E,     &
                    Block_Source_S,     &
                    Block_Source_Si

USE Poseidon_IO_Module, &
            ONLY :  OUTPUT_POSEIDON_SOURCES_1D

use mpi








IMPLICIT NONE

CONTAINS


!+101+##########################################################################!
!                                                                               !
!                           Poseidon_Input_Sources                              !
!                                                                               !
!###############################################################################!
SUBROUTINE Poseidon_Input_Sources(  ProcID, ProcID_Theta, ProcID_Phi,                   &
                                    Local_E, Local_S, Local_Si,                         &
                                    Local_RE_Dim, Local_TE_Dim, Local_PE_Dim,           &
                                    Local_RQ_Dim, Local_TQ_Dim, Local_PQ_Dim,           &
                                    Input_R_Quad, Input_T_Quad, Input_P_Quad,           &
                                    Left_Limit, Right_Limit                             )

INTEGER, INTENT(IN)                                                         ::  ProcID,        &
                                                                                ProcID_Theta,  &
                                                                                ProcID_Phi


REAL(KIND = idp), INTENT(IN), DIMENSION(    1:Local_RQ_Dim*Local_TQ_Dim*Local_PQ_Dim,       &
                                            0:Local_RE_Dim-1,                               &
                                            0:Local_TE_Dim-1,                               &
                                            0:Local_PE_Dim-1  )             ::  Local_E,    &
                                                                                Local_S

REAL(KIND = idp), INTENT(IN), DIMENSION(    1:Local_RQ_Dim*Local_TQ_Dim*Local_PQ_Dim,       &
                                            0:Local_RE_Dim-1,                               &
                                            0:Local_TE_Dim-1,                               &
                                            0:Local_PE_Dim-1,                               &
                                            1:DOMAIN_DIM                )   :: Local_Si





INTEGER, INTENT(IN)                                                     ::  Local_RE_Dim,   &
                                                                            Local_TE_Dim,   &
                                                                            Local_PE_Dim,   &
                                                                            Local_RQ_Dim,   &
                                                                            Local_TQ_Dim,   &
                                                                            Local_PQ_Dim


REAL(KIND = idp), DIMENSION(1:Local_RQ_Dim), INTENT(IN)                 ::  Input_R_Quad
REAL(KIND = idp), DIMENSION(1:Local_TQ_Dim), INTENT(IN)                 ::  Input_T_Quad
REAL(KIND = idp), DIMENSION(1:Local_PQ_Dim), INTENT(IN)                 ::  Input_P_Quad

REAL(KIND = idp), INTENT(IN)                                            ::  Left_Limit,     &
                                                                            Right_Limit


REAL(KIND = idp), DIMENSION(    1:Local_RQ_Dim*Local_TQ_Dim*Local_PQ_Dim,       &
                                            0:Local_RE_Dim-1,                               &
                                            0:Local_TE_Dim-1,                               &
                                            0:Local_PE_Dim-1,                               &
                                            1:DOMAIN_DIM                )   ::      Tmp_Si

INTEGER                                                                 :: re, rd, pe, te





CALL Poseidon_CFA_Block_Share(  ProcID, ProcID_Theta, ProcID_Phi,                   &
                                Local_E, Local_S, Local_Si,                         &
                                Local_RE_Dim, Local_TE_Dim, Local_PE_Dim,           &
                                Local_RQ_Dim, Local_TQ_Dim, Local_PQ_Dim,           &
                                Input_R_Quad, Input_T_Quad, Input_P_Quad,           &
                                Left_Limit, Right_Limit,                            &
                                NUM_R_ELEMENTS, NUM_T_ELEMENTS, NUM_P_ELEMENTS,     &
                                drlocs, dtlocs, dplocs,                             &
                                Block_Source_E, Block_Source_S, Block_Source_Si     )




END SUBROUTINE Poseidon_Input_Sources


!+101+##########################################################################!
!                                                                               !
!                           Poseidon_Input_Sources                              !
!                                                                               !
!###############################################################################!
SUBROUTINE Poseidon_Input_MCL_Sources(  Local_E, Local_S, Local_Si,                         &
                                        Local_RE_Dim, Local_TE_Dim, Local_PE_Dim,           &
                                        Local_RQ_Dim, Local_TQ_Dim, Local_PQ_Dim,           &
                                        Input_R_Quad, Input_T_Quad, Input_P_Quad,           &
                                        Left_Limit, Right_Limit                             )



REAL(KIND = idp), INTENT(IN), DIMENSION(    1:Local_RQ_Dim*Local_TQ_Dim*Local_PQ_Dim,       &
                                            0:Local_RE_Dim-1,                               &
                                            0:Local_TE_Dim-1,                               &
                                            0:Local_PE_Dim-1  )             ::  Local_E,    &
                                                                                Local_S

REAL(KIND = idp), INTENT(IN), DIMENSION(    1:Local_RQ_Dim*Local_TQ_Dim*Local_PQ_Dim,       &
                                            0:Local_RE_Dim-1,                               &
                                            0:Local_TE_Dim-1,                               &
                                            0:Local_PE_Dim-1,                               &
                                            1:DOMAIN_DIM                )   :: Local_Si





INTEGER, INTENT(IN)                                                     ::  Local_RE_Dim,   &
                                                                            Local_TE_Dim,   &
                                                                            Local_PE_Dim,   &
                                                                            Local_RQ_Dim,   &
                                                                            Local_TQ_Dim,   &
                                                                            Local_PQ_Dim


REAL(KIND = idp), DIMENSION(1:Local_RQ_Dim), INTENT(IN)                 ::  Input_R_Quad
REAL(KIND = idp), DIMENSION(1:Local_TQ_Dim), INTENT(IN)                 ::  Input_T_Quad
REAL(KIND = idp), DIMENSION(1:Local_PQ_Dim), INTENT(IN)                 ::  Input_P_Quad

REAL(KIND = idp), INTENT(IN)                                            ::  Left_Limit,     &
                                                                            Right_Limit


REAL(KIND = idp), DIMENSION(    1:Local_RQ_Dim*Local_TQ_Dim*Local_PQ_Dim,       &
                                            0:Local_RE_Dim-1,                               &
                                            0:Local_TE_Dim-1,                               &
                                            0:Local_PE_Dim-1,                               &
                                            1:DOMAIN_DIM                )   ::      Tmp_Si

INTEGER                                                                 :: re, rd, pe, te, td, pd

INTEGER                                                                 :: Here





DO pe = 0,Local_PE_Dim-1
DO te = 0,Local_TE_Dim-1
DO re = 0,Local_RE_Dim-1

DO rd = 1,Local_RQ_Dim
DO td = 1,Local_TQ_Dim
DO pd = 1,Local_PQ_Dim



    Here = (rd-1) * Local_PQ_Dim * Local_TQ_Dim       &
         + (td-1) * Local_PQ_Dim                     &
         + pd
    Block_Source_E(rd,td,pd,re,te,pe) = Local_E(Here,re,te,pe)

END DO
END DO
END DO

END DO
END DO
END DO




END SUBROUTINE Poseidon_Input_MCL_Sources



!+101+##########################################################################!
!                                                                               !
!                           Poseidon_Input_Sources                              !
!                                                                               !
!###############################################################################!
SUBROUTINE Poseidon_Input_Sources_Serial(   Local_E, Local_S, Local_Si,                         &
                                            Local_RE_Dim, Local_TE_Dim, Local_PE_Dim,           &
                                            Local_RQ_Dim, Local_TQ_Dim, Local_PQ_Dim,           &
                                            Input_R_Quad, Input_T_Quad, Input_P_Quad,           &
                                            Left_Limit, Right_Limit                             )



REAL(KIND = idp), INTENT(IN), DIMENSION(    1:Local_RQ_Dim*Local_TQ_Dim*Local_PQ_Dim,       &
                                            0:Local_RE_Dim-1,                               &
                                            0:Local_TE_Dim-1,                               &
                                            0:Local_PE_Dim-1  )             ::  Local_E,    &
                                                                                Local_S

REAL(KIND = idp), INTENT(IN), DIMENSION(    1:Local_RQ_Dim*Local_TQ_Dim*Local_PQ_Dim,       &
                                            0:Local_RE_Dim-1,                               &
                                            0:Local_TE_Dim-1,                               &
                                            0:Local_PE_Dim-1,                               &
                                            1:DOMAIN_DIM                )   :: Local_Si





INTEGER, INTENT(IN)                                                     ::  Local_RE_Dim,   &
                                                                            Local_TE_Dim,   &
                                                                            Local_PE_Dim,   &
                                                                            Local_RQ_Dim,   &
                                                                            Local_TQ_Dim,   &
                                                                            Local_PQ_Dim


REAL(KIND = idp), DIMENSION(1:Local_RQ_Dim), INTENT(IN)                 ::  Input_R_Quad
REAL(KIND = idp), DIMENSION(1:Local_TQ_Dim), INTENT(IN)                 ::  Input_T_Quad
REAL(KIND = idp), DIMENSION(1:Local_PQ_Dim), INTENT(IN)                 ::  Input_P_Quad

REAL(KIND = idp), INTENT(IN)                                            ::  Left_Limit,     &
                                                                            Right_Limit


REAL(KIND = idp), DIMENSION(    1:Local_RQ_Dim*Local_TQ_Dim*Local_PQ_Dim,       &
                                            0:Local_RE_Dim-1,                               &
                                            0:Local_TE_Dim-1,                               &
                                            0:Local_PE_Dim-1,                               &
                                            1:DOMAIN_DIM                )   ::      Tmp_Si

INTEGER                                                                 :: re, rd, pe, te, td, pd

INTEGER                                                                 :: Here





DO pe = 0,Local_PE_Dim-1
DO te = 0,Local_TE_Dim-1
DO re = 0,Local_RE_Dim-1

DO rd = 1,Local_RQ_Dim
DO td = 1,Local_TQ_Dim
DO pd = 1,Local_PQ_Dim



    Here = (rd-1) * Local_PQ_Dim * Local_TQ_Dim       &
         + (td-1) * Local_PQ_Dim                     &
         + pd
    Block_Source_E(rd,td,pd,re,te,pe) = Local_E(Here,re,te,pe)
    Block_Source_S(rd,td,pd,re,te,pe) = Local_S(Here,re,te,pe)
    Block_Source_Si(rd,td,pd,re,te,pe,:) = Local_Si(Here,re,te,pe,:)

END DO
END DO
END DO

END DO
END DO
END DO




END SUBROUTINE Poseidon_Input_Sources_Serial








END MODULE Source_Input_Module
