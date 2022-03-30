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

USE Source_Intput_Poisson,  &
            ONLY :  Poseidon_Input_Sources_Poisson_3D


USE Timer_Routines_Module, &
            ONLY :  TimerStart,                     &
                    TimerStop

USE Timer_Variables_Module, &
            ONLY :  Timer_GR_SourceInput,         &
                    Timer_GR_SourceInput_PartA,        &
                    Timer_GR_SourceInput_PartB

USE Quadrature_Mapping_Functions, &
            ONLY :  Quad_Map


use mpi








IMPLICIT NONE

INTERFACE Poseidon_Input_Sources

    PROCEDURE   Poseidon_Input_Sources_CFA,         &
                Poseidon_Input_Sources_Poisson_3D


END INTERFACE Poseidon_Input_Sources

CONTAINS


!+101+##########################################################################!
!                                                                               !
!                           Poseidon_Input_Sources                              !
!                                                                               !
!###############################################################################!
SUBROUTINE Poseidon_Input_Sources_CFA(  ProcID, ProcID_Theta, ProcID_Phi,                   &
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


CALL TimerStart(Timer_GR_SourceInput)


CALL Poseidon_CFA_Block_Share(  ProcID, ProcID_Theta, ProcID_Phi,                   &
                                Local_E, Local_S, Local_Si,                         &
                                Local_RE_Dim, Local_TE_Dim, Local_PE_Dim,           &
                                Local_RQ_Dim, Local_TQ_Dim, Local_PQ_Dim,           &
                                Input_R_Quad, Input_T_Quad, Input_P_Quad,           &
                                Left_Limit, Right_Limit,                            &
                                NUM_R_ELEMENTS, NUM_T_ELEMENTS, NUM_P_ELEMENTS,     &
                                drlocs, dtlocs, dplocs,                             &
                                Block_Source_E, Block_Source_S, Block_Source_Si     )

CALL TimerStop(Timer_GR_SourceInput)


END SUBROUTINE Poseidon_Input_Sources_CFA










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



INTEGER                                                                 :: re, rd, pe, te, td, pd

INTEGER                                                                 :: Here, There



CALL TimerStart(Timer_GR_SourceInput)

DO pe = 0,Local_PE_Dim-1
DO te = 0,Local_TE_Dim-1
DO re = 0,Local_RE_Dim-1

DO rd = 1,Local_RQ_Dim
DO td = 1,Local_TQ_Dim
DO pd = 1,Local_PQ_Dim



    Here = (rd-1) * Local_PQ_Dim * Local_TQ_Dim       &
         + (td-1) * Local_PQ_Dim                     &
         + pd

    There = Quad_Map(rd,td,pd)
    Block_Source_E(There,re,te,pe) = Local_E(Here,re,te,pe)

END DO
END DO
END DO

END DO
END DO
END DO


CALL TimerStop(Timer_GR_SourceInput)



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



INTEGER                                                                 :: re, rd, pe, te, td, pd

INTEGER                                                                 :: Here, There


CALL TimerStart(Timer_GR_SourceInput)


DO pe = 0,Local_PE_Dim-1
DO te = 0,Local_TE_Dim-1
DO re = 0,Local_RE_Dim-1

DO rd = 1,Local_RQ_Dim
DO td = 1,Local_TQ_Dim
DO pd = 1,Local_PQ_Dim



    Here = (rd-1) * Local_PQ_Dim * Local_TQ_Dim       &
         + (td-1) * Local_PQ_Dim                     &
         + pd

    There = Quad_Map(rd,td,pd)
    Block_Source_E(There,re,te,pe) = Local_E(Here,re,te,pe)
    Block_Source_S(There,re,te,pe) = Local_S(Here,re,te,pe)
    Block_Source_Si(There,re,te,pe,:) = Local_Si(Here,re,te,pe,:)

END DO
END DO
END DO

END DO
END DO
END DO

CALL TimerStop(Timer_GR_SourceInput)


END SUBROUTINE Poseidon_Input_Sources_Serial

















!+101+##########################################################################!
!                                                                               !
!                           Poseidon_Input_Sources                              !
!                                                                               !
!###############################################################################!
SUBROUTINE Poseidon_XCFC_Input_Sources1(    Local_E, Local_Si,                          &
                                            Local_RE_Dim, Local_TE_Dim, Local_PE_Dim,   &
                                            Local_RQ_Dim, Local_TQ_Dim, Local_PQ_Dim,   &
                                            Input_R_Quad, Input_T_Quad, Input_P_Quad,   &
                                            Left_Limit, Right_Limit                     )



REAL(KIND = idp), INTENT(IN), DIMENSION(    1:Local_RQ_Dim*Local_TQ_Dim*Local_PQ_Dim,       &
                                            0:Local_RE_Dim-1,                               &
                                            0:Local_TE_Dim-1,                               &
                                            0:Local_PE_Dim-1  )             ::  Local_E

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

INTEGER                                                                 :: re, rd, pe, te, td, pd
INTEGER                                                                 :: Here, There


CALL TimerStart(Timer_GR_SourceInput)
CALL TimerStart(Timer_GR_SourceInput_PartA)


DO pe = 0,Local_PE_Dim-1
DO te = 0,Local_TE_Dim-1
DO re = 0,Local_RE_Dim-1

DO rd = 1,Local_RQ_Dim
DO td = 1,Local_TQ_Dim
DO pd = 1,Local_PQ_Dim



    Here = (rd-1) * Local_PQ_Dim * Local_TQ_Dim       &
         + (td-1) * Local_PQ_Dim                     &
         + pd

    There = Quad_Map(rd,td,pd)
    Block_Source_E(There,re,te,pe) = Local_E(Here,re,te,pe)
    Block_Source_Si(There,re,te,pe,:) = Local_Si(Here,re,te,pe,:)

END DO
END DO
END DO

END DO
END DO
END DO


CALL TimerStop(Timer_GR_SourceInput)
CALL TimerStop(Timer_GR_SourceInput_PartA)


END SUBROUTINE Poseidon_XCFC_Input_Sources1



!+101+##########################################################################!
!                                                                               !
!                           Poseidon_Input_Sources                              !
!                                                                               !
!###############################################################################!
SUBROUTINE Poseidon_XCFC_Input_Sources2(    Local_S,                                    &
                                            Local_RE_Dim, Local_TE_Dim, Local_PE_Dim,   &
                                            Local_RQ_Dim, Local_TQ_Dim, Local_PQ_Dim,   &
                                            Input_R_Quad, Input_T_Quad, Input_P_Quad,   &
                                            Left_Limit, Right_Limit                     )



REAL(KIND = idp), INTENT(IN), DIMENSION(    1:Local_RQ_Dim*Local_TQ_Dim*Local_PQ_Dim,       &
                                            0:Local_RE_Dim-1,                               &
                                            0:Local_TE_Dim-1,                               &
                                            0:Local_PE_Dim-1  )             ::  Local_S



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

INTEGER                                                                 :: re, rd, pe, te, td, pd
INTEGER                                                                 :: Here, There


CALL TimerStart(Timer_GR_SourceInput)
CALL TimerStart(Timer_GR_SourceInput_PartB)



DO pe = 0,Local_PE_Dim-1
DO te = 0,Local_TE_Dim-1
DO re = 0,Local_RE_Dim-1

DO rd = 1,Local_RQ_Dim
DO td = 1,Local_TQ_Dim
DO pd = 1,Local_PQ_Dim

    Here = (rd-1) * Local_PQ_Dim * Local_TQ_Dim       &
         + (td-1) * Local_PQ_Dim                     &
         + pd

    There = Quad_Map(rd,td,pd)
    Block_Source_S(There,re,te,pe) = Local_S(Here,re,te,pe)

END DO
END DO
END DO

END DO
END DO
END DO


CALL TimerStop(Timer_GR_SourceInput)
CALL TimerStop(Timer_GR_SourceInput_PartB)


END SUBROUTINE Poseidon_XCFC_Input_Sources2











END MODULE Source_Input_Module
