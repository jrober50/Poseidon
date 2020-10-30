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

USE Variables_IO, &
            ONLY :  WRITE_SOURCES_FLAG

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


!Tmp_Si = 0.0_idp
!PRINT*,"Si zeroed in Poseidon_Input_Sources, z.Poseidon_Source_Module.F90"

PRINT*,"Before OUTPUT_POSEIDON_SOURCES_1D"

CALL OUTPUT_POSEIDON_SOURCES_1D(Local_E, Local_S, Local_Si,                         &
                                Local_RE_Dim, Local_TE_Dim, Local_PE_Dim,           &
                                Local_RQ_Dim, Local_TQ_Dim, Local_PQ_Dim,           &
                                Input_R_Quad, Input_T_Quad, Input_P_Quad,           &
                                Left_Limit, Right_Limit                             )
    


! Poseidon_CFA_Block_Share takes source data and redistributes it into Poseidon's preferred decomposition.

PRINT*,"Before Poseidon_CFA_Block_Share"
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







END MODULE Source_Input_Module
