   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE CHIMERA_Params_Read_Module                                                   !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
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
            ONLY : idp


USE CHIMERA_Parameters, &
            ONLY :  CHIMERA_R_ELEMS,        &
                    CHIMERA_C_ELEMS,        &
                    CHIMERA_T_ELEMS,        &
                    CHIMERA_P_ELEMS,        &
                    CHIMERA_R_INPUT_NODES,  &
                    CHIMERA_T_INPUT_NODES,  &
                    CHIMERA_P_INPUT_NODES,  &
                    CHIMERA_LEFT_LIMIT,     &
                    CHIMERA_RIGHT_LIMIT,    &
                    CHIMERA_INNER_RADIUS,   &
                    CHIMERA_CORE_RADIUS,    &
                    CHIMERA_OUTER_RADIUS,   &
                    CHIMERA_DIMENSION,      &
                    CHIMERA_PROCS,          &
                    CHIMERA_y_PROCS,        &
                    CHIMERA_z_PROCS,        &
                    CHIMERA_MESH_TYPE,      &
                    CHIMERA_SOLVER_TYPE,    &
                    RHO_O,                  &
                    POWER_A,                &
                    SELFSIM_T,              &
                    SELFSIM_KAPPA,          &
                    SELFSIM_GAMMA


IMPLICIT NONE

CONTAINS



 !+101+############################################################################!
!                                                                                   !
!                   UNPACK_PARAMETERS                                               !
!                                                                                   !
 !#################################################################################!
SUBROUTINE UNPACK_CHIMERA_PARAMETERS()

INTEGER                                 ::  CHIMERA_read


INTEGER                                 ::  iskipp
INTEGER                                 ::  istat

INTEGER                                 ::  RE, CE, TE, PE
INTEGER                                 ::  MT, ST
INTEGER                                 ::  RQ, TQ, PQ
INTEGER                                 ::  DIM
INTEGER                                 ::  PROCS, yPROCS, zPROCS
REAL(KIND = idp)                        ::  LL, RL, IR, CR, OR
REAL(KIND = idp)                        ::  DEN
INTEGER                                 ::  PWA
REAL(KIND = idp)                        ::  SST, SSK, SSG







CHIMERA_read = 42

OPEN(UNIT=CHIMERA_read, FILE='Params/CHIMERA_Params.d', STATUS='NEW', IOSTAT=istat)
IF ( istat .NE. 0 ) THEN
    OPEN(UNIT=CHIMERA_read, FILE='Params/CHIMERA_Params.d', STATUS='OLD', IOSTAT=istat)
END IF

CALL READ_CHIMERA_PARAMETERS( CHIMERA_read,                     &
                                RE, CE, TE, PE, RQ, TQ, PQ,     &
                                DIM, PROCS, yPROCS, zPROCS,     &
                                LL, RL, IR, CR, OR, MT, ST,     &
                                DEN, PWA,                       &
                                SST, SSK, SSG                   )

CLOSE(UNIT=CHIMERA_read,STATUS='keep',IOSTAT=istat)








CHIMERA_R_ELEMS             = RE
CHIMERA_C_ELEMS             = CE
CHIMERA_T_ELEMS             = TE
CHIMERA_P_ELEMS             = PE
CHIMERA_R_INPUT_NODES       = RQ
CHIMERA_T_INPUT_NODES       = TQ
CHIMERA_P_INPUT_NODES       = PQ
CHIMERA_LEFT_LIMIT          = LL
CHIMERA_RIGHT_LIMIT         = RL
CHIMERA_MESH_TYPE           = MT
CHIMERA_INNER_RADIUS        = IR
CHIMERA_CORE_RADIUS         = CR
CHIMERA_OUTER_RADIUS        = OR
CHIMERA_DIMENSION           = DIM
CHIMERA_PROCS               = PROCS
CHIMERA_y_PROCS             = yPROCS
CHIMERA_z_PROCS             = zPROCS

CHIMERA_SOLVER_TYPE         = ST

RHO_O                       = DEN
POWER_A                     = PWA

SELFSIM_T                   = SST
SELFSIM_KAPPA               = SSK
SELFSIM_GAMMA               = SSG





END SUBROUTINE UNPACK_CHIMERA_PARAMETERS






 !+102+############################################################################!
!                                                                                   !
!                  READ_CHIMERA_PARAMETERS                                          !
!                                                                                   !
 !#################################################################################!
SUBROUTINE READ_CHIMERA_PARAMETERS( nreadp,                             &
                                    RE, CE, TE, PE, RQ, TQ, PQ,         &
                                    DIM, PROCS, yPROCS, zPROCS,         &
                                    LL, RL, IR, CR, OR, MT, ST,         &
                                    DEN, PWA,                           &
                                    SST, SSK, SSG                       )

INTEGER, INTENT(IN)                                 ::  nreadp
INTEGER,INTENT(INOUT)                               ::  RE, CE, TE, PE
INTEGER,INTENT(INOUT)                               ::  RQ, TQ, PQ
INTEGER,INTENT(INOUT)                               ::  DIM
INTEGER,INTENT(INOUT)                               ::  MT, ST
INTEGER,INTENT(INOUT)                               ::  PROCS, yPROCS, zPROCS
REAL(KIND = idp), INTENT(INOUT)                     ::  LL, RL, IR, CR, OR
REAL(KIND = idp), INTENT(INOUT)                     ::  DEN
INTEGER,INTENT(INOUT)                               ::  PWA
REAL(KIND = idp), INTENT(INOUT)                     ::  SST, SSK, SSG

INTEGER                             :: iskipp
INTEGER                             :: istat



CHARACTER(LEN=6)                    :: Param_type
CHARACTER(LEN=20)                   :: Param_name
CHARACTER(LEN=128)                  :: line



101 FORMAT (a128)
111 FORMAT (10x,i10)
121 FORMAT (10x,f2.1)
131 FORMAT (10x,E14.0)



REWIND(nreadp)
READ: DO

READ(nreadp,101,END=5000) line

Param_type = line(1:6)
Param_name = line(40:60)

IF ( Param_type == 'DIM   ' ) THEN
    READ (line, 111) DIM
    CYCLE
END IF

IF ( Param_type == 'RE    ' ) THEN
    READ (line, 111) RE
    CYCLE
END IF

IF ( Param_type == 'CE    ' ) THEN
    READ (line, 111) CE
    CYCLE
END IF


IF ( Param_type == 'TE    ' ) THEN
    READ (line, 111) TE
    CYCLE
END IF

IF ( Param_type == 'PE    ' ) THEN
    READ (line, 111) PE
    CYCLE
END IF




IF ( Param_type == 'RQ    ' ) THEN
    READ (line, 111) RQ
    CYCLE
END IF

IF ( Param_type == 'TQ    ' ) THEN
    READ (line, 111) TQ
    CYCLE
END IF

IF ( Param_type == 'PQ    ' ) THEN
    READ (line, 111) PQ
    CYCLE
END IF





IF ( Param_type == 'PROC  ' ) THEN
    READ (line, 111) PROCS
    CYCLE
END IF


IF ( Param_type == 'yPROC ' ) THEN
    READ (line, 111) yPROCS
    CYCLE
END IF


IF ( Param_type == 'zPROC ' ) THEN
    READ (line, 111) zPROCS
    CYCLE
END IF




IF ( Param_type == 'LL    ' ) THEN
    READ (line, 131) LL
    CYCLE
END IF


IF ( Param_type == 'RL    ' ) THEN
    READ (line, 131) RL
    CYCLE
END IF

IF ( Param_type == 'IR    ' ) THEN
    READ (line, 131) IR
    CYCLE
END IF


IF ( Param_type == 'CR    ' ) THEN
    READ (line, 131) CR
    CYCLE
END IF

IF ( Param_type == 'OR    ' ) THEN
    READ (line, 131) OR
    CYCLE
END IF


IF ( Param_type == 'MT    ' ) THEN
    READ (line, 111) MT
    CYCLE
END IF


IF ( Param_type == 'ST    ' ) THEN
    READ (line, 111) ST
    CYCLE
END IF

IF ( Param_type == 'DEN   ' ) THEN
     READ (line, 131) DEN
     CYCLE
END IF

IF ( Param_type == 'PWA   ' ) THEN
     READ (line, 111) PWA
     CYCLE
END IF



IF ( Param_type == 'SST   ' ) THEN
     READ (line, 131) SST
     CYCLE
END IF

IF ( Param_type == 'SSK   ' ) THEN
     READ (line, 131) SSK
     CYCLE
END IF

IF ( Param_type == 'SSG   ' ) THEN
     READ (line, 131) SSG
     CYCLE
END IF



END DO READ


 5000 RETURN
END SUBROUTINE READ_CHIMERA_PARAMETERS












END MODULE CHIMERA_Params_Read_Module
