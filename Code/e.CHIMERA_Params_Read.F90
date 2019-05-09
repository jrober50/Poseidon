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

INTEGER,            DIMENSION(1:14)     ::  Int_Params
REAL(KIND = idp),   DIMENSION(1:9)      ::  Real_Params
INTEGER,            DIMENSION(1:9)      ::  Real_Params_Flags



CHIMERA_read = 42

OPEN(UNIT=CHIMERA_read, FILE='Params/Driver_Params.d', STATUS='NEW', IOSTAT=istat)
IF ( istat .NE. 0 ) THEN
    OPEN(UNIT=CHIMERA_read, FILE='Params/Driver_Params.d', STATUS='OLD', IOSTAT=istat)
END IF


! Set Initial Values !
Int_Params = -1
Real_Params_Flags = -1

CALL READ_CHIMERA_PARAMETERS( CHIMERA_read,         &
                              Int_Params,           &
                              Real_Params,          &
                              Real_Params_Flags     )

CLOSE(UNIT=CHIMERA_read,STATUS='keep',IOSTAT=istat)



! Integers
CHIMERA_R_ELEMS             = Int_Params(1) ! RE
CHIMERA_C_ELEMS             = Int_Params(2) ! CE
CHIMERA_T_ELEMS             = Int_Params(3) ! TE
CHIMERA_P_ELEMS             = Int_Params(4) ! PE
CHIMERA_R_INPUT_NODES       = Int_Params(5) ! RQ
CHIMERA_T_INPUT_NODES       = Int_Params(6) ! TQ
CHIMERA_P_INPUT_NODES       = Int_Params(7) ! PQ

CHIMERA_MESH_TYPE           = Int_Params(8) ! MT
IF ( Int_Params(13) .NE. -1 ) THEN
    CHIMERA_DIMENSION       = Int_Params(9) ! DIM
END IF
CHIMERA_PROCS               = Int_Params(10) ! PROCS
CHIMERA_y_PROCS             = Int_Params(11) ! yPROCS
CHIMERA_z_PROCS             = Int_Params(12) ! zPROCS

IF ( Int_Params(13) .NE. -1 ) THEN
    CHIMERA_SOLVER_TYPE     = Int_Params(13) ! ST
END IF

POWER_A                     = Int_Params(14) ! PWA


! Reals
CHIMERA_LEFT_LIMIT          = Real_Params(1) ! LL
CHIMERA_RIGHT_LIMIT         = Real_Params(2) ! RL
CHIMERA_INNER_RADIUS        = Real_Params(3) ! IR
CHIMERA_CORE_RADIUS         = Real_Params(4) ! CR
CHIMERA_OUTER_RADIUS        = Real_Params(5) ! OR

IF ( Real_Params_Flags(6) .NE.  -1 ) THEN
    RHO_O                   = Real_Params(6) ! DEN
END IF


SELFSIM_T                   = Real_Params(7) ! SST
SELFSIM_KAPPA               = Real_Params(8) ! SSK
SELFSIM_GAMMA               = Real_Params(9) ! SSG





END SUBROUTINE UNPACK_CHIMERA_PARAMETERS






 !+102+############################################################################!
!                                                                                   !
!                  READ_CHIMERA_PARAMETERS                                          !
!                                                                                   !
 !#################################################################################!
SUBROUTINE READ_CHIMERA_PARAMETERS( nreadp,                     &
                                    Int_Params,                 &
                                    Real_Params,                &
                                    Real_Params_Flags           )

INTEGER, INTENT(IN)                                 ::  nreadp

INTEGER,            DIMENSION(1:14),INTENT(INOUT)   ::  Int_Params
REAL(KIND = idp),   DIMENSION(1:9), INTENT(INOUT)   ::  Real_Params
INTEGER,            DIMENSION(1:9), INTENT(INOUT)   ::  Real_Params_Flags


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
    READ (line, 111) INT_PARAMS(9)
    CYCLE
END IF

IF ( Param_type == 'RE    ' ) THEN
    READ (line, 111) INT_PARAMS(1)
    CYCLE
END IF

IF ( Param_type == 'CE    ' ) THEN
    READ (line, 111) INT_PARAMS(2)
    CYCLE
END IF


IF ( Param_type == 'TE    ' ) THEN
    READ (line, 111) INT_PARAMS(3)
    CYCLE
END IF

IF ( Param_type == 'PE    ' ) THEN
    READ (line, 111) INT_PARAMS(4)
    CYCLE
END IF




IF ( Param_type == 'RQ    ' ) THEN
    READ (line, 111) INT_PARAMS(5)
    CYCLE
END IF

IF ( Param_type == 'TQ    ' ) THEN
    READ (line, 111) INT_PARAMS(6)
    CYCLE
END IF

IF ( Param_type == 'PQ    ' ) THEN
    READ (line, 111) INT_PARAMS(7)
    CYCLE
END IF





IF ( Param_type == 'PROC  ' ) THEN
    READ (line, 111) INT_PARAMS(10)
    CYCLE
END IF


IF ( Param_type == 'yPROC ' ) THEN
    READ (line, 111) INT_PARAMS(11)
    CYCLE
END IF


IF ( Param_type == 'zPROC ' ) THEN
    READ (line, 111) INT_PARAMS(12)
    CYCLE
END IF




IF ( Param_type == 'LL    ' ) THEN
    READ (line, 131) REAL_PARAMS(1)
    REAL_PARAMS_FLAGS(1) = 1
    CYCLE
END IF


IF ( Param_type == 'RL    ' ) THEN
    READ (line, 131) REAL_PARAMS(2)
    REAL_PARAMS_FLAGS(2) = 1
    CYCLE
END IF

IF ( Param_type == 'IR    ' ) THEN
    READ (line, 131) REAL_PARAMS(3)
    REAL_PARAMS_FLAGS(3) = 1
    CYCLE
END IF


IF ( Param_type == 'CR    ' ) THEN
    READ (line, 131) REAL_PARAMS(4)
    REAL_PARAMS_FLAGS(4) = 1
    CYCLE
END IF

IF ( Param_type == 'OR    ' ) THEN
    READ (line, 131) REAL_PARAMS(5)
    REAL_PARAMS_FLAGS(5) = 1
    CYCLE
END IF


IF ( Param_type == 'MT    ' ) THEN
    READ (line, 111) INT_PARAMS(8)
    CYCLE
END IF


IF ( Param_type == 'ST    ' ) THEN
    READ (line, 111) INT_PARAMS(13)
    CYCLE
END IF

IF ( Param_type == 'DEN   ' ) THEN
    READ (line, 131) REAL_PARAMS(6)
    REAL_PARAMS_FLAGS(6) = 1
    CYCLE
END IF

IF ( Param_type == 'PWA   ' ) THEN
    READ (line, 111) INT_PARAMS(14)
    CYCLE
END IF



IF ( Param_type == 'SST   ' ) THEN
    READ (line, 131) REAL_PARAMS(7)
    REAL_PARAMS_FLAGS(7) = 1
    CYCLE
END IF

IF ( Param_type == 'SSK   ' ) THEN
    READ (line, 131) REAL_PARAMS(8)
    REAL_PARAMS_FLAGS(8) = 1
    CYCLE
END IF

IF ( Param_type == 'SSG   ' ) THEN
    READ (line, 131) REAL_PARAMS(9)
    REAL_PARAMS_FLAGS(9) = 1
    CYCLE
END IF



END DO READ


 5000 RETURN
END SUBROUTINE READ_CHIMERA_PARAMETERS












END MODULE CHIMERA_Params_Read_Module
