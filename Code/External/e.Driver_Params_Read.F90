   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE DRIVER_Params_Read_Module                                                   !##!
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


USE DRIVER_Parameters, &
            ONLY :  DRIVER_R_ELEMS,             &
                    DRIVER_C_ELEMS,             &
                    DRIVER_T_ELEMS,             &
                    DRIVER_P_ELEMS,             &
                    DRIVER_R_INPUT_NODES,       &
                    DRIVER_T_INPUT_NODES,       &
                    DRIVER_P_INPUT_NODES,       &
                    DRIVER_LEFT_LIMIT,          &
                    DRIVER_RIGHT_LIMIT,         &
                    DRIVER_INNER_RADIUS,        &
                    DRIVER_CORE_RADIUS,         &
                    DRIVER_OUTER_RADIUS,        &
                    DRIVER_TEST_NUMBER,         &
                    DRIVER_FIRST_GUESS_FLAG,    &
                    DRIVER_SUBSEQUENT_GUESS_FLAG,   &
                    DRIVER_DIMENSION,           &
                    DRIVER_PROCS,               &
                    DRIVER_y_PROCS,             &
                    DRIVER_z_PROCS,             &
                    DRIVER_MESH_TYPE,           &
                    DRIVER_Zoom,                &
                    DRIVER_SOLVER_TYPE,         &
                    SOLVER_MODE,                &
                    CFA_EQs_Flag_Vector,        &
                    RESULTS_OUTPUT_FLAG,        &
                    SOURCE_OUTPUT_FLAG,         &
                    RUN_REPORT_FLAG,            &
                    FRAME_REPORT_FLAG,          &
                    RHO_O,                      &
                    POWER_A,                    &
                    SELFSIM_START_T,            &
                    SELFSIM_END_T,              &
                    SELFSIM_NUM_FRAMES,         &
                    SELFSIM_KAPPA,              &
                    SELFSIM_GAMMA,              &
                    SELFSIM_ECC,                &
                    SELFSIM_V_SWITCH,           &
                    CHIMERA_START_FRAME,        &
                    CHIMERA_END_FRAME,          &
                    OUTPUT_PRIMATIVES_FLAG


IMPLICIT NONE

CONTAINS



 !+101+############################################################################!
!                                                                                   !
!                   UNPACK_PARAMETERS                                               !
!                                                                                   !
 !#################################################################################!
SUBROUTINE UNPACK_DRIVER_PARAMETERS()

INTEGER                                             ::  DRIVER_read


INTEGER                                             ::  iskipp
INTEGER                                             ::  istat

INTEGER                                             ::  Num_Int_Params  = 31
INTEGER                                             ::  Num_Real_Params = 12

INTEGER,            DIMENSION(:), ALLOCATABLE     ::  Int_Params
REAL(KIND = idp),   DIMENSION(:), ALLOCATABLE    ::  Real_Params
INTEGER,            DIMENSION(:), ALLOCATABLE    ::  Real_Params_Flags


ALLOCATE( Int_Params(1:Num_Int_Params) )
ALLOCATE( Real_Params(1:Num_Real_Params) )
ALLOCATE( Real_Params_Flags(1:Num_Real_Params) )

DRIVER_read = 42

OPEN(UNIT=DRIVER_read, FILE='Params/Driver_Params.d', STATUS='NEW', IOSTAT=istat)
IF ( istat .NE. 0 ) THEN
    OPEN(UNIT=DRIVER_read, FILE='Params/Driver_Params.d', STATUS='OLD', IOSTAT=istat)
END IF


! Set Initial Values !
Int_Params = -1
Real_Params_Flags = -1

CALL READ_DRIVER_PARAMETERS( DRIVER_read,           &
                              Int_Params,           &
                              Real_Params,          &
                              Real_Params_Flags,    &
                              Num_Int_Params,       &
                              Num_Real_Params       )

CLOSE(UNIT=DRIVER_read,STATUS='keep',IOSTAT=istat)




! Integers
DRIVER_DIMENSION            = Int_Params(1)     ! DIM
DRIVER_R_ELEMS              = Int_Params(2)     ! RE
DRIVER_C_ELEMS              = Int_Params(3)     ! CE
DRIVER_T_ELEMS              = Int_Params(4)     ! TE
DRIVER_P_ELEMS              = Int_Params(5)     ! PE
DRIVER_R_INPUT_NODES        = Int_Params(6)     ! RQ
DRIVER_T_INPUT_NODES        = Int_Params(7)     ! TQ
DRIVER_P_INPUT_NODES        = Int_Params(8)     ! PQ

DRIVER_MESH_TYPE            = Int_Params(9)     ! MT

DRIVER_PROCS                = Int_Params(10)    ! PROCS
DRIVER_y_PROCS              = Int_Params(11)    ! yPROCS
DRIVER_z_PROCS              = Int_Params(12)    ! zPROCS

IF ( Int_Params(13) .NE. -1 ) THEN
    DRIVER_TEST_NUMBER      = Int_Params(13)    ! DTN (Defaults is 3)
END IF

SOURCE_OUTPUT_FLAG          = Int_Params(14)    ! SOF
RESULTS_OUTPUT_FLAG         = Int_Params(15)    ! ROF
RUN_REPORT_FLAG             = Int_Params(16)    ! RRF
FRAME_REPORT_FLAG           = Int_Params(17)    ! FRF

CHIMERA_START_FRAME         = Int_Params(18)    ! CSF
CHIMERA_END_FRAME           = Int_Params(19)    ! CEF


POWER_A                     = Int_Params(20)    ! PWA

SELFSIM_NUM_FRAMES          = Int_Params(21)    ! YNF
SELFSIM_V_SWITCH            = Int_Params(22)    ! SSV

DRIVER_FIRST_GUESS_FLAG     = Int_Params(23)    ! FGF
DRIVER_SUBSEQUENT_GUESS_FLAG = Int_Params(24)   ! SGF

! Reals
DRIVER_LEFT_LIMIT          = Real_Params(1)     !   LL
DRIVER_RIGHT_LIMIT         = Real_Params(2)     !   RL
DRIVER_INNER_RADIUS        = Real_Params(3)     !   IR
DRIVER_CORE_RADIUS         = Real_Params(4)     !   CR
DRIVER_OUTER_RADIUS        = Real_Params(5)     !   OR

IF ( Real_Params_Flags(6) .NE.  -1 ) THEN
    RHO_O                   = Real_Params(6)    !   DEN
END IF


SELFSIM_START_T             = Real_Params(7)    !   YST
SELFSIM_END_T               = Real_Params(8)    !   YET
SELFSIM_KAPPA               = Real_Params(9)    !   SSK
SELFSIM_GAMMA               = Real_Params(10)   !   SSG
SELFSIM_ECC                 = Real_Params(11)   !   SSE
DRIVER_Zoom                 = Real_Params(12)   !   GZ

OUTPUT_PRIMATIVES_FLAG       = Int_Params(25)    !   OPF

SOLVER_MODE                 = Int_Params(26)    !  DSM

CFA_EQS_Flag_Vector(1)       = Int_Params(27)    !  EQF1
CFA_EQS_Flag_Vector(2)       = Int_Params(28)    !  EQF2
CFA_EQS_Flag_Vector(3)       = Int_Params(29)    !  EQF3
CFA_EQS_Flag_Vector(4)       = Int_Params(30)    !  EQF4
CFA_EQS_Flag_Vector(5)       = Int_Params(31)    !  EQF5



END SUBROUTINE UNPACK_DRIVER_PARAMETERS






 !+102+############################################################################!
!                                                                                   !
!                  READ_DRIVER_PARAMETERS                                          !
!                                                                                   !
 !#################################################################################!
SUBROUTINE READ_DRIVER_PARAMETERS( nreadp,                     &
                                    Int_Params,                 &
                                    Real_Params,                &
                                    Real_Params_Flags,          &
                                    Num_Int_Params,             &
                                    Num_Real_Params             )

INTEGER, INTENT(IN)                                 ::  nreadp

INTEGER, INTENT(IN)                                 ::  Num_Int_Params
INTEGER, INTENT(IN)                                 ::  Num_Real_Params

INTEGER,            DIMENSION(1:Num_Int_Params), INTENT(INOUT)  ::  Int_Params
REAL(KIND = idp),   DIMENSION(1:Num_Real_Params), INTENT(INOUT) ::  Real_Params
INTEGER,            DIMENSION(1:Num_Real_Params), INTENT(INOUT) ::  Real_Params_Flags


INTEGER                             :: iskipp
INTEGER                             :: istat



CHARACTER(LEN=6)                    :: Param_type
CHARACTER(LEN=20)                   :: Param_name
CHARACTER(LEN=128)                  :: line



101 FORMAT (a128)
111 FORMAT (10x,i10)
121 FORMAT (10x,F9.5)
131 FORMAT (10x,E14.0)
141 FORMAT (6x,F17.16)


REWIND(nreadp)
READ: DO

READ(nreadp,101,END=5000) line

Param_type = line(1:6)
Param_name = line(40:60)


IF ( Param_type == 'DIM   ' ) THEN
    READ (line, 111) INT_PARAMS(1)
    CYCLE
END IF

IF ( Param_type == 'RE    ' ) THEN
    READ (line, 111) INT_PARAMS(2)
    CYCLE
END IF

IF ( Param_type == 'CE    ' ) THEN
    READ (line, 111) INT_PARAMS(3)
    CYCLE
END IF


IF ( Param_type == 'TE    ' ) THEN
    READ (line, 111) INT_PARAMS(4)
    CYCLE
END IF

IF ( Param_type == 'PE    ' ) THEN
    READ (line, 111) INT_PARAMS(5)
    CYCLE
END IF




IF ( Param_type == 'RQ    ' ) THEN
    READ (line, 111) INT_PARAMS(6)
    CYCLE
END IF

IF ( Param_type == 'TQ    ' ) THEN
    READ (line, 111) INT_PARAMS(7)
    CYCLE
END IF

IF ( Param_type == 'PQ    ' ) THEN
    READ (line, 111) INT_PARAMS(8)
    CYCLE
END IF

IF ( Param_type == 'MT    ' ) THEN
    READ (line, 111) INT_PARAMS(9)
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



IF ( Param_type == 'DTN   ' ) THEN
    READ (line, 111) INT_PARAMS(13)
    CYCLE
END IF


IF ( Param_type == 'SOF ' ) THEN
    READ (line, 111) INT_PARAMS(14)
    CYCLE
END IF

IF ( Param_type == 'SOF ' ) THEN
    READ (line, 111) INT_PARAMS(15)
    CYCLE
END IF

IF ( Param_type == 'RRF ' ) THEN
    READ (line, 111) INT_PARAMS(16)
    CYCLE
END IF

IF ( Param_type == 'FRF ' ) THEN
    READ (line, 111) INT_PARAMS(17)
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



! CHIMERA Test Parameters

IF ( Param_type == 'CSF   ' ) THEN
    READ (line, 111) INT_PARAMS(18)
    CYCLE
END IF


IF ( Param_type == 'CEF   ' ) THEN
    READ (line, 111) INT_PARAMS(19)
    CYCLE
END IF



!  Spherically Symmetric Test Parameters

IF ( Param_type == 'DEN   ' ) THEN
    READ (line, 131) REAL_PARAMS(6)
    REAL_PARAMS_FLAGS(6) = 1
    CYCLE
END IF

IF ( Param_type == 'PWA   ' ) THEN
    READ (line, 111) INT_PARAMS(20)
    CYCLE
END IF



! Self Similar Test Parameters

IF ( Param_type == 'YST   ' ) THEN
    READ (line, 131) REAL_PARAMS(7)
    REAL_PARAMS_FLAGS(7) = 1
    CYCLE
END IF

IF ( Param_type == 'YET   ' ) THEN
    READ (line, 131) REAL_PARAMS(8)
    REAL_PARAMS_FLAGS(8) = 1
    CYCLE
END IF

IF ( Param_type == 'YNF   ' ) THEN
    READ (line, 111) INT_PARAMS(21)
    CYCLE
END IF


IF ( Param_type == 'SSK   ' ) THEN
    READ (line, 131) REAL_PARAMS(9)
    REAL_PARAMS_FLAGS(9) = 1
    CYCLE
END IF

IF ( Param_type == 'SSG   ' ) THEN
    READ (line, 131) REAL_PARAMS(10)
    REAL_PARAMS_FLAGS(10) = 1
    CYCLE
END IF

IF ( Param_type == 'SSE   ' ) THEN
    READ (line, 121) REAL_PARAMS(11)
    REAL_PARAMS_FLAGS(11) = 1
    CYCLE
END IF


IF ( Param_type == 'SSV   ' ) THEN
    READ (line, 111) INT_PARAMS(22)
    CYCLE
END IF

IF ( Param_type == 'FGF   ' ) THEN
    READ (line, 111) INT_PARAMS(23)
    CYCLE
END IF


IF ( Param_type == 'SGF   ' ) THEN
    READ (line, 111) INT_PARAMS(24)
    CYCLE
END IF

IF ( Param_type == 'GZ    ' ) THEN
    READ (line, 141) REAL_PARAMS(12)
    CYCLE
END IF


IF ( Param_type == 'OPF   ' ) THEN
    READ (line, 111) INT_PARAMS(25)
    CYCLE
END IF


IF ( Param_type == 'DSM   ' ) THEN
    READ ( line,111) INT_PARAMS(26)
    CYCLE
ENDIF

IF ( Param_type == 'EQF1  ' ) THEN
    READ ( line,111) INT_PARAMS(27)
    CYCLE
ENDIF

IF ( Param_type == 'EQF2  ' ) THEN
    READ ( line,111) INT_PARAMS(28)
    CYCLE
ENDIF

IF ( Param_type == 'EQF3  ' ) THEN
    READ ( line,111) INT_PARAMS(29)
    CYCLE
ENDIF

IF ( Param_type == 'EQF4  ' ) THEN
    READ ( line,111) INT_PARAMS(30)
    CYCLE
ENDIF

IF ( Param_type == 'EQF5  ' ) THEN
    READ ( line,111) INT_PARAMS(31)
    CYCLE
ENDIF


END DO READ


 5000 RETURN
END SUBROUTINE READ_DRIVER_PARAMETERS












END MODULE DRIVER_Params_Read_Module
