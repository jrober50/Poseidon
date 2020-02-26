   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_Parameter_Read_Module                                               !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   UNPACK_POSEIDON_PARAMETERS                                          !##!
!##!    +102+    READ_POSEIDON_PARAMETERS                                           !##!
!##!                                                                                !##!
!##!    +201+    WRITE_CFA_COEFFICIENTS                                             !##!
!##!    +202+    READ_CFA_COEFFICIENTS                                              !##!
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



USE Poseidon_Parameters, &
            ONLY :  DOMAIN_DIM,                 &
                    DEGREE,                     &
                    L_LIMIT,                    &
                    NUM_CFA_VARS,               &
                    NUM_R_ELEMS_PER_SHELL,      &
                    NUM_R_ELEMS_PER_SUBSHELL,   &
                    NUM_SHELLS,                 &
                    NUM_SUBSHELLS,              &
                    NUM_BLOCKS,                 &
                    NUM_SUBSHELLS_PER_SHELL,    &
                    NUM_BLOCKS_PER_SHELL,       &
                    NUM_BLOCK_THETA_ROWS,       &
                    NUM_BLOCK_PHI_COLUMNS,      &
                    nPROCS_POSEIDON,            &
                    NUM_R_ELEMS_PER_BLOCK,      &
                    NUM_T_ELEMS_PER_BLOCK,      &
                    NUM_P_ELEMS_PER_BLOCK,      &
                    NUM_R_QUAD_POINTS,          &
                    NUM_T_QUAD_POINTS,          &
                    NUM_P_QUAD_POINTS,          &
                    NUM_QUAD_DOF,               &
                    R_COARSEN_FACTOR,           &
                    T_COARSEN_FACTOR,           &
                    P_COARSEN_FACTOR,           &
                    MAX_ITERATIONS,             &
                    CONVERGENCE_CRITERIA,       &
                    OUTPUT_MATRIX_FLAG,         &
                    OUTPUT_RHS_VECTOR_FLAG,     &
                    WRITE_TIMETABLE_FLAG,       &
                    WRITE_REPORT_FLAG,          &
                    ITER_REPORT_NUM_SAMPLES,    &
                    WRITE_RESULTS_FLAG,         &
                    NEW_PETSC_SOLVER_FLAG


USE Poseidon_Variables_Module, &
            ONLY :  Coefficient_Vector,         &
                    ULM_LENGTH,                 &
                    NUM_R_ELEMENTS,             &
                    PROB_DIM

IMPLICIT NONE

CONTAINS



 !+101+############################################################################!
!                                                                                   !
!                   UNPACK_PARAMETERS                                               !
!                                                                                   !
 !#################################################################################!
SUBROUTINE UNPACK_POSEIDON_PARAMETERS()

INTEGER                                         ::  POSEIDON_read


INTEGER                                         ::  iskipp
INTEGER                                         ::  istat

INTEGER                                         ::  NUM_INT_PARAMS
INTEGER                                         ::  NUM_REAL_PARAMS

INTEGER,          DIMENSION(:), ALLOCATABLE     :: INT_PARAMS
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE     :: REAL_PARAMS



NUM_INT_PARAMS  = 27
NUM_REAL_PARAMS = 1

ALLOCATE( INT_PARAMS(1:NUM_INT_PARAMS) )
ALLOCATE( REAL_PARAMS(1:NUM_REAL_PARAMS) )

INT_PARAMS = -1

POSEIDON_read = 13
OPEN(UNIT=POSEIDON_read, FILE='Params/Poseidon_Params.d', STATUS='NEW', IOSTAT=istat)
IF ( istat .NE. 0 ) THEN
    OPEN(UNIT=POSEIDON_read, FILE='Params/Poseidon_Params.d', STATUS='OLD', IOSTAT=istat)
END IF



CALL READ_POSEIDON_PARAMETERS( POSEIDON_read, NUM_INT_PARAMS, NUM_REAL_PARAMS,      &
                                              INT_PARAMS, REAL_PARAMS               )


CLOSE(UNIT=POSEIDON_read,STATUS='keep',IOSTAT=istat)





DOMAIN_DIM                  = INT_PARAMS(1)
DEGREE                      = INT_PARAMS(2)
L_LIMIT                     = INT_PARAMS(3)

nPROCS_POSEIDON             = INT_PARAMS(4)

NUM_R_ELEMS_PER_SHELL       = INT_PARAMS(5)
NUM_R_ELEMS_PER_SUBSHELL    = INT_PARAMS(6)
NUM_SHELLS                  = INT_PARAMS(7)
NUM_SUBSHELLS_PER_SHELL     = INT_PARAMS(8)

NUM_BLOCKS_PER_SHELL        = INT_PARAMS(9)
NUM_BLOCK_THETA_ROWS        = INT_PARAMS(10)
NUM_BLOCK_PHI_COLUMNS       = INT_PARAMS(11)


NUM_R_ELEMS_PER_BLOCK       = INT_PARAMS(5)
NUM_T_ELEMS_PER_BLOCK       = INT_PARAMS(12)
NUM_P_ELEMS_PER_BLOCK       = INT_PARAMS(13)

NUM_R_QUAD_POINTS           = INT_PARAMS(14)
NUM_T_QUAD_POINTS           = INT_PARAMS(15)
NUM_P_QUAD_POINTS           = INT_PARAMS(16)


IF ( INT_PARAMS(17) .NE. -1 ) THEN
    R_COARSEN_FACTOR            = INT_PARAMS(17) ! Default = 1
END IF
IF ( INT_PARAMS(18) .NE. -1 ) THEN
    T_COARSEN_FACTOR            = INT_PARAMS(18) ! Default = 1
END IF
IF ( INT_PARAMS(19) .NE. -1 ) THEN
    P_COARSEN_FACTOR            = INT_PARAMS(19) ! Default = 1
END IF



MAX_ITERATIONS                  = INT_PARAMS(20)



IF ( INT_PARAMS(21) .NE. -1 ) THEN
    WRITE_TIMETABLE_FLAG        = INT_PARAMS(21)    ! Default = 0, Off
END IF
IF ( INT_PARAMS(22) .NE. -1 ) THEN
    WRITE_REPORT_FLAG           = INT_PARAMS(22)    ! Default = 0, Off
END IF
IF ( INT_PARAMS(23) .NE. -1 ) THEN
    ITER_REPORT_NUM_SAMPLES     = INT_PARAMS(23)    ! Default = 20
END IF
IF ( INT_PARAMS(24) .NE. -1 ) THEN
    WRITE_RESULTS_FLAG          = INT_PARAMS(24)    ! Default = 0, Off
END IF


NEW_PETSC_SOLVER_FLAG           = INT_PARAMS(25)


OUTPUT_MATRIX_FLAG              = INT_PARAMS(26)
OUTPUT_RHS_VECTOR_FLAG          = INT_PARAMS(27)



CONVERGENCE_CRITERIA            = REAL_PARAMS(1)




NUM_QUAD_DOF    = NUM_R_QUAD_POINTS   &
                * NUM_T_QUAD_POINTS   &
                * NUM_P_QUAD_POINTS

NUM_BLOCKS      = NUM_SHELLS*NUM_BLOCKS_PER_SHELL
NUM_SUBSHELLS   = NUM_SHELLS*NUM_SUBSHELLS_PER_SHELL






END SUBROUTINE UNPACK_POSEIDON_PARAMETERS










 !+102+############################################################################!
!                                                                                   !
!                  READ_POSEIDON_PARAMETERS                                         !
!                                                                                   !
 !#################################################################################!
SUBROUTINE READ_POSEIDON_PARAMETERS( nreadp, NUM_INT_PARAMS, NUM_REAL_PARAMS,       &
                                             INT_PARAMS,     REAL_PARAMS            )

INTEGER,                                        INTENT(IN)          ::  nreadp
INTEGER,                                        INTENT(IN)          ::  NUM_INT_PARAMS
INTEGER,                                        INTENT(IN)          ::  NUM_REAL_PARAMS

INTEGER,          DIMENSION(1:NUM_INT_PARAMS),  INTENT(INOUT)       ::  INT_PARAMS
REAL(KIND = idp), DIMENSION(1:NUM_REAL_PARAMS), INTENT(INOUT)       ::  REAL_PARAMS




INTEGER                             :: iskipp
INTEGER                             :: istat



CHARACTER(LEN=6)                    :: Param_type
CHARACTER(LEN=20)                   :: Param_name
CHARACTER(LEN=128)                  :: line



101 FORMAT (a128)
111 FORMAT (10x,i10)
131 FORMAT (10x,E14.0)

REWIND(nreadp)
READ: DO

READ(nreadp,101,END=5000) line

Param_type = line(1:6)
Param_name = line(40:60)

IF ( Param_type == 'DIM' ) THEN
    READ (line, 111) INT_PARAMS(1)
    CYCLE
END IF

IF ( Param_type == 'DEGREE' ) THEN
    READ (line, 111) INT_PARAMS(2)
    CYCLE
END IF

IF ( Param_type == 'LLIMIT' ) THEN
    READ (line, 111) INT_PARAMS(3)
    CYCLE
END IF



IF ( Param_type == 'PPROC ' ) THEN
    READ (line, 111) INT_PARAMS(4)
    CYCLE
END IF






IF ( Param_type == 'NSHELL' ) THEN
    READ (line, 111) INT_PARAMS(7)
    CYCLE
END IF

IF ( Param_type == 'NSSHEL' ) THEN
    READ (line, 111) INT_PARAMS(8)
    CYCLE
END IF







IF ( Param_type == 'NBPSHL' ) THEN
    READ (line, 111) INT_PARAMS(9)
    CYCLE
END IF

IF ( Param_type == 'NBTROW' ) THEN
    READ (line, 111) INT_PARAMS(10)
    CYCLE
END IF


IF ( Param_type == 'NBPCOL' ) THEN
    READ (line, 111) INT_PARAMS(11)
    CYCLE
END IF







IF ( Param_type == 'NREPS ' ) THEN
    READ (line, 111) INT_PARAMS(5)
    CYCLE
END IF

IF ( Param_type == 'NREPSS' ) THEN
    READ (line, 111) INT_PARAMS(6)
    CYCLE
END IF
IF ( Param_type == 'NTEPB ' ) THEN
    READ (line, 111) INT_PARAMS(12)
    CYCLE
END IF

IF ( Param_type == 'NPEPB ' ) THEN
    READ (line, 111) INT_PARAMS(13)
    CYCLE
END IF




IF ( Param_type == 'PRQ   ' ) THEN
    READ (line, 111) INT_PARAMS(14)
    CYCLE
END IF

IF ( Param_type == 'PTQ   ' ) THEN
    READ (line, 111) INT_PARAMS(15)
    CYCLE
END IF

IF ( Param_type == 'PPQ   ' ) THEN
    READ (line, 111) INT_PARAMS(16)
    CYCLE
END IF




IF ( Param_type == 'RCF   ' ) THEN
    READ (line, 111) INT_PARAMS(17)
    CYCLE
END IF

IF ( Param_type == 'TCF   ' ) THEN
    READ (line, 111) INT_PARAMS(18)
    CYCLE
END IF

IF ( Param_type == 'PCF   ' ) THEN
    READ (line, 111) INT_PARAMS(19)
    CYCLE
END IF


IF ( Param_type == 'MI    ' ) THEN
    READ (line, 111) INT_PARAMS(20)
    CYCLE
END IF




IF ( Param_type == 'WRTTT' ) THEN
    READ (line, 111) INT_PARAMS(21)
    CYCLE
END IF

IF ( Param_type == 'WRTIR' ) THEN
    READ (line, 111) INT_PARAMS(22)
    CYCLE
END IF

IF ( Param_type == 'IRNS ' ) THEN
    READ (line, 111) INT_PARAMS(23)
    CYCLE
END IF

IF ( Param_type == 'WRTRS' ) THEN
    READ (line, 111) INT_PARAMS(24)
    CYCLE
END IF



IF ( Param_type == 'CC    ' ) THEN
    READ (line, 131) REAL_PARAMS(1)
    CYCLE
END IF


IF ( Param_type == 'NPS   ' ) THEN
    READ (line, 111) INT_PARAMS(25)
    CYCLE
END IF

IF ( Param_type == 'OMF   ' ) THEN
    READ (line, 111) INT_PARAMS(26)
    CYCLE
END IF

IF ( Param_type == 'ORF   ' ) THEN
    READ (line, 111) INT_PARAMS(27)
    CYCLE
END IF




END DO READ

 5000 RETURN
END SUBROUTINE READ_POSEIDON_PARAMETERS











 !+201+############################################################################!
!                                                                                   !
!                    WRITE_CFA_COEFFICIENTS                                         !
!                                                                                   !
 !#################################################################################!
SUBROUTINE WRITE_CFA_COEFFICIENTS( )


INTEGER                                             :: i
INTEGER                                             ::  re, d, l, m, u, here

INTEGER                                             :: Coeffs_Write


INTEGER                                             :: iskipp
INTEGER                                             :: istat

Coeffs_Write = 13

OPEN(UNIT=Coeffs_Write, FILE='Params/CFA_Coeffs.coefs', ACTION='WRITE', STATUS='REPLACE',IOSTAT=istat)

!DO i = 0,PROB_DIM-1
!    WRITE(Coeffs_Write,'(2ES24.17)') REAL(Coefficient_Vector(i)), AIMAG(Coefficient_Vector(i))
!END DO

DO re = 0,NUM_R_ELEMENTS-1
    DO d = 0,DEGREE
        DO l = 0,L_LIMIT
            DO m = -l,l
                here = (re*DEGREE+d)*ULM_LENGTH+(l*(l+1)+m)*NUM_CFA_VARS
                WRITE(Coeffs_Write,'(A2,I1,A3,I2)')"L=",l,",M=",l
                WRITE(Coeffs_Write,'(2ES24.17)')Coefficient_Vector(here:here+4)
            END DO
            WRITE(Coeffs_Write,'(/ /)')
        END DO
    END DO
END DO

CLOSE(UNIT=Coeffs_Write,STATUS='keep',IOSTAT=istat)





END SUBROUTINE WRITE_CFA_COEFFICIENTS






 !+202+############################################################################!
!                                                                                   !
!                     READ_CFA_COEFFICIENTS                                         !
!                                                                                   !
 !#################################################################################!
SUBROUTINE READ_CFA_COEFFICIENTS( )


INTEGER                                          :: i

INTEGER                                 :: Coeffs_read



INTEGER                                 :: iskipp
INTEGER                                 :: istat

REAL(KIND=idp), DIMENSION(0:1,0:PROB_DIM-1)        :: TMP

Coeffs_Read = 13

OPEN(UNIT=Coeffs_Read, FILE='Params/CFA_Coeffs.coefs', STATUS='NEW', IOSTAT=istat)
IF ( istat .NE. 0 ) THEN
    OPEN(UNIT=Coeffs_READ, FILE='Params/CFA_Coeffs.coefs', STATUS='OLD', IOSTAT=istat)
END IF


!READ(Coeffs_Read,*),TMP
DO i = 0,PROB_DIM-1
   READ(Coeffs_Read,'(2ES24.17)')  TMP(0,i),TMP(1,i)
END DO

CLOSE(UNIT=Coeffs_Read,STATUS='keep',IOSTAT=istat)


!Coefficient_Vector = 0.0_idp
Coefficient_Vector(:) = CMPLX(TMP(0,:), TMP(1,:), KIND =idp)



END SUBROUTINE READ_CFA_COEFFICIENTS












END MODULE Poseidon_Parameter_Read_Module
