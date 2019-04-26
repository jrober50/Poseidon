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



USE Poseidon_Parameters, &
            ONLY :  DOMAIN_DIM,             &
                    DEGREE,                 &
                    L_LIMIT,                &
                    NUM_R_ELEMS_PER_SHELL,  &
                    NUM_R_ELEMS_PER_SUBSHELL,&
                    NUM_SHELLS,             &
                    NUM_SUBSHELLS,          &
                    NUM_BLOCKS,             &
                    NUM_SUBSHELLS_PER_SHELL,&
                    NUM_BLOCKS_PER_SHELL,   &
                    NUM_BLOCK_THETA_ROWS,   &
                    NUM_BLOCK_PHI_COLUMNS,  &
                    nPROCS_POSEIDON,        &
                    NUM_R_ELEMS_PER_BLOCK,  &
                    NUM_T_ELEMS_PER_BLOCK,  &
                    NUM_P_ELEMS_PER_BLOCK,  &
                    NUM_R_QUAD_POINTS,      &
                    NUM_T_QUAD_POINTS,      &
                    NUM_P_QUAD_POINTS,      &
                    NUM_QUAD_DOF,           &
                    R_COARSEN_FACTOR,       &
                    T_COARSEN_FACTOR,       &
                    P_COARSEN_FACTOR,       &
                    MAX_ITERATIONS,         &
                    CONVERGENCE_CRITERIA

USE Global_Variables_And_Parameters, &
            ONLY :  Coefficient_Vector,    &
                    PROB_DIM

IMPLICIT NONE

CONTAINS



 !+101+############################################################################!
!                                                                                   !
!                   UNPACK_PARAMETERS                                               !
!                                                                                   !
 !#################################################################################!
SUBROUTINE UNPACK_POSEIDON_PARAMETERS()

INTEGER                                 :: POSEIDON_read


INTEGER                                 :: iskipp
INTEGER                                 :: istat


INTEGER                                 :: DIMN, DEG, LLIM
INTEGER                                 :: PPROC, NSHELL, NSSHEL
INTEGER                                 :: NBPSHL, NBTROW, NBPCOL
INTEGER                                 :: NREPS, NREPSS
INTEGER                                 :: NTEPB, NPEPB
INTEGER                                 :: PRQ, PTQ, PPQ
INTEGER                                 :: RCF, TCF, PCF

INTEGER                                 :: MI
REAL(KIND = idp)                        :: CC




POSEIDON_read = 13



OPEN(UNIT=POSEIDON_read, FILE='Params/Poseidon_Params.d', STATUS='NEW', IOSTAT=istat)
IF ( istat .NE. 0 ) THEN
    OPEN(UNIT=POSEIDON_read, FILE='Params/Poseidon_Params.d', STATUS='OLD', IOSTAT=istat)
END IF
CALL READ_POSEIDON_PARAMETERS(POSEIDON_read, DIMN, DEG, LLIM, PPROC, NSHELL, NSSHEL,          &
                                             NBPSHL, NBTROW, NBPCOL, NREPS, NREPSS,     &
                                             NTEPB, NPEPB, PRQ, PTQ, PPQ, RCF, TCF, PCF,  &
                                             MI, CC  )
CLOSE(UNIT=POSEIDON_read,STATUS='keep',IOSTAT=istat)








DOMAIN_DIM                  = DIMN
DEGREE                      = DEG
L_LIMIT                     = LLIM

NUM_R_ELEMS_PER_SHELL       = NREPS
NUM_R_ELEMS_PER_SUBSHELL    = NREPSS
NUM_SHELLS                  = NSHELL
NUM_SUBSHELLS_PER_SHELL     = NSSHEL

NUM_BLOCKS_PER_SHELL        = NBPSHL
NUM_BLOCK_THETA_ROWS        = NBTROW
NUM_BLOCK_PHI_COLUMNS       = NBPCOL

nPROCS_POSEIDON             = PPROC

NUM_R_ELEMS_PER_BLOCK       = NREPS
NUM_T_ELEMS_PER_BLOCK       = NTEPB
NUM_P_ELEMS_PER_BLOCK       = NPEPB

NUM_R_QUAD_POINTS           = PRQ
NUM_T_QUAD_POINTS           = PTQ
NUM_P_QUAD_POINTS           = PPQ

R_COARSEN_FACTOR            = RCF
T_COARSEN_FACTOR            = TCF
P_COARSEN_FACTOR            = PCF


MAX_ITERATIONS              = MI
CONVERGENCE_CRITERIA        = CC

NUM_QUAD_DOF    = NUM_R_QUAD_POINTS   &
                * NUM_T_QUAD_POINTS   &
                * NUM_P_QUAD_POINTS

NUM_BLOCKS      = NUM_SHELLS*NUM_BLOCKS_PER_SHELL
NUM_SUBSHELLS   = NUM_SHELLS*NUM_SUBSHELLS_PER_SHELL



!PRINT*,"1",DOMAIN_DIM, DEGREE, L_LIMIT
!PRINT*,"2",NUM_R_ELEMS_PER_SHELL,NUM_R_ELEMS_PER_SUBSHELL, NUM_SHELLS,NUM_SUBSHELLS_PER_SHELL
!PRINT*,"3",NUM_BLOCKS_PER_SHELL,NUM_BLOCK_THETA_ROWS,NUM_BLOCK_PHI_COLUMNS
!PRINT*,"4",nPROCS_POSEIDON
!PRINT*,"5",NUM_R_ELEMS_PER_BLOCK,NUM_T_ELEMS_PER_BLOCK,NUM_P_ELEMS_PER_BLOCK
!PRINT*,"6",NUM_R_QUAD_POINTS,NUM_T_QUAD_POINTS,NUM_P_QUAD_POINTS
!PRINT*,"7",R_COARSEN_FACTOR,T_COARSEN_FACTOR,P_COARSEN_FACTOR
!PRINT*,"8",CC, MI


END SUBROUTINE UNPACK_POSEIDON_PARAMETERS










 !+102+############################################################################!
!                                                                                   !
!                  READ_POSEIDON_PARAMETERS                                         !
!                                                                                   !
 !#################################################################################!
SUBROUTINE READ_POSEIDON_PARAMETERS( nreadp, DIMN, DEG, LLIM, PPROC, NSHELL, NSSHEL,  &
                                     NBPSHL, NBTROW, NBPCOL, NREPS, NREPSS,           &
                                     NTEPB, NPEPB, PRQ, PTQ, PPQ, RCF, TCF, PCF,      &
                                     MI, CC  )

INTEGER, INTENT(IN)                                     :: nreadp
INTEGER, INTENT(INOUT)                                  :: DIMN, DEG, LLIM
INTEGER, INTENT(INOUT)                                  :: PPROC, NSHELL, NSSHEL
INTEGER, INTENT(INOUT)                                  :: NBPSHL, NBTROW, NBPCOL
INTEGER, INTENT(INOUT)                                  :: NREPS, NREPSS
INTEGER, INTENT(INOUT)                                  :: NTEPB, NPEPB
INTEGER, INTENT(INOUT)                                  :: PRQ, PTQ, PPQ
INTEGER, INTENT(INOUT)                                  :: RCF, TCF, PCF

INTEGER, INTENT(INOUT)                                  :: MI
REAL(KIND = idp), INTENT(INOUT)                         :: CC


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
    READ (line, 111) DIMN
    CYCLE
END IF

IF ( Param_type == 'DEGREE' ) THEN
    READ (line, 111) DEG
    CYCLE
END IF

IF ( Param_type == 'LLIMIT' ) THEN
    READ (line, 111) LLIM
    CYCLE
END IF



IF ( Param_type == 'PPROC ' ) THEN
    READ (line, 111) PPROC
    CYCLE
END IF






IF ( Param_type == 'NSHELL' ) THEN
    READ (line, 111) NSHELL
    CYCLE
END IF

IF ( Param_type == 'NSSHEL' ) THEN
    READ (line, 111) NSSHEL
    CYCLE
END IF







IF ( Param_type == 'NBPSHL' ) THEN
    READ (line, 111) NBPSHL
    CYCLE
END IF

IF ( Param_type == 'NBTROW' ) THEN
    READ (line, 111) NBTROW
    CYCLE
END IF


IF ( Param_type == 'NBPCOL' ) THEN
    READ (line, 111) NBPCOL
    CYCLE
END IF







IF ( Param_type == 'NREPS ' ) THEN
    READ (line, 111) NREPS
    CYCLE
END IF

IF ( Param_type == 'NREPSS' ) THEN
    READ (line, 111) NREPSS
    CYCLE
END IF
IF ( Param_type == 'NTEPB ' ) THEN
    READ (line, 111) NTEPB
    CYCLE
END IF

IF ( Param_type == 'NPEPB ' ) THEN
    READ (line, 111) NPEPB
    CYCLE
END IF




IF ( Param_type == 'PRQ   ' ) THEN
    READ (line, 111) PRQ
    CYCLE
END IF

IF ( Param_type == 'PTQ   ' ) THEN
    READ (line, 111) PTQ
    CYCLE
END IF

IF ( Param_type == 'PPQ   ' ) THEN
    READ (line, 111) PPQ
    CYCLE
END IF




IF ( Param_type == 'RCF   ' ) THEN
    READ (line, 111) RCF
    CYCLE
END IF

IF ( Param_type == 'TCF   ' ) THEN
    READ (line, 111) TCF
    CYCLE
END IF

IF ( Param_type == 'PCF   ' ) THEN
    READ (line, 111) PCF
    CYCLE
END IF


IF ( Param_type == 'MI    ' ) THEN
    READ (line, 111) MI
    CYCLE
END IF

IF ( Param_type == 'CC    ' ) THEN
    READ (line, 131) CC
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


INTEGER                                          :: i

INTEGER                                 :: Coeffs_Write


INTEGER                                 :: iskipp
INTEGER                                 :: istat

Coeffs_Write = 13

OPEN(UNIT=Coeffs_Write, FILE='Params/CFA_Coeffs.c', ACTION='WRITE', STATUS='REPLACE',IOSTAT=istat)

DO i = 0,PROB_DIM-1

    WRITE(Coeffs_Write,'(2ES24.17)') REAL(Coefficient_Vector(i)), AIMAG(Coefficient_Vector(i))

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

OPEN(UNIT=Coeffs_Read, FILE='Params/CFA_Coeffs.c', STATUS='NEW', IOSTAT=istat)
IF ( istat .NE. 0 ) THEN
    OPEN(UNIT=Coeffs_READ, FILE='Params/CFA_Coeffs.c', STATUS='OLD', IOSTAT=istat)
END IF


!READ(Coeffs_Read,*),TMP
DO i = 0,PROB_DIM-1
   READ(Coeffs_Read,'(2ES24.17)')  TMP(0,i),TMP(1,i)
END DO

CLOSE(UNIT=Coeffs_Read,STATUS='keep',IOSTAT=istat)


!Coefficient_Vector = 0.0_idp
!Coefficient_Vector(:) = CMPLX(TMP(0,:), TMP(1,:), KIND =idp)



END SUBROUTINE READ_CFA_COEFFICIENTS












END MODULE Poseidon_Parameter_Read_Module
