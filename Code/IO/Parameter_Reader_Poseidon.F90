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
USE Poseidon_Kinds_Module, &
                    ONLY : idp

USE Poseidon_Parameters, &
                    ONLY :  DEGREE,                 &
                            L_LIMIT,                &
                            Domain_Dim,             &
                            Num_CFA_Vars,           &
                            Max_Iterations,         &
                            Convergence_Criteria,   &
                            Method_Flag,            &
                            New_Petsc_Solver_Flag


USE Variables_Derived,  &
                    ONLY :  Prob_Dim,               &
                            ULM_Length

USE Variables_Mesh, &
                    ONLY :  NUM_R_ELEMENTS,         &
                            R_COARSEN_FACTOR,       &
                            T_COARSEN_FACTOR,       &
                            P_COARSEN_FACTOR


USE Variables_IO, &
                    ONLY :  Write_Results_R_Samps,      &
                            Write_Results_T_Samps,      &
                            Write_Results_P_Samps,      &
                            Write_Flags,                &
                            Report_Flags,               &
                            Iter_Report_Num_Samples
                            


USE Variables_MPI, &
                    ONLY :  nProcs_Poseidon,        &
                            Num_Block_Phi_Columns,  &
                            Num_Block_Theta_Rows,   &
                            Num_Shells,             &
                            Num_SubShells,          &
                            Num_SubShells_Per_Shell,&
                            Num_Blocks,             &
                            Num_Blocks_Per_Shell,   &
                            Num_R_Elems_Per_Block,  &
                            Num_T_Elems_Per_Block,  &
                            Num_P_Elems_Per_Block,  &
                            Num_R_Elems_Per_Shell,  &
                            Num_R_Elems_Per_SubShell
                            

USE Variables_Quadrature, &
                    ONLY :  Num_R_Quad_Points,      &
                            Num_T_Quad_Points,      &
                            Num_P_Quad_Points,      &
                            Num_Quad_DOF

IMPLICIT NONE

CONTAINS



 !+101+############################################################################!
!                                                                                   !
!                   UNPACK_PARAMETERS                                               !
!                                                                                   !
 !#################################################################################!
SUBROUTINE UNPACK_POSEIDON_PARAMETERS()

INTEGER                                         ::  POSEIDON_read

INTEGER                                         ::  istat

INTEGER                                         ::  NUM_INT_PARAMS
INTEGER                                         ::  NUM_REAL_PARAMS

INTEGER,          DIMENSION(:), ALLOCATABLE     :: INT_PARAMS
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE     :: REAL_PARAMS



NUM_INT_PARAMS  = 34
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
    Report_Flags(5)        = INT_PARAMS(21)    ! Default = 0, Off
END IF
IF ( INT_PARAMS(22) .NE. -1 ) THEN
    Report_Flags(3)           = INT_PARAMS(22)    ! Default = 0, Off
END IF
IF ( INT_PARAMS(23) .NE. -1 ) THEN
    ITER_REPORT_NUM_SAMPLES     = INT_PARAMS(23)    ! Default = 20
END IF
IF ( INT_PARAMS(24) .NE. -1 ) THEN
    Write_Flags(5)          = INT_PARAMS(24)    ! Default = 0, Off
END IF


NEW_PETSC_SOLVER_FLAG       = INT_PARAMS(25)

IF ( INT_PARAMS(26) .NE. -1 ) THEN
    Write_Flags(1)          = INT_PARAMS(26)
END IF
IF ( INT_PARAMS(27) .NE. -1 ) THEN
    Write_Flags(2)          = INT_PARAMS(27)
END IF
IF ( INT_PARAMS(33) .NE. -1 ) THEN
    Write_Flags(3)          = INT_PARAMS(33)
END IF



WRITE_RESULTS_R_SAMPS           = INT_PARAMS(28)
WRITE_RESULTS_T_SAMPS           = INT_PARAMS(29)
WRITE_RESULTS_P_SAMPS           = INT_PARAMS(30)

IF ( INT_PARAMS(31) .NE. -1 ) THEN
    Report_Flags(4)             = INT_PARAMS(31)    ! Default = 0, Off
END IF

IF ( INT_PARAMS(32) .NE. -1 ) THEN
    Write_Flags(4)              = INT_PARAMS(32)    ! Deafult = 0, Off
END IF

CONVERGENCE_CRITERIA            = REAL_PARAMS(1)

IF ( INT_PARAMS(34) .NE. -1 ) THEN
    Method_Flag                 = INT_PARAMS(34)
END IF

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

! DOMAIN_DIM
IF ( Param_type == 'DIM' ) THEN
    READ (line, 111) INT_PARAMS(1)
    CYCLE
END IF

! DEGREE
IF ( Param_type == 'DEGREE' ) THEN
    READ (line, 111) INT_PARAMS(2)
    CYCLE
END IF

! L_LIMIT
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






! Number of RE per Shell
IF ( Param_type == 'NREPS ' ) THEN
    READ (line, 111) INT_PARAMS(5)
    CYCLE
END IF

! Number of RE per Subshell
IF ( Param_type == 'NREPSS' ) THEN
    READ (line, 111) INT_PARAMS(6)
    CYCLE
END IF

! Number of TE per Block
IF ( Param_type == 'NTEPB ' ) THEN
    READ (line, 111) INT_PARAMS(12)
    CYCLE
END IF

! Number of PE per Block
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

!   Maximum Iterations
IF ( Param_type == 'MI    ' ) THEN
    READ (line, 111) INT_PARAMS(20)
    CYCLE
END IF



! Write Timetable Report Flag
IF ( Param_type == 'WRTTT' ) THEN
    READ (line, 111) INT_PARAMS(21)
    CYCLE
END IF

! Write Iteration Report Flage
IF ( Param_type == 'WRTIR' ) THEN
    READ (line, 111) INT_PARAMS(22)
    CYCLE
END IF

! Number of Samples in each Iteration Report
IF ( Param_type == 'IRNS ' ) THEN
    READ (line, 111) INT_PARAMS(23)
    CYCLE
END IF

! Write Results to File Flag
IF ( Param_type == 'WRTRS' ) THEN
    READ (line, 111) INT_PARAMS(24)
    CYCLE
END IF


! Convergence Critera
IF ( Param_type == 'CC    ' ) THEN
    READ (line, 131) REAL_PARAMS(1)
    CYCLE
END IF


IF ( Param_type == 'NPS   ' ) THEN
    READ (line, 111) INT_PARAMS(25)
    CYCLE
END IF

! Write Jacobian to file
IF ( Param_type == 'OMF   ' ) THEN
    READ (line, 111) INT_PARAMS(26)
    CYCLE
END IF

! Write RHS Vector to file
IF ( Param_type == 'ORF   ' ) THEN
    READ (line, 111) INT_PARAMS(27)
    CYCLE
END IF

! Write Update Vector to file
IF ( Param_type == 'OUF   ' ) THEN
    READ (line, 111) INT_PARAMS(33)
    CYCLE
END IF



! WRITE_RESULTS_R_SAMPS
IF ( Param_type == 'RSMPS ' ) THEN
    READ (line, 111) INT_PARAMS(28)
    CYCLE
END IF

! WRITE_RESULTS_T_SAMPS
IF ( Param_type == 'TSMPS ' ) THEN
    READ (line, 111) INT_PARAMS(29)
    CYCLE
END IF

! WRITE_RESULTS_P_SAMPS
IF ( Param_type == 'PSMPS ' ) THEN
    READ (line, 111) INT_PARAMS(30)
    CYCLE
END IF

! Write Sources
IF ( Param_type == 'WRTSC ' ) THEN
    READ (line, 111) INT_PARAMS(32)
    CYCLE
END IF

! 
IF ( Param_type == 'OSTF ' ) THEN
    READ (line, 111) INT_PARAMS(31)
    CYCLE
END IF

! Method_Flag
IF ( Param_type == 'STF  ' ) THEN
    READ (line, 111) INT_PARAMS(34)
    CYCLE
END IF



END DO READ

 5000 RETURN
END SUBROUTINE READ_POSEIDON_PARAMETERS
















END MODULE Poseidon_Parameter_Read_Module
