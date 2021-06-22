   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE IO_FP_Linear_System                                                          !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
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
USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Poseidon_Parameters, &
            ONLY :  DEGREE,                     &
                    L_LIMIT

USE Poseidon_IO_Parameters, &
            ONLY :  Poseidon_LinSys_Dir

    
USE Variables_Quadrature, &
            ONLY :  Num_R_Quad_Points,          &
                    Num_T_Quad_Points,          &
                    Num_P_Quad_Points,          &
                    Num_TP_Quad_Points,         &
                    Int_R_Locations,            &
                    Int_T_Locations,            &
                    Int_P_Locations,            &
                    Int_R_Weights,              &
                    Int_T_Weights,              &
                    Int_P_Weights


USE Variables_Mesh, &
            ONLY :  Num_R_Elements,             &
                    Num_T_Elements,             &
                    Num_P_Elements,             &
                    rlocs,                      &
                    tlocs,                      &
                    plocs

USE Variables_Derived, &
            ONLY :  Num_R_Nodes,                &
                    ULM_Length,                 &
                    LM_Length,                  &
                    Var_Dim,                    &
                    Beta_Prob_Dim

USE Variables_Tables, &
            ONLY :  LPT_LPT,                    &
                    Ylm_Values,             &
                    Ylm_dt_Values,          &
                    Ylm_dp_Values,          &
                    Ylm_CC_Values,          &
                    Ylm_CC_dt_Values,       &
                    Ylm_CC_dp_Values

USE Variables_IO, &
            ONLY :  File_Suffix

USE Functions_Quadrature, &
            ONLY :  Initialize_LG_Quadrature

USE Variables_FP, &
            ONLY :  Matrix_Format,              &
                    Laplace_Matrix_Beta,        &
                    CFA_EQ_Flags,               &
                    CFA_EQ_Map,                 &
                    CFA_Mat_Map

USE IO_File_Routines_Module, &
            ONLY :  Open_New_File


IMPLICIT NONE

CONTAINS









!+102+###########################################################################!
!                                                                                !
!                   OUTPUT_JACOBIAN_MATRIX                                       !
!                                                                                !
!################################################################################!
SUBROUTINE OUTPUT_LAPLACE_Beta( Matrix, row_in, col_in, flag )

INTEGER,                        INTENT(IN)              ::  row_in, col_in
COMPLEX(idp), DIMENSION(1:Row_In,1:Col_in), INTENT(IN)                ::  Matrix


CHARACTER(LEN = 1 ), INTENT(IN),OPTIONAL                ::  flag

CHARACTER(LEN = 300)                                     ::  FILE_NAME
CHARACTER(LEN = 300)                                     ::  FILE_NAMEb
CHARACTER(LEN = 40)                                     ::  fmt


INTEGER                                                 ::  rows, cols

INTEGER                                                 ::  FILE_ID
INTEGEr                                                 ::  i,j

100 FORMAT (A,A,A,A,A,A)
101 FORMAT (A,A,A,A)

fmt = '(ES24.16E3,SP,ES24.16E3,"i")'

cols = size( Laplace_Matrix_Beta, 1 )
rows = size( Laplace_Matrix_Beta, 2 )


IF ( present(flag) ) THEN
    WRITE(FILE_NAMEb,100) Poseidon_LinSys_Dir,"Beta_Laplace_Dim_",trim(File_Suffix),"_",flag,".out"
ELSE
    WRITE(FILE_NAMEb,101) Poseidon_LinSys_Dir,"Beta_Laplace_Dim_",trim(File_Suffix),".out"
END IF
CALL OPEN_NEW_FILE( trim(FILE_NAMEb), FILE_ID, 300 )
WRITE(FILE_ID,*) Cols, Rows, NUM_R_Elements, DEGREE, L_LIMIT
CLOSE(FILE_ID)


IF ( present(flag) ) THEN
    WRITE(FILE_NAME,100) Poseidon_LinSys_Dir,"Beta_Laplace_",trim(File_Suffix),"_",flag,".out"
ELSE
    WRITE(FILE_NAME,101) Poseidon_LinSys_Dir,"Beta_Laplace_",trim(File_Suffix),".out"
END IF

PRINT*,"Writing Shift Matrix to file : ",File_Name


CALL OPEN_NEW_FILE( trim(FILE_NAME), FILE_ID, 300 )
DO i = 1,rows
    DO j = 1, cols
        WRITE(FILE_ID,fmt) Matrix(j,i)
    END DO
END DO
CLOSE( FILE_ID )



END SUBROUTINE OUTPUT_LAPLACE_Beta



!+101+###########################################################################!
!                                                                                !
!                   OUTPUT_JACOBIAN_MATRIX                                       !
!                                                                                !
!################################################################################!
SUBROUTINE OUTPUT_LAPLACE( Matrix, row_in, col_in, flag )

INTEGER,                        INTENT(IN)              ::  row_in, col_in
COMPLEX(idp), DIMENSION(0:Row_In-1,0:Col_in-1), INTENT(IN)                ::  Matrix


CHARACTER(LEN = 1 ), INTENT(IN)                          ::  flag

CHARACTER(LEN = 300)                                     ::  FILE_NAME
CHARACTER(LEN = 300)                                     ::  FILE_NAMEb
CHARACTER(LEN = 40)                                     ::  fmt


INTEGER                                                 ::  rows, cols

INTEGER                                                 ::  FILE_ID
INTEGEr                                                 ::  i,j

100 FORMAT (A,A,A,A,A,A)

fmt = '(ES24.16E3,SP,ES24.16E3,"i")'

cols = size( Matrix, 1 )
rows = size( Matrix, 2 )



WRITE(FILE_NAMEb,100) Poseidon_LinSys_Dir,"Laplace_Dim_",trim(File_Suffix),"_",flag,".out"
CALL OPEN_NEW_FILE( trim(FILE_NAMEb), FILE_ID, 300 )
WRITE(FILE_ID,*) Cols, Rows, NUM_R_Elements, DEGREE, L_LIMIT
CLOSE(FILE_ID)



WRITE(FILE_NAME,100) Poseidon_LinSys_Dir,"Laplace_",trim(File_Suffix),"_",flag,".out"
CALL OPEN_NEW_FILE( trim(FILE_NAME), FILE_ID, 300 )
DO i = 0,row_in-1
    DO j = 0, col_in-1
        WRITE(FILE_ID,fmt) Matrix(j,i)
    END DO
END DO
CLOSE( FILE_ID )



END SUBROUTINE OUTPUT_LAPLACE









!+202+###########################################################################!
!                                                                                !
!                   OUTPUT_JACOBIAN_MATRIX                                       !
!                                                                                !
!################################################################################!
SUBROUTINE Output_Source_Beta( Source, row_in, flag )

INTEGER,                        INTENT(IN)                  ::  row_in
COMPLEX(idp), DIMENSION(1:Row_In), INTENT(IN)               ::  Source


CHARACTER(LEN = 1 ), INTENT(IN), OPTIONAL                   ::  flag

CHARACTER(LEN = 300)                                        ::  FILE_NAME
CHARACTER(LEN = 300)                                        ::  FILE_NAMEb
CHARACTER(LEN = 40)                                         ::  fmt


INTEGER                                                     ::  rows, cols

INTEGER                                                     ::  FILE_ID
INTEGEr                                                     ::  i

100 FORMAT (A,A,A,A,A,A)
101 FORMAT (A,A,A,A)

fmt = '(ES24.16E3,SP,ES24.16E3,"i")'


rows = size( Source, 1 )

IF ( present(flag) ) THEN
    WRITE(FILE_NAMEb,100) Poseidon_LinSys_Dir,"Beta_Source_Dim_",trim(File_Suffix),"_",flag,".out"
ELSE
    WRITE(FILE_NAMEb,101) Poseidon_LinSys_Dir,"Beta_Source_Dim_",trim(File_Suffix),".out"
END IF


CALL OPEN_NEW_FILE( trim(FILE_NAMEb), FILE_ID, 300 )
WRITE(FILE_ID,*) Cols, Rows, NUM_R_Elements, DEGREE, L_LIMIT
CLOSE(FILE_ID)



IF ( present(flag) ) THEN
    WRITE(FILE_NAME,100) Poseidon_LinSys_Dir,"Beta_Source_",trim(File_Suffix),"_",flag,".out"
ELSE
    WRITE(FILE_NAME,100) Poseidon_LinSys_Dir,"Beta_Source_",trim(File_Suffix),".out"
END IF

CALL OPEN_NEW_FILE( trim(FILE_NAME), FILE_ID, 300 )
DO i = 1,rows
    WRITE(FILE_ID,fmt) Source(i)
END DO
CLOSE( FILE_ID )



END SUBROUTINE Output_Source_Beta









!+202+###########################################################################!
!                                                                                !
!                   OUTPUT_JACOBIAN_MATRIX                                       !
!                                                                                !
!################################################################################!
SUBROUTINE Output_Coeffs_Beta( Source, row_in, flag )

INTEGER,                        INTENT(IN)                  ::  row_in
COMPLEX(idp), DIMENSION(1:Row_In), INTENT(IN)               ::  Source


CHARACTER(LEN = 1 ), INTENT(IN), OPTIONAL                   ::  flag

CHARACTER(LEN = 300)                                        ::  FILE_NAME
CHARACTER(LEN = 300)                                        ::  FILE_NAMEb
CHARACTER(LEN = 40)                                         ::  fmt


INTEGER                                                     ::  rows, cols

INTEGER                                                     ::  FILE_ID
INTEGEr                                                     ::  i

100 FORMAT (A,A,A,A,A,A)
101 FORMAT (A,A,A,A)

fmt = '(ES24.16E3,SP,ES24.16E3,"i")'


rows = size( Source, 1 )

IF ( present(flag) ) THEN
    WRITE(FILE_NAMEb,100) Poseidon_LinSys_Dir,"Beta_Coeffs_Dim_",trim(File_Suffix),"_",flag,".out"
ELSE
    WRITE(FILE_NAMEb,101) Poseidon_LinSys_Dir,"Beta_Coeffs_Dim_",trim(File_Suffix),".out"
END IF


CALL OPEN_NEW_FILE( trim(FILE_NAMEb), FILE_ID, 300 )
WRITE(FILE_ID,*) Cols, Rows, NUM_R_Elements, DEGREE, L_LIMIT
CLOSE(FILE_ID)



IF ( present(flag) ) THEN
    WRITE(FILE_NAME,100) Poseidon_LinSys_Dir,"Beta_Coeffs_",trim(File_Suffix),"_",flag,".out"
ELSE
    WRITE(FILE_NAME,100) Poseidon_LinSys_Dir,"Beta_Coeffs_",trim(File_Suffix),".out"
END IF


PRINT*,"Writing Shift Vector Coefficients to file : ",File_Name


CALL OPEN_NEW_FILE( trim(FILE_NAME), FILE_ID, 300 )
DO i = 1,rows
    WRITE(FILE_ID,fmt) Source(i)
END DO
CLOSE( FILE_ID )



END SUBROUTINE Output_Coeffs_Beta






END MODULE IO_FP_Linear_System

