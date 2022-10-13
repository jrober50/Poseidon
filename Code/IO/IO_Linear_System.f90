   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE IO_Linear_System                                                             !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
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
            ONLY : idp

USE Poseidon_Numbers_Module, &
            ONLY : pi


USE Poseidon_Units_Module, &
            ONLY :  Centimeter

USE Poseidon_IO_Parameters, &
            ONLY :  Poseidon_LinSys_Dir,    &
                    Poseidon_Objects_Dir

USE Poseidon_Parameters, &
            ONLY :  DEGREE,                 &
                    L_LIMIT,                &
                    CUR_ITERATION,          &
                    Poseidon_Frame


USE Variables_Mesh, &
            ONLY :  NUM_R_ELEMENTS,         &
                    rlocs

USE Variables_Derived, &
            ONLY :  Block_Prob_Dim,         &
                    SubShell_Prob_Dim,      &
                    Elem_Prob_Dim_Sqr,      &
                    Num_Off_Diagonals

USE Variables_FEM_Module, &
            ONLY :  FEM_Node_xlocs

USE Poseidon_File_Routines_Module, &
            ONLY :  Open_New_File,                  &
                    Open_Existing_File_Rewind

USE Functions_Quadrature, &
            ONLY :  Initialize_LGL_Quadrature_Locations


IMPLICIT NONE


CHARACTER(LEN = 20), PARAMETER    :: Filename_Format_A = "(A,A)"
CHARACTER(LEN = 20), PARAMETER    :: Filename_Format_B = "(A,A,I2.2,A,I2.2,A)"

CONTAINS

!+403+###########################################################################!
!                                                                                !
!                   OUTPUT_JACOBIAN_MATRIX                                       !
!                                                                                !
!################################################################################!
SUBROUTINE OUTPUT_LAPLACE_MATRIX( Matrix )

COMPLEX(idp), DIMENSION(0:2*NUM_OFF_DIAGONALS, 0:SUBSHELL_PROB_DIM-1), INTENT(IN) :: Matrix

CHARACTER(LEN = 57)                                     ::  FILE_NAME
CHARACTER(LEN = 61)                                     ::  FILE_NAMEb
CHARACTER(LEN = 40)                                     ::  fmt


INTEGER                                                 ::  FILE_ID
INTEGEr                                                 ::  re,e

100 FORMAT (A,A,I2.2,A,I2.2,A)

fmt = '(ES24.16E3,SP,ES24.16E3,"i")'




WRITE(FILE_NAMEb,100) Poseidon_LinSys_Dir,"LAP_MAT_DIM_F",Poseidon_Frame,"_I",CUR_ITERATION,".out"
CALL OPEN_NEW_FILE( FILE_NAMEb, FILE_ID, 300 )
WRITE(FILE_ID,*) Num_R_Elements,NUM_OFF_DIAGONALS, DEGREE, L_LIMIT
CLOSE(FILE_ID)


WRITE(FILE_NAME,Filename_Format_B) Poseidon_LinSys_Dir,"LAP_MAT_F",Poseidon_Frame,"_I",CUR_ITERATION,".out"
CALL OPEN_NEW_FILE( FILE_NAME, FILE_ID, 300 )


PRINT*,"Writing Laplace Matrix to file: ", File_Name


DO re = 0,SUBSHELL_PROB_DIM-1
    DO e = 0,2*NUM_OFF_DIAGONALS
        WRITE(FILE_ID,fmt) Matrix(e,re)
    END DO
END DO

CLOSE( FILE_ID )

END SUBROUTINE OUTPUT_LAPLACE_MATRIX


!+403+###########################################################################!
!                                                                                !
!                   OUTPUT_JACOBIAN_MATRIX                                       !
!                                                                                !
!################################################################################!
SUBROUTINE OUTPUT_JACOBIAN_MATRIX( Matrix )

COMPLEX(idp), DIMENSION(0:ELEM_PROB_DIM_SQR-1 ,0:Num_R_Elements-1), INTENT(IN) :: Matrix

CHARACTER(LEN = 57)                                     ::  FILE_NAME
CHARACTER(LEN = 61)                                     ::  FILE_NAMEb
CHARACTER(LEN = 40)                                     ::  fmt


INTEGER                                                 ::  FILE_ID
INTEGEr                                                 ::  re,e

100 FORMAT (A,A,I2.2,A,I2.2,A)

fmt = '(ES24.16E3,SP,ES24.16E3,"i")'


WRITE(FILE_NAMEb,100) Poseidon_LinSys_Dir,"STF_MAT_DIM_F",Poseidon_Frame,"_I",CUR_ITERATION,".out"
CALL OPEN_NEW_FILE( FILE_NAMEb, FILE_ID, 400 )
WRITE(FILE_ID,*) Num_R_Elements,DEGREE,L_LIMIT
CLOSE(FILE_ID)


WRITE(FILE_NAME,Filename_Format_B) Poseidon_LinSys_Dir,"STF_MAT_F",Poseidon_Frame,"_I",CUR_ITERATION,".out"
CALL OPEN_NEW_FILE( FILE_NAME, FILE_ID, 400 )

PRINT*,"Writing Jacobian Matrix to file: ", File_Name


DO re = 0,Num_R_Elements-1
    DO e = 0,ELEM_PROB_DIM_SQR-1
        WRITE(FILE_ID,fmt) Matrix(e,re)
!        WRITE(*,*) Matrix(e,re)
    END DO
!    WRITE(*,*)" "
!    WRITE(*,*)" "
!    WRITE(*,*)" "
END DO

CLOSE(FILE_ID)

END SUBROUTINE OUTPUT_JACOBIAN_MATRIX




!+404+###########################################################################!
!                                                                                !
!                   OUTPUT_RHS_VECTOR                                       !
!                                                                                !
!################################################################################!
SUBROUTINE OUTPUT_RHS_VECTOR( Vector )

COMPLEX(idp), DIMENSION(0:Block_Prob_Dim-1), INTENT(IN) ::  Vector

CHARACTER(LEN = 57)                                     ::  FILE_NAME
CHARACTER(LEN = 40)                                     ::  fmt

INTEGER                                                 ::  FILE_ID
INTEGER                                                 ::  i

101 FORMAT (I5.5," ",I2.2," ",I2.2)

fmt = '(ES24.16E3,SP,ES24.16E3,"i")'
!fmt = '(F16.10,SP,F16.10,"i")'


WRITE(FILE_NAME,Filename_Format_B) Poseidon_LinSys_Dir,"RHS_VEC_F",Poseidon_Frame,"_I",CUR_ITERATION,".out"
CALL OPEN_NEW_FILE( FILE_NAME, FILE_ID, 500 )


PRINT*,"Writing RHS Vector to file: ", File_Name

WRITE(FILE_ID, 101)NUM_R_ELEMENTS,DEGREE,L_LIMIT
DO i = 0,Block_PROB_DIM-1
    WRITE(FILE_ID,TRIM(fmt)) Vector(i)
END DO

CLOSE(FILE_ID)


END SUBROUTINE OUTPUT_RHS_VECTOR






!+404+###########################################################################!
!                                                                                !
!                   OUTPUT_RHS_VECTOR_Parts                                      !
!                                                                                !
!################################################################################!
SUBROUTINE OUTPUT_RHS_VECTOR_Parts(Laplace, Source)

COMPLEX(KIND = idp), DIMENSION(0:SUBSHELL_PROB_DIM-1), INTENT(IN)       :: Laplace
COMPLEX(KIND = idp), DIMENSION(0:SUBSHELL_PROB_DIM-1), INTENT(IN)       :: Source


CHARACTER(LEN = 57)                                     ::  FILE_NAME
CHARACTER(LEN = 40)                                     ::  fmt

INTEGER                                                 ::  FILE_ID
INTEGER                                                 ::  i

101 FORMAT (I5.5," ",I2.2," ",I2.2)

fmt = '(ES24.16E3,SP,ES24.16E3,"i")'
!fmt = '(F16.10,SP,F16.10,"i")'

CALL OUTPUT_NODE_MESH()

WRITE(FILE_NAME,Filename_Format_B) Poseidon_LinSys_Dir,"RHS_LAP_F",Poseidon_Frame,"_I",CUR_ITERATION,".out"
CALL OPEN_NEW_FILE( FILE_NAME, FILE_ID, 600 )


WRITE(FILE_ID, 101)NUM_R_ELEMENTS,DEGREE,L_LIMIT
DO i = 0,Block_PROB_DIM-1
    WRITE(FILE_ID,TRIM(fmt)) Laplace(i)
END DO

CLOSE(FILE_ID)


WRITE(FILE_NAME,Filename_Format_B) Poseidon_LinSys_Dir,"RHS_SRC_F",Poseidon_Frame,"_I",CUR_ITERATION,".out"
CALL OPEN_NEW_FILE( FILE_NAME, FILE_ID, 600 )


WRITE(FILE_ID, 101)NUM_R_ELEMENTS,DEGREE,L_LIMIT
DO i = 0,Block_PROB_DIM-1
    WRITE(FILE_ID,TRIM(fmt)) Source(i)
END DO
CLOSE(FILE_ID)




END SUBROUTINE OUTPUT_RHS_VECTOR_Parts






!+404+###########################################################################!
!                                                                                !
!                   OUTPUT_UPDATE_VECTOR                                         !
!                                                                                !
!################################################################################!
SUBROUTINE OUTPUT_UPDATE_VECTOR( Vector )

COMPLEX(idp), DIMENSION(0:Block_Prob_Dim-1), INTENT(IN) ::  Vector

CHARACTER(LEN = 57)                                     ::  FILE_NAME
CHARACTER(LEN = 40)                                     ::  fmt

INTEGER                                                 ::  FILE_ID
INTEGER                                                 ::  i

101 FORMAT (I5.5," ",I2.2," ",I2.2)

fmt = '(ES24.16E3,SP,ES24.16E3,"i")'
!fmt = '(F16.10,SP,F16.10,"i")'


WRITE(FILE_NAME,Filename_Format_B) Poseidon_LinSys_Dir,"UPD_VEC_F",Poseidon_Frame,"_I",CUR_ITERATION,".out"
CALL OPEN_NEW_FILE( FILE_NAME, FILE_ID, 700 )


WRITE(FILE_ID, 101)NUM_R_ELEMENTS,DEGREE,L_LIMIT
DO i = 0,Block_PROB_DIM-1
    WRITE(FILE_ID,TRIM(fmt)) Vector(i)
END DO

CLOSE(FILE_ID)

END SUBROUTINE OUTPUT_UPDATE_VECTOR


!+404+###########################################################################!
!                                                                                !
!                   OUTPUT_COEFFICIENT_VECTOR_MATLAB                             !
!                                                                                !
!################################################################################!
SUBROUTINE OUTPUT_COEFFICIENT_VECTOR_MATLAB( Vector )

COMPLEX(idp), DIMENSION(0:Block_Prob_Dim-1), INTENT(IN) ::  Vector

CHARACTER(LEN = 57)                                     ::  FILE_NAME
CHARACTER(LEN = 40)                                     ::  fmt

INTEGER                                                 ::  FILE_ID
INTEGER                                                 ::  i

101 FORMAT (I5.5," ",I2.2," ",I2.2)

fmt = '(ES24.16E3,SP,ES24.16E3,"i")'
!fmt = '(F16.10,SP,F16.10,"i")'


WRITE(FILE_NAME,Filename_Format_B) Poseidon_Objects_Dir,"COEFF_VEC_MATLAB_F",Poseidon_Frame,"_I",CUR_ITERATION,".out"
CALL OPEN_NEW_FILE( FILE_NAME, FILE_ID, 800 )


WRITE(FILE_ID, 101)NUM_R_ELEMENTS,DEGREE,L_LIMIT
DO i = 0,Block_PROB_DIM-1
    WRITE(FILE_ID,TRIM(fmt)) Vector(i)
END DO

CLOSE(FILE_ID)


END SUBROUTINE OUTPUT_COEFFICIENT_VECTOR_MATLAB

!+404+###########################################################################!
!                                                                                !
!                   OUTPUT_COEFFICIENT_VECTOR_FORTRAN                            !
!                                                                                !
!################################################################################!
SUBROUTINE OUTPUT_COEFFICIENT_VECTOR_FORTRAN( Vector )

COMPLEX(idp), DIMENSION(0:Block_Prob_Dim-1), INTENT(IN) ::  Vector

CHARACTER(LEN = 57)                                     ::  FILE_NAME
CHARACTER(LEN = 40)                                     ::  fmt

INTEGER                                                 ::  FILE_ID

fmt = '(ES24.16E3,SP,ES24.16E3,"i")'
!fmt = '(F16.10,SP,F16.10,"i")'

IF ( .FALSE. ) THEN

    WRITE(FILE_NAME,Filename_Format_B) Poseidon_Objects_Dir,"COEFF_VEC_F",Poseidon_Frame,"_I",CUR_ITERATION,".out"
    CALL OPEN_NEW_FILE( FILE_NAME, FILE_ID, 800 )

    WRITE(FILE_ID,*) Vector


    CLOSE(FILE_ID)

END IF

END SUBROUTINE OUTPUT_COEFFICIENT_VECTOR_FORTRAN



!+404+###########################################################################!
!                                                                                !
!                   READ_COEFFICIENT_VECTOR                                      !
!                                                                                !
!################################################################################!
SUBROUTINE READ_COEFFICIENT_VECTOR(Frame_Num, Iter_Num)

INTEGER, INTENT(IN)                                     ::  Frame_Num, Iter_Num

CHARACTER(LEN = 57)                                     ::  FILE_NAME
CHARACTER(LEN = 40)                                     ::  fmt

INTEGER                                                 ::  FILE_ID

COMPLEX(KIND = idp), DIMENSION(0:Block_PROB_DIM-1 )     ::  Test
INTEGER                                                 ::  istat


fmt = '(ES24.16E3,SP,ES24.16E3,"i")'
!fmt = '(F16.10,SP,F16.10,"i")'

WRITE(FILE_NAME,Filename_Format_B)"OUTPUT/Poseidon_Objects/COEFF_VEC_F",Frame_Num,"_I",Iter_Num,".out"
CALL Open_Existing_File_Rewind( FILE_NAME, FILE_ID, istat )

READ(FILE_ID,*)Test


CLOSE(FILE_ID)

END SUBROUTINE READ_COEFFICIENT_VECTOR









!+404+###########################################################################!
!                                                                                !
!                   OUTPUT_NODE_MESH                                             !
!                                                                                !
!################################################################################!
SUBROUTINE OUTPUT_NODE_MESH()

CHARACTER(LEN = 57)                                     ::  FILE_NAME



INTEGER                                                 ::  FILE_ID
INTEGER                                                 ::  re, rd

REAL(KIND = idp), DIMENSION(0:DEGREE)                   ::  CUR_R_LOCS
REAL(KIND = idp)                                        ::  deltar_overtwo



WRITE(FILE_NAME,'(A,A)') Poseidon_LinSys_Dir,"Nodal_Mesh.out"
CALL OPEN_NEW_FILE( FILE_NAME, FILE_ID, 300 )


deltar_overtwo = (rlocs(1) - rlocs(0))/2.0_idp
CUR_R_LOCS(:) = deltar_overtwo * (FEM_Node_xlocs(:)+1.0_idp) + rlocs(0)


WRITE(FILE_ID, '(ES24.16E3)') CUR_R_LOCS(0)/centimeter
DO re = 0,NUM_R_ELEMENTS-1

    deltar_overtwo = (rlocs(re + 1) - rlocs(re))/2.0_idp
    CUR_R_LOCS(:) = deltar_overtwo * (FEM_Node_xlocs(:)+1.0_idp) + rlocs(re)

    DO rd = 1,DEGREE

        WRITE(FILE_ID,'(ES24.16E3)' ) CUR_R_LOCS(rd)/centimeter

    END DO
END DO


CLOSE(FILE_ID)

END SUBROUTINE OUTPUT_NODE_MESH












END MODULE IO_Linear_System
