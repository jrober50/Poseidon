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
                    Num_Off_Diagonals,      &
                    LM_Length

USE Variables_FEM_Module, &
            ONLY :  FEM_Node_xlocs

USE Poseidon_File_Routines_Module, &
            ONLY :  Open_New_File,                  &
                    Open_Existing_File_Rewind

USE Functions_Quadrature, &
            ONLY :  Initialize_LGL_Quadrature_Locations

USE Maps_Fixed_Point, &
            ONLY :  FP_Beta_Array_Map
            
USE Variables_Matrices, &
            ONLY :  dMB_Matrix_Banded,      &
                    iMB_Bandwidth
                    
USE Variables_IO, &
            ONLY :  File_Suffix

IMPLICIT NONE


CHARACTER(LEN = 20), PARAMETER    :: Filename_Format_A = "(A,A)"
CHARACTER(LEN = 20), PARAMETER    :: Filename_Format_B = "(A,A,I2.2,A,I2.2,A)"

CONTAINS

 !+101+############################################################!
!                                                                   !
!          Output_Laplace_Matrix                                    !
!                                                                   !
 !#################################################################!
SUBROUTINE Output_Laplace_Matrix( Matrix )

REAL(idp), DIMENSION(0:2*NUM_OFF_DIAGONALS, 0:SUBSHELL_PROB_DIM-1), INTENT(IN) :: Matrix

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

END SUBROUTINE Output_Laplace_Matrix


 !+102+############################################################!
!                                                                   !
!          Output_Jacobian_Matrix                                   !
!                                                                   !
 !#################################################################!
SUBROUTINE Output_Jacobian_Matrix( Matrix )

REAL(idp), DIMENSION(0:ELEM_PROB_DIM_SQR-1 ,0:Num_R_Elements-1), INTENT(IN) :: Matrix

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

END SUBROUTINE Output_Jacobian_Matrix




 !+201+############################################################!
!                                                                   !
!          Output_RHS_Vector                                        !
!                                                                   !
 !#################################################################!
SUBROUTINE Output_RHS_Vector( Vector )

REAL(idp), DIMENSION(0:Block_Prob_Dim-1), INTENT(IN) ::  Vector

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


END SUBROUTINE Output_RHS_Vector






 !+202+############################################################!
!                                                                   !
!          Output_RHS_Vector_Parts                                  !
!                                                                   !
 !#################################################################!
SUBROUTINE Output_RHS_Vector_Parts(Laplace, Source)

REAL(KIND = idp), DIMENSION(0:SUBSHELL_PROB_DIM-1), INTENT(IN)       :: Laplace
REAL(KIND = idp), DIMENSION(0:SUBSHELL_PROB_DIM-1), INTENT(IN)       :: Source


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




END SUBROUTINE Output_RHS_Vector_Parts






 !+203+############################################################!
!                                                                   !
!          Output_Update_Vector                                     !
!                                                                   !
 !#################################################################!
SUBROUTINE Output_Update_Vector( Vector )

REAL(idp), DIMENSION(0:Block_Prob_Dim-1), INTENT(IN) ::  Vector

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

END SUBROUTINE Output_Update_Vector


 !+204+############################################################!
!                                                                   !
!          Output_Coefficient_Vector_Matlab                         !
!                                                                   !
 !#################################################################!
SUBROUTINE Output_Coefficient_Vector_Matlab( Vector )

REAL(idp), DIMENSION(0:Block_Prob_Dim-1), INTENT(IN) ::  Vector

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


END SUBROUTINE Output_Coefficient_Vector_Matlab



 !+205+############################################################!
!                                                                   !
!          Output_Coefficient_Vector_Fortran                        !
!                                                                   !
 !#################################################################!
SUBROUTINE Output_Coefficient_Vector_Fortran( Vector )

REAL(idp), DIMENSION(0:Block_Prob_Dim-1), INTENT(IN) ::  Vector

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

END SUBROUTINE Output_Coefficient_Vector_Fortran



 !+301+############################################################!
!                                                                   !
!          Read_Coefficient_Vector                                  !
!                                                                   !
 !#################################################################!
SUBROUTINE Read_Coefficient_Vector(Frame_Num, Iter_Num)

INTEGER, INTENT(IN)                                     ::  Frame_Num, Iter_Num

CHARACTER(LEN = 57)                                     ::  FILE_NAME
CHARACTER(LEN = 40)                                     ::  fmt

INTEGER                                                 ::  FILE_ID

REAL(KIND = idp), DIMENSION(0:Block_PROB_DIM-1 )     ::  Test
INTEGER                                                 ::  istat


fmt = '(ES24.16E3,SP,ES24.16E3,"i")'
!fmt = '(F16.10,SP,F16.10,"i")'

WRITE(FILE_NAME,Filename_Format_B)"OUTPUT/Poseidon_Objects/COEFF_VEC_F",Frame_Num,"_I",Iter_Num,".out"
CALL Open_Existing_File_Rewind( FILE_NAME, FILE_ID, istat )

READ(FILE_ID,*)Test


CLOSE(FILE_ID)

END SUBROUTINE Read_Coefficient_Vector




 !+401+############################################################!
!                                                                   !
!          Output_Node_Mesh                                         !
!                                                                   !
 !#################################################################!
SUBROUTINE Output_Node_Mesh()

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

END SUBROUTINE Output_Node_Mesh













 !+501+############################################################!
!                                                                   !
!          Output_Matrix_TypeA                                      !
!                                                                   !
 !#################################################################!
SUBROUTINE Output_Matrix_TypeA()

END SUBROUTINE Output_Matrix_TypeA




 !+502+############################################################!
!                                                                   !
!          Output_Matrix_TypeB                                      !
!                                                                   !
 !#################################################################!
SUBROUTINE Output_Matrix_TypeB()


CHARACTER(LEN = 500),   DIMENSION(2)    ::  Filenames
INTEGER,                DIMENSION(2)    ::  FileIDs

INTEGER                                 ::  re
INTEGER                                 ::  d,  ui, lm
INTEGER                                 ::  dp, uj, lpmp
INTEGER                                 ::  row, col
INTEGER                                 ::  row_bnd
REAL(idp)                               ::  val

!   Dimension and Location Files
WRITE(Filenames(1),'(A,A,A,A)') Poseidon_LinSys_Dir,"Matrix_TypeB_Dims_",TRIM(File_Suffix),".out"
WRITE(Filenames(2),'(A,A,A,A)') Poseidon_LinSys_Dir,"Matrix_TypeB_",TRIM(File_Suffix),".out"

DO ui = 1,2
    CALL OPEN_NEW_FILE( Filenames(ui), FileIDs(ui),250)
END DO

WRITE(FileIDs(1),'(I4.4,1X,I1.1,1XI3.3)')Num_R_Elements, Degree, LM_Length

DO re = 0,Num_R_Elements-1
DO d  = 0,Degree
DO ui = 1,3
DO lm = 1,LM_Length
DO dp = 0,Degree
DO uj = 1,3
DO lpmp = 1,LM_Length
    row_bnd = iMB_Bandwidth + FP_Beta_Array_Map(re,dp,ui,lpmp)
    row = FP_Beta_Array_Map(re,dp,ui,lpmp)
    col = FP_Beta_Array_Map(re,d,uj,lm)
    val = dMB_Matrix_Banded(Row_bnd-Col, Col)
    WRITE(FileIDs(2),'(I5.5,1X,I5.5,1X,ES22.15)') row, col, val
END DO ! lpmp
END DO ! uj
END DO ! dp
END DO ! lm
END DO ! ui
END DO ! d
END DO ! re

END SUBROUTINE Output_Matrix_TypeB




END MODULE IO_Linear_System
