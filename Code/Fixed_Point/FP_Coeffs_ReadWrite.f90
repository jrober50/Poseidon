  !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE FP_Coeffs_ReadWrite                                                           !##!
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

USE Poseidon_IO_Parameters, &
            ONLY :  Poseidon_Coeffs_Dir

USE Variables_Derived,   &
            ONLY :  Prob_Dim,               &
                    LM_Length,              &
                    Num_R_Nodes

USE Variables_FP,   &
            ONLY :  FP_Coeff_Vector,        &
                    FP_Coeff_Vector_Beta

USE Variables_IO, &
            ONLY :  File_Suffix

USE IO_File_Routines_Module, &
            ONLY :  Open_New_File,          &
                    Open_Existing_File

USE Poseidon_Parameters, &
            ONLY :  DEGREE,                 &
                    L_LIMIT

USE Variables_Mesh, &
            ONLY :  Num_R_Elements

USE FP_Functions_Mapping, &
            ONLY :  FP_Beta_Array_Map,          &
                    FP_FEM_Node_Map,            &
                    FP_LM_Map

IMPLICIT NONE

CONTAINS






!+101+###########################################################################!
!                                                                                !
!           Output_FP_Coeffs                                              !
!                                                                                !
!################################################################################!
SUBROUTINE Output_FP_Coeffs()



CHARACTER(LEN = 100)                                        ::  FileName
INTEGER                                                     ::  FileID
INTEGER                                                     ::  ui,lq,rq

116 FORMAT (A,A,A,A)

WRITE(FileName,116) Poseidon_Coeffs_Dir,"Coeffs_",TRIM(File_Suffix),".out"

CALL OPEN_NEW_FILE( FileName, FileID, 342 )



DO rq = 1,Num_R_Nodes
DO ui = 1,5
DO lq = 1,LM_Length
    WRITE(FileID,'(2ES24.17)') REAL(FP_Coeff_Vector_A(rq,lq,ui),idp),AIMAG(FP_Coeff_Vector_A(rq,lq,ui))
END DO
END DO
END DO


CLOSE( Unit = FileID )



END SUBROUTINE Output_FP_Coeffs








!+201+###########################################################################!
!                                                                                !
!           Output_FP_Coeffs                                              !
!                                                                                !
!################################################################################!
SUBROUTINE ReadIn_FP_Coeffs( RE_In, D_In, L_In, Tail_In )



INTEGER, INTENT(IN)                                         ::  RE_In
INTEGER, INTENT(IN)                                         ::  D_In
INTEGER, INTENT(IN)                                         ::  L_In
CHARACTER(LEN=1), INTENT(IN), Optional                      ::  Tail_In

REAL(idp), DIMENSION(:,:), ALLOCATABLE                      ::  tmp
INTEGER                                                     ::  Dim
INTEGER                                                     ::  LDim
INTEGER                                                     ::  i, Here, There

CHARACTER(LEN = 100)                                        ::  FileName
INTEGER                                                     ::  FileID
CHARACTER(LEN=20)                                           ::  Suffix
INTEGER                                                     ::  istat

INTEGER                                                     ::  RE, d, ui, lm

116 FORMAT (A,A,A,A)



IF ( L_In > L_Limit ) THEN

    PRINT*,"Poseidon can not downscale the spherical harmonic expansion."
    PRINT*,"Poseidon is stopping. "
    STOP

ELSE IF ( D_In .NE. Degree ) Then

    PRINT*,"The FEM degree of the requested coefficient file does",  &
            " not match that of the current expansion"
    PRINT*,"Poseidon is stopping. "
    STOP

ELSE IF ( RE_In .NE. Num_R_Elements ) THEN

    PRINT*,"The number of radial elements of the requested coefficient file does",  &
            " not match that of the current expansion"
    PRINT*,"Poseidon is stopping. "
    STOP

ELSE

    LDim =(L_In+1)*(L_In+1)
    Dim = 5*(L_In + 1)*(L_In + 1)*(D_In*RE_In + 1)

    ALLOCATE( TMP(Dim,2) )


    WRITE(Suffix,'(A,I4.4,A,I2.2,A,I2.2)')"RE",RE_In,"_D",D_In,"_L",L_IN
    IF ( PRESENT(Tail_In) ) THEN
        WRITE(Suffix,'(A,A,A)') TRIM(Suffix),"_",Tail_In

    END IF
    WRITE(FileName,116) Poseidon_Coeffs_Dir,"Coeffs_",TRIM(Suffix),".out"

    CALL OPEN_EXISTING_FILE( FileName, FileID, istat )
    PRINT*,Filename


    DO i = 1,Dim
       READ(FileID,'(2ES24.17)') TMP(i,1), TMP(i,2)
!        PRINT*,Tmp(i,1),TMp(i,2)
    END DO

    CLOSE(UNIT=FileID)



    FP_Coeff_Vector = 0.0_idp

    DO RE = 0,RE_In-1
    DO d  = 0,D_In
    DO ui = 1,5
    DO lm = 1,LDim

        Here = (RE*D_In + d)*5*LDim     &
             +  (ui - 1)*LDim           &
             +  lm

        FP_Coeff_Vector(RE*D_In+d+1,lm,ui) = CMPLX(TMP(Here,1), TMP(Here,2),KIND = idp)

    END DO ! lm
    END DO ! ui
    END DO ! d
    END DO ! RE




    DO ui = 1,3
    DO re = 0,Num_R_Elements -1
    DO d = 0,Degree
    DO lm = 1,LM_Length

        Here  = FP_Beta_Array_Map(re,d,ui,lm)
        There = FP_FEM_Node_Map(re,d)

        FP_Coeff_Vector_Beta(Here) = FP_Coeff_Vector(There,lm,ui+2)

    END DO  ! l
    END DO ! d
    END DO ! re
    END DO ! ui




END IF


END SUBROUTINE ReadIn_FP_Coeffs















END MODULE FP_Coeffs_ReadWrite

