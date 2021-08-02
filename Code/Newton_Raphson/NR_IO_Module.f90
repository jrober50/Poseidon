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
                            Num_CFA_Vars

USE Variables_Mesh, &
                    ONLY :  Num_R_Elements

USE Variables_Derived,  &
                    ONLY :  Prob_Dim,               &
                            ULM_Length

USE Variables_NR, &
                    ONLY :  NR_Coeff_Vector

IMPLICIT NONE

CONTAINS







 !+201+############################################################################!
!                                                                                   !
!                    WRITE_CFA_COEFFICIENTS                                         !
!                                                                                   !
 !#################################################################################!
SUBROUTINE WRITE_CFA_COEFFICIENTS( )


INTEGER                                             ::  re, d, l, m, here
INTEGER                                             ::  Coeffs_Write
INTEGER                                             ::  istat

Coeffs_Write = 13

OPEN(UNIT=Coeffs_Write, FILE='Params/CFA_Coeffs.coefs', ACTION='WRITE', STATUS='REPLACE',IOSTAT=istat)

!DO i = 0,PROB_DIM-1
!    WRITE(Coeffs_Write,'(2ES24.17)') REAL(NR_Coeff_Vector(i)), AIMAG(NR_Coeff_Vector(i))
!END DO

DO re = 0,NUM_R_ELEMENTS-1
    DO d = 0,DEGREE
        DO l = 0,L_LIMIT
            DO m = -l,l
                here = (re*DEGREE+d)*ULM_LENGTH+(l*(l+1)+m)*NUM_CFA_VARS
                WRITE(Coeffs_Write,'(A2,I1,A3,I2)')"L=",l,",M=",l
                WRITE(Coeffs_Write,'(2ES24.17)')NR_Coeff_Vector(here:here+4)
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


INTEGER                                         :: i
INTEGER                                         :: Coeffs_read
INTEGER                                         :: istat
REAL(KIND=idp), DIMENSION(0:1,0:PROB_DIM-1)     :: TMP



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


!NR_Coeff_Vector = 0.0_idp
NR_Coeff_Vector(:) = CMPLX(TMP(0,:), TMP(1,:), KIND =idp)



END SUBROUTINE READ_CFA_COEFFICIENTS












END MODULE Poseidon_Parameter_Read_Module

