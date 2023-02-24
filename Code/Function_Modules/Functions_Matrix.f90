   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Functions_Matrix                                                             !##!
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




CONTAINS



!+401+############################################################!
!
!     MVMULT_CCS: Multiply a vector by a matrix in CCS format
!
!#################################################################!
FUNCTION MVMULT_CCS(N, NNZ, ELEM_VAL, COL_PTR, ROW_IND, VECT)


INTEGER, INTENT(IN)                                         :: N, NNZ
INTEGER, DIMENSION(0:N), INTENT(IN)                         :: COL_PTR
INTEGER, DIMENSION(0:NNZ - 1),INTENT(IN)                    :: ROW_IND

REAL(KIND = idp), DIMENSION(0:NNZ - 1), INTENT(IN)       :: ELEM_VAL
REAL(KIND = idp), DIMENSION(0:N-1), INTENT(IN)           :: VECT


REAL(KIND = idp), DIMENSION(0:N-1)                          :: MVMULT_CCS

INTEGER                                                     :: i, j

MVMULT_CCS = 0.0_idp

DO i = 0,N-1

    PRINT*,i,COL_PTR(i),COL_PTR(i+1)-1
    DO j = COL_PTR(i),COL_PTR(i+1)-1
        
        MVMULT_CCS(ROW_IND(j)) =  MVMULT_CCS(ROW_IND(j)) + ELEM_VAL(j)*VECT(i)

    END DO


END DO



END FUNCTION MVMULT_CCS




!+402+############################################################!
!
!     MVMULT_FULL: Multiply a vector by a matrix in CCS format
!
!#################################################################!
PURE FUNCTION MVMULT_FULL(A, V, N, M)


INTEGER, INTENT(IN)                                         :: N, M
REAL(KIND = idp), INTENT(IN), DIMENSION(1:M)             :: V
REAL(KIND = idp), INTENT(IN), DIMENSION(1:N,1:M)         :: A


REAL(KIND = idp), DIMENSION(1:N)                            :: MVMULT_FULL


REAL(KIND = idp), DIMENSION(1:N)                            :: Sol

INTEGER                                                     :: i,j



Do i = 1,N

    Sol(i) = 0.D0

    Do j = 1,M

    Sol(i) = Sol(i) + A(i,j)*V(j)

    END DO

END DO

MVMULT_FULL = Sol

END FUNCTION MVMULT_FULL




!+402+############################################################!
!
!     MVMULT_FULL: Multiply a vector by a matrix in CCS format
!
!#################################################################!
SUBROUTINE MVMULT_FULL_SUB(A, V, N, M)


INTEGER, INTENT(IN)                                         :: N, M
REAL(KIND = idp), INTENT(IN), DIMENSION(1:M)             :: V
REAL(KIND = idp), INTENT(IN), DIMENSION(1:N,1:M)         :: A


REAL(KIND = idp), DIMENSION(1:N)                            :: Sol

INTEGER                                                     :: i,j



Do i = 1,N

    Sol(i) = 0.D0

    Do j = 1,M
        PRINT*,i,j,REAL(A(i,j),idp),REAL(V(j),idp),REAL(A(i,j)*V(j),idp),REAL(Sol(i) + A(i,j)*V(j),idp)
        Sol(i) = Sol(i) + A(i,j)*V(j)
    
    END DO

END DO


END SUBROUTINE MVMULT_FULL_SUB







END MODULE Functions_Matrix
