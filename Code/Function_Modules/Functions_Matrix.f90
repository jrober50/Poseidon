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

COMPLEX(KIND = idp), DIMENSION(0:NNZ - 1), INTENT(IN)       :: ELEM_VAL
COMPLEX(KIND = idp), DIMENSION(0:N-1), INTENT(IN)           :: VECT


COMPLEX(KIND = idp), DIMENSION(0:N-1)                          :: MVMULT_CCS

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
COMPLEX(KIND = idp), INTENT(IN), DIMENSION(1:M)             :: V
COMPLEX(KIND = idp), INTENT(IN), DIMENSION(1:N,1:M)         :: A


COMPLEX(KIND = idp), DIMENSION(1:N)                            :: MVMULT_FULL


COMPLEX(KIND = idp), DIMENSION(1:N)                            :: Sol

INTEGER                                                     :: i,j



Do i = 1,N

    Sol(i) = 0.D0

    Do j = 1,M

    Sol(i) = Sol(i) + A(i,j)*V(j)

    END DO

END DO

MVMULT_FULL = Sol

END FUNCTION MVMULT_FULL






 !+101+################################################!
!                                                       !
!          Matrix_CCS_MVMult                            !
!                                                       !
 !#####################################################!
SUBROUTINE Matrix_CCS_MVMult( nCol,     &
                              nNZ,      &
                              Row_Ind,  &
                              Col_Ptr,  &
                              Val_Ary,  &
                              x_vec,    &
                              b_vec     )

INTEGER,                                INTENT(IN)  ::  nCol
INTEGER,                                INTENT(IN)  ::  nNZ
INTEGER,        DIMENSION(0:nNZ-1),     INTENT(IN)  ::  Row_Ind
INTEGER,        DIMENSION(0:nCol),      INTENT(IN)  ::  Col_Ptr
COMPLEX(idp),   DIMENSION(0:nNZ-1),     INTENT(IN)  ::  Val_Ary

COMPLEX(idp),   DIMENSION(0:nCol-1),    INTENT(IN)  ::  x_Vec
COMPLEX(idp),   DIMENSION(0:nCol-1),    INTENT(OUT) ::  b_Vec

INTEGER                                             ::  i
INTEGER                                             ::  j
INTEGER                                             ::  index
INTEGER                                             ::  nE_in_Col


b_Vec = 0.0_idp

DO j = 0,nCol-1
nE_in_Col = Col_Ptr(j+1)-Col_Ptr(j)
DO i = 0,nE_in_Col-1

    index = Col_Ptr(j)+i

    b_Vec(row_ind(index)) = b_Vec(row_ind(index))       &
                          + Val_Ary(index)*x_Vec(j)

END DO ! j Loop
END DO ! i Loop

END SUBROUTINE Matrix_CCS_MVMult



 !+101+################################################!
!                                                       !
!          Matrix_CCS_MVMult                            !
!                                                       !
 !#####################################################!
SUBROUTINE Matrix_CCS_MtransVMult(  nCol,     &
                                    nNZ,      &
                                    Row_Ind,  &
                                    Col_Ptr,  &
                                    Val_Ary,  &
                                    x_vec,    &
                                    b_vec     )

INTEGER,                                INTENT(IN)  ::  nCol
INTEGER,                                INTENT(IN)  ::  nNZ
INTEGER,        DIMENSION(0:nNZ-1),     INTENT(IN)  ::  Row_Ind
INTEGER,        DIMENSION(0:nCol),      INTENT(IN)  ::  Col_Ptr
COMPLEX(idp),   DIMENSION(0:nNZ-1),     INTENT(IN)  ::  Val_Ary

COMPLEX(idp),   DIMENSION(0:nCol-1),    INTENT(IN)  ::  x_Vec
COMPLEX(idp),   DIMENSION(0:nCol-1),    INTENT(OUT) ::  b_Vec

INTEGER                                             ::  i
INTEGER                                             ::  j
INTEGER                                             ::  index
INTEGER                                             ::  nE_in_Col


b_Vec = 0.0_idp

DO j = 0,nCol-1
nE_in_Col = Col_Ptr(j+1)-Col_Ptr(j)
DO i = 0,nE_in_Col-1

    index = Col_Ptr(j)+i
    
    b_Vec(j) = b_Vec(j) + Val_Ary(index)*x_Vec(row_ind(index))

END DO ! j Loop
END DO ! i Loop


END SUBROUTINE Matrix_CCS_MtransVMult






END MODULE Functions_Matrix
