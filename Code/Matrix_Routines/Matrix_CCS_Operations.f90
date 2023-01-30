   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Matrix_CCS_Operations_Module                                                 !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
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

IMPLICIT NONE



CONTAINS
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




END MODULE Matrix_CCS_Operations_Module
