   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE  Poisson_Back_Substitution_Module                                 !##!
!##!                                                                         !##!
!##!_________________________________________________________________________!##!
!##!                                                                         !##!
!##!                                                                         !##!
!##!=========================================================================!##!
!##!                                                                         !##!
!##!    Contains:                                                            !##!
!##!                                                                         !##!
!##!                                                                         !##!
!###############################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !##########################################################################!


!*D*================================!
!                                   !
!           Dependencies            !
!                                   !
!===================================!
USE Poseidon_Kinds_Module, &
            ONLY :  idp

IMPLICIT NONE


CONTAINS



!+501+#################################################################
!
!   CCS_Back_Substitution - Special Backward substitution call developled for Poseidon.
!
!       This function performs row-oriented back substituion on L_Transpose, where
!       L is a lower triangular matrix from Cholesky factorization.  Poseidon only
!       stores L, and L is stored using the compressed column storage format. Because
!       of the way this format stores data we can access L_Transpose by treating the
!       L in compressed column format as L_Transpose in compressed row format.
!
!#######################################################################
SUBROUTINE CCS_Back_Substitution(N, NNZ, ELEM_VAL, COL_PTR, ROW_IND, b)

INTEGER, INTENT(IN)                                                 :: N, NNZ

INTEGER, DIMENSION(0:N), INTENT(IN)                                 :: COL_PTR
INTEGER, DIMENSION(0:NNZ-1), INTENT(IN)                             :: ROW_IND
REAL(KIND = idp), DIMENSION(0:NNZ-1), INTENT(IN)                    :: ELEM_VAL

COMPLEX(KIND = idp), DIMENSION(0:N-1), INTENT(INOUT)                :: b



INTEGER                                                             :: i, j, here
REAL(KIND = idp)                                                    :: TMP







DO i = N - 1,0,-1


    TMP = 0             ! To store summation
    here = i + 1        ! To keep track of matrix j value, as j index marks
                        !   equivalent location in CCS elem_val vector.




    DO j = COL_PTR(i)+1, COL_PTR(i+1)-1

        !
        !   Perform summation.  ELEM_VAL(j) => A_(here,i)
        !                                   => A_(i, here)^T
        !       which is what we have stored. I.E. We have stored
        !       a lower triangular matrix, L, in column compressed format.
        !       Here we are performing backward substitution on L transpose.
        !       Since this is a row-wise operation and rows transpose into columns
        !       which perform a switch to do row-wise operation on L transpose by
        !       actually doing column-wise math on L.
        !
        TMP = TMP + ELEM_VAL(j) * b(here)


        here = here + 1

    END DO


    b(i) = ( b(i) - TMP ) / ELEM_VAL(COL_PTR(i))




END DO





END SUBROUTINE CCS_Back_Substitution








END MODULE  Poisson_Back_Substitution_Module
