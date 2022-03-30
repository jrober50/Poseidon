   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Poisson_Forward_Substitution_Module                                   !##!
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



!+502+#################################################################
!
!   CCS_Forward_Substitution - This subroutine performs column-wise forward substitution.
!
!#######################################################################
SUBROUTINE CCS_Forward_Substitution(N, NNZ, ELEM_VAL, COL_PTR, ROW_IND, b)

INTEGER, INTENT(IN)                                                 :: N, NNZ

INTEGER, DIMENSION(0:N), INTENT(IN)                                 :: COL_PTR
INTEGER, DIMENSION(0:NNZ-1), INTENT(IN)                             :: ROW_IND
REAL(KIND = idp), DIMENSION(0:NNZ-1), INTENT(IN)                    :: ELEM_VAL

COMPLEX(KIND = idp), DIMENSION(0:N-1), INTENT(INOUT)                :: b


INTEGER                                                             :: i, j, here





DO j = 0,N-2

    b(j) = b(j)/ELEM_VAL(COL_PTR(j))

    here = j + 1

    DO i = COL_PTR(j) + 1, COL_PTR(j+1) - 1



        b(here) = b(here) - b(j)*ELEM_VAL(i)
        here = here + 1

    END DO

END DO


b(N-1) = b(N-1)/ELEM_VAL(NNZ-1)




END SUBROUTINE CCS_Forward_Substitution




END MODULE Poisson_Forward_Substitution_Module
