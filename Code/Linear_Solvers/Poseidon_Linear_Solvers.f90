   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Linear_Solvers_And_Preconditioners                                           !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains the linear solver functions, and preconditioners for the solvers.  !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   PRECOND_CONJ_GRAD_FULL                                              !##!
!##!    +102+   PRECOND_CONJ_GRAD_CCS                                               !##!
!##!                                                                                !##!
!##!    +201+   SSOR_CONDITIONING                                                   !##!
!##!    +202+   JACOBI_CONDITIONING                                                 !##!
!##!    +203+   JACOBI_CONDITIONING_CCS                                             !##!
!##!    +204+   JACOBI_VECTOR_CONDITIONER_CCS                                       !##!
!##!    +205+   JACOBI_VECTOR_CONDITIONER_FULL                                      !##!
!##!    +206+   SSOR_VECTOR_CONDITIONER_FULL                                        !##!
!##!    +207+   SSOR_VECTOR_CONDITIONER_CCS        ! Doesn't Seem to Do Anything !  !##!
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

USE Variables_Derived, &
        ONLY :  Num_R_Nodes




IMPLICIT NONE 





!*F&S*==========================================!
!                                               !
!           Functions & Subroutines             !
!                                               !
!===============================================!
CONTAINS



 !+101+########################################################!
!                                                               !
!          Jacobi_Conditioning                                  !
!                                                               !
 !#############################################################!
SUBROUTINE Jacobi_Conditioning(A,b)

REAL(idp),   DIMENSION(0:NUM_R_NODES-1, 0:NUM_R_NODES-1),    INTENT(INOUT)   :: A
REAL(idp),   DIMENSION(0:NUM_R_NODES -1),                    INTENT(INOUT)   :: b



INTEGER                                                                         :: i

REAL(idp),   DIMENSION(0:NUM_R_NODES -1)                                     :: TMP_VEC
REAL(idp),   DIMENSION(0:NUM_R_NODES-1, 0:NUM_R_NODES-1)                     :: CONDITIONER
REAL(idp),   DIMENSION(0:NUM_R_NODES-1, 0:NUM_R_NODES-1)                     :: Tmp_Mat
REAL(idp)                                                                    :: One
REAL(idp)                                                                    :: Zero

CONDITIONER = 0.0_idp
One = 1.0_idp
Zero = 0.0_idp

DO i = 0,NUM_R_NODES-1
!    PRINT*,"A(i,i)",A(i,i)
    CONDITIONER(i,i) = 1.0_idp/A(i,i)
END DO


CALL ZGEMM('N','N',NUM_R_NODES, NUM_R_NODES, NUM_R_NODES, One, CONDITIONER, NUM_R_NODES, &
                                                A, NUM_R_NODES, Zero,TMP_MAT, NUM_R_NODES )

CALL ZGEMV('N', NUM_R_NODES, NUM_R_NODES, One, CONDITIONER, NUM_R_NODES, b, 1, Zero, TMP_VEC, 1)



b = TMP_VEC
A = TMP_MAT


END SUBROUTINE Jacobi_Conditioning




 !+101+########################################################!
!                                                               !
!          Jacobi_Conditioning_Beta                             !
!                                                               !
 !#############################################################!
SUBROUTINE Jacobi_Conditioning_Beta(A,b, rows, cols )

REAL(idp), DIMENSION(1:Rows, 1:Cols),    INTENT(INOUT)       :: A
REAL(idp), DIMENSION(1:Rows),            INTENT(INOUT)       :: b

INTEGER,                                    INTENT(IN)          :: rows
INTEGER,                                    INTENT(IN)          :: cols

INTEGER                                                         :: i

REAL(idp), DIMENSION(1:cols)                                 :: TMP_VEC
REAL(idp), DIMENSION(1:cols, 1:rows)                         :: CONDITIONER, TMP_MAT
REAL(idp)                                                    :: One, Zero

CONDITIONER = 0.0_idp
One = 1.0_idp
Zero = 0.0_idp

DO i = 1,rows
    CONDITIONER(i,i) = 1.0_idp/A(i,i)
!    PRINT*,CONDITIONER(i,i)
END DO


CALL ZGEMM('N','N',rows, rows, rows, One, CONDITIONER, rows, &
                        A, rows, Zero,TMP_MAT, rows )

CALL ZGEMV('N', rows, rows, One, CONDITIONER, rows, b, 1, Zero, TMP_VEC, 1)



b = TMP_VEC
A = TMP_MAT



END SUBROUTINE Jacobi_Conditioning_Beta








END MODULE Linear_Solvers_And_Preconditioners
