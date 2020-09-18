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


USE Poseidon_Constants_Module, &
        ONLY :  idp, eps

USE Poseidon_Parameters, &
        ONLY :  DEGREE



USE Poseidon_Variables_Module, &
        ONLY :  NUM_R_NODES,        &
                NUM_R_ELEMENTS


USE Poseidon_Matrix_Functions_Module, &
        ONLY :  MVMULT_CCS



IMPLICIT NONE 





!*F&S*==========================================!
!                                               !
!           Functions & Subroutines             !
!                                               !
!===============================================!
CONTAINS
















!+101+#################################################################
!
!   PRECOND_CONJ_GRAD_FULL - Preconditioned Conjugate Gradient method for solving linear systems
!                               using the full matrix (i.e. not sparsity formating)
!
!#######################################################################
SUBROUTINE PRECOND_CONJ_GRAD_FULL(WORK_MAT, WORK_VEC)


COMPLEX(KIND = idp), DIMENSION(0:NUM_R_NODES-1), INTENT(INOUT)                     :: WORK_VEC
COMPLEX(KIND = idp), DIMENSION(0:NUM_R_NODES-1,0:NUM_R_NODES-1), INTENT(IN)        :: WORK_MAT


COMPLEX(KIND = idp), DIMENSION(0:NUM_R_NODES-1, 0:NUM_R_NODES-1)       :: A
COMPLEX(KIND = idp), DIMENSION(0:NUM_R_NODES-1)                        :: b, Z
COMPLEX(KIND = idp), DIMENSION(0:NUM_R_NODES-1)                        :: X

INTEGER                                             :: i, j, k, ierr, MAX_ITER

REAL(KIND = idp)                                    :: alpha, omega

COMPLEX(KIND = idp), DIMENSION(0:NUM_R_NODES-1)            ::  P, Q, R
REAL(KIND = idp)                                        ::  rho, rhoo


MAX_ITER = 1000*NUM_R_NODES
omega = 1.3


!!! INITIAL GUESS AT SOL.  IF NO GUESS AVAILABLE USE X = 0.0 !!!
X = 0.0_idp




A = WORK_MAT
b = WORK_VEC


R = b - MATMUL(A,x)


CALL JACOBI_VECTOR_CONDITIONER_FULL(NUM_R_NODES, A, R, Z)
!CALL SSOR_VECTOR_CONDITIONER_FULL(A, R, Z, omega)


P = Z

rho = DOT_PRODUCT(R,Z)



q = MATMUL(A,p)

alpha = rho/ DOT_PRODUCT(p,q)

X = X + alpha * P
R = R - alpha * Q





k = 1
DO WHILE (( abs(rho) .GE. eps*eps) .AND. (k .LE. MAX_ITER))


    CALL JACOBI_VECTOR_CONDITIONER_FULL(NUM_R_NODES, A, R, Z)
    !CALL SSOR_VECTOR_CONDITIONER_FULL(A, R, Z, omega)



    rhoo = rho
    rho = DOT_PRODUCT(Z,R)
    P = Z + (rho/rhoo)*P
    Q = MATMUL(A,P)
    alpha = rho/DOT_PRODUCT(P,Q)
    X = X + alpha * P
    R = R - alpha * Q


    IF ( k > MAX_ITER) THEN
        PRINT*,"DID NOT CONVERGE IN ", k," STEPS."
    END IF
    k = k + 1

END DO


PRINT*,"CONJ_GRAD_CCS took ",k," iterations"
WORK_VEC = X


END SUBROUTINE PRECOND_CONJ_GRAD_FULL











!+102+##################################################################
!
!   PRECOND_CONJ_GRAD_CCS - Conjugate Gradient method for solving linear systems
!                       stored in Compressed Column Storage Format.
!
!#######################################################################
SUBROUTINE PRECOND_CONJ_GRAD_CCS(N, NNZ, ELEM_VAL, COL_PTR, ROW_IND, B, flag, GUESS, PRECOND_TYPE)

INTEGER, INTENT(IN)                                         :: N, NNZ, flag, PRECOND_TYPE

INTEGER, DIMENSION(0:N), INTENT(IN)                         :: COL_PTR
INTEGER, DIMENSION(0:NNZ-1), INTENT(IN)                     :: ROW_IND

COMPLEX(KIND = idp), DIMENSION(0:NNZ-1), INTENT(IN)            :: ELEM_VAL
!REAL(KIND = idp), DIMENSION(0:N-1), INTENT(INOUT)          :: B, GUESS
COMPLEX(KIND = idp), DIMENSION(0:N-1), INTENT(INOUT)        :: B, GUESS






INTEGER                                                     :: k, MAX_ITER

REAL(KIND = idp)                                            :: rho, rhoo, rhooo, alpha, tol, omega
COMPLEX(KIND = idp), DIMENSION(0:N-1)                          :: X, P, Q, R, Z









MAX_ITER = 10*N
omega = 1.3

!!! INITIAL GUESS AT SOL.  IF NO GUESS AVAILABLE USE X = 0.0 !!!
IF (flag == 0) THEN

    X = 0.0_idp

ELSE IF (flag == 1) THEN

    X = GUESS

END IF


R = B - MVMULT_CCS(N, NNZ, ELEM_VAL, COL_PTR, ROW_IND, X)   ! NEEDED WHEN INITIAL GUESS != 0.0
rho = 1

tol = eps*eps
k = 1
DO WHILE ((MAXVAL(ABS(R)) .GE. eps) .AND. ( abs(rho) .GE. tol)  .AND. (k .LE. MAX_ITER))!

    IF (PRECOND_TYPE == 0) THEN
        CALL JACOBI_VECTOR_CONDITIONER_CCS(N, NNZ, ELEM_VAL, COL_PTR, ROW_IND, R, Z)
    ELSE IF (PRECOND_TYPE == 1) THEN
        CALL SSOR_VECTOR_CONDITIONER_CCS(N, NNZ, ELEM_VAL, COL_PTR, ROW_IND, R, Z, omega)
    END IF

    rhoo = rho
    rho = DOT_PRODUCT(Z, R)

    IF (k == 1) THEN
        P = Z
    ELSE
        P = Z + (rho/rhoo) * P
    END IF

    Q = MVMULT_CCS(N, NNZ, ELEM_VAL, COL_PTR, ROW_IND, P)

    alpha = rho/DOT_PRODUCT(P,Q)
    X = X + alpha * P
    R = R - alpha * Q




    IF ( k >= MAX_ITER) THEN
        PRINT*,"DID NOT CONVERGE IN ", k," STEPS."
    END IF
    k = k + 1


END DO
!PRINT*,"PRECOND_CONJ_GRAD_CCS took ",k," iterations"

B = X





END SUBROUTINE PRECOND_CONJ_GRAD_CCS














!+201+##########################################################################################!
!                                                                                               !
!                                           SSOR_CONDITIONING                                   !
!                                                                                               !
!###############################################################################################!
SUBROUTINE SSOR_CONDITIONING(A, B, omega)


REAL(KIND = idp), INTENT(IN)                                                        :: omega
COMPLEX(KIND = idp), DIMENSION(0:NUM_R_NODES -1), INTENT(INOUT)                        :: B
COMPLEX(KIND = idp), DIMENSION(0:NUM_R_NODES -1, 0:NUM_R_NODES-1), INTENT(INOUT)       :: A


INTEGER                                                                         :: i,j,k, INFO, LWORK
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)                                     :: WORK
INTEGER, DIMENSION(0:NUM_R_NODES-1)                                             :: IPIV

REAL(KIND = idp), DIMENSION(0:NUM_R_NODES -1)                                  :: TMP_VEC
REAL(KIND = idp), DIMENSION(0:NUM_R_NODES-1, 0:NUM_R_NODES-1)                  :: D, L, U, P, C

LWORK = 5*NUM_R_NODES
ALLOCATE(WORK(1:LWORK))

WORK = 0.0_idp
IPIV= 0


D = 0.0;
L = 0.0;
U = 0.0;


DO i = 0,NUM_R_NODES-1



    D(i,i) = 1/A(i,i)
    L(i,i) = A(i,i)
    U(i,i) = A(i,i)

    DO j = 1, DEGREE

        IF (i+j <= NUM_R_NODES-1) THEN
            L(i+j,i) = A(i+j,i)*omega
            U(i,i+j) = A(i,i+j)*omega
        END IF

    END DO


END DO


CALL DGEMM('N','N',NUM_R_NODES, NUM_R_NODES, NUM_R_NODES, 1.0_idp, L, NUM_R_NODES, D, NUM_R_NODES, 0.0_idp, C, NUM_R_NODES )

CALL DGEMM('N','N',NUM_R_NODES, NUM_R_NODES, NUM_R_NODES, 1.0_idp, C, NUM_R_NODES, U, NUM_R_NODES, 0.0_idp, P, NUM_R_NODES )






CALL DGETRF(NUM_R_NODES, NUM_R_NODES, P, NUM_R_NODES, IPIV, INFO)
CALL DGETRI(NUM_R_NODES, P, NUM_R_NODES, IPIV, WORK, LWORK, INFO)



CALL DGEMM('N','N',NUM_R_NODES, NUM_R_NODES, NUM_R_NODES, 1.0_idp, P, NUM_R_NODES, A, NUM_R_NODES, 0.0_idp,C, NUM_R_NODES )


CALL DGEMV('N', NUM_R_NODES, NUM_R_NODES, 1.0_idp, P, NUM_R_NODES, B, 1, 0.0_idp, TMP_VEC, 1)


B = TMP_VEC
A = C



END SUBROUTINE SSOR_CONDITIONING











!+202+##########################################################################################!
!                                                                                               !
!                                           JACOBI CONDITIONING                                 !
!                                                                                               !
!###############################################################################################!
SUBROUTINE JACOBI_CONDITIONING(A,b)




COMPLEX(KIND = idp), DIMENSION(0:NUM_R_NODES -1), INTENT(INOUT)                        :: b
COMPLEX(KIND = idp), DIMENSION(0:NUM_R_NODES -1, 0:NUM_R_NODES-1), INTENT(INOUT)       :: A


INTEGER                                                                         :: i, j

COMPLEX(KIND = idp), DIMENSION(0:NUM_R_NODES -1)                                  :: TMP_VEC
COMPLEX(KIND = idp), DIMENSION(0:NUM_R_NODES-1, 0:NUM_R_NODES-1)                  :: CONDITIONER, TMP_MAT
COMPLEX(KIND = idp)                                                               :: One, Zero

CONDITIONER = 0.0_idp
One = 1.0_idp
Zero = 0.0_idp

DO i = 0,NUM_R_NODES-1
    CONDITIONER(i,i) = 1.0_idp/A(i,i)
END DO


CALL ZGEMM('N','N',NUM_R_NODES, NUM_R_NODES, NUM_R_NODES, One, CONDITIONER, NUM_R_NODES, &
                                                A, NUM_R_NODES, Zero,TMP_MAT, NUM_R_NODES )

CALL ZGEMV('N', NUM_R_NODES, NUM_R_NODES, One, CONDITIONER, NUM_R_NODES, b, 1, Zero, TMP_VEC, 1)



b = TMP_VEC
A = TMP_MAT



END SUBROUTINE JACOBI_CONDITIONING















!+203+#########################################################################################!
!                                                                                               !
!                                       JACOBI CONDITIONING CCS                                 !
!                                                                                               !
!###############################################################################################!
SUBROUTINE JACOBI_CONDITIONING_CCS(N, NNZ, ELEM_VAL, COL_PTR, ROW_IND, WORK_VEC)


INTEGER, INTENT(IN)                                         :: N, NNZ
INTEGER, DIMENSION(0:N), INTENT(IN)                         :: COL_PTR
INTEGER, DIMENSION(0:NNZ-1), INTENT(IN)                     :: ROW_IND


COMPLEX(KIND = idp), DIMENSION(0:NNZ - 1), INTENT(INOUT)       :: ELEM_VAL
COMPLEX(KIND = idp), DIMENSION(0:N-1), INTENT(INOUT)           :: WORK_VEC


INTEGER                                                     :: i, j, COL_HERE
REAL(KIND = idp), DIMENSION(0:N-1)                          :: CONDITIONER



!!!! CALCULATE THE CONDITIONER !!!
CONDITIONER(0) = 1.0_idp/ELEM_VAL(0)
DO i = 0,NUM_R_ELEMENTS-1

    DO j = 1,DEGREE

        COL_HERE = i*DEGREE + j
        CONDITIONER(COL_HERE) = 1.0_idp/ELEM_VAL(COL_PTR(COL_HERE)+j)

    END DO

END DO


!!! APPLY THE CONDITIONER TO THE MATRIX !!!
DO i = 0,NNZ-1
    ELEM_VAL(i) = ELEM_VAL(i)*CONDITIONER(ROW_IND(i))
END DO




!!! APPLY THE CONDITIONER TO THE SOURCE VECTOR !!!
DO i = 0,N-1
    WORK_VEC(i) = WORK_VEC(i) * CONDITIONER(i)
END DO







END SUBROUTINE JACOBI_CONDITIONING_CCS













!+204+##########################################################################################!
!                                                                                               !
!                                       JACOBI_VECTOR_CONDITIONER_CCS                                 !
!                                                                                               !
!###############################################################################################!
SUBROUTINE JACOBI_VECTOR_CONDITIONER_CCS(N, NNZ, ELEM_VAL, COL_PTR, ROW_IND, WORK_VEC, Z)


INTEGER, INTENT(IN)                                         :: N, NNZ
INTEGER, DIMENSION(0:N), INTENT(IN)                         :: COL_PTR
INTEGER, DIMENSION(0:NNZ-1), INTENT(IN)                     :: ROW_IND


COMPLEX(KIND = idp), DIMENSION(0:NNZ - 1), INTENT(IN)          :: ELEM_VAL
COMPLEX(KIND = idp), DIMENSION(0:N-1), INTENT(IN)              :: WORK_VEC
COMPLEX(KIND = idp), DIMENSION(0:N-1), INTENT(INOUT)           :: Z


INTEGER                                                     :: i, j, COL_HERE
REAL(KIND = idp), DIMENSION(0:N-1)                          :: CONDITIONER, CND





!!!! CALCULATE THE CONDITIONER !!!
Z(0) = WORK_VEC(0)/ELEM_VAL(0)
DO i = 0,NUM_R_ELEMENTS-1

    DO j = 1,DEGREE

        COL_HERE = i*DEGREE + j

        Z(COL_HERE) = WORK_VEC(COL_HERE)/ELEM_VAL(COL_PTR(COL_HERE)+j)

    END DO

END DO





END SUBROUTINE JACOBI_VECTOR_CONDITIONER_CCS





!+205+#########################################################################################!
!                                                                                               !
!                                       JACOBI_VECTOR_CONDITIONER_FULL                                !
!                                                                                               !
!###############################################################################################!
SUBROUTINE JACOBI_VECTOR_CONDITIONER_FULL(N, A, b, Z)


INTEGER, INTENT(IN)                                         :: N

COMPLEX(KIND = idp), DIMENSION(0:N - 1), INTENT(IN)                :: b
COMPLEX(KIND = idp), DIMENSION(0:N - 1, 0:N - 1), INTENT(IN)       :: A

COMPLEX(KIND = idp), DIMENSION(0:N - 1), INTENT(INOUT)             :: Z




INTEGER                                                         :: i




!!! APPLY THE CONDITIONER TO THE SOURCE VECTOR !!!
DO i = 0,N-1


    Z(i) = b(i)/A(i,i)


END DO


END SUBROUTINE JACOBI_VECTOR_CONDITIONER_FULL















!+206+#########################################################################################!
!                                                                                               !
!                                           SSOR_CONDITIONING                                   !
!                                                                                               !
!###############################################################################################!
SUBROUTINE SSOR_VECTOR_CONDITIONER_FULL(A, B, Z, omega)


REAL(KIND = idp), INTENT(IN)                                                        :: omega
COMPLEX(KIND = idp), DIMENSION(0:NUM_R_NODES -1), INTENT(IN)                           :: B
COMPLEX(KIND = idp), DIMENSION(0:NUM_R_NODES -1, 0:NUM_R_NODES-1), INTENT(IN)          :: A

COMPLEX(KIND = idp), DIMENSION(0:NUM_R_NODES -1), INTENT(INOUT)                        :: Z


INTEGER                                                                         :: i,j,k, INFO, LWORK
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)                                     :: WORK
INTEGER, DIMENSION(0:NUM_R_NODES-1)                                             :: IPIV

REAL(KIND = idp), DIMENSION(0:NUM_R_NODES -1)                                  :: TMP_VEC
REAL(KIND = idp), DIMENSION(0:NUM_R_NODES-1, 0:NUM_R_NODES-1)                  :: invD, L, P, C

LWORK = 5*NUM_R_NODES
ALLOCATE(WORK(1:LWORK))

WORK = 0.0_idp
IPIV= 0



DO i = 0,NUM_R_NODES-1


    invD(i,i) = 1/A(i,i)
    L(i,i) = A(i,i)

    DO j = 1, DEGREE

        IF (i+j <= NUM_R_NODES-1) THEN
            L(i+j,i) = A(i+j,i)*omega
        END IF

    END DO


END DO






CALL DGEMM('N','N',NUM_R_NODES, NUM_R_NODES, NUM_R_NODES, 1.0_idp, L, NUM_R_NODES, invD, NUM_R_NODES, 0.0_idp, C, NUM_R_NODES )

CALL DGEMM('N','T',NUM_R_NODES, NUM_R_NODES, NUM_R_NODES, 1.0_idp, C, NUM_R_NODES, L, NUM_R_NODES, 0.0_idp, P, NUM_R_NODES )





Z = B
CALL DGESV(NUM_R_NODES, 1, P, NUM_R_NODES, IPIV, Z, NUM_R_NODES, INFO)
IF (INFO > 0) THEN
    print*,"DGESV has failed with INFO = ",INFO
END IF








!!! INVERT THE MATRIX !!!!
!CALL DGETRF(NUM_R_NODES, NUM_R_NODES, P, NUM_R_NODES, IPIV, INFO)
!CALL DGETRI(NUM_R_NODES, P, NUM_R_NODES, IPIV, WORK, LWORK, INFO)

!CALL DGEMV('N', NUM_R_NODES, NUM_R_NODES, 1.0_idp, P, NUM_R_NODES, B, 1, 0.0_idp, Z, 1)







END SUBROUTINE SSOR_VECTOR_CONDITIONER_FULL









!+207+#########################################################################################!
!                                                                                               !
!                                           SSOR_CONDITIONING                                   !
!                                                                                               !
!###############################################################################################!
SUBROUTINE SSOR_VECTOR_CONDITIONER_CCS(N, NNZ, ELEM_VAL, COL_PTR, ROW_IND, WORK_VEC, Z, omega)


INTEGER, INTENT(IN)                                         :: N, NNZ
INTEGER, DIMENSION(0:N), INTENT(IN)                         :: COL_PTR
INTEGER, DIMENSION(0:NNZ-1), INTENT(IN)                     :: ROW_IND


REAL(KIND = idp), INTENT(IN)                                :: omega
COMPLEX(KIND = idp), DIMENSION(0:NNZ - 1), INTENT(IN)          :: ELEM_VAL
COMPLEX(KIND = idp), DIMENSION(0:N-1), INTENT(IN)              :: WORK_VEC
COMPLEX(KIND = idp), DIMENSION(0:N-1), INTENT(INOUT)           :: Z



INTEGER                                                     :: i,j,k


REAL(KIND = idp), DIMENSION(0:N - 1)                        :: TMP_VEC
REAL(KIND = idp), DIMENSION(0:NNZ - 1)                      :: LinvD, L, P, C



DO i = 0,NUM_R_ELEMENTS - 1


    DO j = COL_PTR(i),COL_PTR(i+1)-1






    END DO


END DO










END SUBROUTINE SSOR_VECTOR_CONDITIONER_CCS









END MODULE Linear_Solvers_And_Preconditioners
