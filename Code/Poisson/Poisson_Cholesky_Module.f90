   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poisson_Cholesky_Module                                                      !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains subroutines related to Cholesky Factorization.                     !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   Cholesky_Factorization                                              !##!
!##!                                                                                !##!
!##!    +201+   Sym_Fact_Initialize                                                 !##!
!##!    +202+   Sym_Fact_Matrix_Swap                                                !##!
!##!    +203+   Sym_Fact_Parent                                                     !##!
!##!                                                                                !##!
!##!    +301+   Chol_Factorizer                                                     !##!
!##!                                                                                !##!
!##!    +401+   Chol_CMOD                                                           !##!
!##!    +402+   Chol_CDIV                                                           !##!
!##!                                                                                !##!
!##!    +501+   CCS_Back_Substitution                                               !##!
!##!    +502+   CCS_Forward_Substitution                                            !##!
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

USE Poseidon_Parameters, &
            ONLY :  Degree,             &
                    L_Limit,            &
                    Verbose_Flag



USE Variables_Poisson, &
            ONLY :  Matrix_Cholesky_Factorized_Flag,    &
                    First_Column_Storage,               &
                    Last_Column_Storage

USE Variables_BC, &
            ONLY :  INNER_Poisson_BC_TYPE,              &
                    OUTER_Poisson_BC_TYPE


USE Variables_Mesh, &
            ONLY :  Num_R_Elements

USE Variables_Derived, &
            ONLY :  Num_R_Nodes

USE Variables_Poisson, &
            ONLY :  STF_NNZ,            &
                    STF_Mat_Integrals,  &
                    STF_Elem_Val,       &
                    STF_Row_Ind,        &
                    STF_Col_Ptr
        



IMPLICIT NONE






!*F&S*==========================================!
!                                               !
!           Functions & Subroutines             !
!                                               !
!===============================================!
CONTAINS

!+101+#################################################################
!
!   Cholesky_Factorization - This subroutine performs cholesky factorization
!       on the stiffness matrix.  To do this a few steps must be taken to
!       faciliate the factorization, as well as allow other parts of Poseidon
!       to function properly.
!
!#######################################################################
SUBROUTINE Cholesky_Factorization()


INTEGER                                                             ::  l, i

INTEGER                                                             ::  OLD_NNZ, NEW_NNZ



INTEGER(KIND = 1), ALLOCATABLE, DIMENSION(:,:)                      :: LogMap

INTEGER, ALLOCATABLE, DIMENSION(:)                                  :: NEW_ROW_IND
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)                       :: NEW_ELEM_VAL



IF ( Verbose_Flag ) THEN
    PRINT*," - Begining Cholesky Factorization."
END IF




!
!   Determine the number of non-zero elements in the origional stiffness matrix
!
OLD_NNZ = STF_NNZ









                      !                                     !
                     !!                                     !!
                    !!!     Dirichlet Boundary Conditions   !!!
                     !!                                     !!
                      !                                     !


!   Dirichlet BCs modify the stiffness matrix so we modify it now.
!   But to apply the BCs we will need values from the original matrix,
!   so those values are stored before we modify the matrix.
!
IF ( INNER_Poisson_BC_TYPE .EQ. "D") THEN

    !
    !   This is part of imposing a Dirichlet BC on the FEM system
    !
    STF_ELEM_VAL(0,:) = 1.0_idp


    DO i = 1, DEGREE

        !
        !   Save the needed values first.
        !
        First_Column_Storage(i,:) = STF_ELEM_VAL(i, :)


        !
        !   Now its safe to modify the matrix.
        !
        STF_ELEM_VAL(i,:) = 0.0_idp
        STF_ELEM_VAL(STF_COL_PTR(i),:) = 0.0_idp

    END DO




END IF







IF ( Outer_Poisson_BC_TYPE .EQ. "D") THEN


    DO i = 1, DEGREE

        !
        ! Save the needed values first
        !
        Last_Column_Storage(i,:) = STF_ELEM_VAL(STF_COL_PTR(NUM_R_NODES) - 1 - i, :)

        !
        !   Now its safe to modify the matrix
        !
        STF_ELEM_VAL(OLD_NNZ - 1 - i, :) = 0.0_idp
        STF_ELEM_VAL(STF_COL_PTR(NUM_R_NODES - i) - 1, :) = 0.0_idp


    END DO

    STF_ELEM_VAL(OLD_NNZ-1,:) = 1.0_idp


END IF







                              !                                     !
                             !!                                     !!
                            !!!         Symbolic Factorization      !!!
                             !!                                     !!
                              !                                     !
!
!   During Cholesky Factorization, some zero-valued elements in the stiffnes
!   matrix will become non-zero.  As this code only allocates space to store
!   the non-zero values of the origional matrix, this presents a problem, as
!   there is no place to store the new ones. To make sparse cholesky factorization
!   possible, we need to reallocate space for the new non-zero elements.
!
!   This first step in this process is to determine where these non-zero values
!   will appear. This is what occurs in the following lines of code.
!





!
!   First, a logic map is given space.  This map will use 1,0, and -1 to map the
!   locations of the origional zeros (0), and non-zeros (1) and the locations
!   where a non-zero will appear (-1).
!
ALLOCATE(LogMap(0:NUM_R_NODES-1,0:NUM_R_NODES-1) )




!
!   The process of mapping the matrix occurs in this routine.  It fills the logic map
!   using the 1,0,-1 notation above by a process called Symbolic Factorization. More
!   detials on how this works can be found inside the function.
!
CALL Sym_Fact_Initialize(LogMap, NEW_NNZ, NUM_R_NODES, OLD_NNZ, STF_COL_PTR, STF_ROW_IND)





!
!   Now create space to store the factorized stiffness matrix.
!
ALLOCATE( NEW_ROW_IND(0:NEW_NNZ-1), NEW_ELEM_VAL(0:NEW_NNZ-1,0:L_LIMIT))





!
!   This function transfers the values from the stiffness matrix from the old data space
!   to the newly prepared data space.
!
DO l = 0,L_LIMIT

    CALL Sym_Fact_Matrix_Swap(LogMap, NUM_R_NODES, OLD_NNZ, NEW_NNZ, STF_COL_PTR,           &
                            STF_ROW_IND, STF_ELEM_VAL, NEW_ROW_IND, NEW_ELEM_VAL(:,l))

END DO




!
!   Now that the matrix has been transfered to its new home, we no longer
!   need the logic map.
!
DEALLOCATE(LogMap)











                              !                                     !
                             !!                                     !!
                            !!!        Cholesky Factorization       !!!
                             !!                                     !!
                              !                                     !


!
!   Here the actual factorization of the matrix occurs.  The subroutine called
!   performs Cholesky factorization which produces a lower triangular matrix, L
!   such that
!               A = L * L_Transpose
!
!   where A is the origional stiffness matrix.  This allows one to solve the linear
!   system using backwards and forwards substitution.
!
DO l = 0,L_LIMIT


    CALL Chol_Factorizer(NUM_R_NODES, NEW_NNZ, STF_COL_PTR, NEW_ROW_IND, NEW_ELEM_VAL(:,l))


END DO








                              !                                     !
                             !!                                     !!
                            !!!           Matrix Transfer           !!!
                             !!                                     !!
                              !                                     !

!
!   Now that we no longer need any values from the origional stiffness matrix,
!   we destroy it.
!
DEALLOCATE(STF_ELEM_VAL, STF_ROW_IND)


!
!   Store new number of non-zero values in the global variable STF_NNZ
!
STF_NNZ =  STF_COL_PTR(NUM_R_NODES)


!                                                                       !
!   Reallocate global space for the new factorized stiffness matrix     !
!                                                                       !
ALLOCATE(STF_ELEM_VAL(0:STF_NNZ-1, 0:L_LIMIT), STF_ROW_IND(0:STF_NNZ-1))


!
!   Transfer the new values from function variables to the global variables.
!
STF_ELEM_VAL = NEW_ELEM_VAL
STF_ROW_IND = NEW_ROW_IND


!
!   Destory the old storage for the new stiffness matrix
!
DEALLOCATE(NEW_ELEM_VAL, NEW_ROW_IND)













END SUBROUTINE CHOLESKY_FACTORIZATION





!+201+######################################################################################!
!                                                                                           !
!                                       Sym_Fact_Initialize                                  !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Sym_Fact_Initialize(LogMap,  NEW_NNZ, N, OLD_NNZ, COL_PTR, ROW_IND)

INTEGER, INTENT(IN)                                                     :: N, OLD_NNZ
INTEGER, DIMENSION(0:N), INTENT(IN)                                     :: COL_PTR
INTEGER, DIMENSION(0:OLD_NNZ-1), INTENT(IN)                             :: ROW_IND

INTEGER, INTENT(INOUT)                                                  :: NEW_NNZ
INTEGER(KIND = 1), DIMENSION(0:N-1,0:N-1), INTENT(INOUT)                :: LogMap

INTEGER                                                                 :: k, i, j, row_here
INTEGER                                                                 :: PARENT


!PRINT*,"LOGMAP_INIT FIRST COUNT"

!
!   Map nonzeros to logic map
!   1 -> nonzero value
!   0 -> zero value
!
LogMap = 0
DO i = 0,N-1

    DO j = COL_PTR(i), COL_PTR(i+1)-1

        row_here = ROW_IND(j)

        LogMap(row_here,i) = 1

    END DO

END DO




!PRINT*,"LOGMAP_INIT PARENT COUNT"
!
!   Use logical map to determine parentage of
!   matrix elements, and mark where new non-zero
!   values will appear during factorization.
!
DO k = 0, N-1

    PARENT = Sym_Fact_Parent(k, N, LogMap)

    DO i = k+1,N-1


        IF (LogMap(i,k) .NE. 0) THEN

            LogMap(i,PARENT) = MAX(2*LogMap(i,PARENT) - 1, -1)

        END IF

    END DO


END DO




!PRINT*,"LOGMAP_INIT NEW_NNZ COUNT"
!
! Count number of nonzeros including and below the diagonal.
!   Only count below the diagonal as we will only need to store the
!   the elements below the diagonal as the matrix is symmetric. A
!   second benefit to this is we seek to construct a lower triangular
!   matrix through cholesky factorization, which by definition only
!   has values below the diagonal.
!
NEW_NNZ = 0
DO i = 0,N-1

    NEW_NNZ = NEW_NNZ + COUNT(LogMap(i,i:N-1) .NE. 0)

END DO


END SUBROUTINE Sym_Fact_Initialize





!+202+######################################################################################!
!                                                                                           !
!                                       Sym_Fact_Matrix_Swap                                    !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Sym_Fact_Matrix_Swap(LogMap, N, OLD_NNZ, NEW_NNZ, COL_PTR, OROW_IND, OELEM_VAL, ROW_IND, ELEM_VAL)

INTEGER, INTENT(IN)                                             :: N, OLD_NNZ, NEW_NNZ
INTEGER, DIMENSION(0:N), INTENT(INOUT)                          :: COL_PTR
INTEGER, DIMENSION(0:OLD_NNZ-1), INTENT(IN)                     :: OROW_IND
REAL(KIND = idp), DIMENSION(0:OLD_NNZ-1), INTENT(IN)            :: OELEM_VAL

INTEGER, DIMENSION(0:NEW_NNZ-1), INTENT(INOUT)                  :: ROW_IND
REAL(KIND = idp), DIMENSION(0:NEW_NNZ-1), INTENT(INOUT)         :: ELEM_VAL

INTEGER(KIND = 1), DIMENSION(0:N-1,0:N-1), INTENT(INOUT)        :: LogMap

INTEGER                                                         :: i, j
INTEGER                                                         :: counter, here, fhere





ELEM_VAL = 0
COL_PTR = 0
ROW_IND = 0


here = 0
fhere = 0
DO i = 0,N-1

    counter = 0


    DO j = 0,i-1


        !
        !   Here we count the non-zeros above the diagonal.  As
        !   the SCCS format only stores values below the diagonal
        !   we do not count these towards the new non-zero total,
        !   but we must shift our point of reference (here) to the next
        !   location.
        !
        if (LogMap(j,i) .EQ. 1) THEN

            here = here + 1
        END IF

    END DO



    DO j = i,N-1

        !
        !   Now we are at and below the diagonal.  Here we count
        !   every non-zero (indicated by a non-zero value of the LogMap).
        !
        IF (LogMap(j,i) .NE. 0 ) THEN
            counter = counter + 1
        END IF

    END DO



    !
    !   This updates the Column pointer vector by adding the total number of
    !   non-zeros of the current row to all the vector values after the
    !   the row we are investigating currently.
    !
    COL_PTR(i+1:N) = COL_PTR(i+1:N) + counter



    j = i
    DO WHILE (counter > 0)



        IF (LogMap(j,i) .EQ. 1) THEN
            !
            !   IF LogMap(j,i) == 1 then that means this location already
            !   had a non-zero value to begin with.  Therefore we copy the value
            !   from the old matrix, to the new, and indicate its row appropriotely
            !

            ELEM_VAL(fhere) = OELEM_VAL(here)
            ROW_IND(fhere) = j



            here = here + 1
            fhere = fhere + 1
            counter = counter - 1

        ELSE IF (LogMap(j,i) .EQ. -1) THEN
            !
            !   If LogMap(j,i) == -1 then this location is currently a zero
            !   but will become non-zero during the factorization process
            !   so we need to create space in the new matrix for it, initialize
            !   the value to zero, and indicate its row location appropriotely.
            !
            ELEM_VAL(fhere) = 0.0_idp
            ROW_IND(fhere) = j

            fhere = fhere + 1
            counter = counter - 1

        END IF


        j = j + 1

    END DO



END DO





END SUBROUTINE Sym_Fact_Matrix_Swap























!+203+######################################################################################!
!                                                                                           !
!                 Sym_Fact_Parent - function that returns parent of input col                !
!                                                                                           !
!###########################################################################################!
FUNCTION Sym_Fact_Parent(j, N, LogMap)

INTEGER                                                     :: Sym_Fact_Parent

INTEGER, INTENT(IN)                                         :: j, N
INTEGER(KIND = 1), DIMENSION(0:N-1,0:N-1), INTENT(IN)       :: LogMap

INTEGER                                                     :: here_j



here_j = j+1

IF ( here_j < N - 1 ) THEN

    DO WHILE( LogMap(here_j, j) .EQ. 0 .AND. here_j < N-1)

        here_j = here_j + 1

    END DO

END IF




IF (here_j >= N-1) THEN

    Sym_Fact_Parent = j

ELSE

    Sym_Fact_Parent = here_j

END IF





END FUNCTION Sym_Fact_Parent













!+301+######################################################################################!
!                                                                                           !
!                       Chol_Factorizer - Naive Column - Cholesky Factorization                !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Chol_Factorizer(N, NNZ, COL_PTR, ROW_IND, ELEM_VAL)

INTEGER, INTENT(IN)                                         :: N, NNZ
INTEGER, DIMENSION(0:N), INTENT(IN)                         :: COL_PTR
INTEGER, DIMENSION(0:NNZ-1), INTENT(IN)                     :: ROW_IND
REAL(KIND = idp), DIMENSION(0:NNZ-1), INTENT(INOUT)         :: ELEM_VAL

INTEGER                                                     :: i,j,k


DO j = 0,N-1


    DO k  = 0,j-1


        DO i = COL_PTR(k), COL_PTR(k+1)-1


            IF (ROW_IND(i) .EQ. j) THEN

                CALL Chol_CMOD(N, NNZ, COL_PTR, ROW_IND, ELEM_VAL, i, k)


            END IF

        END DO

    END DO

    call Chol_CDIV(N, NNZ, COL_PTR, ROW_IND, ELEM_VAL, j)


END DO




END SUBROUTINE Chol_Factorizer









!+401+##############################################################################!
!                                                                                   !
!         Chol_CMOD                                                                 !
!                                                                                   !
!###################################################################################!
SUBROUTINE Chol_CMOD(N, NNZ, COL_PTR, ROW_IND, ELEM_VAL, j, k)

INTEGER, INTENT(IN)                                         ::  N, NNZ, j, k

INTEGER, DIMENSION(0:N),INTENT(IN)                          ::  COL_PTR
INTEGER, DIMENSION(0:NNZ-1), INTENT(IN)                     ::  ROW_IND
REAL(KIND = idp), DIMENSION(0:NNZ-1), INTENT(INOUT)         ::  ELEM_VAL

INTEGER                                                     ::  i, ij_here, j_row
REAL(KIND = idp)                                            ::  A_jk




A_jk = ELEM_VAL(j)





IF (A_jk .NE. 0.0_idp ) THEN


    j_row = ROW_IND(j)
    ij_here = COL_PTR(ROW_IND(j))



    DO i = j, COL_PTR(k+1)-1




        DO WHILE(ROW_IND(ij_here) .NE. ROW_IND(i) .AND. ij_here < COL_PTR(j_row + 1))

            ij_here = ij_here+1

        END DO



        ELEM_VAL(ij_here) = ELEM_VAL(ij_here) - ELEM_VAL(i)*A_jk




    END DO
END IF





END SUBROUTINE Chol_CMOD







!+402+##############################################################################!
!                                                                                   !
!         Chol_CDIV                                                                 !
!                                                                                   !
!###################################################################################!
SUBROUTINE Chol_CDIV(N, NNZ, COL_PTR, ROW_IND, ELEM_VAL, j)

INTEGER, INTENT(IN)                                             ::  N, NNZ, j

INTEGER, DIMENSION(0:N),INTENT(IN)                              ::  COL_PTR
INTEGER, DIMENSION(0:NNZ-1), INTENT(IN)                         ::  ROW_IND
REAL(KIND = idp), DIMENSION(0:NNZ-1), INTENT(INOUT)             ::  ELEM_VAL

INTEGER                                                         :: i, here



IF ( ROW_IND(COL_PTR(j)) .EQ. j) THEN  !!! DIAGONAL CHECK !!!
    here = COL_PTR(j)



    ELEM_VAL(here) = SQRT(ELEM_VAL(here))


    DO i = COL_PTR(j)+1, COL_PTR(j+1)-1


        ELEM_VAL(i) = ELEM_VAL(i)/ELEM_VAL(here)


    END DO

END IF





END SUBROUTINE Chol_CDIV






END MODULE  Poisson_Cholesky_Module
