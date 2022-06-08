   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_Cholesky_Module                                                     !##!
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

USE Poseidon_Message_Routines_Module, &
        ONLY :  Init_Message

USE Poseidon_Parameters, &
        ONLY :  DEGREE,             &
                L_LIMIT,            &
                NUM_CFA_EQs,        &
                Verbose_Flag


USE Variables_Derived, &
        ONLY :  Num_R_Nodes

USE Variables_Mesh, &
        ONLY :  Num_R_Elements

USE Variables_BC, &
        ONLY :  INNER_CFA_BC_VALUES,        &
                OUTER_CFA_BC_VALUES,        &
                INNER_CFA_BC_TYPE,          &
                OUTER_CFA_BC_TYPE

USE Variables_FP, &
        ONLY :  First_Column_Storage,   &
                Last_Column_Storage,    &
                Laplace_NNZ,            &
                Factored_NNZ,           &
                Laplace_Matrix_Full,    &
                Laplace_Matrix_VAL,     &
                Laplace_Matrix_ROW,     &
                Laplace_Matrix_COL,     &
                Laplace_Factored_VAL,   &
                Laplace_Factored_ROW,   &
                Laplace_Factored_COL,   &
                Num_Matrices

USE Timer_Routines_Module, &
        ONLY :  TimerStart,             &
                TimerStop

USE Timer_Variables_Module, &
        ONLY :  Timer_XCFC_Matrix_Cholesky

USE Flags_Initialization_Module, &
        ONLY :  lPF_Init_Matrices_Flags,    &
                iPF_Init_Matrices_Type_A_Cholesky

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


INTEGER                                                         ::  l, i, ui
INTEGER                                                         ::  OLD_NNZ, NEW_NNZ

INTEGER(KIND = 1), ALLOCATABLE, DIMENSION(:,:,:)                :: LogMap

INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                          :: NEW_ROW_IND
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)              :: NEW_ELEM_VAL



IF ( Verbose_Flag ) CALL Init_Message('Factorizing scalar Laplace matrix using Cholesky factorization.')
CALL TimerStart(Timer_XCFC_Matrix_Cholesky)



!
!   Determine the number of non-zero elements in the origional stiffness matrix
!
OLD_NNZ = NUM_R_ELEMENTS*(DEGREE + 1)*(DEGREE + 1) - NUM_R_ELEMENTS + 1





                      !                                     !
                     !!                                     !!
                    !!!     Dirichlet Boundary Conditions   !!!
                     !!                                     !!
                      !                                     !

DO ui = 1,Num_Matrices
    !   Dirichlet BCs modify the stiffness matrix so we modify it now.
    !   But to apply the BCs we will need values from the original matrix,
    !   so those values are stored before we modify the matrix.
    !
    IF (INNER_CFA_BC_TYPE(ui) == "D") THEN


        !
        !   This is part of imposing a Dirichlet BC on the FEM system
        !
        Laplace_Factored_VAL(0,:,ui) = 1.0_idp


        !
        !   Save the needed values first.
        !
        First_Column_Storage(1:DEGREE,:,ui) = Laplace_Factored_VAL(1:DEGREE,:,ui)


        !
        !   Now its safe to modify the matrix.
        !
        DO l = 0,L_LIMIT
            Laplace_Factored_VAL(1:Degree,l,ui) = 0.0_idp
            Laplace_Factored_VAL(Laplace_Factored_COL(1:DEGREE,l),l,ui) = 0.0_idp
        END DO

    END IF








    IF (OUTER_CFA_BC_TYPE(ui) == "D") THEN


        !
        ! Save the needed values first
        !
        DO l = 0,L_LIMIT
            DO i = 0,Degree
                Last_Column_Storage(i,l,ui) = Laplace_Factored_VAL(Laplace_Factored_COL(NUM_R_NODES,l) - 1 - i, l,ui)
            END DO
        END DO


        DO i = 1, DEGREE

            !
            !   Now its safe to modify the matrix
            !
            DO l = 0,L_LIMIT
                Laplace_Factored_VAL(Laplace_NNZ - 1 - i, l,ui) = 0.0_idp
                Laplace_Factored_VAL(Laplace_Factored_COL(NUM_R_NODES - i,l) - 1, l,ui) = 0.0_idp
            END DO

        END DO

        Laplace_Factored_VAL(Laplace_NNZ-1,:,ui) = 1.0_idp

    END IF
END DO








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
ALLOCATE(LogMap(0:NUM_R_NODES-1,0:NUM_R_NODES-1,1:Num_Matrices) )




!
!   The process of mapping the matrix occurs in this routine.  It fills the logic map
!   using the 1,0,-1 notation above by a process called Symbolic Factorization. More
!   detials on how this works can be found inside the function.
!
DO ui = 1,Num_Matrices
    CALL Sym_Fact_Initialize(LogMap(:,:,ui),                &
                             NEW_NNZ,               &
                             NUM_R_NODES,           &
                             Laplace_NNZ,           &
                             Laplace_Factored_COL,  &
                             Laplace_Factored_ROW   )
END DO



!
!   Now create space to store the factorized stiffness matrix.
!
ALLOCATE( NEW_ROW_IND(0:NEW_NNZ-1,0:L_LIMIT,1:Num_Matrices),   &
          NEW_ELEM_VAL(0:NEW_NNZ-1,0:L_LIMIT,1:Num_Matrices)    )





!
!   This function transfers the values from the stiffness matrix from the old data space
!   to the newly prepared data space.
!
DO ui = 1,Num_Matrices
    DO l = 0,L_LIMIT

        CALL Sym_Fact_Matrix_Swap(  LogMap(:,:,ui),                 &
                                    NUM_R_NODES,                    &
                                    OLD_NNZ,                        &
                                    NEW_NNZ,                        &
                                    Laplace_Factored_COL(:,l),      &
                                    Laplace_Factored_ROW(:,l),      &
                                    Laplace_Factored_VAL(:,l,ui),   &
                                    NEW_ROW_IND(:,l,ui),            &
                                    NEW_ELEM_VAL(:,l,ui)            )


    END DO
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
DO ui = 1,NUM_MATRICES
    DO l = 0,L_LIMIT


        CALL Chol_Factorizer(NUM_R_NODES, NEW_NNZ, Laplace_Factored_COL(:,l), NEW_ROW_IND(:,L,ui), NEW_ELEM_VAL(:,l,ui))


    END DO
END DO

!PRINT*,"After Chol_Factorizer"
!PRINT*,NEW_ELEM_VAL





                              !                                     !
                             !!                                     !!
                            !!!           Matrix Transfer           !!!
                             !!                                     !!
                              !                                     !

!
!   Now that we no longer need any values from the origional stiffness matrix,
!   we destroy it.
!
DEALLOCATE(Laplace_Factored_VAL, Laplace_Factored_ROW)

!
!   Store new number of non-zero values in the global variable Factored_NNZ
!
Factored_NNZ =  New_NNZ


!                                                                       !
!   Reallocate global space for the new factorized stiffness matrix     !
!                                                                       !
ALLOCATE( Laplace_Factored_VAL(0:Factored_NNZ-1, 0:L_LIMIT, 1:Num_Matrices)  )
ALLOCATE( Laplace_Factored_ROW(0:Factored_NNZ-1, 0:L_LIMIT))


!
!   Transfer the new values from function variables to the global variables.
!
DO ui = 1,Num_Matrices
    Laplace_Factored_VAL(:,:,ui) = NEW_ELEM_VAL(:,:,ui)
    Laplace_Factored_ROW(:,:) = NEW_ROW_IND(:,:,ui)
END DO


!
!   Destory the old storage for the new stiffness matrix
!
DEALLOCATE(NEW_ELEM_VAL, NEW_ROW_IND)

CALL TimerStop(Timer_XCFC_Matrix_Cholesky)
lPF_Init_Matrices_Flags(iPF_Init_Matrices_Type_A_Cholesky) = .TRUE.

!PRINT*,"STOPing at end of CHOLESKY_FACTORIZATION"
!STOP

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
COMPLEX(KIND = idp), DIMENSION(0:OLD_NNZ-1), INTENT(IN)         :: OELEM_VAL

INTEGER, DIMENSION(0:NEW_NNZ-1), INTENT(INOUT)                  :: ROW_IND
COMPLEX(KIND = idp), DIMENSION(0:NEW_NNZ-1), INTENT(INOUT)      :: ELEM_VAL

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
COMPLEX(KIND = idp), DIMENSION(0:NNZ-1), INTENT(INOUT)         :: ELEM_VAL

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
COMPLEX(KIND = idp), DIMENSION(0:NNZ-1), INTENT(INOUT)      ::  ELEM_VAL

INTEGER                                                     ::  i, ij_here, j_row
COMPLEX(KIND = idp)                                            ::  A_jk




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
COMPLEX(KIND = idp), DIMENSION(0:NNZ-1), INTENT(INOUT)          ::  ELEM_VAL

INTEGER                                                         :: i, here



IF ( ROW_IND(COL_PTR(j)) .EQ. j) THEN  !!! DIAGONAL CHECK !!!
    here = COL_PTR(j)


    ELEM_VAL(here) = SQRT(ELEM_VAL(here))


    DO i = COL_PTR(j)+1, COL_PTR(j+1)-1


        ELEM_VAL(i) = ELEM_VAL(i)/ELEM_VAL(here)


    END DO

END IF





END SUBROUTINE Chol_CDIV












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
COMPLEX(KIND = idp), DIMENSION(0:NNZ-1), INTENT(IN)                 :: ELEM_VAL

COMPLEX(KIND = idp), DIMENSION(0:N-1), INTENT(INOUT)                :: b



INTEGER                                                             :: i, j, here
COMPLEX(KIND = idp)                                                 :: TMP







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









!+502+#################################################################
!
!   CCS_Forward_Substitution - This subroutine performs column-wise forward substitution.
!
!#######################################################################
SUBROUTINE CCS_Forward_Substitution(N, NNZ, ELEM_VAL, COL_PTR, ROW_IND, b)

INTEGER, INTENT(IN)                                                 :: N, NNZ

INTEGER, DIMENSION(0:N), INTENT(IN)                                 :: COL_PTR
INTEGER, DIMENSION(0:NNZ-1), INTENT(IN)                             :: ROW_IND
COMPLEX(KIND = idp), DIMENSION(0:NNZ-1), INTENT(IN)                 :: ELEM_VAL

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













END MODULE Poseidon_Cholesky_Module
