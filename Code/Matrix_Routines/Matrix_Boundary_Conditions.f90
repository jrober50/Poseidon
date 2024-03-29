   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Matrix_Boundary_Condition_Routines                                           !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!   Contains the subroutines used to calculate and impliment boundary conditions !##!
!##!    on the system being solved.                                                 !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   DIRICHLET_BC                                                        !##!
!##!    +101+   NEUMANN_BC                                                          !##!
!##!                                                                                !##!
!##!    +201+   DIRICHLET_BC_CSS                                                    !##!
!##!    +202+   NEUMANN_BC_CSS                                                      !##!
!##!                                                                                !##!
!##!    +301+   DIRICHLET_BC_CHOL                                                   !##!
!##!                                                                                !##!
!##!    +401+   BC_Integral                                                         !##!
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
                    
USE Poseidon_Numbers_Module, &
            ONLY :  pi


USE Poseidon_Parameters, &
            ONLY :  Degree,                     &
                    L_Limit,                    &
                    Num_Eqs

USE Variables_Matrices, &
            ONLY :  dMA_First_Col_Storage,       &
                    dMA_Last_Col_Storage,        &
                    dMB_First_Col_Storage,      &
                    dMB_Last_Col_Storage

USE Variables_Mesh, &
            ONLY :  Num_R_Elements,             &
                    R_INNER,                    &
                    R_OUTER

USE Variables_Derived, &
            ONLY :  Num_R_Nodes,                &
                    iVB_Prob_Dim,               &
                    iVB_Elem_Prob_Dim,          &
                    LM_Length

USE Variables_BC, &
            ONLY :  Inner_BC_Values,        &
                    Outer_BC_Values,        &
                    Inner_BC_Type,          &
                    Outer_BC_Type

USE Maps_Fixed_Point, &
            ONLY :  FP_Beta_Array_Map


IMPLICIT NONE

!*F&S*==========================================!
!                                               !
!           Functions & Subroutines             !
!                                               !
!===============================================!
CONTAINS




 !+101+############################################################################!
!                                                                                   !
!       Dirichlet_BC                                                                !
!                                                                                   !
!-----------------------------------------------------------------------------------!
!                                                                                   !
!   Input                                                                           !
!                                                                                   !
 !#################################################################################!
SUBROUTINE DIRICHLET_BC(WORK_MAT, WORK_VEC, L, M, ui)

REAL(idp),  DIMENSION(1:NUM_R_NODES,1:NUM_R_NODES), INTENT(INOUT)   :: WORK_MAT
REAL(idp),  DIMENSION(1:NUM_R_NODES),               INTENT(INOUT)   :: WORK_VEC
INTEGER,                                            INTENT(IN)      :: L, M, ui

INTEGER                             :: i, shift, uj
REAL(idp)                           :: BC_Value



IF (Inner_BC_Type(ui) == "D") THEN

    BC_Value =  2.0_idp*sqrt(pi)*Inner_BC_Values(ui)
    
    WORK_VEC(1) = BC_Value
    DO i = 2,DEGREE+1
        WORK_VEC(i) = WORK_VEC(i) - WORK_MAT(1,i)*BC_Value
    END DO

    WORK_MAT(1,:) = 0.0_idp
    WORK_MAT(:,1) = 0.0_idp
    WORK_MAT(1,1) = 1.0_idp

END IF




shift = 0
IF (NUM_R_ELEMENTS .EQ. 1 ) THEN
    shift = 1
END IF

IF (Outer_BC_Type(ui)  == "D") THEN

    IF ( ( L == 0 ) .AND. ( M == 0 )  ) THEN
        BC_Value = 2.0_idp*sqrt(pi)*Outer_BC_Values(ui)
    ELSE
        BC_Value = 0.0_idp
    END IF

    WORK_VEC(NUM_R_NODES ) = BC_Value

    DO i = 1,DEGREE-shift


        WORK_VEC(NUM_R_NODES-i) = WORK_VEC(NUM_R_NODES-i)       &
                                  - WORK_MAT(NUM_R_NODES, NUM_R_NODES-i)*BC_Value
                                    

    END DO


    WORK_MAT(NUM_R_NODES,:) = 0.0_idp
    WORK_MAT(:,NUM_R_NODES) = 0.0_idp
    WORK_MAT(NUM_R_NODES,NUM_R_NODES) = 1.0_idp

END IF



END SUBROUTINE DIRICHLET_BC




 !+102+############################################################################!
!                                                                                   !
!       Neumann_BC                                                                !
!                                                                                   !
!-----------------------------------------------------------------------------------!
!                                                                                   !
!   Input                                                                           !
!                                                                                   !
 !#################################################################################!
SUBROUTINE NEUMANN_BC(L_VALUE, WORK_VEC)

INTEGER,                                    INTENT(IN)          ::  L_VALUE
REAL(idp),  DIMENSION(0:NUM_R_NODES - 1),   INTENT(INOUT)       ::  WORK_VEC


INTEGER                 :: ui






DO ui = 1,Num_Eqs
    IF (    L_VALUE == 0    )   THEN


        IF (Inner_BC_Type(ui) == "N") THEN

            WORK_VEC(0) = WORK_VEC(0) - sqrt(4.0_idp*pi)*Inner_BC_Values(ui)

        END IF


        IF (Outer_BC_Type(ui) == "N") THEN

            WORK_VEC(NUM_R_NODES-1) = WORK_VEC(NUM_R_NODES-1) + R_OUTER*R_OUTER*Outer_BC_Values(ui)

        END IF



    END IF
END DO

END SUBROUTINE NEUMANN_BC















 !+201+############################################################################!
!                                                                                   !
!       Dirichlet_BC_CCS                                                            !
!                                                                                   !
!-----------------------------------------------------------------------------------!
!                                                                                   !
!   Input                                                                           !
!                                                                                   !
 !#################################################################################!
SUBROUTINE DIRICHLET_BC_CCS(N, NNZ, L, M, ELEM_VAL, COL_PTR, ROW_IND, WORK_VEC)


INTEGER,                                INTENT(IN)                  ::  N, NNZ, L, M

INTEGER, DIMENSION(0:N),                INTENT(IN)                  ::  COL_PTR
INTEGER, DIMENSION(0:NNZ-1),            INTENT(IN)                  ::  ROW_IND

REAL(idp), DIMENSION(0:N - 1),          INTENT(INOUT)               ::  WORK_VEC


REAL(idp), DIMENSION(0:NNZ-1),          INTENT(INOUT)               ::  ELEM_VAL


INTEGER                                                             :: i, shift, ui
REAL(idp)                                                           :: BC_Value





DO ui = 1,Num_Eqs
    IF (Inner_BC_Type(ui) == "D") THEN



        BC_Value = sqrt(4.0_idp*pi)*Inner_BC_Values(ui)



        !!! MODIFY SOURCE VECTOR !!!


        WORK_VEC(0) = BC_Value

        DO i = COL_PTR(0)+1,COL_PTR(1)-1

            WORK_VEC(i) = WORK_VEC(i) - ELEM_VAL(i)*BC_Value

        END DO



        !!! MODIFY MATRIX !!!

        ELEM_VAL(0) = 1.0_idp

        DO i = 1,DEGREE

            ELEM_VAL(i) = 0.0_idp
            ELEM_VAL(COL_PTR(i)) = 0.0_idp


        END DO




    END IF









    shift = 0
    IF (NUM_R_ELEMENTS .EQ. 1 ) THEN
        shift = 1
    END IF

    IF (Outer_BC_Type(ui)  == "D") THEN



        BC_Value = sqrt(4.0_idp*pi)*Outer_BC_Values(ui)

        !!! MODIFY SRC VECTOR !!!
        WORK_VEC(NUM_R_NODES - 1) = BC_Value

        DO i = 1,DEGREE-shift

            WORK_VEC(NUM_R_NODES - 1 - i) = WORK_VEC(NUM_R_NODES - 1 - i)                           &
                                            - ELEM_VAL(COL_PTR(N) - 1 -i)*BC_Value

        END DO



        !!! MODIFY MATRIX !!!
        DO i = 1,DEGREE

            ELEM_VAL(NNZ - 1 - i) = 0.0_idp
            ELEM_VAL(COL_PTR(NUM_R_NODES - i )-1) = 0.0_idp


        END DO

        ELEM_VAL(NNZ-1) = 1.0_idp


    END IF


END DO







END SUBROUTINE DIRICHLET_BC_CCS


















 !+202+############################################################################!
!                                                                                   !
!       Neumann_BC_CCS                                                              !
!                                                                                   !
!-----------------------------------------------------------------------------------!
!                                                                                   !
!   Input                                                                           !
!                                                                                   !
 !#################################################################################!
SUBROUTINE NEUMANN_BC_CCS(N, NNZ, L, M, ELEM_VAL, COL_PTR, ROW_IND, WORK_VEC)


INTEGER,                                INTENT(IN)                  ::  N, NNZ, L, M

INTEGER,    DIMENSION(0:N),             INTENT(IN)                  ::  COL_PTR
INTEGER,    DIMENSION(0:NNZ-1),         INTENT(IN)                  ::  ROW_IND

REAL(idp),  DIMENSION(0:N - 1),         INTENT(INOUT)               ::  WORK_VEC


REAL(idp),  DIMENSION(0:NNZ-1),         INTENT(INOUT)               ::  ELEM_VAL


INTEGER                                                             ::  ui

REAL(KIND = idp)                                                    ::  BC_Enc_Mass,    &
                                                                        Shift_Value




DO ui = 1,Num_Eqs


    IF (    Inner_BC_Type(ui) .EQ. "N"   ) THEN

        BC_Enc_Mass = Inner_BC_Values(ui)

         !                                                   !
        !!   We are asssuming a uniform value on the shell   !!
        !!   so only the L,M = 0 vector is effected          !!
         !                                                   !
        IF (    L .EQ. 0    ) THEN

            !                           !
            !   Calulcate Shift Value   !
            !                           !

            Shift_Value = sqrt(4.0_idp*pi)*BC_Enc_Mass
            WORK_VEC(0) = WORK_VEC(0) - Shift_Value

        END IF
    END IF




    IF (    Outer_BC_Type(ui) .EQ. "N"   ) THEN


        BC_Enc_Mass = Outer_BC_Values(ui)


         !                                                   !
        !!   We are asssuming a uniform value on the shell   !!
        !!   so only the L,M = 0 vector is effected.         !!
         !                                                   !
        IF (    L .EQ. 0    ) THEN

            !                           !
            !   Calulcate Shift Value   !
            !                           !

            Shift_Value = sqrt(4.0_idp*pi)*BC_Enc_Mass


            WORK_VEC(N-1) = WORK_VEC(N-1) + Shift_Value

        END IF
    END IF

END DO ! ui







END SUBROUTINE NEUMANN_BC_CCS













 !+301+############################################################################!
!                                                                                   !
!       Dirichlet_BC_CHOL                                                           !
!                                                                                   !
!-----------------------------------------------------------------------------------!
!                                                                                   !
!   Input                                                                           !
!                                                                                   !
 !#################################################################################!
SUBROUTINE DIRICHLET_BC_CHOL(N, NNZ, L, M, COL_PTR, ROW_IND, WORK_VEC,ui)


INTEGER,                                INTENT(IN)                  ::  N, NNZ, L, M

INTEGER, DIMENSION(0:N),                INTENT(IN)                  ::  COL_PTR
INTEGER, DIMENSION(0:NNZ-1),            INTENT(IN)                  ::  ROW_IND

REAL(idp), DIMENSION(0:N - 1),          INTENT(INOUT)               ::  WORK_VEC
INTEGER,                                INTENT(IN)                  ::  ui



INTEGER                                                             :: i, shift
REAL(idp)                                                           :: BC_Value






IF (Inner_BC_Type(ui) == "D") THEN





    BC_Value = sqrt(4.0_idp*pi)*Inner_BC_Values(ui)



    !!! MODIFY SOURCE VECTOR !!!


    WORK_VEC(0) = BC_Value

    DO i = 1,DEGREE

        WORK_VEC(i) = WORK_VEC(i) - dMA_First_Col_Storage(i,L)*BC_Value

    END DO
    !!! MODIFY MATRIX !!!

    !!! ALREADY DONE IN CHOLESKY FACTORIZATION !!


END IF




shift = 0
IF (NUM_R_ELEMENTS .EQ. 1 ) THEN
    shift = 1
END IF

IF (Outer_BC_Type(ui)  == "D") THEN


    IF ( ( L == 0 ) .AND. ( M == 0 )  ) THEN        
        BC_Value = sqrt(4.0_idp*pi)*Outer_BC_Values(ui)
    ELSE
        BC_Value = 0.0_idp
    END IF

    !!! MODIFY SRC VECTOR !!!
    WORK_VEC(NUM_R_NODES - 1) = BC_Value


    DO i = DEGREE-shift,1,-1

        WORK_VEC(NUM_R_NODES - 1 - i) = WORK_VEC(NUM_R_NODES - 1 - i)                           &
                                        - dMA_Last_Col_Storage(i,L)*BC_Value



    END DO

    !!! MODIFY MATRIX !!!

    !!! ALREADY DONE IN CHOLESKY FACTORIZATION !!!




END IF



END SUBROUTINE DIRICHLET_BC_CHOL













!+101+############################################################################!
!                                                                                   !
!       Dirichlet_BC                                                                !
!                                                                                   !
!-----------------------------------------------------------------------------------!
!                                                                                   !
!   Input                                                                           !
!                                                                                   !
 !#################################################################################!
SUBROUTINE DIRICHLET_BC_Beta(WORK_MAT, WORK_VEC)

REAL(idp), DIMENSION(1:iVB_Prob_Dim),                   INTENT(INOUT)   :: WORK_VEC

REAL(idp), DIMENSION(1:iVB_Prob_Dim,1:iVB_Prob_Dim),    INTENT(INOUT)   :: WORK_MAT




INTEGER                 :: i, shift, uj, ui, m, l, Here


REAL(idp)                                                               :: BC_Value


DO ui = 3,1,-1
    uj = ui + 2

    IF (Inner_BC_Type(ui) == "D") THEN

        BC_Value =  2.0_idp*sqrt(pi)*Inner_BC_Values(uj)


        WORK_VEC(1) = BC_Value
        DO i = 2,DEGREE+1
            WORK_VEC(i) = WORK_VEC(i) - WORK_MAT(1,i)*BC_Value
        END DO


        WORK_MAT(1,:) = 0.0_idp
        WORK_MAT(:,1) = 0.0_idp
        WORK_MAT(1,1) = 1.0_idp

    END IF










    shift = 0
    IF (NUM_R_ELEMENTS .EQ. 1 ) THEN
        shift = 1
    END IF

    IF (Outer_BC_Type(ui)  == "D") THEN

        DO l = 0,L_LIMIT
            DO m = -l,l


                IF ( ( L == 0 ) .AND. ( M == 0 )  ) THEN
                    BC_Value = 2.0_idp*sqrt(pi)*Outer_BC_Values(uj)
                ELSE
                    BC_Value = 0.0_idp
                END IF

                Here = FP_Beta_Array_Map(Num_R_Elements-1,Degree,ui,l,m)


                DO i = 0,iVB_Elem_Prob_Dim-1
                    Work_Vec(iVB_Prob_Dim-i) = Work_Vec(iVB_Prob_Dim-i)       &
                                              - Work_Mat(iVB_Prob_Dim-i,Here)*BC_Value
    
!                    PRINT*,iVB_Prob_Dim-i,                         &
!                            Work_Vec(iVB_Prob_Dim-i),              &
!                            BC_Value,                               &
!                            Work_Mat(iVB_Prob_Dim-i,Here)

                END DO

                ! Set the BC Value in the Coefficient Vector
                ! l,m = 0  is set to the BC
                ! l,m != 0 is set to zero.
                WORK_VEC( Here ) = BC_Value

                ! Transform the Stiffness Matrix
                WORK_MAT(Here,:)    = 0.0_idp
                WORK_MAT(:,Here)    = 0.0_idp
                WORK_MAT(Here,Here) = 1.0_idp

            END DO !
        END DO ! l Loop

    END IF


END DO ! i Loop


END SUBROUTINE DIRICHLET_BC_Beta












 !+201+############################################################################!
!                                                                                   !
!       Dirichlet_BC_Beta_Banded                                                     !
!                                                                                   !
!-----------------------------------------------------------------------------------!
!                                                                                   !
!   Input                                                                           !
!                                                                                   !
 !#################################################################################!
SUBROUTINE DIRICHLET_BC_Beta_Banded( N, WORK_VEC )


INTEGER,                            INTENT(IN)                  ::  N
REAL(idp), DIMENSION(1:N),          INTENT(INOUT)               ::  WORK_VEC

REAL(idp)                                                       :: BC_Value

INTEGER                                                         :: ui
INTEGER                                                         :: lm, d, Row
INTEGER                                                         :: Shift




DO ui = 3,5
    
    IF (Inner_BC_Type(ui) == "D") THEN

        BC_Value = sqrt(4.0_idp*pi)*Inner_BC_Values(ui)

        DO d =  1,DEGREE

            Row = FP_Beta_Array_Map(Num_R_Elements-1, d, ui-2, 0, 0)
            Work_Vec(Row) = Work_Vec(Row) - dMB_First_Col_Storage(0,d,ui-2)*BC_Value
        
        END DO  ! d


        Row = (ui - 3) * LM_Length + 1
        WORK_VEC(Row) = BC_Value

        !!! MODIFY MATRIX !!!

        !!! ALREADY DONE IN Initalization !!

    END IF



    shift = 0
    IF (NUM_R_ELEMENTS .EQ. 1 ) THEN
        shift = 1
    END IF

    

    IF (Outer_BC_Type(ui)  == "D") THEN
    

        DO lm = 1,LM_Length

            IF ( lm == 1 ) THEN
                BC_Value = sqrt(4.0_idp*pi)*Outer_BC_Values(ui)
            ELSE
                BC_Value = 0.0_idp
            END IF
            


            DO d = DEGREE-shift,0,-1
                
                Row = FP_Beta_Array_Map(Num_R_Elements-1,d,ui-2,lm)
                
                Work_Vec(Row) = Work_Vec(Row) - dMB_Last_Col_Storage(lm,d,ui-2)*BC_Value

            END DO  ! d

        
            Row = FP_Beta_Array_Map(Num_R_Elements-1,Degree,ui-2,lm)
            WORK_VEC(Row) = BC_Value

        END DO  ! lm

        !!! MODIFY MATRIX !!!

        !!! ALREADY DONE IN INITIALIZATION !!!


    END IF

END DO ! ui Loop


END SUBROUTINE DIRICHLET_BC_Beta_Banded









END MODULE Matrix_Boundary_Condition_Routines
