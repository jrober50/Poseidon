   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE FP_Functions_BC                                                              !##!
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
            ONLY :  DOMAIN_DIM,                 &
                    DEGREE,                     &
                    L_LIMIT,                    &
                    NUM_CFA_EQs

USE Variables_FP, &
            ONLY :  First_Column_Storage,       &
                    Last_Column_Storage,        &
                    CFA_EQ_Map

USE Variables_Mesh, &
            ONLY :  Num_R_Elements,             &
                    R_INNER,                    &
                    R_OUTER

USE Variables_Derived, &
            ONLY :  Num_R_Nodes,                &
                    Beta_Prob_Dim,              &
                    Beta_Elem_Prob_Dim,         &
                    LM_Length

USE Variables_BC, &
            ONLY :  INNER_CFA_BC_VALUES,        &
                    OUTER_CFA_BC_VALUES,        &
                    INNER_CFA_BC_TYPE,          &
                    OUTER_CFA_BC_TYPE



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

INTEGER, INTENT(IN)                                                                     :: L, M, ui

!REAL(KIND = idp), DIMENSION(0:NUM_R_NODES - 1), INTENT(INOUT)                           :: WORK_VEC
COMPLEX(KIND = idp), DIMENSION(1:NUM_R_NODES), INTENT(INOUT)                            :: WORK_VEC

COMPLEX(KIND = idp), DIMENSION(1:NUM_R_NODES,1:NUM_R_NODES), INTENT(INOUT)         :: WORK_MAT




INTEGER                 :: i, shift, uj


COMPLEX(KIND = idp)                                                         :: BC_Value


uj = CFA_EQ_Map(ui)
IF (INNER_CFA_BC_TYPE(uj) == "D") THEN

    BC_Value =  2.0_idp*sqrt(pi)*INNER_CFA_BC_VALUES(uj)
    

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

IF (OUTER_CFA_BC_TYPE(uj)  == "D") THEN

    IF ( ( L == 0 ) .AND. ( M == 0 )  ) THEN
        BC_Value = 2.0_idp*sqrt(pi)*OUTER_CFA_BC_VALUES(uj)
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

INTEGER,                                            INTENT(IN)                      ::  L_VALUE
COMPLEX(KIND = idp), DIMENSION(0:NUM_R_NODES - 1),  INTENT(INOUT)                   ::  WORK_VEC


INTEGER                 :: i,j, ui






DO ui = 1,NUM_CFA_EQs
    IF (    L_VALUE == 0    )   THEN


        IF (INNER_CFA_BC_TYPE(ui) == "N") THEN

            WORK_VEC(0) = WORK_VEC(0) - sqrt(4.0_idp*pi)*INNER_CFA_BC_VALUES(ui)

        END IF


        IF (OUTER_CFA_BC_TYPE(ui) == "N") THEN

            WORK_VEC(NUM_R_NODES-1) = WORK_VEC(NUM_R_NODES-1) + R_OUTER*R_OUTER*OUTER_CFA_BC_VALUES(ui)

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


INTEGER,                                        INTENT(IN)                  ::  N, NNZ, L, M

INTEGER, DIMENSION(0:N),                        INTENT(IN)                  ::  COL_PTR
INTEGER, DIMENSION(0:NNZ-1),                    INTENT(IN)                  ::  ROW_IND

COMPLEX(KIND = idp), DIMENSION(0:N - 1),        INTENT(INOUT)               ::  WORK_VEC


COMPLEX(KIND = idp), DIMENSION(0:NNZ-1),           INTENT(INOUT)               ::  ELEM_VAL


INTEGER                                                                     :: i, shift, ui
COMPLEX(KIND = idp)                                                         :: BC_Value





DO ui = 1,NUM_CFA_EQs
    IF (INNER_CFA_BC_TYPE(ui) == "D") THEN



        BC_Value = sqrt(4.0_idp*pi)*INNER_CFA_BC_VALUES(ui)



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

    IF (OUTER_CFA_BC_TYPE(ui)  == "D") THEN



        BC_Value = sqrt(4.0_idp*pi)*OUTER_CFA_BC_VALUES(ui)

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


INTEGER,                                        INTENT(IN)                  ::  N, NNZ, L, M

INTEGER, DIMENSION(0:N),                        INTENT(IN)                  ::  COL_PTR
INTEGER, DIMENSION(0:NNZ-1),                    INTENT(IN)                  ::  ROW_IND

COMPLEX(KIND = idp), DIMENSION(0:N - 1),        INTENT(INOUT)               ::  WORK_VEC


COMPLEX(KIND = idp), DIMENSION(0:NNZ-1),           INTENT(INOUT)               ::  ELEM_VAL


INTEGER                                                                     ::  i, shift, ui

REAL(KIND = idp)                                                            ::  BC_Enc_Mass,    &
                                                                                Shift_Value




DO ui = 1,NUM_CFA_EQs


    IF (    INNER_CFA_BC_TYPE(ui) .EQ. "N"   ) THEN

        BC_Enc_Mass = INNER_CFA_BC_VALUES(ui)

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




    IF (    OUTER_CFA_BC_TYPE(ui) .EQ. "N"   ) THEN


        BC_Enc_Mass = OUTER_CFA_BC_VALUES(ui)


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
SUBROUTINE DIRICHLET_BC_CHOL(N, NNZ, L, M, COL_PTR, ROW_IND, WORK_VEC)


INTEGER,                                        INTENT(IN)                  ::  N, NNZ, L, M

INTEGER, DIMENSION(0:N),                        INTENT(IN)                  ::  COL_PTR
INTEGER, DIMENSION(0:NNZ-1),                    INTENT(IN)                  ::  ROW_IND

COMPLEX(KIND = idp), DIMENSION(0:N - 1),        INTENT(INOUT)               ::  WORK_VEC




INTEGER                                                                     :: i, shift, ui
COMPLEX(KIND = idp)                                                         :: BC_Value




DO ui = 1,NUM_CFA_EQs

    IF (INNER_CFA_BC_TYPE(ui) == "D") THEN





        BC_Value = sqrt(4.0_idp*pi)*INNER_CFA_BC_VALUES(ui)



        !!! MODIFY SOURCE VECTOR !!!


        WORK_VEC(0) = BC_Value

        DO i = 1,DEGREE

            WORK_VEC(i) = WORK_VEC(i) - First_Column_Storage(i ,L,ui)*BC_Value

        END DO
        !!! MODIFY MATRIX !!!

        !!! ALREADY DONE IN CHOLESKY FACTORIZATION !!


    END IF




    shift = 0
    IF (NUM_R_ELEMENTS .EQ. 1 ) THEN
        shift = 1
    END IF

    IF (OUTER_CFA_BC_TYPE(ui)  == "D") THEN


        IF ( ( L == 0 ) .AND. ( M == 0 )  ) THEN
            BC_Value = sqrt(4.0_idp*pi)*OUTER_CFA_BC_VALUES(ui)
        ELSE
            BC_Value = 0.0_idp
        END IF

        !!! MODIFY SRC VECTOR !!!
        WORK_VEC(NUM_R_NODES - 1) = BC_Value


        DO i = DEGREE-shift,1,-1


            WORK_VEC(NUM_R_NODES - 1 - i) = WORK_VEC(NUM_R_NODES - 1 - i)                           &
                                            - Last_Column_Storage(i,L,ui)*BC_Value



        END DO





        !!! MODIFY MATRIX !!!

        !!! ALREADY DONE IN CHOLESKY FACTORIZATION !!!




    END IF

END DO



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

COMPLEX(KIND = idp), DIMENSION(1:Beta_Prob_Dim), INTENT(INOUT)                  :: WORK_VEC

COMPLEX(KIND = idp), DIMENSION(1:Beta_Prob_Dim,1:Beta_Prob_Dim), INTENT(INOUT)  :: WORK_MAT




INTEGER                 :: i, shift, uj, ui, m, l, d, Here


COMPLEX(KIND = idp)                                                         :: BC_Value


DO ui = 3,1,-1
    uj = ui + 2

    IF (INNER_CFA_BC_TYPE(ui) == "D") THEN

        BC_Value =  2.0_idp*sqrt(pi)*INNER_CFA_BC_VALUES(uj)


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

    IF (OUTER_CFA_BC_TYPE(ui)  == "D") THEN

        DO l = 0,L_LIMIT
            DO m = -l,l


                IF ( ( L == 0 ) .AND. ( M == 0 )  ) THEN
                    BC_Value = 2.0_idp*sqrt(pi)*OUTER_CFA_BC_VALUES(uj)
                ELSE
                    BC_Value = 0.0_idp
                END IF


                Here = (Num_R_Nodes-1) * 3 * LM_Length  &
                        + (ui - 1) * LM_Length          &
                        + l*(l+1) + m + 1

                DO i = 0,Beta_Elem_Prob_Dim-1
                    Work_Vec(Beta_Prob_Dim-i) = Work_Vec(Beta_Prob_Dim-i)       &
                                              - Work_Mat(Beta_Prob_Dim-i,Here)*BC_Value
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







END MODULE FP_Functions_BC
