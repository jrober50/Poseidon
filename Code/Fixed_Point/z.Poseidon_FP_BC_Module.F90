   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_FP_BC_Module                                                        !##!
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
USE Poseidon_Constants_Module, &
                        ONLY :  idp,                &
                                pi


USE Poseidon_Parameters, &
                                ONLY :  DOMAIN_DIM,                 &
                                        DEGREE,                     &
                                        L_LIMIT,                    &
                                        NUM_CFA_EQs

USE DRIVER_Parameters, &
                        ONLY :  Potential_Solution


USE Poseidon_FP_Variables_Module, &
                        ONLY :  First_Column_Storage,               &
                                Last_Column_Storage,                &
                                CFA_EQ_Map




USE Poseidon_Variables_Module, &
                        ONLY :  NUM_R_ELEMENTS,             &
                                NUM_R_NODES,                &
                                INNER_CFA_BC_VALUES,        &
                                OUTER_CFA_BC_VALUES,        &
                                INNER_CFA_BC_TYPE,          &
                                OUTER_CFA_BC_TYPE,          &
                                R_INNER, R_OUTER


USE Poseidon_Math_Functions_Module, &
                        ONLY :  Spherical_Harmonic

USE Poseidon_Quadrature_Module, &
                        ONLY :  Initialize_LG_Quadrature,       &
                                Initialize_LGL_Quadrature

USE Poseidon_Mapping_Functions_Module, &
                                ONLY :  Map_From_X_Space


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

        !!! ALREADY DONE IN CHOLESKY FACTORIZATION !!!




    END IF









    shift = 0
    IF (NUM_R_ELEMENTS .EQ. 1 ) THEN
        shift = 1
    END IF

    IF (OUTER_CFA_BC_TYPE(ui)  == "D") THEN



        BC_Value = sqrt(4.0_idp*pi)*OUTER_CFA_BC_VALUES(ui)

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





















 !+401+####################################################################!
!                                                                           !
!       BC_Integral - Boundary Condition Integral
!                                                                           !
 !#########################################################################!
FUNCTION BC_Integral(R, L, M)


COMPLEX(KIND = idp)                                :: BC_Integral



REAL(KIND = idp), INTENT(IN)                    :: R
INTEGER, INTENT(IN)                             :: L, M

INTEGER                                         :: te, td, pe, pd


INTEGER                                         ::  T_Degree, P_Degree

INTEGER                                         ::  T_SUB_ELEMENTS, P_SUB_ELEMENTS


REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)     ::  P_xlocs, P_locs, P_weights, &
                                                    T_xlocs, T_locs, T_weights, &
                                                    Sub_tlocs, Sub_plocs

REAL(KIND = idp)                                ::  deltasubt, deltasubp, deltat, deltap


COMPLEX(KIND = idp)                             ::  TMP_INTEGRAL_VALUE


T_SUB_ELEMENTS = 32
P_SUB_ELEMENTS = 32

T_Degree = 5
P_Degree = 5





ALLOCATE(P_xlocs(0:P_Degree), P_locs(0:P_Degree), P_weights(0:P_Degree))
ALLOCATE(T_xlocs(0:T_Degree), T_locs(0:T_Degree), T_weights(0:T_Degree))
ALLOCATE(Sub_tlocs(0:T_SUB_ELEMENTS),Sub_plocs(0:P_SUB_ELEMENTS))

CALL Initialize_LGL_Quadrature(P_Degree, P_xlocs, P_weights)
CALL Initialize_LGL_Quadrature(T_Degree, T_xlocs, T_weights)

deltasubt = (pi)/REAL(T_SUB_ELEMENTS)
deltasubp = (2*pi)/REAL(P_SUB_ELEMENTS)

DO te = 0, T_SUB_ELEMENTS
    Sub_tlocs(te) = 0.0_idp + te*deltasubt
END DO



DO pe = 0,P_SUB_ELEMENTS
    Sub_plocs(pe) = 0.0_idp + pe*deltasubp
END DO




TMP_INTEGRAL_VALUE = 0.0_idp
DO te = 0, T_SUB_ELEMENTS - 1

    T_locs = Map_From_X_Space(Sub_tlocs(te), Sub_tlocs(te+1), T_xlocs)
    deltat = Sub_tlocs(te+1)-Sub_tlocs(te)

    DO pe = 0, P_SUB_ELEMENTS - 1

        P_locs = Map_From_X_Space(Sub_plocs(pe), Sub_plocs(pe+1), P_xlocs)
        deltap = Sub_plocs(pe+1) - Sub_plocs(pe)

        DO td = 0, T_Degree

            DO pd = 0, P_Degree

                TMP_INTEGRAL_VALUE = TMP_INTEGRAL_VALUE                                                     &
                                    + T_weights(td)*P_weights(pd)                                           &
                                    * (-1.0_idp)**M*Spherical_Harmonic(L, -M, T_locs(td), P_locs(pd))   &
                                    * sin(T_locs(td))*Potential_Solution(R, T_locs(td), P_locs(pd))          &
                                    * deltat/2.0_idp*deltap/2.0_idp





            END DO

        END DO

    END DO

END DO




BC_Integral = TMP_INTEGRAL_VALUE




END FUNCTION BC_Integral






END MODULE Poseidon_FP_BC_Module
