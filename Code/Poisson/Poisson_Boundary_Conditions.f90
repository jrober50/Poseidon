   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poisson_Boundary_Conditions_Module                                           !##!
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
            ONLY :  Degree

USE Variables_Mesh, &
            ONLY :  Num_R_Elements,         &
                    R_Inner,                &
                    R_Outer

USE Variables_Derived, &
            ONLY :  LM_LENGTH,              &
                    Num_R_Nodes

USE Variables_Poisson, &
            ONLY :  First_Column_Storage,   &
                    Last_Column_Storage

USE Variables_BC, &
            ONLY :  INNER_Poisson_BC_Type,  &
                    INNER_Poisson_BC_Value, &
                    OUTER_Poisson_BC_Type,  &
                    OUTER_Poisson_BC_Value



IMPLICIT NONE

!*F&S*==========================================!
!                                               !
!           Functions & Subroutines             !
!                                               !
!===============================================!
CONTAINS




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




INTEGER                                                                     :: i, shift
COMPLEX(KIND = idp)                                                         :: BC_Value






IF (INNER_Poisson_BC_Type == "D") THEN


          !                                                                                       !
         !!   For a uniform boundary condition on a sphere the only spherical harmonic expansion  !!
        !!!   coefficient to be effected will be the l,m = 0,0 coefficient due to the symmetry    !!!
        !!!   of the condition, and the orthogonality of the spherical harmonics functions.       !!!
         !!   The integral over the theta, and phi also produces a sqrt(4*pi) scalling factor.    !!
          !                                                                                       !
    IF (    L == 0  )   THEN


        BC_Value = sqrt(4.0_idp*pi)*INNER_Poisson_BC_Value

    ELSE

        BC_Value = 0.0_idp

    END IF






    !!! MODIFY SOURCE VECTOR !!!


    WORK_VEC(0) = BC_Value

    DO i = 1,DEGREE

        WORK_VEC(i) = WORK_VEC(i) - First_Column_Storage(i , L)*BC_Value

    END DO



    !!! MODIFY MATRIX !!!

    !!! ALREADY DONE IN CHOLESKY FACTORIZATION !!!




END IF









shift = 0
IF (NUM_R_ELEMENTS .EQ. 1 ) THEN
    shift = 1
END IF

IF (Outer_Poisson_BC_Type  == "D") THEN

    PRINT*,"In Outer_BC Dirichlet"

          !                                                                                       !
         !!   For a uniform boundary condition on a sphere the only spherical harmonic expansion  !!
        !!!   coefficient to be effected will be the l,m = 0,0 coefficient due to the symmetry    !!!
        !!!   of the condition, and the orthogonality of the spherical harmonics functions.       !!!
         !!   The integral over the theta, and phi also produces a sqrt(4*pi) scalling factor.    !!
          !                                                                                       !
    IF (    L == 0  )   THEN
        
        BC_Value = sqrt(4.0_idp*pi)*Outer_Poisson_BC_Value

    ELSE

        BC_Value = 0.0_idp

    END IF




    !!! MODIFY SRC VECTOR !!!
    WORK_VEC(NUM_R_NODES - 1) = BC_Value


    DO i = DEGREE-shift,1,-1


        WORK_VEC(NUM_R_NODES - 1 - i) = WORK_VEC(NUM_R_NODES - 1 - i)                           &
                                        - Last_Column_Storage(i,L)*BC_Value


    END DO





    !!! MODIFY MATRIX !!!

    !!! ALREADY DONE IN CHOLESKY FACTORIZATION !!!




END IF





END SUBROUTINE DIRICHLET_BC_CHOL















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


REAL(KIND = idp), DIMENSION(0:NNZ-1),           INTENT(INOUT)               ::  ELEM_VAL

REAL(KIND = idp)                                                            ::  BC_Enc_Mass,    &
                                                                                Shift_Value







IF (  INNER_Poisson_BC_Type .EQ. "N"   ) THEN



    BC_Enc_Mass = Outer_Poisson_BC_Value



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










IF (  OUTER_Poisson_BC_Type .EQ. "N"   ) THEN


    BC_Enc_Mass = Outer_Poisson_BC_Value


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









END SUBROUTINE NEUMANN_BC_CCS
































END MODULE Poisson_Boundary_Conditions_Module

