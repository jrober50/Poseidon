   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Poisson_Matrix_Routines                                               !##!
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

USE Poseidon_Parameters, &
            ONLY :  Degree,             &
                    L_Limit

USE Variables_Mesh, &
            ONLY :  Num_R_Elements,     &
                    Num_T_Elements,     &
                    Num_P_Elements,     &
                    rlocs

USE Variables_Poisson, &
            ONLY :  STF_NNZ,            &
                    STF_Mat_Integrals,  &
                    STF_Elem_Val,       &
                    STF_Row_Ind,        &
                    STF_Col_Ptr


USE Functions_Math, &
            ONLY :  Lagrange_Poly,      &
                    Lagrange_Poly_Deriv

USE Functions_Quadrature, &
            ONLY :  Initialize_LGL_Quadrature

IMPLICIT NONE


CONTAINS



 !+103+################################################################!
!                                                                       !
!   Initialize_Stiffness_Matrix                                         !
!                                                                       !
 !#####################################################################!
SUBROUTINE Initialize_Stiffness_Matrix()


INTEGER                                             ::  i, j, k, OrdofQuad, HERE


REAL(KIND = idp), DIMENSION(0:DEGREE)               ::  LagP, DLagP, xlocP, weightP
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)         ::  xlocQ, weightQ



STF_MAT_Integrals = 0.0_idp


OrdofQuad = DEGREE + 10

ALLOCATE(xlocQ(0:OrdofQuad), weightQ(0:OrdofQuad))

!!! Generate nodes for Lagrange Polynomials !!!
CALL Initialize_LGL_Quadrature(DEGREE, xlocP, weightP)

!!! Generate nodes and weights for Quadrature !!!
CALL Initialize_LGL_Quadrature(OrdofQuad, xlocQ, weightQ)


DO k = 0,OrdofQuad

    !!! Generate Lagrange Polynomials for xlocQ(k) using xlocP for nodes !!!
    LagP = Lagrange_Poly(xlocQ(k), DEGREE, xlocP)

    !!! Generate First Derivatives of Lagrange Polynomials for xlocQ(k) using xlocP for nodes !!!
    DLagP = Lagrange_Poly_Deriv(xlocQ(k), DEGREE, xlocP)

    DO j = 0,DEGREE
    DO i = 0,DEGREE


        !!! Integral of LagP*LagP in x-space. Mapped to r-space in GENERATE_STF_MAT
        STF_MAT_Integrals(0,i,j) = STF_MAT_Integrals(0,i,j) + weightQ(k)*LagP(i)*LagP(j)


        !!! Int of r^2 * DLagP * DLagP dr in x-space.
        !!! Mapping to x-space creates three integral
        !!! 1) Int DLagP*DLagP
        STF_MAT_Integrals(1,i,j) = STF_MAT_Integrals(1,i,j) + weightQ(k)*DLagP(i)*DLagP(j)

        !!! 2) Int x*DLagP*DLagP
        STF_MAT_Integrals(2,i,j) = STF_MAT_Integrals(2,i,j) + weightQ(k)*xlocQ(k)*DLagP(i)*DLagP(j)

        !!! 3) Int x^2*DLagP*DLagP
        STF_MAT_Integrals(3,i,j) = STF_MAT_Integrals(3,i,j) + weightQ(k)*xlocQ(k)*xlocQ(k)*DLagP(i)*DLagP(j)

        !!! In GENERATE_STF_MAT these integrals are mapped back to r-space for given a given range of r.

    END DO ! i
    END DO ! j
END DO ! k





!
!   Because the structure of the stiffness matrix is known, we can
!       initialize the ROW_IND, and COL_PTR arrays now, and will
!       fill ELEM_VAL array when GENERATE_STF_MAT_CCS is called.
!

  !                             !
 !!                             !!
!!!    COL_PTR INITIALIZATION   !!!
 !!                             !!
  !                             !

STF_COL_PTR(0) = 0
STF_COL_PTR(1) = STF_COL_PTR(0) + (DEGREE+1)
HERE = 2

DO i = 1,NUM_R_ELEMENTS-1

    DO j = 1,DEGREE - 1

        STF_COL_PTR(Here) = STF_COL_PTR(Here - 1) +  (DEGREE + 1)

        Here = Here + 1

    END DO

    STF_COL_PTR(Here) = STF_COL_PTR(Here - 1) + (2*DEGREE + 1)

    Here = Here + 1

END DO

DO i = 1, DEGREE

    STF_COL_PTR(Here) = STF_COL_PTR(Here - 1) + (DEGREE + 1)

    Here = Here + 1

END DO



  !                             !
 !!                             !!
!!!    ROW_IND INITIALIZATION   !!!
 !!                             !!
  !                             !

Here = 0
DO i = 0, NUM_R_ELEMENTS - 1

DO j = 0,DEGREE - 1
DO k = 0, DEGREE

    STF_ROW_IND(Here) = i*DEGREE + k
    Here = Here + 1

END DO  ! k
END DO  ! j

DO k = 0,DEGREE - 1

    STF_ROW_IND(Here) = i*DEGREE + k
    Here = Here + 1

END DO  ! k

END DO  ! i


STF_ROW_IND(Here) = DEGREE * NUM_R_ELEMENTS



CALL Generate_Stiffness_Matrix()



END SUBROUTINE Initialize_Stiffness_Matrix








!+202+######################################################################!
!                                                                           !
!   Generate_Stiffness_Matrix -  Generates the stiffness matrix in the  !
!                                      Compress Column Storage (CCS) format !
!                                                                           !
!###########################################################################!
SUBROUTINE Generate_Stiffness_Matrix()


INTEGER                     ::  i, j, h, k, Here

REAL(KIND = idp)                                ::  deltar, rplusr



STF_ELEM_VAL = 0.0_idp

Here = 0
DO j = 0,L_LIMIT

    Here = 0
    DO i = 0,NUM_R_ELEMENTS-1

        deltar = rlocs(i+1) - rlocs(i)
        rplusr = rlocs(i+1) + rlocs(i)

        DO h = 0, DEGREE
        DO k = 0, DEGREE

                STF_ELEM_VAL(Here, j)= STF_ELEM_VAL(Here, j) &
                                            + j*(j+1)*(deltar/2.0_idp)*STF_MAT_Integrals(0,h,k) &
                                            + (     (rplusr*rplusr)/(2.0_idp*deltar)            &
                                                        * STF_MAT_Integrals(1,h,k)              &
                                                    + rplusr*STF_MAT_Integrals(2,h,k)           &
                                                    + (deltar/2.0_idp)*STF_MAT_Integrals(3,h,k) )

                Here = Here + 1


        END DO
        END DO

        Here = Here - 1     !!! Take one step back due to overlap of first value
                            !!! in next element with last in current element

    END DO

END DO





END SUBROUTINE Generate_Stiffness_Matrix








END MODULE Poisson_Matrix_Routines
