   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE FP_Functions_Coeffs_Module                                            !##!
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

USE Parameters_Variable_Indices, &
            ONLY :  iU_CF,              &
                    iU_LF,              &
                    iU_S1,              &
                    iU_S2,              &
                    iU_S3,              &
                    iVB_S

USE Variables_Derived, &
            ONLY :  LM_Length,          &
                    Var_Dim,            &
                    Prob_Dim,           &
                    Num_R_Nodes

USE Variables_FP,  &
            ONLY :  FP_Coeff_Vector_A,  &
                    FP_Coeff_Vector_B


IMPLICIT NONE


CONTAINS
!+101+##########################################################################!
!                                                                               !
!                  Coeff_To_Vector                                  !
!                                                                               !
!###############################################################################!
SUBROUTINE Coeff_To_Vector( Vector )

COMPLEX(idp), DIMENSION(Prob_Dim), INTENT(INOUT) :: Vector

INTEGER                                         :: ui

DO ui = iU_CF,iU_LF

    CALL CoeffA_To_Vector(Vector, ui)

END DO

DO ui = iU_S1,iU_S3

    CALL CoeffB_To_Vector(Vector, ui, iVB_S)

END DO

END SUBROUTINE Coeff_To_Vector




!+101+##########################################################################!
!                                                                               !
!                  Vector_To_Coeff                                  !
!                                                                               !
!###############################################################################!
SUBROUTINE Vector_To_Coeff( Vector )

COMPLEX(idp), DIMENSION(Prob_Dim), INTENT(INOUT) :: Vector

INTEGER                                         :: ui

DO ui = iU_CF,iU_LF

    CALL Vector_To_CoeffA(Vector, ui)

END DO

DO ui = iU_S1,iU_S3

    CALL Vector_To_CoeffB(Vector, ui, iVB_S)

END DO

END SUBROUTINE Vector_To_Coeff




!+101+##########################################################################!
!                                                                               !
!                  CoeffA_To_Vector                                  !
!                                                                               !
!###############################################################################!
SUBROUTINE CoeffA_To_Vector( Vector, iU )

COMPLEX(idp), DIMENSION(Prob_Dim), INTENT(INOUT) :: Vector
INTEGER                         , INTENT(IN)    :: iU

INTEGER                                         :: LM_loc, Here, There


DO LM_loc = 1,LM_Length
    here = (iU-1)*Var_Dim + (lm_loc-1)*Num_R_Nodes + 1
    there = (iU-1)*Var_Dim + lm_loc*Num_R_Nodes
    Vector(here:there) = FP_Coeff_Vector_A(:,lm_loc,iU)
END DO


END SUBROUTINE CoeffA_To_Vector


!+101+##########################################################################!
!                                                                               !
!                  CoeffB_To_Vector                                  !
!                                                                               !
!###############################################################################!
SUBROUTINE CoeffB_To_Vector( Vector, iU, iVB )

COMPLEX(idp), DIMENSION(Prob_Dim), INTENT(INOUT) :: Vector
INTEGER                          , INTENT(IN)    :: iU, iVB

INTEGER                                         :: HereA, ThereA
INTEGER                                         :: HereB, ThereB

INTEGER                                     ::  Offset

Offset = 3*iVB

HereA  = (iU-1)*Var_Dim + 1
ThereA = iU * Var_Dim

HereB  = (iU-Offset)*Var_Dim + 1
ThereB = (iU-Offset + 1)*Var_Dim


Vector(HereA:ThereA) = FP_Coeff_Vector_B(HereB:ThereB,iVB)


END SUBROUTINE CoeffB_To_Vector




!+202+##########################################################################!
!                                                                               !
!                    Vector_To_CoeffA                                  !
!                                                                               !
!###############################################################################!
SUBROUTINE Vector_To_CoeffA( Vector, iU )

COMPLEX(idp), DIMENSION(Prob_Dim), INTENT(IN)        :: Vector
INTEGER                         , INTENT(IN)        :: iU

INTEGER                                             :: LM_loc, Here, There

DO LM_loc = 1,LM_Length
    here = (iU-1)*Var_Dim + (lm_loc-1)*Num_R_Nodes + 1
    there = (iU-1)*Var_Dim + lm_loc*Num_R_Nodes
    FP_Coeff_Vector_A(:,lm_loc,iU) = Vector(here:There)
END DO


END SUBROUTINE Vector_To_CoeffA


!+101+##########################################################################!
!                                                                               !
!                 Vector_To_CoeffB                                  !
!                                                                               !
!###############################################################################!
SUBROUTINE Vector_To_CoeffB( Vector, iU, iVB )

COMPLEX(idp), DIMENSION(Prob_Dim), INTENT(IN)    :: Vector
INTEGER                         , INTENT(IN)    :: iU, iVB

INTEGER                                         :: HereA, ThereA
INTEGER                                         :: HereB, ThereB

INTEGER                                     ::  Offset

Offset = 3*iVB

HereA  = (iU-1)*Var_Dim + 1
ThereA = iU * Var_Dim

HereB  = (iU-Offset)*Var_Dim + 1
ThereB = (iU-Offset + 1)*Var_Dim

FP_Coeff_Vector_B(HereB:ThereB,iVB) = Vector(HereA:ThereA)


END SUBROUTINE Vector_To_CoeffB





END MODULE FP_Functions_Coeffs_Module
