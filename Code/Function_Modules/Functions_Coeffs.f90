   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Functions_Coeffs_Module                                               !##!
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

USE Variables_Vectors,  &
            ONLY :  cVA_Coeff_Vector,  &
                    cVB_Coeff_Vector


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

    CALL Coeff_To_Vector_TypeA(Vector, ui)

END DO

DO ui = iU_S1,iU_S3

    CALL Coeff_To_Vector_TypeB(Vector, ui, iVB_S)

END DO

END SUBROUTINE Coeff_To_Vector




!+102+##########################################################################!
!                                                                               !
!                  Vector_To_Coeff                                  !
!                                                                               !
!###############################################################################!
SUBROUTINE Vector_To_Coeff( Vector )

COMPLEX(idp), DIMENSION(Prob_Dim), INTENT(INOUT) :: Vector

INTEGER                                         :: ui

DO ui = iU_CF,iU_LF

    CALL Vector_To_Coeff_TypeA(Vector, ui)

END DO

DO ui = iU_S1,iU_S3

    CALL Vector_To_Coeff_TypeB(Vector, ui, iVB_S)

END DO

END SUBROUTINE Vector_To_Coeff




 !+201+####################################################!
!                                                           !
!          Coeff_To_Vector_TypeA                            !
!                                                           !
 !#########################################################!
SUBROUTINE Coeff_To_Vector_TypeA( Vector, iU )

COMPLEX(idp), DIMENSION(Prob_Dim), INTENT(INOUT) :: Vector
INTEGER                         , INTENT(IN)    :: iU

INTEGER                                         :: LM_loc, Here, There


DO LM_loc = 1,LM_Length
    here = (iU-1)*Var_Dim + (lm_loc-1)*Num_R_Nodes + 1
    there = (iU-1)*Var_Dim + lm_loc*Num_R_Nodes
    Vector(here:there) = cVA_Coeff_Vector(:,lm_loc,iU)
END DO


END SUBROUTINE Coeff_To_Vector_TypeA


 !+202+####################################################!
!                                                           !
!          Coeff_To_Vector_TypeB                            !
!                                                           !
 !#########################################################!
SUBROUTINE Coeff_To_Vector_TypeB( Vector, iU, iVB )

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


Vector(HereA:ThereA) = cVB_Coeff_Vector(HereB:ThereB,iVB)


END SUBROUTINE Coeff_To_Vector_TypeB




 !+301+####################################################!
!                                                           !
!                    Vector_To_Coeff_TypeA                  !
!                                                           !
 !#########################################################!
SUBROUTINE Vector_To_Coeff_TypeA( Vector, iU )

COMPLEX(idp), DIMENSION(Prob_Dim), INTENT(IN)        :: Vector
INTEGER                         , INTENT(IN)        :: iU

INTEGER                                             :: LM_loc, Here, There

DO LM_loc = 1,LM_Length
    here = (iU-1)*Var_Dim + (lm_loc-1)*Num_R_Nodes + 1
    there = (iU-1)*Var_Dim + lm_loc*Num_R_Nodes
    cVA_Coeff_Vector(:,lm_loc,iU) = Vector(here:There)
END DO


END SUBROUTINE Vector_To_Coeff_TypeA


 !+302+####################################################!
!                                                           !
!                    Vector_To_Coeff_TypeB                  !
!                                                           !
 !#########################################################!
SUBROUTINE Vector_To_Coeff_TypeB( Vector, iU, iVB )

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

cVB_Coeff_Vector(HereB:ThereB,iVB) = Vector(HereA:ThereA)


END SUBROUTINE Vector_To_Coeff_TypeB







!+101+##########################################################################!
!                                                                               !
!                  Coeff_To_Vector_TypeA                                  !
!                                                                               !
!###############################################################################!
SUBROUTINE Coeff_To_Vector_TypeA_SV( Vector, iU )

COMPLEX(idp), DIMENSION(Var_Dim), INTENT(INOUT) :: Vector
INTEGER                         , INTENT(IN)    :: iU

INTEGER                                         :: LM_loc, Here, There

DO LM_loc = 1,LM_Length
    Here  = (lm_loc-1)*Num_R_Nodes + 1
    There = lm_loc*Num_R_Nodes

    Vector(here:there) = cVA_Coeff_Vector(:,lm_loc,iU)
END DO


END SUBROUTINE Coeff_To_Vector_TypeA_SV




!+202+##########################################################################!
!                                                                               !
!                    Vector_To_Coeff_TypeA                                  !
!                                                                               !
!###############################################################################!
SUBROUTINE Vector_To_Coeff_TypeA_SV( Vector, iU )

COMPLEX(idp), DIMENSION(Var_Dim), INTENT(IN)        :: Vector
INTEGER                         , INTENT(IN)        :: iU

INTEGER                                             :: LM_loc, Here, There

DO LM_loc = 1,LM_Length
    Here = (lm_loc-1)*Num_R_Nodes + 1
    There = lm_loc*Num_R_Nodes
    cVA_Coeff_Vector(:,lm_loc,iU) = Vector(here:There)
END DO


END SUBROUTINE Vector_To_Coeff_TypeA_SV



END MODULE Functions_Coeffs_Module
