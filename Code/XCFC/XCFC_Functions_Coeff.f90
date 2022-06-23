   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE XCFC_Functions_Coeff_Module                                           !##!
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
                    Num_R_Nodes

USE Variables_Vectors,  &
            ONLY :  cVA_Coeff_Vector,  &
                    cVB_Coeff_Vector


IMPLICIT NONE


CONTAINS
!+101+##########################################################################!
!                                                                               !
!                  Coeff_To_Vector                                 !
!                                                                               !
!###############################################################################!
SUBROUTINE Coeff_To_Vector( Vector )

COMPLEX(idp), DIMENSION(Var_Dim), INTENT(INOUT) :: Vector

INTEGER                                         :: ui

DO ui = iU_CF,iU_LF

    CALL Coeff_To_Vector_TypeA(Vector, ui)

END DO

DO ui = iU_S1,iU_S3

    CALL Coeff_To_Vector_TypeB(Vector, ui, iVB_S)

END DO

END SUBROUTINE Coeff_To_Vector




!+101+##########################################################################!
!                                                                               !
!                  Vector_To_Coeff                                  !
!                                                                               !
!###############################################################################!
SUBROUTINE Vector_To_Coeff( Vector )

COMPLEX(idp), DIMENSION(Var_Dim), INTENT(INOUT) :: Vector

INTEGER                                         :: ui

DO ui = iU_CF,iU_LF

    CALL Vector_To_Coeff_TypeA(Vector, ui)

END DO

DO ui = iU_S1,iU_S3

    CALL Vector_To_Coeff_TypeB(Vector, ui, iVB_S)

END DO

END SUBROUTINE Vector_To_Coeff




!+101+##########################################################################!
!                                                                               !
!                  Coeff_To_Vector_TypeA                                  !
!                                                                               !
!###############################################################################!
SUBROUTINE Coeff_To_Vector_TypeA( Vector, iU )

COMPLEX(idp), DIMENSION(Var_Dim), INTENT(INOUT) :: Vector
INTEGER                         , INTENT(IN)    :: iU

INTEGER                                         :: LM_loc, Here, There

DO LM_loc = 1,LM_Length
    Here  = (lm_loc-1)*Num_R_Nodes + 1
    There = lm_loc*Num_R_Nodes

    Vector(here:there) = cVA_Coeff_Vector(:,lm_loc,iU)
END DO


END SUBROUTINE Coeff_To_Vector_TypeA


!+101+##########################################################################!
!                                                                               !
!                  Coeff_To_Vector_TypeB                                  !
!                                                                               !
!###############################################################################!
SUBROUTINE Coeff_To_Vector_TypeB( Vector, iU, iVB )

COMPLEX(idp), DIMENSION(Var_Dim), INTENT(INOUT) :: Vector
INTEGER                         , INTENT(IN)    :: iU, iVB



Vector(:) = cVB_Coeff_Vector(:,iVB)


END SUBROUTINE Coeff_To_Vector_TypeB




!+202+##########################################################################!
!                                                                               !
!                    Vector_To_Coeff_TypeA                                  !
!                                                                               !
!###############################################################################!
SUBROUTINE Vector_To_Coeff_TypeA( Vector, iU )

COMPLEX(idp), DIMENSION(Var_Dim), INTENT(IN)        :: Vector
INTEGER                         , INTENT(IN)        :: iU

INTEGER                                             :: LM_loc, Here, There

DO LM_loc = 1,LM_Length
    Here = (lm_loc-1)*Num_R_Nodes + 1
    There = lm_loc*Num_R_Nodes
    cVA_Coeff_Vector(:,lm_loc,iU) = Vector(here:There)
END DO


END SUBROUTINE Vector_To_Coeff_TypeA




!+101+##########################################################################!
!                                                                               !
!                  Coeff_To_Vector_TypeB                                        !
!                                                                               !
!###############################################################################!
SUBROUTINE Vector_To_Coeff_TypeB( Vector, iU, iVB )

COMPLEX(idp), DIMENSION(Var_Dim), INTENT(IN)    :: Vector
INTEGER                         , INTENT(IN)    :: iU, iVB


cVB_Coeff_Vector(:,iVB) = Vector(:)


END SUBROUTINE Vector_To_Coeff_TypeB



END MODULE XCFC_Functions_Coeff_Module

