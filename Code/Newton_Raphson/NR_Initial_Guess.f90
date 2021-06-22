   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE NR_Initial_Guess_Module                                                      !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   Initialize_Flat_Space_Guess_Values                                  !##!
!##!    +102+   Initialize_Special_Guess_Values                                     !##!
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

USE Units_Module, &
            ONLY :  C_Square

USE Poseidon_Parameters, &
            ONLY :  DOMAIN_DIM,             &
                    DEGREE,                 &
                    L_LIMIT,                &
                    Verbose_Flag
USE Variables_NR, &
            ONLY :  NR_Coeff_Vector

USE Variables_Mesh, &
            ONLY :  NUM_R_ELEMENTS

USE Functions_Mapping, &
            ONLY :  Map_To_X_Space

USE NR_Mapping_Functions, &
            ONLY :  NR_Array_Map

USE Functions_Quadrature, &
            ONLY :  Initialize_LGL_Quadrature_Locations

USE Functions_Math, &
            ONLY :  Lagrange_Poly


IMPLICIT NONE

CONTAINS



!+101+###########################################################################!
!                                                                                !
!                  Input_FP_Guess                                                !
!                                                                                !
!################################################################################!
SUBROUTINE NR_Input_Guess(  Psi_Guess,                                  &
                            AlphaPsi_Guess,                             &
                            Beta_Guess,                                 &
                            Input_RE, Input_TE, Input_PE,               &
                            Input_RQ, Input_TQ, Input_PQ,               &
                            Input_R_Quad, Input_T_Quad, Input_P_Quad,   &
                            Left_Limit, Right_Limit                     )


REAL(idp), DIMENSION(1:Input_RQ*Input_TQ*Input_PQ,1:Input_RE,1:Input_TE,1:Input_PE),        &
                                                                        INTENT(IN)  :: Psi_Guess

REAL(idp), DIMENSION(1:Input_RQ*Input_TQ*Input_PQ,1:Input_RE,1:Input_TE,1:Input_PE),        &
INTENT(IN)  :: AlphaPsi_Guess

REAL(idp), DIMENSION(1:Input_RQ*Input_TQ*Input_PQ,1:Input_RE,1:Input_TE,1:Input_PE,1:3),        &
INTENT(IN)  :: Beta_Guess

INTEGER, INTENT(IN)                                                                 :: Input_RE
INTEGER, INTENT(IN)                                                                 :: Input_TE
INTEGER, INTENT(IN)                                                                 :: Input_PE

INTEGER, INTENT(IN)                                                                 :: Input_RQ
INTEGER, INTENT(IN)                                                                 :: Input_TQ
INTEGER, INTENT(IN)                                                                 :: Input_PQ

REAL(idp), DIMENSION(1:Input_RQ), INTENT(IN)                                        :: Input_R_Quad
REAL(idp), DIMENSION(1:Input_TQ), INTENT(IN)                                        :: Input_T_Quad
REAL(idp), DIMENSION(1:Input_PQ), INTENT(IN)                                        :: Input_P_Quad

REAL(idp), INTENT(IN)                                                               :: Left_Limit
REAL(idp), INTENT(IN)                                                               :: Right_Limit


INTEGER                                                                             :: Num_Input_DOF

INTEGER                                                                             ::  re, rq, d

REAL(KIND = idp), DIMENSION(0:DEGREE)                                               ::  Node_Locs
REAL(idp), DIMENSION(1:Input_RQ)                                                    ::  Lagrange_Poly_Value
REAL(idp), DIMENSION(1:5)                                                           ::  Tmp_Value
INTEGER                                                                             ::  Here
REAL(idp), DIMENSION(1:Input_RQ)                                                    ::  Mapped_R_Quad

Node_Locs = Initialize_LGL_Quadrature_Locations(DEGREE)

Num_Input_DOF = Input_RQ*Input_TQ*Input_PQ



Mapped_R_Quad = Map_To_X_Space( Left_Limit, Right_Limit, Input_R_Quad )


NR_Coeff_Vector = 0.0_idp
DO re = 1,Input_RE
    DO d = 0,Degree
        Lagrange_Poly_Value = Lagrange_Poly(Node_Locs(d), Input_RQ-1, Mapped_R_Quad)

        Tmp_Value = 0.0_idp
        DO rq = 1,Input_RQ

            Here =  (rq - 1)*Input_TQ*Input_PQ

            Tmp_Value(1) = Tmp_Value(1)          &
                      + Psi_Guess(rq,re,1,1)*Lagrange_Poly_Value(rq)

            Tmp_Value(2) = Tmp_Value(2)           &
                      + AlphaPsi_Guess(rq,re,1,1)*Lagrange_Poly_Value(rq)

            Tmp_Value(3) = Tmp_Value(3)         &
                        + Beta_Guess(rq,re,1,1,1)*Lagrange_Poly_Value(rq)

        END DO ! rq

        Here = NR_Array_Map(re-1,d,1,0)
        NR_Coeff_Vector(Here) = 2.0_idp*Sqrt(pi)*Tmp_Value(1)

        Here = NR_Array_Map(re-1,d,2,0)
        NR_Coeff_Vector(Here) = 2.0_idp*Sqrt(pi)*Tmp_Value(2)

        Here = NR_Array_Map(re-1,d,3,0)
        NR_Coeff_Vector(Here) = 2.0_idp*Sqrt(pi)*Tmp_Value(3)

        
    END DO ! d
END DO ! re




END SUBROUTINE NR_Input_Guess













END MODULE NR_Initial_Guess_Module

