   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE FP_Initial_Guess_Module                                                      !##!
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


USE Variables_Mesh, &
            ONLY :  NUM_R_ELEMENTS

USE Variables_FP, &
            ONLY :  FP_Coeff_Vector,            &
                    FP_Coeff_Vector_Beta,       &
                    CFA_EQ_Flags



USE Functions_Mapping, &
            ONLY :  Map_From_X_Space,   &
                    Map_To_X_Space

USE Functions_Quadrature, &
            ONLY :  Initialize_LGL_Quadrature_Locations

USE Functions_Math, &
            ONLY :  Lagrange_Poly

USE FP_Functions_Mapping, &
            ONLY :  FP_FEM_Node_Map

IMPLICIT NONE

CONTAINS

!+101+###########################################################################!
!                                                                                !
!                  Init_FP_Guess_Flat          !
!                                                                                !
!################################################################################!
SUBROUTINE Init_FP_Guess_Flat()


INTEGER                                     ::  beta_i, j

REAL(KIND = idp)                            ::  Beta_Start,     &
                                                delta_Beta



INTEGER                                     ::  CUR_PSI_LOC,        &
                                                CUR_ALPHPSI_LOC,    &
                                                CUR_BETA_LOC


INTEGER                                     :: re, rd, Here, i


IF ( Verbose_Flag ) THEN
    PRINT*,"Initializing Flat Space Guess"
END IF
!
!   Empty Space Initial Guess
!



Beta_Start = 0.0_idp
!Beta_Start = 1.0E-12
!Beta_Start = 0.05000_idp*2.0_idp*sqrt(pi)


delta_Beta = 0.0_idp
!delta_Beta = 0.0001_idp
!delta_Beta = Beta_Start/VAR_DIM

!
!DO i = 1,2
!    IF ( CFA_EQ_Flags(i) == 1 ) THEN
!        FP_Coeff_Vector(:,:,i) = 1.0_idp * 2.0_idp * sqrt(pi)
!    ELSE
!        FP_Coeff_Vector(:,:,i) = 0.0_idp
!    END IF
!END DO
!IF ( L_LIMIT == 0 ) THEN
!    FP_Coeff_Vector(:,:,1:2) = 1.0_idp * 2.0_idp * sqrt(pi)
!ELSE
!    FP_Coeff_Vector(:,:,1:2) = 1.0_idp
!END IF

FP_Coeff_Vector(:,:,1:2) = 1.0_idp * 2.0_idp * sqrt(pi)
FP_Coeff_Vector(:,:,3:5) = 0.0_idp






!DO re = 0,NUM_R_ELEMENTS - 1
!    DO rd = 0, Degree
!        ! 2 sqrt(pi) is Ylm normalization factor
!
!        Here = FP_FEM_Node_Map(re, rd)
!
!        FP_Coeff_Vector(Here,0,1) = 1.0_idp * 2.0_idp * sqrt(pi)
!        FP_Coeff_Vector(Here,0,2) = 1.0_idp * 2.0_idp * sqrt(pi)
!
!    END DO
!END DO


FP_Coeff_Vector_Beta = 0.0_idp



END SUBROUTINE Init_FP_Guess_Flat










!+101+###########################################################################!
!                                                                                !
!                  Init_FP_Guess_Flat          !
!                                                                                !
!################################################################################!
SUBROUTINE Init_FP_Guess_Informed(  RE, TE, PE,         &
                                    RQ, TQ, PQ,         &
                                    rlocs,              &
                                    Si                  )

INTEGER, INTENT(IN)                                                     ::  RE, TE, PE
INTEGER, INTENT(IN)                                                     ::  RQ, TQ, PQ

REAL(idp), INTENT(IN), DIMENSION(0:RE)                                  ::  rlocs
REAL(idp), INTENT(IN), DIMENSION(1:RQ*TQ*PQ, 1:RE, 1:TE, 1:PE, 1:3)     ::  Si

PRINT*,"The subroutine 'Init_FP_Guess_Informed' is incomplete."
PRINT*,"Please try another subroutine. "
PRINT*,"Poseidon is STOPing the program. "
STOP

IF ( Verbose_Flag ) THEN
    PRINT*,"Initializing Informed Guess"
END IF





END SUBROUTINE Init_FP_Guess_Informed





!+101+###########################################################################!
!                                                                                !
!                  Input_FP_Guess                                                !
!                                                                                !
!################################################################################!
SUBROUTINE Input_FP_Guess(  Psi_Guess,                                  &
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

INTEGER                                                                             ::  re, rq, rqb

REAL(KIND = idp), DIMENSION(0:DEGREE)                                               ::  Node_Locs
REAL(idp), DIMENSION(1:Input_RQ)                                                    ::  Lagrange_Poly_Value
REAL(idp), DIMENSION(1:5)                                                           ::  Tmp_Value
INTEGER                                                                             ::  Here
REAL(idp), DIMENSION(1:Input_RQ)                                                    ::  Mapped_R_Quad

Node_Locs = Initialize_LGL_Quadrature_Locations(DEGREE)

Num_Input_DOF = Input_RQ*Input_TQ*Input_PQ



Mapped_R_Quad = Map_To_X_Space( Left_Limit, Right_Limit, Input_R_Quad )



DO re = 1,Input_RE
    DO rqb = 0,Degree
        Lagrange_Poly_Value = Lagrange_Poly(Node_Locs(rqb), Input_RQ-1, Mapped_R_Quad)

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

        Here = (re-1)*Degree + rqb + 1
        FP_Coeff_Vector(Here,1,:) = 2.0_idp*Sqrt(pi)*Tmp_Value(:)
        
    END DO ! rqb
END DO ! re

FP_Coeff_Vector(:,:,3:5) = 0.0_idp
FP_Coeff_Vector_Beta = 0.0_idp



END SUBROUTINE Input_FP_Guess













END MODULE FP_Initial_Guess_Module
