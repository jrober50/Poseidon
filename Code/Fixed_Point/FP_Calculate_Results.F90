   !#################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE FP_Functions_Results                                                         !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   Calc_3D_Values_At_Location                                          !##!
!##!    +102+   Calc_1D_CFA_Values                                                  !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !#################################################################################!



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
            ONLY :  DEGREE,                 &
                    L_LIMIT,                &
                    NUM_CFA_VARS


USE Variables_Mesh, &
            ONLY :  Num_R_Elements,         &
                    rlocs

USE Variables_Derived, &
            ONLY :  LM_Length,              &
                    ULM_Length

USE Variables_Functions, &
            ONLY : Matrix_Location

USE Variables_FP, &
            ONLY :  FP_Coeff_Vector

USE Variables_Tables, &
            ONLY :  M_Values


USE Functions_Mapping,    &
            ONLY :  Map_To_X_Space

USE Functions_Quadrature, &
            ONLY :  Initialize_LGL_Quadrature_Locations,    &
                    Initialize_LGL_Quadrature

USE Functions_Math, &
            ONLY :  Lagrange_Poly,          &
                    Spherical_Harmonic

USE FP_Functions_Mapping, &
            ONLY :  FP_FEM_Node_Map,    &
                    FP_LM_Map,          &
                    FP_Array_Map

IMPLICIT NONE





CONTAINS





!+101+###########################################################################!
!                                                                                !
!                  Calc_3D_Values_At_Location          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_FP_Values_At_Location( r, theta, phi, Return_Psi, Return_AlphaPsi,  &
                                        Return_Beta1, Return_Beta2, Return_Beta3    )


REAL(KIND = idp), INTENT(IN)                                ::  r, theta, phi
REAL(KIND = idp), INTENT(INOUT)                             ::  Return_Psi,         &
                                                                Return_AlphaPsi,    &
                                                                Return_Beta1,       &
                                                                Return_Beta2,       &
                                                                Return_Beta3



COMPLEX(KIND = idp), DIMENSION(1:5)                         ::  Tmp_U_Value


REAL(KIND = idp)                                                ::  r_tmp
REAL(KIND = idp), DIMENSION(0:DEGREE)                           ::  LagP

INTEGER                                                         ::  re, l, m, d, u


INTEGER                                                         :: Current_Location
INTEGER                                                         :: Loc_RED, Loc_LM

COMPLEX(KIND = idp)                                             ::  TMP_VALUE_A
REAL(KIND = idp), DIMENSION(0:DEGREE)                           ::  xlocP, weightP


Tmp_U_Value = 0.0_idp


IF ( r == rlocs(0) ) THEN

    DO u = 1,NUM_CFA_VARS
    DO l = 0,L_Limit
    DO m = -M_VALUES(l),M_VALUES(l)

        Loc_RED = FP_FEM_Node_Map(0,0)
        Loc_LM  = FP_LM_Map(l,m)
        Tmp_U_Value(u) = Tmp_U_Value(u)                         &
                        + FP_Coeff_Vector(Loc_RED,Loc_LM,u)     &
                        * Spherical_Harmonic(l,m,theta,phi)


    END DO ! m Loop
    END DO  ! l Loop
    END DO  ! u Loop



ELSE

    CALL Initialize_LGL_Quadrature(DEGREE,xlocP,weightP)

    DO re = 0,NUM_R_ELEMENTS-1

        IF ( r > rlocs(re) .AND. r <= rlocs(re+1) ) THEN

            r_tmp = Map_To_X_Space(rlocs(re),rlocs(re+1),r)
            LagP = Lagrange_Poly(r_tmp,DEGREE,xlocP)

            DO u = 1,NUM_CFA_VARS
            DO l = 0,L_Limit
            DO m = -M_VALUES(l),M_VALUES(l)
            DO d = 0,DEGREE

                Loc_RED = FP_FEM_Node_Map(re,d)
                Loc_LM  = FP_LM_Map(l,m)

!                IF ( u == 1 ) THEN
!                IF ( (re >= 118) .and. (re <= 122) ) THEN
!!                    PRINT*,Loc_RED,Loc_LM,u,FP_Coeff_Vector(Loc_RED,Loc_LM,u)
!                    PRINT*,Tmp_U_Value(u),                        &
!                            FP_Coeff_Vector(Loc_RED,Loc_LM,u),     &
!                            Spherical_Harmonic(l,m,theta,phi),     &
!                            LagP(d)
!                END IF
!                END IF
                Tmp_U_Value(u) = Tmp_U_Value(u)                         &
                                + FP_Coeff_Vector(Loc_RED,Loc_LM,u)     &
                                * Spherical_Harmonic(l,m,theta,phi)     &
                                * LagP(d)

            END DO  !   d Loop
            END DO  !   m Loop
            END DO  !   l Loop
            END DO  !   u Loop

!            IF ( (re >= 118) .and. (re <= 122) ) THEN
!!                    PRINT*,Loc_RED,Loc_LM,u,FP_Coeff_Vector(Loc_RED,Loc_LM,u)
!                PRINT*,Tmp_U_Value(1), RE
!            END IF

            EXIT
        END IF

    END DO

    IF ( r > rlocs(NUM_R_ELEMENTS) ) THEN

        DO u = 1,NUM_CFA_VARS
        DO l = 0,L_Limit
        DO m = -M_VALUES(l),M_VALUES(l)


            Loc_RED = FP_FEM_Node_Map(Num_R_Elements-1,Degree)
            Loc_LM  = FP_LM_Map(l,m)
            Tmp_U_Value(u) = Tmp_U_Value(u)                         &
                            + FP_Coeff_Vector(Loc_RED,Loc_LM,u)     &
                            * Spherical_Harmonic(l,m,theta,phi)


        END DO  !   m Loop
        END DO  !   l Loop
        END DO ! u Loop
    END IF

END IF



Return_Psi      = REAL(Tmp_U_Value(1), KIND = idp)
Return_AlphaPsi = REAL(Tmp_U_Value(2), KIND = idp)!/REAL(Tmp_U_Value(1), KIND = idp)
Return_Beta1    = REAL(Tmp_U_Value(3), KIND = idp)
Return_Beta2    = REAL(Tmp_U_Value(4), KIND = idp)
Return_Beta3    = REAL(Tmp_U_Value(5), KIND = idp)


END SUBROUTINE Calc_FP_Values_At_Location









!+102+###########################################################################!
!                                                                                !
!                  Calc_1D_CFA_Values_FP          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_1D_CFA_Values_FP(   Num_RE_Input, Num_RQ_Input, RQ_Input,   &
                                    Left_Limit, Right_Limit,                &
                                    CFA_Lapse, CFA_ConFactor, CFA_Shift     )



INTEGER, INTENT(IN)                                         ::  Num_RE_Input,   &
                                                                Num_RQ_Input

REAL(KIND = idp), DIMENSION(1:Num_RQ_Input), INTENT(IN)     ::  RQ_Input
REAL(KIND = idp), INTENT(IN)                                ::  Left_Limit,     &
                                                                Right_Limit


REAL(KIND = idp), DIMENSION(1:NUM_RQ_Input,1:NUM_RE_Input, 1, 1), INTENT(OUT) ::  CFA_Lapse
REAL(KIND = idp), DIMENSION(1:NUM_RQ_Input,1:NUM_RE_Input, 1, 1), INTENT(OUT) ::  CFA_ConFactor
REAL(KIND = idp), DIMENSION(1:NUM_RQ_Input,1:NUM_RE_Input, 1, 1), INTENT(OUT) ::  CFA_Shift



INTEGER                                                         ::  re, x, u, d

REAL(KIND = idp)                                                ::  Quad_Span
REAL(KIND = idp), DIMENSION(0:DEGREE)                           ::  Local_Locations
REAL(KIND = idp), DIMENSION(0:DEGREE)                           ::  LagP
REAL(KIND = idp), DIMENSION(1:Num_RQ_Input)                     ::  CUR_X_LOCS
COMPLEX(KIND = idp), DIMENSION(1:3)                             ::  TMP_U_Value
INTEGER                                                         ::  Current_Location

Quad_Span = Right_Limit - Left_Limit
Local_Locations = Initialize_LGL_Quadrature_Locations(DEGREE)

DO re = 0,NUM_R_ELEMENTS-1

    CUR_X_LOCS = 2.0_idp * ( RQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp

    DO x = 1,Num_RQ_Input
    
        LagP = Lagrange_Poly(CUR_X_LOCS(x),DEGREE,Local_Locations)
        Tmp_U_Value = 0.0_idp

        DO u = 1,3
            DO d = 0,DEGREE

                Current_Location = FP_FEM_Node_Map(re,d)
                Tmp_U_Value(u) = Tmp_U_Value(u) + FP_Coeff_Vector(Current_Location,1,u)  &
                                                * LagP(d) * Spherical_Harmonic(0,0,pi,pi/2.0_idp)
            END DO ! d Loop
        END DO ! u Loop

        CFA_ConFactor(x,re+1,1,1) = REAL(Tmp_U_Value(1), KIND = idp)
        CFA_Lapse(x,re+1,1,1)     = REAL(Tmp_U_Value(2), KIND = idp)        &
                                  / REAL(Tmp_U_Value(1), KIND = idp)
        CFA_Shift(x,re+1,1,1)     = REAL(Tmp_U_Value(3), KIND = idp)
    END DO ! x Loop
END DO ! re Loop


END SUBROUTINE Calc_1D_CFA_Values_FP











END MODULE FP_Functions_Results

