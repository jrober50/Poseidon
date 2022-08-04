    !#################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Return_Functions_FP                                                         !##!
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
            ONLY :  Degree,                 &
                    L_Limit

USE Parameters_Variable_Indices, &
            ONLY :  iU_CF,                  &
                    iU_LF,                  &
                    iU_S1,                  &
                    iU_S2,                  &
                    iU_S3,                  &
                    iVB_S

USE Variables_Mesh, &
            ONLY :  Num_R_Elements,         &
                    rlocs

USE Variables_Derived, &
            ONLY :  LM_Length,              &
                    ULM_Length

USE Variables_Vectors, &
            ONLY :  cVA_Coeff_Vector,      &
                    cVB_Coeff_Vector

USE Functions_Quadrature, &
            ONLY :  Initialize_LGL_Quadrature_Locations,    &
                    Initialize_LGL_Quadrature

USE Functions_Math, &
            ONLY :  Lagrange_Poly,          &
                    Spherical_Harmonic

USE Maps_Fixed_Point, &
            ONLY :  FP_Array_Map_TypeB, &
                    FP_Array_Map

USE Maps_Domain, &
            ONLY :  Map_To_lm,          &
                    Map_To_FEM_Node


USE Maps_X_Space,    &
            ONLY :  Map_To_X_Space


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

INTEGER                                                         ::  re

REAL(KIND = idp), DIMENSION(0:DEGREE)                           ::  xlocP, weightP


Tmp_U_Value = 0.0_idp


IF ( r .LE. rlocs(0) ) THEN


    LagP    = 0.0_idp
    LagP(0) = 1.0_idp

    CALL Calc_Values_Here_All( 0, theta, phi, LagP, Tmp_U_Value )


ELSE

    CALL Initialize_LGL_Quadrature(DEGREE,xlocP,weightP)

    DO re = 0,NUM_R_ELEMENTS-1
    IF ( r > rlocs(re) .AND. r <= rlocs(re+1) ) THEN


        r_tmp = Map_To_X_Space(rlocs(re),rlocs(re+1),r)
        LagP = Lagrange_Poly(r_tmp,DEGREE,xlocP)

        
        CALL Calc_Values_Here_All( re, theta, phi, LagP, Tmp_U_Value )

        EXIT
    END IF
    END DO

    IF ( r > rlocs(NUM_R_ELEMENTS) ) THEN

        LagP         = 0.0_idp
        LagP(DEGREE) = 1.0_idp

        CALL Calc_Values_Here_All( Num_R_Elements-1, theta, phi, LagP, Tmp_U_Value )

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

        DO u = 1,2
            DO d = 0,DEGREE

                Current_Location = Map_To_FEM_Node(re,d)
                Tmp_U_Value(u) = Tmp_U_Value(u) + cVA_Coeff_Vector(Current_Location,1,u)  &
                                                * LagP(d) * Spherical_Harmonic(0,0,pi,pi/2.0_idp)
            END DO ! d Loop
        END DO ! u Loop

        DO d = 0,DEGREE

            Current_Location = FP_Array_Map_TypeB(iU_S1,iVB_S,re,d,0)
            Tmp_U_Value(u) = Tmp_U_Value(u) + cVB_Coeff_Vector(Current_Location,iVB_S)  &
                                            * LagP(d) * Spherical_Harmonic(0,0,pi,pi/2.0_idp)
        END DO ! d Loop



        CFA_ConFactor(x,re+1,1,1) = REAL(Tmp_U_Value(1), KIND = idp)
        CFA_Lapse(x,re+1,1,1)     = REAL(Tmp_U_Value(2), KIND = idp)        &
                                  / REAL(Tmp_U_Value(1), KIND = idp)
        CFA_Shift(x,re+1,1,1)     = REAL(Tmp_U_Value(3), KIND = idp)
    END DO ! x Loop
END DO ! re Loop


END SUBROUTINE Calc_1D_CFA_Values_FP











 !+501+########################################################!
!                                                               !
!          Calc_Var_At_Location_Type_A                          !
!                                                               !
 !#############################################################!
FUNCTION Calc_Var_At_Location_Type_A( r, theta, phi, iU)

REAL(idp)                                           ::  Calc_Var_At_Location_Type_A
REAL(idp),  INTENT(IN)                              ::  r, theta, phi
INTEGER,    INTENT(IN)                              ::  iU



COMPLEX(KIND = idp)                                 ::  Tmp_U_Value


REAL(KIND = idp)                                    ::  r_tmp
REAL(KIND = idp), DIMENSION(0:DEGREE)               ::  LagP

INTEGER                                             ::  re

REAL(KIND = idp), DIMENSION(0:DEGREE)               ::  xlocP, weightP


Tmp_U_Value = 0.0_idp


IF ( r .LE. rlocs(0) ) THEN

    LagP    = 0.0_idp
    LagP(0) = 1.0_idp

    Tmp_U_Value = Calc_Values_Here_Type_A( 0, theta, phi, LagP, iU )


ELSE

    CALL Initialize_LGL_Quadrature(DEGREE,xlocP,weightP)

    DO re = 0,NUM_R_ELEMENTS-1
    
    IF ( r > rlocs(re) .AND. r <= rlocs(re+1) ) THEN


        r_tmp = Map_To_X_Space(rlocs(re),rlocs(re+1),r)
        LagP = Lagrange_Poly(r_tmp,DEGREE,xlocP)

        
        Tmp_U_Value = Calc_Values_Here_Type_A( re, theta, phi, LagP, iU )

        EXIT
    END IF
    END DO

    IF ( r > rlocs(NUM_R_ELEMENTS) ) THEN

        LagP         = 0.0_idp
        LagP(DEGREE) = 1.0_idp

        Tmp_U_Value = Calc_Values_Here_Type_A( Num_R_Elements-1, theta, phi, LagP, iU )

    END IF

END IF



Calc_Var_At_Location_Type_A = REAL(Tmp_U_Value, KIND = idp)



END FUNCTION Calc_Var_At_Location_Type_A




 !+502+########################################################!
!                                                               !
!          Calc_Var_At_Location_Type_B                          !
!                                                               !
 !#############################################################!
FUNCTION Calc_Var_At_Location_Type_B( r, theta, phi, iU, iVB  )


REAL(idp)                                           ::  Calc_Var_At_Location_Type_B
REAL(idp),  INTENT(IN)                              ::  r, theta, phi
INTEGER,    INTENT(IN)                              ::  iU, iVB



COMPLEX(KIND = idp)                                 ::  Tmp_U_Value


REAL(KIND = idp)                                    ::  r_tmp
REAL(KIND = idp), DIMENSION(0:DEGREE)               ::  LagP

INTEGER                                             ::  re

REAL(KIND = idp), DIMENSION(0:DEGREE)               ::  xlocP, weightP


Tmp_U_Value = 0.0_idp


IF ( r .LE. rlocs(0) ) THEN


    LagP    = 0.0_idp
    LagP(0) = 1.0_idp

    Tmp_U_Value = Calc_Values_Here_Type_B( 0, theta, phi, LagP, iU, iVB )


ELSE

    CALL Initialize_LGL_Quadrature(DEGREE,xlocP,weightP)

    DO re = 0,NUM_R_ELEMENTS-1
    IF ( r > rlocs(re) .AND. r <= rlocs(re+1) ) THEN


        r_tmp = Map_To_X_Space(rlocs(re),rlocs(re+1),r)
        LagP = Lagrange_Poly(r_tmp,DEGREE,xlocP)

        Tmp_U_Value = Calc_Values_Here_Type_B( re, theta, phi, LagP, iU, iVB )


        EXIT
    END IF
    END DO

    IF ( r > rlocs(NUM_R_ELEMENTS) ) THEN

        LagP         = 0.0_idp
        LagP(DEGREE) = 1.0_idp

        Tmp_U_Value = Calc_Values_Here_Type_B( Num_R_Elements-1, theta, phi, LagP, iU, iVB )

    END IF

END IF



Calc_Var_At_Location_Type_B = REAL(Tmp_U_Value, KIND = idp)


END FUNCTION Calc_Var_At_Location_Type_B







!+102+###########################################################################!
!                                                                                !
!                  Calc_1D_CFA_Values_FP          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Values_Here_All( re, theta, phi, LagP, Tmp_U_Value )

INTEGER,                        INTENT(IN)                      ::  re
REAL(idp),                      INTENT(IN)                      ::  theta, phi
REAL(idp), DIMENSION(0:DEGREE), INTENT(INOUT)                   ::  LagP
COMPLEX(idp), DIMENSION(1:5),   INTENT(INOUT)                   ::  Tmp_U_Value


INTEGER                                                         ::  l, m, d, u
INTEGER                                                         ::  Loc_RED, Loc_LM


Tmp_U_Value = 0.0_idp
DO u = iU_CF,iU_LF
DO l = 0,L_Limit
DO m = -l,l
DO d = 0,DEGREE

    Loc_RED = Map_To_FEM_Node(re,d)
    Loc_LM  = Map_to_lm(l,m)


    Tmp_U_Value(u) = Tmp_U_Value(u)                         &
                    + cVA_Coeff_Vector(Loc_RED,Loc_LM,u)   &
                    * Spherical_Harmonic(l,m,theta,phi)     &
                    * LagP(d)


!    IF ( u == iU_CF) THEN
!        PRINT*,"A",Loc_RED,Loc_LM,re,d
!        PRINT*,cVA_Coeff_Vector(Loc_RED,Loc_LM,u),  &
!                Spherical_Harmonic(l,m,theta,phi),     &
!                LagP(d)
!    END IF


END DO  !   d Loop
END DO  !   m Loop
END DO  !   l Loop
END DO  !   u Loop


DO u = iU_S1,iU_S3
DO l = 0,L_Limit
DO m = -l,l
DO d = 0,DEGREE

    Loc_RED = FP_Array_Map_TypeB(u,iVB_S,re,d,l,m)
    Tmp_U_Value(u) = Tmp_U_Value(u)                         &
                    + cVB_Coeff_Vector(Loc_RED,iVB_S)     &
                    * Spherical_Harmonic(l,m,theta,phi)     &
                    * LagP(d)



END DO  !   d Loop
END DO  !   m Loop
END DO  !   l Loop
END DO  !   u Loop


END SUBROUTINE Calc_Values_Here_All





!+102+###########################################################################!
!                                                                                !
!                  Calc_1D_CFA_Values_FP          !
!                                                                                !
!################################################################################!
FUNCTION Calc_Values_Here_Type_A( re, theta, phi, LagP, iU )

REAL(idp)                                               ::  Calc_Values_Here_Type_A
INTEGER,                        INTENT(IN)              ::  re
REAL(idp),                      INTENT(IN)              ::  theta, phi
REAL(idp), DIMENSION(0:DEGREE), INTENT(IN)              ::  LagP
INTEGER,                        INTENT(IN)              ::  iU


COMPLEX(idp)                                            ::  Tmp_U_Value


INTEGER                                                 ::  l, m, d
INTEGER                                                 ::  Loc_RED, Loc_LM


Tmp_U_Value = 0.0_idp
DO l = 0,L_Limit
DO m = -l,l
DO d = 0,DEGREE

    Loc_RED = Map_To_FEM_Node(re,d)
    Loc_LM  = Map_to_lm(l,m)

    Tmp_U_Value = Tmp_U_Value                           &
                + cVA_Coeff_Vector(Loc_RED,Loc_LM,iU)  &
                * Spherical_Harmonic(l,m,theta,phi)     &
                * LagP(d)
!    IF ( iU == iU_CF) THEN
!        PRINT*,"B",Loc_RED,Loc_LM,re,d
!        PRINT*,cVA_Coeff_Vector(Loc_RED,Loc_LM,iU),  &
!                Spherical_Harmonic(l,m,theta,phi),     &
!                LagP(d)
!    END IF
END DO  !   d Loop
END DO  !   m Loop
END DO  !   l Loop


Calc_Values_Here_Type_A = REAL( Tmp_U_Value, KIND = idp )

END FUNCTION Calc_Values_Here_Type_A




!+603+###########################################################################!
!                                                                                !
!                  Calc_1D_CFA_Values_FP          !
!                                                                                !
!################################################################################!
FUNCTION Calc_Values_Here_Type_B( re, theta, phi, LagP, iU, iVB )


REAL(idp)                                                       ::  Calc_Values_Here_Type_B
INTEGER,                        INTENT(IN)                      ::  re
REAL(idp),                      INTENT(IN)                      ::  theta, phi
REAL(idp), DIMENSION(0:DEGREE), INTENT(IN)                      ::  LagP
INTEGER,                        INTENT(IN)                      ::  iU
INTEGER,                        INTENT(IN)                      ::  iVB


COMPLEX(idp)                                                    ::  Tmp_U_Value
INTEGER                                                         ::  l, m, d
INTEGER                                                         ::  Loc_RED


Tmp_U_Value = 0.0_idp
DO l = 0,L_Limit
DO m = -l,l
DO d = 0,DEGREE

    Loc_RED = FP_Array_Map_TypeB(iU,iVB,re,d,l,m)
    Tmp_U_Value = Tmp_U_Value                           &
                + cVB_Coeff_Vector(Loc_RED,iVB)        &
                * Spherical_Harmonic(l,m,theta,phi)     &
                * LagP(d)


END DO  !   d Loop
END DO  !   m Loop
END DO  !   l Loop


Calc_Values_Here_Type_B = REAL( Tmp_U_Value, KIND = idp)

END FUNCTION Calc_Values_Here_Type_B







END MODULE Return_Functions_FP

