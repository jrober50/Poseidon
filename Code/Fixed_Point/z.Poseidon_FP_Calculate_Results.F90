   !#################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_FP_Calculate_Results_Module                                         !##!
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
USE Poseidon_Constants_Module, &
            ONLY :  idp, pi


USE Poseidon_Parameters, &
            ONLY :  DEGREE,                 &
                    L_LIMIT,                &
                    NUM_CFA_VARS

USE Poseidon_Variables_Module,    &
            ONLY :  NUM_R_ELEMENTS,         &
                    ULM_LENGTH,             &
                    LM_LENGTH,              &
                    Coefficient_Vector,     &
                    M_VALUES,               &
                    rlocs,                  &
                    Matrix_Location

USE Poseidon_FP_Variables_Module,           &
            ONLY :  FP_Coeff_Vector

USE Poseidon_Mapping_Functions_Module,    &
            ONLY :  Map_To_X_Space

USE Poseidon_Quadrature_Module, &
            ONLY :  Initialize_LGL_Quadrature_Locations,    &
                    Initialize_LGL_Quadrature

USE Poseidon_Math_Functions_Module, &
            ONLY :  Lagrange_Poly,          &
                    Spherical_Harmonic

USE Poseidon_FP_Mapping_Functions_Module, &
            ONLY :  FP_Vector_Map,  &
                    FP_LM_Map

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
        

                Loc_RED = FP_Vector_Map(0,0)
                LOC_LM  = FP_LM_Map(l,m)
                Tmp_U_Value(u) = Tmp_U_Value(u)                         &
                                + FP_Coeff_Vector(Loc_RED,LOC_LM,u)     &
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


                            Loc_RED = FP_Vector_Map(re,d)
                            LOC_LM  = FP_LM_Map(l,m)
                            Tmp_U_Value(u) = Tmp_U_Value(u)                         &
                                            + FP_Coeff_Vector(Loc_RED,LOC_LM,u)     &
                                            * Spherical_Harmonic(l,m,theta,phi)     &
                                            * LagP(d)

                        END DO  !   d Loop
                    END DO  !   m Loop
                END DO  !   l Loop
            END DO ! u Loop

            EXIT
        END IF

    END DO

    IF ( r > rlocs(NUM_R_ELEMENTS) ) THEN

        DO u = 1,NUM_CFA_VARS
            DO l = 0,L_Limit
                DO m = -M_VALUES(l),M_VALUES(l)

                    Loc_RED = FP_Vector_Map(NUM_R_ELEMENTS-1,DEGREE)
                    LOC_LM  = FP_LM_Map(l,m)
                    Tmp_U_Value(u) = Tmp_U_Value(u)                         &
                                    + FP_Coeff_Vector(Loc_RED,LOC_LM,u)     &
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






END MODULE Poseidon_FP_Calculate_Results_Module

