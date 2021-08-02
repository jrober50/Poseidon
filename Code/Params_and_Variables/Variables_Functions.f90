   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Variables_Functions                                                          !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!


USE Poseidon_Kinds_Module, &
            ONLY : idp


IMPLICIT NONE




ABSTRACT INTERFACE
    PURE FUNCTION Matrix_Location_Pointer(ui, l, m, re, d)
        INTEGER                                             ::  Matrix_Location_Pointer
        INTEGER, INTENT(IN)                                 ::  ui, l, m, re, d
    END FUNCTION Matrix_Location_Pointer
END INTERFACE

PROCEDURE(Matrix_Location_Pointer), POINTER                 ::   Matrix_Location => NULL()






ABSTRACT INTERFACE
    PURE FUNCTION LM_Location_Pointer(l, m)
        INTEGER                                             ::  LM_Location_Pointer
        INTEGER, INTENT(IN)                                 ::  l, m
    END FUNCTION LM_Location_Pointer
END INTERFACE

PROCEDURE(LM_Location_Pointer), POINTER                     ::   LM_Location => NULL()









ABSTRACT INTERFACE
    SUBROUTINE Calc_At_Location_Pointer(r, theta, phi,                          &
                                        Return_Psi, Return_AlphaPsi,            &
                                        Return_Beta1, Return_Beta2, Return_Beta3 )
        REAL(KIND = KIND(1.D0)), INTENT(IN)            ::  r
        REAL(KIND = KIND(1.D0)), INTENT(IN)            ::  theta
        REAL(KIND = KIND(1.D0)), INTENT(IN)            ::  phi
        REAL(KIND = KIND(1.D0)), INTENT(INOUT)         ::  Return_Psi
        REAL(KIND = KIND(1.D0)), INTENT(INOUT)         ::  Return_AlphaPsi
        REAL(KIND = KIND(1.D0)), INTENT(INOUT)         ::  Return_Beta1
        REAL(KIND = KIND(1.D0)), INTENT(INOUT)         ::  Return_Beta2
        REAL(KIND = KIND(1.D0)), INTENT(INOUT)         ::  Return_Beta3

    END SUBROUTINE Calc_At_Location_Pointer
END INTERFACE

PROCEDURE(Calc_At_Location_Pointer), POINTER    ::   Calc_3D_Values_At_Location => NULL()




ABSTRACT INTERFACE
    SUBROUTINE Calc_1D_Array_Pointer(   Num_RE_Input, Num_RQ_Input, RQ_Input,   &
                                        Left_Limit, Right_Limit,                &
                                        CFA_Lapse, CFA_ConFactor, CFA_Shift     )
        INTEGER, INTENT(IN)                                         ::  Num_RE_Input,   &
                                                                        Num_RQ_Input

        REAL(KIND = KIND(1.D0)), DIMENSION(1:Num_RQ_Input), INTENT(IN)     ::  RQ_Input
        REAL(KIND = KIND(1.D0)), INTENT(IN)                                ::  Left_Limit,     &
                                                                                Right_Limit


        REAL(KIND = KIND(1.D0)), DIMENSION(1:NUM_RQ_Input,1:NUM_RE_Input, 1, 1), INTENT(OUT) ::  CFA_Lapse
        REAL(KIND = KIND(1.D0)), DIMENSION(1:NUM_RQ_Input,1:NUM_RE_Input, 1, 1), INTENT(OUT) ::  CFA_ConFactor
        REAL(KIND = KIND(1.D0)), DIMENSION(1:NUM_RQ_Input,1:NUM_RE_Input, 1, 1), INTENT(OUT) ::  CFA_Shift
    END SUBROUTINE Calc_1D_Array_Pointer
END INTERFACE

PROCEDURE(Calc_1D_Array_Pointer), POINTER    ::   Calc_1D_CFA_Values => NULL()





ABSTRACT INTERFACE
    FUNCTION Source_Function_Pointer(r, theta, phi)
        REAL(KIND = KIND(1.D0))                         ::  Source_Function_Pointer
        REAL(KIND = KIND(1.D0)), INTENT(IN)             ::  r, theta, phi
    END FUNCTION Source_Function_Pointer
END INTERFACE

PROCEDURE(Source_Function_Pointer), POINTER             ::  Source_Function => NULL()








ABSTRACT INTERFACE
    FUNCTION Solution_Function_Pointer(r, theta, phi)
        REAL(KIND = KIND(1.D0))                         ::  Solution_Function_Pointer
        REAL(KIND = KIND(1.D0)), INTENT(IN)             ::  r, theta, phi
    END FUNCTION Solution_Function_Pointer
END INTERFACE
PROCEDURE(Solution_Function_Pointer), POINTER           ::   Potential_Solution => NULL()








ABSTRACT INTERFACE
    FUNCTION Shift_Function_Pointer(r, r_locs, num_r_elems)
        REAL(KIND = KIND(1.D0))                         ::  Shift_Function_Pointer
        REAL(KIND = KIND(1.D0)), INTENT(IN)             ::  r
        INTEGER, INTENT(IN)                             ::  num_r_elems
        REAL(KIND = KIND(1.D0)), DIMENSION(0:num_r_elems), INTENT(IN) :: r_locs
    END FUNCTION Shift_Function_Pointer
END INTERFACE
PROCEDURE(Shift_Function_Pointer), POINTER              ::   Shift_Solution => NULL()

END MODULE Variables_Functions

