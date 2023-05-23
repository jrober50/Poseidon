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
