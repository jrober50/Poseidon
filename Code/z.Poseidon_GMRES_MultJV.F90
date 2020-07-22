   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE GMRES_MultJV_Module                                                          !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!

USE Poseidon_Constants_Module, &
                ONLY :  idp

USE GMRES_Functions_Module, &
                ONLY :  CNORM2

USE GMRES_Master_Build_Module, &
                ONLY :  Calc_Equation

USE Poseidon_Variables_Module, &
                ONLY :  Prob_Dim,           &
                        Coefficient_Vector


IMPLICIT NONE

CONTAINS



 !+101+############################################################################!
!                                                                                   !
!                       GMRES_MultJV                                                !
!                                                                                   !
 !#################################################################################!
SUBROUTINE GMRES_MultJV(JV, V, Coeff_Vec, Coeff_Length)

COMPLEX(KIND = idp), DIMENSION(0:Coeff_Length-1), INTENT(IN)    :: V
COMPLEX(KIND = idp), DIMENSION(0:Coeff_Length-1), INTENT(IN)    :: Coeff_Vec
INTEGER, INTENT(IN)                                             :: Coeff_Length

COMPLEX(KIND = idp), DIMENSION(0:Coeff_Length-1), INTENT(OUT)   ::  JV

COMPLEX(KIND = idp), DIMENSION(0:Coeff_Length-1)                :: Eq_a
COMPLEX(KIND = idp), DIMENSION(0:Coeff_Length-1)                :: Eq_b
COMPLEX(KIND = idp), DIMENSION(0:Coeff_Length-1)                :: Coeff_Perturbed

INTEGER                                                         :: Pert_Flag = 5
COMPLEX(KIND = idp)                                             :: Perturbation


CALL Calculate_Perturbation(V, Coeff_Vec, Coeff_Length, Pert_Flag, Perturbation)

Coeff_Perturbed = Coeff_Vec + Perturbation*V

CALL Calc_Equation( Eq_a, Coeff_Vec, Prob_Dim )
CALL Calc_Equation( Eq_b, Coeff_Perturbed, Prob_Dim )

!JV = -(1.0_idp/Perturbation) * ( Equation(x + Perturbation*V) - Equation(x) )
!JV = ( Eq_a - Eq_b )/Perturbation
JV = 0.0_idp

END SUBROUTINE GMRES_MultJV







 !+102+############################################################################!
!                                                                                   !
!                       Calculate_Perturbation                                      !
!                                                                                   !
 !#################################################################################!
SUBROUTINE Calculate_Perturbation(V, x, Length, Flag, Pert)

COMPLEX(KIND = idp), DIMENSION(:), INTENT(IN)       :: V
COMPLEX(KIND = idp), DIMENSION(:), INTENT(IN)       :: x
INTEGER, INTENT(IN)                                 :: Length
INTEGER, INTENT(IN)                                 :: Flag

COMPLEX(KIND = idp), INTENT(OUT)                    :: Pert

INTEGER                                             :: i
REAL(KIND = idp)                                    :: eps
COMPLEX(KIND = idp)                                 :: tmp
INTEGER                                             :: dim

COMPLEX(KIND = idp)                                 :: normx
COMPLEX(KIND = idp)                                 :: normv
COMPLEX(KIND = idp)                                 :: b, a
COMPLEX(KIND = idp)                                 :: xv

IF ( Flag == 1  ) THEN

    Pert = 1E-2

ELSEIF ( Flag == 2 ) THEN

    dim = size(x)
    eps = epsilon(eps)

    tmp = SUM(sqrt(eps)*(1+x(:)))

    IF ( cnorm2(V) > eps ) THEN
        Pert = (1.0_idp/(dim * cnorm2(V)))*tmp
    ELSE
        Pert = tmp/dim
    END IF

ELSE IF (Flag == 3) THEN

    PRINT*,"cnorm2(V)",cnorm2(V)
    Pert = cnorm2(x)/cnorm2(V)*10E-6


ELSE IF (Flag == 4 ) THEN

    normv = cnorm2(V)
    IF ( normv == 0.0_idp ) THEN
        normv = 1.0_idp
    END IF
    dim = size(V)
    b = sqrt(epsilon(eps))
    
    PRINT*,dim, normv,SUM(b*x),b
    Pert = 1.0_idp/(dim*normv) * SUM(b*x) + b

ELSE IF (Flag == 5 ) THEN

    normx = cnorm2(x)
    normv = cnorm2(V)
    IF ( normv == 0.0_idp ) THEN
        normv = 1.0_idp
    END IF
    
    Pert = sqrt((1+normx)*epsilon(eps))/normv


ELSE IF (Flag == 6 ) THEN

    normx = cnorm2(x)
    normv = cnorm2(V)
    xv = x(1)*v(1)+x(2)*v(2)
    b = sqrt(epsilon(eps))
    a = 1.0_idp

!    Pert = b/normv * max(xv,normx)

END IF


END SUBROUTINE Calculate_Perturbation


END MODULE GMRES_MultJV_Module
