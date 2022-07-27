   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Variables_BC                                                                 !##!
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




CHARACTER(LEN = 1), DIMENSION(1:5)              :: INNER_CFA_BC_TYPE
REAL(KIND = idp), DIMENSION(1:5)                :: INNER_CFA_BC_VALUES


CHARACTER(LEN = 1), DIMENSION(1:5)              :: OUTER_CFA_BC_TYPE
REAL(KIND = idp), DIMENSION(1:5)                :: OUTER_CFA_BC_VALUES




CHARACTER(LEN = 1)                              :: INNER_Poisson_BC_TYPE
REAL(KIND = idp)                                :: INNER_Poisson_BC_VALUE

CHARACTER(LEN = 1)                              :: OUTER_Poisson_BC_TYPE
REAL(KIND = idp)                                :: OUTER_Poisson_BC_VALUE


END MODULE Variables_BC



