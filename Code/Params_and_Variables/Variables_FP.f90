   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Variables_FP                                                                 !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
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
            ONLY : idp



IMPLICIT NONE


COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:,:)        ::  FP_Update_Vector
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:,:)        ::  FP_Laplace_Vector
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:,:)        ::  FP_Residual_Vector

COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:)            ::  FP_Laplace_Vector_Beta
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:)            ::  FP_Residual_Vector_Beta

COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:)            ::  FP_Laplace_Vector_X
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:)            ::  FP_Residual_Vector_X






! Cholesky Factorizatio Flag
INTEGER, PUBLIC                                         ::  FP_Anderson_M = 3
INTEGER, PUBLIC, PARAMETER                              ::  FP_Anderson_M_Default = 3



INTEGER, PARAMETER                  :: N_FPTT         = 17
CHARACTER(LEN=31), DIMENSION(N_FPTT):: FPTT_Names = [   'Initalize Fixed Point Matrices ',   &
                                                        '2                              ',   &
                                                        'Calculate FP Source Vector     ',   &
                                                        'Solve Laplacian System(s)      ',   &
                                                        'Solve MVL System(s)            ',   &
                                                        '6                              ',   &
                                                        '7                              ',   &
                                                        'Average Fixed Point Iteration  ',   &
                                                        '9                              ',   &
                                                        '10                             ',   &
                                                        'Cholesky Factorization         ',   &
                                                        'Factorize Shift Matrix         ',   &
                                                        '13                             ',   &
                                                        '14                             ',   &
                                                        '15                             ',   &
                                                        'Solve Laplacian Systems        ',   &
                                                        'Solve MVL Systems              '   ]



END MODULE Variables_FP 

