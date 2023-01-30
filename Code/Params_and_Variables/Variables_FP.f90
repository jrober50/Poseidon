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

LOGICAL                                                 ::  FP_Diagnostics_Flag = .FALSE.
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:,:)        ::  FP_Update_Vector
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:,:)        ::  FP_Laplace_Vector
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:,:)        ::  FP_Residual_Vector

INTEGER,        ALLOCATABLE,    DIMENSION(:)            ::  FP_Iteration_Log

COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:)          ::  FP_Iter_Matrix_Storage
COMPLEX(idp),   ALLOCATABLE,    DIMENSION(:,:)          ::  FP_Iter_Load_Storage

REAL(idp),      ALLOCATABLE,    DIMENSION(:,:,:,:)      ::  Resid_Norms
REAL(idp),      ALLOCATABLE,    DIMENSION(:,:,:)        ::  Update_Norms

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

