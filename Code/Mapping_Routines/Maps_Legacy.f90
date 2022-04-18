   !#################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Maps_Legacy                                                                  !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +201+   CFA_Matrix_Map                                                      !##!
!##!    +202+   CFA_1D_Matrix_Map                                                   !##!
!##!    +203+   CFA_2D_Matrix_Map                                                   !##!
!##!    +204+   CFA_3D_Matrix_Map                                                   !##!
!##!    +205+   CFA_1D_LM_Map                                                       !##!
!##!    +206+   CFA_2D_LM_Map                                                       !##!
!##!    +207+   CFA_3D_LM_Map                                                       !##!
!##!                                                                                !##!
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
            ONLY : idp

USE Poseidon_Numbers_Module, &
            ONLY : pi

USE Poseidon_Parameters, &
            ONLY :  DOMAIN_DIM,             &
                    DEGREE,                 &
                    L_LIMIT,                &
                    NUM_CFA_VARS

USE Variables_Derived, &
            ONLY :  ULM_LENGTH,             &
                    LM_LENGTH

IMPLICIT NONE


CONTAINS






!+201+###########################################################################!
!                                                                                !
!                  CFA_Matrix_Map                                                  !
!                                                                                !
!################################################################################!
PURE FUNCTION CFA_Matrix_Map( ui, l, m, re, d )

INTEGER                                     :: CFA_Matrix_Map

INTEGER, INTENT(IN)                         :: ui, l, m, re, d





CFA_Matrix_Map = (re*DEGREE + d)* ((2+DOMAIN_DIM)*(L_LIMIT+1)*(L_LIMIT+1))  &
                  + (l*(l+1) + m)*(2+DOMAIN_DIM)                            &
                  + (ui - 1)


END FUNCTION CFA_Matrix_Map





!+104+##########################################################################!
!                                                                               !
!                  CFA_1D_Matrix_Map                                            !
!                                                                               !
!###############################################################################!
PURE FUNCTION CFA_1D_Matrix_Map( ui, l ,m, re, d )

INTEGER                                     :: CFA_1D_Matrix_Map
INTEGER, INTENT(IN)                         :: ui, l, m, re, d


CFA_1D_Matrix_Map = (re*DEGREE + d)*NUM_CFA_VARS                &
                  + (ui - 1)




END FUNCTION CFA_1D_Matrix_Map







!+103+##########################################################################!
!                                                                               !
!                  CFA_2D_Matrix_Map                                            !
!                                                                               !
!###############################################################################!
PURE FUNCTION CFA_2D_Matrix_Map( ui, l, m, re, d )

INTEGER                                     :: CFA_2D_Matrix_Map
INTEGER, INTENT(IN)                         :: ui, l, m, re, d


CFA_2D_Matrix_Map = (re*DEGREE + d) * ULM_LENGTH            &
                  + l * NUM_CFA_VARS                        &
                  + (ui - 1)


END FUNCTION CFA_2D_Matrix_Map








!+104+##########################################################################!
!                                                                               !
!                  CFA_3D_Matrix_Map                                            !
!                                                                               !
!###############################################################################!
PURE FUNCTION CFA_3D_Matrix_Map( ui, l, m, re, d )

INTEGER                                     :: CFA_3D_Matrix_Map
INTEGER, INTENT(IN)                         :: ui, l, m, re, d


CFA_3D_Matrix_Map = (re*DEGREE + d) * ULM_LENGTH            &
                  + (ui - 1) * LM_LENGTH                    &
                  + (l*(l+1) + m)


!CFA_3D_Matrix_Map = (re*DEGREE + d) * ULM_LENGTH            &
!                  + (l*(l+1) + m) * NUM_CFA_VARS            &
!                  + (ui - 1)

END FUNCTION CFA_3D_Matrix_Map








!+105+##########################################################################!
!                                                                               !
!                  CFA_1D_LM_Map                                            !
!                                                                               !
!###############################################################################!
PURE FUNCTION CFA_1D_LM_Map( l , m )

INTEGER                                     :: CFA_1D_LM_Map
INTEGER, INTENT(IN)                         :: l, m

CFA_1D_LM_Map = 0

END FUNCTION CFA_1D_LM_Map





!+106+##########################################################################!
!                                                                               !
!                  CFA_2D_LM_Map                                            !
!                                                                               !
!###############################################################################!
PURE FUNCTION CFA_2D_LM_Map( l , m )

INTEGER                                     :: CFA_2D_LM_Map
INTEGER, INTENT(IN)                         :: l, m

CFA_2D_LM_Map = l + 1

END FUNCTION CFA_2D_LM_Map



!+107+##########################################################################!
!                                                                               !
!                  CFA_3D_LM_Map                                            !
!                                                                               !
!###############################################################################!
PURE FUNCTION CFA_3D_LM_Map( l , m )

INTEGER                                     :: CFA_3D_LM_Map
INTEGER, INTENT(IN)                         :: l, m

CFA_3D_LM_Map = l*(l+1) + m + 1

END FUNCTION CFA_3D_LM_Map






!+102+##########################################################################!
!                                                                               !
!                  CFA_3D_Matrix_Map                                            !
!                                                                               !
!###############################################################################!
PURE FUNCTION CFA_ALL_Matrix_Map( ui, lm_loc, re, d )

INTEGER                                     :: CFA_ALL_Matrix_Map
INTEGER, INTENT(IN)                         :: ui, lm_loc, re, d



CFA_ALL_Matrix_Map = (re*DEGREE + d) * ULM_LENGTH           &
                   + (ui - 1) * LM_LENGTH                   &
                   + lm_loc


END FUNCTION CFA_ALL_Matrix_Map








END MODULE Maps_Legacy

