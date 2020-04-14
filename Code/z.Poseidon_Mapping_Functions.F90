   !#################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_Mapping_Functions_Module                                            !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   Map_To_X_Space                                                      !##!
!##!    +102+   Map_From_X_Space                                                    !##!
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
USE Poseidon_Constants_Module, &
            ONLY :  idp, pi, eps


USE Poseidon_Parameters, &
            ONLY :  DOMAIN_DIM,             &
                    DEGREE,                 &
                    L_LIMIT,                &
                    NUM_CFA_VARS

USE Poseidon_Variables_Module, &
            ONLY :  ULM_LENGTH,             &
                    LM_LENGTH

IMPLICIT NONE


CONTAINS




!+101+##########################################################!
!                                                               !
!      Map_To_X_Space - maps r value between ra, and rb to x    !
!                   space such that x in [-1,1].                !
!                                                               !
!###############################################################!
PURE ELEMENTAL FUNCTION Map_To_X_Space(ra, rb, r)

REAL(KIND = idp)                            ::  Map_To_X_Space
REAL(KIND = idp), intent(in)                ::  ra, rb
REAL(KIND = idp), intent(in)                ::  r

Map_To_X_Space = (2.0_idp*(r - ra))/(rb - ra) - 1.0_idp


END FUNCTION Map_To_X_Space







!+102+##########################################################!
!                                                               !
!      Map_From_X_Space - maps x value between -1, and 1 to r   !
!                   space such that r in [ra,rb].               !
!                                                               !
!###############################################################!
PURE ELEMENTAL FUNCTION Map_From_X_Space(ra, rb, x)

REAL(KIND = idp)                            ::  Map_From_X_Space
REAL(KIND = idp), intent(in)                ::  ra, rb
REAL(KIND = idp), intent(in)                ::  x

Map_From_X_Space = ((rb - ra) * (x+1.0_idp))/2.0_idp  + ra


END FUNCTION Map_From_X_Space





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


CFA_1D_Matrix_Map = NUM_CFA_VARS*(re*DEGREE + d)            &
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

CFA_2D_LM_Map = l

END FUNCTION CFA_2D_LM_Map



!+107+##########################################################################!
!                                                                               !
!                  CFA_3D_LM_Map                                            !
!                                                                               !
!###############################################################################!
PURE FUNCTION CFA_3D_LM_Map( l , m )

INTEGER                                     :: CFA_3D_LM_Map
INTEGER, INTENT(IN)                         :: l, m

CFA_3D_LM_Map = l*(l+1) + m

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









END MODULE Poseidon_Mapping_Functions_Module

