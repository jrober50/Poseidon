   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Functions_Translation_Matrix_Module                                   !##!
!##!                                                                         !##!
!##!_________________________________________________________________________!##!
!##!                                                                         !##!
!##!                                                                         !##!
!##!=========================================================================!##!
!##!                                                                         !##!
!##!    Contains:                                                            !##!
!##!                                                                         !##!
!##!                                                                         !##!
!###############################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !##########################################################################!


!*D*================================!
!                                   !
!           Dependencies            !
!                                   !
!===================================!
USE Poseidon_Kinds_Module, &
            ONLY : idp


USE Maps_X_Space, &
            ONLY :  Map_To_X_Space

USE Functions_Math, &
            ONLY :  Lagrange_Poly

IMPLICIT NONE


CONTAINS



FUNCTION Create_Translation_Matrix( Source_NQ,          &
                                    Source_xL,          &
                                    Source_RQ_xlocs,    &
                                    Source_TQ_xlocs,    &
                                    Source_PQ_xlocs,    &
                                    Source_DOF,         &
                                    Dest_NQ,            &
                                    Dest_xL,            &
                                    Dest_RQ_xlocs,      &
                                    Dest_TQ_xlocs,      &
                                    Dest_PQ_xlocs,      &
                                    Dest_DOF            ) RESULT( TransMat )


INTEGER,    DIMENSION(3),               INTENT(IN)      ::  Source_NQ
REAL(idp),  DIMENSION(2),               INTENT(IN)      ::  Source_xL
REAL(idp),  DIMENSION(1:Source_NQ(1)),  INTENT(IN)      ::  Source_RQ_xlocs
REAL(idp),  DIMENSION(1:Source_NQ(2)),  INTENT(IN)      ::  Source_TQ_xlocs
REAL(idp),  DIMENSION(1:Source_NQ(3)),  INTENT(IN)      ::  Source_PQ_xlocs
INTEGER,                                INTENT(IN)      ::  Source_DOF

INTEGER,    DIMENSION(3),               INTENT(IN)      ::  Dest_NQ
REAL(idp),  DIMENSION(2),               INTENT(IN)      ::  Dest_xL
REAL(idp),  DIMENSION(1:Dest_NQ(1)),    INTENT(IN)      ::  Dest_RQ_xlocs
REAL(idp),  DIMENSION(1:Dest_NQ(2)),    INTENT(IN)      ::  Dest_TQ_xlocs
REAL(idp),  DIMENSION(1:Dest_NQ(3)),    INTENT(IN)      ::  Dest_PQ_xlocs
INTEGER,                                INTENT(IN)      ::  Dest_DOF

REAL(idp),  DIMENSION(1:Source_DOF,1:Dest_DOF)          ::  TransMat

INTEGER                                             ::  Dest_R
INTEGER                                             ::  Dest_T
INTEGER                                             ::  Dest_P

INTEGER                                             ::  Source_T
INTEGER                                             ::  Source_P

INTEGER                                             ::  Here
INTEGER                                             ::  There
INTEGER                                             ::  Dest_Here

REAL(idp),  DIMENSION( 1:Source_NQ(1) )              ::  Scaled_R_Quad
REAL(idp),  DIMENSION( 1:Source_NQ(2) )              ::  Scaled_T_Quad
REAL(idp),  DIMENSION( 1:Source_NQ(3) )              ::  Scaled_P_Quad

REAL(idp), DIMENSION(:,:), ALLOCATABLE              ::  R_Lag_Poly_Values
REAL(idp), DIMENSION(:,:), ALLOCATABLE              ::  T_Lag_Poly_Values
REAL(idp), DIMENSION(:,:), ALLOCATABLE              ::  P_Lag_Poly_Values



ALLOCATE( R_Lag_Poly_Values(1:Source_NQ(1),1:Dest_NQ(1)) )
ALLOCATE( T_Lag_Poly_Values(1:Source_NQ(2),1:Dest_NQ(2)) )
ALLOCATE( P_Lag_Poly_Values(1:Source_NQ(3),1:Dest_NQ(3)) )


Scaled_R_Quad = Map_To_X_Space(Source_xL(1),Source_xL(2),Source_RQ_xlocs)
Scaled_T_Quad = Map_To_X_Space(Source_xL(1),Source_xL(2),Source_TQ_xlocs)
Scaled_P_Quad = Map_To_X_Space(Source_xL(1),Source_xL(2),Source_PQ_xlocs)




DO Dest_R = 1,Dest_NQ(1)
 R_Lag_Poly_Values(:,Dest_R) = Lagrange_Poly(Dest_RQ_xlocs(Dest_R),  &
                                             Source_NQ(1)-1,         &
                                             Scaled_R_Quad           )

END DO

DO Dest_T = 1,Dest_NQ(2)
 T_Lag_Poly_Values(:,Dest_T) = Lagrange_Poly(Dest_TQ_xlocs(Dest_T),  &
                                             Source_NQ(2)-1,         &
                                             Scaled_T_Quad           )

END DO

DO Dest_P = 1,Dest_NQ(3)
 P_Lag_Poly_Values(:,Dest_P) = Lagrange_Poly(Dest_PQ_xlocs(Dest_P),  &
                                             Source_NQ(3)-1,         &
                                             Scaled_P_Quad           )

END DO



DO Dest_P = 1,Dest_NQ(3)
DO Dest_T = 1,Dest_NQ(2)
DO Dest_R = 1,Dest_NQ(1)

Dest_Here = (Dest_P-1) * Dest_NQ(2) * Dest_NQ(1)      &
          + (Dest_T-1) * Dest_NQ(1)                   &
          + Dest_R

DO Source_P = 1,Source_NQ(3)
DO Source_T = 1,Source_NQ(2)

    Here = (Source_P-1) * Source_NQ(2) * Source_NQ(1)   &
         + (Source_T-1) * Source_NQ(1)

    There = Here + Source_NQ(1)

    TransMat(Here+1:There, Dest_Here)  =                    &
                R_Lag_Poly_Values(1:Source_NQ(1),Dest_R)    &
              * T_Lag_Poly_Values(Source_T,Source_T)        &
              * P_Lag_Poly_Values(Source_P,Source_P)

END DO  !   Source_T Loop
END DO  !   Source_P Loop
END DO  !   Dest_R Loop
END DO  !   Dest_T Loop
END DO  !   Dest_P Loop



DEALLOCATE( R_Lag_Poly_Values )
DEALLOCATE( T_Lag_Poly_Values )
DEALLOCATE( P_Lag_Poly_Values )



END FUNCTION Create_Translation_Matrix




END MODULE Functions_Translation_Matrix_Module
