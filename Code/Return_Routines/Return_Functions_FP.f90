   !#################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Return_Functions_FP                                                         !##!
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
USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Poseidon_Numbers_Module, &
            ONLY :  pi

USE Poseidon_Parameters, &
            ONLY :  DEGREE,                 &
                    L_LIMIT,                &
                    NUM_CFA_VARS

USE Parameters_Variable_Indices, &
            ONLY :  iU_CF,                  &
                    iU_LF,                  &
                    iU_S1,                  &
                    iU_S2,                  &
                    iU_S3,                  &
                    iVB_S

USE Variables_Mesh, &
            ONLY :  Num_R_Elements,         &
                    rlocs

USE Variables_Derived, &
            ONLY :  LM_Length,              &
                    ULM_Length

USE Variables_FP, &
            ONLY :  FP_Coeff_Vector_A,      &
                    FP_Coeff_Vector_B

USE Functions_Quadrature, &
            ONLY :  Initialize_LGL_Quadrature_Locations,    &
                    Initialize_LGL_Quadrature

USE Functions_Math, &
            ONLY :  Lagrange_Poly,          &
                    Spherical_Harmonic

USE Maps_Fixed_Point, &
            ONLY :  FP_Array_Map_TypeB, &
                    FP_Array_Map

USE Maps_Domain, &
            ONLY :  Map_To_lm,          &
                    Map_To_FEM_Node


USE Maps_X_Space,    &
            ONLY :  Map_To_X_Space


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

INTEGER                                                         ::  re

REAL(KIND = idp), DIMENSION(0:DEGREE)                           ::  xlocP, weightP


Tmp_U_Value = 0.0_idp


IF ( r == rlocs(0) ) THEN


    LagP    = 0.0_idp
    LagP(0) = 1.0_idp

    CALL Calc_Values_Here_All( 0, theta, phi, LagP, Tmp_U_Value )


ELSE

    CALL Initialize_LGL_Quadrature(DEGREE,xlocP,weightP)

    DO re = 0,NUM_R_ELEMENTS-1
    IF ( r > rlocs(re) .AND. r <= rlocs(re+1) ) THEN


        r_tmp = Map_To_X_Space(rlocs(re),rlocs(re+1),r)
        LagP = Lagrange_Poly(r_tmp,DEGREE,xlocP)

        
        CALL Calc_Values_Here_All( re, theta, phi, LagP, Tmp_U_Value )

        EXIT
    END IF
    END DO

    IF ( r > rlocs(NUM_R_ELEMENTS) ) THEN

        LagP         = 0.0_idp
        LagP(DEGREE) = 1.0_idp

        CALL Calc_Values_Here_All( Num_R_Elements-1, theta, phi, LagP, Tmp_U_Value )

    END IF

END IF



Return_Psi      = REAL(Tmp_U_Value(1), KIND = idp)
Return_AlphaPsi = REAL(Tmp_U_Value(2), KIND = idp)!/REAL(Tmp_U_Value(1), KIND = idp)
Return_Beta1    = REAL(Tmp_U_Value(3), KIND = idp)
Return_Beta2    = REAL(Tmp_U_Value(4), KIND = idp)
Return_Beta3    = REAL(Tmp_U_Value(5), KIND = idp)


END SUBROUTINE Calc_FP_Values_At_Location












!+102+###########################################################################!
!                                                                                !
!                  Calc_1D_CFA_Values_FP          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_1D_CFA_Values_FP(   Num_RE_Input, Num_RQ_Input, RQ_Input,   &
                                    Left_Limit, Right_Limit,                &
                                    CFA_Lapse, CFA_ConFactor, CFA_Shift     )



INTEGER, INTENT(IN)                                         ::  Num_RE_Input,   &
                                                                Num_RQ_Input

REAL(KIND = idp), DIMENSION(1:Num_RQ_Input), INTENT(IN)     ::  RQ_Input
REAL(KIND = idp), INTENT(IN)                                ::  Left_Limit,     &
                                                                Right_Limit


REAL(KIND = idp), DIMENSION(1:NUM_RQ_Input,1:NUM_RE_Input, 1, 1), INTENT(OUT) ::  CFA_Lapse
REAL(KIND = idp), DIMENSION(1:NUM_RQ_Input,1:NUM_RE_Input, 1, 1), INTENT(OUT) ::  CFA_ConFactor
REAL(KIND = idp), DIMENSION(1:NUM_RQ_Input,1:NUM_RE_Input, 1, 1), INTENT(OUT) ::  CFA_Shift



INTEGER                                                         ::  re, x, u, d

REAL(KIND = idp)                                                ::  Quad_Span
REAL(KIND = idp), DIMENSION(0:DEGREE)                           ::  Local_Locations
REAL(KIND = idp), DIMENSION(0:DEGREE)                           ::  LagP
REAL(KIND = idp), DIMENSION(1:Num_RQ_Input)                     ::  CUR_X_LOCS
COMPLEX(KIND = idp), DIMENSION(1:3)                             ::  TMP_U_Value
INTEGER                                                         ::  Current_Location

Quad_Span = Right_Limit - Left_Limit
Local_Locations = Initialize_LGL_Quadrature_Locations(DEGREE)

DO re = 0,NUM_R_ELEMENTS-1

    CUR_X_LOCS = 2.0_idp * ( RQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp

    DO x = 1,Num_RQ_Input
    
        LagP = Lagrange_Poly(CUR_X_LOCS(x),DEGREE,Local_Locations)
        Tmp_U_Value = 0.0_idp

        DO u = 1,2
            DO d = 0,DEGREE

                Current_Location = Map_To_FEM_Node(re,d)
                Tmp_U_Value(u) = Tmp_U_Value(u) + FP_Coeff_Vector_A(Current_Location,1,u)  &
                                                * LagP(d) * Spherical_Harmonic(0,0,pi,pi/2.0_idp)
            END DO ! d Loop
        END DO ! u Loop

        DO d = 0,DEGREE

            Current_Location = FP_Array_Map_TypeB(iU_S1,iVB_S,re,d,0)
            Tmp_U_Value(u) = Tmp_U_Value(u) + FP_Coeff_Vector_B(Current_Location,iVB_S)  &
                                            * LagP(d) * Spherical_Harmonic(0,0,pi,pi/2.0_idp)
        END DO ! d Loop



        CFA_ConFactor(x,re+1,1,1) = REAL(Tmp_U_Value(1), KIND = idp)
        CFA_Lapse(x,re+1,1,1)     = REAL(Tmp_U_Value(2), KIND = idp)        &
                                  / REAL(Tmp_U_Value(1), KIND = idp)
        CFA_Shift(x,re+1,1,1)     = REAL(Tmp_U_Value(3), KIND = idp)
    END DO ! x Loop
END DO ! re Loop


END SUBROUTINE Calc_1D_CFA_Values_FP











 !+501+########################################################!
!                                                               !
!          Calc_Var_At_Location_Type_A                          !
!                                                               !
 !#############################################################!
FUNCTION Calc_Var_At_Location_Type_A( r, theta, phi, iU)

REAL(idp)                                           ::  Calc_Var_At_Location_Type_A
REAL(idp),  INTENT(IN)                              ::  r, theta, phi
INTEGER,    INTENT(IN)                              ::  iU



COMPLEX(KIND = idp)                                 ::  Tmp_U_Value


REAL(KIND = idp)                                    ::  r_tmp
REAL(KIND = idp), DIMENSION(0:DEGREE)               ::  LagP

INTEGER                                             ::  re

REAL(KIND = idp), DIMENSION(0:DEGREE)               ::  xlocP, weightP


Tmp_U_Value = 0.0_idp


IF ( r == rlocs(0) ) THEN

    LagP    = 0.0_idp
    LagP(0) = 1.0_idp

    Tmp_U_Value = Calc_Values_Here_Type_A( 0, theta, phi, LagP, iU )


ELSE

    CALL Initialize_LGL_Quadrature(DEGREE,xlocP,weightP)

    DO re = 0,NUM_R_ELEMENTS-1
    
    IF ( r > rlocs(re) .AND. r <= rlocs(re+1) ) THEN


        r_tmp = Map_To_X_Space(rlocs(re),rlocs(re+1),r)
        LagP = Lagrange_Poly(r_tmp,DEGREE,xlocP)

        
        Tmp_U_Value = Calc_Values_Here_Type_A( re, theta, phi, LagP, iU )

        EXIT
    END IF
    END DO

    IF ( r > rlocs(NUM_R_ELEMENTS) ) THEN

        LagP         = 0.0_idp
        LagP(DEGREE) = 1.0_idp

        Tmp_U_Value = Calc_Values_Here_Type_A( Num_R_Elements-1, theta, phi, LagP, iU )

    END IF

END IF



Calc_Var_At_Location_Type_A = REAL(Tmp_U_Value, KIND = idp)



END FUNCTION Calc_Var_At_Location_Type_A




 !+502+########################################################!
!                                                               !
!          Calc_Var_At_Location_Type_B                          !
!                                                               !
 !#############################################################!
FUNCTION Calc_Var_At_Location_Type_B( r, theta, phi, iU, iVB  )


REAL(idp)                                           ::  Calc_Var_At_Location_Type_B
REAL(idp),  INTENT(IN)                              ::  r, theta, phi
INTEGER,    INTENT(IN)                              ::  iU, iVB



COMPLEX(KIND = idp)                                 ::  Tmp_U_Value


REAL(KIND = idp)                                    ::  r_tmp
REAL(KIND = idp), DIMENSION(0:DEGREE)               ::  LagP

INTEGER                                             ::  re

REAL(KIND = idp), DIMENSION(0:DEGREE)               ::  xlocP, weightP


Tmp_U_Value = 0.0_idp


IF ( r == rlocs(0) ) THEN


    LagP    = 0.0_idp
    LagP(0) = 1.0_idp

    Tmp_U_Value = Calc_Values_Here_Type_B( 0, theta, phi, LagP, iU, iVB )


ELSE

    CALL Initialize_LGL_Quadrature(DEGREE,xlocP,weightP)

    DO re = 0,NUM_R_ELEMENTS-1
    IF ( r > rlocs(re) .AND. r <= rlocs(re+1) ) THEN


        r_tmp = Map_To_X_Space(rlocs(re),rlocs(re+1),r)
        LagP = Lagrange_Poly(r_tmp,DEGREE,xlocP)

        Tmp_U_Value = Calc_Values_Here_Type_B( re, theta, phi, LagP, iU, iVB )


        EXIT
    END IF
    END DO

    IF ( r > rlocs(NUM_R_ELEMENTS) ) THEN

        LagP         = 0.0_idp
        LagP(DEGREE) = 1.0_idp

        Tmp_U_Value = Calc_Values_Here_Type_B( Num_R_Elements-1, theta, phi, LagP, iU, iVB )

    END IF

END IF



Calc_Var_At_Location_Type_B = REAL(Tmp_U_Value, KIND = idp)


END FUNCTION Calc_Var_At_Location_Type_B







!+102+###########################################################################!
!                                                                                !
!                  Calc_1D_CFA_Values_FP          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Values_Here_All( re, theta, phi, LagP, Tmp_U_Value )

INTEGER,                        INTENT(IN)                      ::  re
REAL(idp),                      INTENT(IN)                      ::  theta, phi
REAL(idp), DIMENSION(0:DEGREE), INTENT(INOUT)                   ::  LagP
COMPLEX(idp), DIMENSION(1:5),   INTENT(INOUT)                   ::  Tmp_U_Value


INTEGER                                                         ::  l, m, d, u
INTEGER                                                         ::  Loc_RED, Loc_LM


Tmp_U_Value = 0.0_idp
DO u = iU_CF,iU_LF
DO l = 0,L_Limit
DO m = -l,l
DO d = 0,DEGREE

    Loc_RED = Map_To_FEM_Node(re,d)
    Loc_LM  = Map_to_lm(l,m)


    Tmp_U_Value(u) = Tmp_U_Value(u)                         &
                    + FP_Coeff_Vector_A(Loc_RED,Loc_LM,u)   &
                    * Spherical_Harmonic(l,m,theta,phi)     &
                    * LagP(d)


!    IF ( u == iU_CF) THEN
!        PRINT*,"A",Loc_RED,Loc_LM,re,d
!        PRINT*,FP_Coeff_Vector_A(Loc_RED,Loc_LM,u),  &
!                Spherical_Harmonic(l,m,theta,phi),     &
!                LagP(d)
!    END IF


END DO  !   d Loop
END DO  !   m Loop
END DO  !   l Loop
END DO  !   u Loop


DO u = iU_S1,iU_S3
DO l = 0,L_Limit
DO m = -l,l
DO d = 0,DEGREE

    Loc_RED = FP_Array_Map_TypeB(u,iVB_S,re,d,l,m)
    Tmp_U_Value(u) = Tmp_U_Value(u)                         &
                    + FP_Coeff_Vector_B(Loc_RED,iVB_S)     &
                    * Spherical_Harmonic(l,m,theta,phi)     &
                    * LagP(d)



END DO  !   d Loop
END DO  !   m Loop
END DO  !   l Loop
END DO  !   u Loop


END SUBROUTINE Calc_Values_Here_All





!+102+###########################################################################!
!                                                                                !
!                  Calc_1D_CFA_Values_FP          !
!                                                                                !
!################################################################################!
FUNCTION Calc_Values_Here_Type_A( re, theta, phi, LagP, iU )

REAL(idp)                                               ::  Calc_Values_Here_Type_A
INTEGER,                        INTENT(IN)              ::  re
REAL(idp),                      INTENT(IN)              ::  theta, phi
REAL(idp), DIMENSION(0:DEGREE), INTENT(IN)              ::  LagP
INTEGER,                        INTENT(IN)              ::  iU


COMPLEX(idp)                                            ::  Tmp_U_Value


INTEGER                                                 ::  l, m, d
INTEGER                                                 ::  Loc_RED, Loc_LM


Tmp_U_Value = 0.0_idp
DO l = 0,L_Limit
DO m = -l,l
DO d = 0,DEGREE

    Loc_RED = Map_To_FEM_Node(re,d)
    Loc_LM  = Map_to_lm(l,m)

    Tmp_U_Value = Tmp_U_Value                           &
                + FP_Coeff_Vector_A(Loc_RED,Loc_LM,iU)  &
                * Spherical_Harmonic(l,m,theta,phi)     &
                * LagP(d)
!    IF ( iU == iU_CF) THEN
!        PRINT*,"B",Loc_RED,Loc_LM,re,d
!        PRINT*,FP_Coeff_Vector_A(Loc_RED,Loc_LM,iU),  &
!                Spherical_Harmonic(l,m,theta,phi),     &
!                LagP(d)
!    END IF
END DO  !   d Loop
END DO  !   m Loop
END DO  !   l Loop


Calc_Values_Here_Type_A = REAL( Tmp_U_Value, KIND = idp )

END FUNCTION Calc_Values_Here_Type_A




!+603+###########################################################################!
!                                                                                !
!                  Calc_1D_CFA_Values_FP          !
!                                                                                !
!################################################################################!
FUNCTION Calc_Values_Here_Type_B( re, theta, phi, LagP, iU, iVB )


REAL(idp)                                                       ::  Calc_Values_Here_Type_B
INTEGER,                        INTENT(IN)                      ::  re
REAL(idp),                      INTENT(IN)                      ::  theta, phi
REAL(idp), DIMENSION(0:DEGREE), INTENT(IN)                      ::  LagP
INTEGER,                        INTENT(IN)                      ::  iU
INTEGER,                        INTENT(IN)                      ::  iVB


COMPLEX(idp)                                                    ::  Tmp_U_Value
INTEGER                                                         ::  l, m, d
INTEGER                                                         ::  Loc_RED, Loc_LM


Tmp_U_Value = 0.0_idp
DO l = 0,L_Limit
DO m = -l,l
DO d = 0,DEGREE

    Loc_RED = FP_Array_Map_TypeB(iU,iVB,re,d,l,m)
    Tmp_U_Value = Tmp_U_Value                           &
                + FP_Coeff_Vector_B(Loc_RED,iVB)        &
                * Spherical_Harmonic(l,m,theta,phi)     &
                * LagP(d)


END DO  !   d Loop
END DO  !   m Loop
END DO  !   l Loop


Calc_Values_Here_Type_B = REAL( Tmp_U_Value, KIND = idp)

END FUNCTION Calc_Values_Here_Type_B







!!+203+######################################################################################!
!!                                                                                           !
!!       Poseidon_Return_ConFactor - Uses the results from the routine,                      !
!!                                   Poseidon_XCFC_Run_Part2(), to calculate the value       !
!!                                   of the Conformal Factor at the desired locations.       !
!!                                                                                           !
!!###########################################################################################!
!SUBROUTINE Poseidon_Return_Native_ExtrinsicCurvature(  NE, NQ,                 &
!                                                       RQ_Input,               &
!                                                       TQ_Input,               &
!                                                       PQ_Input,               &
!                                                       Left_Limit,             &
!                                                       Right_Limit,            &
!                                                       Return_Kij              )
!
!INTEGER,    DIMENSION(3),                                   INTENT(IN)  ::  NE, NQ
!REAL(idp),  DIMENSION(NQ(1)),                               INTENT(IN)  ::  RQ_Input
!REAL(idp),  DIMENSION(NQ(2)),                               INTENT(IN)  ::  TQ_Input
!REAL(idp),  DIMENSION(NQ(3)),                               INTENT(IN)  ::  PQ_Input
!REAL(idp),                                                  INTENT(IN)  ::  Left_Limit
!REAL(idp),                                                  INTENT(IN)  ::  Right_Limit
!
!REAL(idp),  DIMENSION(NQ(1)*NQ(2)*NQ(3),NE(1),NE(2),NE(3),1:6), INTENT(OUT) ::  Return_Kij
!
!
!INTEGER                                                         ::  re, te, pe
!INTEGER                                                         ::  rd, td, pd, tpd
!INTEGER                                                         ::  d, lm
!INTEGER                                                         ::  Here, There
!INTEGER                                                         ::  i, j
!
!REAL(idp)                                                       ::  Quad_Span
!REAL(idp), DIMENSION(0:DEGREE)                                  ::  Local_Locations
!REAL(idp), DIMENSION(0:DEGREE)                                  ::  LagP
!REAL(idp), DIMENSION(1:NQ(1))                                   ::  CUR_R_LOCS
!REAL(idp), DIMENSION(1:NQ(1))                                   ::  Cur_RX_Locs
!REAL(idp), DIMENSION(1:NQ(2))                                   ::  CUR_T_LOCS
!REAL(idp), DIMENSION(1:NQ(2))                                   ::  Cur_TX_Locs
!REAL(idp), DIMENSION(1:NQ(3))                                   ::  CUR_P_LOCS
!REAL(idp), DIMENSION(1:NQ(3))                                   ::  Cur_PX_Locs
!
!REAL(idp)                                                       ::  DROT
!REAL(idp)                                                       ::  DTOT
!REAL(idp)                                                       ::  DPOT
!
!REAL(idp), DIMENSION(3)                                         ::  gamma
!REAL(idp), DIMENSION(3,3,3)                                     ::  Christoffel
!COMPLEX(idp), DIMENSION(3)                                      ::  TMP_Val
!COMPLEX(idp), DIMENSION(3,3)                                    ::  TMP_Drv
!
!COMPLEX(idp), DIMENSION(4)                                      ::  Reusable_Vals
!COMPLEX(idp), DIMENSION(6)                                      ::  Tmp_A
!
!
!COMPLEX(idp)                                                    ::  Trace(2)
!
!INTEGER                                                         ::  Current_Location
!INTEGER                                                         ::  iVB
!INTEGER, DIMENSION(3)                                           ::  iU
!
!
!iVB   = iVB_X
!iU(1) = iU_X1
!iU(2) = iU_X2
!iU(3) = iU_X3
!
!Quad_Span = Right_Limit - Left_Limit
!
!Local_Locations = Initialize_LGL_Quadrature_Locations(DEGREE)
!CUR_RX_LOCS = 2.0_idp * ( RQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp
!CUR_TX_LOCS = 2.0_idp * ( TQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp
!CUR_PX_LOCS = 2.0_idp * ( PQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp
!
!! nComp order : ij = 11, 12, 13, 22, 23, 33
!Gamma(1) = 1.0_idp
!Christoffel = 0.0_idp
!
!
!DO pe = 1,NE(3)-1
!DO te = 1,NE(2)-1
!DO re = 1,NE(1)-1
!
!DROT = 0.5_idp * (rlocs(re) - rlocs(re-1))
!DTOT = 0.5_idp * (tlocs(te) - tlocs(te-1))
!DPOT = 0.5_idp * (plocs(pe) - plocs(pe-1))
!
!Cur_R_Locs(:) = DROT * (CUR_RX_LOCS(:)+1.0_idp) + rlocs(re)
!Cur_T_Locs(:) = DTOT * (CUR_TX_LOCS(:)+1.0_idp) + tlocs(te)
!Cur_P_Locs(:) = DPOT * (CUR_PX_LOCS(:)+1.0_idp) + plocs(pe)
!
!DO pd = 1,NQ(3)
!DO td = 1,NQ(2)
!DO rd = 1,NQ(1)
!
!    tpd = Map_To_tpd(td,pd)
!    LagP = Lagrange_Poly(CUR_RX_LOCS(rd),DEGREE,Local_Locations)
!
!    TMP_Val = 0.0_idp
!    TMP_Drv = 0.0_idp
!    DO i = 1,3
!    DO d  = 0,DEGREE
!        Here  = FP_Array_Map_TypeB(iU(i),iVB,re-1,d,1)
!        There = FP_Array_Map_TypeB(iU(i),iVB,re-1,d,LM_Length)
!
!        TMP_Val(i) = TMP_Val(i)                                  &
!                + SUM( FP_Coeff_Vector_B( Here:There, iVB )      &
!                        * Ylm_Values( :, tpd, te-1, pe-1 )   )       &
!                * Lagrange_Poly_Table( d, rd, 0 )
!
!
!        TMP_Drv(1,i) = TMP_Drv(1,i)                              &
!                   + SUM( FP_Coeff_Vector_B( Here:There, iVB )   &
!                         * Ylm_Values( :, tpd, te-1, pe-1 )     )    &
!                   * Lagrange_Poly_Table( d, rd, 1 )             &
!                   / DROT
!
!
!        TMP_Drv(2,i) = TMP_Drv(2,i)                              &
!                   + SUM( FP_Coeff_Vector_B( Here:There, iVB )   &
!                         * Ylm_dt_Values( :, tpd, te-1, pe-1)   )    &
!                   * Lagrange_Poly_Table( d, rd, 0)
!
!        TMP_Drv(3,i) = TMP_Drv(3,i)                              &
!                   + SUM( FP_Coeff_Vector_B( Here:There, iVB )   &
!                         * Ylm_dp_Values( :, tpd, te-1, pe-1)   )    &
!                   * Lagrange_Poly_Table( d, rd, 0)
!
!    END DO  ! d
!    END DO  ! i
!
!
!    Gamma(2) = 1.0_idp/(Cur_R_Locs(rd)*Cur_R_Locs(rd))
!    Gamma(3) = Gamma(2) * 1.0_idp/( DSIN(Cur_T_Locs(td))*DSIN(Cur_T_Locs(td)) )
!
!!    PRINT*,"r",re,rd,Cur_R_Locs(rd)
!!    PRINT*,"Sin, Theta",DSIN(Cur_T_Locs(td)),Cur_T_locs(Td)
!!    PRINT*,"Gamma",1.0_idp/Gamma
!
!    Christoffel(1,2,2) = -Cur_R_Locs(rd)
!    Christoffel(1,3,3) = -Cur_R_Locs(rd)*DSIN(Cur_T_Locs(td))*DSIN(Cur_T_Locs(td))
!    
!    Christoffel(2,1,2) = 1.0_idp/Cur_R_Locs(rd)
!    Christoffel(2,2,1) = 1.0_idp/Cur_R_Locs(rd)
!    Christoffel(2,3,3) = -DSIN(Cur_T_Locs(td))*DCOS(Cur_P_Locs(pd))
!
!    Christoffel(3,3,1) = 1.0_idp/Cur_R_Locs(rd)
!    Christoffel(3,1,3) = 1.0_idp/Cur_R_Locs(rd)
!    Christoffel(3,3,2) = 1.0_idp/DTAN(CUR_T_LOCS(td))
!    Christoffel(3,2,3) = 1.0_idp/DTAN(CUR_T_LOCS(td))
!
!
!    Reusable_Vals(1) = 2.0_idp/3.0_idp*(Tmp_Drv(1,1)+Tmp_Drv(2,2)+Tmp_Drv(3,3))
!    Reusable_Vals(2) = 2.0_idp/3.0_idp*(Christoffel(1,1,1)+Christoffel(2,2,1)+Christoffel(3,3,1))
!    Reusable_Vals(3) = 2.0_idp/3.0_idp*(Christoffel(1,1,2)+Christoffel(2,2,2)+Christoffel(3,3,2))
!    Reusable_Vals(4) = 2.0_idp/3.0_idp*(Christoffel(1,1,3)+Christoffel(2,2,3)+Christoffel(3,3,3))
!
!    
!    Here = rd + (td-1)*NQ(1) + (pd-1)*NQ(1)*NQ(2)
!
!
!    ! Ahat^11
!    Tmp_A(1) = Gamma(1)                                                         &
!             * ( 2.0_idp * Tmp_Drv(1,1) - Reusable_Vals(1)                      &
!               +(2.0_idp * Christoffel(1,1,1) - Reusable_Vals(2) )*Tmp_Val(1)   &
!               +(2.0_idp * Christoffel(1,1,2) - Reusable_Vals(3) )*Tmp_Val(2)   &
!               +(2.0_idp * Christoffel(1,1,3) - Reusable_Vals(4) )*Tmp_Val(3)   )
!
!    Return_Kij(Here,re,te,pe,1) = Tmp_A(1)/(Gamma(1)*Gamma(1))
!
!
!    ! Ahat^12
!    i=1
!    j=2
!    Tmp_A(2) = Gamma(i)*Tmp_Drv(i,j)             &
!             + Gamma(j)*Tmp_Drv(j,i)             &
!             +( Gamma(i)*Christoffel(j,i,1)      &
!              + Gamma(j)*Christoffel(i,j,1)      &
!              )*Tmp_Val(1)                       &
!             +( Gamma(i)*Christoffel(j,i,2)      &
!              + Gamma(j)*Christoffel(i,j,2)      &
!              )*Tmp_Val(2)                       &
!             +( Gamma(i)*Christoffel(j,i,3)      &
!              + Gamma(j)*Christoffel(i,j,3)      &
!              )*Tmp_Val(3)
!     Return_Kij(Here,re,te,pe,2) = Tmp_A(2)/(Gamma(1)*Gamma(2))
!
!    ! Ahat^13
!    i=1
!    j=3
!    Tmp_A(3) = Gamma(i)*Tmp_Drv(i,j)             &
!             + Gamma(j)*Tmp_Drv(j,i)             &
!             +( Gamma(i)*Christoffel(j,i,1)      &
!              + Gamma(j)*Christoffel(i,j,1)      &
!              )*Tmp_Val(1)                       &
!             +( Gamma(i)*Christoffel(j,i,2)      &
!              + Gamma(j)*Christoffel(i,j,2)      &
!              )*Tmp_Val(2)                       &
!             +( Gamma(i)*Christoffel(j,i,3)      &
!              + Gamma(j)*Christoffel(i,j,3)      &
!              )*Tmp_Val(3)
!     Return_Kij(Here,re,te,pe,3) = Tmp_A(3)/(Gamma(1)*Gamma(3))
!
!
!    ! Ahat^22
!    Tmp_A(4) = Gamma(2)                                                             &
!             * ( 2.0_idp * Tmp_Drv(2,2) - Reusable_Vals(1)                          &
!               +(2.0_idp * Christoffel(2,2,1) - Reusable_Vals(2) ) * Tmp_Val(1)     &
!               +(2.0_idp * Christoffel(2,2,2) - Reusable_Vals(3) ) * Tmp_Val(2)     &
!               +(2.0_idp * Christoffel(2,2,3) - Reusable_Vals(4) ) * Tmp_Val(3)     )
!     Return_Kij(Here,re,te,pe,4) = Tmp_A(4)/(Gamma(2)*Gamma(2))
!
!
!    ! Ahat^23
!    i=2
!    j=3
!    Tmp_A(5) = Gamma(i)*Tmp_Drv(i,j)             &
!             + Gamma(j)*Tmp_Drv(j,i)             &
!             +( Gamma(i)*Christoffel(j,i,1)      &
!              + Gamma(j)*Christoffel(i,j,1)      &
!              )*Tmp_Val(1)                       &
!             +( Gamma(i)*Christoffel(j,i,2)      &
!              + Gamma(j)*Christoffel(i,j,2)      &
!              )*Tmp_Val(2)                       &
!             +( Gamma(i)*Christoffel(j,i,3)      &
!              + Gamma(j)*Christoffel(i,j,3)      &
!              )*Tmp_Val(3)
!     Return_Kij(Here,re,te,pe,5) = Tmp_A(5)/(Gamma(2)*Gamma(3))
!
!
!    ! Ahat^33
!    Tmp_A(6) = Gamma(3)                                                             &
!             * ( 2.0_idp * Tmp_Drv(3,3) - Reusable_Vals(1)                          &
!               +(2.0_idp * Christoffel(3,3,1) - Reusable_Vals(2) ) * Tmp_Val(1)     &
!               +(2.0_idp * Christoffel(3,3,2) - Reusable_Vals(3) ) * Tmp_Val(2)     &
!               +(2.0_idp * Christoffel(3,3,3) - Reusable_Vals(4) ) * Tmp_Val(3)     )
!     Return_Kij(Here,re,te,pe,6) = Tmp_A(6)/(Gamma(3)*Gamma(3))
!
!
!!    Return_Kij(Here,re,te,pe,1) =
!!    PRINT*,TMP_A(1),TMP_A(4),TMP_A(6)
!
!!    PRINT*,Return_Kij(Here,re,te,pe,1),          &
!!            Return_Kij(Here,re,te,pe,4),          &
!!            Return_Kij(Here,re,te,pe,6)
!
!
!    Trace(1) = REAL(Tmp_A(1)/Gamma(1)+Tmp_A(4)/Gamma(2)+Tmp_A(6)/Gamma(3),KIND = idp)
!
!    Trace(2) = Trace(1)                                     &
!             /REAL(3.0_idp*(  Reusable_Vals(1)                  &
!                        + Reusable_Vals(2)*Tmp_Val(1)       &
!                        + Reusable_Vals(3)*Tmp_Val(2)       &
!                        + Reusable_Vals(4)*Tmp_Val(3)  ), Kind = idp    )
!
!
!
!
!
!END DO ! rd Loop
!END DO ! td Loop
!END DO ! pd Loop
!END DO ! re Loop
!END DO ! te Loop
!END DO ! pe Loop
!
!
!
!
!END SUBROUTINE Poseidon_Return_Native_ExtrinsicCurvature


END MODULE Return_Functions_FP

