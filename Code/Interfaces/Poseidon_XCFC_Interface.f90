   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_XCFC_Interface_Module                                               !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains the top level subroutines needed to run the Poseidon in XCFC mode. !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+                                                                       !##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!

USE Poseidon_Kinds_Module, &
            ONLY : idp

USE Poseidon_Parameters, &
            ONLY :  DEGREE

USE Parameters_Variable_Indices, &
            ONLY :  iVB_X,                      &
                    iVB_S,                      &
                    iU_CF,                      &
                    iU_LF,                      &
                    iU_S1,                      &
                    iU_S2,                      &
                    iU_S3,                      &
                    iU_X1,                      &
                    iU_X2,                      &
                    iU_X3


USE Variables_Tables, &
            ONLY :  Ylm_Values,                 &
                    Ylm_dt_Values,              &
                    Ylm_dp_Values,              &
                    Lagrange_Poly_Table


USE Variables_Derived, &
            ONLY :  LM_LENGTH,          &
                    Beta_Prob_Dim

USE Variables_FP, &
            ONLY :  FP_Coeff_Vector_A,  &
                    FP_Coeff_Vector_B

USE Variables_Mesh, &
        ONLY :  rlocs,              &
                tlocs,              &
                plocs

USE Variables_Mesh, &
        ONLY :  rlocs,              &
                tlocs,              &
                plocs

USE FP_Functions_Mapping, &
        ONLY :  FP_Array_Map_TypeB

USE Functions_Domain_Maps, &
        ONLY :  Map_To_tpd


USE Functions_Domain_Maps, &
            ONLY :  Map_To_tpd,     &
                    Map_To_FEM_Node

USE Functions_Quadrature, &
            ONLY :  Initialize_LGL_Quadrature_Locations

USE Functions_Math, &
            ONLY :  Lagrange_Poly


USE XCFC_Method_Module, &
            ONLY :  XCFC_Method_Part1,      &
                    XCFC_Method_Part2


CONTAINS


!+101+######################################################################################!
!                                                                                           !
!       Poseidon_XCFC_Run_Part1 - Performs the steps to calculate the XCFC Conformal        !
!                                 Factor.  Upon completion of this call, the routine        !
!                                 Poseidon_Return_ConFactor() can be used to acquire        !
!                                 values of the conformal factor at desired locations.      !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_XCFC_Run_Part1()


CALL XCFC_Method_Part1()


END SUBROUTINE Poseidon_XCFC_Run_Part1





!+102+######################################################################################!
!                                                                                           !
!       Poseidon_XCFC_Run_Part2 - Performs the steps to calculate the XCFC Lapse Function,  !
!                                 and Shift Vector.  Upon completion of this call, the      !
!                                 routine, Poseidon_Return_Lapse() and                      !
!                                 Poseidon_Return_Shift() can be used to acquire            !
!                                 values of the those variables at desired locations.       !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_XCFC_Run_Part2()


CALL XCFC_Method_Part2()


END SUBROUTINE Poseidon_XCFC_Run_Part2







!+201+######################################################################################!
!                                                                                           !
!       Poseidon_Return_ConFactor - Uses the results from the routine,                      !
!                                   Poseidon_XCFC_Run_Part1(), to calculate the value       !
!                                   of the Conformal Factor at the desired locations.       !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Return_ConFactor(   NE, NQ,                 &
                                        RQ_Input,           &
                                        TQ_Input,           &
                                        PQ_Input,           &
                                        Left_Limit,             &
                                        Right_Limit,            &
                                        Return_ConFactor        )


INTEGER,    DIMENSION(3),                                   INTENT(IN)  ::  NE, NQ
REAL(idp),  DIMENSION(NQ(1)),                               INTENT(IN)  ::  RQ_Input
REAL(idp),  DIMENSION(NQ(2)),                               INTENT(IN)  ::  TQ_Input
REAL(idp),  DIMENSION(NQ(3)),                               INTENT(IN)  ::  PQ_Input
REAL(idp),                                                  INTENT(IN)  ::  Left_Limit
REAL(idp),                                                  INTENT(IN)  ::  Right_Limit

REAL(idp),  DIMENSION(NQ(1)*NQ(2)*NQ(3),NE(1),NE(2),NE(3)), INTENT(OUT) ::  Return_Confactor



INTEGER                                                         ::  re, te, pe
INTEGER                                                         ::  rd, td, pd, tpd
INTEGER                                                         ::  d, lm, Here

REAL(KIND = idp)                                                ::  Quad_Span
REAL(KIND = idp), DIMENSION(0:DEGREE)                           ::  Local_Locations
REAL(KIND = idp), DIMENSION(0:DEGREE)                           ::  LagP
REAL(KIND = idp), DIMENSION(1:NQ(1))                            ::  CUR_X_LOCS
COMPLEX(KIND = idp)                                             ::  TMP_U_Value
INTEGER                                                         ::  Current_Location

Quad_Span = Right_Limit - Left_Limit

Local_Locations = Initialize_LGL_Quadrature_Locations(DEGREE)
CUR_X_LOCS = 2.0_idp * ( RQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp


DO pe = 1,NE(3)
DO te = 1,NE(2)
DO re = 1,NE(1)


    DO pd = 1,NQ(3)
    DO td = 1,NQ(2)
    DO rd = 1,NQ(1)

        tpd = Map_To_tpd(td,pd)
        LagP = Lagrange_Poly(CUR_X_LOCS(rd),DEGREE,Local_Locations)
        Tmp_U_Value = 0.0_idp

        DO lm = 1,LM_Length
        DO d = 0,DEGREE

            Current_Location = Map_To_FEM_Node(re-1,d)
            Tmp_U_Value = Tmp_U_Value + FP_Coeff_Vector_A(Current_Location,lm,iU_CF)  &
                                      * LagP(d) * Ylm_Values( lm, tpd, te-1, pe-1 )

        END DO ! d Loop
        END DO ! lm Loop

        Here = rd + (td-1)*NQ(1) + (pd-1)*NQ(1)*NQ(2)
        Return_ConFactor(Here,re,te,pe) = REAL(Tmp_U_Value, KIND = idp)

    END DO ! rd Loop
    END DO ! td Loop
    END DO ! pd Loop

END DO ! re Loop
END DO ! te Loop
END DO ! pe Loop


END SUBROUTINE Poseidon_Return_ConFactor






!+202+######################################################################################!
!                                                                                           !
!       Poseidon_Return_ConFactor - Uses the results from the routine,                      !
!                                   Poseidon_XCFC_Run_Part2(), to calculate the value       !
!                                   of the Conformal Factor at the desired locations.       !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Return_Lapse(   NE, NQ,                 &
                                    RQ_Input,           &
                                    TQ_Input,           &
                                    PQ_Input,           &
                                    Left_Limit,             &
                                    Right_Limit,            &
                                    Return_Lapse            )


INTEGER,    DIMENSION(3),                                   INTENT(IN)  ::  NE, NQ
REAL(idp),  DIMENSION(NQ(1)),                               INTENT(IN)  ::  RQ_Input
REAL(idp),  DIMENSION(NQ(2)),                               INTENT(IN)  ::  TQ_Input
REAL(idp),  DIMENSION(NQ(3)),                               INTENT(IN)  ::  PQ_Input
REAL(idp),                                                  INTENT(IN)  ::  Left_Limit
REAL(idp),                                                  INTENT(IN)  ::  Right_Limit

REAL(idp),  DIMENSION(NQ(1)*NQ(2)*NQ(3),NE(1),NE(2),NE(3)), INTENT(OUT) ::  Return_Lapse



INTEGER                                                         ::  re, te, pe
INTEGER                                                         ::  rd, td, pd, tpd
INTEGER                                                         ::  d, lm, Here

REAL(KIND = idp)                                                ::  Quad_Span
REAL(KIND = idp), DIMENSION(0:DEGREE)                           ::  Local_Locations
REAL(KIND = idp), DIMENSION(0:DEGREE)                           ::  LagP
REAL(KIND = idp), DIMENSION(1:NQ(1) )                           ::  CUR_X_LOCS
COMPLEX(KIND = idp), DIMENSION(1:2)                             ::  TMP_U_Value
INTEGER                                                         ::  Current_Location

Quad_Span = Right_Limit - Left_Limit

Local_Locations = Initialize_LGL_Quadrature_Locations(DEGREE)
CUR_X_LOCS = 2.0_idp * ( RQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp


DO pe = 1,NE(3)
DO te = 1,NE(2)
DO re = 1,NE(1)


    DO pd = 1,NQ(3)
    DO td = 1,NQ(2)
    DO rd = 1,NQ(1)

        tpd = Map_To_tpd(td,pd)
        LagP = Lagrange_Poly(CUR_X_LOCS(rd),DEGREE,Local_Locations)
        Tmp_U_Value = 0.0_idp

        DO lm = 1,LM_Length
        DO d = 0,DEGREE

            Current_Location = Map_To_FEM_Node(re-1,d)
            Tmp_U_Value(1) = Tmp_U_Value(1) + FP_Coeff_Vector_A(Current_Location,lm,iU_CF)  &
                                            * LagP(d) * Ylm_Values( lm, tpd, te-1, pe-1 )

            Tmp_U_Value(2) = Tmp_U_Value(2) + FP_Coeff_Vector_A(Current_Location,lm,iU_LF)  &
                                            * LagP(d) * Ylm_Values( lm, tpd, te-1, pe-1 )

        END DO ! d Loop
        END DO ! lm Loop

        Here = rd + (td-1)*NQ(1) + (pd-1)*NQ(1)*NQ(2)
        Return_Lapse(Here,re,te,pe) = REAL(Tmp_U_Value(2), KIND = idp)     &
                                    / REAL(Tmp_U_Value(1), KIND = idp)

    END DO ! rd Loop
    END DO ! td Loop
    END DO ! pd Loop

END DO ! re Loop
END DO ! te Loop
END DO ! pe Loop




END SUBROUTINE Poseidon_Return_Lapse





!+203+######################################################################################!
!                                                                                           !
!       Poseidon_Return_ConFactor - Uses the results from the routine,                      !
!                                   Poseidon_XCFC_Run_Part2(), to calculate the value       !
!                                   of the Conformal Factor at the desired locations.       !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Return_Shift(   NE, NQ,                 &
                                    RQ_Input,           &
                                    TQ_Input,           &
                                    PQ_Input,           &
                                    Left_Limit,             &
                                    Right_Limit,            &
                                    Return_Shift            )

INTEGER,    DIMENSION(3),                                   INTENT(IN)  ::  NE, NQ
REAL(idp),  DIMENSION(NQ(1)),                               INTENT(IN)  ::  RQ_Input
REAL(idp),  DIMENSION(NQ(2)),                               INTENT(IN)  ::  TQ_Input
REAL(idp),  DIMENSION(NQ(3)),                               INTENT(IN)  ::  PQ_Input
REAL(idp),                                                  INTENT(IN)  ::  Left_Limit
REAL(idp),                                                  INTENT(IN)  ::  Right_Limit

REAL(idp),  DIMENSION(NQ(1)*NQ(2)*NQ(3),NE(1),NE(2),NE(3),1:3), INTENT(OUT) ::  Return_Shift


INTEGER                                                         ::  re, te, pe
INTEGER                                                         ::  rd, td, pd, tpd
INTEGER                                                         ::  d, lm, Here, iU

REAL(KIND = idp)                                                ::  Quad_Span
REAL(KIND = idp), DIMENSION(0:DEGREE)                           ::  Local_Locations
REAL(KIND = idp), DIMENSION(0:DEGREE)                           ::  LagP
REAL(KIND = idp), DIMENSION(1:NQ(1))                            ::  CUR_X_LOCS
COMPLEX(KIND = idp)                                             ::  TMP_U_Value
INTEGER                                                         ::  Current_Location

Quad_Span = Right_Limit - Left_Limit

Local_Locations = Initialize_LGL_Quadrature_Locations(DEGREE)
CUR_X_LOCS = 2.0_idp * ( RQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp


DO iU = iU_S1, iU_S3
DO pe = 1,NE(3)
DO te = 1,NE(2)
DO re = 1,NE(1)
DO pd = 1,NQ(3)
DO td = 1,NQ(2)
DO rd = 1,NQ(1)

    tpd = Map_To_tpd(td,pd)
    LagP = Lagrange_Poly(CUR_X_LOCS(rd),DEGREE,Local_Locations)
    Tmp_U_Value = 0.0_idp

    DO lm = 1,LM_Length
    DO d = 0,DEGREE

        Current_Location = FP_Array_Map_TypeB(iU,iVB_S,re-1,d,lm)
        Tmp_U_Value = Tmp_U_Value                                       &
                    + FP_Coeff_Vector_B(Current_Location,iVB_S)         &
                    * LagP(d) * Ylm_Values( lm, tpd, te-1, pe-1 )


    END DO ! d Loop
    END DO ! lm Loop
    Here = rd + (td-1)*NQ(1) + (pd-1)*NQ(1)*NQ(2)
    Return_Shift(Here,re,te,pe,iU-2) = REAL(Tmp_U_Value, KIND = idp)

END DO ! rd Loop
END DO ! td Loop
END DO ! pd Loop
END DO ! re Loop
END DO ! te Loop
END DO ! pe Loop
END DO ! u Loop



END SUBROUTINE Poseidon_Return_Shift






!+203+######################################################################################!
!                                                                                           !
!       Poseidon_Return_ConFactor - Uses the results from the routine,                      !
!                                   Poseidon_XCFC_Run_Part2(), to calculate the value       !
!                                   of the Conformal Factor at the desired locations.       !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Return_ExtrinsicCurvature(  NE, NQ,                 &
                                                RQ_Input,               &
                                                TQ_Input,               &
                                                PQ_Input,               &
                                                Left_Limit,             &
                                                Right_Limit,            &
                                                Return_Kij              )

INTEGER,    DIMENSION(3),                                   INTENT(IN)  ::  NE, NQ
REAL(idp),  DIMENSION(NQ(1)),                               INTENT(IN)  ::  RQ_Input
REAL(idp),  DIMENSION(NQ(2)),                               INTENT(IN)  ::  TQ_Input
REAL(idp),  DIMENSION(NQ(3)),                               INTENT(IN)  ::  PQ_Input
REAL(idp),                                                  INTENT(IN)  ::  Left_Limit
REAL(idp),                                                  INTENT(IN)  ::  Right_Limit

REAL(idp),  DIMENSION(NQ(1)*NQ(2)*NQ(3),NE(1),NE(2),NE(3),1:6), INTENT(OUT) ::  Return_Kij


INTEGER                                                         ::  re, te, pe
INTEGER                                                         ::  rd, td, pd, tpd
INTEGER                                                         ::  d, lm
INTEGER                                                         ::  Here, There
INTEGER                                                         ::  i, j

REAL(idp)                                                       ::  Quad_Span
REAL(idp), DIMENSION(0:DEGREE)                                  ::  Local_Locations
REAL(idp), DIMENSION(0:DEGREE)                                  ::  LagP
REAL(idp), DIMENSION(1:NQ(1))                                   ::  CUR_R_LOCS
REAL(idp), DIMENSION(1:NQ(1))                                   ::  Cur_RX_Locs
REAL(idp), DIMENSION(1:NQ(2))                                   ::  CUR_T_LOCS
REAL(idp), DIMENSION(1:NQ(2))                                   ::  Cur_TX_Locs
REAL(idp), DIMENSION(1:NQ(3))                                   ::  CUR_P_LOCS
REAL(idp), DIMENSION(1:NQ(3))                                   ::  Cur_PX_Locs

REAL(idp)                                                       ::  DROT
REAL(idp)                                                       ::  DTOT
REAL(idp)                                                       ::  DPOT

REAL(idp), DIMENSION(3)                                         ::  gamma
REAL(idp), DIMENSION(3,3,3)                                     ::  Christoffel
COMPLEX(idp), DIMENSION(3)                                      ::  TMP_Val
COMPLEX(idp), DIMENSION(3,3)                                    ::  TMP_Drv

COMPLEX(idp), DIMENSION(4)                                      ::  Reusable_Vals
COMPLEX(idp), DIMENSION(6)                                      ::  Tmp_A


COMPLEX(idp)                                                    ::  Trace(2)

INTEGER                                                         ::  Current_Location
INTEGER                                                         ::  iVB
INTEGER, DIMENSION(3)                                           ::  iU


iVB   = iVB_X
iU(1) = iU_X1
iU(2) = iU_X2
iU(3) = iU_X3

Quad_Span = Right_Limit - Left_Limit

Local_Locations = Initialize_LGL_Quadrature_Locations(DEGREE)
CUR_RX_LOCS = 2.0_idp * ( RQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp
CUR_TX_LOCS = 2.0_idp * ( TQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp
CUR_PX_LOCS = 2.0_idp * ( PQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp

! nComp order : ij = 11, 12, 13, 22, 23, 33
Gamma(1) = 1.0_idp
Christoffel = 0.0_idp


DO pe = 1,NE(3)
DO te = 1,NE(2)
DO re = 1,NE(1)

DROT = 0.5_idp * (rlocs(re) - rlocs(re-1))
DTOT = 0.5_idp * (tlocs(te) - tlocs(te-1))
DPOT = 0.5_idp * (plocs(pe) - plocs(pe-1))

Cur_R_Locs(:) = DROT * (CUR_RX_LOCS(:)+1.0_idp) + rlocs(re)
Cur_T_Locs(:) = DTOT * (CUR_TX_LOCS(:)+1.0_idp) + tlocs(te)
Cur_P_Locs(:) = DPOT * (CUR_PX_LOCS(:)+1.0_idp) + plocs(pe)

DO pd = 1,NQ(3)
DO td = 1,NQ(2)
DO rd = 1,NQ(1)

    tpd = Map_To_tpd(td,pd)
    LagP = Lagrange_Poly(CUR_RX_LOCS(rd),DEGREE,Local_Locations)

    TMP_Val = 0.0_idp
    TMP_Drv = 0.0_idp
    DO i = 1,3
    DO d  = 0,DEGREE
        Here  = FP_Array_Map_TypeB(iU(i),iVB,re-1,d,1)
        There = FP_Array_Map_TypeB(iU(i),iVB,re-1,d,LM_Length)

        TMP_Val(i) = TMP_Val(i)                                  &
                + SUM( FP_Coeff_Vector_B( Here:There, iVB )      &
                        * Ylm_Values( :, tpd, te-1, pe-1 )   )       &
                * Lagrange_Poly_Table( d, rd, 0 )


        TMP_Drv(1,i) = TMP_Drv(1,i)                              &
                   + SUM( FP_Coeff_Vector_B( Here:There, iVB )   &
                         * Ylm_Values( :, tpd, te-1, pe-1 )     )    &
                   * Lagrange_Poly_Table( d, rd, 1 )             &
                   / DROT


        TMP_Drv(2,i) = TMP_Drv(2,i)                              &
                   + SUM( FP_Coeff_Vector_B( Here:There, iVB )   &
                         * Ylm_dt_Values( :, tpd, te-1, pe-1)   )    &
                   * Lagrange_Poly_Table( d, rd, 0)

        TMP_Drv(3,i) = TMP_Drv(3,i)                              &
                   + SUM( FP_Coeff_Vector_B( Here:There, iVB )   &
                         * Ylm_dp_Values( :, tpd, te-1, pe-1)   )    &
                   * Lagrange_Poly_Table( d, rd, 0)

    END DO  ! d
    END DO  ! i


    Gamma(2) = 1.0_idp/(Cur_R_Locs(rd)*Cur_R_Locs(rd))
    Gamma(3) = Gamma(2) * 1.0_idp/( DSIN(Cur_T_Locs(td))*DSIN(Cur_T_Locs(td)) )

!    PRINT*,"r",re,rd,Cur_R_Locs(rd)
!    PRINT*,"Sin, Theta",DSIN(Cur_T_Locs(td)),Cur_T_locs(Td)
!    PRINT*,"Gamma",1.0_idp/Gamma

    Christoffel(1,2,2) = -Cur_R_Locs(rd)
    Christoffel(1,3,3) = -Cur_R_Locs(rd)*DSIN(Cur_T_Locs(td))*DSIN(Cur_T_Locs(td))
    
    Christoffel(2,1,2) = 1.0_idp/Cur_R_Locs(rd)
    Christoffel(2,2,1) = 1.0_idp/Cur_R_Locs(rd)
    Christoffel(2,3,3) = -DSIN(Cur_T_Locs(td))*DCOS(Cur_P_Locs(pd))

    Christoffel(3,3,1) = 1.0_idp/Cur_R_Locs(rd)
    Christoffel(3,1,3) = 1.0_idp/Cur_R_Locs(rd)
    Christoffel(3,3,2) = 1.0_idp/DTAN(CUR_T_LOCS(td))
    Christoffel(3,2,3) = 1.0_idp/DTAN(CUR_T_LOCS(td))


    Reusable_Vals(1) = 2.0_idp/3.0_idp*(Tmp_Drv(1,1)+Tmp_Drv(2,2)+Tmp_Drv(3,3))
    Reusable_Vals(2) = 2.0_idp/3.0_idp*(Christoffel(1,1,1)+Christoffel(2,2,1)+Christoffel(3,3,1))
    Reusable_Vals(3) = 2.0_idp/3.0_idp*(Christoffel(1,1,2)+Christoffel(2,2,2)+Christoffel(3,3,2))
    Reusable_Vals(4) = 2.0_idp/3.0_idp*(Christoffel(1,1,3)+Christoffel(2,2,3)+Christoffel(3,3,3))

    
    Here = rd + (td-1)*NQ(1) + (pd-1)*NQ(1)*NQ(2)


    ! Ahat^11
    Tmp_A(1) = 2.0_idp * Gamma(1) * Tmp_Drv(1,1) + Reusable_Vals(1)                  &
             +(2.0_idp * Gamma(1)*Christoffel(1,1,1) - Reusable_Vals(2) )*Tmp_Val(1) &
             +(2.0_idp * Gamma(1)*Christoffel(1,1,2) - Reusable_Vals(3) )*Tmp_Val(2) &
             +(2.0_idp * Gamma(1)*Christoffel(1,1,3) - Reusable_Vals(4) )*Tmp_Val(3)

    Return_Kij(Here,re,te,pe,1) = Tmp_A(1)/(Gamma(1)*Gamma(1))


    ! Ahat^12
    i=1
    j=2
    Tmp_A(2) = Gamma(i)*Tmp_Drv(i,j)             &
             + Gamma(j)*Tmp_Drv(j,i)             &
             +( Gamma(i)*Christoffel(j,i,1)      &
              + Gamma(j)*Christoffel(i,j,1)      &
              )*Tmp_Val(1)                       &
             +( Gamma(i)*Christoffel(j,i,2)      &
              + Gamma(j)*Christoffel(i,j,2)      &
              )*Tmp_Val(2)                       &
             +( Gamma(i)*Christoffel(j,i,3)      &
              + Gamma(j)*Christoffel(i,j,3)      &
              )*Tmp_Val(3)
     Return_Kij(Here,re,te,pe,2) = Tmp_A(2)/(Gamma(1)*Gamma(2))

    ! Ahat^13
    i=1
    j=3
    Tmp_A(3) = Gamma(i)*Tmp_Drv(i,j)             &
             + Gamma(j)*Tmp_Drv(j,i)             &
             +( Gamma(i)*Christoffel(j,i,1)      &
              + Gamma(j)*Christoffel(i,j,1)      &
              )*Tmp_Val(1)                       &
             +( Gamma(i)*Christoffel(j,i,2)      &
              + Gamma(j)*Christoffel(i,j,2)      &
              )*Tmp_Val(2)                       &
             +( Gamma(i)*Christoffel(j,i,3)      &
              + Gamma(j)*Christoffel(i,j,3)      &
              )*Tmp_Val(3)
     Return_Kij(Here,re,te,pe,3) = Tmp_A(3)/(Gamma(1)*Gamma(3))


    ! Ahat^22
    Tmp_A(4) = 2.0_idp * Gamma(2) * Tmp_Drv(2,2) + Reusable_Vals(1)                      &
             +(2.0_idp * Gamma(2)*Christoffel(2,2,1) - Reusable_Vals(2) )*Tmp_Val(1)     &
             +(2.0_idp * Gamma(2)*Christoffel(2,2,2) - Reusable_Vals(3) )*Tmp_Val(2)     &
             +(2.0_idp * Gamma(2)*Christoffel(2,2,3) - Reusable_Vals(4) )*Tmp_Val(3)
     Return_Kij(Here,re,te,pe,4) = Tmp_A(4)/(Gamma(2)*Gamma(2))


    ! Ahat^23
    i=2
    j=3
    Tmp_A(5) = Gamma(i)*Tmp_Drv(i,j)             &
             + Gamma(j)*Tmp_Drv(j,i)             &
             +( Gamma(i)*Christoffel(j,i,1)      &
              + Gamma(j)*Christoffel(i,j,1)      &
              )*Tmp_Val(1)                       &
             +( Gamma(i)*Christoffel(j,i,2)      &
              + Gamma(j)*Christoffel(i,j,2)      &
              )*Tmp_Val(2)                       &
             +( Gamma(i)*Christoffel(j,i,3)      &
              + Gamma(j)*Christoffel(i,j,3)      &
              )*Tmp_Val(3)
     Return_Kij(Here,re,te,pe,5) = Tmp_A(5)/(Gamma(2)*Gamma(3))


    ! Ahat^33
    Tmp_A(6) = 2.0_idp * Gamma(3) * Tmp_Drv(3,3) + Reusable_Vals(1)                      &
             +(2.0_idp * Gamma(3)*Christoffel(3,3,1) - Reusable_Vals(2) )*Tmp_Val(1)     &
             +(2.0_idp * Gamma(3)*Christoffel(3,3,2) - Reusable_Vals(3) )*Tmp_Val(2)     &
             +(2.0_idp * Gamma(3)*Christoffel(3,3,3) - Reusable_Vals(4) )*Tmp_Val(3)
     Return_Kij(Here,re,te,pe,6) = Tmp_A(6)/(Gamma(3)*Gamma(3))


!    Return_Kij(Here,re,te,pe,1) =
!    PRINT*,TMP_A(1),TMP_A(4),TMP_A(6)

!    PRINT*,Return_Kij(Here,re,te,pe,1),          &
!            Return_Kij(Here,re,te,pe,4),          &
!            Return_Kij(Here,re,te,pe,6)


    Trace(1) = REAL(Tmp_A(1)/Gamma(1)+Tmp_A(4)/Gamma(2)+Tmp_A(6)/Gamma(3),KIND = idp)

    Trace(2) = Trace(1)                                     &
             /REAL(3.0_idp*(  Reusable_Vals(1)                  &
                        + Reusable_Vals(2)*Tmp_Val(1)       &
                        + Reusable_Vals(3)*Tmp_Val(2)       &
                        + Reusable_Vals(4)*Tmp_Val(3)  ), Kind = idp    )

    PRINT*,Trace(1), Trace(2)

!    PRINT*,Trace(1),                                     &
!            REAL(3.0_idp*(  Reusable_Vals(1)                  &
!               + Reusable_Vals(2)*Tmp_Val(1)       &
!               + Reusable_Vals(3)*Tmp_Val(2)       &
!               + Reusable_Vals(4)*Tmp_Val(3)  ), Kind = idp    ), Trace(2)




END DO ! rd Loop
END DO ! td Loop
END DO ! pd Loop
END DO ! re Loop
END DO ! te Loop
END DO ! pe Loop




END SUBROUTINE Poseidon_Return_ExtrinsicCurvature

END MODULE Poseidon_XCFC_Interface_Module
