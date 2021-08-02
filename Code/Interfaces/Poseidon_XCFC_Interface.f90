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


USE Variables_Tables, &
        ONLY :  Ylm_Values,                 &
                Lagrange_Poly_Table

USE Variables_Derived, &
        ONLY :  LM_LENGTH

USE Variables_FP, &
        ONLY :  FP_Coeff_Vector


USE FP_Functions_Mapping, &
        ONLY :  FP_FEM_Node_Map,    &
                FP_tpd_Map

USE Functions_Quadrature, &
        ONLY :  Initialize_LGL_Quadrature_Locations

USE Functions_Math, &
        ONLY :  Lagrange_Poly


USE XCFC_Main_Module, &
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
INTEGER                                                         ::  u, d, lm, Here

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

        tpd = FP_tpd_Map(td,pd)
        LagP = Lagrange_Poly(CUR_X_LOCS(rd),DEGREE,Local_Locations)
        Tmp_U_Value = 0.0_idp

        DO lm = 0,LM_Length
        DO d = 0,DEGREE

            Current_Location = FP_FEM_Node_Map(re,d)
            Tmp_U_Value = Tmp_U_Value + FP_Coeff_Vector(Current_Location,lm,1)  &
                                      * LagP(d) * Ylm_Values( lm, tpd, te, pe )
        END DO ! d Loop
        END DO ! lm Loop

        Here = rd + (te-1)*NQ(1) + (pe-1)*NQ(1)*NQ(2)
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
INTEGER                                                         ::  u, d, lm, Here

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

        tpd = FP_tpd_Map(td,pd)
        LagP = Lagrange_Poly(CUR_X_LOCS(rd),DEGREE,Local_Locations)
        Tmp_U_Value = 0.0_idp

        DO lm = 0,LM_Length
        DO d = 0,DEGREE

            Current_Location = FP_FEM_Node_Map(re,d)
            Tmp_U_Value(1) = Tmp_U_Value(1) + FP_Coeff_Vector(Current_Location,lm,1)  &
                                            * LagP(d) * Ylm_Values( lm, tpd, te, pe )

            Tmp_U_Value(2) = Tmp_U_Value(2) + FP_Coeff_Vector(Current_Location,lm,2)  &
                                            * LagP(d) * Ylm_Values( lm, tpd, te, pe )

        END DO ! d Loop
        END DO ! lm Loop

        Here = rd + (te-1)*NQ(1) + (pe-1)*NQ(1)*NQ(2)
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
INTEGER                                                         ::  u, d, lm, Here

REAL(KIND = idp)                                                ::  Quad_Span
REAL(KIND = idp), DIMENSION(0:DEGREE)                           ::  Local_Locations
REAL(KIND = idp), DIMENSION(0:DEGREE)                           ::  LagP
REAL(KIND = idp), DIMENSION(1:NQ(1))                            ::  CUR_X_LOCS
COMPLEX(KIND = idp), DIMENSION(1:3)                             ::  TMP_U_Value
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

        tpd = FP_tpd_Map(td,pd)
        LagP = Lagrange_Poly(CUR_X_LOCS(rd),DEGREE,Local_Locations)
        Tmp_U_Value = 0.0_idp

        DO lm = 0,LM_Length
        DO d = 0,DEGREE

            Current_Location = FP_FEM_Node_Map(re,d)
            Tmp_U_Value(1) = Tmp_U_Value(1) + FP_Coeff_Vector(Current_Location,lm,3)  &
                                            * LagP(d) * Ylm_Values( lm, tpd, te, pe )

            Tmp_U_Value(2) = Tmp_U_Value(2) + FP_Coeff_Vector(Current_Location,lm,4)  &
                                            * LagP(d) * Ylm_Values( lm, tpd, te, pe )

            Tmp_U_Value(3) = Tmp_U_Value(2) + FP_Coeff_Vector(Current_Location,lm,5)  &
                                            * LagP(d) * Ylm_Values( lm, tpd, te, pe )

        END DO ! d Loop
        END DO ! lm Loop

        Here = rd + (te-1)*NQ(1) + (pe-1)*NQ(1)*NQ(2)
        Return_Shift(Here,re,te,pe,1) = REAL(Tmp_U_Value(1), KIND = idp)
        Return_Shift(Here,re,te,pe,2) = REAL(Tmp_U_Value(2), KIND = idp)
        Return_Shift(Here,re,te,pe,3) = REAL(Tmp_U_Value(3), KIND = idp)

    END DO ! rd Loop
    END DO ! td Loop
    END DO ! pd Loop

END DO ! re Loop
END DO ! te Loop
END DO ! pe Loop




END SUBROUTINE Poseidon_Return_Shift





END MODULE Poseidon_XCFC_Interface_Module
