   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Return_Functions_Native                                               !##!
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

USE Poseidon_Parameters, &
            ONLY :  Degree,                     &
                    L_Limit

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
            ONLY :  Lagrange_Poly_Table,        &
                    Slm_Elem_Values,            &
                    Nlm_Values

USE Variables_Derived, &
            ONLY :  LM_Length,                  &
                    LM_Short_Length

USE Variables_Vectors, &
            ONLY :  dVA_Coeff_Vector,      &
                    dVB_Coeff_Vector

USE Variables_Mesh, &
            ONLY :  rlocs,              &
                    tlocs,              &
                    plocs

USE Variables_AMReX_Source, &
            ONLY :  iLeaf,                &
                    iTrunk

USE Variables_FEM_Module, &
            ONLY :  FEM_Node_xlocs

USE Maps_Fixed_Point, &
            ONLY :  FP_Array_Map_TypeB

USE Maps_Quadrature, &
            ONLY :  Map_To_tpd,             &
                    Quad_Map

USE Maps_Domain, &
            ONLY :  Map_To_FEM_Node,        &
                    FEM_Elem_Map

USE Functions_Quadrature, &
            ONLY :  Initialize_LGL_Quadrature_Locations

USE Functions_Math, &
            ONLY :  Lagrange_Poly
            
USE Initialization_Tables_Slm, &
            ONLY :  Initialize_Am_Table,            &
                    Initialize_Plm_Table,           &
                    Initialize_Slm_Table_on_Elem



IMPLICIT NONE


CONTAINS



!+101+########################################################!
!                                                               !
!       Poseidon_Return_AMReX_Type_A                            !
!                                                               !
!#############################################################!
SUBROUTINE Poseidon_Return_Native_Type_A(   iU,                     &
                                            NE, NQ,                 &
                                            RQ_Input,               &
                                            TQ_Input,               &
                                            PQ_Input,               &
                                            Left_Limit,             &
                                            Right_Limit,            &
                                            Return_Variable        )


INTEGER,                                                    INTENT(IN)  ::  iU
INTEGER,    DIMENSION(3),                                   INTENT(IN)  ::  NE, NQ
REAL(idp),  DIMENSION(NQ(1)),                               INTENT(IN)  ::  RQ_Input
REAL(idp),  DIMENSION(NQ(2)),                               INTENT(IN)  ::  TQ_Input
REAL(idp),  DIMENSION(NQ(3)),                               INTENT(IN)  ::  PQ_Input
REAL(idp),                                                  INTENT(IN)  ::  Left_Limit
REAL(idp),                                                  INTENT(IN)  ::  Right_Limit

REAL(idp),  DIMENSION(NQ(1)*NQ(2)*NQ(3),NE(1),NE(2),NE(3)), INTENT(OUT) ::  Return_Variable



INTEGER                                                         ::  re, te, pe
INTEGER                                                         ::  rd, td, pd, tpd
INTEGER                                                         ::  d, lm, Here

REAL(KIND = idp)                                                ::  Quad_Span
REAL(KIND = idp), DIMENSION(0:DEGREE)                           ::  LagP
REAL(KIND = idp), DIMENSION(1:NQ(1))                            ::  CUR_X_LOCS
REAL(KIND = idp)                                                ::  Tmp_U_Value
INTEGER                                                         ::  Current_Location


REAL(idp),  DIMENSION(1:NQ(2),1:LM_Short_Length,0:NE(2)-1)      ::  Plm_Table
REAL(idp),  DIMENSION(1:NQ(3),1:LM_Length,0:NE(3)-1)            ::  Am_Table
REAL(idp),  DIMENSION(1:LM_Length, 1:NQ(2)*NQ(3) )              ::  Slm_Elem_Table

Quad_Span = Right_Limit - Left_Limit

CUR_X_LOCS = 2.0_idp * ( RQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp

! Initialize Am Table
CALL Initialize_Am_Table(   NQ(3),                      &
                            PQ_Input,                   &
                            L_Limit,                    &
                            NE(3),                      &
                            [0, NE(3)-1],               &
                            plocs,                      &
                            Am_Table                    )

! Initialize Plm Table
CALL Initialize_Plm_Table(  NQ(2),                      &
                            TQ_Input,                   &
                            L_Limit,                    &
                            LM_Short_Length,            &
                            NE(2),                      &
                            [0, nE(2)-1],               &
                            tlocs,                      &
                            Plm_Table                   )
                            
                            
DO pe = 1,NE(3)
DO te = 1,NE(2)
DO re = 1,NE(1)

CALL Initialize_Slm_Table_on_Elem(  te, pe,             &
                                    NQ(2), NQ(3),       &
                                    NE,                 &
                                    [1,1,1],            &
                                    Plm_Table,          &
                                    Am_Table,           &
                                    Slm_Elem_Table      )

DO pd = 1,NQ(3)
DO td = 1,NQ(2)
DO rd = 1,NQ(1)

    tpd = Map_To_tpd(td,pd)
    LagP = Lagrange_Poly(CUR_X_LOCS(rd),DEGREE,FEM_Node_xlocs)
    Tmp_U_Value = 0.0_idp

    DO lm = 1,LM_Length
    DO d = 0,DEGREE

        Current_Location = Map_To_FEM_Node(re-1,d)
        Tmp_U_Value = Tmp_U_Value + dVA_Coeff_Vector(Current_Location,lm,iU)  &
                               * LagP(d) * Slm_Elem_Values( lm, tpd )

    END DO ! d Loop
    END DO ! lm Loop

    
    Here = Quad_Map(rd, td, pd, NQ(1), NQ(2),NQ(3))
    Return_Variable(Here,re,te,pe) = REAL(Tmp_U_Value, KIND = idp)


END DO ! rd Loop
END DO ! td Loop
END DO ! pd Loop

END DO ! re Loop
END DO ! te Loop
END DO ! pe Loop

END SUBROUTINE Poseidon_Return_Native_Type_A








 !+201+########################################################!
!                                                               !
!       Poseidon_Return_AMReX_Type_B                            !
!                                                               !
 !#############################################################!
SUBROUTINE Poseidon_Return_Native_Type_B(  iU, iVB,             &
                                           NE, NQ,              &
                                           RQ_Input,            &
                                           TQ_Input,            &
                                           PQ_Input,            &
                                           Left_Limit,          &
                                           Right_Limit,         &
                                           Return_Variable         )

INTEGER,                                                    INTENT(IN)  ::  iU, iVB
INTEGER,    DIMENSION(3),                                   INTENT(IN)  ::  NE, NQ
REAL(idp),  DIMENSION(NQ(1)),                               INTENT(IN)  ::  RQ_Input
REAL(idp),  DIMENSION(NQ(2)),                               INTENT(IN)  ::  TQ_Input
REAL(idp),  DIMENSION(NQ(3)),                               INTENT(IN)  ::  PQ_Input
REAL(idp),                                                  INTENT(IN)  ::  Left_Limit
REAL(idp),                                                  INTENT(IN)  ::  Right_Limit

REAL(idp),  DIMENSION(NQ(1)*NQ(2)*NQ(3),NE(1),NE(2),NE(3)), INTENT(OUT) ::  Return_Variable


INTEGER                                                         ::  re, te, pe
INTEGER                                                         ::  rd, td, pd, tpd
INTEGER                                                         ::  d, lm, Here

REAL(KIND = idp)                                                ::  Quad_Span
REAL(KIND = idp), DIMENSION(0:DEGREE)                           ::  LagP
REAL(KIND = idp), DIMENSION(1:NQ(1))                            ::  CUR_X_LOCS
REAL(KIND = idp)                                                ::  TMP_U_Value
INTEGER                                                         ::  Current_Location


REAL(idp),  DIMENSION(1:NQ(2),1:LM_Short_Length,0:NE(2)-1)      ::  Plm_Table
REAL(idp),  DIMENSION(1:NQ(3),1:LM_Length,0:NE(3)-1)            ::  Am_Table
REAL(idp),  DIMENSION(1:LM_Length, 1:NQ(2)*NQ(3) )              ::  Slm_Elem_Table

Quad_Span = Right_Limit - Left_Limit

CUR_X_LOCS = 2.0_idp * ( RQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp

! Initialize Am Table
CALL Initialize_Am_Table(   NQ(3),                      &
                            PQ_Input,                   &
                            L_Limit,                    &
                            NE(3),                      &
                            [0, NE(3)-1],               &
                            plocs,                      &
                            Am_Table                    )

! Initialize Plm Table
CALL Initialize_Plm_Table(  NQ(2),                      &
                            TQ_Input,                   &
                            L_Limit,                    &
                            LM_Short_Length,            &
                            NE(2),                      &
                            [0, nE(2)-1],               &
                            tlocs,                      &
                            Plm_Table                   )
DO pe = 1,NE(3)
DO te = 1,NE(2)
DO re = 1,NE(1)

CALL Initialize_Slm_Table_on_Elem(  te, pe,             &
                                    NQ(2), NQ(3),       &
                                    NE,                 &
                                    [1,1,1],            &
                                    Plm_Table,          &
                                    Am_Table,           &
                                    Slm_Elem_Table      )
DO pd = 1,NQ(3)
DO td = 1,NQ(2)
DO rd = 1,NQ(1)

    tpd = Map_To_tpd(td,pd)
    LagP = Lagrange_Poly(CUR_X_LOCS(rd),DEGREE,FEM_Node_xlocs)
    Tmp_U_Value = 0.0_idp

    DO lm = 1,LM_Length
    DO d = 0,DEGREE

        Current_Location = FP_Array_Map_TypeB(iU,iVB,re-1,d,lm)
        Tmp_U_Value = Tmp_U_Value                                       &
                    + dVB_Coeff_Vector(Current_Location,iVB_S)         &
                    * LagP(d) * Slm_Elem_Values( lm, tpd )


    END DO ! d Loop
    END DO ! lm Loop
    
    Here = Quad_Map(rd, td, pd, NQ(1), NQ(2),NQ(3))
    Return_Variable(Here,re,te,pe) = REAL(Tmp_U_Value, KIND = idp)

END DO ! rd Loop
END DO ! td Loop
END DO ! pd Loop
END DO ! re Loop
END DO ! te Loop
END DO ! pe Loop






END SUBROUTINE Poseidon_Return_Native_Type_B







END MODULE Return_Functions_Native
