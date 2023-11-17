  !##########################################################################!
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Poseidon_Return_Routines_Kij                                           !##!
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
            ONLY :  idp

USE Poseidon_Parameters, &
            ONLY :  Degree,                      &
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
            ONLY :  Slm_Elem_Values,        &
                    Slm_Elem_dt_Values,     &
                    Slm_Elem_dp_Values,              &
                    Lagrange_Poly_Table,        &
                    Level_DX

USE Variables_Mesh, &
           ONLY :  Num_R_Elements,         &
                   Num_T_Elements,         &
                   Num_P_Elements

USE Variables_Derived, &
            ONLY :  LM_Length,              &
                    LM_Short_Length

USE Variables_Vectors, &
            ONLY :  dVA_Coeff_Vector,       &
                    dVB_Coeff_Vector

USE Variables_Mesh, &
           ONLY :  rlocs,              &
                   tlocs,              &
                   plocs,              &
                   iNE_Base

USE Parameters_AMReX, &
           ONLY :  iLeaf,                &
                   iTrunk

USE Variables_FEM_Module, &
            ONLY :  FEM_Node_xlocs

USE Maps_Fixed_Point, &
            ONLY :  FP_Array_Map_TypeB

USE Maps_Quadrature, &
            ONLY :  Map_To_tpd

USE Maps_Domain, &
            ONLY :  Map_To_FEM_Node,        &
                   FEM_Elem_Map

USE Functions_Quadrature, &
            ONLY :  Initialize_LGL_Quadrature_Locations

USE Functions_Math, &
            ONLY :  Lagrange_Poly,                   &
                    Lagrange_Poly_Deriv

USE Initialization_Tables_Slm, &
            ONLY :  Initialize_Am_Tables,            &
                    Initialize_Plm_Tables,           &
                    Initialize_Slm_Tables_on_Elem


USE Variables_Interface, &
            ONLY :  Caller_NQ,                      &
                    Caller_Quad_DOF,                      &
                    Caller_xL,                      &
                    Caller_RQ_xlocs,                &
                    Caller_TQ_xlocs,                &
                    Caller_PQ_xlocs


#ifdef POSEIDON_AMREX_FLAG
use amrex_base_module

USE amrex_box_module,   &
           ONLY :  amrex_box

USE amrex_boxarray_module, &
           ONLY :  amrex_boxarray

use amrex_fort_module, &
           ONLY :  amrex_spacedim

USE amrex_multifab_module,  &
           ONLY :  amrex_multifab,         &
                   amrex_multifab_build,   &
                   amrex_imultifab_build,  &
                   amrex_imultifab_destroy

USE Variables_AMReX_Core, &
            ONLY :  AMReX_Num_Levels,        &
                    AMReX_Max_Grid_Size

USE Poseidon_AMReX_MakeFineMask_Module, &
            ONLY :  AMReX_MakeFineMask

#endif

USE Return_Functions_Native, &
            ONLY :  Poseidon_Return_Native_Type_A,  &
                    Poseidon_Return_Native_Type_B

USE Return_Functions_AMReX, &
            ONLY :  Poseidon_Return_AMReX_Type_A,  &
                    Poseidon_Return_AMReX_Type_B



IMPLICIT NONE



INTEGER                 ::  iU_K11 = 6
INTEGER                 ::  iU_K12 = 7
INTEGER                 ::  iU_K13 = 8
INTEGER                 ::  iU_K22 = 9
INTEGER                 ::  iU_K23 = 10
INTEGER                 ::  iU_K33 = 11



CONTAINS
!+101+######################################################################################!
!                                                                                           !
!       Poseidon_Return_ConFactor - Uses the results from the routine,                      !
!                                   Poseidon_XCFC_Run_Part1(), to calculate the value       !
!                                   of the Conformal Factor at the desired locations.       !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Return_Kij_Native(  NE, NQ,                 &
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

REAL(idp),  DIMENSION(NQ(1)*NQ(2)*NQ(3),                                &
                      0:NE(1)-1,0:NE(2)-1,0:NE(3)-1,                    &
                      1:6),                                 INTENT(OUT) ::  Return_Kij


INTEGER                                                         ::  re, te, pe
INTEGER                                                         ::  rd, td, pd, tpd
INTEGER                                                         ::  d
INTEGER                                                         ::  Here, There
INTEGER                                                         ::  i, j, l, m, lm

REAL(idp)                                                       ::  Quad_Span
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
REAL(idp), DIMENSION(3)                                      ::  TMP_Val
REAL(idp), DIMENSION(3,3)                                    ::  TMP_Drv

REAL(idp), DIMENSION(4)                                      ::  Reusable_Vals
REAL(idp), DIMENSION(6)                                      ::  Tmp_A
REAL(idp)                                                    ::  Tmp_Val_Psi

REAL(idp)                                                    ::  Trace(2)

INTEGER                                                         ::  iVB
INTEGER, DIMENSION(3)                                           ::  iU



INTEGER                                                         ::  Long_LM
INTEGER                                                         ::  Short_LM

REAL(idp), DIMENSION(:,:,:), ALLOCATABLE                        ::  Caller_LPT


REAL(idp), DIMENSION(1:NQ(3),-L_Limit:L_Limit,0:NE(3)-1)        ::  Am_Table
REAL(idp), DIMENSION(1:NQ(3),-L_Limit:L_Limit,0:NE(3)-1)        ::  Am_dp_Table

REAL(idp), DIMENSION(1:NQ(2),0:L_Limit,0:NE(2)-1)               ::  Plm_Table
REAL(idp), DIMENSION(1:NQ(2),0:L_Limit,0:NE(2)-1)               ::  Plm_dt_Table

REAL(idp), DIMENSION( 1:LM_Length, 1:NQ(2)*NQ(3))               ::  Slm_Elem_Table
REAL(idp), DIMENSION( 1:LM_Length, 1:NQ(2)*NQ(3))               ::  Slm_Elem_dt_Table
REAL(idp), DIMENSION( 1:LM_Length, 1:NQ(2)*NQ(3))               ::  Slm_Elem_dp_Table

ALLOCATE( Caller_LPT(0:1,0:DEGREE,1:NQ(1)))

iVB   = iVB_X
iU(1) = iU_X1
iU(2) = iU_X2
iU(3) = iU_X3

Quad_Span = Right_Limit - Left_Limit

CUR_RX_LOCS = 2.0_idp * ( RQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp
CUR_TX_LOCS = 2.0_idp * ( TQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp
CUR_PX_LOCS = 2.0_idp * ( PQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp

! nComp order : ij = 11, 12, 13, 22, 23, 33
Gamma(1) = 1.0_idp
Christoffel = 0.0_idp


CALL Initialize_Am_Tables(  NQ(3),                  &
                            CUR_PX_LOCS,            &
                            L_Limit,                &
                            NE(3),                  &
                            [0, NE(3)],             &
                            plocs,                  &
                            Am_Table,               &
                            Am_dp_Table             )

CALL Initialize_Plm_Tables( NQ(2),                  &
                            CUR_TX_LOCS,            &
                            L_Limit,                &
                            LM_Short_Length,        &
                            NE(2),                  &
                            [0, NE(2)],             &
                            tlocs,                  &
                            Plm_Table,              &
                            Plm_dt_Table            )

DO pe = 0,NE(3)-1
DO te = 0,NE(2)-1

CALL Initialize_Slm_Tables_on_Elem( te, pe,             &
                                   NQ(2), NQ(3),        &
                                   NE,                  &
                                   [0,0,0],             &
                                   Plm_Table,            &
                                   Plm_dt_Table,         &
                                   Am_Table,             &
                                   Am_dp_Table,          &
                                   Slm_Elem_Table,       &
                                   Slm_Elem_dt_Table,    &
                                   Slm_Elem_dp_Table     )

DO re = 0,NE(1)-1

DROT = 0.5_idp * (rlocs(re+1) - rlocs(re))
DTOT = 0.5_idp * (tlocs(te+1) - tlocs(te))
DPOT = 0.5_idp * (plocs(pe+1) - plocs(pe))

Cur_R_Locs(:) = DROT * (CUR_RX_LOCS(:)+1.0_idp) + rlocs(re)
Cur_T_Locs(:) = DTOT * (CUR_TX_LOCS(:)+1.0_idp) + tlocs(te)
Cur_P_Locs(:) = DPOT * (CUR_PX_LOCS(:)+1.0_idp) + plocs(pe)


DO rd = 1,NQ(1)
    Caller_LPT(0,:,rd)=Lagrange_Poly(CUR_RX_LOCS(rd),DEGREE,FEM_Node_xlocs)
    Caller_LPT(1,:,rd)=Lagrange_Poly_Deriv(CUR_RX_LOCS(rd),DEGREE,FEM_Node_xlocs)
END DO



DO pd = 1,NQ(3)
DO td = 1,NQ(2)
DO rd = 1,NQ(1)

    tpd = Map_To_tpd(td,pd)

    TMP_Val = 0.0_idp
    TMP_Drv = 0.0_idp
    
    TMP_Val_Psi = 0.0_idp
    DO lm = 1,LM_Length
    DO d  = 0,DEGREE

        Here  = Map_To_FEM_Node(RE,d)

        TMP_Val_Psi = TMP_Val_Psi                           &
                    + dVA_Coeff_Vector(Here,lm,iU_CF)       &
                    * Caller_LPT(0,d,rd)                    &
                    * Slm_Elem_Table( lm, tpd )
        
    END DO
    END DO
    
    
    DO i = 1,3
    DO d  = 0,DEGREE
        Here  = FP_Array_Map_TypeB(iU(i),iVB,re,d,1)
        There = FP_Array_Map_TypeB(iU(i),iVB,re,d,LM_Length)

        TMP_Val(i) = TMP_Val(i)                                  &
                + SUM( dVB_Coeff_Vector( Here:There, iVB )      &
                        * Slm_Elem_Values(:,tpd)   )       &
                * Caller_LPT(0,d,rd)


        TMP_Drv(1,i) = TMP_Drv(1,i)                              &
                   + SUM( dVB_Coeff_Vector( Here:There, iVB )   &
                         * Slm_Elem_Values(:,tpd)     )    &
                   * Caller_LPT(1,d,rd)             &
                   / DROT


        TMP_Drv(2,i) = TMP_Drv(2,i)                              &
                   + SUM( dVB_Coeff_Vector( Here:There, iVB )   &
                         * Slm_Elem_dt_Values(:,tpd)   )    &
                   * Caller_LPT(0,d,rd)

        TMP_Drv(3,i) = TMP_Drv(3,i)                              &
                   + SUM( dVB_Coeff_Vector( Here:There, iVB )   &
                         * Slm_Elem_dp_Values(:,tpd)   )    &
                   * Caller_LPT(0,d,rd)

    END DO  ! d
    END DO  ! i


    Gamma(2) = 1.0_idp/(Cur_R_Locs(rd)*Cur_R_Locs(rd))
    Gamma(3) = Gamma(2) * 1.0_idp/( DSIN(Cur_T_Locs(td))*DSIN(Cur_T_Locs(td)) )

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
    Tmp_A(1) = Gamma(1)                                                         &
             * ( 2.0_idp * Tmp_Drv(1,1) - Reusable_Vals(1)                      &
               +(2.0_idp * Christoffel(1,1,1) - Reusable_Vals(2) )*Tmp_Val(1)   &
               +(2.0_idp * Christoffel(1,1,2) - Reusable_Vals(3) )*Tmp_Val(2)   &
               +(2.0_idp * Christoffel(1,1,3) - Reusable_Vals(4) )*Tmp_Val(3)   )

    Return_Kij(Here,re,te,pe,1) = REAL(Tmp_A(1)/(Gamma(1)*Gamma(1)*Tmp_Val_Psi*Tmp_Val_Psi), KIND = idp)
    

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
     Return_Kij(Here,re,te,pe,2) = REAL(Tmp_A(2)/(Gamma(1)*Gamma(2)*Tmp_Val_Psi*Tmp_Val_Psi), KIND = idp)

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
     Return_Kij(Here,re,te,pe,3) = REAL(Tmp_A(3)/(Gamma(1)*Gamma(3)*Tmp_Val_Psi*Tmp_Val_Psi), KIND = idp)


    ! Ahat^22
    Tmp_A(4) = Gamma(2)                                                             &
             * ( 2.0_idp * Tmp_Drv(2,2) - Reusable_Vals(1)                          &
               +(2.0_idp * Christoffel(2,2,1) - Reusable_Vals(2) ) * Tmp_Val(1)     &
               +(2.0_idp * Christoffel(2,2,2) - Reusable_Vals(3) ) * Tmp_Val(2)     &
               +(2.0_idp * Christoffel(2,2,3) - Reusable_Vals(4) ) * Tmp_Val(3)     )
     Return_Kij(Here,re,te,pe,4) = REAL(Tmp_A(4)/(Gamma(2)*Gamma(2)*Tmp_Val_Psi*Tmp_Val_Psi), KIND = idp)


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
     Return_Kij(Here,re,te,pe,5) = REAL(Tmp_A(5)/(Gamma(2)*Gamma(3)*Tmp_Val_Psi*Tmp_Val_Psi), KIND = idp)


    ! Ahat^33
    Tmp_A(6) = Gamma(3)                                                             &
             * ( 2.0_idp * Tmp_Drv(3,3) - Reusable_Vals(1)                          &
               +(2.0_idp * Christoffel(3,3,1) - Reusable_Vals(2) ) * Tmp_Val(1)     &
               +(2.0_idp * Christoffel(3,3,2) - Reusable_Vals(3) ) * Tmp_Val(2)     &
               +(2.0_idp * Christoffel(3,3,3) - Reusable_Vals(4) ) * Tmp_Val(3)     )
     Return_Kij(Here,re,te,pe,6) = REAL(Tmp_A(6)/(Gamma(3)*Gamma(3)*Tmp_Val_Psi*Tmp_Val_Psi), KIND = idp)



    Trace(1) = REAL(Tmp_A(1)/Gamma(1)+Tmp_A(4)/Gamma(2)+Tmp_A(6)/Gamma(3),KIND = idp)

    Trace(2) = Trace(1)                                     &
             /REAL(3.0_idp*(  Reusable_Vals(1)                  &
                        + Reusable_Vals(2)*Tmp_Val(1)       &
                        + Reusable_Vals(3)*Tmp_Val(2)       &
                        + Reusable_Vals(4)*Tmp_Val(3)  ), Kind = idp    )





END DO ! rd Loop
END DO ! td Loop
END DO ! pd Loop
END DO ! re Loop
END DO ! te Loop
END DO ! pe Loop

END SUBROUTINE Poseidon_Return_Kij_Native


!+101+######################################################################################!
!                                                                                           !
!       Poseidon_Return_ConFactor - Uses the results from the routine,                      !
!                                   Poseidon_XCFC_Run_Part1(), to calculate the value       !
!                                   of the Conformal Factor at the desired locations.       !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Return_Kij_Native_Caller( Return_Kij       )


REAL(idp),  DIMENSION( Caller_Quad_DOF,                      &
                      Num_R_Elements,                   &
                      Num_T_Elements,                   &
                      Num_P_Elements,                   &
                      1:6           ),  INTENT(OUT)    ::  Return_Kij




INTEGER, DIMENSION(3)                                                   ::  NE

NE = [ Num_R_Elements, Num_T_Elements, Num_P_Elements ]


CALL Poseidon_Return_Kij_Native(  NE,                   &
                                  Caller_NQ,            &
                                  Caller_RQ_xlocs,      &
                                  Caller_TQ_xlocs,      &
                                  Caller_PQ_xlocs,      &
                                  Caller_xL(1),         &
                                  Caller_xL(2),         &
                                  Return_Kij            )


END SUBROUTINE Poseidon_Return_Kij_Native_Caller




#ifdef POSEIDON_AMREX_FLAG
!+101+######################################################################################!
!                                                                                           !
!       Poseidon_Return_ConFactor_AMReX                                                     !
!           - Uses the results from the routine, Poseidon_XCFC_Run_Part1(), to calculate    !
!             the value of the Conformal Factor at the desired locations, and fill an       !
!             AMReX multifab with the data.                                                 !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Return_Kij_AMReX( NQ,                     &
                                   RQ_Input,               &
                                   TQ_Input,               &
                                   PQ_Input,               &
                                   Left_Limit,             &
                                   Right_Limit,            &
                                   nLevels,                &
                                   MF_Results,             &
                                   FillGhostCells_Option   )


INTEGER,    DIMENSION(3),                       INTENT(IN)      ::  NQ
REAL(idp),  DIMENSION(NQ(1)),                   INTENT(IN)      ::  RQ_Input
REAL(idp),  DIMENSION(NQ(2)),                   INTENT(IN)      ::  TQ_Input
REAL(idp),  DIMENSION(NQ(3)),                   INTENT(IN)      ::  PQ_Input
REAL(idp),                                      INTENT(IN)      ::  Left_Limit
REAL(idp),                                      INTENT(IN)      ::  Right_Limit

INTEGER,                                        INTENT(IN)      ::  nLevels
TYPE(amrex_multifab),                           INTENT(INOUT)   ::  MF_Results(0:nLevels-1)

LOGICAL,                        OPTIONAL,   INTENT(IN)      ::  FillGhostCells_Option


 INTEGER                                                 ::  iRE
 INTEGER                                                 ::  re, te, pe
 INTEGER                                                 ::  rd, td, pd, tpd
 INTEGER                                                 ::  d
 INTEGER                                                 ::  i, j, lvl, lm
 INTEGER                                                 ::  nComp
 INTEGER                                                 ::  Here, There

 TYPE(amrex_mfiter)                                              ::  mfi
 TYPE(amrex_box)                                                 ::  Box
 TYPE(amrex_imultifab)                                           ::  Level_Mask
 INTEGER, DIMENSION(3)                                           ::  iEL, iEU
 INTEGER,    CONTIGUOUS, POINTER                                 ::  Mask_PTR(:,:,:,:)
 REAL(idp),  CONTIGUOUS, POINTER                                 ::  Result_PTR(:,:,:,:)



 REAL(idp)                                                       ::  Quad_Span
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
 REAL(idp), DIMENSION(3)                                      ::  TMP_Val
 REAL(idp), DIMENSION(3,3)                                    ::  TMP_Drv

 REAL(idp), DIMENSION(4)                                      ::  Reusable_Vals
 REAL(idp), DIMENSION(6)                                      ::  Tmp_A
  REAL(idp)                                                       ::  TMP_Val_Psi

 INTEGER                                                         ::  iVB
 INTEGER, DIMENSION(3)                                           ::  iU

 INTEGER                                                         ::  Num_DOF
 
 REAL(idp), DIMENSION(:,:,:), ALLOCATABLE                        ::  Caller_LPT
 
INTEGER, DIMENSION(3)                                       ::  iEL_A, iEU_A
LOGICAL                                                     ::  FillGhostCells
INTEGER,        DIMENSION(1:3)                              ::  nGhost_Vec


INTEGER,    DIMENSION(3)                                    ::  iNE
REAL(idp),  DIMENSION(0:AMReX_Max_Grid_Size(2)-1)           ::  tlocs_subarray
REAL(idp),  DIMENSION(0:AMReX_Max_Grid_Size(3)-1)           ::  plocs_subarray

REAL(idp),  DIMENSION(1:NQ(2),1:LM_Short_Length,0:AMReX_Max_Grid_Size(2)-1) ::  Plm_Table
REAL(idp),  DIMENSION(1:NQ(2),1:LM_Short_Length,0:AMReX_Max_Grid_Size(2)-1) ::  Plm_dt_Table

REAL(idp),  DIMENSION(1:NQ(3),1:LM_Length,0:AMReX_Max_Grid_Size(3)-1)       ::  Am_Table
REAL(idp),  DIMENSION(1:NQ(3),1:LM_Length,0:AMReX_Max_Grid_Size(3)-1)       ::  Am_dp_Table

REAL(idp),  DIMENSION(1:LM_Length, 1:NQ(2)*NQ(3) )                          ::  Slm_Elem_Table
REAL(idp),  DIMENSION(1:LM_Length, 1:NQ(2)*NQ(3) )                          ::  Slm_Elem_dt_Table
REAL(idp),  DIMENSION(1:LM_Length, 1:NQ(2)*NQ(3) )                          ::  Slm_Elem_dp_Table


IF ( PRESENT(FillGhostCells_Option) ) THEN
    FillGhostCells = FillGhostCells_Option
ELSE
    FillGhostCells = .FALSE.
END IF


 iVB   = iVB_X
 iU(1) = iU_X1
 iU(2) = iU_X2
 iU(3) = iU_X3

 Quad_Span = Right_Limit - Left_Limit

 CUR_RX_LOCS = 2.0_idp * ( RQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp
 CUR_TX_LOCS = 2.0_idp * ( TQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp
 CUR_PX_LOCS = 2.0_idp * ( PQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp

 Num_DOF = NQ(1)*NQ(2)*NQ(3)

 ! nComp order : ij = 11, 12, 13, 22, 23, 33
 Gamma(1) = 1.0_idp
 Christoffel = 0.0_idp

ALLOCATE( Caller_LPT(0:1,0:DEGREE,1:NQ(1)))

 DO lvl = nLevels-1,0,-1

    IF ( FillGhostCells ) THEN
        nGhost_Vec = MF_Results(lvl)%nghostvect()
    ELSE
        nGhost_Vec = 0
    END IF

     DROT = 0.5_idp * Level_dx(lvl,1)
     DTOT = 0.5_idp * Level_dx(lvl,2)
     DPOT = 0.5_idp * Level_dx(lvl,3)

     !
     !   MakeFineMask
     !
     IF ( lvl < AMReX_Num_Levels-1 ) THEN
         CALL AMReX_MakeFineMask(  Level_Mask,               &
                                   MF_Results(lvl)%ba,       &
                                   MF_Results(lvl)%dm,       &
                                   nGhost_Vec,               &
                                   MF_Results(lvl+1)%ba,     &
                                   iLeaf, iTrunk            )
     ELSE
         ! Create Level_Mask all equal to 1
         CALL amrex_imultifab_build( Level_Mask,             &
                                     MF_Results(lvl)%ba,     &
                                     MF_Results(lvl)%dm,     &
                                     1,                      &  ! ncomp = 1
                                     0                       )  ! nghost = 0
         CALL Level_Mask%SetVal(iLeaf)
     END IF


     CALL amrex_mfiter_build(mfi, MF_Results(lvl), tiling = .true. )

     DO WHILE(mfi%next())

        Result_PTR => MF_Results(lvl)%dataPtr(mfi)
        Mask_PTR   => Level_Mask%dataPtr(mfi)

        Box = mfi%tilebox()
        nComp =  MF_Results(lvl)%ncomp()

        iEL_A = Box%lo
        iEU_A = Box%hi
        
        iEL = iEL_A-nGhost_Vec
        iEU = iEU_A+nGhost_Vec

        IF ( ANY( iEL < 0 ) ) THEN
            ! Reflecting Conditions
            iEL = iEL_A
        END IF

        IF ( ANY( iEU .GE. (2**lvl)*iNE_Base(1) ) ) THEN
            iEU = iEU_A
        END IF
        
        iNE = iEU-iEL+1
        
        
        DO te = iEL(2),iEU(2)
            tlocs_subarray(te-iEL(2)) = Level_dx(lvl,2)*te
        END DO
        DO pe = iEL(3),iEU(3)
            plocs_subarray(pe-iEL(3)) = Level_dx(lvl,3)*pe
        END DO
        
        
        ! Initialize Am Table
        CALL Initialize_Am_Tables(  NQ(3),                      &
                                    PQ_Input,                   &
                                    L_Limit,                    &
                                    iNE(3),                     &
                                    [iEL(3), iEU(3)],           &
                                    plocs_subarray(0:iNE(3)-1), &
                                    Am_Table,                   &
                                    Am_dp_Table                 )

        ! Initialize Plm Table
        CALL Initialize_Plm_Tables( NQ(2),                      &
                                    TQ_Input,                   &
                                    L_Limit,                    &
                                    LM_Short_Length,            &
                                    iNE(2),                     &
                                    [iEL(2), iEU(2)],           &
                                    tlocs_subarray(0:iNE(2)-1), &
                                    Plm_Table,                  &
                                    Plm_dt_Table                )

        DO re = iEL(1),iEU(1)
        DO te = iEL(2),iEU(2)
        DO pe = iEL(3),iEU(3)


        IF ( Mask_PTR(RE,TE,PE,1) == iLeaf ) THEN
             iRE = FEM_Elem_Map(re,lvl)

             CALL Initialize_Slm_Tables_on_Elem( te, pe,             &
                                                NQ(2), NQ(3),           &
                                                iNE,                &
                                                iEL,                &
                                                Plm_Table,          &
                                                Plm_dt_Table,       &
                                                Am_Table,           &
                                                Am_dp_Table,        &
                                                Slm_Elem_Table,     &
                                                Slm_Elem_dt_Table,  &
                                                Slm_Elem_dp_Table   )
     
            Cur_R_Locs(:) = DROT * (CUR_RX_LOCS(:)+1.0_idp + 2.0_idp*re)
            Cur_T_Locs(:) = DTOT * (CUR_TX_LOCS(:)+1.0_idp + 2.0_idp*te)
            Cur_P_Locs(:) = DPOT * (CUR_PX_LOCS(:)+1.0_idp + 2.0_idp*pe)


            DO rd = 1,NQ(1)
                Caller_LPT(0,:,rd)=Lagrange_Poly(CUR_RX_LOCS(rd),DEGREE,FEM_Node_xlocs)
                Caller_LPT(1,:,rd)=Lagrange_Poly_Deriv(CUR_RX_LOCS(rd),DEGREE,FEM_Node_xlocs)
            END DO
            
             DO pd = 1,NQ(3)
             DO td = 1,NQ(2)
             DO rd = 1,NQ(1)

                 tpd = Map_To_tpd(td,pd)
                 Tmp_Val = 0.0_idp
                 Tmp_Drv = 0.0_idp

                 TMP_Val_Psi = 0.0_idp
                 DO lm = 1,LM_Length
                 DO d  = 0,DEGREE

                     Here  = Map_To_FEM_Node(RE-1,d)


                     TMP_Val_Psi = TMP_Val_Psi                      &
                             + dVA_Coeff_Vector(Here,lm,iU_CF)      &
                             * Caller_LPT(0,d,rd)                   &
                             * Slm_Elem_Table( lm, tpd )
                    
                 END DO
                 END DO




                 DO i  = 1,3
                 DO d  = 0,DEGREE

                     Here  = FP_Array_Map_TypeB(iU(i),iVB,iRE,d,1)                      ! Check iRE index
                     There = FP_Array_Map_TypeB(iU(i),iVB,iRE,d,LM_Length)


                     TMP_Val(i) = TMP_Val(i)                                  &
                             + SUM( dVB_Coeff_Vector( Here:There, iVB )      &
                                     * Slm_Elem_Table( :, tpd )   )       &
                             * Caller_LPT(0,d,rd)


                     TMP_Drv(1,i) = TMP_Drv(1,i)                              &
                                + SUM( dVB_Coeff_Vector( Here:There, iVB )   &
                                      * Slm_Elem_Table( :, tpd  )     )    &
                                * Caller_LPT(1,d,rd)             &
                                / DROT


                     TMP_Drv(2,i) = TMP_Drv(2,i)                              &
                                + SUM( dVB_Coeff_Vector( Here:There, iVB )   &
                                      * Slm_Elem_dt_Table( :, tpd )   )    &
                                * Caller_LPT(0,d,rd)

                     TMP_Drv(3,i) = TMP_Drv(3,i)                              &
                                + SUM( dVB_Coeff_Vector( Here:There, iVB )   &
                                      * Slm_Elem_dp_Table( :, tpd)   )    &
                                * Caller_LPT(0,d,rd)

                 END DO ! d Loop
                 END DO ! i Loop

                 Gamma(2) = 1.0_idp/(Cur_R_Locs(rd)*Cur_R_Locs(rd))
                 Gamma(3) = Gamma(2) * 1.0_idp/( DSIN(Cur_T_Locs(td))*DSIN(Cur_T_Locs(td)) )

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



                 


                 ! Ahat^11
                 Tmp_A(1) = Gamma(1)                                                         &
                          * ( 2.0_idp * Tmp_Drv(1,1) - Reusable_Vals(1)                      &
                            +(2.0_idp * Christoffel(1,1,1) - Reusable_Vals(2) )*Tmp_Val(1)   &
                            +(2.0_idp * Christoffel(1,1,2) - Reusable_Vals(3) )*Tmp_Val(2)   &
                            +(2.0_idp * Christoffel(1,1,3) - Reusable_Vals(4) )*Tmp_Val(3)   )

                 


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
                  


                 ! Ahat^22
                 Tmp_A(4) = Gamma(2)                                                             &
                          * ( 2.0_idp * Tmp_Drv(2,2) - Reusable_Vals(1)                          &
                            +(2.0_idp * Christoffel(2,2,1) - Reusable_Vals(2) ) * Tmp_Val(1)     &
                            +(2.0_idp * Christoffel(2,2,2) - Reusable_Vals(3) ) * Tmp_Val(2)     &
                            +(2.0_idp * Christoffel(2,2,3) - Reusable_Vals(4) ) * Tmp_Val(3)     )
                  


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
                 


                 ! Ahat^33
                 Tmp_A(6) = Gamma(3)                                                             &
                          * ( 2.0_idp * Tmp_Drv(3,3) - Reusable_Vals(1)                          &
                            +(2.0_idp * Christoffel(3,3,1) - Reusable_Vals(2) ) * Tmp_Val(1)     &
                            +(2.0_idp * Christoffel(3,3,2) - Reusable_Vals(3) ) * Tmp_Val(2)     &
                            +(2.0_idp * Christoffel(3,3,3) - Reusable_Vals(4) ) * Tmp_Val(3)     )



                 
                 Result_PTR(re,te,pe,AMReX_nCOMP_Map( iU_K11, rd, td, pd, NQ )) = Tmp_A(1)/(Gamma(1)*Gamma(1)*Tmp_Val_Psi*Tmp_Val_Psi)
                 Result_PTR(re,te,pe,AMReX_nCOMP_Map( iU_K12, rd, td, pd, NQ )) = Tmp_A(2)/(Gamma(1)*Gamma(2)*Tmp_Val_Psi*Tmp_Val_Psi)
                 Result_PTR(re,te,pe,AMReX_nCOMP_Map( iU_K13, rd, td, pd, NQ )) = Tmp_A(3)/(Gamma(1)*Gamma(3)*Tmp_Val_Psi*Tmp_Val_Psi)
                 Result_PTR(re,te,pe,AMReX_nCOMP_Map( iU_K22, rd, td, pd, NQ )) = Tmp_A(4)/(Gamma(2)*Gamma(2)*Tmp_Val_Psi*Tmp_Val_Psi)
                 Result_PTR(re,te,pe,AMReX_nCOMP_Map( iU_K23, rd, td, pd, NQ )) = Tmp_A(5)/(Gamma(2)*Gamma(3)*Tmp_Val_Psi*Tmp_Val_Psi)
                 Result_PTR(re,te,pe,AMReX_nCOMP_Map( iU_K33, rd, td, pd, NQ )) = Tmp_A(6)/(Gamma(3)*Gamma(3)*Tmp_Val_Psi*Tmp_Val_Psi)

             END DO ! rd Loop
             END DO ! td Loop
             END DO ! pd Loop



         END IF

         END DO ! pe
         END DO ! te
         END DO ! re

     END DO

     CALL amrex_mfiter_destroy(mfi)
     CALL amrex_imultifab_destroy( Level_Mask )

 END DO ! lvl

END SUBROUTINE Poseidon_Return_Kij_AMReX






!+102+######################################################################################!
!                                                                                           !
!       Poseidon_Return_ConFactor_AMReX_Caller                                              !
!           - Uses the results from the routine, Poseidon_XCFC_Run_Part1(), to calculate    !
!             the value of the Conformal Factor at the locations set during initialization, !
!             and fill an AMReX multifab with the data.                                     !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Return_Kij_AMReX_Caller( MF_Results )

TYPE(amrex_multifab),   INTENT(INOUT)           ::  MF_Results(0:AMReX_Num_Levels-1)



CALL Poseidon_Return_Kij_AMReX( Caller_NQ,                      &
                                Caller_RQ_xlocs,                &
                                Caller_TQ_xlocs,                &
                                Caller_PQ_xlocs,                &
                                Caller_xL(1),                   &
                                Caller_xL(2),                   &
                                AMReX_Num_Levels,               &
                                MF_Results                      )


END SUBROUTINE Poseidon_Return_Kij_AMReX_Caller


#else

!+101+######################################################################################!
!                                                                                           !
!       Poseidon_Return_ConFactor_AMReX                                                     !
!           - Uses the results from the routine, Poseidon_XCFC_Run_Part1(), to calculate    !
!             the value of the Conformal Factor at the desired locations, and fill an       !
!             AMReX multifab with the data.                                                 !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Return_Kij_AMReX()
END SUBROUTINE Poseidon_Return_Kij_AMReX


!+102+######################################################################################!
!                                                                                           !
!       Poseidon_Return_ConFactor_AMReX_Caller                                              !
!           - Uses the results from the routine, Poseidon_XCFC_Run_Part1(), to calculate    !
!             the value of the Conformal Factor at the locations set during initialization, !
!             and fill an AMReX multifab with the data.                                     !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Return_Kij_AMReX_Caller( A_Difference )
INTEGER, INTENT(IN)         :: A_Difference
END SUBROUTINE Poseidon_Return_Kij_AMReX_Caller

#endif












 !+202+########################################################!
!                                                               !
!       Poseidon_Return_AMReX_Extrinsic_Curvature               !
!                                                               !
 !#############################################################!
PURE FUNCTION AMReX_nCOMP_Map( iU, rd, td, pd, NQ )

INTEGER, INTENT(IN)                         ::  iU
INTEGER, INTENT(IN)                         ::  rd
INTEGER, INTENT(IN)                         ::  td
INTEGER, INTENT(IN)                         ::  pd
INTEGER, DIMENSION(3),  INTENT(IN)          ::  NQ

INTEGER                                     ::  AMReX_nCOMP_Map

INTEGER                                     ::  Here
INTEGER                                     ::  Num_QP

! CF = 1
! LF = 2
! S1 = 3
! S2 = 4
! S3 = 5
! K11 = 6
! K12 = 7
! K13 = 8
! K22 = 9
! K23 = 10
! K33 = 11

Here = (pd-1)*NQ(1)*NQ(2)   &
     + (td-1)*NQ(1)         &
     + rd
Num_QP = NQ(1)*NQ(2)*NQ(3)

AMReX_nCOMP_Map = (iU-1)*Num_QP + Here


END FUNCTION AMReX_nCOMP_Map








END MODULE Poseidon_Return_Routines_Kij
