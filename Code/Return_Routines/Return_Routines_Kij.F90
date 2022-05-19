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
           ONLY :   Ylm_Values,                 &
                    Ylm_dt_Values,              &
                    Ylm_dp_Values,              &
                    Ylm_Elem_Values,        &
                    Ylm_Elem_dt_Values,     &
                    Ylm_Elem_dp_Values,              &
                    Lagrange_Poly_Table,        &
                    Level_DX

USE Variables_Mesh, &
           ONLY :  Num_R_Elements,         &
                   Num_T_Elements,         &
                   Num_P_Elements

USE Variables_Derived, &
           ONLY :  LM_LENGTH

USE Variables_FP, &
           ONLY :  FP_Coeff_Vector_A,  &
                   FP_Coeff_Vector_B

USE Variables_Mesh, &
           ONLY :  rlocs,              &
                   tlocs,              &
                   plocs

USE Variables_AMReX_Source, &
           ONLY :  iLeaf,                &
                   iTrunk

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
           ONLY :  Lagrange_Poly

USE Initialization_Tables, &
           ONLY :  Initialize_Normed_Legendre_Tables_On_Level,     &
                   Initialize_Ylm_Tables_On_Elem


USE Variables_Interface, &
           ONLY :  Caller_nLevels,                 &
                   Caller_NQ,                      &
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
           ONLY :  AMReX_Num_Levels

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

REAL(idp),  DIMENSION(NQ(1)*NQ(2)*NQ(3),NE(1),NE(2),NE(3),1:6), INTENT(OUT) ::  Return_Kij


INTEGER                                                         ::  re, te, pe
INTEGER                                                         ::  rd, td, pd, tpd
INTEGER                                                         ::  d
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


DO pe = 1,NE(3)-1
DO te = 1,NE(2)-1
DO re = 1,NE(1)-1

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
    Tmp_A(1) = Gamma(1)                                                         &
             * ( 2.0_idp * Tmp_Drv(1,1) - Reusable_Vals(1)                      &
               +(2.0_idp * Christoffel(1,1,1) - Reusable_Vals(2) )*Tmp_Val(1)   &
               +(2.0_idp * Christoffel(1,1,2) - Reusable_Vals(3) )*Tmp_Val(2)   &
               +(2.0_idp * Christoffel(1,1,3) - Reusable_Vals(4) )*Tmp_Val(3)   )

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
    Tmp_A(4) = Gamma(2)                                                             &
             * ( 2.0_idp * Tmp_Drv(2,2) - Reusable_Vals(1)                          &
               +(2.0_idp * Christoffel(2,2,1) - Reusable_Vals(2) ) * Tmp_Val(1)     &
               +(2.0_idp * Christoffel(2,2,2) - Reusable_Vals(3) ) * Tmp_Val(2)     &
               +(2.0_idp * Christoffel(2,2,3) - Reusable_Vals(4) ) * Tmp_Val(3)     )
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
    Tmp_A(6) = Gamma(3)                                                             &
             * ( 2.0_idp * Tmp_Drv(3,3) - Reusable_Vals(1)                          &
               +(2.0_idp * Christoffel(3,3,1) - Reusable_Vals(2) ) * Tmp_Val(1)     &
               +(2.0_idp * Christoffel(3,3,2) - Reusable_Vals(3) ) * Tmp_Val(2)     &
               +(2.0_idp * Christoffel(3,3,3) - Reusable_Vals(4) ) * Tmp_Val(3)     )
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
                                   MF_Results              )


INTEGER,    DIMENSION(3),                       INTENT(IN)      ::  NQ
REAL(idp),  DIMENSION(NQ(1)),                   INTENT(IN)      ::  RQ_Input
REAL(idp),  DIMENSION(NQ(2)),                   INTENT(IN)      ::  TQ_Input
REAL(idp),  DIMENSION(NQ(3)),                   INTENT(IN)      ::  PQ_Input
REAL(idp),                                      INTENT(IN)      ::  Left_Limit
REAL(idp),                                      INTENT(IN)      ::  Right_Limit

INTEGER,                                        INTENT(IN)      ::  nLevels
TYPE(amrex_multifab),                           INTENT(INOUT)   ::  MF_Results(0:nLevels-1)



 INTEGER                                                 ::  iRE
 INTEGER                                                 ::  re, te, pe
 INTEGER                                                 ::  rd, td, pd, tpd
 INTEGER                                                 ::  d
 INTEGER                                                 ::  i, j, lvl
 INTEGER                                                 ::  nComp
 INTEGER                                                 ::  Here, There

 TYPE(amrex_mfiter)                                              ::  mfi
 TYPE(amrex_box)                                                 ::  Box
 TYPE(amrex_imultifab)                                           ::  Level_Mask
 INTEGER, DIMENSION(3)                                           ::  iEL, iEU
 INTEGER,    CONTIGUOUS, POINTER                                 ::  Mask_PTR(:,:,:,:)
 REAL(idp),  CONTIGUOUS, POINTER                                 ::  Result_PTR(:,:,:,:)



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

 INTEGER                                                         ::  iVB
 INTEGER, DIMENSION(3)                                           ::  iU

 INTEGER                                                         ::  Num_DOF


 iVB   = iVB_X
 iU(1) = iU_X1
 iU(2) = iU_X2
 iU(3) = iU_X3

 Quad_Span = Right_Limit - Left_Limit

 Local_Locations = Initialize_LGL_Quadrature_Locations(DEGREE)
 CUR_RX_LOCS = 2.0_idp * ( RQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp
 CUR_TX_LOCS = 2.0_idp * ( TQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp
 CUR_PX_LOCS = 2.0_idp * ( PQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp

 Num_DOF = NQ(1)*NQ(2)*NQ(3)

 ! nComp order : ij = 11, 12, 13, 22, 23, 33
 Gamma(1) = 1.0_idp
 Christoffel = 0.0_idp


 DO lvl = nLevels-1,0,-1


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

     CALL Initialize_Normed_Legendre_Tables_on_Level( iEU, iEL, lvl )

     CALL amrex_mfiter_build(mfi, MF_Results(lvl), tiling = .true. )

     DO WHILE(mfi%next())

         Result_PTR => MF_Results(lvl)%dataPtr(mfi)
         Mask_PTR   => Level_Mask%dataPtr(mfi)

         Box = mfi%tilebox()
         nComp =  MF_Results(lvl)%ncomp()

         iEL = Box%lo
         iEU = Box%hi

         DO re = iEL(1),iEU(1)
         DO te = iEL(2),iEU(2)
         DO pe = iEL(3),iEU(3)


         IF ( Mask_PTR(RE,TE,PE,1) == iLeaf ) THEN
             iRE = FEM_Elem_Map(re,lvl)

             CALL Initialize_Ylm_Tables_on_Elem( te, pe, iEL, lvl )
     
             Cur_R_Locs(:) = DROT * (CUR_RX_LOCS(:)+1.0_idp + 2.0_idp*re)
             Cur_T_Locs(:) = DTOT * (CUR_TX_LOCS(:)+1.0_idp + 2.0_idp*te)
             Cur_P_Locs(:) = DPOT * (CUR_PX_LOCS(:)+1.0_idp + 2.0_idp*pe)



             DO pd = 1,NQ(3)
             DO td = 1,NQ(2)
             DO rd = 1,NQ(1)

                 tpd = Map_To_tpd(td,pd)
                 LagP = Lagrange_Poly(CUR_RX_LOCS(rd),DEGREE,Local_Locations)
                 Tmp_Val = 0.0_idp
                 Tmp_Drv = 0.0_idp

                 DO i  = 1,3
                 DO d  = 0,DEGREE

                     Here  = FP_Array_Map_TypeB(iU(i),iVB,iRE,d,1)
                     There = FP_Array_Map_TypeB(iU(i),iVB,iRE,d,LM_Length)


                     TMP_Val(i) = TMP_Val(i)                                  &
                             + SUM( FP_Coeff_Vector_B( Here:There, iVB )      &
                                     * Ylm_Elem_Values( :, tpd )   )       &
                             * Lagrange_Poly_Table( d, rd, 0 )


                     TMP_Drv(1,i) = TMP_Drv(1,i)                              &
                                + SUM( FP_Coeff_Vector_B( Here:There, iVB )   &
                                      * Ylm_Elem_Values( :, tpd  )     )    &
                                * Lagrange_Poly_Table( d, rd, 1 )             &
                                / DROT


                     TMP_Drv(2,i) = TMP_Drv(2,i)                              &
                                + SUM( FP_Coeff_Vector_B( Here:There, iVB )   &
                                      * Ylm_Elem_dt_Values( :, tpd )   )    &
                                * Lagrange_Poly_Table( d, rd, 0)

                     TMP_Drv(3,i) = TMP_Drv(3,i)                              &
                                + SUM( FP_Coeff_Vector_B( Here:There, iVB )   &
                                      * Ylm_Elem_dp_Values( :, tpd)   )    &
                                * Lagrange_Poly_Table( d, rd, 0)

                 END DO ! d Loop
                 END DO ! i Loop

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



                 
                 Result_PTR(re,te,pe,AMReX_nCOMP_Map( iU_K11, rd, td, pd, NQ )) = Tmp_A(1)/(Gamma(1)*Gamma(1))
                 Result_PTR(re,te,pe,AMReX_nCOMP_Map( iU_K12, rd, td, pd, NQ )) = Tmp_A(2)/(Gamma(1)*Gamma(2))
                 Result_PTR(re,te,pe,AMReX_nCOMP_Map( iU_K13, rd, td, pd, NQ )) = Tmp_A(3)/(Gamma(1)*Gamma(3))
                 Result_PTR(re,te,pe,AMReX_nCOMP_Map( iU_K22, rd, td, pd, NQ )) = Tmp_A(4)/(Gamma(2)*Gamma(2))
                 Result_PTR(re,te,pe,AMReX_nCOMP_Map( iU_K23, rd, td, pd, NQ )) = Tmp_A(5)/(Gamma(2)*Gamma(3))
                 Result_PTR(re,te,pe,AMReX_nCOMP_Map( iU_K33, rd, td, pd, NQ )) = Tmp_A(6)/(Gamma(3)*Gamma(3))

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

TYPE(amrex_multifab),   INTENT(INOUT)           ::  MF_Results(0:Caller_nLevels-1)



CALL Poseidon_Return_Kij_AMReX( Caller_NQ,                      &
                                Caller_RQ_xlocs,                &
                                Caller_TQ_xlocs,                &
                                Caller_PQ_xlocs,                &
                                Caller_xL(1),                   &
                                Caller_xL(2),                   &
                                Caller_nLevels,                 &
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
