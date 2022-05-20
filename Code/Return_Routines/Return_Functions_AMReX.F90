   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Return_Functions_AMReX                                          	     !##!
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
            ONLY :  Ylm_Elem_Values,        &
                    Ylm_Elem_dt_Values,     &
                    Ylm_Elem_dp_Values,              &
                    Lagrange_Poly_Table,        &
                    Level_DX


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

USE Maps_X_Space, &
            ONLY :  Map_To_X_Space

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
                    Caller_xL,                      &
                    Caller_RQ_xlocs,                &
                    Caller_TQ_xlocs,                &
                    Caller_PQ_xlocs

USE Variables_Quadrature, &
            ONLY :  Num_R_Quad_Points,      &
                    Num_T_Quad_Points,      &
                    Num_P_Quad_Points,      &
                    Num_TP_Quad_Points,     &
                    Local_Quad_DOF,         &
                    xLeftLimit,             &
                    xRightLimit,            &
                    Int_R_Locations,        &
                    Int_T_Locations,        &
                    Int_P_Locations,        &
                    Int_R_Weights,          &
                    Int_T_Weights,          &
                    Int_P_Weights,          &
                    Int_TP_Weights


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

IMPLICIT NONE




CONTAINS


#ifdef POSEIDON_AMREX_FLAG
 !+201+########################################################!
!                                                               !
!       Poseidon_Return_AMReX_Type_A                            !
!                                                               !
 !#############################################################!
SUBROUTINE Poseidon_Return_AMReX_Type_A(iU,                     &
                                        NQ,                     &
                                        RQ_Input,               &
                                        TQ_Input,               &
                                        PQ_Input,               &
                                        Left_Limit,             &
                                        Right_Limit,            &
                                        nLevels,                &
                                        MF_Results              )


INTEGER,                                    INTENT(IN)  ::  iU
INTEGER,    DIMENSION(3),                   INTENT(IN)  ::  NQ
REAL(idp),  DIMENSION(NQ(1)),               INTENT(IN)  ::  RQ_Input
REAL(idp),  DIMENSION(NQ(2)),               INTENT(IN)  ::  TQ_Input
REAL(idp),  DIMENSION(NQ(3)),               INTENT(IN)  ::  PQ_Input
REAL(idp),                                  INTENT(IN)  ::  Left_Limit
REAL(idp),                                  INTENT(IN)  ::  Right_Limit

INTEGER,                                    INTENT(IN)     ::  nLevels
TYPE(amrex_multifab),                       INTENT(INOUT)  ::  MF_Results(0:nLevels-1)


INTEGER                                                 ::  DOF_Input
INTEGER                                                 ::  Output_Here

INTEGER                                                 ::  iRE
INTEGER                                                 ::  re, te, pe
INTEGER                                                 ::  rd, td, pd, tpd
INTEGER                                                 ::  d, lm, Here

REAL(KIND = idp)                                        ::  Quad_Span
REAL(KIND = idp), DIMENSION(0:DEGREE)                   ::  Local_Locations
REAL(KIND = idp), DIMENSION(0:DEGREE)                   ::  LagP
REAL(KIND = idp), DIMENSION(1:NQ(1))                    ::  CUR_X_LOCS
COMPLEX(KIND = idp)                                     ::  TMP_U_Value
INTEGER                                                 ::  Current_Location


INTEGER                                                 ::  lvl
TYPE(amrex_mfiter)                              ::  mfi
TYPE(amrex_box)                                 ::  Box
TYPE(amrex_imultifab)                           ::  Level_Mask
INTEGER                                                 ::  nComp

INTEGER,    CONTIGUOUS, POINTER                 ::  Mask_PTR(:,:,:,:)
REAL(idp),  CONTIGUOUS, POINTER                 ::  Result_PTR(:,:,:,:)
REAL(idp),  DIMENSION(1:Local_Quad_DOF)         ::  Var_Holder_Elem
REAL(idp),  DIMENSION(:,:), ALLOCATABLE         ::  TransMat


INTEGER, DIMENSION(3)                           ::  iEL, iEU

Quad_Span = Right_Limit - Left_Limit

Local_Locations = Initialize_LGL_Quadrature_Locations(DEGREE)
CUR_X_LOCS = 2.0_idp * ( RQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp

DOF_Input = NQ(1)*NQ(2)*NQ(3)

Allocate( TransMat(1:Local_Quad_DOF,1:DOF_Input))

TransMat = Create_Translation_Matrix( [Num_R_Quad_Points, Num_T_Quad_Points, Num_P_Quad_Points ],            &
                                        [xLeftLimit, xRightLimit ],            &
                                        Int_R_Locations,      &
                                        Int_R_Locations,      &
                                        Int_R_Locations,      &
                                        Local_Quad_DOF,         &
                                        NQ,          &
                                        [Left_Limit, Right_Limit],          &
                                        RQ_Input,    &
                                        TQ_Input,    &
                                        PQ_Input,    &
                                        DOF_Input               )



DO lvl = nLevels-1,0,-1


    !
    !   MakeFineMask
    !
    IF ( lvl < AMReX_Num_Levels-1 ) THEN
        CALL AMReX_MakeFineMask(  Level_Mask,               &
                                  MF_Results(lvl)%ba,        &
                                  MF_Results(lvl)%dm,        &
                                  MF_Results(lvl+1)%ba,      &
                                  iLeaf, iTrunk            )
    ELSE
        ! Create Level_Mask all equal to 1
        CALL amrex_imultifab_build( Level_Mask,             &
                                    MF_Results(lvl)%ba,      &
                                    MF_Results(lvl)%dm,      &
                                    1,                      &
                                    0                       )
        CALL Level_Mask%SetVal(iLeaf)
    END IF


    CALL amrex_mfiter_build(mfi, MF_Results(lvl), tiling = .true. )

    DO WHILE(mfi%next())

        Result_PTR => MF_Results(lvl)%dataPtr(mfi)
        Mask_PTR   => Level_Mask%dataPtr(mfi)

        Box = mfi%tilebox()
        nComp =  MF_Results(lvl)%ncomp()

        iEL = Box%lo
        iEU = Box%hi

        CALL Initialize_Normed_Legendre_Tables_on_Level( iEU, iEL, lvl )

        DO re = iEL(1),iEU(1)
        DO te = iEL(2),iEU(2)
        DO pe = iEL(3),iEU(3)
            
            IF ( Mask_PTR(RE,TE,PE,1) == iLeaf ) THEN
            iRE = FEM_Elem_Map(re,lvl)
            CALL Initialize_Ylm_Tables_on_Elem( te, pe, iEL, lvl )

            DO pd = 1,Num_P_Quad_Points
            DO td = 1,NUM_T_QUAD_POINTS
            DO rd = 1,NUM_R_QUAD_POINTS

                tpd = Map_To_tpd(td,pd)
                LagP = Lagrange_Poly(Int_R_Locations(rd),DEGREE,Local_Locations)
                Tmp_U_Value = 0.0_idp

                
                DO d = 0,DEGREE
                    Current_Location = Map_To_FEM_Node(iRE,d)
                    Tmp_U_Value = Tmp_U_Value                                    &
                                + SUM( FP_Coeff_Vector_A(Current_Location,:,iU)  &
                                        * Ylm_Elem_Values( :, tpd )            ) &
                                * LagP(d)

                END DO ! d Loop

                Here = (pd-1)*NQ(2)*NQ(1)           &
                    + (td-1)*NQ(1)                  &
                    + rd
                Var_Holder_Elem(Here) = Tmp_U_Value
    
            END DO ! rd
            END DO ! td
            END DO ! pd



            DO Output_Here = 1,DOF_Input

                Here = (iU-1)*DOF_Input+Output_Here

                Result_PTR(re,te,pe,Here) = DOT_PRODUCT( TransMat(:,Output_Here), &
                                                         Var_Holder_Elem(:)         )

            END DO

            END IF

        END DO ! pe
        END DO ! te
        END DO ! re

    END DO

    CALL amrex_mfiter_destroy(mfi)
    CALL amrex_imultifab_destroy( Level_Mask )

END DO ! lvl

END SUBROUTINE Poseidon_Return_AMReX_Type_A




 !+202+########################################################!
!                                                               !
!       Poseidon_Return_AMReX_Type_B                            !
!                                                               !
 !#############################################################!
SUBROUTINE Poseidon_Return_AMReX_Type_B(iU,                     &
                                        iVB,                    &
                                        NQ,                     &
                                        RQ_Input,               &
                                        TQ_Input,               &
                                        PQ_Input,               &
                                        Left_Limit,             &
                                        Right_Limit,            &
                                        nLevels,                &
                                        MF_Results              )


INTEGER,                                    INTENT(IN)  ::  iU
INTEGER,                                    INTENT(IN)  ::  iVB
INTEGER,    DIMENSION(3),                   INTENT(IN)  ::  NQ
REAL(idp),  DIMENSION(NQ(1)),               INTENT(IN)  ::  RQ_Input
REAL(idp),  DIMENSION(NQ(2)),               INTENT(IN)  ::  TQ_Input
REAL(idp),  DIMENSION(NQ(3)),               INTENT(IN)  ::  PQ_Input
REAL(idp),                                  INTENT(IN)  ::  Left_Limit
REAL(idp),                                  INTENT(IN)  ::  Right_Limit

INTEGER,                                    INTENT(IN)  ::  nLevels
TYPE(amrex_multifab),                       INTENT(INOUT)  ::  MF_Results(0:nLevels-1)

INTEGER                                                 ::  DOF_Input
INTEGER                                                 ::  Output_Here
INTEGER                                                 ::  iRE
INTEGER                                                 ::  re, te, pe
INTEGER                                                 ::  rd, td, pd, tpd
INTEGER                                                 ::  d
INTEGER                                                 ::  lvl

REAL(KIND = idp)                                        ::  Quad_Span
REAL(KIND = idp), DIMENSION(0:DEGREE)                   ::  Local_Locations
REAL(KIND = idp), DIMENSION(0:DEGREE)                   ::  LagP
REAL(KIND = idp), DIMENSION(1:NQ(1))                    ::  CUR_X_LOCS
COMPLEX(KIND = idp)                                     ::  TMP_U_Value
INTEGER                                                 ::  Here, There

TYPE(amrex_mfiter)                                      ::  mfi
TYPE(amrex_box)                                         ::  Box
TYPE(amrex_imultifab)                                   ::  Level_Mask
INTEGER, DIMENSION(3)                                   ::  iEL, iEU
INTEGER                                                 ::  nComp

INTEGER,    CONTIGUOUS, POINTER                         ::  Mask_PTR(:,:,:,:)
REAL(idp),  CONTIGUOUS, POINTER                         ::  Result_PTR(:,:,:,:)
REAL(idp),  DIMENSION(1:Local_Quad_DOF)         ::  Var_Holder_Elem
REAL(idp),  DIMENSION(:,:), ALLOCATABLE         ::  TransMat

Quad_Span = Right_Limit - Left_Limit

Local_Locations = Initialize_LGL_Quadrature_Locations(DEGREE)
CUR_X_LOCS = 2.0_idp * ( RQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp

DOF_Input = NQ(1)*NQ(2)*NQ(3)

Allocate( TransMat(1:Local_Quad_DOF,1:DOF_Input))

TransMat = Create_Translation_Matrix( [Num_R_Quad_Points, Num_T_Quad_Points, Num_P_Quad_Points ],            &
                                        [xLeftLimit, xRightLimit ],            &
                                        Int_R_Locations,      &
                                        Int_R_Locations,      &
                                        Int_R_Locations,      &
                                        Local_Quad_DOF,         &
                                        NQ,          &
                                        [Left_Limit, Right_Limit],          &
                                        RQ_Input,    &
                                        TQ_Input,    &
                                        PQ_Input,    &
                                        DOF_Input               )




DO lvl = nLevels-1,0,-1

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
                DO pd = 1,NQ(3)
                DO td = 1,NQ(2)
                DO rd = 1,NQ(1)

                    tpd = Map_To_tpd(td,pd)
                    LagP = Lagrange_Poly(CUR_X_LOCS(rd),DEGREE,Local_Locations)
                    Tmp_U_Value = 0.0_idp

                    DO d = 0,DEGREE

                        Here  = FP_Array_Map_TypeB(iU,iVB,iRE,d,1)
                        There = FP_Array_Map_TypeB(iU,iVB,iRE,d,LM_Length)
                        
                        Tmp_U_Value = Tmp_U_Value                                   &
                                    + SUM( FP_Coeff_Vector_B(Here:There,iVB)        &
                                            * Ylm_Elem_Values( :, tpd ) )   &
                                    * LagP(d)
                                    

                    END DO ! d Loop

                    Here = (pd-1)*NQ(2)*NQ(1)           &
                        + (td-1)*NQ(1)                  &
                        + rd
                    Var_Holder_Elem(Here) = Tmp_U_Value

                END DO ! rd
                END DO ! td
                END DO ! pd



                DO Output_Here = 1,DOF_Input

                    Here = (iU-1)*DOF_Input+Output_Here

                    Result_PTR(re,te,pe,Here) = DOT_PRODUCT( TransMat(:,Output_Here), &
                                                             Var_Holder_Elem(:)         )

                END DO

            END IF

        END DO ! pe
        END DO ! te
        END DO ! re

    END DO

    CALL amrex_mfiter_destroy(mfi)
    CALL amrex_imultifab_destroy( Level_Mask )

END DO ! lvl

END SUBROUTINE Poseidon_Return_AMReX_Type_B






#else



 !+201+########################################################!
!                                                               !
!       Poseidon_Return_AMReX_Type_A                            !
!                                                               !
 !#############################################################!
SUBROUTINE Poseidon_Return_AMReX_Type_A(  )

END SUBROUTINE Poseidon_Return_AMReX_Type_A


 !+202+########################################################!
!                                                               !
!       Poseidon_Return_AMReX_Type_B                            !
!                                                               !
 !#############################################################!
SUBROUTINE Poseidon_Return_AMReX_Type_B()

END SUBROUTINE Poseidon_Return_AMReX_Type_B



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

REAL(idp),  DIMENSION(Source_DOF,Dest_DOF)              ::  TransMat

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
                      R_Lag_Poly_Values(1:Source_NQ(1),Dest_R)      &
                    * T_Lag_Poly_Values(Source_T,Source_T)          &
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





END MODULE Return_Functions_AMReX
