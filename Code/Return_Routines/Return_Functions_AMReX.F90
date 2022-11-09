   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Return_Functions_AMReX                                                !##!
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
            ONLY :  Ylm_Elem_Values,            &
                    Ylm_Elem_dt_Values,         &
                    Ylm_Elem_dp_Values,         &
                    Lagrange_Poly_Table,        &
                    Level_DX


USE Variables_Derived, &
            ONLY :  LM_LENGTH

USE Variables_Vectors, &
            ONLY :  cVA_Coeff_Vector,      &
                    cVB_Coeff_Vector

USE Variables_Mesh, &
            ONLY :  rlocs,              &
                    tlocs,              &
                    plocs,              &
                    iNE_Base


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

USE Maps_X_Space, &
            ONLY :  Map_To_X_Space

USE Functions_Quadrature, &
            ONLY :  Initialize_LGL_Quadrature_Locations

USE Functions_Math, &
            ONLY :  Lagrange_Poly

USE Functions_Translation_Matrix_Module, &
            ONLY :  Create_Translation_Matrix

USE Initialization_Tables, &
            ONLY :  Initialize_Normed_Legendre_Tables_On_Level,     &
                    Initialize_Ylm_Tables_On_Elem


USE Variables_Interface, &
            ONLY :  Caller_NQ,                      &
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
                                        MF_Results,             &
                                        FillGhostCells_Option   )


INTEGER,                                    INTENT(IN)      ::  iU
INTEGER,    DIMENSION(3),                   INTENT(IN)      ::  NQ
REAL(idp),  DIMENSION(NQ(1)),               INTENT(IN)      ::  RQ_Input
REAL(idp),  DIMENSION(NQ(2)),               INTENT(IN)      ::  TQ_Input
REAL(idp),  DIMENSION(NQ(3)),               INTENT(IN)      ::  PQ_Input
REAL(idp),                                  INTENT(IN)      ::  Left_Limit
REAL(idp),                                  INTENT(IN)      ::  Right_Limit

INTEGER,                                    INTENT(IN)      ::  nLevels
TYPE(amrex_multifab),                       INTENT(INOUT)   ::  MF_Results(0:nLevels-1)

LOGICAL,                        OPTIONAL,   INTENT(IN)      ::  FillGhostCells_Option


INTEGER                                                     ::  Num_DOF
INTEGER                                                     ::  Output_Here

INTEGER                                                     ::  iRE
INTEGER                                                     ::  re, te, pe
INTEGER                                                     ::  rd, td, pd, tpd
INTEGER                                                     ::  d, Here

REAL(KIND = idp)                                            ::  Quad_Span
REAL(KIND = idp), DIMENSION(0:DEGREE)                       ::  LagP
REAL(KIND = idp), DIMENSION(1:NQ(1))                        ::  Cur_RX_Locs
COMPLEX(KIND = idp)                                         ::  TMP_U_Value
INTEGER                                                     ::  Current_Location

INTEGER,        DIMENSION(1:3)                              ::  nGhost_Vec


INTEGER                                                     ::  lvl
TYPE(amrex_mfiter)                                          ::  mfi
TYPE(amrex_box)                                             ::  Box
TYPE(amrex_imultifab)                                       ::  Level_Mask
INTEGER                                                     ::  nComp

INTEGER,    CONTIGUOUS, POINTER                             ::  Mask_PTR(:,:,:,:)
REAL(idp),  CONTIGUOUS, POINTER                             ::  Result_PTR(:,:,:,:)
REAL(idp),  DIMENSION(1:Local_Quad_DOF)                     ::  Var_Holder_Elem
REAL(idp),  DIMENSION(:,:), ALLOCATABLE                     ::  Translation_Matrix


INTEGER, DIMENSION(3)                                       ::  iEL, iEU
INTEGER, DIMENSION(3)                                       ::  iEL_A, iEU_A
LOGICAL                                                     ::  FillGhostCells



IF ( PRESENT(FillGhostCells_Option) ) THEN
    FillGhostCells = FillGhostCells_Option
ELSE
    FillGhostCells = .FALSE.
END IF


Quad_Span = Right_Limit - Left_Limit

Cur_RX_Locs = 2.0_idp * ( RQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp

Num_DOF = NQ(1)*NQ(2)*NQ(3)

Allocate( Translation_Matrix(1:Local_Quad_DOF,1:Num_DOF))

Translation_Matrix = Create_Translation_Matrix( [Num_R_Quad_Points, Num_T_Quad_Points, Num_P_Quad_Points ],            &
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
                                        Num_DOF               )



DO lvl = nLevels-1,0,-1


    IF ( FillGhostCells ) THEN
        nGhost_Vec = MF_Results(lvl)%nghostvect()
    ELSE
        nGhost_Vec = 0
    END IF
    
!    PRINT*,"In Return Function, nGhost",nGhost_Vec(1),FillGhostCells
    
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
                                    nGhost_Vec(1)           )
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
                    LagP = Lagrange_Poly(Int_R_Locations(rd),DEGREE,FEM_Node_xlocs)
                    Tmp_U_Value = 0.0_idp

                    
                    DO d = 0,DEGREE
                        Current_Location = Map_To_FEM_Node(iRE,d)
                        Tmp_U_Value = Tmp_U_Value                                    &
                                    + SUM( cVA_Coeff_Vector(Current_Location,:,iU)  &
                                            * Ylm_Elem_Values( :, tpd )            ) &
                                    * LagP(d)

                    END DO ! d Loop

                    Here = Quad_Map(rd, td, pd, NQ(1), NQ(2),NQ(3))
                    Var_Holder_Elem(Here) = Tmp_U_Value
                END DO ! rd
                END DO ! td
                END DO ! pd



                DO Output_Here = 1,Num_DOF

                    Here = (iU-1)*Num_DOF+Output_Here

                    Result_PTR(re,te,pe,Here) = DOT_PRODUCT( Translation_Matrix(:,Output_Here), &
                                                             Var_Holder_Elem(:)         )

                END DO ! Output_Here

            END IF !  Mask_PTR(RE,TE,PE,1) == iLeaf

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

INTEGER                                                 ::  Num_DOF
INTEGER                                                 ::  Output_Here
INTEGER                                                 ::  iRE
INTEGER                                                 ::  re, te, pe
INTEGER                                                 ::  rd, td, pd, tpd
INTEGER                                                 ::  d
INTEGER                                                 ::  lvl

REAL(KIND = idp)                                        ::  Quad_Span
REAL(KIND = idp), DIMENSION(0:DEGREE)                   ::  LagP
REAL(KIND = idp), DIMENSION(1:NQ(1))                    ::  Cur_RX_Locs
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
REAL(idp),  DIMENSION(:,:), ALLOCATABLE         ::  Translation_Matrix

Quad_Span = Right_Limit - Left_Limit

Cur_RX_Locs = 2.0_idp * ( RQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp

Num_DOF = NQ(1)*NQ(2)*NQ(3)

Allocate( Translation_Matrix(1:Local_Quad_DOF,1:Num_DOF))

Translation_Matrix = Create_Translation_Matrix( [Num_R_Quad_Points, Num_T_Quad_Points, Num_P_Quad_Points ],            &
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
                                        Num_DOF               )




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
                    LagP = Lagrange_Poly(Cur_RX_Locs(rd),DEGREE,FEM_Node_xlocs)
                    Tmp_U_Value = 0.0_idp

                    DO d = 0,DEGREE

                        Here  = FP_Array_Map_TypeB(iU,iVB,iRE,d,1)
                        There = FP_Array_Map_TypeB(iU,iVB,iRE,d,LM_Length)
                        
                        Tmp_U_Value = Tmp_U_Value                                   &
                                    + SUM( cVB_Coeff_Vector(Here:There,iVB)        &
                                            * Ylm_Elem_Values( :, tpd ) )   &
                                    * LagP(d)
                                    

                    END DO ! d Loop

                    Here = Quad_Map(rd, td, pd, NQ(1), NQ(2),NQ(3))
                    Var_Holder_Elem(Here) = Tmp_U_Value

                END DO ! rd
                END DO ! td
                END DO ! pd



                DO Output_Here = 1,Num_DOF

                    Here = (iU-1)*Num_DOF+Output_Here

                    Result_PTR(re,te,pe,Here) = DOT_PRODUCT( Translation_Matrix(:,Output_Here), &
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

Here = Quad_Map(rd, td, pd, NQ(1), NQ(2),NQ(3))
Num_QP = NQ(1)*NQ(2)*NQ(3)

AMReX_nCOMP_Map = (iU-1)*Num_QP + Here


END FUNCTION AMReX_nCOMP_Map










END MODULE Return_Functions_AMReX
