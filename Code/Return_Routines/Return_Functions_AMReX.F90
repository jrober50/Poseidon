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
            ONLY :  Slm_Elem_Values,            &
                    Slm_Elem_dt_Values,         &
                    Slm_Elem_dp_Values,         &
                    Nlm_Values,                 &
                    Lagrange_Poly_Table,        &
                    Level_DX


USE Variables_Derived, &
            ONLY :  LM_Length,                  &
                    LM_Short_Length

USE Variables_Vectors, &
            ONLY :  dVA_Coeff_Vector,      &
                    dVB_Coeff_Vector

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

USE Initialization_Tables_Slm, &
            ONLY :  Initialize_Am_Table,            &
                    Initialize_Plm_Table,           &
                    Initialize_Slm_Table_on_Elem



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
            ONLY :  AMReX_Num_Levels,       &
                    AMReX_Max_Grid_Size

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

INTEGER                                                     ::  iRE
INTEGER                                                     ::  re, te, pe
INTEGER                                                     ::  rd, td, pd, tpd
INTEGER                                                     ::  d, Here
INTEGER,        DIMENSION(3)                                ::  iNE

REAL(idp)                                                   ::  Quad_Span
REAL(idp),      DIMENSION(0:DEGREE)                         ::  LagP
REAL(idp),      DIMENSION(1:NQ(1))                          ::  Cur_RX_Locs
REAL(idp)                                                   ::  TMP_U_Value
INTEGER                                                     ::  Current_Location

INTEGER,        DIMENSION(1:3)                              ::  nGhost_Vec


INTEGER                                                     ::  lvl
TYPE(amrex_mfiter)                                          ::  mfi
TYPE(amrex_box)                                             ::  Box
TYPE(amrex_imultifab)                                       ::  Level_Mask
INTEGER                                                     ::  nComp

INTEGER,    CONTIGUOUS, POINTER                             ::  Mask_PTR(:,:,:,:)
REAL(idp),  CONTIGUOUS, POINTER                             ::  Result_PTR(:,:,:,:)


INTEGER, DIMENSION(3)                                       ::  iEL, iEU
INTEGER, DIMENSION(3)                                       ::  iEL_A, iEU_A
LOGICAL                                                     ::  FillGhostCells

REAL(idp),  DIMENSION(0:AMReX_Max_Grid_Size(2)-1)           ::  tlocs_subarray
REAL(idp),  DIMENSION(0:AMReX_Max_Grid_Size(3)-1)           ::  plocs_subarray
    
REAL(idp),  DIMENSION(1:NQ(2),1:LM_Short_Length,0:AMReX_Max_Grid_Size(2)-1) ::  Plm_Table
REAL(idp),  DIMENSION(1:NQ(3),1:LM_Length,0:AMReX_Max_Grid_Size(3)-1)       ::  Am_Table
REAL(idp),  DIMENSION(1:LM_Length, 1:NQ(2)*NQ(3) )                          ::  Slm_Elem_Table


IF ( PRESENT(FillGhostCells_Option) ) THEN
    FillGhostCells = FillGhostCells_Option
ELSE
    FillGhostCells = .FALSE.
END IF


Quad_Span = Right_Limit - Left_Limit

Cur_RX_Locs = 2.0_idp * ( RQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp

Num_DOF = NQ(1)*NQ(2)*NQ(3)



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
                                  MF_Results(lvl)%ba,       &
                                  MF_Results(lvl)%dm,       &
                                  nGhost_Vec,               &
                                  MF_Results(lvl+1)%ba,     &
                                  iLeaf, iTrunk             )
    ELSE
        ! Create Level_Mask all equal to 1
        CALL amrex_imultifab_build( Level_Mask,             &
                                    MF_Results(lvl)%ba,     &
                                    MF_Results(lvl)%dm,     &
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
        
        iNE = iEU-iEL+1
        
        
        
        
        DO te = iEL(2),iEU(2)
            tlocs_subarray(te-iEL(2)) = Level_dx(lvl,2)*te
        END DO
        DO pe = iEL(3),iEU(3)
            plocs_subarray(pe-iEL(3)) = Level_dx(lvl,3)*pe
        END DO
        
        
        ! Initialize Am Table
        CALL Initialize_Am_Table(   NQ(3),                      &
                                    PQ_Input,                   &
                                    L_Limit,                    &
                                    iNE(3),                     &
                                    [iEL(3), iEU(3)],           &
                                    plocs_subarray(0:iNE(3)-1), &
                                    Am_Table                    )

        ! Initialize Plm Table
        CALL Initialize_Plm_Table(  NQ(2),                      &
                                    TQ_Input,                   &
                                    L_Limit,                    &
                                    LM_Short_Length,            &
                                    iNE(2),                     &
                                    [iEL(2), iEU(2)],           &
                                    tlocs_subarray(0:iNE(2)-1), &
                                    Plm_Table                   )

        DO re = iEL(1),iEU(1)
        DO te = iEL(2),iEU(2)
        DO pe = iEL(3),iEU(3)
            
            IF ( Mask_PTR(RE,TE,PE,1) == iLeaf ) THEN
                iRE = FEM_Elem_Map(re,lvl)
                CALL Initialize_Slm_Table_on_Elem(  te, pe,             &
                                                    NQ(2), NQ(3),       &
                                                    iNE,                &
                                                    iEL,                &
                                                    Plm_Table,          &
                                                    Am_Table,           &
                                                    Slm_Elem_Table      )
  

                DO pd = 1,NQ(3)
                DO td = 1,NQ(2)
                DO rd = 1,NQ(1)

                    tpd = Map_To_tpd(td,pd)
                    LagP = Lagrange_Poly(RQ_Input(rd),DEGREE,FEM_Node_xlocs)
                    Tmp_U_Value = 0.0_idp

                    
                    DO d = 0,DEGREE
                        Current_Location = Map_To_FEM_Node(iRE,d)
                        Tmp_U_Value = Tmp_U_Value                                    &
                                    + SUM( dVA_Coeff_Vector(Current_Location,:,iU)  &
                                            * Slm_Elem_Values( :, tpd )            ) &
                                    * LagP(d)

                    END DO ! d Loop

                    Here = Quad_Map(rd, td, pd, NQ(1), NQ(2),NQ(3))
                    Result_PTR(re,te,pe,Here) = Tmp_U_Value
                END DO ! rd
                END DO ! td
                END DO ! pd

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
                                        MF_Results,             &
                                        FillGhostCells_Option   )


INTEGER,                                    INTENT(IN)      ::  iU
INTEGER,                                    INTENT(IN)      ::  iVB
INTEGER,    DIMENSION(3),                   INTENT(IN)      ::  NQ
REAL(idp),  DIMENSION(NQ(1)),               INTENT(IN)      ::  RQ_Input
REAL(idp),  DIMENSION(NQ(2)),               INTENT(IN)      ::  TQ_Input
REAL(idp),  DIMENSION(NQ(3)),               INTENT(IN)      ::  PQ_Input
REAL(idp),                                  INTENT(IN)      ::  Left_Limit
REAL(idp),                                  INTENT(IN)      ::  Right_Limit

INTEGER,                                    INTENT(IN)      ::  nLevels
TYPE(amrex_multifab),                       INTENT(INOUT)   ::  MF_Results(0:nLevels-1)

LOGICAL,                        OPTIONAL,   INTENT(IN)      ::  FillGhostCells_Option


INTEGER                                                 ::  Num_DOF
INTEGER                                                 ::  iRE
INTEGER                                                 ::  re, te, pe
INTEGER                                                 ::  rd, td, pd, tpd
INTEGER                                                 ::  d
INTEGER                                                 ::  lvl
INTEGER,    DIMENSION(3)                                ::  iNE

REAL(idp)                                               ::  Quad_Span
REAL(idp),  DIMENSION(0:DEGREE)                         ::  LagP
REAL(idp),  DIMENSION(1:NQ(1))                          ::  Cur_RX_Locs
REAL(idp)                                               ::  TMP_U_Value
INTEGER                                                 ::  Here, There

TYPE(amrex_mfiter)                                      ::  mfi
TYPE(amrex_box)                                         ::  Box
TYPE(amrex_imultifab)                                   ::  Level_Mask
INTEGER,    DIMENSION(3)                                ::  iEL, iEU
INTEGER                                                 ::  nComp

INTEGER,    CONTIGUOUS, POINTER                         ::  Mask_PTR(:,:,:,:)
REAL(idp),  CONTIGUOUS, POINTER                         ::  Result_PTR(:,:,:,:)

LOGICAL                                                 ::  FillGhostCells
INTEGER,    DIMENSION(1:3)                              ::  nGhost_Vec

REAL(idp),  DIMENSION(0:AMReX_Max_Grid_Size(2)-1)       ::  tlocs_subarray
REAL(idp),  DIMENSION(0:AMReX_Max_Grid_Size(3)-1)       ::  plocs_subarray

REAL(idp),  DIMENSION(1:NQ(2),1:LM_Short_Length,0:AMReX_Max_Grid_Size(2)-1) ::  Plm_Table
REAL(idp),  DIMENSION(1:NQ(3),1:LM_Length,0:AMReX_Max_Grid_Size(3)-1)       ::  Am_Table
REAL(idp),  DIMENSION(1:LM_Length, 1:NQ(2)*NQ(3) )                          ::  Slm_Elem_Table

tlocs_subarray = 0.0_idp
plocs_subarray = 0.0_idp

Quad_Span = Right_Limit - Left_Limit

Cur_RX_Locs = 2.0_idp * ( RQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp

Num_DOF = NQ(1)*NQ(2)*NQ(3)


IF ( PRESENT(FillGhostCells_Option) ) THEN
    FillGhostCells = FillGhostCells_Option
ELSE
    FillGhostCells = .FALSE.
END IF




DO lvl = nLevels-1,0,-1

    IF ( FillGhostCells ) THEN
        nGhost_Vec = MF_Results(lvl)%nghostvect()
    ELSE
        nGhost_Vec = 0
    END IF


    !
    !   MakeFineMask
    !
    IF ( lvl < AMReX_Num_Levels-1 ) THEN
        CALL AMReX_MakeFineMask(  Level_Mask,               &
                                  MF_Results(lvl)%ba,       &
                                  MF_Results(lvl)%dm,       &
                                  nGhost_Vec,               &
                                  MF_Results(lvl+1)%ba,     &
                                  iLeaf, iTrunk              )
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

        iEL = Box%lo
        iEU = Box%hi
        
        iNE = iEU-iEL+1
        
        DO te = iEL(2),iEU(2)
            tlocs_subarray(te-iEL(2)) = Level_dx(lvl,2)*te
        END DO
        DO pe = iEL(3),iEU(3)
            plocs_subarray(pe-iEL(3)) = Level_dx(lvl,3)*pe
        END DO
        
        
        ! Initialize Am Table
        CALL Initialize_Am_Table(   NQ(3),                      &
                                    PQ_Input,                   &
                                    L_Limit,                    &
                                    iNE(3),                     &
                                    [iEL(3), iEU(3)],           &
                                    plocs_subarray(0:iNE(3)-1), &
                                    Am_Table                    )

        ! Initialize Plm Table
        CALL Initialize_Plm_Table(  NQ(2),                      &
                                    TQ_Input,                   &
                                    L_Limit,                    &
                                    LM_Short_Length,            &
                                    iNE(2),                     &
                                    [iEL(2), iEU(2)],           &
                                    tlocs_subarray(0:iNE(2)-1), &
                                    Plm_Table                   )

        DO re = iEL(1),iEU(1)
        DO te = iEL(2),iEU(2)
        DO pe = iEL(3),iEU(3)


            IF ( Mask_PTR(RE,TE,PE,1) == iLeaf ) THEN



                iRE = FEM_Elem_Map(re,lvl)

                CALL Initialize_Slm_Table_on_Elem(  te, pe,             &
                                                    NQ(2), NQ(3),       &
                                                    iNE,                &
                                                    iEL,                &
                                                    Plm_Table,          &
                                                    Am_Table,           &
                                                    Slm_Elem_Table      )
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
                                    + SUM( dVB_Coeff_Vector(Here:There,iVB)        &
                                            * Slm_Elem_Table( :, tpd ) )   &
                                    * LagP(d)
                                    

                    END DO ! d Loop

                    Here = Quad_Map(rd, td, pd, NQ(1), NQ(2),NQ(3))
                    Result_PTR(re,te,pe,Here) = Tmp_U_Value

                END DO ! rd
                END DO ! td
                END DO ! pd

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
