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
            ONLY :  cVA_Coeff_Vector,       &
                    cVB_Coeff_Vector

USE Variables_Mesh, &
            ONLY :  rlocs,                  &
                    tlocs,                  &
                    plocs,                  &
                    iNE_Base


USE Variables_AMReX_Source, &
            ONLY :  iLeaf,                  &
                    iTrunk,                 &
                    iCovered,               &
                    iNotCovered,            &
                    iOutside,               &
                    iInterior

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
            ONLY :  Map_To_X_Space,         &
                    Map_From_X_Space

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
    
USE amrex_amrcore_module, &
            ONLY:   amrex_geom

USE amrex_multifab_module,  &
            ONLY :  amrex_multifab,         &
                    amrex_multifab_build,   &
                    amrex_imultifab_build,  &
                    amrex_imultifab_destroy

USE Variables_AMReX_Core, &
            ONLY :  AMReX_Num_Levels

USE Poseidon_AMReX_MakeFineMask_Module, &
            ONLY :  AMReX_MakeFineMask
            
USE Poseidon_AMReX_BuildMask_Module, &
            ONLY :  AMReX_BuildMask
            
USE Variables_AMReX_Source, &
            ONLY :  Source_PTR,         &
                    Mask_PTR,           &
                    Ghost_PTR

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

REAL(idp)                                                   ::  Quad_Span
REAL(idp),      DIMENSION(0:DEGREE)                         ::  LagP
REAL(idp),      DIMENSION(1:NQ(1))                          ::  Cur_RX_Locs
REAL(idp),      DIMENSION(1:NQ(2))                          ::  Cur_TX_Locs
REAL(idp),      DIMENSION(1:NQ(3))                          ::  Cur_PX_Locs
INTEGER                                                     ::  Current_Location

INTEGER,        DIMENSION(1:3)                              ::  nGhost_Vec


INTEGER                                                     ::  lvl
TYPE(amrex_mfiter)                                          ::  mfi
TYPE(amrex_box)                                             ::  Box
TYPE(amrex_imultifab)                                       ::  Level_Mask
INTEGER                                                     ::  nComp

TYPE(amrex_imultifab)                                       ::  Ghost_Mask

REAL(idp),  CONTIGUOUS, POINTER                             ::  Result_PTR(:,:,:,:)
REAL(idp),      DIMENSION(1:Local_Quad_DOF)                 ::  Var_Holder
REAL(idp),      DIMENSION(:,:), ALLOCATABLE                 ::  Translation_Matrix

INTEGER,        DIMENSION(1:3)                              ::  iE
INTEGER,        DIMENSION(1:3)                              ::  iEL, iEU
INTEGER,        DIMENSION(1:3)                              ::  iEL_A, iEU_A
LOGICAL                                                     ::  FillGhostCells



IF ( PRESENT(FillGhostCells_Option) ) THEN
    FillGhostCells = FillGhostCells_Option
ELSE
    FillGhostCells = .FALSE.
END IF


Quad_Span = Right_Limit - Left_Limit

Cur_RX_Locs = 2.0_idp * ( RQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp
Cur_TX_Locs = 2.0_idp * ( TQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp
Cur_PX_Locs = 2.0_idp * ( PQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp


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



    !
    !   BuildMask
    !
    CALL amrex_imultifab_build( Ghost_Mask,             &
                                MF_Results(lvl)%ba,     &
                                MF_Results(lvl)%dm,     &
                                1,                      &
                                nGhost_Vec(1)           )
                                
    CALL AMReX_BuildMask( Ghost_Mask,       &
                          amrex_geom(lvl),  &
                          iCovered,         &
                          iNotCovered,      &
                          iOutside,         &
                          iInterior         )


    CALL amrex_mfiter_build(mfi, MF_Results(lvl), tiling = .true. )

    DO WHILE(mfi%next())

        Result_PTR => MF_Results(lvl)%dataPtr(mfi)
        Mask_PTR   => Level_Mask%dataPtr(mfi)
        Ghost_PTR  => Ghost_Mask%dataPtr(mfi)

        Box = mfi%tilebox()
        nComp =  MF_Results(lvl)%ncomp()

        iEL_A = Box%lo
        iEU_A = Box%hi
        
        iEL = iEL_A-nGhost_Vec
        iEU = iEU_A+nGhost_Vec


        CALL Initialize_Normed_Legendre_Tables_on_Level( iEU, iEL, lvl )

        ! Fill Leaf Elements
        DO re = iEL(1),iEU(1)
        DO te = iEL(2),iEU(2)
        DO pe = iEL(3),iEU(3)
        
!            PRINT*,lvl,re,te,pe,Ghost_PTR(re,te,pe,1)
            
            IF ( Mask_PTR(RE,TE,PE,1) == iLeaf ) THEN

            iE = [re, te, pe]
            
            IF (     ( Ghost_PTR(re,te,pe,1) == iInterior )     &
                .OR. ( Ghost_PTR(re,te,pe,1) == iCovered  )     ) THEN
                
                CALL Poseidon_Valid_Type_A( iE, iEL, NQ,    &
                                            Cur_RX_Locs,    &
                                            Cur_TX_Locs,    &
                                            Cur_PX_Locs,    &
                                            lvl,            &
                                            iU,             &
                                            Var_Holder      )
                                            
                
                
            ELSE IF ( Ghost_PTR(re,te,pe,1) == iNotCovered) THEN
            
                CALL Poseidon_NotCovered_Type_A(  iE, iEL, NQ,    &
                                                  Cur_RX_Locs,    &
                                                  Cur_TX_Locs,    &
                                                  Cur_PX_Locs,    &
                                                  lvl,            &
                                                  iU,             &
                                                  Var_Holder      )
            
            
            ELSE IF ( Ghost_PTR(re,te,pe,1) == iOutside ) THEN
            
                CALL Poseidon_Outside_Type_A( iE, iEL, NQ,    &
                                              Cur_RX_Locs,    &
                                              Cur_TX_Locs,    &
                                              Cur_PX_Locs,    &
                                              lvl,            &
                                              iU,             &
                                              Var_Holder      )
        
            END IF !  Ghost_PTR
            END IF !  Mask_PTR(RE,TE,PE,1) == iLeaf

        END DO ! pe
        END DO ! te
        END DO ! re
        
    END DO

    CALL amrex_mfiter_destroy(mfi)
    CALL amrex_imultifab_destroy( Level_Mask )

END DO ! lvl


!STOP "In Poseidon_Return_AMReX_Type_A"

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

LOGICAL,                        OPTIONAL,   INTENT(IN)      ::  FillGhostCells_Option


INTEGER                                                 ::  Num_DOF
INTEGER                                                 ::  Output_Here
INTEGER                                                 ::  iRE
INTEGER                                                 ::  re, te, pe
INTEGER                                                 ::  rd, td, pd, tpd
INTEGER                                                 ::  d
INTEGER                                                 ::  lvl
INTEGER                                                 ::  Here

REAL(idp)                                               ::  Quad_Span
REAL(idp),      DIMENSION(0:DEGREE)                     ::  LagP
REAL(idp),      DIMENSION(1:NQ(1))                      ::  Cur_RX_Locs
REAL(idp),      DIMENSION(1:NQ(2))                      ::  Cur_TX_Locs
REAL(idp),      DIMENSION(1:NQ(3))                      ::  Cur_PX_Locs


TYPE(amrex_mfiter)                                      ::  mfi
TYPE(amrex_box)                                         ::  Box
TYPE(amrex_imultifab)                                   ::  Level_Mask
INTEGER,        DIMENSION(1:3)                          ::  iEL, iEU, iE
INTEGER                                                 ::  nComp

INTEGER,    CONTIGUOUS, POINTER                         ::  Mask_PTR(:,:,:,:)
REAL(idp),  CONTIGUOUS, POINTER                         ::  Result_PTR(:,:,:,:)
REAL(idp),      DIMENSION(1:Local_Quad_DOF)             ::  Var_Holder


COMPLEX(idp)                                            ::  TMP_U_Value

LOGICAL                                                 ::  FillGhostCells
INTEGER,        DIMENSION(1:3)                          ::  nGhost_Vec

Quad_Span = Right_Limit - Left_Limit

Cur_RX_Locs = 2.0_idp * ( RQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp
Cur_TX_Locs = 2.0_idp * ( TQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp
Cur_PX_Locs = 2.0_idp * ( PQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp

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

    CALL Initialize_Normed_Legendre_Tables_on_Level( iEU, iEL, lvl )

    CALL amrex_mfiter_build(mfi, MF_Results(lvl), tiling = .true. )

    DO WHILE(mfi%next())

        Result_PTR => MF_Results(lvl)%dataPtr(mfi)
        Mask_PTR   => Level_Mask%dataPtr(mfi)

        Box   = mfi%tilebox()
        nComp = MF_Results(lvl)%ncomp()

        iEL = Box%lo
        iEU = Box%hi

        DO re = iEL(1),iEU(1)
        DO te = iEL(2),iEU(2)
        DO pe = iEL(3),iEU(3)

            

            IF ( Mask_PTR(RE,TE,PE,1) == iLeaf ) THEN
            
            iE = [re,te,pe]
            
            IF (     ( Ghost_PTR(re,te,pe,1) == iInterior )     &
                .OR. ( Ghost_PTR(re,te,pe,1) == iCovered  )     ) THEN
                
                CALL Poseidon_Valid_Type_B( iE, iEL, NQ,    &
                                            Cur_RX_Locs,    &
                                            Cur_TX_Locs,    &
                                            Cur_PX_Locs,    &
                                            lvl,            &
                                            iU, iVB,        &
                                            Var_Holder      )
                
                
            ELSE IF ( Ghost_PTR(re,te,pe,1) == iNotCovered) THEN
            
                CALL Poseidon_NotCovered_Type_B( iE, iEL, NQ,    &
                                                 Cur_RX_Locs,    &
                                                 Cur_TX_Locs,    &
                                                 Cur_PX_Locs,    &
                                                 lvl,            &
                                                 iU, iVB,        &
                                                 Var_Holder      )
            
            
            ELSE IF ( Ghost_PTR(re,te,pe,1) == iOutside ) THEN
            
                CALL Poseidon_Outside_Type_B( iE, iEL, NQ,    &
                                              Cur_RX_Locs,    &
                                              Cur_TX_Locs,    &
                                              Cur_PX_Locs,    &
                                              lvl,            &
                                              iU, iVB,        &
                                              Var_Holder      )
        
            END IF !  Ghost_PTR
            END IF !  Mask_PTR(RE,TE,PE,1) == iLeaf

        
            DO Output_Here = 1,Num_DOF
                Here = (iU-1)*Num_DOF+Output_Here
                Result_PTR(re,te,pe,:) = Var_Holder(:)
            END DO ! Output_Here
        
        
        
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











 !+201+########################################################!
!                                                               !
!          Poseidon_Valid_Type_A                                !
!                                                               !
 !#############################################################!
SUBROUTINE Poseidon_Valid_Type_A( iE, iEL, NQ,  &
                                  Cur_RX_Locs,  &
                                  Cur_TX_Locs,  &
                                  Cur_PX_Locs,  &
                                  lvl,          &
                                  iU,           &
                                  Var_Holder    )

INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  iE
INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  NQ
INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  iEL
REAL(idp),  DIMENSION(1:NQ(1)),             INTENT(IN)      ::  Cur_RX_Locs
REAL(idp),  DIMENSION(1:NQ(2)),             INTENT(IN)      ::  Cur_TX_Locs
REAL(idp),  DIMENSION(1:NQ(3)),             INTENT(IN)      ::  Cur_PX_Locs
INTEGER,                                    INTENT(IN)      ::  lvl
INTEGER,                                    INTENT(IN)      ::  iU
REAL(idp),  DIMENSION(1:NQ(1)*NQ(2)*NQ(3)), INTENT(OUT)     ::  Var_Holder

INTEGER                                                     ::  iRE
INTEGER                                                     ::  pd, td, rd
INTEGER                                                     ::  tpd
INTEGER                                                     ::  d
INTEGER                                                     ::  Current_Location
INTEGER                                                     ::  Here
INTEGER                                                     ::  Output_Here

REAL(idp)                                                   ::  Tmp_U_Value

REAL(idp),  DIMENSION(0:DEGREE)                             ::  LagP

iRE = FEM_Elem_Map(iE(1),lvl)
CALL Initialize_Ylm_Tables_on_Elem( iE(2), iE(3), iEL, lvl )

DO pd = 1,NQ(3)
DO td = 1,NQ(2)
DO rd = 1,NQ(1)

    tpd = Map_To_tpd(td,pd)
    LagP = Lagrange_Poly(Cur_RX_Locs(rd),DEGREE,FEM_Node_xlocs)
    Tmp_U_Value = 0.0_idp


    DO d = 0,DEGREE
        Current_Location = Map_To_FEM_Node(iRE,d)
        Tmp_U_Value = Tmp_U_Value                                    &
                    + SUM( cVA_Coeff_Vector(Current_Location,:,iU)   &
                            * Ylm_Elem_Values( :, tpd )            ) &
                    * LagP(d)

    END DO ! d Loop

    Here = Quad_Map(rd, td, pd, NQ(1), NQ(2),NQ(3))
    Var_Holder(Here) = Tmp_U_Value
END DO ! rd
END DO ! td
END DO ! pd







END SUBROUTINE Poseidon_Valid_Type_A




 !+201+########################################################!
!                                                               !
!          Poseidon_Valid_Type_B                                !
!                                                               !
 !#############################################################!
SUBROUTINE Poseidon_Valid_Type_B( iE, iEL, NQ,    &
                                  Cur_RX_Locs,    &
                                  Cur_TX_Locs,  &
                                  Cur_PX_Locs,  &
                                  lvl,            &
                                  iU, iVB,        &
                                  Var_Holder      )

INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  iE
INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  NQ
INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  iEL
REAL(idp),  DIMENSION(1:NQ(1)),             INTENT(IN)      ::  Cur_RX_Locs
REAL(idp),  DIMENSION(1:NQ(2)),             INTENT(IN)      ::  Cur_TX_Locs
REAL(idp),  DIMENSION(1:NQ(3)),             INTENT(IN)      ::  Cur_PX_Locs
INTEGER,                                    INTENT(IN)      ::  lvl
INTEGER,                                    INTENT(IN)      ::  iU, iVB
REAL(idp),  DIMENSION(1:NQ(1)*NQ(2)*NQ(3)), INTENT(OUT)     ::  Var_Holder

INTEGER                                                     ::  iRE
INTEGER                                                     ::  pd, td, rd
INTEGER                                                     ::  tpd
INTEGER                                                     ::  d

INTEGER                                                     ::  Here, There

REAL(idp)                                                   ::  Tmp_U_Value

REAL(idp),  DIMENSION(0:DEGREE)                             ::  LagP

iRE = FEM_Elem_Map(iE(1),lvl)
CALL Initialize_Ylm_Tables_on_Elem( iE(2), iE(3), iEL, lvl )


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
    Var_Holder(Here) = Tmp_U_Value

END DO ! rd
END DO ! td
END DO ! pd



END SUBROUTINE Poseidon_Valid_Type_B







 !+201+########################################################!
!                                                               !
!          Poseidon_NotCovered_Type_A                           !
!                                                               !
 !#############################################################!
SUBROUTINE Poseidon_NotCovered_Type_A( iE, iEL, NQ,     &
                                       Cur_RX_Locs,     &
                                       Cur_TX_Locs,     &
                                       Cur_PX_Locs,     &
                                       lvl,             &
                                       iU,              &
                                       Var_Holder       )

INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  iE
INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  NQ
INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  iEL
REAL(idp),  DIMENSION(1:NQ(1)),             INTENT(IN)      ::  Cur_RX_Locs
REAL(idp),  DIMENSION(1:NQ(2)),             INTENT(IN)      ::  Cur_TX_Locs
REAL(idp),  DIMENSION(1:NQ(3)),             INTENT(IN)      ::  Cur_PX_Locs
INTEGER,                                    INTENT(IN)      ::  lvl
INTEGER,                                    INTENT(IN)      ::  iU
REAL(idp),  DIMENSION(1:NQ(1)*NQ(2)*NQ(3)), INTENT(OUT)     ::  Var_Holder


INTEGER,    DIMENSION(1:3)                                  ::  iE_Coarse
REAL(idp),  DIMENSION(0:Degree)                             ::  Coarse_xLocs
INTEGER                                                     ::  Num_Coarse_Locs
REAL(idp),  DIMENSION(:),   ALLOCATABLE                     ::  Coarse_Holder

INTEGER                                                     ::  pd, td, rd, tpd
INTEGER                                                     ::  d
INTEGER                                                     ::  iRE
INTEGER                                                     ::  Here
INTEGER                                                     ::  Current_Location

REAL(idp),  DIMENSION(0:DEGREE)                             ::  LagP
REAL(idp),  DIMENSION(0:DEGREE)                             ::  LagP_x
!REAL(idp),  DIMENSION(0:DEGREE)                             ::  LagP_y
!REAL(idp),  DIMENSION(0:DEGREE)                             ::  LagP_z

COMPLEX(idp)                                                ::  TMP_U_Value


Num_Coarse_Locs = (Degree+1)*NQ(2)*NQ(3)
ALLOCATE( Coarse_Holder(1:Num_Coarse_Locs) )


iE_Coarse = iE/2
Coarse_xLocs = Map_From_X_Space( -1.0_idp, 0.0_idp, FEM_Node_xlocs)


iRE = FEM_Elem_Map(iE_Coarse(1),lvl)
CALL Initialize_Ylm_Tables_on_Elem( iE(2), iE(3), iEL, lvl )


DO pd = 1,NQ(3)
DO td = 1,NQ(2)
DO rd = 1,Degree+1

    tpd = Map_To_tpd(td,pd)
    LagP = Lagrange_Poly(Coarse_xLocs(rd-1),DEGREE,FEM_Node_xlocs)
    Tmp_U_Value = 0.0_idp

    DO d = 0,DEGREE
    
        Current_Location = Map_To_FEM_Node(iRE,d)
        Tmp_U_Value = Tmp_U_Value                                    &
                    + SUM( cVA_Coeff_Vector(Current_Location,:,iU)  &
                            * Ylm_Elem_Values( :, tpd )            ) &
                    * LagP(d)

    END DO ! d Loop

    Here = Quad_Map(rd, td, pd, Degree+1, NQ(2),NQ(3))
    Coarse_Holder(Here) = Tmp_U_Value

END DO ! rd
END DO ! td
END DO ! pd



DO pd = 1,NQ(3)
DO td = 1,NQ(2)
DO rd = 1,NQ(1)

    tpd = Map_To_tpd(td,pd)
    LagP_x = Lagrange_Poly(Cur_RX_Locs(rd),DEGREE,FEM_Node_xlocs)
!    LagP_y = Lagrange_Poly(Cur_TX_Locs(td),DEGREE,FEM_Node_xlocs)
!    LagP_z = Lagrange_Poly(Cur_PX_Locs(pd),DEGREE,FEM_Node_xlocs)
    Tmp_U_Value = 0.0_idp

    DO d = 0,DEGREE
    
        Here = Quad_Map(rd, td, pd, Degree+1, NQ(2),NQ(3))
        
        Tmp_U_Value = Tmp_U_Value           &
                    + Coarse_Holder(Here)   &
                    * LagP_x(d)

    END DO ! d Loop

    Here = Quad_Map(rd, td, pd, NQ(1), NQ(2),NQ(3))
    Var_Holder(Here) = Tmp_U_Value

END DO ! rd
END DO ! td
END DO ! pd



END SUBROUTINE Poseidon_NotCovered_Type_A



 !+201+########################################################!
!                                                               !
!          Poseidon_NotCovered_Type_B                           !
!                                                               !
 !#############################################################!
SUBROUTINE Poseidon_NotCovered_Type_B( iE, iEL, NQ,     &
                                       Cur_RX_Locs,     &
                                       Cur_TX_Locs,     &
                                       Cur_PX_Locs,     &
                                       lvl,             &
                                       iU, iVB,         &
                                       Var_Holder       )

INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  iE
INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  NQ
INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  iEL
REAL(idp),  DIMENSION(1:NQ(1)),             INTENT(IN)      ::  Cur_RX_Locs
REAL(idp),  DIMENSION(1:NQ(2)),             INTENT(IN)      ::  Cur_TX_Locs
REAL(idp),  DIMENSION(1:NQ(3)),             INTENT(IN)      ::  Cur_PX_Locs
INTEGER,                                    INTENT(IN)      ::  lvl
INTEGER,                                    INTENT(IN)      ::  iU, iVB
REAL(idp),  DIMENSION(1:NQ(1)*NQ(2)*NQ(3)), INTENT(OUT)     ::  Var_Holder

INTEGER,    DIMENSION(1:3)                                  ::  iE_Coarse
REAL(idp),  DIMENSION(0:Degree)                             ::  Coarse_xLocs
INTEGER                                                     ::  Num_Coarse_Locs
REAL(idp),  DIMENSION(:),   ALLOCATABLE                     ::  Coarse_Holder

INTEGER                                                     ::  iRE
INTEGER                                                     ::  pd, td, rd, tpd
INTEGER                                                     ::  d
INTEGER                                                     ::  Here, There

REAL(idp),  DIMENSION(0:DEGREE)                             ::  LagP
REAL(idp),  DIMENSION(0:DEGREE)                             ::  LagP_x

COMPLEX(idp)                                                ::  TMP_U_Value

Num_Coarse_Locs = (Degree+1)*NQ(2)*NQ(3)
ALLOCATE( Coarse_Holder(1:Num_Coarse_Locs) )


iE_Coarse = iE/2
Coarse_xLocs = Map_From_X_Space( -1.0_idp, 0.0_idp, FEM_Node_xlocs)


iRE = FEM_Elem_Map(iE_Coarse(1),lvl)
CALL Initialize_Ylm_Tables_on_Elem( iE(2), iE(3), iEL, lvl )


DO pd = 1,NQ(3)
DO td = 1,NQ(2)
DO rd = 1,Degree+1

    tpd = Map_To_tpd(td,pd)
    LagP = Lagrange_Poly(Coarse_xLocs(rd-1),DEGREE,FEM_Node_xlocs)
    Tmp_U_Value = 0.0_idp
    
    DO d = 0,DEGREE
    
        Here  = FP_Array_Map_TypeB(iU,iVB,iRE,d,1)
        There = FP_Array_Map_TypeB(iU,iVB,iRE,d,LM_Length)
        
        Tmp_U_Value = Tmp_U_Value                                   &
                    + SUM( cVB_Coeff_Vector(Here:There,iVB)        &
                            * Ylm_Elem_Values( :, tpd ) )   &
                    * LagP(d)

    END DO ! d Loop

    Here = Quad_Map(rd, td, pd, Degree+1, NQ(2),NQ(3))
    Coarse_Holder(Here) = Tmp_U_Value

END DO ! rd
END DO ! td
END DO ! pd



DO pd = 1,NQ(3)
DO td = 1,NQ(2)
DO rd = 1,NQ(1)

    tpd = Map_To_tpd(td,pd)
    LagP_x = Lagrange_Poly(Cur_RX_Locs(rd),DEGREE,FEM_Node_xlocs)
!    LagP_y = Lagrange_Poly(Cur_TX_Locs(td),DEGREE,FEM_Node_xlocs)
!    LagP_z = Lagrange_Poly(Cur_PX_Locs(pd),DEGREE,FEM_Node_xlocs)
    Tmp_U_Value = 0.0_idp

    DO d = 0,DEGREE
    
        Here = Quad_Map(rd, td, pd, Degree+1, NQ(2),NQ(3))
        
        Tmp_U_Value = Tmp_U_Value           &
                    + Coarse_Holder(Here)   &
                    * LagP_x(d)

    END DO ! d Loop

    Here = Quad_Map(rd, td, pd, NQ(1), NQ(2),NQ(3))
    Var_Holder(Here) = Tmp_U_Value

END DO ! rd
END DO ! td
END DO ! pd


END SUBROUTINE Poseidon_NotCovered_Type_B














 !+201+########################################################!
!                                                               !
!          Poseidon_Outside_Type_A                              !
!                                                               !
 !#############################################################!
SUBROUTINE Poseidon_Outside_Type_A( iE, iEL, NQ,    &
                                    Cur_RX_Locs,    &
                                    Cur_TX_Locs,    &
                                    Cur_PX_Locs,    &
                                    lvl,            &
                                    iU,             &
                                    Var_Holder      )

INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  iE
INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  NQ
INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  iEL
REAL(idp),  DIMENSION(1:NQ(1)),             INTENT(IN)      ::  Cur_RX_Locs
REAL(idp),  DIMENSION(1:NQ(2)),             INTENT(IN)      ::  Cur_TX_Locs
REAL(idp),  DIMENSION(1:NQ(3)),             INTENT(IN)      ::  Cur_PX_Locs
INTEGER,                                    INTENT(IN)      ::  lvl
INTEGER,                                    INTENT(IN)      ::  iU
REAL(idp),  DIMENSION(1:NQ(1)*NQ(2)*NQ(3)), INTENT(OUT)     ::  Var_Holder

Var_Holder = 0.0_idp

END SUBROUTINE Poseidon_Outside_Type_A



 !+201+########################################################!
!                                                               !
!          Poseidon_Outside_Type_B                              !
!                                                               !
 !#############################################################!
SUBROUTINE Poseidon_Outside_Type_B( iE, iEL, NQ,    &
                                    Cur_RX_Locs,    &
                                    Cur_TX_Locs,    &
                                    Cur_PX_Locs,    &
                                    lvl,            &
                                    iU, iVB,        &
                                    Var_Holder      )

INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  iE
INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  NQ
INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  iEL
REAL(idp),  DIMENSION(1:NQ(1)),             INTENT(IN)      ::  Cur_RX_Locs
REAL(idp),  DIMENSION(1:NQ(2)),             INTENT(IN)      ::  Cur_TX_Locs
REAL(idp),  DIMENSION(1:NQ(3)),             INTENT(IN)      ::  Cur_PX_Locs
INTEGER,                                    INTENT(IN)      ::  lvl
INTEGER,                                    INTENT(IN)      ::  iU, iVB
REAL(idp),  DIMENSION(1:NQ(1)*NQ(2)*NQ(3)), INTENT(OUT)     ::  Var_Holder

Var_Holder = 0.0_idp

END SUBROUTINE Poseidon_Outside_Type_B



END MODULE Return_Functions_AMReX
