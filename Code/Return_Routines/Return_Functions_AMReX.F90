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
            ONLY :  idp

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
            ONLY :  Level_DX

USE Variables_Derived, &
            ONLY :  LM_Length,                  &
                    LM_Short_Length

USE Variables_Vectors, &
            ONLY :  dVA_Coeff_Vector,      &
                    dVB_Coeff_Vector


USE Variables_Mesh, &
            ONLY :  rlocs,                  &
                    tlocs,                  &
                    plocs,                  &
                    iNE_Base

USE Variables_FEM_Module, &
            ONLY :  FEM_Node_xlocs

USE Maps_Fixed_Point, &
            ONLY :  FP_Array_Map_TypeB

USE Maps_Domain, &
            ONLY :  Map_To_FEM_Node,        &
                    FEM_Elem_Map

USE Maps_X_Space, &
            ONLY :  Map_From_X_Space

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

USE Maps_Quadrature, &
            ONLY :  Map_To_tpd,             &
                    Quad_Map

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
    
USE amrex_amrcore_module, &
            ONLY:   amrex_geom,             &
                    amrex_get_numlevels

USE amrex_multifab_module,  &
            ONLY :  amrex_multifab,         &
                    amrex_multifab_build,   &
                    amrex_imultifab_build,  &
                    amrex_imultifab_destroy

USE Variables_AMReX_Core, &
            ONLY :  AMReX_Num_Levels,       &
                    AMReX_Max_Grid_Size,    &
                    Source_PTR,             &
                    Mask_PTR,               &
                    Ghost_PTR

USE Parameters_AMReX, &
            ONLY :  iLeaf,                  &
                    iTrunk,                 &
                    iCovered,               &
                    iNotCovered,            &
                    iOutside,               &
                    iInterior

USE Poseidon_AMReX_MakeFineMask_Module, &
            ONLY :  AMReX_MakeFineMask
            
USE Poseidon_AMReX_BuildMask_Module, &
            ONLY :  AMReX_BuildMask


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
REAL(idp),      DIMENSION(1:NQ(2))                          ::  Cur_TX_Locs
REAL(idp),      DIMENSION(1:NQ(3))                          ::  Cur_PX_Locs

INTEGER                                                     ::  Current_Location

INTEGER,        DIMENSION(1:3)                              ::  nGhost_Vec

INTEGER                                                     ::  Output_Here

INTEGER                                                     ::  lvl
TYPE(amrex_mfiter)                                          ::  mfi
TYPE(amrex_box)                                             ::  Box
TYPE(amrex_imultifab)                                       ::  Level_Mask
INTEGER                                                     ::  nComp


INTEGER,    CONTIGUOUS, POINTER                             ::  Mask_PTR(:,:,:,:)
REAL(idp),  CONTIGUOUS, POINTER                             ::  Result_PTR(:,:,:,:)
TYPE(amrex_imultifab)                                       ::  Ghost_Mask

REAL(idp),      DIMENSION(:),   ALLOCATABLE                 ::  Var_Holder

INTEGER,        DIMENSION(1:3)                              ::  iE
INTEGER,        DIMENSION(1:3)                              ::  iEL, iEU
INTEGER,        DIMENSION(1:3)                              ::  iEL_A, iEU_A
LOGICAL                                                     ::  FillGhostCells

REAL(idp),  DIMENSION(0:AMReX_Max_Grid_Size(2))             ::  tlocs_subarray
REAL(idp),  DIMENSION(0:AMReX_Max_Grid_Size(3))             ::  plocs_subarray
    
REAL(idp),  DIMENSION(1:NQ(2),1:LM_Short_Length,0:AMReX_Max_Grid_Size(2)-1) ::  Plm_Table
REAL(idp),  DIMENSION(1:NQ(3),1:LM_Length,0:AMReX_Max_Grid_Size(3)-1)       ::  Am_Table


INTEGER     :: mLevels


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
ALLOCATE( Var_Holder(Num_DOF) )


mLevels = amrex_get_numlevels()

DO lvl = mLevels-1,0,-1

    IF ( FillGhostCells ) THEN
        nGhost_Vec = MF_Results(lvl)%nghostvect()
    ELSE
        nGhost_Vec = 0
    END IF
    
    !
    !   MakeFineMask
    !
    IF ( lvl < mLevels-1 ) THEN
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

        IF ( ANY( iEL < 0 ) ) THEN
            ! Reflecting Conditions
            iEL = iEL_A
        END IF

        IF ( ANY( iEU .GE. (2**lvl)*iNE_Base(1) ) ) THEN
            iEU = iEU_A
        END IF
        
        iNE = iEU-iEL+1
        
        
        tlocs_subarray = 0.0_idp
        plocs_subarray = 0.0_idp
        DO te = iEL(2),iEU(2)+1
            tlocs_subarray(te-iEL(2)) = Level_dx(lvl,2)*te
        END DO
        DO pe = iEL(3),iEU(3)+1
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
!        Plm_Table(:,1,:) = 1.0_idp
        CALL Initialize_Plm_Table(  NQ(2),                      &
                                    TQ_Input,                   &
                                    L_Limit,                    &
                                    LM_Short_Length,            &
                                    iNE(2),                     &
                                    [iEL(2), iEU(2)],           &
                                    tlocs_subarray(0:iNE(2)),   &
                                    Plm_Table                   )


        ! Fill Leaf Elements
        DO re = iEL(1),iEU(1)
        DO te = iEL(2),iEU(2)
        DO pe = iEL(3),iEU(3)
        
            IF ( Mask_PTR(RE,TE,PE,1) == iLeaf ) THEN

            iE = [re, te, pe]
            
            IF (     ( Ghost_PTR(re,te,pe,1) == iInterior )     &
                .OR. ( Ghost_PTR(re,te,pe,1) == iCovered  )     ) THEN
                
                CALL Poseidon_Valid_Type_A( iE, iEL, iNE, NQ,    &
                                            Cur_RX_Locs,    &
                                            Cur_TX_Locs,    &
                                            Cur_PX_Locs,    &
                                            lvl,            &
                                            iU,             &
                                            Am_Table,       &
                                            Plm_Table,      &
                                            Var_Holder      )
                                        
            ELSE IF ( Ghost_PTR(re,te,pe,1) == iNotCovered) THEN
                CALL Poseidon_NotCovered_Type_A(  iE, iEL, iNE, NQ,    &
                                                  Cur_RX_Locs,    &
                                                  Cur_TX_Locs,    &
                                                  Cur_PX_Locs,    &
                                                  lvl,            &
                                                  iU,             &
                                                  Am_Table,       &
                                                  Plm_Table,      &
                                                  Var_Holder      )
            
            
            ELSE IF ( Ghost_PTR(re,te,pe,1) == iOutside ) THEN
                CALL Poseidon_Outside_Type_A( iE, iEL, iNE, NQ,    &
                                              Cur_RX_Locs,    &
                                              Cur_TX_Locs,    &
                                              Cur_PX_Locs,    &
                                              lvl,            &
                                              iU,             &
                                              Am_Table,       &
                                              Plm_Table,      &
                                              Var_Holder      )
        
            END IF !  Ghost_PTR
            
            
            DO Output_Here = 1,Num_DOF
                Here = (iU-1)*Num_DOF+Output_Here
                Result_PTR(re,te,pe,Here) = Var_Holder(Output_Here)
            END DO ! Output_Here

            END IF !  Mask_PTR(RE,TE,PE,1) == iLeaf

        END DO ! pe
        END DO ! te
        END DO ! re
        
    END DO

    CALL amrex_mfiter_destroy(mfi)
    CALL amrex_imultifab_destroy( Level_Mask )

END DO ! lvl

DEALLOCATE( Var_Holder )

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
INTEGER,    DIMENSION(3)                                ::  iE

INTEGER                                                 ::  Output_Here
REAL(idp),      DIMENSION(:), ALLOCATABLE               ::  Var_Holder


REAL(idp)                                               ::  Quad_Span
REAL(idp),  DIMENSION(0:DEGREE)                         ::  LagP
REAL(idp),  DIMENSION(1:NQ(1))                          ::  Cur_RX_Locs
REAL(idp),  DIMENSION(1:NQ(2))                          ::  Cur_TX_Locs
REAL(idp),  DIMENSION(1:NQ(3))                          ::  Cur_PX_Locs
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

tlocs_subarray = 0.0_idp
plocs_subarray = 0.0_idp


Quad_Span = Right_Limit - Left_Limit

Cur_RX_Locs = 2.0_idp * ( RQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp
Cur_TX_Locs = 2.0_idp * ( TQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp
Cur_PX_Locs = 2.0_idp * ( PQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp

Num_DOF = NQ(1)*NQ(2)*NQ(3)

ALLOCATE( Var_Holder(1:Num_DOF) )

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

        Box   = mfi%tilebox()
        nComp = MF_Results(lvl)%ncomp()

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
                                    Cur_PX_Locs,                &
                                    L_Limit,                    &
                                    iNE(3),                     &
                                    [iEL(3), iEU(3)],           &
                                    plocs_subarray(0:iNE(3)-1), &
                                    Am_Table                    )

        ! Initialize Plm Table
        CALL Initialize_Plm_Table(  NQ(2),                      &
                                    Cur_TX_Locs,                &
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
            
            iE = [re,te,pe]
            
            IF (     ( Ghost_PTR(re,te,pe,1) == iInterior )     &
                .OR. ( Ghost_PTR(re,te,pe,1) == iCovered  )     ) THEN
                
                
                CALL Poseidon_Valid_Type_B( iE, iEL, iNE, NQ,    &
                                            Cur_RX_Locs,    &
                                            Cur_TX_Locs,    &
                                            Cur_PX_Locs,    &
                                            lvl,            &
                                            iU, iVB,        &
                                            Am_Table,       &
                                            Plm_Table,      &
                                            Var_Holder      )
                
                
            ELSE IF ( Ghost_PTR(re,te,pe,1) == iNotCovered) THEN
            
                CALL Poseidon_NotCovered_Type_B( iE, iEL, iNE, NQ,    &
                                                 Cur_RX_Locs,    &
                                                 Cur_TX_Locs,    &
                                                 Cur_PX_Locs,    &
                                                 lvl,            &
                                                 iU, iVB,        &
                                                Am_Table,       &
                                                Plm_Table,      &
                                                 Var_Holder      )
            
            
            ELSE IF ( Ghost_PTR(re,te,pe,1) == iOutside ) THEN
            
                CALL Poseidon_Outside_Type_B( iE, iEL, iNE, NQ,    &
                                              Cur_RX_Locs,    &
                                              Cur_TX_Locs,    &
                                              Cur_PX_Locs,    &
                                              lvl,            &
                                              iU, iVB,        &
                                              Am_Table,       &
                                              Plm_Table,      &
                                              Var_Holder      )
        
            END IF !  Ghost_PTR
            
            
            DO Output_Here = 1,Num_DOF
                Here = (iU-1)*Num_DOF+Output_Here
                Result_PTR(re,te,pe,Here) = Var_Holder(Output_Here)
                

            END DO ! Output_Here
            
            END IF !  Mask_PTR(RE,TE,PE,1) == iLeaf

        END DO ! pe
        END DO ! te
        END DO ! re

    END DO

    CALL amrex_mfiter_destroy(mfi)
    CALL amrex_imultifab_destroy( Level_Mask )

END DO ! lvl



DEALLOCATE(Var_Holder)

END SUBROUTINE Poseidon_Return_AMReX_Type_B








 !+201+########################################################!
!                                                               !
!          Poseidon_Valid_Type_A                                !
!                                                               !
 !#############################################################!
SUBROUTINE Poseidon_Valid_Type_A( iE, iEL, iNE, NQ,  &
                                  Cur_RX_Locs,  &
                                  Cur_TX_Locs,  &
                                  Cur_PX_Locs,  &
                                  lvl,          &
                                  iU,           &
                                  Am_Table,       &
                                  Plm_Table,      &
                                  Var_Holder    )

INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  iE
INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  NQ
INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  iEL
INTEGER,    DIMENSION(3),                   INTENT(IN)      ::  iNE
REAL(idp),  DIMENSION(1:NQ(1)),             INTENT(IN)      ::  Cur_RX_Locs
REAL(idp),  DIMENSION(1:NQ(2)),             INTENT(IN)      ::  Cur_TX_Locs
REAL(idp),  DIMENSION(1:NQ(3)),             INTENT(IN)      ::  Cur_PX_Locs
INTEGER,                                    INTENT(IN)      ::  lvl
INTEGER,                                    INTENT(IN)      ::  iU

REAL(idp),  DIMENSION(1:NQ(2),                              &
                      1:LM_Short_Length,                    &
                      0:AMReX_Max_Grid_Size(2)-1),          &
                                            INTENT(IN)      ::  Plm_Table
                                            
REAL(idp),  DIMENSION(1:NQ(3),                              &
                      1:LM_Length,                          &
                      0:AMReX_Max_Grid_Size(3)-1)           ::  Am_Table
                      
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
REAL(idp),  DIMENSION(1:LM_Length, 1:NQ(2)*NQ(3) )          ::  Slm_Elem_Table



iRE = FEM_Elem_Map(iE(1),lvl)
CALL Initialize_Slm_Table_on_Elem(  iE(2), iE(3),       &
                                    NQ(2), NQ(3),       &
                                    iNE,                &
                                    iEL,                &
                                    Plm_Table,          &
                                    Am_Table,           &
                                    Slm_Elem_Table      )

Var_Holder = 0.0_idp
DO pd = 1,NQ(3)
DO td = 1,NQ(2)
DO rd = 1,NQ(1)

    tpd = Map_To_tpd(td,pd,NQ(3))
    LagP = Lagrange_Poly(Cur_RX_Locs(rd),DEGREE,FEM_Node_xlocs)
    Tmp_U_Value = 0.0_idp

    
    DO d = 0,DEGREE
        Current_Location = Map_To_FEM_Node(iRE,d)
        Tmp_U_Value = Tmp_U_Value                                    &
                    + SUM( dVA_Coeff_Vector(Current_Location,:,iU)  &
                            * Slm_Elem_Table( :, tpd )            ) &
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
SUBROUTINE Poseidon_Valid_Type_B( iE, iEL, iNE, NQ,    &
                                  Cur_RX_Locs,    &
                                  Cur_TX_Locs,  &
                                  Cur_PX_Locs,  &
                                  lvl,            &
                                  iU, iVB,        &
                                  Am_Table,       &
                                  Plm_Table,      &
                                  Var_Holder      )

INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  iE
INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  NQ
INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  iEL
INTEGER,    DIMENSION(3),                   INTENT(IN)      ::  iNE
REAL(idp),  DIMENSION(1:NQ(1)),             INTENT(IN)      ::  Cur_RX_Locs
REAL(idp),  DIMENSION(1:NQ(2)),             INTENT(IN)      ::  Cur_TX_Locs
REAL(idp),  DIMENSION(1:NQ(3)),             INTENT(IN)      ::  Cur_PX_Locs
INTEGER,                                    INTENT(IN)      ::  lvl
INTEGER,                                    INTENT(IN)      ::  iU, iVB
REAL(idp),  DIMENSION(1:NQ(2),                              &
                      1:LM_Short_Length,                    &
                      0:AMReX_Max_Grid_Size(2)-1),          &
                                            INTENT(IN)      ::  Plm_Table
                                            
REAL(idp),  DIMENSION(1:NQ(3),                              &
                      1:LM_Length,                          &
                      0:AMReX_Max_Grid_Size(3)-1)           ::  Am_Table
REAL(idp),  DIMENSION(1:NQ(1)*NQ(2)*NQ(3)), INTENT(OUT)     ::  Var_Holder

INTEGER                                                     ::  iRE
INTEGER                                                     ::  pd, td, rd
INTEGER                                                     ::  tpd
INTEGER                                                     ::  d

INTEGER                                                     ::  Here, There

REAL(idp)                                                   ::  Tmp_U_Value

REAL(idp),  DIMENSION(0:DEGREE)                             ::  LagP
REAL(idp),  DIMENSION(1:LM_Length, 1:NQ(2)*NQ(3) )          ::  Slm_Elem_Table

iRE = FEM_Elem_Map(iE(1),lvl)

CALL Initialize_Slm_Table_on_Elem(  iE(2), iE(3),       &
                                    NQ(2), NQ(3),       &
                                    iNE,                &
                                    iEL,                &
                                    Plm_Table,          &
                                    Am_Table,           &
                                    Slm_Elem_Table      )
DO pd = 1,NQ(3)
DO td = 1,NQ(2)
DO rd = 1,NQ(1)

    tpd = Map_To_tpd(td,pd,NQ(3))
    LagP = Lagrange_Poly(Cur_RX_Locs(rd),DEGREE,FEM_Node_xlocs)
    Tmp_U_Value = 0.0_idp

    DO d = 0,DEGREE

        Here  = FP_Array_Map_TypeB(iU,iVB,iRE,d,1)
        There = FP_Array_Map_TypeB(iU,iVB,iRE,d,LM_Length)
        
        Tmp_U_Value = Tmp_U_Value                               &
                    + SUM( dVB_Coeff_Vector(Here:There,iVB)     &
                            * Slm_Elem_Table( :, tpd ) )        &
                    * LagP(d)
                    

    END DO ! d Loop

    Here = Quad_Map(rd, td, pd, NQ(1), NQ(2), NQ(3))
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
SUBROUTINE Poseidon_NotCovered_Type_A( iE, iEL, iNE, NQ,     &
                                       Cur_RX_Locs,     &
                                       Cur_TX_Locs,     &
                                       Cur_PX_Locs,     &
                                       lvl,             &
                                       iU,              &
                                       Am_Table,        &
                                       Plm_Table,       &
                                       Var_Holder       )

INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  iE
INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  NQ
INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  iEL
INTEGER,    DIMENSION(3),                   INTENT(IN)      ::  iNE
REAL(idp),  DIMENSION(1:NQ(1)),             INTENT(IN)      ::  Cur_RX_Locs
REAL(idp),  DIMENSION(1:NQ(2)),             INTENT(IN)      ::  Cur_TX_Locs
REAL(idp),  DIMENSION(1:NQ(3)),             INTENT(IN)      ::  Cur_PX_Locs
INTEGER,                                    INTENT(IN)      ::  lvl
INTEGER,                                    INTENT(IN)      ::  iU
REAL(idp),  DIMENSION(1:NQ(2),                              &
                      1:LM_Short_Length,                    &
                      0:AMReX_Max_Grid_Size(2)-1),          &
                                            INTENT(IN)      ::  Plm_Table
                                            
REAL(idp),  DIMENSION(1:NQ(3),                              &
                      1:LM_Length,                          &
                      0:AMReX_Max_Grid_Size(3)-1)           ::  Am_Table
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

REAL(idp)                                                   ::  TMP_U_Value
REAL(idp),  DIMENSION(1:LM_Length, 1:NQ(2)*NQ(3) )          ::  Slm_Elem_Table


Num_Coarse_Locs = (Degree+1)*NQ(2)*NQ(3)
ALLOCATE( Coarse_Holder(1:Num_Coarse_Locs) )

iE_Coarse = iE/2
Coarse_xLocs = Map_From_X_Space( -1.0_idp, 0.0_idp, FEM_Node_xlocs)



iRE = FEM_Elem_Map(iE_Coarse(1),lvl-1)
CALL Initialize_Slm_Table_on_Elem(  iE(2), iE(3),       &
                                    NQ(2), NQ(3),       &
                                    iNE,                &
                                    iEL,                &
                                    Plm_Table,          &
                                    Am_Table,           &
                                    Slm_Elem_Table      )


DO pd = 1,NQ(3)
DO td = 1,NQ(2)
DO rd = 1,Degree+1

    tpd = Map_To_tpd(td,pd,NQ(3))
    LagP = Lagrange_Poly(Coarse_xLocs(rd-1),DEGREE,FEM_Node_xlocs)
    Tmp_U_Value = 0.0_idp

    DO d = 0,DEGREE
    
        Current_Location = Map_To_FEM_Node(iRE,d)

        Tmp_U_Value = Tmp_U_Value                                    &
                    + SUM( dVA_Coeff_Vector(Current_Location,:,iU)  &
                            * Slm_Elem_Table( :, tpd )            ) &
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

    tpd = Map_To_tpd(td,pd,NQ(3))
    LagP_x = Lagrange_Poly(Cur_RX_Locs(rd),DEGREE,FEM_Node_xlocs)
!    LagP_y = Lagrange_Poly(Cur_TX_Locs(td),DEGREE,FEM_Node_xlocs)
!    LagP_z = Lagrange_Poly(Cur_PX_Locs(pd),DEGREE,FEM_Node_xlocs)
    Tmp_U_Value = 0.0_idp

    DO d = 0,DEGREE
    
        Here = Quad_Map(d+1, td, pd, Degree+1, NQ(2),NQ(3))
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
SUBROUTINE Poseidon_NotCovered_Type_B( iE, iEL, iNE, NQ,     &
                                       Cur_RX_Locs,     &
                                       Cur_TX_Locs,     &
                                       Cur_PX_Locs,     &
                                       lvl,             &
                                       iU, iVB,         &
                                       Am_Table,        &
                                       Plm_Table,       &
                                       Var_Holder       )

INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  iE
INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  NQ
INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  iEL
INTEGER,    DIMENSION(3),                   INTENT(IN)      ::  iNE
REAL(idp),  DIMENSION(1:NQ(1)),             INTENT(IN)      ::  Cur_RX_Locs
REAL(idp),  DIMENSION(1:NQ(2)),             INTENT(IN)      ::  Cur_TX_Locs
REAL(idp),  DIMENSION(1:NQ(3)),             INTENT(IN)      ::  Cur_PX_Locs
INTEGER,                                    INTENT(IN)      ::  lvl
INTEGER,                                    INTENT(IN)      ::  iU, iVB
REAL(idp),  DIMENSION(1:NQ(2),                              &
                      1:LM_Short_Length,                    &
                      0:AMReX_Max_Grid_Size(2)-1),          &
                                            INTENT(IN)      ::  Plm_Table
                                            
REAL(idp),  DIMENSION(1:NQ(3),                              &
                      1:LM_Length,                          &
                      0:AMReX_Max_Grid_Size(3)-1)           ::  Am_Table
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

REAL(idp)                                                   ::  TMP_U_Value
REAL(idp),  DIMENSION(1:LM_Length, 1:NQ(2)*NQ(3) )          ::  Slm_Elem_Table

Num_Coarse_Locs = (Degree+1)*NQ(2)*NQ(3)
ALLOCATE( Coarse_Holder(1:Num_Coarse_Locs) )


iE_Coarse = iE/2
Coarse_xLocs = Map_From_X_Space( -1.0_idp, 0.0_idp, FEM_Node_xlocs)


iRE = FEM_Elem_Map(iE_Coarse(1),lvl-1)

CALL Initialize_Slm_Table_on_Elem(  iE(2), iE(3),       &
                                    NQ(2), NQ(3),       &
                                    iNE,                &
                                    iEL,                &
                                    Plm_Table,          &
                                    Am_Table,           &
                                    Slm_Elem_Table      )

DO pd = 1,NQ(3)
DO td = 1,NQ(2)
DO rd = 1,Degree+1

    tpd = Map_To_tpd(td,pd,NQ(3))
    LagP = Lagrange_Poly(Coarse_xLocs(rd-1),DEGREE,FEM_Node_xlocs)
    Tmp_U_Value = 0.0_idp
    
    DO d = 0,DEGREE
    
        Here  = FP_Array_Map_TypeB(iU,iVB,iRE,d,1)
        There = FP_Array_Map_TypeB(iU,iVB,iRE,d,LM_Length)
        
        Tmp_U_Value = Tmp_U_Value                                   &
                    + SUM( dVB_Coeff_Vector(Here:There,iVB)        &
                            * Slm_Elem_Table( :, tpd ) )   &
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

    tpd = Map_To_tpd(td,pd,NQ(3))
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
SUBROUTINE Poseidon_Outside_Type_A( iE, iEL, iNE, NQ,    &
                                    Cur_RX_Locs,    &
                                    Cur_TX_Locs,    &
                                    Cur_PX_Locs,    &
                                    lvl,            &
                                    iU,             &
                                    Am_Table,       &
                                    Plm_Table,      &
                                    Var_Holder      )

INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  iE
INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  NQ
INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  iEL
INTEGER,    DIMENSION(3),                   INTENT(IN)      ::  iNE
REAL(idp),  DIMENSION(1:NQ(1)),             INTENT(IN)      ::  Cur_RX_Locs
REAL(idp),  DIMENSION(1:NQ(2)),             INTENT(IN)      ::  Cur_TX_Locs
REAL(idp),  DIMENSION(1:NQ(3)),             INTENT(IN)      ::  Cur_PX_Locs
INTEGER,                                    INTENT(IN)      ::  lvl
INTEGER,                                    INTENT(IN)      ::  iU
REAL(idp),  DIMENSION(1:NQ(2),                              &
                      1:LM_Short_Length,                    &
                      0:AMReX_Max_Grid_Size(2)-1),          &
                                            INTENT(IN)      ::  Plm_Table
                                            
REAL(idp),  DIMENSION(1:NQ(3),                              &
                      1:LM_Length,                          &
                      0:AMReX_Max_Grid_Size(3)-1)           ::  Am_Table
REAL(idp),  DIMENSION(1:NQ(1)*NQ(2)*NQ(3)), INTENT(OUT)     ::  Var_Holder

Var_Holder = 0.0_idp

END SUBROUTINE Poseidon_Outside_Type_A



 !+201+########################################################!
!                                                               !
!          Poseidon_Outside_Type_B                              !
!                                                               !
 !#############################################################!
SUBROUTINE Poseidon_Outside_Type_B( iE, iEL, iNE, NQ,    &
                                    Cur_RX_Locs,    &
                                    Cur_TX_Locs,    &
                                    Cur_PX_Locs,    &
                                    lvl,            &
                                    iU, iVB,        &
                                    Am_Table,       &
                                    Plm_Table,      &
                                    Var_Holder      )

INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  iE
INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  NQ
INTEGER,    DIMENSION(1:3),                 INTENT(IN)      ::  iEL
INTEGER,    DIMENSION(3),                   INTENT(IN)      ::  iNE
REAL(idp),  DIMENSION(1:NQ(1)),             INTENT(IN)      ::  Cur_RX_Locs
REAL(idp),  DIMENSION(1:NQ(2)),             INTENT(IN)      ::  Cur_TX_Locs
REAL(idp),  DIMENSION(1:NQ(3)),             INTENT(IN)      ::  Cur_PX_Locs
INTEGER,                                    INTENT(IN)      ::  lvl
INTEGER,                                    INTENT(IN)      ::  iU, iVB
REAL(idp),  DIMENSION(1:NQ(2),                              &
                      1:LM_Short_Length,                    &
                      0:AMReX_Max_Grid_Size(2)-1),          &
                                            INTENT(IN)      ::  Plm_Table
                                            
REAL(idp),  DIMENSION(1:NQ(3),                              &
                      1:LM_Length,                          &
                      0:AMReX_Max_Grid_Size(3)-1)           ::  Am_Table
REAL(idp),  DIMENSION(1:NQ(1)*NQ(2)*NQ(3)), INTENT(OUT)     ::  Var_Holder

Var_Holder = 0.0_idp

END SUBROUTINE Poseidon_Outside_Type_B


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

END MODULE Return_Functions_AMReX
