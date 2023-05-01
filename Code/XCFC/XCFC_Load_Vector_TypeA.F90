  !################################################################################!
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE XCFC_Load_Vector_TypeA_Module                                              !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  !################################################################################!

!*D*================================!
!                                   !
!           Dependencies            !
!                                   !
!===================================!
USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Poseidon_Message_Routines_Module, &
            ONLY :  Run_Message

USE Poseidon_Numbers_Module, &
            ONLY :  pi,                         &
                    TwoPi

USE Poseidon_Units_Module, &
            ONLY :  GR_Source_Scalar

USE Poseidon_Parameters, &
            ONLY :  Degree,                     &
                    L_Limit,                    &
                    Verbose_Flag

USE Parameters_Variable_Indices, &
            ONLY :  iU_CF,                        &
                    iU_LF,                        &
                    iU_S1,                        &
                    iU_S2,                        &
                    iU_S3,                        &
                    iU_X1,                        &
                    iU_X2,                        &
                    iU_X3,                       &
                    iVB_S,                       &
                    iVB_X

USE Poseidon_IO_Parameters, &
            ONLY :  CFA_Var_Names

USE Variables_Quadrature, &
            ONLY :  NUM_R_QUAD_POINTS,          &
                    NUM_T_QUAD_POINTS,          &
                    NUM_P_QUAD_POINTS,          &
                    NUM_TP_QUAD_POINTS,         &
                    INT_R_LOCATIONS,            &
                    INT_T_LOCATIONS,            &
                    INT_P_LOCATIONS

USE Variables_Mesh, &
            ONLY :  NUM_R_ELEMENTS,             &
                    NUM_T_ELEMENTS,             &
                    NUM_P_ELEMENTS,             &
                    rlocs,                      &
                    tlocs,                      &
                    plocs
                 
USE Variables_Source, &
            ONLY :  Block_Source_E,             &
                    Block_Source_S,             &
                    Block_Source_Si

USE Variables_Tables, &
            ONLY :  Plm_Values,                 &
                    Plm_dt_Values,              &
                    Am_Values,                  &
                    Am_dp_Values,               &
                    Slm_Elem_Values,            &
                    Slm_Elem_dt_Values,         &
                    Slm_Elem_dp_Values,         &
                    Level_dx,                   &
                    Level_Ratios,               &
                    Lagrange_Poly_Table

USE Variables_Derived, &
            ONLY :  LM_Length,                  &
                    LM_Short_Length

USE Variables_Vectors, &
            ONLY :  dVA_Load_Vector

USE Functions_Jacobian, &
            ONLY :  Calc_Ahat

USE Maps_Domain, &
            ONLY :  Map_To_FEM_Node,            &
                    FEM_Elem_Map

USE Maps_Quadrature, &
            ONLY :  Map_To_tpd

USE XCFC_Source_Routine_Variables_Module, &
            ONLY :  Cur_R_Locs,         &
                    Cur_T_Locs,         &
                    R_Square,           &
                    Sin_Square,         &
                    RSin_Square,        &
                    TP_Sin_Val,         &
                    TP_Sin_Square,      &
                    TP_Cotan_Val,       &
                    TP_RSin_Square,     &
                    R_Int_Weights,      &
                    TP_Int_Weights,     &
                    Cur_Val_Psi,        &
                    Cur_Val_AlphaPsi,   &
                    Cur_Val_X,          &
                    Cur_Drv_X,          &
                    SourceTerm
                    
USE XCFC_Functions_Physical_Source_Module, &
            ONLY : Get_Physical_Source

USE XCFC_Functions_Calc_Values_Module, &
            ONLY :  Calc_Int_Weights,               &
                    Calc_Int_Weights_AMReX,         &
                    Calc_Val_On_Elem_TypeA,         &
                    Calc_Val_And_Drv_On_Elem_TypeB

USE Initialization_Tables_Slm, &
            ONLY :  Initialize_Am_Tables,            &
                    Initialize_Plm_Tables,           &
                    Initialize_Slm_Tables_on_Elem,  &
                    Initialize_Slm_Tables
                    
USE Timer_Routines_Module, &
            ONLY :  TimerStart,                     &
                    TimerStop

USE Timer_Variables_Module, &
            ONLY :  Timer_Lapse_SourceVector,  &
                    Timer_ConFactor_SourceVector

USE MPI



#ifdef POSEIDON_AMREX_FLAG
use amrex_base_module

USE amrex_box_module,   &
            ONLY:   amrex_box

USE amrex_boxarray_module, &
            ONLY:   amrex_boxarray


USE amrex_multifab_module,  &
            ONLY:   amrex_multifab,         &
                    amrex_multifab_build,   &
                    amrex_imultifab_build

USE Variables_AMReX_Core, &
            ONLY :  MF_Source,              &
                    AMReX_Num_Levels,       &
                    AMReX_Max_Grid_Size,    &
                    Source_PTR,             &
                    Mask_PTR
                    
USE Parameters_AMReX, &
            ONLY :  iTrunk,             &
                    iLeaf

USE Poseidon_AMReX_MakeFineMask_Module, &
            ONLY :  AMReX_MakeFineMask

USE Poseidon_AMReX_Multilayer_Utilities_Module, &
            ONLY :  Find_Coarsest_Parent

#endif




IMPLICIT NONE


CONTAINS

!+102+###########################################################################!
!                                                                                !
!           XCFC_Calc_Load_Vector_TypeA                                        !
!                                                                                !
!################################################################################!
SUBROUTINE XCFC_Calc_Load_Vector_TypeA( iU, iEU, iEL )

INTEGER, INTENT(IN)                     ::  iU
INTEGER, INTENT(IN), DIMENSION(3)       ::  iEU
INTEGER, INTENT(IN), DIMENSION(3)       ::  iEL

CHARACTER(LEN = 300)                    ::  Message


#ifndef POSEIDON_AMREX_FLAG
INTEGER,             DIMENSION(3)       ::  iE
INTEGER                                 ::  re, te, pe
#endif

IF ( Verbose_Flag ) THEN
    WRITE(Message,'(A,A,A)')'Calculating ',TRIM(CFA_Var_Names(iU)),' Load Vector.'
    CALL Run_Message(TRIM(Message))
END IF



IF ( iU == iU_CF ) THEN
    CALL TimerStart( Timer_ConFactor_SourceVector)
ELSEIF ( iU == iU_LF ) THEN
    CALL TimerStart( Timer_Lapse_SourceVector)
END IF




#ifdef POSEIDON_AMREX_FLAG

    CALL XCFC_AMReX_Calc_Load_Vector_TypeA( iU )

#else
    
    dVA_Load_Vector(:,:,iU) = 0.0_idp
    
    CALL Initialize_Slm_Tables()

    DO re = iEL(1),iEU(1)
    DO te = iEL(2),iEU(2)
    DO pe = iEL(3),iEU(3)
        iE = [re,te,pe]
        CALL XCFC_Calc_Load_Vector_On_Element_TypeA( iU, iE, ELo_Opt = [0,0,0])
    END DO ! pe
    END DO ! te
    END DO ! re
#endif



IF ( iU == iU_CF ) THEN
    CALL TimerStop( Timer_ConFactor_SourceVector)
ELSEIF ( iU == iU_LF ) THEN
    CALL TimerStop( Timer_Lapse_SourceVector)
END IF




END SUBROUTINE XCFC_Calc_Load_Vector_TypeA



!+102+###########################################################################!
!                                                                                !
!           XCFC_Calc_Load_Vector_TypeA                                        !
!                                                                                !
!################################################################################!
SUBROUTINE XCFC_Calc_Load_Vector_On_Element_TypeA( iU, iE, Level_Option, iNE_Opt, ELo_Opt )

INTEGER, INTENT(IN)                                 ::  iU
INTEGER, INTENT(IN), DIMENSION(3)                   ::  iE
INTEGER, INTENT(IN),                    OPTIONAL    ::  Level_Option
INTEGER, INTENT(IN), DIMENSION(3),      OPTIONAL    ::  iNE_Opt
INTEGER, INTENT(IN), DIMENSION(3),      OPTIONAL    ::  ELo_Opt

INTEGER                                             ::  rd, tpd, td, pd

INTEGER                                             ::  FEM_Elem
REAL(KIND = idp)                                    ::  DROT,     &
                                                        DTOT

INTEGER                                             ::  Level, i
INTEGER                                             ::  iCE(3)
INTEGER                                             ::  iRE(3)
INTEGER                                             ::  iEOff(3)
INTEGER                                             ::  iNE(3)
INTEGER                                             ::  ELo(3)

IF (Present(Level_Option)) THEN
    Level = Level_Option
ELSE
    Level = 0
END IF

IF (Present(iNE_Opt) ) THEN
    iNE = iNE_Opt
ELSE
    iNE = [Num_R_Elements, Num_T_Elements, Num_P_Elements]
END IF

IF (Present(ELo_Opt) ) THEN
    ELo = ELo_Opt
ELSE
    ELo = [1, 1, 1]
END IF


#ifdef POSEIDON_AMREX_FLAG
IF ( amrex_spacedim == 1 ) THEN
    iEoff(2:3) = 0
ELSEIF ( amrex_spacedim == 2) THEN
    iEoff(2)   = iE(2)
    iEoff(3)   = 0
ELSEIF ( amrex_spacedim == 3 ) THEN
    iEoff(2:3) = iE(2:3)
END IF

DO i = 1,3
    iCE(i) = Find_Coarsest_Parent(iE(i), Level)
    iRE(i) = 2.0_idp*MOD(iE(i),Level_Ratios(Level))
END DO

FEM_Elem = FEM_Elem_Map(iE(1),Level)
DROT = Level_dx(Level,1)/2.0_idp
DTOT = Level_dx(Level,2)/2.0_idp

CUR_R_LOCS(:) = DROT * (Int_R_Locations(:) + 1.0_idp + iE(1)*2.0_idp)
CUR_T_LOCS(:) = DTOT * (Int_T_Locations(:) + 1.0_idp + iEOff(2)*2.0_idp)




#else

FEM_Elem = iE(1)
DROT = 0.5_idp * (rlocs(iE(1)+1) - rlocs(iE(1)))
DTOT = 0.5_idp * (tlocs(iE(2)+1) - tlocs(iE(2)))

CUR_R_LOCS(:) = DROT * (INT_R_LOCATIONS(:)+1.0_idp) + rlocs(iE(1))
CUR_T_LOCS(:) = DTOT * (INT_T_LOCATIONS(:)+1.0_idp) + tlocs(iE(2))



#endif


CALL Initialize_Slm_Tables_on_Elem( iE(2), iE(3),       &
                                    Num_T_Quad_Points,  &
                                    Num_P_Quad_Points,  &
                                    iNE,                &
                                    ELo,                &
                                    Plm_Values,         &
                                    Plm_dt_Values,      &
                                    Am_Values,          &
                                    Am_dp_Values,       &
                                    Slm_Elem_Values,    &
                                    Slm_Elem_dt_Values, &
                                    Slm_Elem_dp_Values  )


R_SQUARE(:) = CUR_R_LOCS(:)*CUR_R_LOCS(:)
DO td = 1,NUM_T_QUAD_POINTS
DO pd = 1,NUM_P_QUAD_POINTS
    tpd = Map_To_tpd(td,pd)
    TP_Sin_Val(tpd)    = DSIN(CUR_T_LOCS(td))
    TP_Cotan_Val(tpd)  = 1.0_idp/DTAN(CUR_T_LOCS(td))
END DO
END DO
TP_Sin_Square(:) = TP_Sin_Val(:)*TP_Sin_Val



DO rd = 1,NUM_R_QUAD_POINTS
    TP_RSIN_SQUARE(:,rd) = R_SQUARE(rd)*TP_SIN_SQUARE(:)
END DO


CALL Calc_Int_Weights( DROT, DTOT,                  &
                       R_Square, TP_Sin_Val,        &
                       R_Int_Weights, TP_Int_Weights )

CALL Calc_XCFC_CurVals_TypeA( iE,       &
                              iU,       &
                              DROT,     &
                              DTOT,     &
                              Level     )

CALL Create_XCFC_Vector_TypeA( iE, iU, Level, FEM_Elem )


END SUBROUTINE XCFC_Calc_Load_Vector_On_Element_TypeA











!+201+##########################################################################!
!                                                                               !
!                  Create_XCFC_Vector_TypeA                                     !
!                                                                               !
!###############################################################################!
SUBROUTINE Create_XCFC_Vector_TypeA( iE, iU, Level, FEM_Elem )


INTEGER, INTENT(IN), DIMENSION(3)                           ::  iE
INTEGER, INTENT(IN)                                         ::  iU
INTEGER, INTENT(IN)                                         ::  Level
INTEGER, INTENT(IN)                                         ::  FEM_Elem

INTEGER                                                     ::  rd, d, lm_loc
INTEGER                                                     ::  Current_i_Location

REAL(idp)                                                   ::  RHS_TMP
INTEGER                                                     ::  iCT

!iCT = 2**(level+1) - mod(iE(1),2**level) - 2
iCT = 0


DO lm_loc = 1,LM_LENGTH
DO d = 0,DEGREE
    
    RHS_TMP = 0.0_idp
    DO rd = 1,NUM_R_QUAD_POINTS


        RHS_TMP =  RHS_TMP                                          &
                 + SUM( SourceTerm( :, rd, iU )                     &
                       * Slm_Elem_Values( lm_loc, :)                &
                       * TP_Int_Weights(:)                     )    &
               * Lagrange_Poly_Table( d, rd, 0)                     &
               * R_Int_Weights(rd)

    END DO  ! rd Loop
    


    Current_i_Location = Map_To_FEM_Node(FEM_Elem,d)
    dVA_Load_Vector(Current_i_Location,lm_loc,iU)          &
        = dVA_Load_Vector(Current_i_Location,lm_loc,iU)    &
        + RHS_TMP


END DO  ! d Loop
END DO  ! lm_loc Loop






END SUBROUTINE Create_XCFC_Vector_TypeA










!+202+###########################################################################!
!                                                                                !
!                  Calc_Current_Values          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_XCFC_CurVals_TypeA( iE, iU, DROT, DTOT, Level )

INTEGER, INTENT(IN), DIMENSION(3)                               ::  iE
INTEGER, INTENT(IN)                                             ::  iU
REAL(KIND = idp), INTENT(IN)                                    ::  DROT, DTOT
INTEGER, INTENT(IN)                                             ::  Level


REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3,3) ::  Ahat_Array
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3 )  ::  f
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points )    ::  AA_Array
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points)     ::  PhysSrc

INTEGER                                                         ::  tpd, rd
INTEGER                                                         ::  i, j


CALL Calc_Int_Weights( DROT, DTOT,                      &
                       R_Square, TP_Sin_Val,               &
                       R_Int_Weights, TP_Int_Weights    )

CALL Calc_Val_On_Elem_TypeA( iE, Cur_Val_Psi, iU_CF, Level )
CALL Calc_Val_On_Elem_TypeA( iE, Cur_Val_AlphaPsi, iU_LF, Level )

CALL Calc_Val_And_Drv_On_Elem_TypeB(iE, DROT,      &
                                    CUR_Val_X(:,:,1),       &
                                    CUR_DRV_X(:,:,:,1),     &
                                    iU_X1, iVB_X,           &
                                    Level                   )

CALL Calc_Val_And_Drv_On_Elem_TypeB( iE, DROT,      &
                                    CUR_Val_X(:,:,2),       &
                                    CUR_DRV_X(:,:,:,2),     &
                                    iU_X2, iVB_X,           &
                                    Level                   )

CALL Calc_Val_And_Drv_On_Elem_TypeB( iE, DROT,      &
                                    CUR_Val_X(:,:,3),       &
                                    CUR_DRV_X(:,:,:,3),     &
                                    iU_X3, iVB_X,           &
                                    Level                   )



CALL Calc_Ahat( Ahat_Array,                             &
                NUM_R_QUAD_POINTS, NUM_TP_QUAD_POINTS,  &
                Cur_R_Locs, R_SQUARE,                   &
                TP_RSIN_SQUARE, TP_COTAN_VAL,           &
                CUR_VAL_X, CUR_DRV_X                  )



DO rd = 1,NUM_R_QUAD_POINTS
DO tpd = 1,NUM_TP_QUAD_POINTS

    ! Calc Flat Metric Terms
    f(tpd,rd,1) = 1.0_idp
    f(tpd,rd,2) = R_Square(rd)
    f(tpd,rd,3) = TP_RSIN_SQUARE(tpd,rd)
    
END DO ! tpd
END DO ! rd


AA_Array = 0.0_idp
DO i = 1,3
DO j = 1,3

    AA_Array(:,:) = AA_Array(:,:)           &
        + f(:,:,i)*f(:,:,j) * (Ahat_Array(:,:,i,j))**2

END DO ! i
END DO ! j




CALL Get_Physical_Source( PhysSrc, iU, iE )


IF ( iU == iU_CF ) THEN


!   XCFC Source
!
    SourceTerm(:,:,iU) = -TwoPi * GR_Source_Scalar / Cur_Val_Psi(:,:)   &
                        * PhysSrc(:,:)                                  &
                      -1.0_idp / ( 8.0_idp * Cur_Val_Psi(:,:)**7)       &
                        * AA_Array(:,:)




!   CFA Source
!
!    SourceTerm(:,:,iU) = -TwoPi * GR_Source_Scalar * Cur_Val_Psi(:,:)**5   &
!                        * PhysSrc(:,:)                                  &
!                      -1.0_idp / ( 8.0_idp * Cur_Val_Psi(:,:)**7)       &
!                        * AA_Array(:,:)


ELSEIF ( iU == iU_LF ) THEN


    SourceTerm(:,:,iU) = TwoPi * GR_Source_Scalar * Cur_Val_AlphaPsi(:,:)   &
                        / (Cur_Val_Psi(:,:)**2) * PhysSrc(:,:)              &
                    + 7.0_idp*Cur_Val_AlphaPsi(:,:)                         &
                        / ( 8.0_idp * Cur_Val_Psi(:,:)**8) * AA_Array(:,:)

ENDIF



END SUBROUTINE Calc_XCFC_CurVals_TypeA






!+102+##########################################################################!
!                                                                               !
!          XCFC_AMReX_Calc_Load_Vector_TypeA                                  !
!                                                                               !
!###############################################################################!
SUBROUTINE XCFC_AMReX_Calc_Load_Vector_TypeA( iU )

INTEGER, INTENT(IN)                             ::  iU

#ifdef POSEIDON_AMREX_FLAG

TYPE(amrex_mfiter)                              ::  mfi
TYPE(amrex_box)                                 ::  Box

TYPE(amrex_imultifab)                           ::  Level_Mask

INTEGER                                         ::  re, te, pe
INTEGER, DIMENSION(3)                           ::  iE
INTEGER, DIMENSION(3)                           ::  iEL
INTEGER, DIMENSION(3)                           ::  iEU
INTEGER, DIMENSION(3)                           ::  iEL_Off
INTEGER, DIMENSION(3)                           ::  iEU_Off
INTEGER                                         ::  nComp
INTEGER                                         ::  lvl

INTEGER, DIMENSION(1:3)                         ::  nGhost_Vec

INTEGER,    DIMENSION(3)                                ::  iNE
REAL(idp),  DIMENSION(0:AMReX_Max_Grid_Size(2))         ::  tlocs_subarray
REAL(idp),  DIMENSION(0:AMReX_Max_Grid_Size(3))         ::  plocs_subarray


nGhost_Vec = 0


dVA_Load_Vector(:,:,iU) = 0.0_idp

DO lvl = AMReX_Num_Levels-1,0,-1




    IF ( lvl < AMReX_Num_Levels-1 ) THEN
        CALL AMReX_MakeFineMask(  Level_Mask,               &
                                  MF_Source(lvl)%ba,        &
                                  MF_Source(lvl)%dm,        &
                                  nGhost_Vec,               &
                                  MF_Source(lvl+1)%ba,      &
                                  iLeaf, iTrunk            )
    ELSE
        ! Create Level_Mask all equal to 1
        CALL amrex_imultifab_build( Level_Mask,             &
                                    MF_Source(lvl)%ba,      &
                                    MF_Source(lvl)%dm,      &
                                    1,                      &
                                    0                       )
        CALL Level_Mask%SetVal(iLeaf)
    END IF



    CALL amrex_mfiter_build(mfi, MF_Source(lvl), tiling = .false. )
    DO WHILE(mfi%next())

        Source_PTR => MF_Source(lvl)%dataPtr(mfi)
        Mask_PTR   => Level_Mask%dataPtr(mfi)


        Box = mfi%tilebox()

        nComp =  MF_Source(lvl)%ncomp()

        iEL = Box%lo
        iEU = Box%hi

        iNE = iEU-iEL+1
        
        IF ( amrex_spacedim == 1 ) THEN
            iEL_off(2:3) = 0
            iEU_off(2:3) = 0
        ELSEIF ( amrex_spacedim == 2) THEN
            iEL_off(2)   = iEL(2)
            iEL_off(3)   = 0
            iEU_off(2)   = iEU(2)
            iEU_off(3)   = 0
        ELSEIF ( amrex_spacedim == 3 ) THEN
            iEL_off(2:3) = iEL(2:3)
            iEU_off(2:3) = iEU(2:3)
        END IF
        
        
        DO te = iEL_Off(2),iEU_Off(2)+1
            tlocs_subarray(te-iEL_Off(2)) = Level_dx(lvl,2)*te
        END DO
        DO pe = iEL_Off(3),iEU_Off(3)+1
            plocs_subarray(pe-iEL_Off(3)) = Level_dx(lvl,3)*pe
        END DO
        
        
        ! Initialize Am Table
        CALL Initialize_Am_Tables(  Num_P_Quad_Points,          &
                                    Int_P_Locations,            &
                                    L_Limit,                    &
                                    iNE(3),                     &
                                    [iEL_Off(3), iEU_Off(3)],   &
                                    plocs_subarray(0:iNE(3)),   &
                                    Am_Values,                  &
                                    Am_dp_Values                )

        ! Initialize Plm Table
        CALL Initialize_Plm_Tables( Num_T_Quad_Points,          &
                                    Int_T_Locations,            &
                                    L_Limit,                    &
                                    LM_Short_Length,            &
                                    iNE(2),                     &
                                    [iEL_Off(2), iEU_Off(2)],   &
                                    tlocs_subarray(0:iNE(2)),   &
                                    Plm_Values,                 &
                                    Plm_dt_Values               )

        DO re = iEL(1),iEU(1)
        DO te = iEL(2),iEU(2)
        DO pe = iEL(3),iEU(3)
            
            IF ( Mask_PTR(RE,TE,PE,1) == iLeaf ) THEN

                iE = [re,te,pe]
                CALL XCFC_Calc_Load_Vector_On_Element_TypeA( iU,        &
                                                             iE,        &
                                                             lvl,       &
                                                             iNE,       &
                                                             iEl        )
            END IF
        END DO ! pe
        END DO ! te
        END DO ! re

    END DO
    
    
    CALL amrex_mfiter_destroy(mfi)
    CALL amrex_imultifab_destroy( Level_Mask )
    Source_PTR => Null()
    Mask_PTR => Null()
    
    

END DO ! lvl

#endif

END SUBROUTINE XCFC_AMReX_Calc_Load_Vector_TypeA



END MODULE XCFC_Load_Vector_TypeA_Module
