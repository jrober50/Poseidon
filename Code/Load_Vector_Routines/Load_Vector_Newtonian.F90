  !################################################################################!
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Load_Vector_Newtonian_Module                                                 !##!
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
            ONLY :  GR_Source_Scalar,           &
                    Grav_Constant_G

USE Poseidon_Parameters, &
            ONLY :  Degree,                     &
                    L_Limit,                    &
                    Verbose_Flag

USE Parameters_Variable_Indices, &
            ONLY :  iU_NP

USE Poseidon_IO_Parameters, &
            ONLY :  CFA_Var_Names

USE Variables_Quadrature, &
            ONLY :  NUM_R_QUAD_POINTS,          &
                    NUM_T_QUAD_POINTS,          &
                    NUM_P_QUAD_POINTS,          &
                    NUM_TP_QUAD_POINTS,         &
                    INT_R_LOCATIONS,            &
                    INT_T_LOCATIONS,            &
                    INT_P_LOCATIONS,            &
                    Int_R_Weights,              &
                    Int_T_Weights,              &
                    Int_P_Weights

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

USE Load_Vector_Variables_Module, &
            ONLY :  Cur_R_Locs,         &
                    Cur_T_Locs,         &
                    R_Square,           &
                    TP_Sin_Val,         &
                    TP_Sin_Square,      &
                    TP_RSin_Square,     &
                    R_Int_Weights,      &
                    TP_Int_Weights,     &
                    SourceTerm
                    
USE Load_Vector_Functions_Physical_Source_Module, &
            ONLY : Get_Physical_Source

USE Load_Vector_Functions_Calc_Values_Module, &
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
            ONLY :  Timer_Newtonian_LoadVector

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

 !+101+############################################################!
!                                                                   !
!          Create_Load_Vector_Newtonian                             !
!                                                                   !
 !#################################################################!
SUBROUTINE Create_Load_Vector_Newtonian( )

CHARACTER(LEN = 300)                    ::  Message

IF ( Verbose_Flag ) THEN
    WRITE(Message,'(A,A,A)')'Calculating ',TRIM(CFA_Var_Names(1)),' Load Vector.'
    CALL Run_Message(TRIM(Message))
END IF



CALL TimerStart( Timer_Newtonian_LoadVector)


#ifdef POSEIDON_AMREX_FLAG
    CALL Create_Load_Vector_Newtonian_AMReX( )
#else
    CALL Create_Load_Vector_Newtonian_Native( )
#endif


CALL TimerStop( Timer_Newtonian_LoadVector)




END SUBROUTINE Create_Load_Vector_Newtonian



 !+203+############################################################!
!                                                                   !
!          Create_Load_Vector_Newtonian_Native                      !
!                                                                   !
 !#################################################################!
SUBROUTINE Create_Load_Vector_Newtonian_Native(  )



INTEGER,    DIMENSION(3)        ::  iEU
INTEGER,    DIMENSION(3)        ::  iEL
INTEGER,    DIMENSION(3)        ::  iE
INTEGER                         ::  re, te, pe


iEL = [0, 0, 0]
iEU = [Num_R_Elements-1,Num_T_Elements-1,Num_P_Elements-1]


dVA_Load_Vector(:,:,iU_NP) = 0.0_idp


CALL Initialize_Slm_Tables()


DO re = iEL(1),iEU(1)
DO te = iEL(2),iEU(2)
DO pe = iEL(3),iEU(3)
    iE = [re,te,pe]
    CALL Calc_Load_Vector_On_Element_Newtonian( iE, ELo_Opt = [0,0,0])
END DO ! pe
END DO ! te
END DO ! re


END SUBROUTINE Create_Load_Vector_Newtonian_Native




 !+203+############################################################!
!                                                                   !
!          Create_Load_Vector_Newtonian_AMReX                       !
!                                                                   !
 !#################################################################!
SUBROUTINE Create_Load_Vector_Newtonian_AMReX( )


#ifdef POSEIDON_AMREX_FLAG

TYPE(amrex_mfiter)                              ::  mfi
TYPE(amrex_box)                                 ::  Box

TYPE(amrex_imultifab)                           ::  Level_Mask

INTEGER                                         ::  re, te, pe
INTEGER,    DIMENSION(3)                        ::  iE
INTEGER,    DIMENSION(3)                        ::  iEL
INTEGER,    DIMENSION(3)                        ::  iEU
INTEGER,    DIMENSION(3)                        ::  iEL_Off
INTEGER,    DIMENSION(3)                        ::  iEU_Off
INTEGER                                         ::  nComp
INTEGER                                         ::  lvl

INTEGER,    DIMENSION(1:3)                      ::  nGhost_Vec

INTEGER,    DIMENSION(3)                        ::  iNE
REAL(idp),  DIMENSION(0:AMReX_Max_Grid_Size(2)) ::  tlocs_subarray
REAL(idp),  DIMENSION(0:AMReX_Max_Grid_Size(3)) ::  plocs_subarray


nGhost_Vec = 0


dVA_Load_Vector(:,:,iU_NP) = 0.0_idp

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
                CALL Calc_Load_Vector_On_Element_Newtonian( iE,        &
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

END SUBROUTINE Create_Load_Vector_Newtonian_AMReX







 !+102+############################################################!
!                                                                   !
!          Calc_Load_Vector_On_Element_Newtonian                    !
!                                                                   !
 !#################################################################!
SUBROUTINE Calc_Load_Vector_On_Element_Newtonian( iE, Level_Option, iNE_Opt, ELo_Opt )

INTEGER, INTENT(IN), DIMENSION(3)                   ::  iE
INTEGER, INTENT(IN),                    OPTIONAL    ::  Level_Option
INTEGER, INTENT(IN), DIMENSION(3),      OPTIONAL    ::  iNE_Opt
INTEGER, INTENT(IN), DIMENSION(3),      OPTIONAL    ::  ELo_Opt

INTEGER                                             ::  rd, tpd, td, pd

INTEGER                                             ::  FEM_Elem
REAL(KIND = idp)                                    ::  DROT,       &
                                                        DTOT,       &
                                                        DPOT

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
DPOT = Level_dx(Level,3)/2.0_idp

Cur_R_Locs(:) = DROT * (Int_R_Locations(:) + 1.0_idp + iE(1)*2.0_idp)
Cur_T_Locs(:) = DTOT * (Int_T_Locations(:) + 1.0_idp + iEOff(2)*2.0_idp)


#else

FEM_Elem = iE(1)
DROT = 0.5_idp * (rlocs(iE(1)+1) - rlocs(iE(1)))
DTOT = 0.5_idp * (tlocs(iE(2)+1) - tlocs(iE(2)))
DPOT = 0.5_idp * (plocs(iE(3)+1) - plocs(iE(3)))

Cur_R_Locs(:) = DROT * (INT_R_LOCATIONS(:)+1.0_idp) + rlocs(iE(1))
Cur_T_Locs(:) = DTOT * (INT_T_LOCATIONS(:)+1.0_idp) + tlocs(iE(2))



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


R_Square(:) = Cur_R_Locs(:)*Cur_R_Locs(:)
DO td = 1,NUM_T_QUAD_POINTS
DO pd = 1,NUM_P_QUAD_POINTS
    tpd = Map_To_tpd(td,pd)
    TP_Sin_Val(tpd)    = DSIN(Cur_T_Locs(td))
END DO
END DO
TP_Sin_Square(:) = TP_Sin_Val(:)*TP_Sin_Val


CALL Calc_Int_Weights( DROT, DTOT, DPOT,            &
                       R_Square, TP_Sin_Val,        &
                       R_Int_Weights, TP_Int_Weights )

CALL Calc_CurVals_Newtonian( iE )


CALL Create_Vector_Newtonian( iE, Level, FEM_Elem )




END SUBROUTINE Calc_Load_Vector_On_Element_Newtonian




!+202+############################################################!
!                                                                   !
!          Calc_CurVals_Newtonian                                  !
!                                                                   !
!#################################################################!
SUBROUTINE Calc_CurVals_Newtonian( iE )

INTEGER, INTENT(IN), DIMENSION(3)                               ::  iE

REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points)     ::  PhysSrc


CALL Get_Physical_Source( PhysSrc, iU_NP, iE )

SourceTerm(:,:,iU_NP) = 4.0 * Pi * Grav_Constant_G * PhysSrc(:,:)

PRINT*,iE(1),SourceTerm(:,:,iU_NP)

END SUBROUTINE Calc_CurVals_Newtonian







 !+201+############################################################!
!                                                                   !
!          Create_Vector_Newtonian                                 !
!                                                                   !
 !#################################################################!
SUBROUTINE Create_Vector_Newtonian( iE, iU, FEM_Elem )


INTEGER, INTENT(IN), DIMENSION(3)                           ::  iE
INTEGER, INTENT(IN)                                         ::  iU
INTEGER, INTENT(IN)                                         ::  FEM_Elem

INTEGER                                                     ::  rd, d, lm_loc
INTEGER                                                     ::  Current_i_Location

REAL(idp)                                                   ::  RHS_TMP


DO lm_loc = 1,LM_LENGTH
DO d = 0,DEGREE
    
    RHS_TMP = 0.0_idp
    DO rd = 1,NUM_R_QUAD_POINTS

        RHS_TMP =  RHS_TMP                                          &
                 + SUM( SourceTerm( :, rd, iU_NP )                     &
                       * Slm_Elem_Values( lm_loc, :)                &
                       * TP_Int_Weights(:)                     )    &
               * Lagrange_Poly_Table( d, rd, 0)                     &
               * R_Int_Weights(rd)

    END DO  ! rd Loop


    Current_i_Location = Map_To_FEM_Node(FEM_Elem,d)
    dVA_Load_Vector(Current_i_Location,lm_loc,iU_NP)          &
        = dVA_Load_Vector(Current_i_Location,lm_loc,iU_NP)    &
        + RHS_TMP


END DO  ! d Loop
END DO  ! lm_loc Loop


END SUBROUTINE Create_Vector_Newtonian








END MODULE Load_Vector_Newtonian_Module

