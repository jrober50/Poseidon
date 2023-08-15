   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE ADM_Mass_Module                                                       !##!
!##!                                                                         !##!
!##!_________________________________________________________________________!##!
!##!                                                                         !##!
!##!                                                                         !##!
!##!=========================================================================!##!
!##!                                                                         !##!
!##!    Contains:                                                            !##!
!##!                                                                         !##!
!##!        +101+                   Calc_ADM_Mass                            !##!
!##!                                                                         !##!
!##!        +201+                   Calc_ADM_Mass_Native                     !##!
!##!        +202+                   Calc_ADM_Mass_AMReX                      !##!
!##!                                                                         !##!
!##!        +301+                   Calc_ADM_Mass_On_Element                 !##!
!##!                                                                         !##!
!##!        +401+                   Calc_Cur_Values                          !##!
!##!        +402+                   Calc_Int_Weights                         !##!
!##!        +403+                   Calc_Int_Source_XCFC                     !##!
!##!        +404+                   Calc_Int_Source_Newtonian                !##!
!##!        +405+                   Calc_Element_ADM_Integral                !##!
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

USE Poseidon_Numbers_Module, &
            ONLY :  Pi

USE Poseidon_Parameters, &
            ONLY :  Degree,                     &
                    L_Limit

USE Poseidon_Units_Module, &
            ONLY :  GR_Source_Scalar,           &
                    C_Square

USE Parameters_Variable_Indices, &
            ONLY :  iU_CF,                      &
                    iU_X1,                      &
                    iU_X2,                      &
                    iU_X3,                      &
                    iVB_X

USE Variables_Vectors, &
            ONLY :  dVB_Coeff_Vector

USE Variables_Quadrature, &
            ONLY :  NUM_R_QUAD_POINTS,          &
                    NUM_T_QUAD_POINTS,          &
                    NUM_P_QUAD_POINTS,          &
                    NUM_TP_QUAD_POINTS,         &
                    INT_R_LOCATIONS,            &
                    INT_T_LOCATIONS,            &
                    INT_P_LOCATIONS,            &
                    INT_R_WEIGHTS,              &
                    INT_T_WEIGHTS,              &
                    INT_P_WEIGHTS,              &
                    INT_TP_WEIGHTS

USE Variables_Mesh, &
            ONLY :  NUM_R_ELEMENTS,             &
                    NUM_T_ELEMENTS,             &
                    NUM_P_ELEMENTS,             &
                    rlocs,                      &
                    tlocs,                      &
                    plocs

USE Variables_Source, &
            ONLY :  Block_Source_E,             &
                    Source_Rho
                    
USE Variables_Derived, &
            ONLY :  LM_Length,                  &
                    LM_Short_Length

USE Maps_Quadrature, &
            ONLY :  Map_To_tpd,                 &
                    Quad_Map

USE Maps_Domain, &
            ONLY :  FEM_Elem_Map

USE XCFC_Functions_Calc_Values_Module, &
            ONLY :  Calc_Val_On_Elem_TypeA,         &
                    Calc_Val_And_Drv_On_Elem_TypeB

USE Functions_Jacobian, &
            ONLY :  Calc_Ahat


USE Variables_Tables, &
            ONLY :  Plm_Values,                 &
                    Plm_dt_Values,              &
                    Am_Values,                  &
                    Am_dp_Values,               &
                    Slm_Elem_Values,            &
                    Slm_Elem_dt_Values,         &
                    Slm_Elem_dp_Values,         &
                    Level_dx

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
            ONLY :  iTrunk,                 &
                    iLeaf

USE Poseidon_AMReX_MakeFineMask_Module, &
            ONLY :  AMReX_MakeFineMask

USE Poseidon_AMReX_Multilayer_Utilities_Module, &
            ONLY :  Find_Coarsest_Parent

#endif


USE Variables_MPI, &
            ONLY :  POSEIDON_COMM_WORLD,        &
                    myID_Poseidon,              &
                    iErr

USE Initialization_Tables_Slm, &
            ONLY :  Initialize_Am_Tables,            &
                    Initialize_Plm_Tables,           &
                    Initialize_Slm_Tables_on_Elem


USE XCFC_Functions_Physical_Source_Module, &
            ONLY : Get_Physical_Source

USE Flags_Core_Module, &
            ONLY :  iPF_Core_Flags,         &
                    iPF_Core_Method_Mode,   &
                    iPF_Core_Method_Newtonian

USE MPI

IMPLICIT NONE




CONTAINS
!+101+##################################################################!
!                                                                       !
!          Calc_ADM_Mass                                                  !
!                                                                       !
!#######################################################################!
SUBROUTINE Calc_ADM_Mass( ADM_Mass )

REAL(idp), INTENT(OUT)                              ::  ADM_Mass


#ifdef POSEIDON_AMREX_FLAG

CALL Calc_ADM_Mass_AMReX(ADM_Mass)

#else

CALL Calc_ADM_Mass_Native(ADM_Mass)

#endif

CALL MPI_ALLREDUCE( MPI_IN_PLACE,       &
                    ADM_Mass,           &
                    1,                  &
                    MPI_DOUBLE,         &
                    MPI_SUM,            &
                    POSEIDON_COMM_WORLD,&
                    ierr                )



END SUBROUTINE Calc_ADM_Mass














!+201+##################################################################!
!                                                                       !
!          Calc_ADM_Mass_Native                                         !
!                                                                       !
!#######################################################################!
SUBROUTINE Calc_ADM_Mass_Native( ADM_Mass )

REAL(idp), INTENT(OUT)                              ::  ADM_Mass



INTEGER                                             ::  re, te, pe
REAL(idp)                                           ::  Int_Val


ADM_Mass = 0.0_idp


DO re = 0,Num_R_Elements-1
DO te = 0,Num_T_Elements-1
DO pe = 0,Num_P_Elements-1


    CALL Calc_ADM_Mass_On_Element( [re, te, pe], Int_Val )

    ADM_Mass = ADM_Mass + Int_Val

END DO ! pe Loop
END DO ! te Loop
END DO ! re Loop


END SUBROUTINE Calc_ADM_Mass_Native



!+202+##################################################################!
!                                                                       !
!          Calc_ADM_Mass_AMReX                                          !
!                                                                       !
!#######################################################################!
SUBROUTINE Calc_ADM_Mass_AMReX( ADM_Mass )

REAL(idp), INTENT(OUT)                          ::  ADM_Mass

#ifdef POSEIDON_AMREX_FLAG

TYPE(amrex_mfiter)                              ::  mfi
TYPE(amrex_box)                                 ::  Box

TYPE(amrex_imultifab)                           ::  Level_Mask

INTEGER                                         ::  re, te, pe
INTEGER, DIMENSION(3)                           ::  iE
INTEGER, DIMENSION(3)                           ::  iEL, iEU
INTEGER                                         ::  nComp
INTEGER                                         ::  lvl

REAL(idp)                                       ::  Int_Val

INTEGER, DIMENSION(1:3)                         ::  nGhost_Vec

INTEGER,    DIMENSION(3)                                ::  iNE
REAL(idp),  DIMENSION(0:AMReX_Max_Grid_Size(2)-1)       ::  tlocs_subarray
REAL(idp),  DIMENSION(0:AMReX_Max_Grid_Size(3)-1)       ::  plocs_subarray



nGhost_Vec = 0

ADM_Mass = 0.0_idp
Int_Val  = 0.0_idp

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
        
        
        DO te = iEL(2),iEU(2)
            tlocs_subarray(te-iEL(2)) = Level_dx(lvl,2)*te
        END DO
        DO pe = iEL(3),iEU(3)
            plocs_subarray(pe-iEL(3)) = Level_dx(lvl,3)*pe
        END DO
        
        
        ! Initialize Am Table
        CALL Initialize_Am_Tables(  Num_P_Quad_Points,          &
                                    Int_P_Locations,            &
                                    L_Limit,                    &
                                    iNE(3),                     &
                                    [iEL(3), iEU(3)],           &
                                    plocs_subarray(0:iNE(3)-1), &
                                    Am_Values,                  &
                                    Am_dp_Values                )

        ! Initialize Plm Table
        CALL Initialize_Plm_Tables( Num_T_Quad_Points,          &
                                    Int_T_Locations,            &
                                    L_Limit,                    &
                                    LM_Short_Length,            &
                                    iNE(2),                     &
                                    [iEL(2), iEU(2)],           &
                                    tlocs_subarray(0:iNE(2)-1), &
                                    Plm_Values,                 &
                                    Plm_dt_Values               )


        DO re = iEL(1),iEU(1)
        DO te = iEL(2),iEU(2)
        DO pe = iEL(3),iEU(3)
            
            IF ( Mask_PTR(RE,TE,PE,1) == iLeaf ) THEN
                CALL Initialize_Slm_Tables_on_Elem( te, pe,             &
                                                    Num_T_Quad_Points,  &
                                                    Num_P_Quad_Points,  &
                                                    iNE,                &
                                                    iEL,                &
                                                    Plm_Values,         &
                                                    Plm_dt_Values,      &
                                                    Am_Values,          &
                                                    Am_dp_Values,       &
                                                    Slm_Elem_Values,    &
                                                    Slm_Elem_dt_Values, &
                                                    Slm_Elem_dp_Values  )
                iE = [re,te,pe]
                CALL Calc_ADM_Mass_On_Element( iE, Int_Val, lvl )

                ADM_Mass = ADM_Mass + Int_Val
            END IF
        END DO ! pe
        END DO ! te
        END DO ! re

    END DO
    CALL amrex_mfiter_destroy(mfi)
END DO ! lvl

#endif

END SUBROUTINE Calc_ADM_Mass_AMReX


!+301+##################################################################!
!                                                                       !
!          Calc_ADM_Mass_On_Element                                     !
!                                                                       !
!#######################################################################!
SUBROUTINE Calc_ADM_Mass_On_Element( iE, Int_Val, Level_Option )

INTEGER,    INTENT(IN), DIMENSION(1:3)              ::  iE
REAL(idp),  INTENT(OUT)                             ::  Int_Val
INTEGER,    INTENT(IN), OPTIONAL                    ::  Level_Option



REAL(idp)                                           ::  DROT, DTOT, DPOT
REAL(idp), DIMENSION(1:Num_R_Quad_Points)           ::  crlocs
REAL(idp), DIMENSION(1:Num_T_Quad_Points)           ::  ctlocs
REAL(idp), DIMENSION(1:Num_R_Quad_Points)           ::  rSquare

REAL(idp), DIMENSION(1:Num_TP_Quad_Points)          ::  TP_Sin_Val
REAL(idp), DIMENSION(1:Num_TP_Quad_Points)          ::  TP_Cotan_Val

REAL(idp), DIMENSION(1:Num_TP_Quad_Points,          &
                     1:Num_R_Quad_Points)           ::  TP_Rsin_Square

REAL(idp), DIMENSION(1:Num_TP_Quad_Points,          &
                     1:Num_R_Quad_Points)           ::  Int_Weights

REAL(idp), DIMENSION(1:Num_TP_Quad_Points,          &
                     1:Num_R_Quad_Points)           ::  Int_Source

INTEGER                                             ::  Level
INTEGER                                             ::  iEOff(3)


IF (Present(Level_Option)) THEN
    Level = Level_Option
ELSE
    Level = 0
END IF


#ifdef POSEIDON_AMREX_FLAG
iEoff(1) = iE(1)
IF ( amrex_spacedim == 1 ) THEN
    iEoff(2:3) = 0
ELSEIF ( amrex_spacedim == 2) THEN
    iEoff(2)   = iE(2)
    iEoff(3)   = 0
ELSEIF ( amrex_spacedim == 3 ) THEN
    iEoff(2:3) = iE(2:3)
END IF
#else
iEoff = iE
#endif



CALL Calc_Cur_Values( iEoff,                &
                      DROT, DTOT, DPOT,     &
                      crlocs, ctlocs,       &
                      rSquare,              &
                      TP_Sin_Val,           &
                      TP_Cotan_Val,         &
                      TP_rSin_Square,       &
                      Level                 )


CALL Calc_Int_Weights( DROT, DTOT,              &
                       rSquare, TP_Sin_Val,     &
                       Int_Weights              )



IF ( iPF_Core_Flags(iPF_Core_Method_Mode) == iPF_Core_Method_Newtonian ) THEN

    CALL Calc_Int_Source_Newtonian( iE, Int_Source  )

ELSE

    CALL Calc_Int_Source_XCFC(  iE, Level,       &
                                DROT, DTOT, DPOT,&
                                crlocs, rSquare, &
                                TP_Sin_Val,      &
                                TP_Cotan_Val,    &
                                TP_rSin_Square,  &
                                Int_Source       )

END IF


CALL Calc_Element_ADM_Integral( Int_Weights,    &
                                Int_Source,     &
                                Int_Val         )







END SUBROUTINE Calc_ADM_Mass_On_Element








!+401+##################################################################!
!                                                                       !
!          Calc_Cur_Values                                              !
!                                                                       !
!#######################################################################!
SUBROUTINE Calc_Cur_Values( iE,                     &
                            DROT, DTOT, DPOT,       &
                            crlocs, ctlocs,         &
                            rSquare,                &
                            TP_Sin_Val,             &
                            TP_Cotan_Val,           &
                            TP_RSin_Square ,        &
                            Level                   )

INTEGER,   INTENT(IN),  DIMENSION(1:3)                      ::  iE

REAL(idp), INTENT(OUT)                                      ::  DROT, DTOT, DPOT
REAL(idp), INTENT(OUT), DIMENSION(1:Num_R_Quad_Points)      ::  crlocs
REAL(idp), INTENT(OUT), DIMENSION(1:Num_T_Quad_Points)      ::  ctlocs
REAL(idp), INTENT(OUT), DIMENSION(1:Num_R_Quad_Points)      ::  rSquare

REAL(idp), INTENT(OUT), DIMENSION(1:Num_TP_Quad_Points)     ::  TP_Sin_Val
REAL(idp), INTENT(OUT), DIMENSION(1:Num_TP_Quad_Points)     ::  TP_Cotan_Val

REAL(idp), INTENT(OUT), DIMENSION(1:Num_TP_Quad_Points,     &
                                  1:Num_R_Quad_Points)      ::  TP_Rsin_Square

INTEGER,   INTENT(IN)                                       ::  Level


INTEGER                                                     ::  FEM_Elem
INTEGER                                                     ::  rd, td, pd, tpd



#ifdef POSEIDON_AMREX_FLAG
    FEM_Elem = FEM_Elem_Map(iE(1),Level)
    DROT = Level_dx(Level,1)/2.0_idp
    DTOT = Level_dx(Level,2)/2.0_idp

    crlocs(:) = DROT * (Int_R_Locations(:) + 1.0_idp + iE(1)*2.0_idp)
    ctlocs(:) = DTOT * (Int_T_Locations(:) + 1.0_idp + iE(2)*2.0_idp)

#else

    FEM_Elem = iE(1)
    DROT = 0.5_idp * (rlocs(iE(1)+1) - rlocs(iE(1)))
    DTOT = 0.5_idp * (tlocs(iE(2)+1) - tlocs(iE(2)))

    crlocs(:) = DROT * (INT_R_LOCATIONS(:)+1.0_idp) + rlocs(iE(1))
    ctlocs(:) = DTOT * (INT_T_LOCATIONS(:)+1.0_idp) + tlocs(iE(2))

#endif



DPOT = 0.5_idp * (plocs(iE(3) + 1) - plocs(iE(3)))


rSquare(:) = crlocs(:)*crlocs(:)

DO td = 1,NUM_T_QUAD_POINTS
DO pd = 1,NUM_P_QUAD_POINTS
    tpd = Map_To_tpd(td,pd)
    TP_Sin_Val(tpd)    = DSIN(ctlocs(td))
    TP_Cotan_Val(tpd)  = 1.0_idp/DTAN(ctlocs(td))
END DO
END DO


DO rd = 1,NUM_R_QUAD_POINTS
    TP_RSIN_SQUARE(:,rd) = rSquare(rd)*TP_Sin_Val(:)*TP_Sin_Val(:)
END DO



END SUBROUTINE Calc_Cur_Values





!+402+##################################################################!
!                                                                       !
!          Calc_Int_Weights                                             !
!                                                                       !
!#######################################################################!
SUBROUTINE Calc_Int_Weights( DROT, DTOT,            &
                             R_Square, Sin_Val,     &
                             Int_Weights            )


REAL(idp), INTENT(IN)                                       ::  DROT, DTOT
REAL(idp), INTENT(IN),    DIMENSION(1:Num_R_Quad_Points)    ::  R_Square
REAL(idp), INTENT(IN),    DIMENSION(1:Num_TP_Quad_Points)   ::  Sin_Val
REAL(idp), INTENT(INOUT), DIMENSION(1:Num_TP_Quad_Points,   &
                                    1:Num_R_Quad_Points)    ::  Int_Weights

INTEGER                                                     ::  tpd, rd, td, pd



DO rd = 1,Num_R_Quad_Points
DO td = 1,NUM_T_QUAD_POINTS
DO pd = 1,NUM_P_QUAD_POINTS

   tpd = Map_To_tpd(td,pd)
   Int_Weights( tpd, rd ) = DROT * R_SQUARE(rd) * INT_R_WEIGHTS(rd)  &
                          * DTOT * SIN_VAL(tpd) * INT_T_WEIGHTS(td)  &
                          * INT_P_WEIGHTS(pd)

END DO ! pd Loop
END DO ! td Loop
END DO ! rd Loop

END SUBROUTINE Calc_Int_Weights



!+403+##################################################################!
!                                                                       !
!          Calc_Int_Source_XCFC                                         !
!                                                                       !
!#######################################################################!
SUBROUTINE Calc_Int_Source_XCFC(iE, Level,         &
                                DROT, DTOT, DPOT,  &
                                rlocs, rSquare,    &
                                TP_Sin_Val,        &
                                TP_Cotan_Val,      &
                                TP_RSin_Square,    &
                                Int_Source         )

INTEGER,    INTENT(IN), DIMENSION(1:3)                      ::  iE
INTEGER,    INTENT(IN)                                      ::  Level
REAL(idp),  INTENT(IN)                                      ::  DROT, DTOT, DPOT

REAL(idp),  INTENT(IN), DIMENSION(1:Num_R_Quad_Points)      ::  rlocs
REAL(idp),  INTENT(IN), DIMENSION(1:Num_R_Quad_Points)      ::  rSquare

REAL(idp),  INTENT(IN), DIMENSION(1:Num_TP_Quad_Points)     ::  TP_Sin_Val
REAL(idp),  INTENT(IN), DIMENSION(1:Num_TP_Quad_Points)     ::  TP_Cotan_Val

REAL(idp),  INTENT(IN), DIMENSION(1:Num_TP_Quad_Points,     &
                                  1:Num_R_Quad_Points)      ::  TP_Rsin_Square

REAL(idp),  INTENT(INOUT), DIMENSION(1:Num_TP_Quad_Points,  &
                                     1:Num_R_Quad_Points)   ::  Int_Source




REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3,3) ::  Ahat_Array
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3 )  ::  f
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points )    ::  AA_Array


REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points )    ::  Cur_Val_Psi
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3 )  ::  Cur_Val_X
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3,3) ::  Cur_Drv_X
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points)     ::  PhysSrc

INTEGER                                                         ::  tpd, rd, td, pd
INTEGER                                                         ::  i, j, HEre




CALL Calc_Val_On_Elem_TypeA( iE, Cur_Val_Psi, iU_CF, Level )

IF ( ALLOCATED(dVB_Coeff_Vector) ) THEN
    CALL Calc_Val_And_Drv_On_Elem_TypeB( iE, DROT,      &
                                         CUR_Val_X(:,:,1),      &
                                         CUR_DRV_X(:,:,:,1),    &
                                         iU_X1, iVB_X,          &
                                         Level                  )

    CALL Calc_Val_And_Drv_On_Elem_TypeB( iE, DROT,      &
                                         CUR_Val_X(:,:,2),      &
                                         CUR_DRV_X(:,:,:,2),    &
                                         iU_X2, iVB_X,          &
                                         Level                  )

    CALL Calc_Val_And_Drv_On_Elem_TypeB( iE, DROT,      &
                                         CUR_Val_X(:,:,3),      &
                                         CUR_DRV_X(:,:,:,3),    &
                                         iU_X3, iVB_X,          &
                                         Level                  )

ELSE

    Cur_Val_X = 0.0_idp
    Cur_Drv_X = 0.0_idp

END IF

CALL Calc_Ahat( Ahat_Array,                             &
                NUM_R_QUAD_POINTS, NUM_TP_QUAD_POINTS,  &
                rlocs, rSquare,                         &
                TP_RSIN_SQUARE, TP_COTAN_VAL,           &
                CUR_VAL_X, CUR_DRV_X                  )





DO rd = 1,NUM_R_QUAD_POINTS
DO tpd = 1,NUM_T_QUAD_POINTS

    ! Calc Flat Metric Terms
    f(tpd,rd,1) = 1.0_idp
    f(tpd,rd,2) = rSquare(rd)
    f(tpd,rd,3) = rSquare(rd) * TP_Sin_Val(tpd) * TP_Sin_Val(tpd)
    
END DO ! tpd
END DO ! rd


AA_Array = 0.0_idp
DO i = 1,3
DO j = 1,3
    AA_Array(:,:) = AA_Array(:,:)           &
        + f(:,:,i)*f(:,:,j) * (Ahat_Array(:,:,i,j))**2

END DO ! i
END DO ! j


CALL Get_Physical_Source( PhysSrc, iU_CF, iE )



DO rd = 1,NUM_R_QUAD_POINTS
DO td = 1,NUM_T_QUAD_POINTS
DO pd = 1,NUM_P_QUAD_POINTS

    tpd = Map_To_tpd(td,pd)
    Here = Quad_Map(rd,td,pd)


    Int_Source(tpd,rd) = GR_Source_Scalar                               &
                            / Cur_Val_Psi(tpd,rd)                       &
                            * PhysSrc(tpd,rd)                           &
                         + AA_Array(tpd,rd)                             &
                            / ( 16.0_idp * pi * Cur_Val_Psi(tpd,rd)**7)
                            

END DO ! pd
END DO ! td
END DO ! rd



END SUBROUTINE Calc_Int_Source_XCFC



!+404+##################################################################!
!                                                                       !
!          Calc_Int_Source_Newtonian                                    !
!                                                                       !
!#######################################################################!
SUBROUTINE Calc_Int_Source_Newtonian(   iE, Int_Source         )

INTEGER,    INTENT(IN),     DIMENSION(1:3)                      ::  iE

REAL(idp),  INTENT(INOUT), DIMENSION(1:Num_TP_Quad_Points,  &
                                     1:Num_R_Quad_Points)       ::  Int_Source


REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points )    ::  Cur_Val_Psi

INTEGER                                                         ::  rd, td, pd
INTEGER                                                         ::  tpd, Here

REAL(idp)                                                       ::  Psi


INTEGER                                                         ::  Level = 0


CALL Calc_Val_On_Elem_TypeA( iE, Cur_Val_Psi, iU_CF, Level )



DO rd = 1,NUM_R_QUAD_POINTS
DO td = 1,NUM_T_QUAD_POINTS
DO pd = 1,NUM_P_QUAD_POINTS

    tpd = Map_To_tpd(td,pd)
    Here = Quad_Map(rd,td,pd)


    Psi = 1.0_idp - Cur_Val_Psi(tpd,rd)/(2.0_idp*C_Square)

    Int_Source(tpd,rd) = GR_Source_Scalar / Psi * Source_Rho(Here,iE(1),iE(2),iE(3))
    

END DO ! pd
END DO ! td
END DO ! rd


END SUBROUTINE Calc_Int_Source_Newtonian




!+405+##################################################################!
!                                                                       !
!          Calc_Element_ADM_Integral                                    !
!                                                                       !
!#######################################################################!
SUBROUTINE Calc_Element_ADM_Integral( Int_Weights,      &
                                      Int_Source,       &
                                      Int_Val           )

REAL(idp),  INTENT(IN), DIMENSION(1:Num_TP_Quad_Points, &
                                  1:Num_R_Quad_Points)  ::  Int_Weights

REAL(idp),  INTENT(IN), DIMENSION(1:Num_TP_Quad_Points, &
                                  1:Num_R_Quad_Points)  ::  Int_Source


REAL(idp), INTENT(OUT)                                  ::  Int_Val

INTEGER                                                 ::  rd
REAL(idp)                                               ::  Tmp_Val




Tmp_Val = 0.0_idp
DO rd = 1,NUM_R_QUAD_POINTS

    Tmp_Val =  Tmp_Val                      &
             + SUM( Int_Source( :, rd )     &
                  * Int_Weights(:, rd )     )


END DO  ! rd Loop
    

Int_Val = Tmp_Val


END SUBROUTINE Calc_Element_ADM_Integral














END MODULE ADM_Mass_Module


