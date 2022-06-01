  !################################################################################!
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE XCFC_Source_Vector_TypeB_Module                                              !##!
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

USE Poseidon_Numbers_Module, &
            ONLY :  pi,                         &
                    TwoPi

USE Poseidon_Units_Module, &
            ONLY :  GR_Source_Scalar

USE Poseidon_Parameters, &
            ONLY :  DEGREE,                     &
                    L_LIMIT,                    &
                    NUM_CFA_EQs,                &
                    NUM_CFA_VARs

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


USE Variables_Quadrature, &
            ONLY :  NUM_R_QUAD_POINTS,          &
                    NUM_T_QUAD_POINTS,          &
                    NUM_P_QUAD_POINTS,          &
                    NUM_TP_QUAD_POINTS,         &
                    INT_R_LOCATIONS,            &
                    INT_T_LOCATIONS,            &
                    INT_P_LOCATIONS

USE Variables_Mesh, &
            ONLY :  Num_P_Elements,             &
                    rlocs,                      &
                    drlocs,                     &
                    tlocs,                      &
                    plocs
                

USE Variables_Tables, &
            ONLY :  Ylm_CC_Values,              &
                    Ylm_Elem_CC_Values,         &
                    Lagrange_Poly_Table,        &
                    Lagpoly_MultiLayer_Table,   &
                    Level_dx,                   &
                    Level_Ratios

USE Variables_Derived, &
            ONLY :  LM_LENGTH

USE Variables_FP, &
            ONLY :  FP_Coeff_Vector_A,           &
                    FP_Coeff_Vector_B,           &
                    FP_Source_Vector_B

USE Functions_Jacobian, &
            ONLY :  Calc_Ahat

USE Maps_Fixed_Point, &
            ONLY :  FP_Array_Map_TypeB

USE Maps_Domain, &
            ONLY :  Map_To_FEM_Node,            &
                    FEM_Elem_Map

USE Maps_Quadrature, &
            ONLY :  Map_To_tpd


USE XCFC_Functions_Calc_Values_Module, &
            ONLY :  Calc_Int_Weights,               &
                    Calc_Val_On_Elem_TypeB,         &
                    Calc_Val_And_Drv_On_Elem_TypeA, &
                    Calc_Val_And_Drv_On_Elem_TypeB

USE XCFC_Functions_Physical_Source_Module, &
            ONLY :  Get_Physical_Source

USE XCFC_Source_Variables_Module, &
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
                    Cur_Drv_Psi,        &
                    Cur_Val_AlphaPsi,   &
                    Cur_Drv_AlphaPsi,   &
                    Cur_Val_X,          &
                    Cur_Drv_X,          &
                    SourceTerm
                    
USE Variables_MPI, &
            ONLY :  myID_Poseidon,      &
                    nPROCS_Poseidon,    &
                    Poseidon_Comm_World

USE Poseidon_MPI_Utilities_Module, &
            ONLY :  STOP_MPI,               &
                    MPI_Master_Print,       &
                    MPI_All_Print

USE Initialization_Tables, &
            ONLY :  Initialize_Normed_Legendre_Tables_On_Level,     &
                    Initialize_Ylm_Tables_On_Elem


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
            ONLY :  MF_Source,              &
                    AMReX_Num_Levels

USE Variables_AMReX_Source, &
            ONLY :  Source_PTR,             &
                    Mask_PTR,               &
                    iLeaf,                &
                    iTrunk

USE Poseidon_AMReX_MakeFineMask_Module, &
            ONLY :  AMReX_MakeFineMask

USE Poseidon_AMReX_Multilayer_Utilities_Module, &
            ONLY :  Find_Coarsest_Parent

#endif

USE MPI




IMPLICIT NONE

REAL(idp)           :: E_Mass

CONTAINS
!+101+###########################################################################!
!                                                                                !
!           XCFC_Calc_Source_Vector_TypeB                                        !
!                                                                                !
!################################################################################!
SUBROUTINE XCFC_Calc_Source_Vector_TypeB( iU, iVB, iEU, iEL )

INTEGER, INTENT(IN), DIMENSION(3)       ::  iU          ! Variable Reference Numbers
INTEGER, INTENT(IN)                     ::  iVB         ! Variable Array Reference Number
INTEGER, INTENT(IN), DIMENSION(3)       ::  iEU         ! Upper Element Triplet
INTEGER, INTENT(IN), DIMENSION(3)       ::  iEL         ! Lower Element Triplet




#ifdef POSEIDON_AMREX_FLAG

    CALL XCFC_AMReX_Calc_Source_Vector_TypeB( iU, iVB )

#else

    INTEGER, DIMENSION(3)           ::  iE
    INTEGER                         ::  re, te, pe

    FP_Source_Vector_B(:,iVB) = 0.0_idp
    DO re = iEL(1),iEU(1)
    DO te = iEL(2),iEU(2)
    DO pe = iEL(3),iEU(3)
        iE = [re,te,pe]
        CALL XCFC_Calc_Source_Vector_On_Element_TypeB( iU, iVB, iE )
    END DO ! pe
    END DO ! te
    END DO ! re

#endif

!IF ( iVB == iVB_X ) THEN
!    PRINT*,"*********************"
!    PRINT*,"E_Mass ",E_Mass
!    PRINT*,"*********************"
!END IF


END SUBROUTINE XCFC_Calc_Source_Vector_TypeB



!+102+###########################################################################!
!                                                                                !
!           XCFC_Calc_Source_Vector_TypeB                                        !
!                                                                                !
!################################################################################!
SUBROUTINE XCFC_Calc_Source_Vector_On_Element_TypeB( iU, iVB, iE, Level_Option )

INTEGER, INTENT(IN), DIMENSION(3)               ::  iU
INTEGER, INTENT(IN)                             ::  iVB
INTEGER, INTENT(IN), DIMENSION(3)               ::  iE
INTEGER, INTENT(IN), OPTIONAL                   ::  Level_Option

INTEGER                                         ::  rd, tpd, td, pd, i

INTEGER                                         ::  FEM_Elem
REAL(KIND = idp)                                ::  DROT,     &
                                                    DTOT

INTEGER                                         ::  Level
INTEGER                                         ::  iCE(3)
INTEGER                                         ::  iRE(3)
INTEGER                                         ::  iEOff(3)






IF (Present(Level_Option)) THEN
    Level = Level_Option
ELSE
    Level = 0
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

DROT = drlocs(FEM_Elem)/2.0_idp
DTOT = Level_dx(Level,2)/2.0_idp

CUR_R_LOCS(:) = DROT * (Int_R_Locations(:) + 1.0_idp) + rlocs(FEM_Elem)
CUR_T_LOCS(:) = DTOT * (Int_T_Locations(:) + 1.0_idp + iEOff(2)*2.0_idp)



#else

FEM_Elem = iE(1)
DROT = 0.5_idp * (rlocs(iE(1)+1) - rlocs(iE(1)))
DTOT = 0.5_idp * (tlocs(iE(2)+1) - tlocs(iE(2)))


CUR_R_LOCS(:) = DROT * (INT_R_LOCATIONS(:)+1.0_idp) + rlocs(iE(1))
CUR_T_LOCS(:) = DTOT * (INT_T_LOCATIONS(:)+1.0_idp) + tlocs(iE(2))


#endif





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


CALL Calc_XCFC_CurVals_TypeB( iU, iVB,          &
                              iE,               &
                              DROT, DTOT,       &
                              Level             )




CALL Create_XCFC_Vector_TypeB( iE, iU, iVB, Level, FEM_Elem )


END SUBROUTINE XCFC_Calc_Source_Vector_On_Element_TypeB





!+201+##########################################################################!
!                                                                               !
!                  Create_XCFC_Vector_TypeB                                     !
!                                                                               !
!###############################################################################!
SUBROUTINE Create_XCFC_Vector_TypeB( iE, iU, iVB, Level, FEM_Elem )

INTEGER, INTENT(IN), DIMENSION(3)                   ::  iE
INTEGER, INTENT(IN), DIMENSION(3)                   ::  iU
INTEGER, INTENT(IN)                                 ::  iVB
INTEGER, INTENT(IN)                                 ::  Level
INTEGER, INTENT(IN)                                 ::  FEM_Elem

INTEGER                                             ::  ui, rd, d, lm_loc
INTEGER                                             ::  Current_i_Location

COMPLEX(KIND = idp)                                 ::  RHS_TMP
INTEGER                                             ::  iCT

REAL(idp)                                           ::  mass_tmp


! replace level with (level - iIRL)?
!iCT = 2**(level+1) - mod(iE(1),2**level) - 2
iCT = 0


DO ui = iU(1),iU(3)
DO lm_loc = 1,LM_LENGTH
DO d = 0,DEGREE


    RHS_TMP = 0.0_idp
    Mass_TMP = 0.0_idp

    DO rd = 1,NUM_R_QUAD_POINTS


#ifdef POSEIDON_AMREX_FLAG

        RHS_TMP = RHS_TMP                                       &
                + SUM( SourceTerm( :, rd, ui )                  &
                    * Ylm_Elem_CC_Values( :, lm_loc )           &
                    * TP_Int_Weights(:)                     )   &
                * Lagpoly_MultiLayer_Table( d, rd, 0, iCT )     &
                * R_Int_Weights(rd)


        MASS_TMP = MASS_TMP                                     &
                + SUM( SourceTerm( :, rd, ui )                  &
                    * Ylm_Elem_CC_Values( :, lm_loc )           &
                    * TP_Int_Weights(:)                     )   &
                * R_Int_Weights(rd)


!        IF ( ui == iU(1) ) THEN
!            PRINT*,level,iE,lm_loc,d,rd
!
!            PRINT*,SourceTerm( :, rd, ui )
!            PRINT*,"++++++++++++++++++++++"
!            PRINT*,Ylm_Elem_CC_Values( :, lm_loc )
!            PRINT*,"======================"
!            PRINT*,TP_Int_Weights(:)
!            PRINT*,"~~~~~~~~~~~~~~~~~~~~~~"

!            PRINT*,SUM( SourceTerm( :, rd, ui )             &
!            * Ylm_Elem_CC_Values( :, lm_loc )               &
!            * TP_Int_Weights(:)                     ),      &
!            Lagpoly_MultiLayer_Table( d, rd, 0, iCT ),      &
!            R_Int_Weights(rd)
!        END IF

#else


        RHS_TMP =  RHS_TMP                                          &
                + SUM( SourceTerm( :, rd, ui )                      &
                       * Ylm_CC_Values( :, lm_loc, iE(2), iE(3))    &
                       * TP_Int_Weights(:)                     )    &
                * Lagrange_Poly_Table(d, rd, 0)                     &
                * R_Int_Weights(rd)


        MASS_TMP = MASS_TMP                                         &
                + SUM( SourceTerm( :, rd, ui )                      &
                       * Ylm_CC_Values( :, lm_loc, iE(2), iE(3))    &
                       * TP_Int_Weights(:)                     )    &
                * R_Int_Weights(rd)

!        IF ( ui == iU(1) ) THEN
!            PRINT*,Level,iE,lm_loc,d,rd
!
!            PRINT*,SourceTerm( :, rd, ui )
!            PRINT*,"++++++++++++++++++++++"
!            PRINT*,Ylm_CC_Values( :, lm_loc, iE(2), iE(3))
!            PRINT*,"======================"
!            PRINT*,TP_Int_Weights(:)
!            PRINT*,"~~~~~~~~~~~~~~~~~~~~~~"


!            PRINT*,SUM( SourceTerm( :, rd, ui )             &
!            * Ylm_CC_Values( :, lm_loc, iE(2), iE(3))       &
!            * TP_Int_Weights(:)                     ),      &
!            Lagrange_Poly_Table(d, rd, 0),                      &
!            R_Int_Weights(rd)
!        END IF



#endif


    END DO  ! rd Loop



    Current_i_Location = FP_Array_Map_TypeB(ui,iVB,     &
                                            FEM_Elem,   &
                                            d, lm_loc   )


    FP_Source_Vector_B(Current_i_Location,iVB)          &
        = FP_Source_Vector_B(Current_i_Location,iVB)    &
        + RHS_TMP


    E_Mass = E_Mass + Mass_TMP

!    PRINT*,Current_i_Location,FP_Source_Vector_B(Current_i_Location,iVB)
END DO  ! d Loop
END DO  ! lm_loc Loop
END DO  ! ui Loop





END SUBROUTINE Create_XCFC_Vector_TypeB










!+202+###########################################################################!
!                                                                                !
!                  Calc_Current_Values          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_XCFC_CurVals_TypeB( iU, iVB, iE,            &
                                    DROT, DTOT,             &
                                    Level                   )


INTEGER, INTENT(IN), DIMENSION(3)                               ::  iU
INTEGER, INTENT(IN)                                             ::  iVB
INTEGER, INTENT(IN), DIMENSION(3)                               ::  iE
REAL(KIND = idp), INTENT(IN)                                    ::  DROT, DTOT
INTEGER, INTENT(IN)                                             ::  Level


REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3,3) ::  Ahat_Array
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3 )  ::  n_Array
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3 )  ::  f
REAL(idp), DIMENSION(Num_TP_Quad_Points, Num_R_Quad_Points,3)   ::  PhysSrc

INTEGER                                                         ::  tpd, rd
INTEGER                                                         ::  ui, i

IF ( iVB == iVB_X ) THEN

    DO rd = 1,NUM_R_QUAD_POINTS
    DO tpd = 1,NUM_TP_QUAD_POINTS

        ! Calc Flat Metric Terms
        f(tpd,rd,1) = 1.0_idp
        f(tpd,rd,2) = R_Square(rd)
        f(tpd,rd,3) = R_Square(rd) * TP_SIN_SQUARE(tpd)
    END DO ! tpd
    END DO ! rd


!    PRINT*,Level, iE
!    PRINT*,"A1"
    DO ui = iU(1),iU(3)
        CALL Get_Physical_Source( PhysSrc(:,:,ui-5), ui-3, iE )
        
!        CALL MPI_Barrier(Poseidon_Comm_World, ierr)
!        DO i = 0,nPROCs_Poseidon-1
!            IF( myID_Poseidon == i ) THEN
!                PRINT*,"myID ",i
!                PRINT*,PhysSrc(:,:,ui-5)
!            END IF
!            CALL MPI_Barrier(Poseidon_Comm_World, ierr)
!        END DO
!        IF ( ui == iU(1) ) THEN
!        DO rd = 1,Num_R_Quad_Points
!            PRINT*,"f(:,:,ui-5)",ui
!            PRINT*,f(:,rd,ui-5)
!            PRINT*,"PhysSrc(:,:,ui-5)"
!            PRINT*,PhysSrc(:,rd,ui-5)
!        END DO
!        END IF
!        PRINT*,"B1"
!        PRINT*,pi
!        PRINT*,GR_Source_Scalar
!        PRINT*,f(:,:,ui-5)
!        PRINT*," "
!        PRINT*,PhysSrc(:,:,ui-5)
!        PRINT*,""
!        PRINT*, 8.0_idp * pi * GR_Source_Scalar * f(:,:,ui-5)*PhysSrc(:,:,ui-5)
        SourceTerm(:,:,ui) = 8.0_idp * pi * GR_Source_Scalar * f(:,:,ui-5)*PhysSrc(:,:,ui-5)

    END DO
!    PRINT*,Level, iE
!    PRINT*,PhysSrc(:,:,iU(1)-5)
   

ELSE IF ( iVB == iVB_S ) THEN


    DO ui = iU(1),iU(3)
        CALL Get_Physical_Source( PhysSrc(:,:,ui-2), ui, iE )
    END DO

    CALL Calc_Val_And_Drv_On_Elem_TypeA( iE, DROT,    &
                                         Cur_Val_Psi,               &
                                         Cur_DRV_Psi,               &
                                         iU_CF,                     &
                                         Level                      )

    CALL Calc_Val_And_Drv_On_Elem_TypeA( iE, DROT,    &
                                         Cur_Val_AlphaPsi,          &
                                         Cur_DRV_AlphaPsi,          &
                                         iU_LF,                     &
                                         Level                      )
     
    CALL Calc_Val_And_Drv_On_Elem_TypeB( iE, DROT,    &
                                         CUR_Val_X(:,:,1),          &
                                         CUR_DRV_X(:,:,:,1),        &
                                         iU_X1, iVB_X ,             &
                                         Level                      )

    CALL Calc_Val_And_Drv_On_Elem_TypeB( iE, DROT,    &
                                         CUR_Val_X(:,:,2),          &
                                         CUR_DRV_X(:,:,:,2),        &
                                         iU_X2, iVB_X,              &
                                         Level                      )

    CALL Calc_Val_And_Drv_On_Elem_TypeB( iE, DROT,    &
                                         CUR_Val_X(:,:,3),          &
                                         CUR_DRV_X(:,:,:,3),        &
                                         iU_X3, iVB_X,              &
                                         Level                      )


    CALL Calc_Ahat( Ahat_Array,                             &
                    NUM_R_QUAD_POINTS, NUM_TP_QUAD_POINTS,  &
                    Cur_R_Locs, R_SQUARE,                   &
                    TP_RSIN_SQUARE, TP_COTAN_VAL,           &
                    CUR_VAL_X, CUR_DRV_X                  )



    DO i = 1,3
        n_Array(:,:,i) = Cur_DRV_AlphaPsi(:,:,i)/(Cur_VAL_Psi(:,:)**7)     &
                       - 7.0_idp * Cur_DRV_Psi(:,:,i)                      &
                         * Cur_Val_AlphaPsi(:,:)/(Cur_Val_psi(:,:)**8)
    END DO



    DO ui = 1,3

    SourceTerm(:,:,iU(ui)) = 16.0_idp * pi * GR_Source_Scalar               &   ! Physical Source
                        * Cur_Val_AlphaPsi(:,:)/(Cur_Val_Psi(:,:)**7)   &
                        * PhysSrc(:,:,ui)                               &
                      + 2.0_idp                                         &   ! Geometry Source
                        * ( N_Array(:,:,1)*Ahat_Array(:,:,ui,1)       &
                            + N_Array(:,:,2)*Ahat_Array(:,:,ui,2)     &
                            + N_Array(:,:,3)*Ahat_Array(:,:,ui,3)     )

    END DO


!    IF ( .TRUE. ) THEN
!    DO rd = 1,Num_R_Quad_Points
!        ui = 1
!        PRINT*,"Cur_Val_AlphaPsi",iE
!        PRINT*,Cur_Val_AlphaPsi(:,rd)
!        PRINT*,"Cur_Val_Psi"
!        pRINT*,Cur_Val_Psi(:,rd)
!        PRINT*,"PhysSrc(:,:,ui)"
!        PRINT*,PhysSrc(:,rd,ui)
!        PRINT*,"SourceTerm"
!        PRINT*,SourceTerm(:,rd,iU(ui))
!    END DO
!    END IF


ELSE

    WRITE(*,'(A)')" Error in Poseidon Subroutine : Calc_XCFC_CurVals_TypeB"
    WRITE(*,'(A)')" Invalid input. "
    WRITE(*,'(A)')" Variable : iVB "
    WRITE(*,'(A,I2.2)')" Value    : ",iVB

END IF



END SUBROUTINE Calc_XCFC_CurVals_TypeB







!+101+##########################################################################!
!                                                                               !
!          XCFC_AMReX_Calc_Source_Vector_TypeB                                  !
!                                                                               !
!###############################################################################!
SUBROUTINE XCFC_AMReX_Calc_Source_Vector_TypeB( iU, iVB )

INTEGER, INTENT(IN), DIMENSION(3)               ::  iU
INTEGER, INTENT(IN)                             ::  iVB

#ifdef POSEIDON_AMREX_FLAG

TYPE(amrex_mfiter)                              ::  mfi
TYPE(amrex_box)                                 ::  Box

TYPE(amrex_imultifab)                           ::  Level_Mask

INTEGER                                         ::  re, te, pe
INTEGER, DIMENSION(3)                           ::  iE
INTEGER, DIMENSION(3)                           ::  iEL, iEU
INTEGER                                         ::  nComp
INTEGER                                         ::  lvl




E_Mass = 0.0_idp

FP_Source_Vector_B(:,iVB) = 0.0_idp
DO lvl = AMReX_Num_Levels-1,0,-1


    !
    !   MakeFineMask
    !
    IF ( lvl < AMReX_Num_Levels-1 ) THEN
        CALL AMReX_MakeFineMask(  Level_Mask,               &
                                  MF_Source(lvl)%ba,        &
                                  MF_Source(lvl)%dm,        &
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


!    PRINT*,"After MakeFineMask"

    !
    !   Build mfiter
    !
!    PRINT*,"Before mfiter build"
    CALL amrex_mfiter_build(mfi, MF_Source(lvl), tiling = .true. )

!    PRINT*,"Before DO WHile"
    DO WHILE(mfi%next())

!        PRINT*,"Before PTRs"
        Source_PTR => MF_Source(lvl)%dataPtr(mfi)
        Mask_PTR   => Level_Mask%dataPtr(mfi)

!        PRINT*,"Before Box and nComp"
        Box = mfi%tilebox()
        nComp =  MF_Source(lvl)%ncomp()

!        PRINT*,"Before lo hi"
        iEL = Box%lo
        iEU = Box%hi

!        PRINT*,"Before Init_Normed_Legendre_Tables"
        !
        !   Initialize Legendre Polynomials on Box
        !
        CALL Initialize_Normed_Legendre_Tables_on_Level( iEU, iEL, lvl )

!        PRINT*,"After Init_Normed_Legendre_Tables"
        DO re = iEL(1),iEU(1)
        DO te = iEL(2),iEU(2)
        DO pe = iEL(3),iEU(3)

            IF ( Mask_PTR(RE,TE,PE,1) == iLeaf ) THEN

                !
                ! Initalize Ylm Table on Elem
                !
!                PRINT*,"Before Initialize_Ylm_Tables_on_Elem",lvl,re,te,pe
                CALL Initialize_Ylm_Tables_on_Elem( te, pe, iEL, lvl )

                iE = [re,te,pe]
!                PRINT*,"Before XCFC_Calc_Source_Vector_On_Element_TypeB"
                CALL XCFC_Calc_Source_Vector_On_Element_TypeB( iU, iVB, iE, lvl )
!                PRINT*,"After XCFC_Calc_Source_Vector_On_Element_TypeB"
            END IF
        END DO ! pe
        END DO ! te
        END DO ! re

    END DO

!    pRINT*,"Before Destroys"
    CALL amrex_mfiter_destroy(mfi)
!    PRINT*,"Between Destroys"
    CALL amrex_imultifab_destroy( Level_Mask )
!    PRINT*,"After Destroys"

END DO ! lvl


#endif

!PRINT*,"E_Mass",E_Mass

END SUBROUTINE XCFC_AMReX_Calc_Source_Vector_TypeB





END MODULE XCFC_Source_Vector_TypeB_Module
