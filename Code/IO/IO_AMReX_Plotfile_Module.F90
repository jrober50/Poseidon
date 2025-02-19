   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE IO_AMReX_Plotfile_Module                                              !##!
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
            
USE Poseidon_Numbers_Module, &
            ONLY :  pi
    
USE Poseidon_Units_Module, &
            ONLY :  C_Square,           &
                    Centimeter
            
USE Poseidon_Parameters, &
            ONLY :  Degree,                     &
                    L_Limit
                    
USE Variables_Mesh, &
            ONLY :  R_Inner,                &
                    R_Outer

USE Variables_Derived, &
            ONLY :  LM_Length,                  &
                    LM_Short_Length

USE Poseidon_IO_Parameters, &
            ONLY :  CFA_ShortVars,          &
                    Kij_ShortVars
            
USE Parameters_Variable_Indices, &
            ONLY :  iU_CF,      &
                    iU_LF,      &
                    iU_S1,      &
                    iU_S2,      &
                    iU_S3,      &
                    iU_X1,      &
                    iU_X2,      &
                    iU_X3,      &
                    iVB_S,      &
                    iVB_X

USE Variables_IO, &
            ONLY :  File_Suffix
            
USE Poseidon_IO_Parameters, &
            ONLY :  Poseidon_Output_Dir
            
USE Variables_FEM_Module, &
            ONLY :  FEM_Node_xlocs

USE Variables_Tables, &
            ONLY :  Level_dx
                    
USE Variables_Interface, &
            ONLY :  Caller_NQ,                  &
                    Caller_RQ_xlocs,            &
                    Caller_TQ_xlocs,            &
                    Caller_PQ_xlocs,            &
                    Caller_XL,                  &
                    Caller_Quad_DOF
                    
USE Variables_Quadrature, &
            ONLY :  Local_Quad_DOF,             &
                    NUM_R_QUAD_POINTS,          &
                    NUM_T_QUAD_POINTS,          &
                    NUM_P_QUAD_POINTS,          &
                    NUM_TP_QUAD_POINTS,         &
                    INT_R_LOCATIONS,            &
                    INT_T_LOCATIONS,            &
                    INT_P_LOCATIONS,            &
                    INT_R_WEIGHTS,              &
                    INT_T_WEIGHTS,              &
                    INT_P_WEIGHTS,              &
                    INT_TP_WEIGHTS,             &
                    xLeftLimit,                 &
                    xRightLimit
                    
USE Variables_Vectors, &
            ONLY :  dVA_Coeff_Vector,          &
                    dVB_Coeff_Vector
            
USE Maps_Fixed_Point, &
            ONLY :  FP_Array_Map_TypeB
            
USE Functions_Math, &
            ONLY :  Lagrange_Poly
            
USE Maps_Domain, &
            ONLY :  Map_To_FEM_Node
            
USE Poseidon_Return_Routines_All, &
            ONLY :  Poseidon_Return_All_AMReX
            
USE Maps_Quadrature, &
            ONLY :  Quad_Map
            
!USE External_Yahil_Profile_Module, &
!            ONLY :  Yahil_Potential_Solution,   &
!                    Yahil_Potential_Solution_Sub
                    
#ifdef POSEIDON_AMREX_FLAG
use amrex_base_module

USE amrex_amrcore_module, &
            ONLY :  amrex_ref_ratio,        &
                    amrex_geom,             &
                    amrex_get_finest_level
            
USE amrex_amr_module, &
            ONLY :  amrex_geom

USE amrex_box_module,   &
            ONLY :  amrex_box

USE amrex_boxarray_module, &
            ONLY :  amrex_boxarray


USE amrex_multifab_module,  &
            ONLY :  amrex_multifab,         &
                    amrex_multifab_build,   &
                    amrex_imultifab_build
                    
USE amrex_multifabutil_module, &
            ONLY :  amrex_average_down
                    
USE amrex_string_module, &
            ONLY :  amrex_string,           &
                    amrex_string_build
    
USE amrex_plotfile_module, &
            ONLY :  amrex_write_plotfile
            
USE Variables_AMReX_Core, &
            ONLY :  MF_Source,                  &
                    AMReX_Num_Levels
            
USE Poseidon_AMReX_MakeFineMask_Module, &
            ONLY :  AMReX_MakeFineMask
            
USE Parameters_AMReX, &
            ONLY :  iLeaf,                &
                    iTrunk

#endif





IMPLICIT NONE

CONTAINS
#ifdef POSEIDON_AMREX_FLAG
 !+101+########################################################!
!                                                               !
!          Write_AMReX_Plotfile                                 !
!                                                               !
 !#############################################################!
SUBROUTINE Write_AMReX_Plotfile( WriteX_Option )

LOGICAL,        INTENT(IN),     OPTIONAL    ::  WriteX_Option

CHARACTER(500)                              ::  IO_Str
CHARACTER(200)                              ::  PlotFileName
INTEGER                                     ::  Level
LOGICAL                                     ::  WriteX

TYPE(amrex_multifab)                        ::  MF_plt(0:AMReX_Num_Levels-1)
TYPE(amrex_multifab)                        ::  MF_sol(0:AMReX_Num_Levels-1)
TYPE(amrex_imultifab)                       ::  MF_Mask

INTEGER                                     ::  NumVars
INTEGER                                     ::  NumVars_Base
INTEGER                                     ::  NumVars_Mesh
INTEGER                                     ::  NumVars_Sol
INTEGER                                     ::  NumVars_Src
INTEGER                                     ::  NumVars_Err

TYPE(amrex_string),     ALLOCATABLE         ::  VarNames(:)

INTEGER                                     ::  iOS, iErr, lvl
INTEGER                                     ::  NumSolVars
INTEGER,        DIMENSION(1:3)              ::  nGhost_Vec

REAL(idp), DIMENSION(2)                     ::  GlobalError


PlotFileName = Poseidon_Output_Dir//"AMReX.plt_"//TRIM(File_Suffix)

nGhost_Vec = 0

NumVars_Base = 1 ! MPI proc
NumVars_Mesh = 6 ! X1_C, X2_C, X3_C, dX1, dX2, dX3
NumVars_Sol  = 11 ! CF, LF, B1, B2, B3, K11, K12, K13, K22, K23, K33
NumVars_Src  = 5  ! E, S, S1, S2, S3
NumVars_Err  = 2  ! CF_Err, CF_Sol

WriteX = .FALSE.
IF( PRESENT( WriteX_Option ) )THEN
    IF ( WriteX_Option )THEN
        WriteX = .TRUE.
        NumVars_Sol = NumVars_Sol + 3
    END IF
END IF


NumVars = NumVars_Base+NumVars_Mesh+NumVars_Sol+NumVars_Src+NumVars_Err

GlobalError = 0.0_idp


!
!   Build VarNames
!


ALLOCATE( VarNames(NumVars) )

CALL amrex_string_build( VarNames( 1 ), 'MPIProcess' )
CALL amrex_string_build( VarNames( 2 ), 'X1_C' )
CALL amrex_string_build( VarNames( 3 ), 'X2_C' )
CALL amrex_string_build( VarNames( 4 ), 'X3_C' )
CALL amrex_string_build( VarNames( 5 ), 'dX1' )
CALL amrex_string_build( VarNames( 6 ), 'dX2' )
CALL amrex_string_build( VarNames( 7 ), 'dX3' )

CALL amrex_string_build( VarNames( 8 ), 'E' )
CALL amrex_string_build( VarNames( 9 ), 'S' )
CALL amrex_string_build( VarNames( 10 ), 'S1' )
CALL amrex_string_build( VarNames( 11 ), 'S2' )
CALL amrex_string_build( VarNames( 12 ), 'S3' )

CALL amrex_string_build( VarNames( 13 ), CFA_ShortVars(iU_CF) )
CALL amrex_string_build( VarNames( 14 ), CFA_ShortVars(iU_LF) )

CALL amrex_string_build( VarNames( 15 ), CFA_ShortVars(iU_S1) )
CALL amrex_string_build( VarNames( 16 ), CFA_ShortVars(iU_S2) )
CALL amrex_string_build( VarNames( 17 ), CFA_ShortVars(iU_S3) )

CALL amrex_string_build( VarNames( 18 ), Kij_ShortVars(1) )
CALL amrex_string_build( VarNames( 19 ), Kij_ShortVars(2) )
CALL amrex_string_build( VarNames( 20 ), Kij_ShortVars(3) )
CALL amrex_string_build( VarNames( 21 ), Kij_ShortVars(4) )
CALL amrex_string_build( VarNames( 22 ), Kij_ShortVars(5) )
CALL amrex_string_build( VarNames( 23 ), Kij_ShortVars(6) )

iOS = 0
IF ( WriteX ) THEN
    CALL amrex_string_build( VarNames( 24 ), CFA_ShortVars(iU_X1) )
    CALL amrex_string_build( VarNames( 25 ), CFA_ShortVars(iU_X2) )
    CALL amrex_string_build( VarNames( 26 ), CFA_ShortVars(iU_X3) )
    iOS = 3
END IF

iErr = 24+iOS
CALL amrex_string_build( VarNames( iErr ), 'CF_Err' )
CALL amrex_string_build( VarNames( iErr+1 ), 'CF_Sol' )




Call BuildSolutionMultfab(MF_Sol, NumSolVars, WriteX )



!
!   Build multifab containing data to be output.
!
DO lvl = 0, AMReX_Num_Levels-1

    !
    !   Build and Initialize MF_plt
    !
    CALL amrex_multifab_build( MF_plt(lvl),             &
                               MF_Source(lvl) % BA,     &
                               MF_Source(lvl) % DM,     &
                               NumVars,                 &
                               0                        )
                               
    CALL MF_plt(lvl) % setVal( 0.0_idp )

    !
    !   MakeFineMask
    !
    IF ( lvl < AMReX_Num_Levels-1 ) THEN
        CALL AMReX_MakeFineMask(  MF_Mask,              &
                                  MF_Source(lvl)%ba,    &
                                  MF_Source(lvl)%dm,    &
                                  nGhost_Vec,           &
                                  MF_Source(lvl+1)%ba,  &
                                  iLeaf, iTrunk         )
                                  
    ELSE
            
        CALL amrex_imultifab_build( MF_Mask,            &
                                    MF_Source(lvl)%ba,  &
                                    MF_Source(lvl)%dm,  &
                                    1,                  &  ! ncomp = 1
                                    nGhost_Vec(1)       )
        CALL MF_Mask%SetVal(iLeaf)
    END IF



    !
    !   Fill MF_plt with data
    !
    CALL Fill_MPI( MF_plt(lvl) )

    CALL Fill_Mesh( lvl, MF_plt(lvl) )

    CALL FillSourceVariables( lvl, MF_plt(lvl), MF_Source(lvl), MF_Mask )

    CALL FillSolutionVariables( lvl, WriteX, NumVars_Sol, MF_plt(lvl), MF_Sol(lvl), MF_Mask )
    
!    CALL FillErrorVariables( lvl, NumVars_Sol, MF_plt(lvl), MF_Sol(lvl), MF_Mask, iErr, GlobalError )

END DO

CALL AverageDown(MF_plt)

CALL amrex_write_plotfile(  PlotFileName,       &
                            AMReX_Num_Levels,   &
                            MF_plt,             &
                            VarNames,           &
                            amrex_geom,         &
                            0.0_idp,            &
                            [1],                &
                            amrex_ref_ratio     )

print*,"Global L2 Error: ",GlobalError(1),GlobalError(2),GlobalError(1)/GlobalError(2)
   
   
DO lvl = 0, AMReX_Num_Levels-1
    CALL amrex_multifab_destroy ( MF_plt(lvl) )
END DO


END SUBROUTINE Write_AMReX_Plotfile







 !+201+########################################################!
!                                                               !
!          Fill_MPI                                             !
!                                                               !
 !#############################################################!
SUBROUTINE Fill_MPI( MF_plt )

TYPE(amrex_multifab),   INTENT(INOUT)       ::  MF_plt

INTEGER                                     ::  iX1, iX2, iX3
INTEGER                                     ::  iEL(3), iEU(3)
TYPE(amrex_box)                             ::  Box
TYPE(amrex_mfiter)                          ::  mfi

REAL(idp),           CONTIGUOUS,     POINTER :: U_plt(:,:,:,:)



CALL amrex_mfiter_build( mfi, MF_plt, tiling = .FALSE. )

DO WHILE( mfi % next() )

    U_plt => MF_plt%DataPtr( mfi )

    Box = mfi%TileBox()

    iEL = Box % lo
    iEU = Box % hi

    DO iX3 = iEL(3), iEU(3)
    DO iX2 = iEL(2), iEU(2)
    DO iX1 = iEL(1), iEU(1)

        U_plt(iX1,iX2,iX3,1) = amrex_parallel_myproc()

    END DO
    END DO
    END DO

END DO

CALL amrex_mfiter_destroy( mfi )

END SUBROUTINE Fill_MPI




 !+202+########################################################!
!                                                               !
!          Fill_Mesh                                            !
!                                                               !
 !#############################################################!
SUBROUTINE Fill_Mesh( lvl, MF_plt )

INTEGER,                INTENT(IN)          ::  lvl
TYPE(amrex_multifab),   INTENT(INOUT)       ::  MF_plt

INTEGER                                     ::  iX1, iX2, iX3
INTEGER                                     ::  iEL(3), iEU(3)
INTEGER                                     ::  Offset(3)
TYPE(amrex_box)                             ::  Box
TYPE(amrex_mfiter)                          ::  mfi



REAL(idp),           CONTIGUOUS,     POINTER :: U_plt(:,:,:,:)




CALL amrex_mfiter_build( mfi, MF_plt, tiling = .FALSE. )

DO WHILE( mfi % next() )

    U_plt => MF_plt%DataPtr( mfi )

    Box = mfi%TileBox()

    iEL = Box % lo
    iEU = Box % hi
    

    
    Offset = 0
    if (amrex_spacedim == 1 ) THEN
        Offset(2:3) = -1
    ELSE IF ( amrex_spacedim == 2 ) THEN
        Offset(3) = -1
    END IF

    DO iX3 = iEL(3), iEU(3)
    DO iX2 = iEL(2), iEU(2)
    DO iX1 = iEL(1), iEU(1)

        U_plt(iX1,iX2,iX3,2) = Level_dx(Lvl,1) * ( 0.5_idp + iX1 + Offset(1) )
        U_plt(iX1,iX2,iX3,3) = Level_dx(Lvl,2) * ( 0.5_idp + iX2 + Offset(2) )
        U_plt(iX1,iX2,iX3,4) = Level_dx(Lvl,3) * ( 0.5_idp + iX3 + Offset(3) )

        U_plt(iX1,iX2,iX3,5) = Level_dx(Lvl,1)
        U_plt(iX1,iX2,iX3,6) = Level_dx(Lvl,2)
        U_plt(iX1,iX2,iX3,7) = Level_dx(Lvl,3)

    END DO
    END DO
    END DO

END DO

CALL amrex_mfiter_destroy( mfi )

END SUBROUTINE Fill_Mesh









 !+301+########################################################!
!                                                               !
!          BuildSolutionMultfab                                 !
!                                                               !
 !#############################################################!
SUBROUTINE BuildSolutionMultfab( MF_Sol, NumSolVars, WriteX )

TYPE(amrex_multifab),       INTENT(OUT)     ::  MF_Sol(0:AMReX_Num_Levels-1)
INTEGER,                    INTENT(OUT)     ::  NumSolVars
LOGICAL,                    INTENT(IN)      ::  WriteX

INTEGER                                     ::  lvl


NumSolVars = 11 ! CF, LF, B1, B2, B3, K11, K12, K13, K22, K23, K33

IF ( WriteX ) THEN
    NumSolVars = 14 ! Above + X1, X2, X3
END IF

!
!   Build multifab containing data to be output.
!
DO lvl = 0, AMReX_Num_Levels-1

    !
    !   Build and Initialize MF_sol
    !
    CALL amrex_multifab_build( MF_Sol(lvl),                 &
                               MF_Source(lvl) % BA,         &
                               MF_Source(lvl) % DM,         &
                               NumSolVars*Local_Quad_DOF,   &
                               0                            )
                               
    CALL MF_Sol(lvl) % setVal( 0.0_idp )



END DO


CALL Poseidon_Return_All_AMReX([Num_R_Quad_Points, Num_T_Quad_Points, Num_P_Quad_Points],            &
                         Int_R_Locations,   &
                         Int_T_Locations,   &
                         Int_P_Locations,   &
                         xLeftLimit,        &
                         xRightLimit,       &
                         AMReX_Num_Levels,  &
                         MF_Sol,            &
                         ReturnX_Option = WriteX )



END SUBROUTINE BuildSolutionMultfab




 !+203+########################################################!
!                                                               !
!          FillSolutionVariables                                !
!                                                               !
 !#############################################################!
SUBROUTINE FillSolutionVariables( lvl, WriteX, NumVars, MF_plt, MF_Sol, MF_Mask )

INTEGER,                INTENT(IN)              ::  lvl
LOGICAL,                INTENT(IN)              ::  WriteX
INTEGER,                INTENT(IN)              ::  NumVars
TYPE(amrex_multifab),   INTENT(INOUT)           ::  MF_plt
TYPE(amrex_multifab),   INTENT(INOUT)           ::  MF_Sol
TYPE(amrex_imultifab),  INTENT(IN)              ::  MF_Mask


INTEGER                                         ::  iEL(3), iEU(3)
TYPE(amrex_box)                                 ::  Box
TYPE(amrex_mfiter)                              ::  mfi

REAL(idp),  DIMENSION(Local_Quad_DOF,NumVars)   ::  U
INTEGER,    DIMENSION(3)                        ::  nE
INTEGER                                         ::  iX1, iX2, iX3, i
INTEGER,    DIMENSION(4)                        ::  LoSol, HiSol

REAL(idp),           CONTIGUOUS,     POINTER    :: U_plt(:,:,:,:)
REAL(idp),           CONTIGUOUS,     POINTER    :: U_Sol(:,:,:,:)
INTEGER,             CONTIGUOUS,     POINTER    :: U_Mask(:,:,:,:)

REAL(idp)                                       ::  DROT,       &
                                                    DTOT,       &
                                                    DPOT
                                                    
INTEGER                                         ::  Here
INTEGER                                         ::  rd, td, pd
REAL(idp)                                       ::  Volume
INTEGER                                         ::  iEOff(3)
REAL(idp)                                       ::  CellAverage(NumVars)
REAL(idp)                                       ::  Cur_R_Locs(Num_R_Quad_Points)
REAL(idp)                                       ::  Cur_T_Locs(Num_T_Quad_Points)
REAL(idp)                                       ::  IntWeights(Local_Quad_DOF)


DROT = Level_dx(Lvl,1)/2.0_idp
DTOT = Level_dx(Lvl,2)/2.0_idp
DPOT = Level_dx(Lvl,3)/2.0_idp


CALL amrex_mfiter_build( mfi, MF_plt, tiling = .false. )

DO WHILE( mfi % next() )

    U_plt => MF_plt%DataPtr( mfi )
    U_Sol => MF_Sol%DataPtr( mfi )
    U_Mask => MF_Mask%DataPtr( mfi )

    Box = mfi%TileBox()

    iEL = Box % lo
    iEU = Box % hi
    
    LoSol = LBOUND( U_Sol )
    HiSol = UBOUND( U_Sol )
  
    DO iX3 = iEL(3), iEU(3)
    DO iX2 = iEL(2), iEU(2)
    DO iX1 = iEL(1), iEU(1)
    
        IF ( amrex_spacedim == 1 ) THEN
            iEoff(2:3) = 0
        ELSEIF ( amrex_spacedim == 2) THEN
            iEoff(2)   = iX2
            iEoff(3)   = 0
        ELSEIF ( amrex_spacedim == 3 ) THEN
            iEoff(2) = iX2
            iEoff(3) = iX3
        END IF
    
        Cur_R_Locs(:) = DROT * (Int_R_Locations(:) + 1.0_idp + iX1*2.0_idp)
        Cur_T_Locs(:) = DTOT * (Int_T_Locations(:) + 1.0_idp + iEOff(2)*2.0_idp)
        
                               
        DO rd = 1,NUM_R_QUAD_POINTS
        DO td = 1,NUM_T_QUAD_POINTS
        DO pd = 1,NUM_P_QUAD_POINTS
           Here = Quad_Map(rd,td,pd)
           IntWeights( Here ) = Cur_R_Locs(rd)              &
                              * Cur_R_Locs(rd)              &
                              * DSIN(CUR_T_LOCS(td))        &
                              * DROT * INT_R_WEIGHTS(rd)    &
                              * DTOT * INT_T_WEIGHTS(td)    &
                              * INT_P_WEIGHTS(pd)
        END DO
        END DO
        END DO
    
        Volume = SUM(IntWeights)
    
        U(1:Local_Quad_DOF,1:NumVars) = RESHAPE(U_sol(iX1,iX2,iX3,LoSol(4):HiSol(4)),[Local_Quad_DOF,NumVars] )
    
        
        DO i = 1,NumVars
            
            CellAverage(i) = SUM(U(:,i)*IntWeights)/Volume
            
        END DO
        
!        if (iX2 == 7) THEN
!            print*,iX1,iX2,CellAverage(1),SUM(U(:,1))/Local_Quad_DOF,U(1,1)
!
!        END IF

!        DO rd = 1,NUM_R_QUAD_POINTS
!        DO td = 1,NUM_T_QUAD_POINTS
!        DO pd = 1,NUM_P_QUAD_POINTS
!           Here = Quad_Map(rd,td,pd)
!           print*,lvl,iX1,iX2,Cur_r_Locs(rd),Cur_T_Locs(td),U(Here,1)
!        END DO
!        END DO
!        END DO
        
        U_plt(iX1,iX2,iX3,13:13+NumVars-1) = CellAverage
    END DO
    END DO
    END DO


END DO

CALL amrex_mfiter_destroy( mfi )

END SUBROUTINE FillSolutionVariables





 !+204+########################################################!
!                                                               !
!          FillSourceVariables                                  !
!                                                               !
 !#############################################################!
SUBROUTINE FillSourceVariables( lvl, MF_plt, MF_Src, MF_Mask )

INTEGER,                INTENT(IN)                  ::  lvl
TYPE(amrex_multifab),   INTENT(INOUT)               ::  MF_plt
TYPE(amrex_multifab),   INTENT(IN)                  ::  MF_Src
TYPE(amrex_imultifab),  INTENT(IN)                  ::  MF_Mask

INTEGER, PARAMETER                                  ::  NumVars_Src = 5

TYPE(amrex_box)                                     ::  Box
TYPE(amrex_mfiter)                                  ::  mfi
    
INTEGER,    DIMENSION(3)                            ::  nE
INTEGER                                             ::  iX1, iX2, iX3, i
INTEGER                                             ::  iEL(3), iEU(3)
INTEGER                                             ::  iEOff(3)
INTEGER,    DIMENSION(4)                            ::  LoSrc, HiSrc

REAL(idp),           CONTIGUOUS,     POINTER        :: U_plt(:,:,:,:)
REAL(idp),           CONTIGUOUS,     POINTER        :: U_Src(:,:,:,:)
INTEGER,             CONTIGUOUS,     POINTER        :: U_Mask(:,:,:,:)

REAL(idp)                                           ::  DROT,       &
                                                        DTOT,       &
                                                        DPOT
                                                                                                        
INTEGER                                             ::  Here
INTEGER                                             ::  rd, td, pd
REAL(idp)                                           ::  Cur_R_Locs(Num_R_Quad_Points)
REAL(idp)                                           ::  Cur_T_Locs(Num_T_Quad_Points)
REAL(idp)                                           ::  IntWeights(Local_Quad_DOF)

REAL(idp)                                           ::  Volume
REAL(idp),  DIMENSION(NumVars_Src)                  ::  CellAverage(NumVars_Src)
REAL(idp),  DIMENSION(Local_Quad_DOF,NumVars_Src)   ::  U
    

DROT = Level_dx(Lvl,1)/2.0_idp
DTOT = Level_dx(Lvl,2)/2.0_idp
DPOT = Level_dx(Lvl,3)/2.0_idp


CALL amrex_mfiter_build( mfi, MF_plt, tiling = .FALSE. )

DO WHILE( mfi % next() )

    U_plt => MF_plt%DataPtr( mfi )
    U_Src => MF_Src%DataPtr( mfi )
    U_Mask => MF_Mask%DataPtr( mfi )
    
    Box = mfi%TileBox()

    iEL = Box % lo
    iEU = Box % hi
    
    LoSrc = LBOUND( U_Src )
    HiSrc = UBOUND( U_Src )
  
    DO iX3 = iEL(3), iEU(3)
    DO iX2 = iEL(2), iEU(2)
    DO iX1 = iEL(1), iEU(1)
    
        IF ( amrex_spacedim == 1 ) THEN
            iEoff(2:3) = 0
        ELSEIF ( amrex_spacedim == 2) THEN
            iEoff(2)   = iX2
            iEoff(3)   = 0
        ELSEIF ( amrex_spacedim == 3 ) THEN
            iEoff(2) = iX2
            iEoff(3) = iX3
        END IF
    
        Cur_R_Locs(:) = DROT * (Int_R_Locations(:) + 1.0_idp + iX1*2.0_idp)
        Cur_T_Locs(:) = DTOT * (Int_T_Locations(:) + 1.0_idp + iEOff(2)*2.0_idp)
        
                               
        DO rd = 1,NUM_R_QUAD_POINTS
        DO td = 1,NUM_T_QUAD_POINTS
        DO pd = 1,NUM_P_QUAD_POINTS
           Here = Quad_Map(rd,td,pd)
           IntWeights( Here ) = Cur_R_Locs(rd)              &
                              * Cur_R_Locs(rd)              &
                              * DSIN(CUR_T_LOCS(td))        &
                              * DROT * INT_R_WEIGHTS(rd)    &
                              * DTOT * INT_T_WEIGHTS(td)    &
                              * INT_P_WEIGHTS(pd)
        END DO
        END DO
        END DO
    
        Volume = SUM(IntWeights)
    
        U(1:Local_Quad_DOF,1:NumVars_Src) = RESHAPE(U_Src(iX1,iX2,iX3,LoSrc(4):HiSrc(4)),[Local_Quad_DOF,NumVars_Src] )
        
        DO i = 1,NumVars_Src
            
            CellAverage(i) = SUM(U(:,i)*IntWeights)/Volume
            
        END DO

        U_plt(iX1,iX2,iX3,8:8+NumVars_Src-1) = CellAverage

    END DO
    END DO
    END DO

END DO

CALL amrex_mfiter_destroy( mfi )

END SUBROUTINE FillSourceVariables







!+203+########################################################!
!                                                               !
!          FillSolutionVariables                                !
!                                                               !
 !#############################################################!
SUBROUTINE FillErrorVariables( lvl, NumVars, MF_plt, MF_Sol, MF_Mask, iErr, GlobalError )

INTEGER,                INTENT(IN)              ::  lvl
INTEGER,                INTENT(IN)              ::  NumVars
TYPE(amrex_multifab),   INTENT(INOUT)           ::  MF_plt
TYPE(amrex_multifab),   INTENT(INOUT)           ::  MF_Sol
TYPE(amrex_imultifab),  INTENT(IN)              ::  MF_Mask
REAL(idp),              INTENT(INOUT)           ::  GlobalError(2)

INTEGER,                INTENT(IN)              ::  iErr

INTEGER                                         ::  iEL(3), iEU(3)
TYPE(amrex_box)                                 ::  Box
TYPE(amrex_mfiter)                              ::  mfi

REAL(idp),  DIMENSION(Local_Quad_DOF,NumVars)   ::  U
INTEGER,    DIMENSION(3)                        ::  nE
INTEGER                                         ::  iX1, iX2, iX3, i
INTEGER,    DIMENSION(4)                        ::  LoSol, HiSol

REAL(idp),           CONTIGUOUS,     POINTER    :: U_plt(:,:,:,:)
REAL(idp),           CONTIGUOUS,     POINTER    :: U_Sol(:,:,:,:)
INTEGER,             CONTIGUOUS,     POINTER    :: U_Mask(:,:,:,:)

REAL(idp)                                       ::  DROT,       &
                                                    DTOT,       &
                                                    DPOT
                                                    
INTEGER                                         ::  Here
INTEGER                                         ::  rd, td, pd
REAL(idp)                                       ::  Volume
REAL(idp)                                       ::  ErrInt
INTEGER                                         ::  iEOff(3)
REAL(idp)                                       ::  Cur_R_Locs(Num_R_Quad_Points)
REAL(idp)                                       ::  Cur_T_Locs(Num_T_Quad_Points)
REAL(idp)                                       ::  Cur_P_Locs(Num_P_Quad_Points)
REAL(idp)                                       ::  IntWeights(Local_Quad_DOF)
REAL(idp)                                       ::  PsiSol(Local_Quad_DOF)
REAL(idp)                                       ::  Error(Local_Quad_DOF)
REAL(idp)                                       ::  Potential
REAL(idp)                                       ::  CellAverage

DROT = Level_dx(Lvl,1)/2.0_idp
DTOT = Level_dx(Lvl,2)/2.0_idp
DPOT = Level_dx(Lvl,3)/2.0_idp


CALL amrex_mfiter_build( mfi, MF_plt, tiling = .false. )

DO WHILE( mfi % next() )

    U_plt => MF_plt%DataPtr( mfi )
    U_Sol => MF_Sol%DataPtr( mfi )
    U_Mask => MF_Mask%DataPtr( mfi )

    Box = mfi%TileBox()

    iEL = Box % lo
    iEU = Box % hi
    
    LoSol = LBOUND( U_Sol )
    HiSol = UBOUND( U_Sol )
  
    DO iX3 = iEL(3), iEU(3)
    DO iX2 = iEL(2), iEU(2)
    DO iX1 = iEL(1), iEU(1)
    
        IF ( amrex_spacedim == 1 ) THEN
            iEoff(2:3) = 0
        ELSEIF ( amrex_spacedim == 2) THEN
            iEoff(2)   = iX2
            iEoff(3)   = 0
        ELSEIF ( amrex_spacedim == 3 ) THEN
            iEoff(2) = iX2
            iEoff(3) = iX3
        END IF
    
        Cur_R_Locs(:) = DROT * (Int_R_Locations(:) + 1.0_idp + iX1*2.0_idp)
        Cur_T_Locs(:) = DTOT * (Int_T_Locations(:) + 1.0_idp + iEOff(2)*2.0_idp)
        Cur_P_Locs(:) = DPOT * (Int_P_Locations(:) + 1.0_idp + iEOff(3)*2.0_idp)
                               
        DO rd = 1,NUM_R_QUAD_POINTS
        DO td = 1,NUM_T_QUAD_POINTS
        DO pd = 1,NUM_P_QUAD_POINTS
           Here = Quad_Map(rd,td,pd)
           IntWeights( Here ) = Cur_R_Locs(rd)              &
                              * Cur_R_Locs(rd)              &
                              * DSIN(CUR_T_LOCS(td))        &
                              * DROT * INT_R_WEIGHTS(rd)    &
                              * DTOT * INT_T_WEIGHTS(td)    &
                              * INT_P_WEIGHTS(pd)
        END DO
        END DO
        END DO
        
        DO rd = 1,NUM_R_QUAD_POINTS
        DO td = 1,NUM_T_QUAD_POINTS
        DO pd = 1,NUM_P_QUAD_POINTS
           Here = Quad_Map(rd,td,pd)
!           Potential = Yahil_Potential_Solution(Cur_R_Locs(rd)*Centimeter,Cur_T_Locs(td),Cur_P_Locs(pd))

            PsiSol(Here) = 1.0_idp - Potential/(2.0_idp*C_Square)
        END DO
        END DO
        END DO
    
        U(1:Local_Quad_DOF,1:NumVars) = RESHAPE(U_sol(iX1,iX2,iX3,LoSol(4):HiSol(4)),[Local_Quad_DOF,NumVars] )
        Error = ABS(ABS(U(:,1) - PsiSol(:)))/ABS(PsiSol(:))
        
        
        Volume = SQRT(SUM(PsiSol*IntWeights)*SUM(PsiSol*IntWeights))
        ErrInt = SQRT(SUM(Error*IntWeights)*SUM(Error*IntWeights))
        
        CellAverage = ErrInt/Volume
        
        
        IF ( U_Mask(iX1,iX2,iX3,1) == iLeaf ) THEN
            GlobalError(1) = GlobalError(1) + ErrInt
            GlobalError(2) = GlobalError(2) + Volume
         END IF

        
        U_plt(iX1,iX2,iX3,iErr) = ErrInt/Volume
        U_plt(iX1,iX2,iX3,iErr+1)= SUM(PsiSol*IntWeights)/Volume
        
    END DO
    END DO
    END DO


END DO

CALL amrex_mfiter_destroy( mfi )

END SUBROUTINE FillErrorVariables














!+301+########################################################!
!                                                               !
!          AverageDown                                          !
!                                                               !
 !#############################################################!
SUBROUTINE AverageDown( MF )
TYPE(amrex_multifab),           INTENT(INOUT)       ::  MF(0:)


INTEGER                                             ::  lvl
INTEGER                                             ::  FinestLevel

FinestLevel = amrex_get_finest_level()

DO lvl = FinestLevel-1,0,-1

    CALL amrex_average_down( MF(lvl+1),             &
                             MF(lvl),               &
                             amrex_geom(lvl+1),     &
                             amrex_geom(lvl),       &
                             1,                 &
                             MF(lvl)%nComp(),       &
                             amrex_ref_ratio(lvl)   )

END DO


END SUBROUTINE AverageDown



#endif





END MODULE IO_AMReX_Plotfile_Module
