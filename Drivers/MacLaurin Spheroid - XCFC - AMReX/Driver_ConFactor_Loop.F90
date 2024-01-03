   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Driver_ConFactor_Loop_Module                                          !##!
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
use amrex_base_module

USE amrex_box_module,   &
            ONLY:   amrex_box
            
USE amrex_multifab_module,  &
            ONLY :  amrex_multifab,         &
                    amrex_multifab_build,   &
                    amrex_imultifab_build,  &
                    amrex_imultifab_destroy
                    
USE amrex_amrcore_module, &
            ONLY:   amrex_geom,             &
                    amrex_get_numlevels


USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag
            
USE Parameters_Variable_Indices, &
            ONLY :  iU_CF

USE Poseidon_Message_Routines_Module, &
            ONLY :  Driver_Init_Message

USE Poseidon_Interface_Source_Input, &
            ONLY :  Poseidon_Input_Sources

USE Variables_External, &
            ONLY :  HCT_Rhoo,       &
                    HCT_Star_Radius

USE Driver_SetSource_Module, &
            ONLY :  Driver_SetSource

USE Poseidon_Interface_Run, &
            ONLY :  Poseidon_Run

USE Poseidon_Interface_Return_Routines, &
            ONLY :  Poseidon_Return_Conformal_Factor

USE Maps_Quadrature, &
            ONLY :  Quad_Map

USE IO_Print_Results, &
            ONLY :  Print_Results

USE Driver_SetGuess_Module, &
            ONLY :  Driver_SetGuess
            
USE Variables_AMReX_Core, &
            ONLY :  AMReX_Num_Levels,       &
                    AMReX_Max_Grid_Size,    &
                    MF_Source,              &
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
                    
USE Variables_Interface, &
            ONLY :  Caller_Quad_DOF,        &
                    Caller_NQ,                      &
                    Caller_xL,                      &
                    Caller_RQ_xlocs,                &
                    Caller_TQ_xlocs,                &
                    Caller_PQ_xlocs
            
USE Driver_InitSource_Module, &
            ONLY :  Driver_InitSource
            
USE Variables_Driver_AMReX, &
            ONLY :  MF_Src_nComps,      &
                    MF_Driver_Source
                    
USE IO_Write_Final_Results, &
            ONLY :  Output_ConFactorLoopData

USE Variables_External, &
            ONLY :  CFLD_Update,            &
                    CFLD_Residual,          &
                    CFLD_Tolerance,         &
                    CFLD_MaxIters,          &
                    CFLD_Iters

USE Variables_Tables, &
            ONLY :  Level_DX
            
USE External_HCT_Solution_Module, &
            ONLY :  HCT_Solution

USE External_IO_Test_Results_Module, &
            ONLY :  Print_MacLaurin_Error

USE Flags_Source_Input_Module, &
            ONLY :  lPF_SI_Flags,       &
                    iPF_SI_MF_Ready


IMPLICIT NONE


INTEGER                                     ::  MF_nVars    = 11

INTEGER                                     ::  MF_nGhosts  = 0
INTEGER, DIMENSION(1:3)                     ::  nGhost_Vec



CONTAINS



 !+101+############################################################!
!                                                                   !
!          Driver_ConFactor_Loop                                    !
!                                                                   !
 !#################################################################!
SUBROUTINE Driver_ConFactor_Loop( )


TYPE(amrex_multifab),   ALLOCATABLE         ::  MF_Old(:)
TYPE(amrex_multifab),   ALLOCATABLE         ::  MF_New(:)



INTEGER                                     ::  MF_nComps

INTEGER                                     ::  Iter


REAL(idp)                                   ::  Max_Difference
REAL(idp)                                   ::  Error_Norms(3)

LOGICAL                                     ::  Flag
INTEGER                                     ::  lvl

IF ( Verbose_Flag ) CALL Driver_Init_Message('Begining the conformal factor loop.')

CFLD_MaxIters   = 500
CFLD_Tolerance  = 6E-15

MF_nComps       = MF_Src_nComps


ALLOCATE( CFLD_Residual(1:CFLD_MaxIters) )
ALLOCATE( CFLD_Update(1:CFLD_MaxIters) )
CFLD_Update   = 0.0_idp
CFLD_Residual = 0.0_idp

ALLOCATE( MF_Old(0:AMReX_Num_Levels-1) )
ALLOCATE( MF_New(0:AMReX_Num_Levels-1) )

DO lvl = 0,AMReX_Num_Levels-1

    nGhost_Vec = MF_Source(lvl)%nghostvect()

    CALL amrex_multifab_build(  MF_Old(lvl),        &
                                MF_Source(Lvl)%BA,  &
                                MF_Source(Lvl)%DM,  &
                                MF_nComps,          &
                                nGhost_Vec(1)       )

    CALL MF_Old(lvl)%SetVal(0.0_idp)

    CALL amrex_multifab_build(  MF_New(lvl),        &
                                MF_Source(Lvl)%BA,  &
                                MF_Source(Lvl)%DM,  &
                                MF_nComps,          &
                                nGhost_Vec(1)       )

    CALL MF_New(lvl)%SetVal(0.0_idp)

END DO



CALL Poseidon_Run()
    
CALL Poseidon_Return_Conformal_Factor( MF_Old )


Iter = 2
Flag = .TRUE.
DO WHILE ( Flag )

!    DEALLOCATE( MF_Driver_Source )
    lPF_SI_Flags(iPF_SI_MF_Ready) = .FALSE.
    CALL Driver_InitSource( )

    CALL Driver_SetSource( )

    CALL Poseidon_Run()

    PRINT*,"Here",Iter
    CALL Poseidon_Return_Conformal_Factor( MF_New )
    PRINT*,"There",Iter
    Max_Difference = Calculate_Max_Difference( MF_New, MF_Old)
    
    Error_Norms = Calculate_Error_Norms( MF_New )


    DO lvl = 0,AMReX_Num_Levels-1
        CALL MF_Old(lvl)%copy(  MF_New(lvl),        &
                                1,                  &
                                1,                  &
                                MF_nComps,          &
                                nGhost_Vec(1)       )
    END DO



    CFLD_Update(Iter)   = Max_Difference
    CFLD_Residual(Iter) = Error_Norms(2)
    
    
    IF ( Iter >= CFLD_MaxIters ) THEN
        Flag = .FALSE.
    ELSEIF ( Max_Difference .LE. CFLD_Tolerance) THEN
        Flag = .FALSE.
    ELSE
        Flag = .TRUE.
        Iter = Iter + 1
    END IF


    IF ( .TRUE. ) THEN
!    IF ( Verbose_Flag ) THEN
        WRITE(*,'(A,I3.3,A,ES22.15,A,ES22.15)')         &
                "End of Conformal Factor Iteration: ",Iter,  &
                 " - Maximum Difference: ",Max_Difference,            &
                 " Tolerance: ", CFLD_Tolerance
    END IF

END DO ! While Iter < Iter_Max


CFLD_Iters = Iter


CALL Output_ConFactorLoopData()

DEALLOCATE( MF_Driver_Source )



IF ( .TRUE. ) THEN
    CALL Print_MacLaurin_Error()
END IF

END SUBROUTINE Driver_ConFactor_Loop





 !+101+############################################################!
!                                                                   !
!        Calculate_Max_Difference                                   !
!                                                                   !
 !#################################################################!
FUNCTION Calculate_Max_Difference( MF_New, MF_Old  )

REAL(idp)                                   ::  Calculate_Max_Difference
TYPE(amrex_multifab),   INTENT(IN)          ::  MF_New(0:AMReX_Num_Levels-1)
TYPE(amrex_multifab),   INTENT(IN)          ::  MF_Old(0:AMReX_Num_Levels-1)

INTEGER,   CONTIGUOUS, POINTER              ::  Mask_PTR(:,:,:,:)
REAL(idp), CONTIGUOUS, POINTER              ::  New_PTR(:,:,:,:)
REAL(idp), CONTIGUOUS, POINTER              ::  Old_PTR(:,:,:,:)


INTEGER                                     ::  lvl
TYPE(amrex_mfiter)                          ::  mfi
TYPE(amrex_box)                             ::  Box
TYPE(amrex_imultifab)                       ::  Level_Mask
TYPE(amrex_imultifab)                       ::  Ghost_Mask

INTEGER,    DIMENSION(1:3)                  ::  iEL, iEU
INTEGER                                     ::  re, te, pe
INTEGER                                     ::  i
INTEGER                                     ::  Here

REAL(idp)                                   ::  Difference
REAL(idp)                                   ::  Max_Difference
REAL(idp)                                   ::  New_Value
REAL(idp)                                   ::  Old_Value



Max_Difference = 0.0_idp
DO lvl = 0,AMReX_Num_Levels-1


    nGhost_Vec = MF_New(lvl)%nghostvect()

    !
    !   MakeFineMask
    !
    IF ( lvl < AMReX_Num_Levels-1 ) THEN
        CALL AMReX_MakeFineMask(  Level_Mask,               &
                                  MF_New(lvl)%ba,           &
                                  MF_New(lvl)%dm,           &
                                  nGhost_Vec,               &
                                  MF_New(lvl+1)%ba,     &
                                  iLeaf, iTrunk             )
    ELSE
        ! Create Level_Mask all equal to 1
        CALL amrex_imultifab_build( Level_Mask,             &
                                    MF_New(lvl)%ba,         &
                                    MF_New(lvl)%dm,         &
                                    1,                      &
                                    nGhost_Vec(1)           )
        CALL Level_Mask%SetVal(iLeaf)
    END IF
    
    !
    !   BuildMask
    !
    CALL amrex_imultifab_build( Ghost_Mask,             &
                                MF_New(lvl)%ba,         &
                                MF_New(lvl)%dm,         &
                                1,                      &
                                nGhost_Vec(1)           )
                                
    CALL AMReX_BuildMask( Ghost_Mask,       &
                          amrex_geom(lvl),  &
                          iCovered,         &
                          iNotCovered,      &
                          iOutside,         &
                          iInterior         )


    CALL amrex_mfiter_build(mfi, MF_New(lvl), tiling = .true. )

    DO WHILE(mfi%next())

        Mask_PTR => Level_Mask%dataPtr(mfi)
        New_PTR => MF_New(lvl)%dataPtr(mfi)
        Old_PTR => MF_Old(lvl)%dataPtr(mfi)

        Box = mfi%tilebox()
        
        iEL = Box%lo
        iEU = Box%hi
                             
        DO re = iEL(1),iEU(1)
        DO te = iEL(2),iEU(2)
        DO pe = iEL(3),iEU(3)
        
            IF ( Mask_PTR(RE,TE,PE,1) == iLeaf ) THEN
                DO i = 1,Caller_Quad_DOF

                
                    Here = (iU_CF-1)*Caller_Quad_DOF+i
                    New_Value = New_PTR(re,te,pe,Here)
                    Old_Value = Old_PTR(re,te,pe,Here)
                    
                    Difference = abs(New_Value - Old_Value)
        !            PRINT*,New_Value,Old_Value,Difference
                    IF ( Difference > Max_Difference ) THEN
                        Max_Difference = Difference
                    END IF
                
                END DO ! i
            END IF
        END DO ! re
        END DO ! te
        END DO ! pe
        
    END DO !     While Loop
END DO !         Lvl Loop

Calculate_Max_Difference = Max_Difference

END FUNCTION Calculate_Max_Difference








 !+101+############################################################!
!                                                                   !
!        Calculate_Error_Norms                                      !
!                                                                   !
 !#################################################################!
FUNCTION Calculate_Error_Norms( MF_Result )

REAL(idp)                                   ::  Calculate_Error_Norms(3)
TYPE(amrex_multifab),   INTENT(IN)          ::  MF_Result(0:AMReX_Num_Levels-1)


INTEGER,   CONTIGUOUS, POINTER              ::  Mask_PTR(:,:,:,:)
REAL(idp), CONTIGUOUS, POINTER              ::  Result_PTR(:,:,:,:)


INTEGER                                     ::  lvl
TYPE(amrex_mfiter)                          ::  mfi
TYPE(amrex_box)                             ::  Box
TYPE(amrex_imultifab)                       ::  Level_Mask
TYPE(amrex_imultifab)                       ::  Ghost_Mask

INTEGER,    DIMENSION(1:3)                  ::  iEL, iEU
INTEGER                                     ::  re, te, pe
INTEGER                                     ::  rd, td, pd
INTEGER                                     ::  i
INTEGER                                     ::  Here

REAL(idp)                                   ::  Dif
REAL(idp)                                   ::  SumDif
REAL(idp)                                   ::  SumSol
REAL(idp)                                   ::  SumDif2
REAL(idp)                                   ::  SumSol2
REAL(idp)                                   ::  MaxDif
REAL(idp)                                   ::  MaxSol

REAL(idp)                                   ::  Sol
REAL(idp)                                   ::  Result

REAL(idp),  DIMENSION(1:Caller_NQ(1))       ::  Cur_R_Locs

REAL(idp)                                   ::  DROT

! Entry-Wise Error Norms

Calculate_Error_Norms = 0.0_idp

MaxDif = 0.0_idp
MaxSol = 0.0_idp
SumDif = 0.0_idp
SumSol = 0.0_idp
SumDif2 = 0.0_idp
SumSol2 = 0.0_idp

DO lvl = 0,AMReX_Num_Levels-1

    DROT = Level_dx(Lvl,1)/(Caller_xL(2)-Caller_xL(1))

    nGhost_Vec = MF_Result(lvl)%nghostvect()

    !
    !   MakeFineMask
    !
    IF ( lvl < AMReX_Num_Levels-1 ) THEN
        CALL AMReX_MakeFineMask(  Level_Mask,               &
                                  MF_Result(lvl)%ba,        &
                                  MF_Result(lvl)%dm,        &
                                  nGhost_Vec,               &
                                  MF_Result(lvl+1)%ba,      &
                                  iLeaf, iTrunk             )
    ELSE
        ! Create Level_Mask all equal to 1
        CALL amrex_imultifab_build( Level_Mask,             &
                                    MF_Result(lvl)%ba,      &
                                    MF_Result(lvl)%dm,      &
                                    1,                      &
                                    nGhost_Vec(1)           )
        CALL Level_Mask%SetVal(iLeaf)
    END IF
    
    !
    !   BuildMask
    !
    CALL amrex_imultifab_build( Ghost_Mask,             &
                                MF_Result(lvl)%ba,      &
                                MF_Result(lvl)%dm,      &
                                1,                      &
                                nGhost_Vec(1)           )
                                
    CALL AMReX_BuildMask( Ghost_Mask,       &
                          amrex_geom(lvl),  &
                          iCovered,         &
                          iNotCovered,      &
                          iOutside,         &
                          iInterior         )


    CALL amrex_mfiter_build(mfi, MF_Result(lvl), tiling = .true. )

    DO WHILE(mfi%next())

        Mask_PTR => Level_Mask%dataPtr(mfi)
        Result_PTR  => MF_Result(lvl)%dataPtr(mfi)

        Box = mfi%tilebox()
        
        iEL = Box%lo
        iEU = Box%hi
                             
        DO re = iEL(1),iEU(1)
        DO te = iEL(2),iEU(2)
        DO pe = iEL(3),iEU(3)
        
            IF ( Mask_PTR(RE,TE,PE,1) == iLeaf ) THEN
        
                CUR_R_LOCS(:) = DROT * (Caller_RQ_xlocs(:) + 1.0_idp        &
                              + re*(Caller_xL(2)-Caller_xL(1)))
        
                DO pd = 1,Caller_NQ(3)
                DO td = 1,Caller_NQ(2)
                DO rd = 1,Caller_NQ(1)

                    Here = Quad_Map(rd, td, pd,     &
                                    Caller_NQ(1),   &
                                    Caller_NQ(2),   &
                                    Caller_NQ(3))

                    Result = Result_PTR(re,te,pe,Here)
                    Sol = HCT_Solution(CUR_R_LOCS(rd))
                    Sol = abs(Sol)
                    Dif = abs(Result - Sol)
                    
                    SumDif = SumDif + Dif
                    SumSol = SumSol + Sol
                    
                    SumDif2 = SumDif2 + Dif*Dif
                    SumSol2 = SumSol2 + Sol*Sol
                    
                    IF (Dif > MaxDif) THEN
                        MaxDif = Dif
                    END IF
                    
                    IF (Sol > MaxSol) THEN
                        MaxSol = Sol
                    END IF

                END DO ! rd
                END DO ! td
                END DO ! pd
                
            END IF
            
            
        END DO ! re
        END DO ! te
        END DO ! pe
        
    END DO !     While Loop
END DO !         Lvl Loop

Calculate_Error_Norms(1) = SumDif/SumSol
Calculate_Error_Norms(2) = sqrt(SumDif2)/sqrt(SumSol2)
IF ( MaxSol .NE. 0.0_idp ) THEN
    Calculate_Error_Norms(3) = MaxDif/MaxSol
ELSE
    Calculate_Error_Norms(3) = 10000.0
END IF

END FUNCTION Calculate_Error_Norms












END MODULE Driver_ConFactor_Loop_Module

