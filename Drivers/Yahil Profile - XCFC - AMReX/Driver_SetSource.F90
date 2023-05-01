   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Driver_SetSource_Module                                              !##!
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
            ONLY :  Verbose_Flag

USE Poseidon_Message_Routines_Module, &
            ONLY :  Driver_Init_Message

USE amrex_base_module

USE amrex_box_module,       ONLY: &
  amrex_box
USE amrex_boxarray_module,  ONLY: &
  amrex_boxarray,         &
  amrex_boxarray_build,   &
  amrex_boxarray_destroy
USE amrex_distromap_module, ONLY: &
  amrex_distromap,       &
  amrex_distromap_build, &
  amrex_distromap_destroy
USE amrex_multifab_module,  ONLY: &
  amrex_multifab, &
  amrex_multifab_build
USE amrex_amrcore_module, ONLY: &
  amrex_amrcore_init, &
  amrex_init_virtual_functions, &
  amrex_init_from_scratch, &
  amrex_ref_ratio, &
  amrex_max_level


USE Poseidon_Units_Module, &
            ONLY :  Grav_Constant_G,    &
                    Speed_of_Light,     &
                    C_Square,           &
                    GR_Source_Scalar,   &
                    Centimeter,         &
                    Second,             &
                    Millisecond,         &
                    Erg,                &
                    Gram,               &
                    E_Units


USE Variables_AMReX_Core, &
            ONLY :  MF_Source,          &
                    AMReX_Num_Levels,   &
                    Source_PTR,         &
                    Mask_PTR

USE Parameters_AMReX, &
            ONLY :  iLeaf,            &
                    iTrunk

USE Variables_Driver_AMReX, &
            ONLY :  MF_Driver_Source,   &
                    MF_Src_nComps,      &
                    MF_Src_nGhost,      &
                    dt, t_old, t_new,   &
                    StepNo,             &
                    nLevels



USE Poseidon_AMReX_MakeFineMask_Module, &
            ONLY :  AMReX_MakeFineMask

USE Variables_MPI, &
            ONLY :  myID_Poseidon,      &
                    Poseidon_Comm_World

USE Poseidon_Interface_Source_Input, &
            ONLY :  Poseidon_Input_Sources


USE Poseidon_MPI_Utilities_Module, &
            ONLY :  STOP_MPI,               &
                    MPI_Master_Print,       &
                    MPI_All_Print

USE Variables_External, &
            ONLY :  SelfSim_T,              &
                    SelfSim_Kappa,          &
                    SelfSim_Gamma,          &
                    Central_E

USE Driver_AMReX_Virtual_Functions_Module, &
            ONLY :  VF_Make_New_Level_From_Scratch, &
                    VF_Make_New_Level_From_Coarse,  &
                    VF_Remake_Level,                &
                    VF_Clear_Level,                 &
                    VF_Error_Estimate

USE Timer_Routines_Module, &
            ONLY :  TimerStart,     &
                    TimerSTop


USE Timer_Variables_Module, &
            ONLY :  Timer_Driver_SetSource_InitTest,        &
                    Timer_Driver_SetSource_SetSource,       &
                    Timer_Driver_SetSource_Scale

USE External_Yahil_Profile_Module, &
            ONLY :  Calc_Yahil_Central_E

USE Variables_Interface, &
            ONLY :  Caller_Quad_DOF

USE MPI

IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!     Driver_SetSource                                                			!
!                                                                               !
!###############################################################################!
SUBROUTINE Driver_SetSource( Yahil_Params,              &
                             nLevels_Input              )

REAL(idp),  INTENT(IN), DIMENSION(3)                    ::  Yahil_Params
INTEGER,    INTENT(IN)                                  ::  nLevels_Input

INTEGER                                                 ::  nVars_Source


REAL(idp)                                               ::  Kappa_wUnits

IF ( Verbose_Flag ) CALL Driver_Init_Message('Creating AMReX source variables.')
IF ( Verbose_Flag ) CALL Driver_Init_Message('Initializing Yahil Source Multifab.')

CALL TimerStart( Timer_Driver_SetSource_InitTest )


SelfSim_T     = Yahil_Params(1)
SelfSim_Kappa = Yahil_Params(2)
SelfSim_Gamma = Yahil_Params(3)


Kappa_wUnits = SelfSim_Kappa*((Erg/Centimeter**3)/(Gram/Centimeter**3)**SelfSim_Gamma)


Central_E = Calc_Yahil_Central_E(SelfSim_T, SelfSim_Kappa, SelfSim_Gamma)


IF ( Verbose_Flag ) THEN
    WRITE(*,'(A)')'------------- Test Parameters ----------------'
    WRITE(*,'(A)')' Source Configuration : Yahil Self-Similar Collapse Profile'
    WRITE(*,'(A,ES12.5,A)') ' - Yahil Time      : ', SelfSim_T,' ms'
    WRITE(*,'(A,ES12.5)')   ' - Kappa           : ', Kappa_wUnits
    WRITE(*,'(A,ES12.5)')   ' - Gamma           : ', SelfSim_Gamma
    WRITE(*,'(A,ES12.5,A)') ' - Central E       : ', Central_E/E_Units," Erg/cm^3"
    WRITE(*,'(/)')
END IF


CALL amrex_init_virtual_functions &
       ( VF_Make_New_Level_From_Scratch, &
         VF_Make_New_Level_From_Coarse, &
         VF_Remake_Level, &
         VF_Clear_Level, &
         VF_Error_Estimate )



nVars_Source    = 5
MF_Src_nComps   = nVars_Source*Caller_Quad_DOF
MF_Src_nGhost   = 0



ALLOCATE( MF_Driver_Source(0:nLevels_Input-1) )
CALL amrex_init_from_scratch( 0.0_idp )

CALL TimerStop( Timer_Driver_SetSource_InitTest )




CALL Poseidon_Input_Sources(MF_Driver_Source,    &   ! Source Multifab
                            MF_Src_nComps        )   ! AMReX Levels





END SUBROUTINE Driver_SetSource



!+101+##########################################################################!
!                                                                               !
!     Driver_SetSource                                                			!
!                                                                               !
!###############################################################################!
SUBROUTINE Driver_CreateSource( Yahil_Params,               &
                                nLevels_Input               )

REAL(idp),  INTENT(IN), DIMENSION(3)                    ::  Yahil_Params
INTEGER,    INTENT(IN)                                  ::  nLevels_Input

INTEGER                                                 ::  nVars_Source


REAL(idp)                                               ::  Kappa_wUnits

IF ( Verbose_Flag ) CALL Driver_Init_Message('Creating AMReX source variables.')
IF ( Verbose_Flag ) CALL Driver_Init_Message('Initializing Yahil Source Multifab.')


CALL TimerStart( Timer_Driver_SetSource_InitTest )


SelfSim_T     = Yahil_Params(1)
SelfSim_Kappa = Yahil_Params(2)
SelfSim_Gamma = Yahil_Params(3)

Kappa_wUnits = SelfSim_Kappa*((Erg/Centimeter**3)/(Gram/Centimeter**3)**SelfSim_Gamma)


Central_E = Calc_Yahil_Central_E(SelfSim_T, SelfSim_Kappa, SelfSim_Gamma)


IF ( Verbose_Flag ) THEN
    WRITE(*,'(A)')'------------- Test Parameters ----------------'
    WRITE(*,'(A)')' Source Configuration : Yahil Self-Similar Collapse Profile'
    WRITE(*,'(A,ES12.5,A)') ' - Yahil Time      : ', SelfSim_T,' ms'
    WRITE(*,'(A,ES12.5)')   ' - Kappa           : ', Kappa_wUnits
    WRITE(*,'(A,ES12.5)')   ' - Gamma           : ', SelfSim_Gamma
    WRITE(*,'(A,ES12.5,A)') ' - Central E       : ', Central_E/E_Units," Erg/cm^3"
    WRITE(*,'(/)')
END IF


CALL amrex_init_virtual_functions &
       ( VF_Make_New_Level_From_Scratch, &
         VF_Make_New_Level_From_Coarse, &
         VF_Remake_Level, &
         VF_Clear_Level, &
         VF_Error_Estimate )



nVars_Source    = 5
MF_Src_nComps   = nVars_Source*Caller_Quad_DOF
MF_Src_nGhost   = 0


ALLOCATE( MF_Driver_Source(0:nLevels_Input-1) )
CALL amrex_init_from_scratch( 0.0_idp )

CALL TimerStop( Timer_Driver_SetSource_InitTest )


CALL Poseidon_Input_Sources(MF_Driver_Source,    &   ! Source Multifab
                            MF_Src_nComps        )   ! AMReX Levels

END SUBROUTINE Driver_CreateSource








!+201+##########################################################################!
!                                                                               !
!     SetSource_Parallel                                                            !
!                                                                               !
!###############################################################################!
SUBROUTINE AMReX_Source_Test()


TYPE(amrex_mfiter)                              ::  mfi
TYPE(amrex_box)                                 ::  Box

TYPE(amrex_imultifab)                           ::  Level_Mask

INTEGER                                         ::  re, te, pe
INTEGER, DIMENSION(3)                           ::  iEL, iEU
INTEGER                                         ::  nComp
INTEGER                                         ::  lvl
INTEGER, DIMENSION(1:3)                         ::  nGhost_Vec

nGhost_Vec = 0

DO lvl = 0,AMReX_Num_Levels-1

    IF ( lvl < AMReX_Num_Levels-1 ) THEN
        CALL AMReX_MakeFineMask(  Level_Mask,               &
                                  MF_Source(lvl)%ba,        &
                                  MF_Source(lvl)%dm,        &
                                  nGhost_Vec,               &
                                  MF_Source(lvl+1)%ba,      &
                                  iLeaf, iTrunk   )
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



        DO re = iEL(1),iEU(1)
        DO te = iEL(2),iEU(2)
        DO pe = iEL(3),iEU(3)

            PRINT*,lvl,re,te,pe
            PRINT*,Source_PTR(re,te,pe,:)

        END DO ! pe
        END DO ! te
        END DO ! re

    END DO
    CALL amrex_mfiter_destroy(mfi)
END DO ! lvl


END SUBROUTINE AMReX_Source_Test








END MODULE Driver_SetSource_Module
