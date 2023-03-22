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


use amrex_base_module

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
            ONLY :  Centimeter,         &
                    Gram

USE Variables_External, &
            ONLY :  MacLaurin_SemiMinor,         &
                    MacLaurin_SemiMajor,         &
                    MacLaurin_Ecc,          &
                    MacLaurin_SphereType,   &
                    MacLaurin_Rho


USE Variables_AMReX_Core, &
            ONLY :  MF_Source,          &
                    AMReX_Num_Levels


USE Variables_Driver_AMReX, &
            ONLY :  MF_Driver_Source,   &
                    MF_Src_nComps,      &
                    MF_Src_nGhost,      &
                    nLevels

USE Variables_AMReX_Source, &
            ONLY :  Source_PTR,         &
                    Mask_PTR,           &
                    iLeaf,            &
                    iTrunk

USE Poseidon_AMReX_MakeFineMask_Module, &
            ONLY :  AMReX_MakeFineMask

USE Poseidon_Interface_Source_Input, &
            ONLY :  Poseidon_Input_Sources

USE Variables_Interface, &
            ONLY :  Caller_Quad_DOF

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


USE MPI

IMPLICIT NONE


CONTAINS



 !+101+####################################################!
!                                                           !
!     Driver_SetSource                                      !
!                                                           !
 !#########################################################!
SUBROUTINE Driver_SetSource( )


INTEGER                                         ::  nVars_Source

IF ( Verbose_Flag )             &
        CALL Driver_Init_Message('Creating AMReX source variables.')
IF ( Verbose_Flag )             &
        CALL Driver_Init_Message('Initializing Homogenous Spheroid Source  Multifab.')


CALL TimerStart( Timer_Driver_SetSource_InitTest )


MacLaurin_SemiMinor  = 250.0_idp*Centimeter
MacLaurin_SemiMajor  = 2500.0_idp*Centimeter
MacLaurin_SphereType = 'O'
MacLaurin_Rho        = 1.0E20_idp*(Gram/Centimeter**3)

CALL amrex_init_virtual_functions &
       ( VF_Make_New_Level_From_Scratch, &
         VF_Make_New_Level_From_Coarse, &
         VF_Remake_Level, &
         VF_Clear_Level, &
         VF_Error_Estimate )


nVars_Source    = 5
MF_Src_nComps   = nVars_Source*Caller_Quad_DOF
MF_Src_nGhost   = 0


ALLOCATE( MF_Driver_Source(0:nLevels-1) )

CALL amrex_init_from_scratch( 0.0_idp )

CALL TimerStop( Timer_Driver_SetSource_InitTest )





CALL Poseidon_Input_Sources( MF_Driver_Source,    &
                             MF_Src_nComps        )






END SUBROUTINE Driver_SetSource









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
