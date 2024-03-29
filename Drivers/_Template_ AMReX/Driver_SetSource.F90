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
                    AMReX_Num_Levels


USE Variables_Driver_AMReX, &
            ONLY :  MF_Driver_Source,   &
                    MF_Src_nComps,      &
                    MF_Src_nGhost,      &
                    dt, t_old, t_new,   &
                    StepNo,             &
                    nLevels

USE Variables_AMReX_Source, &
            ONLY :  Source_PTR,         &
                    Mask_PTR,           &
                    iCoarse,            &
                    iFine

USE Poseidon_AMReX_MakeFineMask_Module, &
            ONLY :  AMReX_MakeFineMask

USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag


USE Variables_MPI, &
            ONLY :  myID_Poseidon,      &
                    Poseidon_Comm_World

USE Source_Input_AMReX, &
            ONLY :  Poseidon_Input_Sources_AMREX


USE Poseidon_MPI_Utilities_Module, &
            ONLY :  STOP_MPI,               &
                    MPI_Master_Print,       &
                    MPI_All_Print

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



!+101+##########################################################################!
!                                                                               !
!     Driver_SetSource                                                			!
!                                                                               !
!###############################################################################!
SUBROUTINE Driver_SetSource( )


IF ( Verbose_Flag ) THEN
    WRITE(*,'(A)')"In Driver, Creating AMReX Source Variables."
END IF



CALL TimerStart( Timer_Driver_SetSource_InitTest )



CALL amrex_init_virtual_functions &
       ( VF_Make_New_Level_From_Scratch, &
         VF_Make_New_Level_From_Coarse, &
         VF_Remake_Level, &
         VF_Clear_Level, &
         VF_Error_Estimate )



!PRINT*,"nLevels",nlevels
ALLOCATE( MF_Driver_Source(0:nLevels-1) )

CALL amrex_init_from_scratch( 0.0_idp )




CALL Poseidon_Input_Sources_AMREX( MF_Driver_Source, MF_Src_nComps, nLevels )


CALL TimerStop( Timer_Driver_SetSource_InitTest )





DEALLOCATE(MF_Driver_Source)



END SUBROUTINE Driver_SetSource
























!+201+##########################################################################!
!                                                                               !
!     SetSource_Parallel                                                            !
!                                                                               !
!###############################################################################!
SUBROUTINE SetSource_Parallel(  nDOF, nVars, nLevels,          &
                                NE_Lower, NE_Upper,             &
                                Local_E, Local_S, Local_Si,     &
                                MF_SRC_Input                      )

INTEGER,    INTENT(IN)                              ::  nDOF, nVars, nLevels
INTEGER,    INTENT(IN), DIMENSION(3)                ::  NE_Lower, NE_Upper
REAL(idp),  INTENT(IN), DIMENSION( nDOF,        &
                                NE_Upper(1),    &
                                NE_Upper(2),    &
                                NE_Upper(3)     )   ::  Local_E

REAL(idp),  INTENT(IN), DIMENSION( nDOF,        &
                                NE_Upper(1),    &
                                NE_Upper(2),    &
                                NE_Upper(3)     )   ::  Local_S

REAL(idp),  INTENT(IN), DIMENSION( nDOF,        &
                                NE_Upper(1),    &
                                NE_Upper(2),    &
                                NE_Upper(3),    &
                                1:3             )   ::  Local_Si


TYPE(amrex_multifab), INTENT(INOUT)                 ::  MF_SRC_Input(0:nLevels-1)

REAL(idp), CONTIGUOUS, POINTER                      :: p(:,:,:,:)

INTEGER                                             :: PE, TE, RE, Var
INTEGER, DIMENSION(3)                               :: ELo, EHi

INTEGER                                             :: Here, There, lvl

TYPE(amrex_mfiter)                                  :: mfi
TYPE(amrex_box)                                     :: Box

!INTEGER                                             :: ierr

PRINT*,"In SetSource_Parallel ", myID_Poseidon




DO lvl = 0,nLevels-1
    PRINT*,"A"
    CALL amrex_mfiter_build(mfi, MF_SRC_Input(lvl), tiling = .false. )
    DO WHILE(mfi%next())

!        CALL MPI_All_Print("In mfiter loop")
        PRINT*,"B"
        p => MF_SRC_Input(lvl)%dataPtr(mfi)
        Box = mfi%tilebox()

        PRINT*,"C"
        ELo = Box%lo
        EHi = Box%hi

        PRINT*,"D"
        DO PE = ELo(3), EHi(3)
        DO TE = ELo(2), EHi(2)
        DO RE = ELo(1), EHi(1)

            Here  = 1
            There = nDOF
            p(RE,TE,PE, Here:There ) = Local_E(:,RE,TE,PE)
            
            Here  = 1*nDOF+1
            There = 2*nDOF
            p(RE,TE,PE, Here:There ) = Local_S(:,RE,TE,PE)

            DO Var = 3,5
                Here  = (Var-1)*nDOF+1
                There = Var*nDOF
                p(RE,TE,PE,Here:There) = Local_Si(:,RE,TE,PE,Var-2)
            END DO

        

        END DO
        END DO
        END DO

    END DO
    CALL amrex_mfiter_destroy(mfi)
END DO ! lvl




!CALL MPI_Barrier(Poseidon_Comm_World, ierr)
!CALL MPI_Master_Print("At the end of SetSource_Parallel")
!CALL STOP_MPI(ierr)


END SUBROUTINE SetSource_Parallel









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
INTEGER, DIMENSION(3)                           ::  iE
INTEGER, DIMENSION(3)                           ::  iEL, iEU
INTEGER                                         ::  nComp
INTEGER                                         ::  lvl


DO lvl = 0,AMReX_Num_Levels-1

    IF ( lvl < AMReX_Num_Levels-1 ) THEN
        CALL AMReX_MakeFineMask(  Level_Mask,               &
                                  MF_Source(lvl)%ba,        &
                                  MF_Source(lvl)%dm,        &
                                  MF_Source(lvl+1)%ba,      &
                                  iCoarse, iFine   )
    ELSE
        ! Create Level_Mask all equal to 1
        CALL amrex_imultifab_build( Level_Mask,             &
                                    MF_Source(lvl)%ba,      &
                                    MF_Source(lvl)%dm,      &
                                    1,                      &
                                    0                       )
        CALL Level_Mask%SetVal(iCoarse)
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
