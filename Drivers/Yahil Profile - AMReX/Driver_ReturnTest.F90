   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Driver_ReturnTest                                                     !##!
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


USE amrex_base_module

USE amrex_box_module, &
            ONLY :  amrex_box

USE amrex_boxarray_module,  &
            ONLY :  amrex_boxarray,                 &
                    amrex_boxarray_build,           &
                    amrex_boxarray_destroy

USE amrex_distromap_module, &
            ONLY :  amrex_distromap,                &
                    amrex_distromap_build,          &
                    amrex_distromap_destroy

USE amrex_multifab_module,  &
            ONLY :  amrex_multifab,                 &
                    amrex_multifab_build

USE amrex_amrcore_module,   &
            ONLY :  amrex_amrcore_init,             &
                    amrex_init_virtual_functions,   &
                    amrex_init_from_scratch,        &
                    amrex_ref_ratio,                &
                    amrex_max_level

USE amrex_geometry_module,  &
            ONLY :  amrex_geometry


USE Poseidon_AMReX_MakeFineMask_Module, &
            ONLY :  AMReX_MakeFineMask

USE Functions_Quadrature, &
            ONLY :  Initialize_LG_Quadrature_Locations, &
                    Initialize_Trapezoid_Quadrature_Locations

USE Variables_AMReX_Source, &
            ONLY :  iCoarse,                &
                    iFine

USE Return_Functions_AMReX, &
            ONLY :  Poseidon_Return_ConFactor_AMReX,    &
                    Poseidon_Return_All_AMReX

IMPLICIT NONE

INTEGER                 ::  iU_K11 = 6
INTEGER                 ::  iU_K12 = 7
INTEGER                 ::  iU_K13 = 8
INTEGER                 ::  iU_K22 = 9
INTEGER                 ::  iU_K23 = 10
INTEGER                 ::  iU_K33 = 11


CONTAINS



 !+101+####################################################!
!                                                           !
!          Return_Test	                                    !
!                                                           !
 !#########################################################!
SUBROUTINE Return_Test( nLevels, NQ, MF_Source )

INTEGER,                INTENT(IN)              ::  nLevels
INTEGER, DIMENSION(3),  INTENT(IN)              ::  NQ
TYPE(amrex_multifab),   INTENT(IN)              ::  MF_Source(0:nLevels-1)



INTEGER                                         ::  MF_Results_nVars    = 11
INTEGER                                         ::  MF_Results_nGhosts  = 0


INTEGER                                         ::  Num_Quad
INTEGER                                         ::  MF_Results_nComps


REAL(idp),  DIMENSION(NQ(1))                    ::  R_Quad
REAL(idp),  DIMENSION(NQ(2))                    ::  T_Quad
REAL(idp),  DIMENSION(NQ(3))                    ::  P_Quad
REAL(idp),  DIMENSION(2)                        ::  xL

TYPE(amrex_multifab),           ALLOCATABLE     ::  MF_Results(:)
TYPE(amrex_boxarray),           ALLOCATABLE     ::  BA_Results(:)
TYPE(amrex_distromap),          ALLOCATABLE     ::  DM_Results(:)
TYPE(amrex_geometry),           ALLOCATABLE     ::  GM_Results(:)

INTEGER                                         ::  lvl



ALLOCATE( MF_Results(0:nLevels-1) )
ALLOCATE( BA_Results(0:nLevels-1) )
ALLOCATE( DM_Results(0:nLevels-1) )
ALLOCATE( GM_Results(0:nLevels-1) )


Num_Quad = NQ(1)*NQ(2)*NQ(3)
MF_Results_nComps = Num_Quad*MF_Results_nVars

DO lvl = 0,nLevels-1

    CALL amrex_multifab_build(  MF_Results(lvl),            &
                                MF_Source(Lvl)%BA,          &
                                MF_Source(Lvl)%DM,          &
                                MF_Results_nComps,          &
                                MF_Results_nGhosts          )

    CALL MF_Results(lvl)%SetVal(0.0_idp)
END DO


xL(1) = -1.0_idp
xL(2) = +1.0_idp
R_Quad = Initialize_LG_Quadrature_Locations(NQ(1))
T_Quad = Initialize_LG_Quadrature_Locations(NQ(2))
P_Quad = Initialize_Trapezoid_Quadrature_Locations(NQ(3))

CALL Poseidon_Return_ConFactor_AMReX( NQ,                   &
                                      R_Quad,               &
                                      T_Quad,               &
                                      P_Quad,               &
                                      xL(1),                &
                                      xL(2),                &
                                      nLevels,              &
                                      MF_Results            )

!CALL Output_Variable( nLevels, NQ, MF_Results, iU_CF )


CALL Poseidon_Return_ALL_AMReX( NQ,                   &
                                R_Quad,               &
                                T_Quad,               &
                                P_Quad,               &
                                xL(1),                &
                                xL(2),                &
                                nLevels,              &
                                MF_Results            )

!CALL Output_All_Variables( nLevels, NQ, MF_Results )

END SUBROUTINE Return_Test










 !+101+####################################################!
!                                                           !
!          Output_ConFactor                                 !
!                                                           !
 !#########################################################!
SUBROUTINE Output_Variable( nLevels, NQ, MF_Results, iU )

INTEGER,                INTENT(IN)              ::  nLevels
INTEGER, DIMENSION(3),  INTENT(IN)              ::  NQ
TYPE(amrex_multifab),   INTENT(IN)              ::  MF_Results(0:nLevels-1)
INTEGER,                INTENT(IN)              ::  iU

INTEGER                                         ::  lvl
TYPE(amrex_mfiter)                              ::  mfi
TYPE(amrex_box)                                 ::  Box
TYPE(amrex_imultifab)                           ::  Level_Mask
INTEGER                                         ::  nComp
INTEGER, DIMENSION(3)                           ::  iEL, iEU
INTEGER,    CONTIGUOUS, POINTER                 ::  Mask_PTR(:,:,:,:)
REAL(idp),  CONTIGUOUS, POINTER                 ::  Results_PTR(:,:,:,:)

INTEGER                                         ::  re, te, pe
INTEGER                                         ::  rd, td, pd
INTEGER                                         ::  Here
INTEGER                                         ::  Num_Quad

Num_Quad = NQ(1)*NQ(2)*NQ(3)


DO lvl = nLevels-1,0,-1

    !
    !   MakeFineMask
    !
    IF ( lvl < nLevels-1 ) THEN
        CALL AMReX_MakeFineMask(  Level_Mask,               &
                                  MF_Results(lvl)%ba,        &
                                  MF_Results(lvl)%dm,        &
                                  MF_Results(lvl+1)%ba,      &
                                  iCoarse, iFine            )
    ELSE
        ! Create Level_Mask all equal to 1
        CALL amrex_imultifab_build( Level_Mask,             &
                                    MF_Results(lvl)%ba,      &
                                    MF_Results(lvl)%dm,      &
                                    1,                      &
                                    0                       )
        CALL Level_Mask%SetVal(iCoarse)
    END IF


    CALL amrex_mfiter_build(mfi, MF_Results(lvl), tiling = .true. )

    DO WHILE(mfi%next())

        Results_PTR => MF_Results(lvl)%dataPtr(mfi)
        Mask_PTR   => Level_Mask%dataPtr(mfi)

        Box = mfi%tilebox()
        nComp =  MF_Results(lvl)%ncomp()

        iEL = Box%lo
        iEU = Box%hi



        DO re = iEL(1),iEU(1)
        DO te = iEL(2),iEU(2)
        DO pe = iEL(3),iEU(3)

            IF ( Mask_PTR(RE,TE,PE,1) == iCoarse ) THEN

            DO pd = 1,NQ(3)
            DO td = 1,NQ(2)
            DO rd = 1,NQ(1)


                Here = (iU-1) * Num_Quad        &
                     + (pd-1)*NQ(1)*NQ(2)       &
                     + (td-1)*NQ(1)             &
                     + rd

                PRINT*,re,te,pe,rd,td,pd,iU,Results_PTR(re,te,pe,Here)

            END DO ! rd Loop
            END DO ! td Loop
            END DO ! pd Loop

            END IF

        END DO ! pe
        END DO ! te
        END DO ! re

    END DO

    CALL amrex_mfiter_destroy(mfi)
    CALL amrex_imultifab_destroy( Level_Mask )

END DO ! lvl


END SUBROUTINE Output_Variable




!+101+####################################################!
!                                                           !
!          Output_ConFactor                                 !
!                                                           !
 !#########################################################!
SUBROUTINE Output_ALL_Variables( nLevels, NQ, MF_Results )

INTEGER,                INTENT(IN)              ::  nLevels
INTEGER, DIMENSION(3),  INTENT(IN)              ::  NQ
TYPE(amrex_multifab),   INTENT(IN)              ::  MF_Results(0:nLevels-1)

INTEGER                                         ::  lvl
TYPE(amrex_mfiter)                              ::  mfi
TYPE(amrex_box)                                 ::  Box
TYPE(amrex_imultifab)                           ::  Level_Mask
INTEGER                                         ::  nComp
INTEGER, DIMENSION(3)                           ::  iEL, iEU
INTEGER,    CONTIGUOUS, POINTER                 ::  Mask_PTR(:,:,:,:)
REAL(idp),  CONTIGUOUS, POINTER                 ::  Results_PTR(:,:,:,:)

INTEGER                                         ::  re, te, pe
INTEGER                                         ::  rd, td, pd
INTEGER                                         ::  Here
INTEGER                                         ::  Num_Quad
INTEGER                                         ::  iU


Num_Quad = NQ(1)*NQ(2)*NQ(3)


DO lvl = nLevels-1,0,-1

    !
    !   MakeFineMask
    !
    IF ( lvl < nLevels-1 ) THEN
        CALL AMReX_MakeFineMask(  Level_Mask,               &
                                  MF_Results(lvl)%ba,        &
                                  MF_Results(lvl)%dm,        &
                                  MF_Results(lvl+1)%ba,      &
                                  iCoarse, iFine            )
    ELSE
        ! Create Level_Mask all equal to 1
        CALL amrex_imultifab_build( Level_Mask,             &
                                    MF_Results(lvl)%ba,      &
                                    MF_Results(lvl)%dm,      &
                                    1,                      &
                                    0                       )
        CALL Level_Mask%SetVal(iCoarse)
    END IF


    CALL amrex_mfiter_build(mfi, MF_Results(lvl), tiling = .true. )

    DO WHILE(mfi%next())

        Results_PTR => MF_Results(lvl)%dataPtr(mfi)
        Mask_PTR   => Level_Mask%dataPtr(mfi)

        Box = mfi%tilebox()
        nComp =  MF_Results(lvl)%ncomp()

        iEL = Box%lo
        iEU = Box%hi



        DO re = iEL(1),iEU(1)
        DO te = iEL(2),iEU(2)
        DO pe = iEL(3),iEU(3)

            IF ( Mask_PTR(RE,TE,PE,1) == iCoarse ) THEN

            DO pd = 1,NQ(3)
            DO td = 1,NQ(2)
            DO rd = 1,NQ(1)




                PRINT*,re,te,pe,rd,td,pd
                DO iU = 1,11
                    Here = (iU-1) * Num_Quad        &
                         + (pd-1)*NQ(1)*NQ(2)       &
                         + (td-1)*NQ(1)             &
                         + rd
                    PRINT*,Results_PTR(re,te,pe,Here)
                END DO
            END DO ! rd Loop
            END DO ! td Loop
            END DO ! pd Loop

            END IF

        END DO ! pe
        END DO ! te
        END DO ! re

    END DO

    CALL amrex_mfiter_destroy(mfi)
    CALL amrex_imultifab_destroy( Level_Mask )

END DO ! lvl



END SUBROUTINE Output_ALL_Variables




END MODULE Driver_ReturnTest
