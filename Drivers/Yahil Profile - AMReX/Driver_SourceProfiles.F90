   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Driver_SourceProfiles_Module                                          !##!
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



USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag





IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!          Load_AMReX_Profile                                           	    !
!                                                                               !
!###############################################################################!
SUBROUTINE Load_AMReX_Profile(  Profile,                &
                                nComps,                 &
                                nGhost,                 &
                                nLevels,                &
                                NE_Lower,               &
                                NE_Upper,               &
                                GM_Src_Input,           &
                                MF_Src_Input            )
                
CHARACTER(LEN=*),                   INTENT(IN)          ::  Profile

INTEGER,                            INTENT(IN)          ::  nComps
INTEGER,                            INTENT(IN)          ::  nGhost
INTEGER,                            INTENT(IN)          ::  nLevels
INTEGER, DIMENSION(3),              INTENT(IN)          ::  NE_Lower
INTEGER, DIMENSION(3),              INTENT(IN)          ::  NE_Upper

TYPE(amrex_geometry),   DIMENSION(0:nLevels-1), INTENT(INOUT)   ::  GM_Src_Input
TYPE(amrex_multifab),   DIMENSION(0:nLevels-1), INTENT(INOUT)   ::  MF_Src_Input


IF ( Verbose_Flag ) THEN
    WRITE(*,'(A,A)')"-Loading AMReX structure profile : ",TRIM(Profile)
END IF



SELECT CASE (Profile)
    CASE( 'SingleLayer' )

        CALL AMReX_Profile_SingleLayer( nComps,                 &
                                        nGhost,                 &
                                        nLevels,                &
                                        NE_Lower,               &
                                        NE_Upper,               &
                                        GM_Src_Input,           &
                                        MF_Src_Input            )

    CASE( 'TwoLayer' )

        CALL AMReX_Profile_TwoLayer(    nComps,                 &
                                        nGhost,                 &
                                        nLevels,                &
                                        NE_Lower,               &
                                        NE_Upper,               &
                                        GM_Src_Input,           &
                                        MF_Src_Input            )


    CASE DEFAULT

        WRITE(*,'(A,A,A)')'AMReX Profile :',TRIM(Profile),' not recognized.'
        WRITE(*,'(A)' )'Poseidon is STOPing.'
        STOP
END SELECT


IF ( Verbose_Flag ) THEN
    WRITE(*,'(A,A)')"-AMReX structure profile loaded."
END IF

END SUBROUTINE Load_AMReX_Profile






!+201+##########################################################################!
!                                                                               !
!           AMReX_Profile_SingleLayer                                           !
!                                                                               !
!###############################################################################!
SUBROUTINE AMReX_Profile_SingleLayer(   nComps,             &
                                        nGhost,             &
                                        nLevels,            &
                                        NE_Lower,           &
                                        NE_Upper,           &
                                        GM_Src_Input,       &
                                        MF_Src_Input        )
                

INTEGER,                            INTENT(IN)              ::  nComps
INTEGER,                            INTENT(IN)              ::  nGhost
INTEGER,                            INTENT(IN)              ::  nLevels
INTEGER, DIMENSION(3),              INTENT(IN)              ::  NE_Lower
INTEGER, DIMENSION(3),              INTENT(IN)              ::  NE_Upper

TYPE(amrex_multifab),   DIMENSION(0:nLevels-1), INTENT(INOUT)   ::  MF_Src_Input
TYPE(amrex_geometry),   DIMENSION(0:nLevels-1), INTENT(INOUT)   ::  GM_Src_Input



INTEGER                                                     ::  lvl
INTEGER,                DIMENSION(3,0:nLevels-1)            ::  MaxGridSize
TYPE(amrex_box),        DIMENSION(0:nLevels-1)              ::  BX_Src_Layer
TYPE(amrex_boxarray),   DIMENSION(0:nLevels-1)              ::  BA_Src_Layer
TYPE(amrex_distromap),  DIMENSION(0:nLevels-1)              ::  DM_Src_Layer



MaxGridSize(:,0) = NE_Upper
BX_Src_Layer(0) = amrex_box(NE_Lower, NE_Upper)

CALL Build_AMReX_Layer( nComps,             &
                        nGhost,             &
                        MaxGridSize(:,0),   &
                        BX_Src_Layer(0),    &
                        BA_Src_Layer(0),    &
                        DM_Src_Layer(0),    &
                        GM_Src_Input(0),    &
                        MF_Src_Input(0)     )



END SUBROUTINE AMReX_Profile_SingleLayer






!+201+##########################################################################!
!                                                                               !
!           AMReX_Profile_TwoLayer                                              !
!                                                                               !
!###############################################################################!
SUBROUTINE AMReX_Profile_TwoLayer(  nComps,                 &
                                    nGhost,                 &
                                    nLevels,                &
                                    NE_Lower,               &
                                    NE_Upper,               &
                                    GM_Src_Input,           &
                                    MF_Src_Input            )
            

INTEGER,                            INTENT(IN)              ::  nComps
INTEGER,                            INTENT(IN)              ::  nGhost
INTEGER,                            INTENT(IN)              ::  nLevels
INTEGER, DIMENSION(3),              INTENT(IN)              ::  NE_Lower
INTEGER, DIMENSION(3),              INTENT(IN)              ::  NE_Upper

TYPE(amrex_geometry),   DIMENSION(0:nLevels-1), INTENT(INOUT)   ::  GM_Src_Input
TYPE(amrex_multifab),   DIMENSION(0:nLevels-1), INTENT(INOUT)   ::  MF_Src_Input

INTEGER                                                     ::  lvl
INTEGER,                DIMENSION(3,0:nLevels-1)            ::  MaxGridSize
TYPE(amrex_box),        DIMENSION(0:nLevels-1)              ::  BX_Src_Layer
TYPE(amrex_boxarray),   DIMENSION(0:nLevels-1)              ::  BA_Src_Layer
TYPE(amrex_distromap),  DIMENSION(0:nLevels-1)              ::  DM_Src_Layer


IF ( nLevels > 1 ) THEN

    ! Level 0 - Coarsest Layer
    MaxGridSize(:,0) = NE_Upper
    BX_Src_Layer(0) = amrex_box(NE_Lower, NE_Upper)

    PRINT*,"Before Layer 0 Build"
    CALL Build_AMReX_Layer( nComps,             &
                            nGhost,             &
                            MaxGridSize(:,0),   &
                            BX_Src_Layer(0),    &
                            BA_Src_Layer(0),    &
                            DM_Src_Layer(0),    &
                            GM_Src_Input(0),    &
                            MF_Src_Input(0)     )




    ! Level 1 - Finest Layer
    MaxGridSize(:,1) = NE_Upper/2
    BX_Src_Layer(1)  = amrex_box(NE_Lower, NE_Upper/2)  ! Box Spans Lower half of 1D space.

    PRINT*,"Before Layer 1 Build"
    CALL Build_AMReX_Layer( nComps,             &
                            nGhost,             &
                            MaxGridSize(:,1),   &
                            BX_Src_Layer(1),    &
                            BA_Src_Layer(1),    &
                            DM_Src_Layer(1),    &
                            GM_Src_Input(1),    &
                            MF_Src_Input(1)     )



ELSE

    WRITE(*,'(A)') "Warning : Fatal Poseidon Error "
    WRITE(*,'(A)') "Insufficient number of levels for the 'TwoLayer' AMReX Profile."
    WRITE(*,'(A)') "Input : nLevels = ",nLevels," < 2."
    WRITE(*,'(A)') "Poseidon is STOPing."
    STOP

END IF


END SUBROUTINE AMReX_Profile_TwoLayer






!+301+##########################################################################!
!                                                                               !
!     Build_AMReX_Layer                                                         !
!                                                                               !
!###############################################################################!
SUBROUTINE Build_AMReX_Layer(   nComps,                 &
                                nGhost,                 &
                                GridSize_Layer,         &
                                BX_Src_Layer,           &
                                BA_Src_Layer,           &
                                DM_Src_Layer,           &
                                GM_Src_Layer,           &
                                MF_Src_Layer            )

INTEGER,                            INTENT(IN)          ::  nComps
INTEGER,                            INTENT(IN)          ::  nGhost
INTEGER, DIMENSION(3),              INTENT(IN)          ::  GridSize_Layer

TYPE(amrex_box),                    INTENT(INOUT)       ::  BX_Src_Layer
TYPE(amrex_boxarray),               INTENT(INOUT)       ::  BA_Src_Layer
TYPE(amrex_distromap),              INTENT(INOUT)       ::  DM_Src_Layer
TYPE(amrex_geometry),               INTENT(INOUT)       ::  GM_Src_Layer
TYPE(amrex_multifab),               INTENT(INOUT)       ::  MF_Src_Layer


PRINT*,"A"
CALL amrex_boxarray_build(  BA_Src_Layer,           &
                            BX_Src_Layer            )

PRINT*,"B"
CALL BA_Src_Layer%maxSize(  GridSize_Layer          )


PRINT*,"C"
PRINT*,GM_Src_Layer%dx
PRINT*,BX_Src_Layer%hi
CALL amrex_geometry_build(  GM_Src_Layer,           &
                            BX_Src_Layer            )
PRINT*,GM_Src_Layer%dx
PRINT*,BX_Src_Layer%hi

PRINT*,"D"
CALL amrex_distromap_build( DM_Src_Layer,           &
                            BA_Src_Layer            )

PRINT*,"E"
CALL amrex_multifab_build(  MF_Src_Layer,           &
                            BA_Src_Layer,           &
                            DM_Src_Layer,           &
                            nComps,                 &
                            nGhost                  )

PRINT*,"F"
CALL MF_Src_Layer % SetVal(0.0_idp)


END SUBROUTINE Build_AMReX_Layer




END MODULE Driver_SourceProfiles_Module
