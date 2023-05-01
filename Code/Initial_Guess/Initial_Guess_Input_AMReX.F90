   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE IG_Input_AMReX_Module                                                        !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
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
            ONLY : idp
        
USE Poseidon_Numbers_Module, &
            ONLY :  pi
            
USE Poseidon_Parameters, &
            ONLY :  Degree
            
USE Parameters_Variable_Indices, &
            ONLY :  iU_CF,          &
                    iU_LF,          &
                    iVB_S,          &
                    iVB_X,          &
                    iU_S1,          &
                    iU_S2,          &
                    iU_S3

USE Maps_Domain, &
            ONLY :  Map_To_FEM_Node,        &
                    FEM_Elem_Map
                    
USE Variables_Vectors, &
            ONLY :  dVA_Coeff_Vector,      &
                    dVB_Coeff_Vector

USE Parameters_AMReX, &
            ONLY :  iLeaf,                &
                    iTrunk

USE Variables_FEM_Module, &
            ONLY :  FEM_Node_xlocs

USE Maps_Fixed_Point, &
            ONLY :  FP_Array_Map_TypeB
                    
                    
#ifdef POSEIDON_AMREX_FLAG
use amrex_base_module

USE amrex_box_module,   &
            ONLY :  amrex_box
            
USE amrex_multifab_module,  &
            ONLY :  amrex_multifab,                 &
                    amrex_imultifab_build,          &
                    amrex_imultifab_destroy
                    
USE Variables_Interface, &
            ONLY :  Caller_NQ,                      &
                    Caller_Quad_DOF,                      &
                    Caller_xL,                      &
                    Caller_RQ_xlocs,                &
                    Caller_TQ_xlocs,                &
                    Caller_PQ_xlocs


USE Functions_Math, &
            ONLY :  Lagrange_Poly

USE Variables_AMReX_Core, &
           ONLY :  AMReX_Num_Levels

USE Poseidon_AMReX_MakeFineMask_Module, &
           ONLY :  AMReX_MakeFineMask

#endif


CONTAINS


#ifdef POSEIDON_AMREX_FLAG
 !+101+########################################################!
!                                                               !
!       IG_Input_XCFC_AMReX                                     !
!                                                               !
 !#############################################################!
SUBROUTINE IG_Input_XCFC_AMReX( NQ,                     &
                                RQ_Input,               &
                                TQ_Input,               &
                                PQ_Input,               &
                                Left_Limit,             &
                                Right_Limit,            &
                                nLevels,                &
                                MF_Guess                )


INTEGER,    DIMENSION(3),       INTENT(IN)      ::  NQ
REAL(idp),  DIMENSION(NQ(1)),   INTENT(IN)      ::  RQ_Input
REAL(idp),  DIMENSION(NQ(2)),   INTENT(IN)      ::  TQ_Input
REAL(idp),  DIMENSION(NQ(3)),   INTENT(IN)      ::  PQ_Input
REAL(idp),                      INTENT(IN)      ::  Left_Limit
REAL(idp),                      INTENT(IN)      ::  Right_Limit

INTEGER,                        INTENT(IN)      ::  nLevels
TYPE(amrex_multifab),           INTENT(IN)      ::  MF_Guess(0:nLevels-1)
                                        

REAL(idp)                                       ::  Quad_Span
REAL(idp),  DIMENSION(1:NQ(1))                  ::  Cur_RX_Locs
REAL(idp),  DIMENSION(1:NQ(2))                  ::  Cur_TX_Locs
REAL(idp),  DIMENSION(1:NQ(3))                  ::  Cur_PX_Locs
        
        
INTEGER,    DIMENSION(1:3)                      ::  nGhost_Vec


INTEGER                                         ::  lvl
TYPE(amrex_mfiter)                              ::  mfi
TYPE(amrex_box)                                 ::  Box
TYPE(amrex_imultifab)                           ::  Level_Mask

INTEGER                                         ::  nComp
INTEGER                                         ::  nQuad

INTEGER,    CONTIGUOUS, POINTER                 ::  Mask_PTR(:,:,:,:)
REAL(idp),  CONTIGUOUS, POINTER                 ::  Guess_PTR(:,:,:,:)

INTEGER                                         ::  re, te, pe
INTEGER, DIMENSION(3)                           ::  iEL, iEU, iE
    
    
nQuad     = NQ(1)*NQ(2)*NQ(3)
Quad_Span = Right_Limit - Left_Limit

Cur_RX_Locs = 2.0_idp * ( RQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp
Cur_TX_Locs = 2.0_idp * ( TQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp
Cur_PX_Locs = 2.0_idp * ( PQ_Input(:) - Left_Limit )/Quad_Span - 1.0_idp

nGhost_Vec = 0


DO lvl = nLevels-1,0,-1

    
    !
    !   MakeFineMask
    !
    IF ( lvl < AMReX_Num_Levels-1 ) THEN
        CALL AMReX_MakeFineMask(  Level_Mask,           &
                                  MF_Guess(lvl)%ba,     &
                                  MF_Guess(lvl)%dm,     &
                                  nGhost_Vec,           &
                                  MF_Guess(lvl+1)%ba,   &
                                  iLeaf, iTrunk         )
    ELSE
        ! Create Level_Mask all equal to 1
        CALL amrex_imultifab_build( Level_Mask,         &
                                    MF_Guess(lvl)%ba,   &
                                    MF_Guess(lvl)%dm,   &
                                    1,                  &
                                    nGhost_Vec(1)       )
        CALL Level_Mask%SetVal(iLeaf)
    END IF



    CALL amrex_mfiter_build(mfi, MF_Guess(lvl), tiling = .true. )

    DO WHILE(mfi%next())

        Guess_PTR => MF_Guess(lvl)%dataPtr(mfi)
        Mask_PTR  => Level_Mask%dataPtr(mfi)

        Box = mfi%tilebox()
        nComp =  MF_Guess(lvl)%ncomp()

        iEL = Box%lo
        iEU = Box%hi

        ! Fill Leaf Elements
        DO re = iEL(1),iEU(1)
        DO te = iEL(2),iEU(2)
        DO pe = iEL(3),iEU(3)

            IF ( Mask_PTR(RE,TE,PE,1) == iLeaf ) THEN

                iE = [re, te, pe]
                
                CALL SetGuess_From_Elem(iE, iEL, NQ,    &
                                        Cur_RX_Locs,    &
                                        Cur_TX_Locs,    &
                                        Cur_PX_Locs,    &
                                        lvl,            &
                                        Guess_PTR       )
                            
                
            END IF !  Mask_PTR(RE,TE,PE,1) == iLeaf

        END DO ! pe
        END DO ! te
        END DO ! re
        
    END DO

    CALL amrex_mfiter_destroy(mfi)
    CALL amrex_imultifab_destroy( Level_Mask )

END DO ! lvl


dVB_Coeff_Vector(:,iVB_X) = 0.0_idp

END SUBROUTINE IG_Input_XCFC_AMReX






!+101+###########################################################################!
!                                                                                !
!               Poseidon_Input_Guess                                             !
!                                                                                !
!################################################################################!
SUBROUTINE IG_Input_XCFC_AMReX_Caller( MF_Guess )


TYPE(amrex_multifab),       INTENT(INOUT)   ::  MF_Guess(0:AMReX_Num_Levels-1)



CALL IG_Input_XCFC_AMReX( Caller_NQ,                    &
                                Caller_RQ_xlocs,        &
                                Caller_TQ_xlocs,        &
                                Caller_PQ_xlocs,        &
                                Caller_xL(1),           &
                                Caller_xL(2),           &
                                AMReX_Num_Levels,       &
                                MF_Guess                )



END SUBROUTINE IG_Input_XCFC_AMReX_Caller




 !+201+########################################################!
!                                                               !
!          SetGuess_From_Elem                                   !
!                                                               !
 !#############################################################!
SUBROUTINE SetGuess_From_Elem(  iE, iEL, NQ,    &
                                Cur_RX_Locs,    &
                                Cur_TX_Locs,    &
                                Cur_PX_Locs,    &
                                lvl,            &
                                Guess_PTR       )
                                
INTEGER,    DIMENSION(1:3),         INTENT(IN)      ::  iE
INTEGER,    DIMENSION(1:3),         INTENT(IN)      ::  iEL
INTEGER,    DIMENSION(1:3),         INTENT(IN)      ::  NQ
REAL(idp),  DIMENSION(1:NQ(1)),     INTENT(IN)      ::  Cur_RX_Locs
REAL(idp),  DIMENSION(1:NQ(2)),     INTENT(IN)      ::  Cur_TX_Locs
REAL(idp),  DIMENSION(1:NQ(3)),     INTENT(IN)      ::  Cur_PX_Locs
INTEGER,                            INTENT(IN)      ::  lvl
REAL(idp),  CONTIGUOUS, POINTER,    INTENT(IN)      ::  Guess_PTR(:,:,:,:)

INTEGER                                             ::  d, lm
INTEGER                                             ::  iRE
INTEGER                                             ::  iU, iVB, iU_Offset
INTEGER                                             ::  Here, There, Somewhere
INTEGER                                             ::  nQuad
REAL(idp),      DIMENSION(0:NQ(1)-1)                ::  LagP
REAL(idp)                                           ::  sqrtfourpi

nQuad = NQ(1)*NQ(2)*NQ(3)
lm = 1
sqrtfourpi = sqrt(4.0_idp*pi)

! 1-D implementation
IF ((iE(2) == 1) .AND. (iE(3) == 1)) THEN


iRE = FEM_Elem_Map(iE(1),lvl)
IF ( iRE == 0 ) THEN
    d = 0
    LagP = Lagrange_Poly(FEM_Node_xlocs(d), NQ(1)-1, Cur_RX_Locs)

    DO iU = iU_CF,iU_LF
        Here  = (iU-1)*nQuad + 1
        There = (iU-1)*nQuad + NQ(1)
        Somewhere = Map_To_FEM_Node(iRE,d)
!        PRINT*,Guess_PTR(iE(1),iE(2),iE(3),Here:There),dVA_Coeff_Vector(Somewhere,lm,iU)
        dVA_Coeff_Vector(Somewhere,lm,iU) = sqrtfourpi*SUM( LagP(:)*Guess_PTR(iE(1),iE(2),iE(3),Here:There) )
!        PRINT*,Guess_PTR(iE(1),iE(2),iE(3),Here:There),dVA_Coeff_Vector(Somewhere,lm,iU)
    END DO ! iU
    
    iVB = iVB_S
    DO iU = iU_S1,iU_S3
        iU_Offset = iU-3*iVB+1
        Here  = (iU-1)*nQuad + 1
        There = (iU-1)*nQuad + NQ(1)
        Somewhere = FP_Array_Map_TypeB(iU,iVB,iRE,d,1)
        dVB_Coeff_Vector(Somewhere,iVB) = sqrtfourpi*SUM( LagP(:)*Guess_PTR(iE(1),iE(2),iE(3),Here:There) )
    END DO ! iU
END IF





! Use Guess Values to Interpolate to FEM Nodes

DO d = 1,Degree
    LagP = Lagrange_Poly(FEM_Node_xlocs(d), NQ(1)-1, Cur_RX_Locs)

    DO iU = iU_CF,iU_LF
        Here  = (iU-1)*nQuad + 1
        There = (iU-1)*nQuad + NQ(1)
        Somewhere = Map_To_FEM_Node(iRE,d)
        dVA_Coeff_Vector(Somewhere,lm,iU) = sqrtfourpi*SUM( LagP(:)*Guess_PTR(iE(1),iE(2),iE(3),Here:There) )
    END DO ! iU
    
    iVB = iVB_S
    DO iU = iU_S1,iU_S3
        iU_Offset = iU-3*iVB+1
        Here  = (iU-1)*nQuad + 1
        There = (iU-1)*nQuad + NQ(1)
        Somewhere = FP_Array_Map_TypeB(iU,iVB,iRE,d,1)
        dVB_Coeff_Vector(Somewhere,iVB) = sqrtfourpi*SUM( LagP(:)*Guess_PTR(iE(1),iE(2),iE(3),Here:There) )
    END DO ! iU
    
END DO

END IF



END SUBROUTINE SetGuess_From_Elem


#else


SUBROUTINE IG_Input_XCFC_AMReX()

END SUBROUTINE IG_Input_XCFC_AMReX



SUBROUTINE IG_Input_XCFC_AMReX_Caller( difference )

INTEGER,    INTENT(IN)          :: difference

END SUBROUTINE IG_Input_XCFC_AMReX_Caller

#endif




END MODULE IG_Input_AMReX_Module

