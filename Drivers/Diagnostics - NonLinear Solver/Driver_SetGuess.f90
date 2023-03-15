   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Driver_SetGuess_Module                                                !##!
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
            
USE Poseidon_Parameters, &
            ONLY :  Degree,                     &
                    L_Limit,                    &
                    Verbose_Flag
                    
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

USE Poseidon_Units_Module, &
            ONLY :  Grav_Constant_G,    &
                    Speed_of_Light,     &
                    C_Square,           &
                    GR_Source_Scalar,   &
                    Centimeter,         &
                    Second,             &
                    Millisecond,         &
                    Erg,                &
                    Gram

USE Poseidon_Interface_Initial_Guess, &
            ONLY :  Poseidon_Initialize_Flat_Guess
            
USE Variables_Mesh, &
            ONLY :  Num_R_Elements,             &
                    Num_T_Elements,             &
                    Num_P_Elements,             &
                    rlocs,                      &
                    drlocs,                     &
                    tlocs,                      &
                    plocs

USE Variables_Functions, &
            ONLY :  Potential_Solution
            
USE Variables_FEM_Module, &
            ONLY :  FEM_Node_xlocs
            
USE Variables_Vectors, &
            ONLY :  cVA_Coeff_Vector,   &
                    cVB_Coeff_Vector
                    
USE Maps_Domain, &
            ONLY :  Map_To_FEM_Node
            
USE Flags_Initial_Guess_Module, &
            ONLY :  lPF_IG_Flags,           &
                    iPF_IG_Set

IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!     Driver_SetGuess                                                			!
!                                                                               !
!###############################################################################!
SUBROUTINE Driver_SetGuess()


CALL Poseidon_Initialize_Flat_Guess()

!CALL SetGuess_Yahil_Newtonian()

END SUBROUTINE Driver_SetGuess




!+101+##########################################################################!
!                                                                               !
!     Driver_SetGuess                                                			!
!                                                                               !
!###############################################################################!
SUBROUTINE SetGuess_Yahil_Newtonian()


INTEGER                                 ::  RE
REAL(idp)                               ::  Potential
REAL(idp)                               ::  DROT
REAL(idp),  DIMENSION(0:Degree)         ::  Cur_R_Locs

REAL(idp)                               ::  CF_Val
REAL(idp)                               ::  LF_Val


INTEGER                                 ::  d
INTEGER                                 ::  Cur_Loc

PRINT*,"In SetGuess_Yahil_Newtonian"

cVA_Coeff_Vector = 0.0_idp
cVB_Coeff_Vector = 0.0_idp

DO RE = 0,Num_R_Elements-1

    DROT = (rlocs(re+1)-rlocs(re))/2.0_idp

    Cur_R_Locs(:) = DROT * (FEM_Node_xlocs(:)+1.0_idp) + rlocs(re)
    
    
    DO d = 0,Degree
        Potential = Potential_Solution(Cur_R_Locs(d)*Centimeter, 0.0_idp, 0.0_idp)
        CF_Val = 1.0_idp - 0.5_idp*Potential/C_Square
        LF_Val = 1.0_idp + 0.5_idp*Potential/C_Square
        
!        PRINT*,re,d,Cur_R_Locs(d),Potential,CF_Val,LF_Val
        Cur_Loc = Map_To_FEM_Node(re,d)

        cVA_Coeff_Vector(Cur_Loc,1,iU_CF) = SQRT(4.0*pi)*CF_Val
        cVA_Coeff_Vector(Cur_Loc,1,iU_LF) = SQRT(4.0*pi)*LF_Val

    END DO

END DO

!STOP "End of SetGuess_Yahil_Newtonian"

lPF_IG_Flags(iPF_IG_Set) = .TRUE.

END SUBROUTINE SetGuess_Yahil_Newtonian





END MODULE Driver_SetGuess_Module
