   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE New_PETSc_Module                                                             !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!================================================================================!##!
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


USE OMP_LIB

USE Poseidon_Constants_Module, &
                        ONLY :  idp,                &
                                pi,                 &
                                TwoPi,              &
                                OneThird,           &
                                TwoThirds,          &
                                FourThirds,         &
                                OneThirtySecond

USE Poseidon_Parameters, &
                        ONLY :  DOMAIN_DIM,                 &
                                DEGREE,                     &
                                L_LIMIT,                    &
                                nPROCS_POSEIDON,            &
                                NUM_R_ELEMS_PER_SHELL,      &
                                NUM_R_ELEMS_PER_SUBSHELL,   &
                                NUM_R_ELEMS_PER_BLOCK,      &
                                NUM_T_ELEMS_PER_BLOCK,      &
                                NUM_P_ELEMS_PER_BLOCK,      &
                                NUM_R_QUAD_POINTS,          &
                                NUM_T_QUAD_POINTS,          &
                                NUM_P_QUAD_POINTS,          &
                                NUM_BLOCK_THETA_ROWS



USE Poseidon_Variables_Module, &
                        ONLY :  NUM_R_ELEMENTS,             &
                                NUM_T_ELEMENTS,             &
                                NUM_P_ELEMENTS,             &
                                NUM_R_NODES,                &
                                INT_R_LOCATIONS,            &
                                INT_T_LOCATIONS,            &
                                INT_P_LOCATIONS,            &
                                INT_R_WEIGHTS,              &
                                INT_T_WEIGHTS,              &
                                INT_P_WEIGHTS,              &
                                rlocs,                      &
                                tlocs,                      &
                                plocs,                      &
                                myShell,                    &
                                myID_Shell,                 &
                                myID_PETSC,                 &
                                SUBSHELL_PROB_DIM,          &
                                ULM_LENGTH,                 &
                                Coefficient_Vector


USE Jacobian_Internal_Functions_Module,  &
                        ONLY :  Initialize_Guess_Values,        &
                                Initialize_Special_Guess_Values



USE CFA_3D_Master_Build_Module,  &
                        ONLY :  Calc_3D_Current_Values,             &
                                Finish_3D_RHS_Vector



#include <petsc/finclude/petscsnes.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscpc.h>


use petscsnes
use petscdm
use petscksp
use petscsys
use petscmat
use petscvec
use petscpc


IMPLICIT NONE




CONTAINS


!+101+###########################################################################!
!                                                                                !
!                  New_PETSc_Routine                                             !
!                                                                                !
!################################################################################!
SUBROUTINE New_PETSc_Routine()


PetscErrorCode                                  ::  jerr
DM                                              ::  da
SNES                                            ::  snes
PC                                              ::  pc
KSP                                             ::  ksp

Vec                                             ::  Coeffs_Vector, Residual_Vector
Mat                                             ::  Jacobian_Mat
PetscInt                                        ::  Radial_Nodes, Stencil_Width, DOF, CTX
PetscInt, DIMENSION(nPROCS_POSEIDON)            ::  lx

PetscInt                                        ::   Start_here, Mid_Here, End_here
PetscInt, DIMENSION(0:SUBSHELL_PROB_DIM-1)      ::   Index_Array
PetscInt                                        ::   i

PetscScalar, DIMENSION(0:SUBSHELL_PROB_DIM-1)   ::   Test_Vec




Radial_Nodes = NUM_R_NODES
Stencil_Width = 1
DOF = 5*(L_LIMIT + 1)*(L_LIMIT + 1)



! Initalize PETSc !
CALL PetscInitialize(PETSC_NULL_CHARACTER, jerr)
IF (jerr .ne. 0 ) THEN
    PRINT*,"Unable to Initialize PETSc"
    STOP
END IF



! Create Nonlinear Solver Context !
CALL SNESCreate( PETSC_COMM_WORLD, snes, jerr )



!  Create DMDA  !
lx(1:nPROCS_POSEIDON)   = (Radial_Nodes - 1)/nPROCS_POSEIDON
lx(nPROCS_POSEIDON)     = lx(nPROCS_POSEIDON) + 1

CALL DMDACreate1d(  PETSC_COMM_WORLD,           &
                    DM_BOUNDARY_NONE,           &
                    Radial_Nodes,               &
                    DOF,                        &
                    Stencil_Width,              &
                    lx,                         &
                    da,                         &
                    jerr                        )

CALL DMSetFromOptions( da, jerr )
CALL DMSetUp( da, jerr )


!  Create Solution and Residual Vectors  !
CALL DMCreateGlobalVector( da, Coeffs_Vector, jerr )
CALL VecDuplicate( Coeffs_Vector, Residual_Vector,jerr )

!  Create Jacobian Matrix  !
CALL DMSetMatType( da, MATSEQAIJ, jerr )
CALL DMCreateMatrix( da, Jacobian_Mat, jerr )






! Set Function/Residual Evaluation Routine
CALL SNESSetFunction(snes, Residual_Vector, Residual_Function, CTX, jerr)

! Set Jacobian Matrix Routine
CALL SNESSetJacobian(snes, Jacobian_Mat, Jacobian_Mat, Jacobian_Function, CTX, jerr)





! Create Linear Solver and Preconditioner for SNES
CALL SNESGetKSP( snes, ksp, jerr )
CALL KSPGetPC( ksp, pc, jerr )

CALL SNESSetFromOptions( snes, jerr )








!                   !
! Set Initial Guess !
!                   !
!CALL Initialize_Guess_Values()
CALL Initialize_Special_Guess_Values()

! Crete Index_Array specifying where values go in the PETSc Vector
start_here=myID_PETSc*NUM_R_ELEMS_PER_SUBSHELL*DEGREE*ULM_LENGTH
end_here=start_here + SUBSHELL_PROB_DIM - 1
Index_Array =  (/(i,i=start_here,end_here,1)/)

! Define bounds within Block_RHS_VECTOR where values going into PETSc exist
start_here=myID_SHELL*NUM_R_ELEMS_PER_SUBSHELL*DEGREE*ULM_LENGTH
end_here=start_here + SUBSHELL_PROB_DIM - 1



CALL VecSetValues(Coeffs_Vector, SUBSHELL_PROB_DIM, Index_Array,              &
                  Coefficient_Vector(start_here:end_here), INSERT_VALUES, jerr)

CALL VecAssemblyBegin(Coeffs_Vector, jerr)
CALL VecAssemblyEnd(Coeffs_Vector, jerr)


!CALL VecGetValues(Coeffs_Vector, SUBSHELL_PROB_DIM, Index_Array, Test_Vec, jerr)








! Have PETSc Solve the Non-linear System !

!CALL SNESSolve( snes, PETSC_NULL_OBJECT, Coeffs, jerr )




CALL VecDestroy( Coeffs_Vector, jerr )
CALL VecDestroy( Residual_Vector, jerr )
CALL MatDestroy( Jacobian_Mat, jerr )
CALL SNESDestroy( snes, jerr )



CALL PetscFinalize( jerr )


END SUBROUTINE New_PETSc_Routine


















!+201+###########################################################################!
!                                                                                !
!                  Residual_Function                                             !
!                                                                                !
!################################################################################!
SUBROUTINE Residual_Function( snes, Coeffs, Residual, CTX, jerr )


SNES                        ::  snes
Vec                         ::  Coeffs, Residual
PetscInt                    ::  CTX
PetscErrorCode              ::  jerr


INTEGER                                                                 ::  Local_re,       &
                                                                            Local_te,       &
                                                                            Local_pe,       &
                                                                            Global_re,      &
                                                                            Global_te,      &
                                                                            Global_pe,      &
                                                                            d, l, m,        &
                                                                            dp, lp, mp,     &
                                                                            rd, td, pd

REAL(KIND = idp),  DIMENSION(1:5,                   &
                             1:NUM_R_QUAD_POINTS,   &
                             1:NUM_T_QUAD_POINTS,   &
                             1:NUM_P_QUAD_POINTS    )                   ::  RHS_TERMS

REAL(KIND = idp), DIMENSION(1:NUM_R_QUAD_POINTS, 1:NUM_T_QUAD_POINTS)   ::  RSIN_SQUARE

REAL(KIND = idp), DIMENSION(1:NUM_R_QUAD_POINTS)                        ::  R_VAL,              &
                                                                            R_SQUARE,           &
                                                                            R_CUBED,             &
                                                                            R_INVERSE


REAL(KIND = idp), DIMENSION(1:NUM_T_QUAD_POINTS)                        ::  SIN_VAL,            &
                                                                            SIN_SQUARE,         &
                                                                            CSC_VAL,            &
                                                                            CSC_SQUARE,         &
                                                                            COTAN_VAL,          &
                                                                            COS_VAL,            &
                                                                            COS_SQUARE


COMPLEX(KIND = idp), DIMENSION(1:NUM_P_QUAD_POINTS)                     ::  PHI_EXP,            &
                                                                            PHI_TWOEXP



REAL(KIND = idp)                                                        ::  TWOOVER_DELTAR,    &
                                                                            deltar_overtwo,     &
                                                                            deltat_overtwo,     &
                                                                            deltap_overtwo



REAL(KIND = idp), DIMENSION( 1:NUM_R_QUAD_POINTS,        &
                             1:NUM_T_QUAD_POINTS,        &
                             1:NUM_P_QUAD_POINTS         )              ::  Int_Factor



REAL(KIND = idp), DIMENSION( 1:NUM_R_QUAD_POINTS,   &
                             1:NUM_T_QUAD_POINTS,   &
                             1:NUM_P_QUAD_POINTS)               ::  CUR_VAL_PSI

REAL(KIND = idp), DIMENSION( 1:NUM_R_QUAD_POINTS,   &
                             1:NUM_T_QUAD_POINTS,   &
                             1:NUM_P_QUAD_POINTS)               ::  CUR_VAL_ALPHAPSI

REAL(KIND = idp), DIMENSION( 1:3,                   &
                             1:NUM_R_QUAD_POINTS,   &
                             1:NUM_T_QUAD_POINTS,   &
                             1:NUM_P_QUAD_POINTS )              ::  CUR_VAL_BETA

REAL(KIND = idp), DIMENSION( 1:3,                   &
                             1:NUM_R_QUAD_POINTS,   &
                             1:NUM_T_QUAD_POINTS,   &
                             1:NUM_P_QUAD_POINTS    )           ::  CUR_DRV_PSI

REAL(KIND = idp), DIMENSION( 1:3,                   &
                             1:NUM_R_QUAD_POINTS,   &
                             1:NUM_T_QUAD_POINTS,   &
                             1:NUM_P_QUAD_POINTS    )           ::  CUR_DRV_ALPHAPSI

REAL(KIND = idp), DIMENSION( 1:3, 1:3,              &
                             1:NUM_R_QUAD_POINTS,   &
                             1:NUM_T_QUAD_POINTS,   &
                             1:NUM_P_QUAD_POINTS    )           ::  CUR_DRV_BETA

REAL(KIND = idp), DIMENSION( 1:3, 1:3, 1:3,         &
                             1:NUM_R_QUAD_POINTS,   &
                             1:NUM_T_QUAD_POINTS,   &
                             1:NUM_P_QUAD_POINTS    )           ::  CUR_DDRV_BETA


REAL(KIND = idp), DIMENSION(1:NUM_R_QUAD_POINTS)                ::  CUR_R_LOCS
REAL(KIND = idp), DIMENSION(1:NUM_T_QUAD_POINTS)                ::  CUR_T_LOCS
REAL(KIND = idp), DIMENSION(1:NUM_P_QUAD_POINTS)                ::  CUR_P_LOCS

INTEGER                                                         ::  Block_T_Begin, Block_P_Begin








Block_T_Begin = MOD(myID_Shell,NUM_BLOCK_THETA_ROWS)*NUM_T_ELEMS_PER_BLOCK
Block_P_Begin = (myID_Shell/NUM_BLOCK_THETA_ROWS)*NUM_P_ELEMS_PER_BLOCK


DO Local_re = 0,NUM_R_ELEMS_PER_BLOCK-1


    Global_re = myShell*NUM_R_ELEMS_PER_SHELL + Local_re


    deltar_overtwo = (rlocs(Global_re + 1) - rlocs(Global_re))/2.0_idp
    TWOOVER_DELTAR = 1.0_idp/deltar_overtwo
    CUR_R_LOCS(:) = deltar_overtwo * (INT_R_LOCATIONS(:)+1.0_idp) + rlocs(Global_re)


    R_SQUARE(:) = CUR_R_LOCS(:)*CUR_R_LOCS(:)
    R_CUBED(:) = R_SQUARE(:)*CUR_R_LOCS(:)


    R_INVERSE(:) = 1.0_idp/CUR_R_LOCS(:)


    !
    !   Move through phi
    !
    DO Local_pe = 0, NUM_P_ELEMS_PER_BLOCK-1



        Global_pe = Block_P_Begin + Local_pe

        deltap_overtwo = (plocs(Global_pe + 1) - plocs(Global_pe))/2.0_idp
        CUR_P_LOCS(:) = deltap_overtwo * (INT_P_LOCATIONS(:)+1.0_idp) + plocs(Global_pe)

        PHI_EXP(:) = EXP( CMPLX(0, -CUR_P_LOCS(:), KIND = idp) )
        PHI_TWOEXP(:) = EXP( CMPLX(0, -2.0_idp*CUR_P_LOCS(:), KIND = idp) )


        DO Local_te = 0, NUM_T_ELEMS_PER_BLOCK-1


            Global_te = Block_T_Begin + Local_te

            deltat_overtwo = (tlocs(Global_te + 1) - tlocs(Global_te))/2.0_idp
            CUR_T_LOCS(:) = deltat_overtwo * (INT_T_LOCATIONS(:)+1.0_idp) + tlocs(Global_te)


            COTAN_VAL(:) = 1.0_idp/DTAN(CUR_T_LOCS(:))
            CSC_VAL(:) = 1.0_idp/DSIN(CUR_T_LOCS(:))
            SIN_VAL(:) = DSIN(CUR_T_LOCS(:))
            COS_VAL(:) = DCOS(CUR_T_LOCS(:))

            COS_SQUARE(:) = COS_VAL(:)*COS_VAL(:)
            SIN_SQUARE(:) = SIN_VAL(:)*SIN_VAL(:)
            CSC_SQUARE(:) = CSC_VAL(:)*CSC_VAL(:)

            DO td = 1,NUM_T_QUAD_POINTS

                RSIN_SQUARE(:,td) = R_SQUARE(:)*SIN_SQUARE(td)

            END DO




            !*!
            !*! Calculate Current Values of CFA Varaiables and their Deriviatives
            !*!
            CALL Calc_3D_Current_Values(Global_re , Local_te  , Local_pe,       &
                                        CUR_VAL_PSI, CUR_DRV_PSI,                       &
                                        CUR_VAL_ALPHAPSI, CUR_DRV_ALPHAPSI,             &
                                        CUR_VAL_BETA, CUR_DRV_BETA, CUR_DDRV_BETA,      &
                                        Int_Factor,                                     &
                                        CUR_R_LOCS, CUR_T_LOCS, CUR_P_LOCS,             &
                                        R_SQUARE,                                       &
                                        SIN_VAL,                                        &
                                        DELTAR_OVERTWO, DELTAT_OVERTWO, DELTAP_OVERTWO  )


            !*!
            !*!  Calculate the Sub-Jacobian and RHS Terms
            !*!
!            CALL Create_Resdiual_Vector( Local_re, Local_te, Local_pe,                   &
!                                         CUR_VAL_PSI, CUR_DRV_PSI,                       &
!                                         CUR_VAL_ALPHAPSI, CUR_DRV_ALPHAPSI,             &
!                                         CUR_VAL_BETA, CUR_DRV_BETA, CUR_DDRV_BETA,      &
!                                         Int_Factor,                                     &
!                                         RHS_TERMS,                                      &
!                                         CUR_R_LOCS, CUR_T_LOCS, CUR_P_LOCS,             &
!                                         R_SQUARE, R_CUBED, R_INVERSE,                   &
!                                         SIN_VAL, COS_VAL, CSC_VAL, COTAN_VAL,           &
!                                         SIN_SQUARE, CSC_SQUARE, RSIN_SQUARE             )




        END DO  ! pe Loop
    END DO  ! te Loop
END DO  ! re Loop



CALL FINISH_3D_RHS_VECTOR()










END SUBROUTINE Residual_Function


















!+202+###########################################################################!
!                                                                                !
!                  Jacobian_Function                                             !
!                                                                                !
!################################################################################!
SUBROUTINE Jacobian_Function( snes, Coeffs, Jacobian_Mat, PMat, CTX, jerr )


SNES                        ::  snes
Vec                         ::  Coeffs
Mat                         ::  Jacobian_Mat, PMat
PetscInt                    ::  CTX
PetscErrorCode              ::  jerr


INTEGER                                                                 ::  Local_re,       &
                                                                            Local_te,       &
                                                                            Local_pe,       &
                                                                            Global_re,      &
                                                                            Global_te,      &
                                                                            Global_pe,      &
                                                                            d, l, m,        &
                                                                            dp, lp, mp,     &
                                                                            rd, td, pd

REAL(KIND = idp),  DIMENSION(1:5,                   &
                             1:NUM_R_QUAD_POINTS,   &
                             1:NUM_T_QUAD_POINTS,   &
                             1:NUM_P_QUAD_POINTS    )                   ::  RHS_TERMS

REAL(KIND = idp), DIMENSION(1:NUM_R_QUAD_POINTS, 1:NUM_T_QUAD_POINTS)   ::  RSIN_SQUARE

REAL(KIND = idp), DIMENSION(1:NUM_R_QUAD_POINTS)                        ::  R_VAL,              &
                                                                            R_SQUARE,           &
                                                                            R_CUBED,             &
                                                                            R_INVERSE


REAL(KIND = idp), DIMENSION(1:NUM_T_QUAD_POINTS)                        ::  SIN_VAL,            &
                                                                            SIN_SQUARE,         &
                                                                            CSC_VAL,            &
                                                                            CSC_SQUARE,         &
                                                                            COTAN_VAL,          &
                                                                            COS_VAL,            &
                                                                            COS_SQUARE


COMPLEX(KIND = idp), DIMENSION(1:NUM_P_QUAD_POINTS)                     ::  PHI_EXP,            &
                                                                            PHI_TWOEXP



REAL(KIND = idp)                                                        ::  TWOOVER_DELTAR,    &
                                                                            deltar_overtwo,     &
                                                                            deltat_overtwo,     &
                                                                            deltap_overtwo



REAL(KIND = idp), DIMENSION( 1:NUM_R_QUAD_POINTS,        &
                             1:NUM_T_QUAD_POINTS,        &
                             1:NUM_P_QUAD_POINTS         )              ::  Int_Factor



REAL(KIND = idp), DIMENSION( 1:NUM_R_QUAD_POINTS,   &
                             1:NUM_T_QUAD_POINTS,   &
                             1:NUM_P_QUAD_POINTS)               ::  CUR_VAL_PSI

REAL(KIND = idp), DIMENSION( 1:NUM_R_QUAD_POINTS,   &
                             1:NUM_T_QUAD_POINTS,   &
                             1:NUM_P_QUAD_POINTS)               ::  CUR_VAL_ALPHAPSI

REAL(KIND = idp), DIMENSION( 1:3,                   &
                             1:NUM_R_QUAD_POINTS,   &
                             1:NUM_T_QUAD_POINTS,   &
                             1:NUM_P_QUAD_POINTS )              ::  CUR_VAL_BETA

REAL(KIND = idp), DIMENSION( 1:3,                   &
                             1:NUM_R_QUAD_POINTS,   &
                             1:NUM_T_QUAD_POINTS,   &
                             1:NUM_P_QUAD_POINTS    )           ::  CUR_DRV_PSI

REAL(KIND = idp), DIMENSION( 1:3,                   &
                             1:NUM_R_QUAD_POINTS,   &
                             1:NUM_T_QUAD_POINTS,   &
                             1:NUM_P_QUAD_POINTS    )           ::  CUR_DRV_ALPHAPSI

REAL(KIND = idp), DIMENSION( 1:3, 1:3,              &
                             1:NUM_R_QUAD_POINTS,   &
                             1:NUM_T_QUAD_POINTS,   &
                             1:NUM_P_QUAD_POINTS    )           ::  CUR_DRV_BETA

REAL(KIND = idp), DIMENSION( 1:3, 1:3, 1:3,         &
                             1:NUM_R_QUAD_POINTS,   &
                             1:NUM_T_QUAD_POINTS,   &
                             1:NUM_P_QUAD_POINTS    )           ::  CUR_DDRV_BETA


REAL(KIND = idp), DIMENSION(1:NUM_R_QUAD_POINTS)                ::  CUR_R_LOCS
REAL(KIND = idp), DIMENSION(1:NUM_T_QUAD_POINTS)                ::  CUR_T_LOCS
REAL(KIND = idp), DIMENSION(1:NUM_P_QUAD_POINTS)                ::  CUR_P_LOCS

INTEGER                                                         ::  Block_T_Begin, Block_P_Begin











END SUBROUTINE Jacobian_Function





END MODULE New_PETSc_Module
