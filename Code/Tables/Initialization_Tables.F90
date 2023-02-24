   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Initialization_Tables                                                        !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +201+   Initialize_Ylm_Tables                                               !##!
!##!    +202+   Initialize_Lagrange_Poly_Tables                                     !##!
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
            ONLY :  idp

USE Poseidon_Numbers_Module, &
            ONLY :  pi, twopi

USE Poseidon_Message_Routines_Module, &
            ONLY :  Init_Message

USE Poseidon_Parameters, &
            ONLY :  Domain_Dim,             &
                    Degree,                 &
                    L_Limit,                &
                    Verbose_Flag

USE Variables_AMReX_Core, &
            ONLY :  AMReX_Max_Grid_Size,          &
                    AMReX_Num_Levels

USE Variables_Quadrature, &
            ONLY :  Num_R_Quad_Points,      &
                    Num_T_Quad_Points,      &
                    Num_P_Quad_Points,      &
                    Int_R_Locations,        &
                    Int_T_Locations,        &
                    Int_P_Locations

USE Variables_Mesh, &
            ONLY :  R_Inner,                &
                    R_Outer,                &
                    Num_R_Elements,         &
                    Num_T_Elements,         &
                    Num_P_Elements,         &
                    rlocs,                  &
                    tlocs,                  &
                    plocs

USE Variables_Derived, &
            ONLY :  LM_Length


USE Variables_Tables, &
            ONLY :  Lagrange_Poly_Table,    &
                    LPT_LPT,                &
                    M_Values,               &
                    Slm_Elem_Values,        &
                    Slm_Elem_dt_Values,     &
                    Slm_Elem_dp_Values,     &
                    Plm_Values,             &
                    Plm_dt_Values,          &
                    Nlm_Values,             &
                    Am_Values,              &
                    AM_dp_Values,           &
                    Level_dx,               &
                    Level_Ratios

USE Variables_FEM_Module, &
            ONLY :  FEM_Node_xlocs

USE Functions_Quadrature, &
            ONLY :  Initialize_LGL_Quadrature_Locations

USE Functions_Math, &
            ONLY :  Lagrange_Poly,          &
                    Lagrange_Poly_Deriv,    &
                    Legendre_Poly,          &
                    Legendre_Poly_Array,    &
                    Norm_Factor,            &
                    Sqrt_Factor

USE Maps_Quadrature, &
            ONLY :  Map_To_tpd

USE Maps_X_Space, &
            ONLY :  Map_From_X_Space

USE Maps_Domain, &
            ONLY :  Map_To_lm,              &
                    Map_To_Short_lm

USE Allocation_Tables, &
            ONLY :  Allocate_Tables
            
USE Initialization_Tables_Slm, &
            ONLY :  Initialize_Slm_Tables

USE Flags_Initialization_Module, &
            ONLY :  lPF_Init_Tables_Flags,  &
                    iPF_Init_Tables_Init

#ifdef POSEIDON_AMREX_FLAG
USE amrex_fort_module, &
            ONLY :  amrex_spacedim
#endif


IMPLICIT NONE

CONTAINS



!+101+##########################################################################!
!                                                                               !
!                       Initialize_Tables                                       !
!                                                                               !
!###############################################################################!
SUBROUTINE Initialize_Tables()



IF ( Verbose_Flag ) CALL Init_Message('Initializing Basis Functions Tables.')



CALL Allocate_Tables()
CALL Initialize_Lagrange_Poly_Tables( Degree, Num_R_Quad_Points )
CALL Initialize_Slm_Tables()
    
    
#ifdef POSEIDON_AMREX_FLAG
    CALL Initialize_Level_Tables()
#endif

lPF_Init_Tables_Flags(iPF_Init_Tables_Init) = .TRUE.


END SUBROUTINE Initialize_Tables






!+202+##########################################################################!
!                                                                               !
!                  Initialize_Lagrange_Poly_Tables                              !
!                                                                               !
!###############################################################################!
SUBROUTINE Initialize_Lagrange_Poly_Tables( Ord, Num_Quad_Points )

INTEGER, INTENT(IN)                         ::  Ord
INTEGER, INTENT(IN)                         ::  Num_Quad_Points


REAL(idp), DIMENSION(0:Ord)          ::  Lagrange_Poly_Values
REAL(idp), DIMENSION(0:Ord)          ::  Lagrange_DRV_Values

INTEGER                                     ::  Eval_Point
INTEGER                                     ::  rd, d, dp,dd



Lagrange_Poly_Table = 0.0_idp

DO Eval_Point = 1,Num_Quad_Points
    
    Lagrange_Poly_Values = Lagrange_Poly(INT_R_Locations(Eval_Point), Ord, FEM_Node_xlocs)
    Lagrange_DRV_Values  = Lagrange_Poly_Deriv(INT_R_Locations(Eval_Point), Ord, FEM_Node_xlocs)

    Lagrange_Poly_Table( :, Eval_Point, 0) = Lagrange_Poly_Values
    Lagrange_Poly_Table( :, Eval_Point, 1) = Lagrange_DRV_Values

END DO


DO dd = 0,1
    DO d = 0,Ord
        DO dp = 0,Ord
            DO rd = 1,Num_Quad_Points
           
                LPT_LPT(rd, d, dp, 0:1, dd) = Lagrange_Poly_Table(d, rd, 0:1)       &
                                            * Lagrange_Poly_Table(dp, rd, dd)
            END DO
        END DO
    END DO
END DO






END SUBROUTINE Initialize_Lagrange_Poly_Tables







 !+203+############################################################!
!                                                                   !
!          Initialize_Level_dx_Table                                !
!                                                                   !
 !#################################################################!
SUBROUTINE Initialize_Level_Tables()

INTEGER                 ::  lvl
REAL(idp)               ::  dr, dt, dp

#ifndef POSEIDON_AMREX_FLAG
INTEGER                 ::  amrex_spacedim = 3
#endif


dr = (R_Outer-R_Inner)/Num_R_Elements
dt = pi/Num_T_Elements
dp = TwoPi/Num_P_Elements


DO lvl = 0,AMReX_Num_Levels-1
    Level_dx(lvl,1) = dr/2.0_idp**lvl
END DO



IF ( amrex_spacedim < 2 ) THEN
    Level_dx(:,2) = dt
ELSE
    DO lvl = 0,AMReX_Num_Levels-1
        Level_dx(lvl,2) = dt/2.0_idp**lvl
    END DO
END IF



IF ( amrex_spacedim < 3 ) THEN
    Level_dx(:,3) = dp
ELSE
    DO lvl = 0,AMReX_Num_Levels-1
        Level_dx(lvl,3) = dp/2.0_idp**lvl
    END DO
END IF


! Level Ratios !
DO lvl = 0,AMReX_Num_Levels-1
    Level_Ratios(lvl) = 2**lvl
END DO



END SUBROUTINE Initialize_Level_Tables
















END MODULE Initialization_Tables
