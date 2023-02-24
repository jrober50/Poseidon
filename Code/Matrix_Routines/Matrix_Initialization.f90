   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Matrix_Initialization_Module                                                 !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
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
            ONLY :  idp

USE Poseidon_Bailout_Module, &
            ONLY :  Poseidon_Bailout
            
USE Poseidon_Numbers_Module, &
            ONLY : pi, TwoPi

USE Poseidon_Message_Routines_Module, &
            ONLY :  Init_Message

USE Poseidon_Parameters, &
            ONLY :  DEGREE,                     &
                    L_LIMIT,                    &
                    Verbose_Flag,               &
                    Domain_Dim

USE Variables_Quadrature, &
            ONLY :  Num_R_Quad_Points

USE Variables_Mesh, &
            ONLY :  Num_R_Elements,             &
                    rlocs

USE Variables_Derived, &
            ONLY :  Num_R_Nodes,                &
                    ULM_Length,                 &
                    LM_Length,                  &
                    LM_Short_Length,            &
                    Var_Dim


USE Variables_Matrices, &
            ONLY :  Matrix_Format,              &
                    Laplace_Matrix_Full,        &
                    Laplace_Matrix_VAL,         &
                    Laplace_Matrix_ROW,         &
                    Laplace_Matrix_COL,         &
                    Laplace_Factored_VAL,       &
                    Laplace_Factored_ROW,       &
                    Laplace_Factored_COL,       &
                    dMB_Matrix_Banded,            &
                    iMB_Bandwidth

USE Variables_FEM_Module, &
            ONLY :  FEM_Node_xlocs

USE Matrix_Cholesky_Factorization_Module,   &
            ONLY :  Cholesky_Factorization

USE Maps_Fixed_Point, &
            ONLY :  FP_Beta_Array_Map

USE Maps_Domain, &
            ONLY :  Map_To_FEM_Node,            &
                    Map_To_LM,                  &
                    Map_To_Short_LM

USE Maps_X_Space, &
            ONLY :  Map_From_X_Space

USE Functions_Math, &
            ONLY :  Lagrange_Poly,              &
                    Lagrange_Poly_Deriv,        &
                    Legendre_Poly,              &
                    Norm_Factor

USE Functions_Quadrature, &
            ONLY :  Initialize_LGL_Quadrature_Locations,    &
                    Initialize_LG_Quadrature,               &
                    Initialize_Trapezoid_Quadrature

USE Initialization_Tables_Slm, &
            ONLY :  Initialize_Nlm_Table,           &
                    Initialize_Am_Tables,           &
                    Initialize_Plm_Tables,          &
                    Initialize_Slm_Tables_on_Elem


USE Flags_Initialization_Module, &
            ONLY :  lPF_Init_Matrices_Flags,            &
                    iPF_Init_Matrices_Type_A,           &
                    iPF_Init_Matrices_Type_B,           &
                    iPF_Init_Matrices_Type_A_Cholesky,  &
                    iPF_Init_Matrices_Type_B_LU


USE Timer_Routines_Module, &
            ONLY :  TimerStart,                     &
                    TimerStop

USE Timer_Variables_Module, &
            ONLY :  Timer_Matrix_Init,              &
                    Timer_Matrix_Radial_Terms,      &
                    Timer_Matrix_Angular_Terms,     &
                    Timer_Matrix_Laplace_Init,      &
                    Timer_Matrix_MVL_Init

USE MPI

IMPLICIT NONE


INTEGER,                                            PRIVATE     ::  Int_R_Deg
INTEGER,                                            PRIVATE     ::  Int_T_Deg
INTEGER,                                            PRIVATE     ::  Int_P_Deg
INTEGER,                                            PRIVATE     ::  Int_TP_Deg


REAL(idp),  ALLOCATABLE,    DIMENSION(:),           PRIVATE     ::  Cur_T_Locs
REAL(idp),  ALLOCATABLE,    DIMENSION(:),           PRIVATE     ::  Cur_P_Locs

REAL(idp),  ALLOCATABLE,    DIMENSION(:),           PRIVATE     ::  Int_R_Locs, Int_R_Weights
REAL(idp),  ALLOCATABLE,    DIMENSION(:),           PRIVATE     ::  Int_T_Locs, Int_T_Weights
REAL(idp),  ALLOCATABLE,    DIMENSION(:),           PRIVATE     ::  Int_P_Locs, Int_P_Weights

REAL(idp),  ALLOCATABLE,    DIMENSION( :,:,:,:,: ), PRIVATE     ::  LP_LP_Table
REAL(idp),  ALLOCATABLE,    DIMENSION( :,:,: ),     PRIVATE     ::  TP_TP_Integrals

REAL(idp),  ALLOCATABLE,    DIMENSION( :,:,: ),     PRIVATE     ::  RR_Factor
REAL(idp),  ALLOCATABLE,    DIMENSION( :,:,: ),     PRIVATE     ::  DRR_Factor
REAL(idp),  ALLOCATABLE,    DIMENSION( :,:,: ),     PRIVATE     ::  RDR_Factor
REAL(idp),  ALLOCATABLE,    DIMENSION( :,:,: ),     PRIVATE     ::  DRDR_Factor

REAL(idp),  ALLOCATABLE,    DIMENSION(:,:),         PRIVATE     ::  Slm_Table
REAL(idp),  ALLOCATABLE,    DIMENSION(:,:),         PRIVATE     ::  Slm_dp_Table
REAL(idp),  ALLOCATABLE,    DIMENSION(:,:),         PRIVATE     ::  Slm_dt_Table
REAL(idp),  ALLOCATABLE,    DIMENSION(:,:,:),       PRIVATE     ::  Plm_Table
REAL(idp),  ALLOCATABLE,    DIMENSION(:,:,:),       PRIVATE     ::  Plm_dt_Table
REAL(idp),  ALLOCATABLE,    DIMENSION(:,:,:),       PRIVATE     ::  Am_Table
REAL(idp),  ALLOCATABLE,    DIMENSION(:,:,:),       PRIVATE     ::  Am_dp_Table
REAL(idp),  ALLOCATABLE,    DIMENSION(:),           PRIVATE     ::  Nlm_Table
REAL(idp),  ALLOCATABLE,    DIMENSION(:),           PRIVATE     ::  TP_Int_Weights

CONTAINS


!+101+##########################################################################!
!                                                                               !
!           Initialize_XCFC_Matrices                                            !
!                                                                               !
!###############################################################################!
SUBROUTINE Initialize_XCFC_Matrices()


IF ( Verbose_Flag ) CALL Init_Message('Beginning Matrix Initialization.')
CALL TimerStart( Timer_Matrix_Init )


! Set number of quadrature points for different dimensions
!
! Radial : Using Gaussian Quadrature. 2n-1 points for exact integration of Polynomials of degree n.
!           Radially this is a product of two Lagrange Polys so n = 2*Degree.
! Phi    : Using Trapezoid rule quadrature. m + 1 for exact integration of exp(-i*m*phi).
!           Max m is defined by the product of exp(i*L_Limit*Phi)*exp(i*L_Limit*Phi) = exp(i*(2*L_Limit)*Phi)
!
Int_R_Deg  = 4*Degree - 1
!INT_R_Deg  = 5
Int_T_Deg  = 15
Int_P_Deg  = 2*L_LIMIT + 1

Int_TP_Deg = Int_T_Deg*Int_P_Deg



CALL Allocate_Matrix_Init_Variables()


CALL Initialize_LG_Quadrature(Int_R_Deg, Int_R_Locs, Int_R_Weights)
CALL Initialize_LG_Quadrature(Int_T_Deg, Int_T_Locs, Int_T_Weights)
!CALL Initialize_LG_Quadrature(Int_P_Deg, Int_P_Locs, Int_P_Weights)
CALL Initialize_Trapezoid_Quadrature(Int_P_Deg, 1, Int_P_Locs, Int_P_Weights)



 !                                                          !
!!   Map Theta and Phi Locations from [-1,1] space to real space.   !!
 !                                                          !
Cur_T_Locs = Map_From_X_Space( 0.0_idp, pi, Int_T_Locs )
Cur_P_Locs = Int_P_Locs

CALL Calculate_Radial_Terms()
CALL Calculate_Angular_Terms()

CALL Calculate_Laplace_Matrix()
CALL Calculate_MVL_Banded()

Call Deallocate_Matrix_Init_Variables()


CALL TimerStop( Timer_Matrix_Init )




END SUBROUTINE Initialize_XCFC_Matrices



!+101+##########################################################################!
!                                                                               !
!           Initialize_XCFC_Matrices                                            !
!                                                                               !
!###############################################################################!
SUBROUTINE Initialize_Poisson_Matrix()


IF ( Verbose_Flag ) CALL Init_Message('Beginning Matrix Initialization.')


! Set number of quadrature points for different dimensions
!
! Radial : Using Gaussian Quadrature. 2n-1 points for exact integration of Polynomials of degree n.
!           Radially this is a product of two Lagrange Polys so n = 2*Degree.
! Phi    : Using Trapezoid rule quadrature. m + 1 for exact integration of exp(-i*m*phi).
!           Max m is defined by the product of exp(i*L_Limit*Phi)*exp(i*L_Limit*Phi) = exp(i*(2*L_Limit)*Phi)
!
Int_R_Deg  = 4*Degree - 1
!INT_R_Deg  = 5
Int_T_Deg  = 15
Int_P_Deg  = 2*L_LIMIT + 1

Int_TP_Deg = Int_T_Deg*Int_P_Deg



CALL Allocate_Matrix_Init_Variables()


CALL Initialize_LG_Quadrature(Int_R_Deg, Int_R_Locs, Int_R_Weights)
!CALL Initialize_LG_Quadrature(Int_T_Deg, Int_T_Locs, Int_T_Weights)
!CALL Initialize_Trapezoid_Quadrature(Int_P_Deg, 1, Int_P_Locs, Int_P_Weights)



 !                                                          !
!!   Map Theta and Phi Locations from [-1,1] space to real space.   !!
 !                                                          !
Cur_T_Locs = Map_From_X_Space( 0.0_idp, pi, Int_T_Locs )
Cur_P_Locs = Int_P_Locs

CALL Calculate_Radial_Terms()
!CALL Calculate_Angular_Terms()


CALL Calculate_Laplace_Matrix()


Call Deallocate_Matrix_Init_Variables()


END SUBROUTINE Initialize_Poisson_Matrix



 !+201+################################################################!
!                                                                       !
!          Allocate_Matrix_Init_Variables                               !
!                                                                       !
 !#####################################################################!
SUBROUTINE Allocate_Matrix_Init_Variables()

ALLOCATE( Cur_T_Locs(1:Int_T_Deg) )
ALLOCATE( Cur_P_Locs(1:Int_P_Deg) )

ALLOCATE( Int_R_Locs(1:Int_R_Deg), Int_R_Weights(1:Int_R_Deg) )
ALLOCATE( Int_T_Locs(1:Int_T_Deg), Int_T_Weights(1:Int_T_Deg) )
ALLOCATE( Int_P_Locs(1:Int_P_Deg), Int_P_Weights(1:Int_P_Deg) )

ALLOCATE( LP_LP_Table( 1:Int_R_Deg, 0:Degree, 0:Degree, 0:1, 0:1) )
ALLOCATE( TP_TP_Integrals( 1:LM_Length, 1:LM_Length, 1:16) )

ALLOCATE( Slm_Table(1:Int_TP_Deg, 1:LM_Length)      )
ALLOCATE( Slm_dt_Table(1:Int_TP_Deg, 1:LM_Length)   )
ALLOCATE( Slm_dp_Table(1:Int_TP_Deg, 1:LM_Length)   )

ALLOCATE( Plm_Table(1:Int_T_Deg, 1:LM_Short_Length,1)      )
ALLOCATE( Plm_dt_Table(1:Int_T_Deg, 1:LM_Short_Length,1)   )

ALLOCATE( Am_Table(1:Int_P_Deg, -L_Limit:L_Limit,1)      )
ALLOCATE( Am_dp_Table(1:Int_P_Deg, -L_Limit:L_Limit,1)   )

ALLOCATE( Nlm_Table(1:LM_Short_Length) )

ALLOCATE( TP_Int_Weights(1:Int_TP_Deg)     )


END SUBROUTINE Allocate_Matrix_Init_Variables



 !+201+################################################################!
!                                                                       !
!          Deallocate_Matrix_Init_Variables                             !
!                                                                       !
 !#####################################################################!
SUBROUTINE Deallocate_Matrix_Init_Variables()

DEALLOCATE( Cur_T_Locs )
DEALLOCATE( Cur_P_Locs )

DEALLOCATE( Int_R_Locs, Int_R_Weights )
DEALLOCATE( Int_T_Locs, Int_T_Weights )
DEALLOCATE( Int_P_Locs, Int_P_Weights )

DEALLOCATE( LP_LP_Table )
DEALLOCATE( TP_TP_Integrals )

DEALLOCATE( Slm_Table      )
DEALLOCATE( Slm_dt_Table   )
DEALLOCATE( Slm_dp_Table   )
DEALLOCATE( Plm_Table      )
DEALLOCATE( Plm_dt_Table   )
DEALLOCATE( Am_Table       )
DEALLOCATE( Am_dp_Table    )
DEALLOCATE( Nlm_Table      )
DEALLOCATE( TP_Int_Weights )



END SUBROUTINE Deallocate_Matrix_Init_Variables



!+201+##########################################################################!
!                                                                               !
!           Calculate_Radial_Terms                                              !
!                                                                               !
!###############################################################################!
SUBROUTINE Calculate_Radial_Terms()

REAL(idp),  DIMENSION(0:Degree)                     ::  Lagrange_Poly_Values
REAL(idp),  DIMENSION(0:Degree)                     ::  Lagrange_DRV_Values
REAL(idp),  DIMENSION(0:Degree, 0:Int_R_Deg, 0:1)   ::  Lagrange_Poly_Table

INTEGER                                             ::  d, dp, dv, rd


CALL TimerStart(Timer_Matrix_Radial_Terms)


DO rd = 1,Int_R_Deg
    
    Lagrange_Poly_Values = Lagrange_Poly(Int_R_Locs(rd), Degree, FEM_Node_xlocs)
    Lagrange_DRV_Values  = Lagrange_Poly_Deriv(Int_R_Locs(rd), Degree, FEM_Node_xlocs)

    Lagrange_Poly_Table(:, rd, 0) = Lagrange_Poly_Values
    Lagrange_Poly_Table(:, rd, 1) = Lagrange_DRV_Values

END DO



DO dv = 0,1
DO dp = 0,Degree
DO d  = 0,Degree
DO rd = 1,Int_R_Deg

    LP_LP_Table(rd, d, dp, 0:1, dv) = Lagrange_Poly_Table(d, rd, 0:1)       &
                                    * Lagrange_Poly_Table(dp, rd, dv)
END DO
END DO
END DO
END DO

CALL TimerStop(Timer_Matrix_Radial_Terms)

END SUBROUTINE Calculate_Radial_Terms










!+202+##########################################################################!
!                                                                               !
!           Calculate_Angular_Terms                                             !
!                                                                               !
!###############################################################################!
SUBROUTINE Calculate_Angular_Terms()


INTEGER                                                         ::  l, m
INTEGER                                                         ::  td, pd, tpd
INTEGER                                                         ::  lm_loc, lpmp_loc
INTEGER                                                         ::  Long_LM, Short_LM

REAL(idp),                  DIMENSION(1:Int_TP_Deg)             ::  Cotan_Val
REAL(idp),                  DIMENSION(1:Int_Tp_Deg)             ::  Sin_Square




CALL TimerStart(Timer_Matrix_Angular_Terms)




DO td = 1,Int_T_Deg
DO pd = 1,Int_P_Deg

    tpd = (td-1)*Int_P_Deg + pd

    Sin_Square(tpd)     = DSIN( Cur_T_Locs(td) )*DSIN( Cur_T_Locs(td) )
    Cotan_Val(tpd)      = 1.0_idp/DTAN( Cur_T_Locs(td) )
    TP_Int_Weights(tpd) = DSIN( Cur_T_Locs(td) )                        &
                        * pi / 2.0_idp * Int_T_Weights(td)              &
                        * Int_P_Weights(pd)

END DO ! pd
END DO ! td

CALL Initialize_Nlm_Table(  L_Limit,                &
                            LM_Short_Length,        &
                            Nlm_Table               )


CALL Initialize_Am_Tables(  Int_P_Deg,              &
                            Int_P_Locs,             &
                            L_Limit,                &
                            1,                      &
                            [0, 0],                 &
                            [0.0_idp, TwoPi],       &
                            Am_Table,               &
                            Am_dp_Table             )


CALL Initialize_Plm_Tables( Int_T_Deg,              &
                            Int_T_Locs,             &
                            L_Limit,                &
                            LM_Short_Length,        &
                            1,                      &
                            [0, 0],                 &
                            [0.0_idp, Pi],          &
                            Plm_Table,              &
                            Plm_dt_Table            )


CALL Initialize_Slm_Tables_on_Elem( 0, 0,                   &
                                    Int_T_Deg,              &
                                    Int_P_Deg,              &
                                    [1, 1, 1],              &
                                    [0, 0, 0],              &
                                    Plm_Table,              &
                                    Plm_dt_Table,           &
                                    Am_Table,               &
                                    Am_dp_Table,            &
                                    Slm_Table,              &
                                    Slm_dt_Table,           &
                                    Slm_dp_Table            )



!PRINT*,"Slm_Table"
!PRINT*,Slm_Table
!PRINT*,"Slm_dt_Table"
!PRINT*,Slm_dt_Table
!PRINT*,"Slm_dp_Table"
!PRINT*,Slm_dp_Table



DO lpmp_loc = 1,LM_Length
DO lm_loc = 1,LM_Length

!
    ! S^lm * S^lpmp
    IF ( lm_loc == lpmp_loc ) THEN
        TP_TP_Integrals( lm_loc, lpmp_loc, 1)=1.0_idp
    ELSE
        TP_TP_Integrals( lm_loc, lpmp_loc, 1)=0.0_idp
    END IF


    ! SUM( TP_dTP_Factor(:,lm_loc,lpmp_loc)  )
    TP_TP_Integrals( lm_loc, lpmp_loc, 2 ) = SUM( Slm_Table( :, lm_loc )          &
                                                * Slm_dt_Table( :, lpmp_loc)   &
                                                * TP_Int_Weights(:)         )



    ! Slm * S^lpmp * Cotan
    TP_TP_Integrals( lm_loc, lpmp_loc, 3 ) = SUM( Slm_Table( :, lm_loc )          &
                                                * Slm_Table( :, lpmp_loc )     &
                                                * TP_Int_Weights(:)         &
                                                * Cotan_Val(:)              )

    ! d Slm/dt * S^lpmp     SUM( dTP_TP_Factor(:,lm_loc,lpmp_loc))
    TP_TP_Integrals( lm_loc, lpmp_loc, 4 ) = SUM( Slm_dt_Table(:, lm_loc )        &
                                                * Slm_Table( :, lpmp_loc )     &
                                                * TP_Int_Weights(:)         )

    ! Slm * d S^lpmp/dp
    TP_TP_Integrals( lm_loc, lpmp_loc, 5 ) = SUM( Slm_Table( :, lm_loc )          &
                                                * Slm_dp_Table( :, lpmp_loc )  &
                                                * TP_Int_Weights(:)         )

    ! d Slm/dp * S^lpmp
    TP_TP_Integrals( lm_loc, lpmp_loc, 6 ) = SUM( Slm_dp_Table( :, lpmp_loc )     &
                                                * Slm_Table( :, lpmp_loc )     &
                                                * TP_Int_Weights(:)         )

    ! SUM( dTP_dTP_Factor(:,lm_loc,lpmp_loc) )
    TP_TP_Integrals( lm_loc, lpmp_loc, 7 ) = SUM( Slm_dt_Table(:, lm_loc )        &
                                                * Slm_dt_Table( :, lpmp_loc )  &
                                                * TP_Int_Weights(:)         )


    ! SUM( TP_TP_Factor(:,lm_loc,lpmp_loc) / Sin_Square(:) )
    TP_TP_Integrals( lm_loc, lpmp_loc, 8 ) = SUM( Slm_Table( :, lm_loc )          &
                                                * Slm_Table( :, lpmp_loc )     &
                                                * TP_Int_Weights(:)         &
                                                / Sin_Square(:)             )

    ! SUM( TP_TP_Factor(:,lm_loc,lpmp_loc) * (1-Cotan_Val(:)**2) )
    TP_TP_Integrals( lm_loc, lpmp_loc, 9 ) = SUM( Slm_Table( :, lm_loc )          &
                                                * Slm_Table( :, lpmp_loc )     &
                                                * TP_Int_Weights(:)         &
                                                * (1-Cotan_Val(:)**2)       )

    ! SUM( dTP_TdP_Factor(:,lm_loc,lpmp_loc) )
    TP_TP_Integrals( lm_loc, lpmp_loc, 10 ) = SUM( Slm_dt_Table(:, lm_loc )       &
                                                 * Slm_dp_Table( :, lpmp_loc ) &
                                                 * TP_Int_Weights(:)        )

    ! SUM( TdP_TP_Factor(:,lm_loc,lpmp_loc) * Cotan_Val(:)   )
    TP_TP_Integrals( lm_loc, lpmp_loc, 11 ) = SUM( Slm_dp_Table( :, lpmp_loc )    &
                                                 * Slm_Table( :, lpmp_loc )    &
                                                 * TP_Int_Weights(:)        &
                                                 * Cotan_Val(:)             )

    ! SUM( TdP_TP_Factor(:,lm_loc,lpmp_loc) / Sin_Square(:)     )
    TP_TP_Integrals( lm_loc, lpmp_loc, 12 ) = SUM( Slm_dp_Table( :, lpmp_loc )    &
                                                 * Slm_Table( :, lpmp_loc )    &
                                                 * TP_Int_Weights(:)        &
                                                 / Sin_Square(:)            )

    ! SUM( TdP_dTP_Factor(:,lm_loc,lpmp_loc)/ Sin_Square(:) )
    TP_TP_Integrals( lm_loc, lpmp_loc, 13 ) = SUM( Slm_dp_Table( :, lpmp_loc )    &
                                                 * Slm_dt_Table( :, lpmp_loc ) &
                                                 * TP_Int_Weights(:)        &
                                                 / Sin_Square(:)            )

    ! SUM( TdP_TP_Factor(:,lm_loc,lpmp_loc) * Cotan_Val(:) / Sin_Square(:)     )
    TP_TP_Integrals( lm_loc, lpmp_loc, 14 ) = SUM( Slm_dp_Table( :, lpmp_loc )    &
                                                 * Slm_Table( :, lpmp_loc )    &
                                                 * TP_Int_Weights(:)        &
                                                 * Cotan_Val(:)             &
                                                 / Sin_Square(:)            )

    ! SUM( TdP_TdP_Factor(:,lm_loc,lpmp_loc)/Sin_Square(:) )
    TP_TP_Integrals( lm_loc, lpmp_loc, 15 ) = SUM( Slm_dp_Table( :, lpmp_loc )    &
                                                 * Slm_dp_Table( :, lpmp_loc ) &
                                                 * TP_Int_Weights(:)        &
                                                 / Sin_Square(:)            )

    ! SUM( dTP_TP_Factor(:,lm_loc,lpmp_loc) * Cotan_Val(:) )
    TP_TP_Integrals( lm_loc, lpmp_loc, 16 ) = SUM( Slm_dt_Table(:, lm_loc )       &
                                                 * Slm_Table( :, lpmp_loc )    &
                                                 * TP_Int_Weights(:)        &
                                                 * Cotan_Val(:)             )


!    PRINT*,TP_TP_Integrals(lm_loc, lpmp_loc,:)

END DO
END DO




!CALL Poseidon_Bailout("End of Calculate_Angular_Terms")

CALL TimerStop(Timer_Matrix_Angular_Terms)


END SUBROUTINE Calculate_Angular_Terms









!+301+##########################################################################!
!                                                                               !
!           Calculate_Laplace_Matrix                                            !
!                                                                               !
!###############################################################################!
SUBROUTINE Calculate_Laplace_Matrix()

INTEGER                                                 ::  l, re, d, dp
INTEGER                                                 ::  i, j, Here

REAL(KIND = idp)                                        ::  DR, TODR
REAL(KIND = idp)                                        ::  L_Lp1

REAL(idp),  DIMENSION(:),   ALLOCATABLE                 ::  CUR_R_LOCS
REAL(idp),  DIMENSION(:),   ALLOCATABLE                 ::  R_SQUARE


CALL TimerStart(Timer_Matrix_Laplace_Init)



ALLOCATE( Cur_R_locs(1:Int_R_Deg)  )
ALLOCATE( R_Square(1:Int_R_Deg)    )




IF ( Matrix_Format == 'Full') THEN


    IF ( Verbose_Flag ) CALL Init_Message('Initializing Laplace Matrix.  Format: Full.')

    DO l = 0,L_LIMIT
    DO re = 0,NUM_R_ELEMENTS-1

        L_Lp1 = REAL( l*(l+1), idp )

        DR   = rlocs(re+1) - rlocs(re)
        TODR = 2.0_idp/DR

        CUR_R_LOCS(:) = (DR/2.0_idp) * (Int_R_Locs(:)+1.0_idp) + rlocs(re)
        R_SQUARE = CUR_R_LOCS**2

        DO dp = 0,DEGREE
        DO d = 0,DEGREE

            i = Map_To_FEM_Node(re,d)
            j = Map_To_FEM_Node(re,dp)

            Laplace_Matrix_Full(i, j, l) = Laplace_Matrix_Full(i, j, l)                 &
                                         + SUM( R_SQUARE(:) * LP_LP_Table(:,d,dp,1,1)   &
                                                * TODR * Int_R_Weights(:)               &
                                               + L_Lp1 * LP_LP_Table(:,d,dp,0,0) / TODR        &
                                                * Int_R_Weights(:)                    )


        END DO  ! dp Loop
        END DO  ! d Loop
    END DO  ! re Loop
    END DO  ! l Loop







ELSEIF ( Matrix_Format == 'CCS') THEN


    IF ( Verbose_Flag ) CALL Init_Message('Initializing Laplace Matrix.  Format: CCS.')


    Laplace_Matrix_COL(0,:) = 0
    Laplace_Matrix_COL(1,:) = Laplace_Matrix_COL(0,:) + (DEGREE+1)
    HERE = 2

    DO re = 1,NUM_R_ELEMENTS-1
        DO d = 1,DEGREE - 1

            Laplace_Matrix_COL(Here,:) = Laplace_Matrix_COL(Here - 1,:) +  (DEGREE + 1)
            Here = Here + 1
        END DO

        Laplace_Matrix_COL(Here,:) = Laplace_Matrix_COL(Here - 1,:) + (2*DEGREE + 1)
        Here = Here + 1

    END DO


    DO d = 1, DEGREE

        Laplace_Matrix_COL(Here,:) = Laplace_Matrix_COL(Here - 1,:) + (DEGREE + 1)
        Here = Here + 1

    END DO



      !                             !
     !!                             !!
    !!!    ROW_IND INITIALIZATION   !!!
     !!                             !!
      !                             !
    Here = 0
    DO re = 0, NUM_R_ELEMENTS - 1
        DO d = 0,DEGREE - 1
        DO dp = 0, DEGREE

            Laplace_Matrix_ROW(Here,:) = re*DEGREE + dp
            Here = Here + 1

        END DO ! dp Loop
        END DO ! d Loop

        DO d = 0,DEGREE - 1

            Laplace_Matrix_ROW(Here,:) = re*DEGREE + d
            Here = Here + 1

        END DO ! d Loop
    END DO ! re Loop
    Laplace_Matrix_ROW(Here,:) = DEGREE * NUM_R_ELEMENTS


    Laplace_Matrix_VAL = 0.0_idp

    DO l = 0,L_LIMIT
        L_Lp1 = REAL( l*(l+1), idp )
        Here = 0
        DO re = 0,NUM_R_ELEMENTS-1

            DR   = rlocs(re+1) - rlocs(re)
            TODR = 2.0_idp/DR
            CUR_R_LOCS(:) = (DR/2.0_idp) * (Int_R_Locs(:)+1.0_idp) + rlocs(re)
            R_SQUARE = CUR_R_LOCS**2

            DO dp = 0,DEGREE
                DO d = 0,DEGREE


                    Laplace_Matrix_VAL(Here,l) = Laplace_Matrix_VAL(Here,l)     &
                                    + SUM( R_SQUARE(:) * LP_LP_Table(:,d,dp,1,1)    &
                                            * TODR * Int_R_Weights(:)               &
                                           + L_Lp1 * LP_LP_Table(:,d,dp,0,0) / TODR &
                                            * Int_R_Weights(:)                      )

                    Here = Here + 1

                END DO  ! dp Loop
            END DO  ! d Loop
            Here = Here - 1 !!! Take one step back due to overlap of first value
                            !!! in next element with last in current element
        END DO  ! re Loop
    END DO  ! l Loop


    Laplace_Factored_VAL = Laplace_Matrix_VAL
    Laplace_Factored_ROW = Laplace_Matrix_ROW
    Laplace_Factored_COL = Laplace_Matrix_COL


END IF


lPF_Init_Matrices_Flags(iPF_Init_Matrices_Type_A_Cholesky) = .FALSE.
lPF_Init_Matrices_Flags(iPF_Init_Matrices_Type_A) = .TRUE.



CALL TimerStop(Timer_Matrix_Laplace_Init)


END SUBROUTINE Calculate_Laplace_Matrix












!+101+###########################################################################!
!                                                                                !
!           Initialize_Laplace_Matrices_Full                                     !
!                                                                                !
!################################################################################!
SUBROUTINE Calculate_MVL_Banded()



INTEGER                                                 ::  l, m

INTEGER                                                 ::  re,rd
INTEGER                                                 ::  d, dp
INTEGER                                                 ::  i, j, ui

REAL(idp)                                               ::  L_Lp1

REAL(idp), ALLOCATABLE, DIMENSION(:)                    ::  Cur_R_Locs
REAL(idp), ALLOCATABLE, DIMENSION(:)                    ::  R_Square

REAL(idp),              DIMENSION(0:DEGREE)             ::  Reusable_Values


REAL(idp)                                               ::  DR, TODR



IF ( Verbose_Flag ) CALL Init_Message('Initializing Modified Vector Laplacian Matrix.  Format: Banded.')
CALL TimerStart(Timer_Matrix_MVL_Init)



dMB_Matrix_Banded = 0.0_idp



ALLOCATE( Cur_R_Locs(1:Int_R_Deg), R_Square(1:Int_R_Deg) )

ALLOCATE( RR_Factor(   1:Int_R_Deg, 0:DEGREE, 0:DEGREE )    )
ALLOCATE( RDR_Factor(  1:Int_R_Deg, 0:DEGREE, 0:DEGREE )    )
ALLOCATE( DRR_Factor(  1:Int_R_Deg, 0:DEGREE, 0:DEGREE )    )
ALLOCATE( DRDR_Factor( 1:Int_R_Deg, 0:DEGREE, 0:DEGREE )    )



DO re = 0,NUM_R_ELEMENTS-1
DO d = 0,DEGREE
DO ui = 1,3
DO l = 0,L_LIMIT


    L_Lp1 = REAL( l*(l+1), idp )

    DR            = rlocs(re+1) - rlocs(re)
    TODR          = 2.0_idp/DR
    CUR_R_LOCS(:) = (DR/2.0_idp) * (Int_R_Locs(:)+1.0_idp) + rlocs(re)
    R_SQUARE      = CUR_R_LOCS**2


    Reusable_Values = 0.0_idp
    DO dp = 0,DEGREE
        DO rd = 1,Int_R_Deg

            Reusable_Values(dp) = Reusable_Values(dp) &
                                - R_SQUARE(rd) * LP_LP_Table(rd,d,dp,1,1)   &
                                    * TODR * Int_R_Weights(rd)              &
                                - L_Lp1 * LP_LP_Table(rd,d,dp,0,0) / TODR   &
                                    * Int_R_Weights(rd)
                                    
        END DO  ! rd Loop
    END DO ! dp Loop


    DO m = -l,l
        j = FP_Beta_Array_Map(re,d,ui,l,m)

        DO dp = 0,Degree
            i = iMB_Bandwidth + FP_Beta_Array_Map(re,dp,ui,l,m)

            
            dMB_Matrix_Banded(i-j,j)                   &
                        = dMB_Matrix_Banded(i-j,j)     &
                          + Reusable_Values(dp)

        END DO ! dp Loop
    END DO  ! m Loop


END DO  ! l Loop
END DO  ! ui Loop
END DO  ! d Loop
END DO  ! re Loop




! At this point Laplace Matrix_Beta contains the discretized scalr Laplace operator for the three
! Beta components.  To test set up three Poisson/Laplace problems and assign one to each
! Beta component.  There is no coupling yet, so this should solve the three independently.



DO re = 0,Num_R_Elements-1

    DR = rlocs(re+1) - rlocs(re)
    TODR = 2.0_idp/DR

    CUR_R_LOCS(:) = (DR/2.0_idp) * (Int_R_Locs(:)+1.0_idp) + rlocs(re)
    R_SQUARE = CUR_R_LOCS**2

    CALL Calc_RR_Values( R_Square, TODR, RR_Factor, DRR_Factor, RDR_Factor, DRDR_Factor )


    DO d = 0,DEGREE
    DO l = 0,L_LIMIT
    DO m = -l,l

        CALL Calc_Beta1_Terms( re, d, l, m,                                         &
                                RR_Factor, dRR_Factor, dRdR_Factor,                 &
                                TP_TP_Integrals,                                    &
                                Cur_R_Locs, R_Square                                )


        CALL Calc_Beta2_Terms( re, d, l, m,                                         &
                                RR_Factor, dRR_Factor, dRdR_Factor, RdR_Factor,     &
                                TP_TP_Integrals,                                    &
                                Cur_R_Locs, R_Square                                )



        CALL Calc_Beta3_Terms( re, d, l, m,                                         &
                                RR_Factor, dRR_Factor, dRdR_Factor, RdR_Factor,     &
                                TP_TP_Integrals,                                    &
                                Cur_R_Locs, R_Square                                )



    END DO  ! mp Loop
    END DO  ! lp Loop
    END DO  ! dp Loop

END DO ! re Loop


DEALLOCATE( Cur_R_Locs, R_Square )
DEALLOCATE( RR_Factor    )
DEALLOCATE( RDR_Factor   )
DEALLOCATE( DRR_Factor   )
DEALLOCATE( DRDR_Factor  )



lPF_Init_Matrices_Flags(iPF_Init_Matrices_Type_B_LU) = .FALSE.
lPF_Init_Matrices_Flags(iPF_Init_Matrices_Type_B)    = .TRUE.

CALL TimerStop(Timer_Matrix_MVL_Init)


END SUBROUTINE Calculate_MVL_Banded







!+403+###########################################################################!
!                                                                                !
!                   CALC_RR_Values                                               !
!                                                                                !
!################################################################################!
SUBROUTINE CALC_RR_Values( R_Square, TODR,          &
                           RR_Factor, DRR_Factor,   &
                           RDR_Factor, DRDR_Factor  )

REAL(idp), DIMENSION(1:Int_R_Deg),                      INTENT(IN)      :: R_Square
REAL(idp),                                                      INTENT(IN)      :: TODR

REAL(idp), DIMENSION(1:Int_R_Deg, 0:DEGREE, 0:DEGREE ), INTENT(INOUT)   :: RR_Factor
REAL(idp), DIMENSION(1:Int_R_Deg, 0:DEGREE, 0:DEGREE ), INTENT(INOUT)   :: DRR_Factor
REAL(idp), DIMENSION(1:Int_R_Deg, 0:DEGREE, 0:DEGREE ), INTENT(INOUT)   :: RDR_Factor
REAL(idp), DIMENSION(1:Int_R_Deg, 0:DEGREE, 0:DEGREE ), INTENT(INOUT)   :: DRDR_Factor


INTEGER                                                                         :: d, dp
REAL(idp), ALLOCATABLE, DIMENSION(:)                                            :: R_Int_Weights



ALLOCATE( R_Int_Weights( 1:Int_R_Deg ) )




R_Int_Weights(:) = R_SQUARE(:) * INT_R_WEIGHTS(:)/TODR
DO d = 0,DEGREE
    DO dp = 0,DEGREE

        RR_Factor(:, d, dp)      = R_Int_Weights(:) * LP_LP_Table(:, d, dp, 0, 0)
        DRR_Factor(:, d, dp)     = R_Int_Weights(:) * LP_LP_Table(:, d, dp, 1, 0) * TODR
        RDR_Factor(:, d, dp)     = R_Int_Weights(:) * LP_LP_Table(:, d, dp, 0, 1) * TODR
        DRDR_Factor(:, d, dp)    = R_Int_Weights(:) * LP_LP_Table(:, d, dp, 1, 1) * TODR * TODR

    END DO
END DO

DEALLOCATE( R_Int_Weights)


END SUBROUTINE CALC_RR_Values






!+403+###########################################################################!
!                                                                                !
!                                                          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Beta1_Terms( re, dp, lp, mp,                                    &
                             RR_Factor, dRR_Factor, dRdR_Factor,                &
                             TP_TP_Integrals,                                   &
                             Cur_R_Locs, R_Square                               )

INTEGER,                                                        INTENT(IN)  :: re, dp, lp, mp

REAL(idp), DIMENSION(1:Int_R_Deg, 0:DEGREE, 0:DEGREE ), INTENT(IN)  :: RR_Factor
REAL(idp), DIMENSION(1:Int_R_Deg, 0:DEGREE, 0:DEGREE ), INTENT(IN)  :: DRR_Factor
REAL(idp), DIMENSION(1:Int_R_Deg, 0:DEGREE, 0:DEGREE ), INTENT(IN)  :: DRDR_Factor

REAL(idp), DIMENSION( 1:LM_Length, 1:LM_Length, 1:16),  INTENT(IN)  :: TP_TP_Integrals


REAL(idp),    DIMENSION( 1:Int_R_Deg ),                 INTENT(IN)  :: Cur_R_Locs
REAL(idp),    DIMENSION( 1:Int_R_Deg ),                 INTENT(IN)  :: R_Square


INTEGER                                                     :: d, ui, uj
INTEGER                                                     :: row, col
INTEGER                                                     :: lm_loc, lpmp_loc



ui       = 1
lpmp_loc = Map_To_lm(lp,mp)
Row      = iMB_Bandwidth + FP_Beta_Array_Map(re,dp,ui,lpmp_loc)




DO d = 0,Degree

    uj = 1

    DO lm_loc = 1,LM_Length
        Col = FP_Beta_Array_Map(re,d,uj,lm_loc)


        dMB_Matrix_Banded(Row-Col, Col) = dMB_Matrix_Banded(Row-Col, Col)               &
                                      - SUM( dRdR_Factor(:, d, dp) )/3.0_idp        &
                                        * TP_TP_Integrals( lm_loc, lpmp_loc, 1)     &
                                      - 8.0_idp/3.0_idp * SUM( RR_Factor(:, d, dp)  &
                                                               / R_Square(:)    )   &
                                        * TP_TP_Integrals( lm_loc, lpmp_loc, 1)

        
    END DO ! lpmp_loc Loop





    uj = 2

    DO lm_loc = 1,LM_Length
        Col = FP_Beta_Array_Map(re,d,uj,lm_loc)

        dMB_Matrix_Banded(Row-Col, Col) = dMB_Matrix_Banded(Row-Col, Col)               &
                                      - SUM( dRR_Factor(:, d, dp) )/3.0_idp         &
                                        * TP_TP_Integrals( lm_loc, lpmp_loc, 2 )    &
                                      - SUM( dRR_Factor(:, d, dp) )/3.0_idp         &
                                        * TP_TP_Integrals( lm_loc, lpmp_loc, 3 )    &
                                      - 2.0_idp * SUM( RR_Factor(:, d, dp)          &
                                                        / CUR_R_LOCS(:)         )   &
                                        * TP_TP_Integrals( lm_loc, lpmp_loc, 4 )    &
                                      - 2.0_idp * SUM( RR_Factor(:, d, dp)          &
                                                        / CUR_R_LOCS(:)         )   &
                                        * TP_TP_Integrals( lm_loc, lpmp_loc, 3 )



    END DO ! lpmp_loc Loop
    



    uj = 3 ! beta^phi

    DO lm_loc = 1,LM_Length
        Col = FP_Beta_Array_Map(re,d,uj,lm_loc)


        dMB_Matrix_Banded(Row-Col, Col) = dMB_Matrix_Banded(Row-Col, Col)               &
                                      - SUM( dRR_Factor(:, d, dp)  )/3.0_idp        &
                                        * TP_TP_Integrals( lm_loc, lpmp_loc, 5 )    &
                                      - 2.0_idp * SUM( RR_Factor(:, d, dp)         &
                                                        / CUR_R_LOCS(:)        )   &
                                        * TP_TP_Integrals( lm_loc, lpmp_loc, 6 )

    END DO ! lpmp_loc Loop

END DO ! dp Loop




END SUBROUTINE Calc_Beta1_Terms





!+403+###########################################################################!
!                                                                                !
!                   Calc_Beta2_Terms                                       !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Beta2_Terms( re, dp, lp, mp,                                    &
                             RR_Factor, dRR_Factor, dRdR_Factor, RdR_Factor,    &
                             TP_TP_Integrals,                                   &
                             Cur_R_Locs, R_Square                               )

INTEGER,                                                        INTENT(IN)  :: re, dp, lp, mp

REAL(idp), DIMENSION(1:Int_R_Deg, 0:DEGREE, 0:DEGREE ), INTENT(IN)  :: RR_Factor
REAL(idp), DIMENSION(1:Int_R_Deg, 0:DEGREE, 0:DEGREE ), INTENT(IN)  :: DRR_Factor
REAL(idp), DIMENSION(1:Int_R_Deg, 0:DEGREE, 0:DEGREE ), INTENT(IN)  :: RdR_Factor
REAL(idp), DIMENSION(1:Int_R_Deg, 0:DEGREE, 0:DEGREE ), INTENT(IN)  :: DRDR_Factor

REAL(idp), DIMENSION( 1:LM_Length, 1:LM_Length, 1:16),  INTENT(IN)  :: TP_TP_Integrals


REAL(idp),    DIMENSION( 1:Int_R_Deg ),                 INTENT(IN)  :: Cur_R_Locs
REAL(idp),    DIMENSION( 1:Int_R_Deg ),                 INTENT(IN)  :: R_Square


INTEGER                                                     :: d, ui, uj
INTEGER                                                     :: row, col
INTEGER                                                     :: lm_loc, lpmp_loc



ui       = 2
lpmp_loc = Map_To_lm(lp,mp)
Row      = iMB_Bandwidth + FP_Beta_Array_Map(re,dp,ui,lpmp_loc)


DO d = 0,Degree

uj = 1
DO lm_loc = 1,LM_Length
    Col = FP_Beta_Array_Map(re,d,uj,lm_loc)
    

    dMB_Matrix_Banded(Row-Col, Col) = dMB_Matrix_Banded(Row-Col, Col)               &
                                  + SUM( RdR_Factor(:, d, dp)                   &
                                         /(3.0_idp*R_Square(:) )       )        &   ! Term 1
                                    * TP_TP_Integrals( lm_loc, lpmp_loc, 4 )    &
                                  - SUM( 8.0_idp * RR_Factor(:, d, dp)          &   ! Term 2
                                         /(3.0_idp*R_Square(:)*Cur_R_Locs(:)))  &
                                    * TP_TP_Integrals( lm_loc, lpmp_loc, 4 )

END DO ! lpmp_loc Loop




uj = 2
DO lm_loc = 1,LM_Length
    Col = FP_Beta_Array_Map(re,d,uj,lm_loc)

    dMB_Matrix_Banded(Row-Col, Col) = dMB_Matrix_Banded(Row-Col, Col)               &
                                  - SUM( RR_Factor(:, d, dp)                    &
                                         /(3.0_idp*R_Square(:) )       )        &
                                    * TP_TP_Integrals( lm_loc, lpmp_loc, 7 )    &
                                  - SUM( RR_Factor(:, d, dp)                    &
                                         /(3.0_idp*R_Square(:) )           )    &
                                    * TP_TP_Integrals( lm_loc, lpmp_loc, 4 )    &
                                  + SUM( 2.0_idp * dRR_Factor(:, d, dp)         &
                                          /Cur_R_Locs(:)                    )   &
                                    * TP_TP_Integrals( lm_loc, lpmp_loc, 1 )    &
                                  - SUM( RR_Factor(:, d, dp)                    &
                                         / (3.0_idp*R_Square(:) )           )   &
                                    * TP_TP_Integrals( lm_loc, lpmp_loc, 8 )    &
                                  + SUM( RR_Factor(:, d, dp)                    &
                                         / R_Square(:)                     )   &
                                    * TP_TP_Integrals( lm_loc, lpmp_loc, 9 )
    

END DO ! lpmp_loc Loop





uj = 3 ! beta^phi

DO lm_loc = 1,LM_Length
    Col = FP_Beta_Array_Map(re,d,uj,lm_loc)

    dMB_Matrix_Banded(Row-Col, Col) = dMB_Matrix_Banded(Row-Col, Col)           &
                                  - SUM( RR_Factor(:, d, dp)                    &
                                         / (3.0_idp*R_Square(:) )           )   &
                                    * TP_TP_Integrals( lm_loc, lpmp_loc, 10 )   &
                                  - SUM ( 2.0_idp * RR_Factor(:, d, dp)         &
                                          / R_Square(:)                     )   &
                                    * TP_TP_Integrals( lm_loc, lpmp_loc, 11 )

END DO ! lpmp_loc Loop

END DO ! dp Loop



END SUBROUTINE Calc_Beta2_Terms













!+403+###########################################################################!
!                                                                                !
!                   Calc_Beta3_Terms                                       !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Beta3_Terms( re, dp, lp, mp,                                    &
                             RR_Factor, dRR_Factor, dRdR_Factor, RdR_Factor,    &
                             TP_TP_Integrals,                                   &
                             Cur_R_Locs, R_Square                               )

INTEGER,                                                        INTENT(IN)    :: re, dp, lp, mp

REAL(idp), DIMENSION(1:Int_R_Deg, 0:DEGREE, 0:DEGREE ), INTENT(IN)    :: RR_Factor
REAL(idp), DIMENSION(1:Int_R_Deg, 0:DEGREE, 0:DEGREE ), INTENT(IN)    :: DRR_Factor
REAL(idp), DIMENSION(1:Int_R_Deg, 0:DEGREE, 0:DEGREE ), INTENT(IN)    :: RDR_Factor
REAL(idp), DIMENSION(1:Int_R_Deg, 0:DEGREE, 0:DEGREE ), INTENT(IN)    :: DRDR_Factor

REAL(idp), DIMENSION( 1:LM_Length, 1:LM_Length, 1:16),  INTENT(IN)    :: TP_TP_Integrals

REAL(idp),    DIMENSION( 1:Int_R_Deg ),                 INTENT(IN)    :: Cur_R_Locs
REAL(idp),    DIMENSION( 1:Int_R_Deg ),                 INTENT(IN)    :: R_Square


INTEGER                                                     :: d, ui, uj
INTEGER                                                     :: row, col
INTEGER                                                     :: lm_loc, lpmp_loc


lpmp_loc = Map_To_lm(lp,mp)

ui       = 3
lpmp_loc = Map_To_lm(lp,mp)
Row      = iMB_Bandwidth + FP_Beta_Array_Map(re,dp,ui,lpmp_loc)



DO d = 0,Degree

uj = 1

DO lm_loc = 1,LM_Length

    Col = FP_Beta_Array_Map(re,d,uj,lm_loc)

    dMB_Matrix_Banded(Row-Col, Col) = dMB_Matrix_Banded(Row-Col, Col)                   &
                                  - SUM( RdR_Factor(:, d, dp)                       &   ! Term 1
                                        /(3.0_idp * R_Square(:) )    )              &
                                    * TP_TP_Integrals( lm_loc, lpmp_loc, 6 )        &
                                  + SUM( 8.0_idp * RR_Factor(:, d, dp)              &   ! Term 2
                                        /(3.0_idp*Cur_R_Locs(:) * R_Square(:) ) )   &
                                    * TP_TP_Integrals( lm_loc, lpmp_loc, 12 )

END DO ! lm_loc Loop



uj = 2

DO lm_loc = 1,LM_Length
    Col = FP_Beta_Array_Map(re,d,uj,lm_loc)
    



    dMB_Matrix_Banded(Row-Col, Col) = dMB_Matrix_Banded(Row-Col, Col)               &
                                  - SUM( RR_Factor(:, d, dp)                    &   ! Term 1
                                         /( 3.0_idp * R_Square(:) )   )         &
                                    * TP_TP_Integrals( lm_loc, lpmp_loc, 13 )   &
                                  + SUM( 8.0_idp * RR_Factor(:, d, dp)          &   ! Term 2
                                         / ( 3.0_idp * R_Square(:) ) )          &
                                    * TP_TP_Integrals( lm_loc, lpmp_loc, 14 )

END DO ! lpmp_loc Loop



uj = 3 ! beta^phi
DO lm_loc = 1,LM_Length
    Col = FP_Beta_Array_Map(re,d,uj,lm_loc)

    dMB_Matrix_Banded(Row-Col, Col) = dMB_Matrix_Banded(Row-Col, Col)           &
                                  - SUM( RR_Factor(:, d, dp )                   &   ! Term 1
                                         /(3.0_idp * R_Square(:) )    )         &
                                    * TP_TP_Integrals( lm_loc, lpmp_loc, 15 )   &
                                  + SUM( 2.0_idp * dRR_Factor(:, d, dp)         &   ! Term 2
                                         / Cur_R_Locs(:)                    )   &
                                    * TP_TP_Integrals( lm_loc, lpmp_loc, 1)     &
                                  + SUM ( 2.0_idp * RR_Factor(:, d, dp)         &   ! Term 3
                                            / R_Square(:)                   )   &
                                    * TP_TP_Integrals( lm_loc, lpmp_loc, 16 )

END DO ! lpmp_loc Loop


END DO ! dp Loop








END SUBROUTINE Calc_Beta3_Terms




END MODULE Matrix_Initialization_Module

