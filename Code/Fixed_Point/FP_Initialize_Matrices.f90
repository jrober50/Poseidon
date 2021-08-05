   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE FP_Intialize_Matrices                                                        !##!
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

USE Poseidon_Numbers_Module, &
            ONLY : pi

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
                    Var_Dim


USE Variables_FP, &
            ONLY :  Matrix_Format,              &
                    Laplace_Matrix_Full,        &
                    Laplace_Matrix_VAL,         &
                    Laplace_Matrix_ROW,         &
                    Laplace_Matrix_COL,         &
                    Laplace_Factored_VAL,       &
                    Laplace_Factored_ROW,       &
                    Laplace_Factored_COL,       &
                    CFA_EQ_Flags,               &
                    CFA_EQ_Map,                 &
                    CFA_Mat_Map,                &
                    MCF_Flag,                   &
                    Beta_MVL_Banded,            &
                    Beta_Factorized_Flag,       &
                    Beta_Bandwidth


USE Poseidon_Cholesky_Module,   &
            ONLY :  Cholesky_Factorization

USE FP_Functions_Laplace_Beta, &
            ONLY :  Initialize_Laplace_Matrices_Beta

USE FP_Beta_Banded, &
            ONLY :  Initialize_Beta_MVL_Banded

USE IO_FP_Linear_System, &
            ONLY :  Output_Laplace

USE FP_Functions_Mapping, &
            ONLY :  FP_Beta_Array_Map,          &
                    FP_FEM_Node_Map,            &
                    FP_LM_Map

USE Functions_Math, &
            ONLY :  Lagrange_Poly,          &
                    Lagrange_Poly_Deriv,    &
                    Legendre_Poly,          &
                    Norm_Factor

USE Functions_Mapping, &
            ONLY :  Map_From_X_Space

USE Functions_Quadrature, &
            ONLY :  Initialize_LGL_Quadrature_Locations,    &
                    Initialize_LG_Quadrature,               &
                    Initialize_Trapezoid_Quadrature


USE MPI

IMPLICIT NONE


INTEGER                                                 ::  Int_R_Deg
INTEGER                                                 ::  Int_T_Deg
INTEGER                                                 ::  Int_P_Deg
INTEGER                                                 ::  Int_TP_Deg


REAL(idp),      ALLOCATABLE,    DIMENSION(:)            ::  Cur_T_Locs
REAL(idp),      ALLOCATABLE,    DIMENSION(:)            ::  Cur_P_Locs

REAL(idp),      ALLOCATABLE,    DIMENSION(:)            ::  Int_R_Locs, Int_R_Weights
REAL(idp),      ALLOCATABLE,    DIMENSION(:)            ::  Int_T_Locs, Int_T_Weights
REAL(idp),      ALLOCATABLE,    DIMENSION(:)            ::  Int_P_Locs, Int_P_Weights
REAL(idp),      ALLOCATABLE,    DIMENSION(:)            ::  Int_TP_Weights

REAL(idp),      ALLOCATABLE,    DIMENSION( :,:,:,:,: )  ::  LP_LP_Table
COMPLEX(idp),   ALLOCATABLE,    DIMENSION( :,:,: )      ::  TP_TP_Integrals

REAL(idp),      ALLOCATABLE,    DIMENSION( :,:,: )      ::  RR_Factor
REAL(idp),      ALLOCATABLE,    DIMENSION( :,:,: )      ::  DRR_Factor
REAL(idp),      ALLOCATABLE,    DIMENSION( :,:,: )      ::  RDR_Factor
REAL(idp),      ALLOCATABLE,    DIMENSION( :,:,: )      ::  DRDR_Factor


COMPLEX(idp),   ALLOCATABLE,    DIMENSION( :,: )        ::  Ylm
COMPLEX(idp),   ALLOCATABLE,    DIMENSION( :,: )        ::  Ylm_dt
COMPLEX(idp),   ALLOCATABLE,    DIMENSION( :,: )        ::  Ylm_dp

COMPLEX(idp),   ALLOCATABLE,    DIMENSION( :,: )        ::  Ylm_CC
COMPLEX(idp),   ALLOCATABLE,    DIMENSION( :,: )        ::  Ylm_CC_dt
COMPLEX(idp),   ALLOCATABLE,    DIMENSION( :,: )        ::  Ylm_CC_dp

CONTAINS


!+101+##########################################################################!
!                                                                               !
!           Initialize_FP_Matrices                                              !
!                                                                               !
!###############################################################################!
SUBROUTINE Initialize_FP_Matrices()



! Set number of quadrature points for different dimensions
!
! Radial : Using Gaussian Quadrature. 2n-1 points for exact integration of Polynomials of degree n.
!           Radially this is a product of two Lagrange Polys so n = 2*Degree.
! Phi    : Using Trapezoid rule quadrature. m + 1 for exact integration of exp(-i*m*phi).
!           Max m is defined by the product of exp(i*L_Limit*Phi)*exp(i*L_Limit*Phi) = exp(i*(2*L_Limit)*Phi)
!
Int_R_Deg  = 4*Degree - 1
Int_T_Deg  = 15
Int_P_Deg  = 2*L_LIMIT + 1

Int_TP_Deg = Int_T_Deg*Int_P_Deg


ALLOCATE( Cur_T_Locs(1:Int_T_Deg) )
ALLOCATE( Cur_P_Locs(1:Int_P_Deg) )

ALLOCATE( Int_R_Locs(1:Int_R_Deg), Int_R_Weights(1:Int_R_Deg) )
ALLOCATE( Int_T_Locs(1:Int_T_Deg), Int_T_Weights(1:Int_T_Deg) )
ALLOCATE( Int_P_Locs(1:Int_P_Deg), Int_P_Weights(1:Int_P_Deg) )
ALLOCATE( Int_TP_Weights(1:Int_TP_Deg) )

ALLOCATE( LP_LP_Table( 1:Int_R_Deg, 0:Degree, 0:Degree, 0:1, 0:1) )
ALLOCATE( TP_TP_Integrals( 1:LM_Length, 1:LM_Length, 1:16) )


ALLOCATE( Ylm(1:Int_TP_Deg,1:LM_Length) )
ALLOCATE( Ylm_dt(1:Int_TP_Deg,1:LM_Length) )
ALLOCATE( Ylm_dp(1:Int_TP_Deg,1:LM_Length) )

ALLOCATE( Ylm_CC(1:Int_TP_Deg,1:LM_Length) )
ALLOCATE( Ylm_CC_dt(1:Int_TP_Deg,1:LM_Length) )
ALLOCATE( Ylm_CC_dp(1:Int_TP_Deg,1:LM_Length) )



CALL Initialize_LG_Quadrature(Int_R_Deg, Int_R_Locs, Int_R_Weights)
CALL Initialize_LG_Quadrature(Int_T_Deg, Int_T_Locs, Int_T_Weights)
CALL Initialize_Trapezoid_Quadrature(Int_P_Deg, Int_P_Locs, Int_P_Weights)



 !                                                          !
!!   Map Theta and Phi Locations from [-1,1] space to real space.   !!
 !                                                          !
Cur_T_Locs = Map_From_X_Space( 0.0_idp, pi, Int_T_Locs )
Cur_P_Locs = Int_P_Locs





CALL Calculate_Radial_Terms()
CALL Calculate_Angular_Terms()


CALL Calculate_Laplace_Matrix()
CALL Calculate_MVL_Banded()

DEALLOCATE( Cur_T_Locs )
DEALLOCATE( Cur_P_Locs )

DEALLOCATE( Int_R_Locs, Int_R_Weights )
DEALLOCATE( Int_T_Locs, Int_T_Weights )
DEALLOCATE( Int_P_Locs, Int_P_Weights )
DEALLOCATE( Int_TP_Weights )

DEALLOCATE( LP_LP_Table )
DEALLOCATE( TP_TP_Integrals )


DEALLOCATE( Ylm    )
DEALLOCATE( Ylm_dt )
DEALLOCATE( Ylm_dp )

DEALLOCATE( Ylm_CC )
DEALLOCATE( Ylm_CC_dt )
DEALLOCATE( Ylm_CC_dp )



END SUBROUTINE Initialize_FP_Matrices








!+201+##########################################################################!
!                                                                               !
!           Calculate_Radial_Terms                                              !
!                                                                               !
!###############################################################################!
SUBROUTINE Calculate_Radial_Terms()

REAL(idp),  DIMENSION(0:Degree)                     ::  Local_Locations

REAL(idp),  DIMENSION(0:Degree)                     ::  Lagrange_Poly_Values
REAL(idp),  DIMENSION(0:Degree)                     ::  Lagrange_DRV_Values
REAL(idp),  DIMENSION(0:Degree, 0:Int_R_Deg, 0:1)   ::  Lagrange_Poly_Table

INTEGER                                             ::  d, dp, dv, rd

Local_Locations = Initialize_LGL_Quadrature_Locations(Degree)

DO rd = 1,Int_R_Deg
    
    Lagrange_Poly_Values = Lagrange_Poly(Int_R_Locs(rd), Degree, Local_Locations)
    Lagrange_DRV_Values  = Lagrange_Poly_Deriv(Int_R_Locs(rd), Degree, Local_Locations)

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



END SUBROUTINE Calculate_Radial_Terms










!+202+##########################################################################!
!                                                                               !
!           Calculate_Angular_Terms                                             !
!                                                                               !
!###############################################################################!
SUBROUTINE Calculate_Angular_Terms()


INTEGER                                                         ::  l, m
INTEGER                                                         ::  td, pd, tpd


REAL(idp), DIMENSION(1:1)                                       ::  Legendre_Poly_Value


REAL(idp)                                                       ::  Norm_Storage
REAL(idp)                                                       ::  REAL_L


REAL(idp), DIMENSION(-L_LIMIT:L_LIMIT)                          ::  M_POWER_TABLE

INTEGER                                                         ::  lm_loc, lpmp_loc

COMPLEX(idp), ALLOCATABLE, DIMENSION(:,:,:)                     ::  Ylm_Table

COMPLEX(idp),    DIMENSION(1:LM_LENGTH)                         ::  Sqrt_Term
REAL(idp),       DIMENSION(1:Int_TP_Deg)                        ::  CSC_Val
REAL(idp),       DIMENSION(1:Int_TP_Deg)                        ::  Cotan_Val
REAL(idp),       DIMENSION(1:Int_Tp_Deg)                        ::  Sin_Square


COMPLEX(idp)                                                    ::  Tmp_Value_A
COMPLEX(idp)                                                    ::  Tmp_Value_B

REAL(idp), ALLOCATABLE, DIMENSION(:)                            ::  TP_Int_Weights

INTEGER,        ALLOCATABLE, DIMENSION(:)                       ::  M_Values


ALLOCATE( Ylm_Table(1:Int_TP_Deg, -L_LIMIT:L_LIMIT, -1:L_LIMIT )    )
ALLOCATE( TP_Int_Weights( 1:Int_TP_Deg)     )
ALLOCATE( M_VALUES(0:L_LIMIT) )

M_POWER_TABLE(0) = 1.0_idp
IF ( L_LIMIT > 0 ) THEN
    DO m = 1, L_LIMIT

        M_POWER_TABLE(m) = -1.0_idp*M_POWER_TABLE(m-1)
        M_POWER_TABLE(-m) = -1.0_idp*M_POWER_TABLE(m-1)
    END DO
END IF


IF ( DOMAIN_DIM == 1 ) THEN
     M_VALUES = 0
ELSE IF ( DOMAIN_DIM == 2 ) THEN
     M_VALUES = 0
ELSE IF ( DOMAIN_DIM == 3 ) THEN
     M_VALUES = (/(l,l=0,L_LIMIT,1)/)
END IF


Sqrt_Term(1) = 0.0_idp
DO l = 1,L_LIMIT

    REAL_L = REAL(l, idp)

    DO m = -M_VALUES(l),M_VALUES(l)

        Tmp_Value_A = COMPLEX( (2.0_idp * REAL_L + 1.0_idp)/(2.0_idp*Real_L - 1.0_idp),0.0_idp )
        Tmp_Value_B = COMPLEX( (l-m)*(l+m),0.0_idp )

        Sqrt_Term(FP_LM_Map(l,m)) = zsqrt( Tmp_Value_A)*zsqrt(Tmp_Value_B)

    END DO ! m Loop
END DO ! l Loop







DO td = 1,Int_T_Deg
DO pd = 1,Int_P_Deg

    tpd = (td-1)*Int_P_Deg + pd

    Sin_Square(tpd)     = DSIN( Cur_T_Locs(td) )*DSIN( Cur_T_Locs(td) )
    CSC_VAL(tpd)        = 1.0_idp/DSIN(Cur_T_Locs(td))
    Cotan_Val(tpd)      = 1.0_idp/DTAN( Cur_T_Locs(td) )
    TP_Int_Weights(tpd) = DSIN( Cur_T_Locs(td) )                        &
                        * pi / 2.0_idp * Int_T_Weights(td)              &
                        * Int_P_Weights(pd)

END DO ! pd
END DO ! td


        
Ylm_Table = 0.0_idp
DO l = 0,L_LIMIT
DO m = -M_VALUES(l),M_VALUES(l)
DO td = 1,Int_T_Deg
DO pd = 1,Int_P_Deg

    tpd = (td-1)*Int_P_Deg + pd
    Legendre_Poly_Value = Legendre_Poly(l, m, 1, Cur_T_Locs(td))

    Ylm_Table(tpd, m, l ) = Norm_Factor(l,m)                            &   ! Normalization Factor
                          * Legendre_Poly_Value(1)                      &   ! Legendre Polynomial
                          * CDEXP(CMPLX(0.0_idp,m*Cur_P_Locs(pd),idp))      ! exp(im phi)

END DO ! pd Loop
END DO ! td Loop
END DO ! m Loop
END DO ! l Loop





DO l = 0,L_LIMIT
DO m = -M_VALUES(l),M_VALUES(l)


    REAL_L = REAL(l, idp)

    lm_loc  = FP_LM_Map(l,m)
    Norm_Storage = Norm_Factor(l,m)


    Ylm(:, lm_loc )    = Ylm_Table(:,m,l)
    Ylm_dt(:, lm_loc ) = REAL_L * Cotan_Val(:) * Ylm_Table(:,m,l)                     &
                        - SQRT_TERM(lm_loc) * CSC_VAL(:) * Ylm_Table(:,m,l-1)
    Ylm_dp(:, lm_loc ) = CMPLX(0,m,idp) * Ylm_Table(:,m,l)


    Ylm_CC( :, lm_loc )    = M_POWER_TABLE(m) * Ylm_Table(:,-m,l)
    Ylm_CC_dt( :, lm_loc ) = REAL_L*Cotan_Val(:) * Ylm_CC( :, lm_loc)          &
                             - SQRT_TERM(lm_loc) * CSC_VAL(:)                       &
                                * M_POWER_TABLE(m) * Ylm_Table(:,-m, l-1)
    Ylm_CC_dp( :, lm_loc ) = CMPLX(0,-m,idp) * Ylm_CC(:, lm_loc)


END DO ! m Loop
END DO ! l Loop







DO lpmp_loc = 1,LM_Length
DO lm_loc = 1,LM_Length

!
    ! Ylm * Y^lpmp
    IF ( lm_loc == lpmp_loc ) THEN
        TP_TP_Integrals( lm_loc, lpmp_loc, 1)=1.0_idp
    ELSE
        TP_TP_Integrals( lm_loc, lpmp_loc, 1)=0.0_idp
    END IF


    ! SUM( TP_dTP_Factor(:,lm_loc,lpmp_loc)  )
    TP_TP_Integrals( lm_loc, lpmp_loc, 2 ) = SUM( Ylm( :, lm_loc )          &
                                                * Ylm_CC_dp( :, lpmp_loc)   &
                                                * TP_Int_Weights(:)         )



    ! Ylm * Y^lpmp * Cotan
    TP_TP_Integrals( lm_loc, lpmp_loc, 2 ) = SUM( Ylm( :, lm_loc )          &
                                                * Ylm_CC( :, lpmp_loc )     &
                                                * TP_Int_Weights(:)         &
                                                * Cotan_Val(:)              )

    ! d Ylm/dt * Y^lpmp     SUM( dTP_TP_Factor(:,lm_loc,lpmp_loc))
    TP_TP_Integrals( lm_loc, lpmp_loc, 3 ) = SUM( Ylm_dt(:, lm_loc )        &
                                                * Ylm_CC( :, lpmp_loc )     &
                                                * TP_Int_Weights(:)         )

    ! Ylm * d Y^lpmp/dp
    TP_TP_Integrals( lm_loc, lpmp_loc, 4 ) = SUM( Ylm( :, lm_loc )          &
                                                * Ylm_CC_dt( :, lpmp_loc )  &
                                                * TP_Int_Weights(:)         )

    ! d Ylm/dp * Y^lpmp
    TP_TP_Integrals( lm_loc, lpmp_loc, 5 ) = SUM( Ylm_dp( :, lpmp_loc )     &
                                                * Ylm_CC( :, lpmp_loc )     &
                                                * TP_Int_Weights(:)         )

    ! SUM( dTP_dTP_Factor(:,lm_loc,lpmp_loc) )
    TP_TP_Integrals( lm_loc, lpmp_loc, 6 ) = SUM( Ylm_dt(:, lm_loc )        &
                                                * Ylm_CC_dt( :, lpmp_loc )  &
                                                * TP_Int_Weights(:)         )

    ! SUM( dTP_TP_Factor(:,lm_loc,lpmp_loc) * Cotan_Val(:) )
    TP_TP_Integrals( lm_loc, lpmp_loc, 7 ) = SUM( Ylm_dt(:, lm_loc )        &
                                                * Ylm_CC( :, lpmp_loc )     &
                                                * TP_Int_Weights(:)         &
                                                * Cotan_Val(:)              )

    ! SUM( TP_TP_Factor(:,lm_loc,lpmp_loc) / Sin_Square(:) )
    TP_TP_Integrals( lm_loc, lpmp_loc, 8 ) = SUM( Ylm( :, lm_loc )          &
                                                * Ylm_CC( :, lpmp_loc )     &
                                                * TP_Int_Weights(:)         &
                                                / Sin_Square(:)             )

    ! SUM( TP_TP_Factor(:,lm_loc,lpmp_loc) * (1-Cotan_Val(:)**2) )
    TP_TP_Integrals( lm_loc, lpmp_loc, 9 ) = SUM( Ylm( :, lm_loc )          &
                                                * Ylm_CC( :, lpmp_loc )     &
                                                * TP_Int_Weights(:)         &
                                                * (1-Cotan_Val(:)**2)       )

    ! SUM( dTP_TdP_Factor(:,lm_loc,lpmp_loc) )
    TP_TP_Integrals( lm_loc, lpmp_loc, 10 ) = SUM( Ylm_dt(:, lm_loc )       &
                                                 * Ylm_CC_dp( :, lpmp_loc ) &
                                                 * TP_Int_Weights(:)        )

    ! SUM( TdP_TP_Factor(:,lm_loc,lpmp_loc) * Cotan_Val(:)   )
    TP_TP_Integrals( lm_loc, lpmp_loc, 11 ) = SUM( Ylm_dp( :, lpmp_loc )    &
                                                 * Ylm_CC( :, lpmp_loc )    &
                                                 * TP_Int_Weights(:)        &
                                                 * Cotan_Val(:)             )

    ! SUM( TdP_TP_Factor(:,lm_loc,lpmp_loc) / Sin_Square(:)     )
    TP_TP_Integrals( lm_loc, lpmp_loc, 12 ) = SUM( Ylm_dp( :, lpmp_loc )    &
                                                 * Ylm_CC( :, lpmp_loc )    &
                                                 * TP_Int_Weights(:)        &
                                                 / Sin_Square(:)            )

    ! SUM( TdP_dTP_Factor(:,lm_loc,lpmp_loc)/ Sin_Square(:) )
    TP_TP_Integrals( lm_loc, lpmp_loc, 13 ) = SUM( Ylm_dp( :, lpmp_loc )    &
                                                 * Ylm_CC_dt( :, lpmp_loc ) &
                                                 * TP_Int_Weights(:)        &
                                                 / Sin_Square(:)            )

    ! SUM( TdP_TP_Factor(:,lm_loc,lpmp_loc) * Cotan_Val(:) / Sin_Square(:)     )
    TP_TP_Integrals( lm_loc, lpmp_loc, 14 ) = SUM( Ylm_dp( :, lpmp_loc )    &
                                                 * Ylm_CC( :, lpmp_loc )    &
                                                 * TP_Int_Weights(:)        &
                                                 * Cotan_Val(:)             &
                                                 / Sin_Square(:)            )

    ! SUM( TdP_TdP_Factor(:,lm_loc,lpmp_loc)/Sin_Square(:) )
    TP_TP_Integrals( lm_loc, lpmp_loc, 15 ) = SUM( Ylm_dp( :, lpmp_loc )    &
                                                 * Ylm_CC_dp( :, lpmp_loc ) &
                                                 * TP_Int_Weights(:)        &
                                                 / Sin_Square(:)            )

    ! SUM( dTP_TP_Factor(:,lm_loc,lpmp_loc) * Cotan_Val(:) )
    TP_TP_Integrals( lm_loc, lpmp_loc, 16 ) = SUM( Ylm_dt(:, lm_loc )       &
                                                 * Ylm_CC( :, lpmp_loc )    &
                                                 * TP_Int_Weights(:)        &
                                                 * Cotan_Val(:)             )


END DO
END DO


DEALLOCATE( Ylm_Table)


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




ALLOCATE( Cur_R_locs(1:Int_R_Deg)  )
ALLOCATE( R_Square(1:Int_R_Deg)    )



IF ( Matrix_Format == 'Full') THEN

    IF ( Verbose_Flag ) THEN
        PRINT*,"-Initializing Laplace Matrix.  Format: Full. "
    END IF


    DO l = 0,L_LIMIT
    DO re = 0,NUM_R_ELEMENTS-1

        L_Lp1 = REAL( l*(l+1), idp )

        DR   = rlocs(re+1) - rlocs(re)
        TODR = 2.0_idp/DR

        CUR_R_LOCS(:) = (DR/2.0_idp) * (Int_R_Locs(:)+1.0_idp) + rlocs(re)
        R_SQUARE = CUR_R_LOCS**2

        DO dp = 0,DEGREE
        DO d = 0,DEGREE

            i = FP_FEM_Node_Map(re,d)
            j = FP_FEM_Node_Map(re,dp)

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

    IF ( Verbose_Flag ) THEN
        PRINT*,"-Initializing Laplace Matrix.  Format: CCS. "
    END IF



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


                    Laplace_Matrix_VAL(Here,l,:) = Laplace_Matrix_VAL(Here,l,:)     &
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


    MCF_Flag = 0
    Laplace_Factored_VAL = Laplace_Matrix_VAL
    Laplace_Factored_ROW = Laplace_Matrix_ROW
    Laplace_Factored_COL = Laplace_Matrix_COL


END IF



!PRINT*,"Row"
!PRINT*,-Laplace_Matrix_VAL
!PRINT*,"col"
!PRINT*,Laplace_Matrix_ROW
!PRiNT*,"Val"
!PRINT*,Laplace_Matrix_COL




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

COMPLEX(idp), DIMENSION(0:DEGREE)                       ::  Reusable_Values


REAL(idp)                                               ::  DR, TODR



Beta_MVL_Banded = 0.0_idp



ALLOCATE( Cur_R_Locs(1:Int_R_Deg), R_Square(1:Int_R_Deg) )

ALLOCATE( RR_Factor(   1:Int_R_Deg, 0:DEGREE, 0:DEGREE )    )
ALLOCATE( RDR_Factor(  1:Int_R_Deg, 0:DEGREE, 0:DEGREE )    )
ALLOCATE( DRR_Factor(  1:Int_R_Deg, 0:DEGREE, 0:DEGREE )    )
ALLOCATE( DRDR_Factor( 1:Int_R_Deg, 0:DEGREE, 0:DEGREE )    )


IF ( Verbose_Flag ) THEN
    PRINT*,"-Initializing Modified Vector Laplacian Matrix.  Format: Banded. "
END IF



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
            i = Beta_Bandwidth + FP_Beta_Array_Map(re,dp,ui,l,m)

            Beta_MVL_Banded(i-j,j)                   &
                        = Beta_MVL_Banded(i-j,j)     &
                          + Reusable_Values(dp)


        END DO ! dp Loop
    END DO  ! m Loop


END DO  ! l Loop
END DO  ! ui Loop
END DO  ! d Loop
END DO  ! re Loop





! At this point Laplace Matrix_Beta contains the discretized Laplace operator for the three
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




!Call Output_Laplace_Beta(Beta_MVL_Banded,Beta_Prob_Dim, Beta_Prob_Dim)
Beta_Factorized_Flag = .FALSE.


DEALLOCATE( Cur_R_Locs, R_Square )

DEALLOCATE( RR_Factor    )
DEALLOCATE( RDR_Factor   )
DEALLOCATE( DRR_Factor   )
DEALLOCATE( DRDR_Factor  )


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
!                   Calc_Beta1_Terms                                       !
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

COMPLEX(idp), DIMENSION( 1:LM_Length, 1:LM_Length, 1:16),       INTENT(IN)  :: TP_TP_Integrals


REAL(idp),    DIMENSION( 1:Int_R_Deg ),                 INTENT(IN)  :: Cur_R_Locs
REAL(idp),    DIMENSION( 1:Int_R_Deg ),                 INTENT(IN)  :: R_Square


INTEGER                                                     :: d, ui, uj
INTEGER                                                     :: row, col
INTEGER                                                     :: lm_loc, lpmp_loc



ui       = 1
lpmp_loc = FP_LM_Map(lp,mp)
Row      = Beta_Bandwidth + FP_Beta_Array_Map(re,dp,ui,lpmp_loc)




DO d = 0,Degree

    uj = 1

    DO lm_loc = 1,LM_Length
        Col = FP_Beta_Array_Map(re,d,uj,lm_loc)


        Beta_MVL_Banded(Row-Col, Col) = Beta_MVL_Banded(Row-Col, Col)               &
                                      - SUM( dRdR_Factor(:, d, dp) )/3.0_idp        &
                                        * TP_TP_Integrals( lm_loc, lpmp_loc, 1)     &
                                      - 8.0_idp/3.0_idp * SUM( RR_Factor(:, d, dp)  &
                                                               / R_Square(:)    )   &
                                        * TP_TP_Integrals( lm_loc, lpmp_loc, 1)

    END DO ! lpmp_loc Loop






    uj = 2

    DO lm_loc = 1,LM_Length
        Col = FP_Beta_Array_Map(re,d,uj,lm_loc)

        Beta_MVL_Banded(Row-Col, Col) = Beta_MVL_Banded(Row-Col, Col)               &
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


        Beta_MVL_Banded(Row-Col, Col) = Beta_MVL_Banded(Row-Col, Col)               &
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

COMPLEX(idp), DIMENSION( 1:LM_Length, 1:LM_Length, 1:16),       INTENT(IN)  :: TP_TP_Integrals


REAL(idp),    DIMENSION( 1:Int_R_Deg ),                 INTENT(IN)  :: Cur_R_Locs
REAL(idp),    DIMENSION( 1:Int_R_Deg ),                 INTENT(IN)  :: R_Square


INTEGER                                                     :: d, ui, uj
INTEGER                                                     :: row, col
INTEGER                                                     :: lm_loc, lpmp_loc



ui       = 2
lpmp_loc = FP_LM_Map(lp,mp)
Row      = Beta_Bandwidth + FP_Beta_Array_Map(re,dp,ui,lpmp_loc)


DO d = 0,Degree

uj = 1
DO lm_loc = 1,LM_Length
    Col = FP_Beta_Array_Map(re,d,uj,lm_loc)
    

    Beta_MVL_Banded(Row-Col, Col) = Beta_MVL_Banded(Row-Col, Col)               &
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

    Beta_MVL_Banded(Row-Col, Col) = Beta_MVL_Banded(Row-Col, Col)               &
                                  - SUM( RR_Factor(:, d, dp)                    &
                                         /(3.0_idp*R_Square(:) )       )        &
                                    * TP_TP_Integrals( lm_loc, lpmp_loc, 7 )    &
                                  - SUM( RR_Factor(:, d, dp)                    &
                                         /(3.0_idp*R_Square(:) )           )    &
                                    * TP_TP_Integrals( lm_loc, lpmp_loc, 8 )    &
                                  + SUM( 2.0_idp * dRR_Factor(:, d, dp)         &
                                          /Cur_R_Locs(:)                    )   &
                                    * TP_TP_Integrals( lm_loc, lpmp_loc, 1 )    &
                                  - SUM( RR_Factor(:, d, dp)                    &
                                         / (3.0_idp*R_Square(:) )           )   &
                                    * TP_TP_Integrals( lm_loc, lpmp_loc, 9 )    &
                                  + SUM( RR_Factor(:, d, dp)                    &
                                         / R_Square(:)                     )   &
                                    * TP_TP_Integrals( lm_loc, lpmp_loc, 10 )
    

END DO ! lpmp_loc Loop





uj = 3 ! beta^phi

DO lm_loc = 1,LM_Length
    Col = FP_Beta_Array_Map(re,d,uj,lm_loc)

    Beta_MVL_Banded(Row-Col, Col) = Beta_MVL_Banded(Row-Col, Col)           &
                                  - SUM( RR_Factor(:, d, dp)                    &
                                         / (3.0_idp*R_Square(:) )           )   &
                                    * TP_TP_Integrals( lm_loc, lpmp_loc, 11 )   &
                                  - SUM ( 2.0_idp * RR_Factor(:, d, dp)         &
                                          / R_Square(:)                     )   &
                                    * TP_TP_Integrals( lm_loc, lpmp_loc, 12 )

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

COMPLEX(idp), DIMENSION( 1:LM_Length, 1:LM_Length, 1:16),       INTENT(IN)    :: TP_TP_Integrals

REAL(idp),    DIMENSION( 1:Int_R_Deg ),                 INTENT(IN)    :: Cur_R_Locs
REAL(idp),    DIMENSION( 1:Int_R_Deg ),                 INTENT(IN)    :: R_Square


INTEGER                                                     :: d, ui, uj
INTEGER                                                     :: row, col
INTEGER                                                     :: lm_loc, lpmp_loc


lpmp_loc = FP_LM_Map(lp,mp)

ui       = 3
lpmp_loc = FP_LM_Map(lp,mp)
Row      = Beta_Bandwidth + FP_Beta_Array_Map(re,dp,ui,lpmp_loc)



DO d = 0,Degree

uj = 1

DO lm_loc = 1,LM_Length

    Col = FP_Beta_Array_Map(re,d,uj,lm_loc)

    Beta_MVL_Banded(Row-Col, Col) = Beta_MVL_Banded(Row-Col, Col)                   &
                                  - SUM( RdR_Factor(:, d, dp)                       &   ! Term 1
                                        /(3.0_idp * R_Square(:) )    )              &
                                    * TP_TP_Integrals( lm_loc, lpmp_loc, 6 )        &
                                  + SUM( 8.0_idp * RR_Factor(:, d, dp)              &   ! Term 2
                                        /(3.0_idp*Cur_R_Locs(:) * R_Square(:) ) )   &
                                    * TP_TP_Integrals( lm_loc, lpmp_loc, 13 )

END DO ! lm_loc Loop



uj = 2

DO lm_loc = 1,LM_Length
    Col = FP_Beta_Array_Map(re,d,uj,lm_loc)
    



    Beta_MVL_Banded(Row-Col, Col) = Beta_MVL_Banded(Row-Col, Col)               &
                                  - SUM( RR_Factor(:, d, dp)                    &   ! Term 1
                                         /( 3.0_idp * R_Square(:) )   )         &
                                    * TP_TP_Integrals( lm_loc, lpmp_loc, 14 )   &
                                  + SUM( 8.0_idp * RR_Factor(:, d, dp)          &   ! Term 2
                                         / ( 3.0_idp * R_Square(:) ) )          &
                                    * TP_TP_Integrals( lm_loc, lpmp_loc, 15 )

END DO ! lpmp_loc Loop



uj = 3 ! beta^phi
DO lm_loc = 1,LM_Length
    Col = FP_Beta_Array_Map(re,d,uj,lm_loc)

    Beta_MVL_Banded(Row-Col, Col) = Beta_MVL_Banded(Row-Col, Col)           &
                                  - SUM( RR_Factor(:, d, dp )                   &   ! Term 1
                                         /(3.0_idp * R_Square(:) )    )         &
                                    * TP_TP_Integrals( lm_loc, lpmp_loc, 16 )   &
                                  + SUM( 2.0_idp * dRR_Factor(:, d, dp)         &   ! Term 2
                                         / Cur_R_Locs(:)                    )   &
                                    * TP_TP_Integrals( lm_loc, lpmp_loc, 1)     &
                                  + SUM ( 2.0_idp * RR_Factor(:, d, dp)         &   ! Term 3
                                            / R_Square(:)                   )   &
                                    * TP_TP_Integrals( lm_loc, lpmp_loc, 7 )

END DO ! lpmp_loc Loop


END DO ! dp Loop








END SUBROUTINE Calc_Beta3_Terms

END MODULE FP_Intialize_Matrices

