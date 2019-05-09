   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Jacobian_Internal_Functions_Module                                           !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!        Contains the functions used to calculate the components of the          !##!
!##!    extended Jacobian matrix as well as the derivative coefficients. These      !##!
!##!    are used to construct the stiffness matrix.                                 !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   Initialize_Guess_Values                                             !##!
!##!    +102+   Initialize_Special_Guess_Values                                     !##!
!##!                                                                                !##!
!##!    +201+   Initialize_Ylm_Table                                                !##!
!##!    +202+   Initialize_Lagrange_Poly_Tables                                     !##!
!##!                                                                                !##!
!##!    +301+   JCBN_mu_FUNCTION_3D_ALL                                             !##!
!##!    +302+   JCBN_BIGK_FUNCTION                                                  !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!




!*D*================================!
!                                   !
!           Dependencies            !
!                                   !
!===================================!
USE Poseidon_Constants_Module, &
                        ONLY : idp, pi,                 &
                               TwoThirds,               &
                               FourThirds,              &
                               Speed_of_Light,          &
                               GR_Source_Scalar


USE Poseidon_Parameters, &
                        ONLY :  DOMAIN_DIM,             &
                                DEGREE,                 &
                                L_LIMIT,                &
                                NUM_BLOCK_THETA_ROWS,   &
                                NUM_T_ELEMS_PER_BLOCK,  &
                                NUM_P_ELEMS_PER_BLOCK,  &
                                NUM_R_QUAD_POINTS,      &
                                NUM_T_QUAD_POINTS,      &
                                NUM_P_QUAD_POINTS

USE CHIMERA_Parameters, &
                        ONLY :  myID,                   &
                                Analytic_Solution,      &
                                Shift_Solution


USE Global_Variables_And_Parameters, &
                        ONLY :  NUM_R_ELEMENTS,         &
                                NUM_T_ELEMENTS,         &
                                NUM_P_ELEMENTS,         &
                                rlocs, tlocs, plocs,    &
                                NUM_R_NODES,            &
                                INT_R_LOCATIONS,        &
                                INT_T_LOCATIONS,        &
                                INT_P_LOCATIONS,        &
                                VAR_DIM,                &
                                PROB_DIM,               &
                                ULM_LENGTH,             &
                                Ylm_Table_Block,        &
                                Ylm_Values,             &
                                Ylm_CC_Values,          &
                                Ylm_dt_Values,          &
                                Ylm_dp_Values,          &
                                Ylm_dtt_Values,         &
                                Ylm_dtp_Values,         &
                                Ylm_dpp_Values,         &
                                LM_LENGTH,              &
                                myID_Shell,             &
                                Lagrange_Poly_Table,    &
                                LPT_LPT,                &
                                Coefficient_Vector,     &
                                LM_Location,            &
                                Matrix_Location,        &
                                M_VALUES,               &
                                R_Inner,                &
                                R_Outer



USE Additional_Functions_Module, &
                        ONLY :  Norm_Factor,                            &
                                Legendre_Poly,                          &
                                Lagrange_Poly,                          &
                                Lagrange_Poly_Deriv,                    &
                                Lagrange_Second_Deriv,                  &
                                Map_From_X_Space,                       &
                                Initialize_LGL_Quadrature_Locations,    &
                                Initialize_LG_Quadrature,               &
                                Initialize_LGL_Quadrature,              &
                                Map_To_X_SPace,                         &
                                Spherical_Harmonic,                     &
                                CFA_Matrix_Map




IMPLICIT NONE


CONTAINS



!+101+###########################################################################!
!                                                                                !
!                  Initialize_Guess_Values          !
!                                                                                !
!################################################################################!
SUBROUTINE Initialize_Guess_Values()


INTEGER                                     ::  beta_i, j

REAL(KIND = idp)                            ::  Beta_Start,     &
                                                delta_Beta



INTEGER                                     ::  CUR_PSI_LOC,        &
                                                CUR_ALPHPSI_LOC,    &
                                                CUR_BETA_LOC


INTEGER                                     :: re, rd



!
!   Empty Space Initial Guess
!
Coefficient_Vector(:) = 0.0_idp


Beta_Start = 0.0_idp
!Beta_Start = 0.05000_idp*2.0_idp*sqrt(pi)


delta_Beta = 0.0_idp
!delta_Beta = 0.0001_idp
!delta_Beta = Beta_Start/VAR_DIM







DO re = 0,NUM_R_ELEMENTS - 1


    DO rd = 0,DEGREE


        ! 2 sqrt(pi) is Ylm normalization factor


        CUR_PSI_LOC = CFA_Matrix_Map( 1, 0, 0, re, rd )
        Coefficient_Vector(CUR_PSI_LOC) = 0.9_idp * 2.0_idp * sqrt(pi)




        CUR_ALPHPSI_LOC = CFA_Matrix_Map( 2, 0, 0, re, rd )
        Coefficient_Vector(CUR_ALPHPSI_LOC) = 0.9_idp * 2.0_idp * sqrt(pi)


        DO beta_i = 1,DOMAIN_DIM

            CUR_BETA_LOC = CFA_Matrix_Map( 2+beta_i, 0, 0, re, rd )
            Coefficient_Vector(CUR_BETA_LOC) = Beta_Start - delta_Beta*(re*DEGREE+rd)

        END DO


    END DO


END DO




END SUBROUTINE Initialize_Guess_Values








!+102+###########################################################################!
!                                                                                !
!                  Initialize_Guess_Values          !
!                                                                                !
!################################################################################!
SUBROUTINE Initialize_Special_Guess_Values()


INTEGER                                 :: re, rd, Here


REAL(KIND = idp), DIMENSION(0:DEGREE)                           :: R_Values


REAL(KIND = idp), DIMENSION(0:DEGREE)                           ::  Local_Locations


INTEGER                                                         ::  CUR_PSI_LOC,    &
                                                                    CUR_ALPHPSI_LOC,&
                                                                    CUR_SHIFT_LOC

INTEGER          :: Matrix_Map

REAL(KIND = idp)  :: csqr

Local_Locations = Initialize_LGL_Quadrature_Locations(DEGREE)


!
!   Empty Space Initial Guess
!
Coefficient_Vector = 0.0_idp


csqr = Speed_of_Light*Speed_of_Light



DO re = 0,NUM_R_ELEMENTS - 1

    R_Values = Map_From_X_Space(rlocs(re), rlocs(re + 1), Local_Locations)

    DO rd = 0,DEGREE
 
        ! 2 sqrt(pi) is Ylm normalization factor


        CUR_PSI_LOC = (re*Degree + rd)*(ULM_LENGTH)

        Coefficient_Vector(CUR_PSI_LOC) = 2.0_idp * sqrt(pi)                                       &
                                        * ( 1.0_idp - 0.5_idp                                          &
                                            * Analytic_Solution(R_Values(rd),0.0_idp,0.0_idp)/csqr   )



        CUR_ALPHPSI_LOC = (re*Degree + rd)*(ULM_LENGTH) + 1

        Coefficient_Vector(CUR_ALPHPSI_LOC) = 2.0_idp * sqrt(pi)                             &
                                        * ( 1.0_idp + 0.5_idp                                   &
                                            * Analytic_Solution(R_Values(rd),0.0_idp,0.0_idp)/csqr  )




        CUR_SHIFT_LOC = (re*Degree + rd)*ULM_LENGTH + 2

        Coefficient_Vector(Cur_Shift_Loc) = 2.0_idp*sqrt(pi)*Shift_Solution(rlocs(re),rlocs,NUM_R_ELEMENTS)


    END DO


END DO



END SUBROUTINE Initialize_Special_Guess_Values












!+201+###########################################################################!
!                                                                                !
!                  Init_Ylm_Table_Block          !
!                                                                                !
!################################################################################!
SUBROUTINE Initialize_Ylm_Table


INTEGER                                                         ::  l, m, te, pe, td, pd


REAL(KIND = idp), DIMENSION(1:NUM_T_QUAD_POINTS)                ::  T_Locations
REAL(KIND = idp), DIMENSION(1:NUM_P_QUAD_POINTS)                ::  P_Locations

REAL(KIND = idp), DIMENSION(1:NUM_T_QUAD_POINTS)                ::  Legendre_Poly_Values


REAL(KIND = idp)                                                ::  Norm_Storage

INTEGER                                                         ::  Block_T_Begin, Block_P_Begin,   &
                                                                    Global_TE, Global_PE


Ylm_Table_Block = 0.0_idp

Block_T_Begin = MOD(myID_Shell,NUM_BLOCK_THETA_ROWS)*NUM_T_ELEMS_PER_BLOCK
Block_P_Begin = (myID_Shell/NUM_BLOCK_THETA_ROWS)*NUM_P_ELEMS_PER_BLOCK

!PRINT*,"ylm T/P_Begin",Block_T_Begin, Block_P_Begin

DO l = 0,L_LIMIT

    DO m = -l,l

        Norm_Storage = Norm_Factor(l,m)


        DO pe = 0,NUM_P_ELEMS_PER_BLOCK-1

             !                                                          !
            !!    Map Phi Locations from [-1,1] space to real space.    !!
             !                                                          !
            Global_pe = Block_P_Begin + pe
            P_Locations = Map_From_X_Space(plocs(Global_pe), plocs(Global_pe + 1), INT_P_LOCATIONS)

!            print*,"Ylm_Block P_locs",P_Locations

            DO te = 0,NUM_T_ELEMS_PER_BLOCK-1

                 !                                                          !
                !!   Map Theta Locations from [-1,1] space to real space.   !!
                 !                                                          !
                Global_te = Block_T_Begin + te
                T_Locations = Map_From_X_Space(tlocs(Global_te), tlocs(Global_te + 1), INT_T_LOCATIONS)

!                PRINT*,"ylm_Block T_locs",T_Locations

                Legendre_Poly_Values = Legendre_Poly(l, m, NUM_T_QUAD_POINTS, T_Locations)



                DO pd = 1,NUM_P_QUAD_POINTS

                    DO td = 1,NUM_T_QUAD_POINTS





                        Ylm_Table_Block(m, l, td, pd, te, pe) = Norm_Storage                            &   ! Normalization Factor
                                                            * Legendre_Poly_Values(td)                  &   ! Legendre Polynomial
                                                            * CDEXP(CMPLX(0.0_idp,m*P_Locations(pd),idp))   ! exp(im phi)



                    END DO

                END DO

            END DO

        END DO

    END DO

END DO



END SUBROUTINE Initialize_Ylm_Table








!+201+###########################################################################!
!                                                                                !
!                  Init_Ylm_Table_Block          !
!                                                                                !
!################################################################################!
SUBROUTINE Initialize_Ylm_Tables


INTEGER                                                         ::  l, m, te, pe, td, pd


REAL(KIND = idp), DIMENSION(1:NUM_T_QUAD_POINTS)                ::  T_Locations
REAL(KIND = idp), DIMENSION(1:NUM_P_QUAD_POINTS)                ::  P_Locations

REAL(KIND = idp), DIMENSION(1:1)                                ::  Legendre_Poly_Value


REAL(KIND = idp)                                                ::  Norm_Storage
REAL(KIND = idp)                                                ::  REAL_L
INTEGER                                                         ::  Block_T_Begin,      &
                                                                    Block_P_Begin,      &
                                                                    Global_TE,          &
                                                                    Global_PE

REAL(KIND = idp), DIMENSION(-L_LIMIT:L_LIMIT)                           :: M_POWER_TABLE

INTEGER                                                         ::  lm_loc

INTEGER, DIMENSION(1:NUM_T_QUAD_POINTS)                         ::  here

COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:,:,:,:)        ::  Ylm_Table

REAL(KIND = idp), DIMENSION(0:LM_LENGTH-1)                      ::  Sqrt_Term
REAL(KIND = idp), DIMENSION(1:NUM_T_QUAD_POINTS)                ::  SIN_VAL, TAN_VAL
REAL(KIND = idp), DIMENSION(1:NUM_T_QUAD_POINTS)                ::  CSC_VAL, COT_VAL


ALLOCATE( Ylm_Table(-L_LIMIT:L_LIMIT, -1:L_LIMIT,                           &
                    1:NUM_T_QUAD_POINTS, 1:NUM_P_QUAD_POINTS,               &
                    0:NUM_T_ELEMS_PER_BLOCK-1, 0:NUM_P_ELEMS_PER_BLOCK-1)   )



Block_T_Begin = MOD(myID_Shell,NUM_BLOCK_THETA_ROWS)*NUM_T_ELEMS_PER_BLOCK
Block_P_Begin = (myID_Shell/NUM_BLOCK_THETA_ROWS)*NUM_P_ELEMS_PER_BLOCK

!PRINT*,"ylm T/P_Begin",Block_T_Begin, Block_P_Begin

M_POWER_TABLE(0) = 1.0_idp
IF ( L_LIMIT > 0 ) THEN
    DO m = 1, L_LIMIT

        M_POWER_TABLE(m) = -1.0_idp*M_POWER_TABLE(m-1)
        M_POWER_TABLE(-m) = -1.0_idp*M_POWER_TABLE(m-1)

    END DO
END IF


Sqrt_Term(0) = 0.0_idp
DO l = 1,L_LIMIT

    REAL_L = REAL(l, idp)

    DO m = -M_VALUES(l),M_VALUES(l)

        Sqrt_Term(LM_Location(l,m)) = SQRT((2*REAL_L + 1)/(2*REAL_L - 1)* REAL( (l-m)*(l+m), idp) )


    END DO

END DO


Ylm_Table = 0.0_idp

DO pe = 0,NUM_P_ELEMS_PER_BLOCK-1

     !                                                          !
    !!    Map Phi Locations from [-1,1] space to real space.    !!
     !                                                          !
    Global_pe = Block_P_Begin + pe
    P_Locations = Map_From_X_Space(plocs(Global_pe), plocs(Global_pe + 1), INT_P_LOCATIONS)


    DO te = 0,NUM_T_ELEMS_PER_BLOCK-1

         !                                                          !
        !!   Map Theta Locations from [-1,1] space to real space.   !!
         !                                                          !
        Global_te = Block_T_Begin + te
        T_Locations = Map_From_X_Space(tlocs(Global_te), tlocs(Global_te + 1), INT_T_LOCATIONS)


        DO pd = 1,NUM_P_QUAD_POINTS

            DO td = 1,NUM_T_QUAD_POINTS

                DO l = 0,L_LIMIT

                    DO m = -M_VALUES(l),M_VALUES(l)

                        Norm_Storage = Norm_Factor(l,m)
                        Legendre_Poly_Value = Legendre_Poly(l, m, 1, T_Locations(td))

                        Ylm_Table(m, l, td, pd, te, pe) = Norm_Storage                            &   ! Normalization Factor
                                                        * Legendre_Poly_Value(1)                  &   ! Legendre Polynomial
                                                        * CDEXP(CMPLX(0.0_idp,m*P_Locations(pd),idp))   ! exp(im phi)


                    END DO

                END DO

            END DO

        END DO

    END DO

END DO



DO pe = 0,NUM_P_ELEMS_PER_BLOCK-1

     !                                                          !
    !!    Map Phi Locations from [-1,1] space to real space.    !!
     !                                                          !
    Global_pe = Block_P_Begin + pe
    P_Locations = Map_From_X_Space(plocs(Global_pe), plocs(Global_pe + 1), INT_P_LOCATIONS)


    DO te = 0,NUM_T_ELEMS_PER_BLOCK-1

         !                                                          !
        !!   Map Theta Locations from [-1,1] space to real space.   !!
         !                                                          !
        Global_te = Block_T_Begin + te
        T_Locations = Map_From_X_Space(tlocs(Global_te), tlocs(Global_te + 1), INT_T_LOCATIONS)

        CSC_VAL(:) = 1.0_idp/DSIN(T_Locations(:))
        COT_VAL(:) = 1.0_idp/DTAN(T_Locations(:))






        DO pd = 1,NUM_P_QUAD_POINTS

            DO td = 1,NUM_T_QUAD_POINTS

                DO l = 0,L_LIMIT

                    REAL_L = REAL(l, idp)

                    DO m = -M_VALUES(l),M_VALUES(l)



                        

                        lm_loc = LM_Location(l,m)
                        Norm_Storage = Norm_Factor(l,m)
                        Legendre_Poly_Value = Legendre_Poly(l, m, 1, T_Locations(td))

                        Ylm_Values(lm_loc, td, pd, te, pe) = Ylm_Table(m,l,td, pd, te, pe)
                        Ylm_CC_Values(lm_loc, td, pd, te, pe) = M_POWER_TABLE(m) * Ylm_Table(-m,l,td, pd, te, pe)

                        Ylm_dt_Values(lm_loc, td, pd, te, pe) = REAL_L*COT_VAL(td)                      &
                                                                * Ylm_Table(m,l,td, pd, te, pe)         &
                                                              - SQRT_TERM(lm_loc) * CSC_VAL(td)         &
                                                                * Ylm_Table(m,l-1,td, pd, te, pe)

                        Ylm_dp_Values(lm_loc, td, pd, te, pe) = CMPLX(0,m,idp)                          &
                                                                * Ylm_Table(m,l,td, pd, te, pe)

                        Ylm_dtt_Values(lm_loc, td, pd, te, pe) = (REAL(m*m - l, idp) * (CSC_VAL(td)*CSC_VAL(td))    &
                                                                  - REAL_L*REAL_L) * Ylm_Table(m,l,td, pd, te, pe)  &
                                                               + SQRT_TERM(lm_loc) *( CSC_VAL(td)*COT_VAL(td) )     &
                                                                  * Ylm_Table(m,l-1,td, pd, te, pe)

                        Ylm_dtp_Values(lm_loc, td, pd, te, pe) = CMPLX(0,m,idp)*Ylm_dt_Values(lm_loc, td, pd, te, pe)

                        Ylm_dpp_Values(lm_loc, td, pd, te, pe) = CMPLX(0,-m*m,idp)*Ylm_Table(m,l,td, pd, te, pe)




                    END DO

                END DO

            END DO

        END DO

    END DO

END DO












END SUBROUTINE Initialize_Ylm_Tables








!+202+###########################################################################!
!                                                                                !
!                  Init_Lagrange_Poly_Tables          !
!                                                                                !
!################################################################################!
SUBROUTINE Initialize_Lagrange_Poly_Tables()


INTEGER                                         ::  Derv_Degree,    &
                                                    Basis_Func,     &
                                                    Eval_Point

REAL(KIND = idp), DIMENSION(0:DEGREE)           ::  Local_Locations
REAL(KIND = idp), DIMENSION(0:DEGREE)           ::  Lagrange_Poly_Values
REAL(KIND = idp), DIMENSION(0:DEGREE)           ::  Lagrange_DRV_Values
REAL(KIND = idp), DIMENSION(0:DEGREE)           ::  Lagrange_DDRV_Values


INTEGER                                         ::  rd, d, dp,dd


Lagrange_Poly_Table = 0.0_idp


Local_Locations = Initialize_LGL_Quadrature_Locations(DEGREE)

DO Eval_Point = 1,NUM_R_QUAD_POINTS

    Lagrange_Poly_Values = Lagrange_Poly(INT_R_Locations(Eval_Point), DEGREE, Local_Locations)
    Lagrange_DRV_Values  = Lagrange_Poly_Deriv(INT_R_Locations(Eval_Point), DEGREE, Local_Locations)
    Lagrange_DDRV_Values = Lagrange_Second_Deriv(INT_R_Locations(Eval_Point), DEGREE, Local_Locations)


    Lagrange_Poly_Table(:, Eval_Point, 0) = Lagrange_Poly_Values
    Lagrange_Poly_Table(:, Eval_Point, 1) = Lagrange_DRV_Values
    Lagrange_Poly_Table(:, Eval_Point, 2) = Lagrange_DDRV_Values


END DO


DO dd = 0,2
    DO d = 0,DEGREE
        DO dp = 0,DEGREE
            DO rd = 1,NUM_R_QUAD_POINTS
           
                LPT_LPT(rd, dp, d,0:1, dd) = Lagrange_Poly_Table(d, rd, 0:1)       &
                                           * Lagrange_Poly_Table(dp, rd, dd)
            END DO 
        END DO
    END DO
END DO



END SUBROUTINE Initialize_Lagrange_Poly_Tables












!+301+###########################################################################!
!                                                                                !
!                                                                                !
!################################################################################!
PURE FUNCTION JCBN_mu_FUNCTION_3D_ALL(rd, td, pd,                                       &
                                        r, R_SQUARE, R_CUBED, R_INVERSE, RSIN_SQUARE,    &
                                        SIN_VAL, SIN_SQUARE, CSC_SQUARE,                &
                                        COS_VAL, COTAN_VAL,                             &
                                        CUR_VAL_BETA, CUR_DRV_BETA                  )

REAL(KIND = idp), DIMENSION(1:3,1:3)      ::  JCBN_mu_FUNCTION_3D_ALL

INTEGER, INTENT(IN)                     ::  rd, td, pd
REAL(KIND = idp), INTENT(IN)            ::  r, R_SQUARE, RSIN_SQUARE, R_CUBED, R_INVERSE,    &
                                            SIN_VAL, SIN_SQUARE, CSC_SQUARE, COS_VAL, COTAN_VAL





REAL(KIND = idp), INTENT(IN), DIMENSION( 1:3,                   &
                                         1:NUM_R_QUAD_POINTS,   &
                                         1:NUM_T_QUAD_POINTS,   &
                                         1:NUM_P_QUAD_POINTS    )          ::  CUR_VAL_BETA


REAL(KIND = idp), INTENT(IN), DIMENSION( 1:3, 1:3,              &
                                         1:NUM_R_QUAD_POINTS,   &
                                         1:NUM_T_QUAD_POINTS,   &
                                         1:NUM_P_QUAD_POINTS    )          ::  CUR_DRV_BETA


REAL(KIND = idp)                    ::  Beta_Drv_Trace





Beta_Drv_Trace = CUR_DRV_BETA(1, 1, rd, td, pd )    &
               + CUR_DRV_BETA(2, 2, rd, td, pd )    &
               + CUR_DRV_BETA(3, 3, rd, td, pd )






JCBN_mu_FUNCTION_3D_ALL(1,1) = 2.0_idp                                      &
                                * CUR_DRV_BETA(1, 1, rd, td, pd )           &
                             - FourThirds/r                                 &
                                * CUR_VAL_BETA(1, rd, td, pd )              &
                             - TwoThirds * COTAN_VAL                        &
                                * CUR_VAL_BETA(2, rd, td, pd )              &
                             - TwoThirds                                    &
                                * Beta_Drv_Trace

JCBN_mu_FUNCTION_3D_ALL(2,1) =  CUR_DRV_BETA(1, 2, rd, td, pd )             &
                             + (1.0_idp/R_SQUARE)                           &
                                * CUR_DRV_BETA(2, 1, rd, td, pd )

JCBN_mu_FUNCTION_3D_ALL(3,1) =  CUR_DRV_BETA(1, 3, rd, td, pd )             &
                             + (CSC_SQUARE/R_SQUARE)                        &
                                * CUR_DRV_BETA(3, 1, rd, td, pd )





JCBN_mu_FUNCTION_3D_ALL(1,2) = JCBN_mu_FUNCTION_3D_ALL(2,1)

JCBN_mu_FUNCTION_3D_ALL(2,2) = FourThirds / R_CUBED                          &
                                * CUR_VAL_BETA(1, rd, td, pd )              &
                             - (TwoThirds/R_SQUARE)*COTAN_VAL               &
                                * CUR_VAL_BETA(2, rd, td, pd )              &
                             + 2 / R_SQUARE                                 &
                                * CUR_DRV_BETA(2, 2, rd, td, pd )           &
                            - TwoThirds / R_SQUARE                          &
                                * Beta_Drv_Trace

JCBN_mu_FUNCTION_3D_ALL(3,2) = (1.0_idp/R_SQUARE)                           &
                                * CUR_DRV_BETA(2, 3, rd, td, pd )           &
                             + (CSC_SQUARE/R_SQUARE)                        &
                                * CUR_DRV_BETA(3, 2, rd, td, pd )



JCBN_mu_FUNCTION_3D_ALL(1,3) = JCBN_mu_FUNCTION_3D_ALL(3,1)

JCBN_mu_FUNCTION_3D_ALL(2,3) = JCBN_mu_FUNCTION_3D_ALL(3,2)

JCBN_mu_FUNCTION_3D_ALL(3,3) = TwoThirds*CSC_SQUARE/( r * R_SQUARE)             &
                                * CUR_VAL_BETA(1, rd, td, pd )                  &
                             +  FourThirds * COTAN_VAL*CSC_SQUARE/R_SQUARE      &
                                * CUR_VAL_BETA(2, rd, td, pd )                  &
                             + 2.0_idp*CSC_SQUARE / R_SQUARE                    &
                                * CUR_DRV_BETA(3, 3, rd, td, pd )               &
                             - TwoThirds / RSIN_SQUARE                          &
                                * Beta_Drv_Trace






END FUNCTION JCBN_mu_FUNCTION_3D_ALL








!+401+###########################################################################!
!                                                                                !
!                                                                                !
!################################################################################!
PURE FUNCTION JCBN_BIGK_FUNCTION(   rd, td, pd,                         &
                                    CUR_VAL_BETA, CUR_DRV_BETA,         &
                                    R_VAL, R_SQUARE, SIN_SQUARE,  CSC_SQUARE,        &
                                    RSIN_SQUARE, COTAN_VAL              )

REAL(KIND = idp)                    ::  JCBN_BIGK_FUNCTION

INTEGER, INTENT(IN)                 ::  rd, td, pd

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:3,                   &
                                         1:NUM_R_QUAD_POINTS,   &
                                         1:NUM_T_QUAD_POINTS,   &
                                         1:NUM_P_QUAD_POINTS    )   ::  CUR_VAL_BETA


REAL(KIND = idp), INTENT(IN), DIMENSION( 1:3, 1:3,              &
                                         1:NUM_R_QUAD_POINTS,   &
                                         1:NUM_T_QUAD_POINTS,   &
                                         1:NUM_P_QUAD_POINTS    )   ::  CUR_DRV_BETA


REAL(KIND = idp), INTENT(IN)                                    ::  R_VAL,                  &
                                                                    R_SQUARE,               &
                                                                    SIN_SQUARE,             &
                                                                    CSC_SQUARE,             &
                                                                    RSIN_SQUARE,            &
                                                                    COTAN_VAL




REAL(KIND = idp)                    ::  FourThirds,                 &
                                        EightThirds




FourThirds = (4.0_idp/3.0_idp)
EightThirds = (8.0_idp/3.0_idp)

JCBN_BIGK_FUNCTION =                                                                                        &
        FourThirds                      * CUR_DRV_BETA(1,1, rd, td, pd) * CUR_DRV_BETA(1,1, rd, td, pd)     &
      + R_SQUARE                        * CUR_DRV_BETA(1,2, rd, td, pd) * CUR_DRV_BETA(1,2, rd, td, pd)     &
      + 1.0_idp/R_SQUARE                * CUR_DRV_BETA(2,1, rd, td, pd) * CUR_DRV_BETA(2,1, rd, td, pd)     &
      + 2.0_idp                         * CUR_DRV_BETA(1,2, rd, td, pd) * CUR_DRV_BETA(2,1, rd, td, pd)     &
      - FourThirds                      * CUR_DRV_BETA(1,1, rd, td, pd) * CUR_DRV_BETA(2,2, rd, td, pd)     &
      + FourThirds                      * CUR_DRV_BETA(2,2, rd, td, pd) * CUR_DRV_BETA(2,2, rd, td, pd)     &
      + SIN_SQUARE                      * CUR_DRV_BETA(2,3, rd, td, pd) * CUR_DRV_BETA(2,3, rd, td, pd)     &
      + CSC_SQUARE                      * CUR_DRV_BETA(3,2, rd, td, pd) * CUR_DRV_BETA(3,2, rd, td, pd)     &
      + 2.0_idp                         * CUR_DRV_BETA(2,3, rd, td, pd) * CUR_DRV_BETA(3,2, rd, td, pd)     &
      - FourThirds                      * CUR_DRV_BETA(2,2, rd, td, pd) * CUR_DRV_BETA(3,3, rd, td, pd)     &
      + FourThirds                      * CUR_DRV_BETA(3,3, rd, td, pd) * CUR_DRV_BETA(3,3, rd, td, pd)     &
      + RSIN_SQUARE                     * CUR_DRV_BETA(1,3, rd, td, pd) * CUR_DRV_BETA(1,3, rd, td, pd)     &
      + 2.0_idp                         * CUR_DRV_BETA(1,3, rd, td, pd) * CUR_DRV_BETA(3,1, rd, td, pd)     &
      + CSC_SQUARE/R_SQUARE             * CUR_DRV_BETA(3,1, rd, td, pd) * CUR_DRV_BETA(3,1, rd, td, pd)     &
      - FourThirds                      * CUR_DRV_BETA(1,1, rd, td, pd) * CUR_DRV_BETA(3,3, rd, td, pd)     &
      +(- EightThirds/R_VAL             * CUR_DRV_BETA(1,1, rd, td, pd)                                     &
      + FourThirds/R_VAL                * CUR_DRV_BETA(2,2, rd, td, pd)                                     &
      + FourThirds/R_VAL                * CUR_DRV_BETA(3,3, rd, td, pd) ) * CUR_VAL_BETA(1, rd, td, pd)     &
      + (- FourThirds*COTAN_VAL         * CUR_DRV_BETA(1,1, rd, td, pd)                                     &
      - FourThirds*COTAN_VAL            * CUR_DRV_BETA(2,2, rd, td, pd)                                     &
      + EightThirds*COTAN_VAL           * CUR_DRV_BETA(3,3, rd, td, pd) ) * CUR_VAL_BETA(2, rd, td, pd)     &
      + FourThirds/R_SQUARE             * CUR_VAL_BETA(1, rd, td, pd) * CUR_VAL_BETA(1, rd, td, pd)         &
      + FourThirds*COTAN_VAL/R_VAL      * CUR_VAL_BETA(1, rd, td, pd) * CUR_VAL_BETA(2, rd, td, pd)         &
      + FourThirds*COTAN_VAL*COTAN_VAL  * CUR_VAL_BETA(2, rd, td, pd) * CUR_VAL_BETA(2, rd, td, pd)




END FUNCTION JCBN_BIGK_FUNCTION








!+401+###########################################################################!
!                                                                                !
!                                                                                !
!################################################################################!
SUBROUTINE  JCBN_BIGK_SUBROUTINE(   rd, td, pd,                         &
                                    CUR_VAL_BETA, CUR_DRV_BETA,         &
                                    R_VAL, R_SQUARE, SIN_SQUARE,  CSC_SQUARE,        &
                                    RSIN_SQUARE, COTAN_VAL              )

REAL(KIND = idp)                    ::  JCBN_BIGK_FUNCTION

INTEGER, INTENT(IN)                 ::  rd, td, pd

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:3,                   &
                                         1:NUM_R_QUAD_POINTS,   &
                                         1:NUM_T_QUAD_POINTS,   &
                                         1:NUM_P_QUAD_POINTS    )   ::  CUR_VAL_BETA


REAL(KIND = idp), INTENT(IN), DIMENSION( 1:3, 1:3,              &
                                         1:NUM_R_QUAD_POINTS,   &
                                         1:NUM_T_QUAD_POINTS,   &
                                         1:NUM_P_QUAD_POINTS    )   ::  CUR_DRV_BETA


REAL(KIND = idp), INTENT(IN)                                    ::  R_VAL,                  &
                                                                    R_SQUARE,               &
                                                                    SIN_SQUARE,             &
                                                                    CSC_SQUARE,             &
                                                                    RSIN_SQUARE,            &
                                                                    COTAN_VAL




REAL(KIND = idp)                    ::  FourThirds,                 &
                                        EightThirds




FourThirds = (4.0_idp/3.0_idp)
EightThirds = (8.0_idp/3.0_idp)

JCBN_BIGK_FUNCTION =                                                                                        &
        FourThirds                      * CUR_DRV_BETA(1,1, rd, td, pd) * CUR_DRV_BETA(1,1, rd, td, pd)     &
      + R_SQUARE                        * CUR_DRV_BETA(1,2, rd, td, pd) * CUR_DRV_BETA(1,2, rd, td, pd)     &
      + 1.0_idp/R_SQUARE                * CUR_DRV_BETA(2,1, rd, td, pd) * CUR_DRV_BETA(2,1, rd, td, pd)     &
      + 2.0_idp                         * CUR_DRV_BETA(1,2, rd, td, pd) * CUR_DRV_BETA(2,1, rd, td, pd)     &
      - FourThirds                      * CUR_DRV_BETA(1,1, rd, td, pd) * CUR_DRV_BETA(2,2, rd, td, pd)     &
      + FourThirds                      * CUR_DRV_BETA(2,2, rd, td, pd) * CUR_DRV_BETA(2,2, rd, td, pd)     &
      + SIN_SQUARE                      * CUR_DRV_BETA(2,3, rd, td, pd) * CUR_DRV_BETA(2,3, rd, td, pd)     &
      + CSC_SQUARE                      * CUR_DRV_BETA(3,2, rd, td, pd) * CUR_DRV_BETA(3,2, rd, td, pd)     &
      + 2.0_idp                         * CUR_DRV_BETA(2,3, rd, td, pd) * CUR_DRV_BETA(3,2, rd, td, pd)     &
      - FourThirds                      * CUR_DRV_BETA(2,2, rd, td, pd) * CUR_DRV_BETA(3,3, rd, td, pd)     &
      + FourThirds                      * CUR_DRV_BETA(3,3, rd, td, pd) * CUR_DRV_BETA(3,3, rd, td, pd)     &
      + RSIN_SQUARE                     * CUR_DRV_BETA(1,3, rd, td, pd) * CUR_DRV_BETA(1,3, rd, td, pd)     &
      + 2.0_idp                         * CUR_DRV_BETA(1,3, rd, td, pd) * CUR_DRV_BETA(3,1, rd, td, pd)     &
      + CSC_SQUARE/R_SQUARE             * CUR_DRV_BETA(3,1, rd, td, pd) * CUR_DRV_BETA(3,1, rd, td, pd)     &
      - FourThirds                      * CUR_DRV_BETA(1,1, rd, td, pd) * CUR_DRV_BETA(3,3, rd, td, pd)     &
      +(- EightThirds/R_VAL             * CUR_DRV_BETA(1,1, rd, td, pd)                                     &
      + FourThirds/R_VAL                * CUR_DRV_BETA(2,2, rd, td, pd)                                     &
      + FourThirds/R_VAL                * CUR_DRV_BETA(3,3, rd, td, pd) ) * CUR_VAL_BETA(1, rd, td, pd)     &
      + (- FourThirds*COTAN_VAL         * CUR_DRV_BETA(1,1, rd, td, pd)                                     &
      - FourThirds*COTAN_VAL            * CUR_DRV_BETA(2,2, rd, td, pd)                                     &
      + EightThirds*COTAN_VAL           * CUR_DRV_BETA(3,3, rd, td, pd) ) * CUR_VAL_BETA(2, rd, td, pd)     &
      + FourThirds/R_SQUARE             * CUR_VAL_BETA(1, rd, td, pd) * CUR_VAL_BETA(1, rd, td, pd)         &
      + FourThirds*COTAN_VAL/R_VAL      * CUR_VAL_BETA(1, rd, td, pd) * CUR_VAL_BETA(2, rd, td, pd)         &
      + FourThirds*COTAN_VAL*COTAN_VAL  * CUR_VAL_BETA(2, rd, td, pd) * CUR_VAL_BETA(2, rd, td, pd)

IF (.TRUE.) THEN
	PRINT*,"BETA1_VAL",CUR_VAL_BETA(1,rd,td,pd)
	PRINT*,"BETA1_DRV",CUR_DRV_BETA(1,1,rd,td,pd)
	PRINT*,"1",FourThirds                      * CUR_DRV_BETA(1,1, rd, td, pd) * CUR_DRV_BETA(1,1, rd, td, pd)     
	PRINT*,"2", R_SQUARE                        * CUR_DRV_BETA(1,2, rd, td, pd) * CUR_DRV_BETA(1,2, rd, td, pd)     
	PRINT*,"3", 1.0_idp/R_SQUARE                * CUR_DRV_BETA(2,1, rd, td, pd) * CUR_DRV_BETA(2,1, rd, td, pd)     
	PRINT*,"4", 2.0_idp                         * CUR_DRV_BETA(1,2, rd, td, pd) * CUR_DRV_BETA(2,1, rd, td, pd)  
	PRINT*,"5",- FourThirds                      * CUR_DRV_BETA(1,1, rd, td, pd) * CUR_DRV_BETA(2,2, rd, td, pd)   
	PRINT*,"6", FourThirds                      * CUR_DRV_BETA(2,2, rd, td, pd) * CUR_DRV_BETA(2,2, rd, td, pd)
	PRINT*,"7",SIN_SQUARE                      * CUR_DRV_BETA(2,3, rd, td, pd) * CUR_DRV_BETA(2,3, rd, td, pd)
	PRINT*,"8",CSC_SQUARE                      * CUR_DRV_BETA(3,2, rd, td, pd) * CUR_DRV_BETA(3,2, rd, td, pd)
	PRINT*,"9",2.0_idp                         * CUR_DRV_BETA(2,3, rd, td, pd) * CUR_DRV_BETA(3,2, rd, td, pd)
	PRINT*,"10",- FourThirds                      * CUR_DRV_BETA(2,2, rd, td, pd) * CUR_DRV_BETA(3,3, rd, td, pd)
	PRINT*,"11",FourThirds                      * CUR_DRV_BETA(3,3, rd, td, pd) * CUR_DRV_BETA(3,3, rd, td, pd)
	PRINT*,"12",RSIN_SQUARE                     * CUR_DRV_BETA(1,3, rd, td, pd) * CUR_DRV_BETA(1,3, rd, td, pd)
	PRINT*,"13",2.0_idp                         * CUR_DRV_BETA(1,3, rd, td, pd) * CUR_DRV_BETA(3,1, rd, td, pd)
	PRINT*,"14",CSC_SQUARE/R_SQUARE             * CUR_DRV_BETA(3,1, rd, td, pd) * CUR_DRV_BETA(3,1, rd, td, pd)
	PRINT*,"15",- FourThirds                      * CUR_DRV_BETA(1,1, rd, td, pd) * CUR_DRV_BETA(3,3, rd, td, pd)   
	PRINT*,"16",- EightThirds/R_VAL             * CUR_DRV_BETA(1,1, rd, td, pd)  * CUR_VAL_BETA(1,rd,td,pd)
	PRINT*,"17", FourThirds/R_VAL                * CUR_DRV_BETA(2,2, rd, td, pd) * CUR_VAL_BETA(1,rd,td,pd)
	PRINT*,"18",FourThirds/R_VAL                * CUR_DRV_BETA(3,3, rd, td, pd)  * CUR_VAL_BETA(1, rd, td, pd)
	PRINT*,"19",- FourThirds*COTAN_VAL         * CUR_DRV_BETA(1,1, rd, td, pd)   * CUR_VAL_BETA(2,rd,td,pd)
	PRINT*,"20",- FourThirds*COTAN_VAL            * CUR_DRV_BETA(2,2, rd, td, pd)    * CUR_VAL_BETA(2,rd,td,pd)
	PRINT*,"21",EightThirds*COTAN_VAL           * CUR_DRV_BETA(3,3, rd, td, pd)  * CUR_VAL_BETA(2, rd, td, pd)     
	PRINT*,"22",FourThirds/R_SQUARE             * CUR_VAL_BETA(1, rd, td, pd) * CUR_VAL_BETA(1, rd, td, pd)         
	PRINT*,"23",FourThirds*COTAN_VAL/R_VAL      * CUR_VAL_BETA(1, rd, td, pd) * CUR_VAL_BETA(2, rd, td, pd)         
	PRINT*,"24",FourThirds*COTAN_VAL*COTAN_VAL  * CUR_VAL_BETA(2, rd, td, pd) * CUR_VAL_BETA(2, rd, td, pd)
	PRINT*,"25",JCBN_BIGK_FUNCTION
END IF

END SUBROUTINE JCBN_BIGK_SUBROUTINE






!+401+###########################################################################!
!                                                                                !
!              Calc_Shift_BC_1D                                                  !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Shift_1D( Shift_Vector,                                       &
                             NUM_R_QUAD, NUM_T_QUAD, NUM_P_QUAD,                 &
                             NUM_R_ELEM, NUM_T_ELEM, NUM_P_ELEM, Dimn,           &
                             Sr, r_locs                                          )


REAL(KIND = idp), DIMENSION( 0:NUM_R_ELEM ),INTENT( OUT )          ::  Shift_Vector

INTEGER,               INTENT( IN )                                ::  NUM_R_QUAD,  &
                                                                       NUM_T_QUAD,  &
                                                                       NUM_P_QUAD,  &
                                                                       NUM_R_ELEM,  &
                                                                       NUM_T_ELEM,  &
                                                                       NUM_P_ELEM,  &
                                                                       Dimn

REAL(KIND = idp), DIMENSION( 1:NUM_R_QUAD*NUM_T_QUAD*NUM_P_QUAD,  &
                             0:NUM_R_ELEM-1,                      &
                             0:NUM_T_ELEM-1,                      &
                             0:NUM_P_ELEM-1,                      &
                             1:Dimn ), INTENT(IN)                  ::  Sr

REAL(KIND = idp), DIMENSION( 0:NUM_R_ELEM ),  INTENT(IN)          ::  r_locs


COMPLEX(KIND = idp)                                               ::  Basis_Funcs,  &
                                                                      Tmp_Psi,      &
                                                                      Tmp_AlphaPsi

INTEGER                                                           ::  Ord,          &
                                                                      Current_Location


INTEGER                                                           ::  i, j, l, m, d, re, reb

REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                       :: x_locs,     &
                                                                     ri_locs,    &
                                                                     wi

REAL(KIND = idp), DIMENSION(:,:), ALLOCATABLE                     :: rij_locs,    &
                                                                     PSI_10,     &
                                                                     Sr_New

REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                       :: AlphaPsi
REAL(KIND = idp)                                                  :: Psi
REAL(KIND = idp)                                                  :: Beta_Tmp
REAL(KIND = idp)                                                  :: Inner_Int

REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                       :: xlocP, yi, LagP
REAL(KIND = idp)                                                  :: x_tmp
Ord = 6






ALLOCATE( LagP(0:DEGREE) )
ALLOCATE( xlocP(0:DEGREE) )
ALLOCATE( yi(0:DEGREE) )

ALLOCATE( x_locs(1:Ord) )
ALLOCATE( ri_locs(1:Ord) )
ALLOCATE( wi(1:Ord)  )
ALLOCATE( AlphaPsi(1:Ord) )

ALLOCATE( rij_locs(1:Ord,1:Ord) )
ALLOCATE( PSI_10(1:Ord,1:Ord)  )
ALLOCATE( Sr_New(1:Ord,1:Ord)  )

CALL Initialize_LG_Quadrature( Ord, x_locs, wi )
CALL Initialize_LGL_Quadrature( DEGREE, xlocP, yi)

Shift_Vector(0) = 0.0_idp
DO re = 0,NUM_R_ELEM-1

   ! Calculate the r locations for the Outer Integral's Quadrature Points !
   ri_locs(:) = (r_locs(re+1)-r_locs(re))/2.0_idp * ( x_locs(:) + 1.0_idp ) + r_locs(re)




   ! Calculate the Alpha Psi values at each of the Outer Integral's Quadrature Points !
   DO i = 1,Ord

      x_tmp = Map_To_X_Space(r_locs(re), r_locs(re+1), ri_locs(i))
      LagP = Lagrange_Poly(x_tmp, DEGREE, xlocP)
      Tmp_AlphaPsi = 0.0_idp

      DO l = 0,L_LIMIT
         DO m = -l,l
            DO d = 0,DEGREE

               Current_Location = Matrix_Location( 2, l, m, re, d )
               Basis_Funcs = Spherical_Harmonic(l,m,0.0_idp,0.0_idp)*LagP(d)

               Tmp_AlphaPsi = Tmp_AlphaPsi                                          &
                            + Coefficient_Vector(Current_Location)                  &
                            * Basis_Funcs


             END DO ! d Loop
          END DO ! m Loop
       END DO ! l Loop

    

       AlphaPsi(i) = REAL( Tmp_AlphaPsi, KIND = idp )

   END DO ! i Loop





   DO i = 1,Ord

      ! Calculate the Quadrature Points for each of the Inner Integrals
      rij_locs(:,i) = (ri_locs(i) - R_Inner)/2.0_idp *(x_locs(:) + 1.0_idp) + R_Inner


      ! Calculate Psi^10 values at each of the Inner Quadrature Points
      DO j = 1,Ord
         DO reb = 0,NUM_R_ELEM - 1
            IF ( (rij_Locs(j,i) > r_locs(reb)) .AND. (rij_Locs(j,i) < rlocs(reb+1)) ) THEN

               x_tmp = Map_To_X_Space(r_locs(reb), r_locs(reb+1), rij_locs(j,i))
               LagP = Lagrange_Poly(x_tmp, DEGREE, xlocP)
               Tmp_Psi = 0.0_idp

               DO l = 0,L_LIMIT
                  DO m = -l,l
                     DO d = 0,DEGREE

                        Current_Location = Matrix_Location( 1, l, m, reb, d )
                        Basis_Funcs = Spherical_Harmonic(l,m,0.0_idp,0.0_idp)*LagP(d)

                        Tmp_Psi = Tmp_Psi                                                    &
                                + Coefficient_Vector(Current_Location+0)                     &
                                * Basis_Funcs


                     END DO ! d Loop
                  END DO ! m Loop
               END DO ! l Loop

               ! Interpolate Sr Value to Quadrature Point !
               Sr_New(j,i) = Sr(1,reb,0,0,1)

            END IF
         END DO ! reb Loop

         Psi = REAL( Tmp_Psi, KIND = idp)
         Psi_10(j,i) = Psi**10

      END DO

   END DO ! i Loop



   ! Do the Outer Integral
   Beta_Tmp = 0.0_idp
   DO i = 1,Ord


     ! Do the Inner Integrals
     Inner_Int = 0.0_idp
     DO j = 1,Ord

        Inner_Int = Inner_Int                                 &
                  + rij_locs(j,i)*rij_locs(j,i)*rij_locs(j,i) &
                  * PSI_10(j,i)                               &
                  * Sr_New(j,i)                               &
                  * wi(j)


!         WRITE(*,'(ES22.15,1X,ES22.15,1X,ES22.5)')rij_locs(j,i)*rij_locs(j,i)*rij_locs(j,i),PSI_10(j,i),Sr_New(j,i)
      END DO ! j Loop

!      WRITE(*,'(I3,1X,I3,1X,ES22.15)') re,i,Inner_Int*(ri_locs(i) - R_Inner)/2.0_idp


      Beta_Tmp = Beta_Tmp                                        &
               + AlphaPsi(i)                                     &
               / (ri_locs(i)*ri_locs(i))                         &
               * Inner_Int                                       &
               * ( ri_locs(i) - R_Inner)/2.0_idp                 &
               * wi(i)

   

   END DO ! i Loop






   IF ( re == 0 ) THEN

      Shift_Vector(re+1) = (3.0_idp/2.0_idp)                        &
                         * r_locs(re+1)                             &
                         * ( r_locs(re+1) - R_locs(re))/2.0_idp        &
                         * Beta_Tmp

   ELSE 

      Shift_Vector(re+1) = r_locs(re+1)/r_locs(re)                  &
                         * Shift_Vector(re)                         &
                         + (3.0_idp/2.0_idp)                        &
                         * r_locs(re+1)                             &
                         * ( r_locs(re+1) - r_locs(re))/2.0_idp        &
                         * Beta_Tmp

   END IF



END DO ! re loop





!IF ( .TRUE. ) THEN
!   PRINT*,"    r                 Shift_Vector_R "
!   DO re = 0,NUM_R_ELEM

!       WRITE(*,'(ES22.15,1X,ES22.15,1X,ES22.15)') r_locs(re), Shift_Vector(re), Shift_Solution(r_locs(re),r_locs,NUM_R_ELEM)
!
!   END DO 
!END IF

END SUBROUTINE Calc_Shift_1D





!+401+###########################################################################!
!                                                                                !
!              Calc_Shift_1D                                                  !
!                                                                                !
!################################################################################!
SUBROUTINE Write_Shift_1D( Shift_Vector, NUM_R_ELEM, R_LOCS, SelfSim_T )

REAL(KIND = idp), DIMENSION( 0:NUM_R_ELEM ),INTENT( IN )          ::  Shift_Vector

INTEGER,               INTENT( IN )                               ::  NUM_R_ELEM

REAL(KIND = idp), DIMENSION( 0:NUM_R_ELEM ),  INTENT(IN)          ::  r_locs
REAL(KIND = idp)                                                  ::  SelfSim_T


INTEGER                                                           ::  re

INTEGER                                                           ::  File_id

CHARACTER( LEN = 43 )                                             ::  FileName

CHARACTER( LEN = 20 )                                             ::  FileDir
CHARACTER( LEN = 13 )                                             ::  FilePre
CHARACTER( LEN = 4 )                                              ::  FileExt


FileDir = "Shift_Vector_Output/"
FilePre = "Shift_Vector_"
FileExt = ".out"

WRITE( FileName, '(A,A,F6.4,A)' ) FileDir,FilePre,SelfSim_T,FileExt


PRINT*,Filename

File_id = 42

OPEN(Unit = File_id, file = FileName)

DO re = 0,NUM_R_ELEM

   WRITE(File_id,'(ES22.15,1X,ES22.15)') r_locs(re), Shift_Vector(re)

END DO


CLOSE(Unit = File_id)

END SUBROUTINE Write_Shift_1D






!+401+###########################################################################!
!                                                                                !
!              Calc_Shift_1D                                                  !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Shift_BC_1D( Shift_Vector_BC,                                    &
                          NUM_R_QUAD, NUM_T_QUAD, NUM_P_QUAD,                 &
                          NUM_R_ELEM, NUM_T_ELEM, NUM_P_ELEM, Dimn,           &
                          Sr, r_locs                                          )


REAL(KIND = idp),      INTENT( OUT )                               ::  Shift_Vector_BC

INTEGER,               INTENT( IN )                                ::  NUM_R_QUAD,  &
                                                                       NUM_T_QUAD,  &
                                                                       NUM_P_QUAD,  &
                                                                       NUM_R_ELEM,  &
                                                                       NUM_T_ELEM,  &
                                                                       NUM_P_ELEM,  &
                                                                       Dimn

REAL(KIND = idp), DIMENSION( 1:NUM_R_QUAD*NUM_T_QUAD*NUM_P_QUAD,  &
                             0:NUM_R_ELEM-1,                      &
                             0:NUM_T_ELEM-1,                      &
                             0:NUM_P_ELEM-1,                      &
                             1:Dimn ), INTENT(IN)                  ::  Sr

REAL(KIND = idp), DIMENSION( 0:NUM_R_ELEM ),  INTENT(IN)          ::  r_locs


COMPLEX(KIND = idp)                                               ::  Basis_Funcs,  &
                                                                      Tmp_Psi,      &
                                                                      Tmp_AlphaPsi

INTEGER                                                           ::  Ord,          &
                                                                      Current_Location


INTEGER                                                           ::  i, j, l, m, d, re

REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                       :: x_locs,     &
                                                                     ri_locs,    &
                                                                     wi

REAL(KIND = idp), DIMENSION(:,:), ALLOCATABLE                     :: rij_locs,    &
                                                                     PSI_10,     &
                                                                     Sr_New

REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                       :: AlphaPsi
REAL(KIND = idp)                                                  :: Psi
REAL(KIND = idp)                                                  :: Beta_Tmp
REAL(KIND = idp)                                                  :: Inner_Int

REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                       :: xlocP, yi, LagP
REAL(KIND = idp)                                                  :: x_tmp
Ord = 100






ALLOCATE( LagP(0:DEGREE) )
ALLOCATE( xlocP(0:DEGREE) )
ALLOCATE( yi(0:DEGREE) )

ALLOCATE( x_locs(1:Ord) )
ALLOCATE( ri_locs(1:Ord) )
ALLOCATE( wi(1:Ord)  )
ALLOCATE( AlphaPsi(1:Ord) )

ALLOCATE( rij_locs(1:Ord,1:Ord) )
ALLOCATE( PSI_10(1:Ord,1:Ord)  )
ALLOCATE( Sr_New(1:Ord,1:Ord)  )

CALL Initialize_LG_Quadrature( Ord, x_locs, wi )
CALL Initialize_LGL_Quadrature( DEGREE, xlocP, yi)


! Calculate the r locations for the Outer Integral's Quadrature Points !
ri_locs(:) = (R_Outer-R_inner)/2.0_idp * ( x_locs(:) + 1.0_idp ) + R_Inner




! Calculate the Alpha Psi values at each of the Outer Integral's Quadrature Points !
DO i = 1,Ord
   DO re = 0,NUM_R_ELEM - 1
       IF ( (ri_locs(i) > r_locs(re)) .AND. (ri_locs(i) <= rlocs(re+1)) ) THEN

             x_tmp = Map_To_X_Space(r_locs(re), r_locs(re+1), ri_locs(i))
             LagP = Lagrange_Poly(x_tmp, DEGREE, xlocP)
             Tmp_AlphaPsi = 0.0_idp

             DO l = 0,L_LIMIT
                DO m = -l,l
                   DO d = 0,DEGREE

                        Current_Location = Matrix_Location( 2, l, m, re, d )
                        Basis_Funcs = Spherical_Harmonic(l,m,0.0_idp,0.0_idp)*LagP(d)

                        Tmp_AlphaPsi = Tmp_AlphaPsi                                          &
                                     + Coefficient_Vector(Current_Location)                  &
                                     * Basis_Funcs


                   END DO ! d Loop
                END DO ! m Loop
             END DO ! l Loop

       END IF
   END DO ! re Loop


   AlphaPsi(i) = REAL( Tmp_AlphaPsi, KIND = idp )

END DO ! i Loop


DO i = 1,Ord

   ! Calculate the Quadrature Points for each of the Inner Integrals
   rij_locs(:,i) = (ri_locs(i) - R_Inner)/2.0_idp *(x_locs(:) + 1.0_idp) + R_Inner


   ! Calculate Psi^10 values at each of the Inner Quadrature Points
   DO j = 1,Ord
      DO re = 0,NUM_R_ELEM - 1
         IF ( (rij_Locs(j,i) > r_locs(re)) .AND. (rij_Locs(j,i) < rlocs(re+1)) ) THEN

            x_tmp = Map_To_X_Space(r_locs(re), r_locs(re+1), rij_locs(j,i))
            LagP = Lagrange_Poly(x_tmp, DEGREE, xlocP)
            Tmp_Psi = 0.0_idp

            DO l = 0,L_LIMIT
               DO m = -l,l
                  DO d = 0,DEGREE

                     Current_Location = Matrix_Location( 1, l, m, re, d )
                     Basis_Funcs = Spherical_Harmonic(l,m,0.0_idp,0.0_idp)*LagP(d)

                     Tmp_Psi = Tmp_Psi                                                    &
                             + Coefficient_Vector(Current_Location+0)                     &
                             * Basis_Funcs


                  END DO ! d Loop
               END DO ! m Loop
            END DO ! l Loop

            ! Interpolate Sr Value to Quadrature Point !
            Sr_New(j,i) = Sr(1,re,0,0,1)

         END IF
      END DO ! reb Loop

      Psi = REAL( Tmp_Psi, KIND = idp)
      Psi_10(j,i) = Psi**10

   END DO

END DO ! i Loop



! Do the Outer Integral
Beta_Tmp = 0.0_idp
DO i = 1,Ord


   ! Do the Inner Integrals
   Inner_Int = 0.0_idp
   DO j = 1,Ord

      Inner_Int = Inner_Int                                 &
                + rij_locs(j,i)*rij_locs(j,i)*rij_locs(j,i) &
                * PSI_10(j,i)                               &
                * Sr_New(j,i)                               &
                * wi(j)

   END DO ! j Loop


   Beta_Tmp = Beta_Tmp                                       &
            + AlphaPsi(i)                                    & 
            / (ri_locs(i)*ri_locs(i)*ri_locs(i)*ri_locs(i))  &
            * Inner_Int                                      &
            * (ri_locs(i) - R_Inner)/2.0_idp                 &
            * wi(i)



END DO ! i Loop

Shift_Vector_BC = (3.0_idp/2.0_idp)                        &
                * 8.0_idp * pi * GR_Source_Scalar          &
                * R_OUTER                                  &
                * Beta_Tmp                                 &
                * (R_Outer - R_Inner)/2.0_idp





END SUBROUTINE Calc_Shift_BC_1D













END MODULE Jacobian_Internal_Functions_Module
