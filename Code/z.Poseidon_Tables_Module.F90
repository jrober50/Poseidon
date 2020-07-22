   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_Tables_Module                                                       !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   Initialize_Ylm_Tables                                               !##!
!##!    +102+   Initialize_Lagrange_Poly_Tables                                     !##!
!##!                                                                                !##!
!##!    +201+   Lagrange_Poly                                                       !##!
!##!    +202+   Lagrange_Poly_Deriv                                                 !##!
!##!    +203+   Lagrange_Second_Deriv                                               !##!
!##!                                                                                !##!
!##!    +301+   Legendre_Poly                                                       !##!
!##!    +302+   Norm_Factor                                                         !##!
!##!    +303+   Factorial                                                           !##!
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
            ONLY :  idp, pi


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

USE Poseidon_Variables_Module, &
            ONLY :  NUM_R_ELEMENTS,         &
                    NUM_T_ELEMENTS,         &
                    NUM_P_ELEMENTS,         &
                    NUM_TP_QUAD_POINTS,     &
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
                    Ylm_dt_Values,          &
                    Ylm_dp_Values,          &
                    Ylm_CC_Values,          &
                    Ylm_CC_dt_Values,       &
                    Ylm_CC_dp_Values,       &
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

USE Poseidon_Quadrature_Module, &
            ONLY :  Initialize_LGL_Quadrature_Locations

USE Poseidon_Mapping_Functions_Module, &
            ONLY :  Map_From_X_Space

USE Poseidon_Math_Functions_Module, &
            ONLY :  Lagrange_Poly,          &
                    Lagrange_Poly_Deriv,    &
                    Lagrange_Second_Deriv,  &
                    Legendre_Poly,          &
                    Norm_Factor,            &
                    Factorial


IMPLICIT NONE

CONTAINS






!+101+##########################################################################!
!                                                                               !
!                  Initialize_Ylm_Tables                                        !
!                                                                               !
!###############################################################################!
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

INTEGER                                                         ::  lm_loc, tpd_loc

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



                        
                        tpd_loc = (td-1)*NUM_P_QUAD_POINTS + pd
                        lm_loc = LM_Location(l,m)
                        Norm_Storage = Norm_Factor(l,m)
                        Legendre_Poly_Value = Legendre_Poly(l, m, 1, T_Locations(td))

                        Ylm_Values(lm_loc, tpd_loc, te, pe) = Ylm_Table(m,l,td, pd, te, pe)

                        Ylm_dt_Values(lm_loc, tpd_loc, te, pe) = REAL_L*COT_VAL(td)                     &
                                                                * Ylm_Table(m,l,td, pd, te, pe)         &
                                                              - SQRT_TERM(lm_loc) * CSC_VAL(td)         &
                                                                * Ylm_Table(m,l-1,td, pd, te, pe)

                        Ylm_dp_Values(lm_loc, tpd_loc, te, pe) = CMPLX(0,m,idp)                         &
                                                                * Ylm_Table(m,l,td, pd, te, pe)



                        Ylm_CC_Values( tpd_loc, lm_loc, te, pe) = M_POWER_TABLE(m) * Ylm_Table(-m,l,td, pd, te, pe)

                        Ylm_CC_dp_Values( tpd_loc, lm_loc, te, pe) = CMPLX(0,-m,idp)                      &
                                                                 * Ylm_CC_Values(tpd_loc, lm_loc, te, pe)

                        
                        Ylm_CC_dt_Values( tpd_loc, lm_loc, te, pe) = REAL_L*COT_VAL(td)                     &
                                                                   * Ylm_CC_Values( tpd_loc, lm_loc, te, pe)  &
                                                                 - SQRT_TERM(lm_loc) * CSC_VAL(td)          &
                                                                   * M_POWER_TABLE(m)                       &  ! Complex Conjugate
                                                                   * Ylm_Table(-m, l-1, td, pd, te, pe)     ! of Y^{l-1}_{m}

                    END DO

                END DO

            END DO

        END DO

    END DO

END DO



DEALLOCATE( Ylm_Table)








END SUBROUTINE Initialize_Ylm_Tables






!+102+##########################################################################!
!                                                                               !
!                  Initialize_Lagrange_Poly_Tables                              !
!                                                                               !
!###############################################################################!
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
           
                LPT_LPT(rd, d, dp, 0:1, dd) = Lagrange_Poly_Table(d, rd, 0:1)       &
                                            * Lagrange_Poly_Table(dp, rd, dd)
            END DO
        END DO
    END DO
END DO





END SUBROUTINE Initialize_Lagrange_Poly_Tables







END MODULE Poseidon_Tables_Module
