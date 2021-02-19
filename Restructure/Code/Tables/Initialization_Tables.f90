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
                    ONLY : idp

USE Poseidon_Numbers_Module, &
                    ONLY : pi

USE Poseidon_Parameters, &
            ONLY :  Domain_Dim,             &
                    Degree,                 &
                    L_Limit

USE Variables_MPI, &
            ONLY :  Num_Block_Theta_Rows,   &
                    Num_T_Elems_Per_Block,  &
                    Num_P_Elems_Per_Block,  &
                    myID_Shell

USE Variables_Quadrature, &
            ONLY :  Num_R_Quad_Points,      &
                    Num_T_Quad_Points,      &
                    Num_P_Quad_Points,      &
                    Int_R_Locations,        &
                    Int_T_Locations,        &
                    Int_P_Locations

USE Variables_Mesh, &
            ONLY :  tlocs,                  &
                    plocs

USE Variables_Derived, &
            ONLY :  LM_Length


USE Variables_Tables, &
            ONLY :  Ylm_Values,             &
                    Ylm_dt_Values,          &
                    Ylm_dp_Values,          &
                    Ylm_CC_Values,          &
                    Ylm_CC_dt_Values,       &
                    Ylm_CC_dp_Values,       &
                    Lagrange_Poly_Table,    &
                    LPT_LPT,              &
                    M_Values

USE Variables_Functions, &
            ONLY :  LM_Location

USE Functions_Quadrature, &
            ONLY :  Initialize_LGL_Quadrature_Locations

USE Functions_Mapping, &
            ONLY :  Map_From_X_Space

USE Functions_Math, &
            ONLY :  Lagrange_Poly,          &
                    Lagrange_Poly_Deriv,    &
                    Legendre_Poly,          &
                    Norm_Factor,            &
                    Factorial

USE Allocation_Tables, &
            ONLY :  Allocate_Tables


IMPLICIT NONE

CONTAINS



!+101+##########################################################################!
!                                                                               !
!                       Initialize_Tables                                       !
!                                                                               !
!###############################################################################!
SUBROUTINE Initialize_Tables()

CALL Allocate_Tables()

CALL Initialize_Ylm_Tables()
CALL Initialize_Lagrange_Poly_Tables( Degree, Num_R_Quad_Points )

END SUBROUTINE Initialize_Tables





!+201+##########################################################################!
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

COMPLEX(KIND = idp), DIMENSION(0:LM_LENGTH-1)                      ::  Sqrt_Term
REAL(KIND = idp), DIMENSION(1:NUM_T_QUAD_POINTS)                ::  SIN_VAL, TAN_VAL
REAL(KIND = idp), DIMENSION(1:NUM_T_QUAD_POINTS)                ::  CSC_VAL, COT_VAL

COMPLEX(Kind = idp)                                             ::  Tmp_Value_A
COMPLEX(Kind = idp)                                             ::  Tmp_Value_B

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


IF ( DOMAIN_DIM == 1 ) THEN
     M_VALUES = 0
ELSE IF ( DOMAIN_DIM == 2 ) THEN
     M_VALUES = 0
ELSE IF ( DOMAIN_DIM == 3 ) THEN
     M_VALUES = (/(l,l=0,L_LIMIT,1)/)
END IF


Sqrt_Term(0) = 0.0_idp
DO l = 1,L_LIMIT

    REAL_L = REAL(l, idp)

    DO m = -M_VALUES(l),M_VALUES(l)

        Tmp_Value_A = COMPLEX( (2.0_idp * REAL_L + 1.0_idp)/(2.0_idp*Real_L - 1.0_idp),0.0_idp )
        Tmp_Value_B = COMPLEX( (l-m)*(l+m),0.0_idp )

        Sqrt_Term(LM_Location(l,m)) = zsqrt( Tmp_Value_A)*zsqrt(Tmp_Value_B)

    END DO ! m Loop
END DO ! l Loop




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


        END DO ! m Loop
        END DO ! l Loop
        END DO ! td Loop
        END DO ! pd Loop
    END DO ! te Loop
END DO ! pe Loop




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

!                PRINT*,t_locations(td),P_Locations(pd),lm_loc, Ylm_dt_Values(lm_loc,tpd_loc,te,pe)

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



            END DO ! m Loop
        END DO ! l Loop
        END DO ! td Loop
        END DO ! pd Loop
    END DO ! te Loop
END DO ! pe Loop



DEALLOCATE( Ylm_Table)








END SUBROUTINE Initialize_Ylm_Tables






!+202+##########################################################################!
!                                                                               !
!                  Initialize_Lagrange_Poly_Tables                              !
!                                                                               !
!###############################################################################!
SUBROUTINE Initialize_Lagrange_Poly_Tables( Ord, Num_Quad_Points )

INTEGER, INTENT(IN)                         ::  Ord
INTEGER, INTENT(IN)                         ::  Num_Quad_Points


REAL(KIND = idp), DIMENSION(0:Ord)          ::  Local_Locations
REAL(KIND = idp), DIMENSION(0:Ord)          ::  Lagrange_Poly_Values
REAL(KIND = idp), DIMENSION(0:Ord)          ::  Lagrange_DRV_Values
REAL(KIND = idp), DIMENSION(0:Ord)          ::  Lagrange_DDRV_Values

INTEGER                                     ::  Eval_Point
INTEGER                                     ::  rd, d, dp,dd


Lagrange_Poly_Table = 0.0_idp


Local_Locations = Initialize_LGL_Quadrature_Locations(Ord)

DO Eval_Point = 1,Num_Quad_Points

    Lagrange_Poly_Values = Lagrange_Poly(INT_R_Locations(Eval_Point), Ord, Local_Locations)
    Lagrange_DRV_Values  = Lagrange_Poly_Deriv(INT_R_Locations(Eval_Point), Ord, Local_Locations)

    Lagrange_Poly_Table(:, Eval_Point, 0) = Lagrange_Poly_Values
    Lagrange_Poly_Table(:, Eval_Point, 1) = Lagrange_DRV_Values

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







END MODULE Initialization_Tables
