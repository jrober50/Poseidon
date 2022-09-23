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
            ONLY :  Ylm_Values,             &
                    Ylm_dt_Values,          &
                    Ylm_dp_Values,          &
                    Ylm_CC_Values,          &
                    Ylm_CC_dt_Values,       &
                    Ylm_CC_dp_Values,       &
                    Lagrange_Poly_Table,    &
                    LPT_LPT,                &
                    M_Values,               &
                    Ylm_Sqrt_Table,         &
                    Ylm_Norm_Table,         &
                    rBT_NormedLegendre,     &
                    Ylm_Elem_Values,        &
                    Ylm_Elem_dt_Values,     &
                    Ylm_Elem_dp_Values,     &
                    Ylm_Elem_CC_Values,     &
                    Level_dx,               &
                    Level_Ratios,            &
                    LagPoly_MultiLayer_Table,&
                    LagPoly_Num_Tables

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
            ONLY :  Map_To_lm

USE Allocation_Tables, &
            ONLY :  Allocate_Tables

USE Flags_Initialization_Module, &
            ONLY :  lPF_Init_Tables_Flags,        &
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
LagPoly_Num_Tables = 2**(AMReX_Num_Levels+1) - 1  ! Sum of power of 2



#ifdef POSEIDON_AMREX_FLAG

    CALL Allocate_Tables()
    CALL Initialize_Lagrange_Poly_Tables( Degree, Num_R_Quad_Points )
    CALL Initialize_Ylm_Norm_Tables_AMReX()
    CALL Initialize_Level_Tables()

#else

    CALL Allocate_Tables()
    CALL Initialize_Lagrange_Poly_Tables( Degree, Num_R_Quad_Points )
    CALL Initialize_Ylm_Tables()
#endif


lPF_Init_Tables_Flags(iPF_Init_Tables_Init) = .TRUE.

END SUBROUTINE Initialize_Tables





!+201+##########################################################################!
!                                                                               !
!                  Initialize_Ylm_Tables                                        !
!                                                                               !
!###############################################################################!
SUBROUTINE Initialize_Ylm_Tables


INTEGER                                                         ::  l, m, te, pe, td, pd


REAL(idp), DIMENSION(1:NUM_T_QUAD_POINTS)                ::  T_Locations
REAL(idp), DIMENSION(1:NUM_P_QUAD_POINTS)                ::  P_Locations

REAL(idp), DIMENSION(1:1)                                ::  Legendre_Poly_Value


REAL(idp)                                                ::  Norm_Storage
REAL(idp)                                                ::  REAL_L

REAL(idp), DIMENSION(-L_LIMIT:L_LIMIT)                   :: M_POWER_TABLE

INTEGER                                                         ::  lm_loc, tpd_loc

COMPLEX(idp), ALLOCATABLE, DIMENSION(:,:,:,:,:,:)        ::  Ylm_Table

COMPLEX(idp), DIMENSION(1:LM_LENGTH)                     ::  Sqrt_Term
REAL(idp), DIMENSION(1:NUM_T_QUAD_POINTS)                ::  CSC_VAL, COT_VAL

COMPLEX(idp)                                             ::  Tmp_Value_A
COMPLEX(idp)                                             ::  Tmp_Value_B

ALLOCATE( Ylm_Table(-L_LIMIT:L_LIMIT, -1:L_LIMIT,                           &
                    1:NUM_T_QUAD_POINTS, 1:NUM_P_QUAD_POINTS,               &
                    0:NUM_T_Elements-1, 0:NUM_P_Elements-1)   )





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

        Tmp_Value_A = CMPLX( (2.0_idp * REAL_L + 1.0_idp)/(2.0_idp*Real_L - 1.0_idp),0.0_idp,idp )
        Tmp_Value_B = CMPLX( (l-m)*(l+m),0.0_idp, idp )

        Sqrt_Term(Map_To_lm(l,m)) = sqrt( Tmp_Value_A)*sqrt(Tmp_Value_B)

    END DO ! m Loop
END DO ! l Loop




Ylm_Table = 0.0_idp
DO pe = 0,Num_P_Elements-1

     !                                                          !
    !!    Map Phi Locations from [-1,1] space to real space.    !!
     !                                                          !
    P_Locations = Map_From_X_Space(plocs(pe), plocs(pe + 1), INT_P_LOCATIONS)
    

    DO te = 0,Num_T_Elements-1

         !                                                          !
        !!   Map Theta Locations from [-1,1] space to real space.   !!
         !                                                          !
        T_Locations = Map_From_X_Space(tlocs(te), tlocs(te + 1), INT_T_LOCATIONS)


        

        DO pd = 1,NUM_P_QUAD_POINTS
        DO td = 1,NUM_T_QUAD_POINTS
        DO l = 0,L_LIMIT
        DO m = -M_VALUES(l),M_VALUES(l)

            Norm_Storage = Norm_Factor(l,m)
            Legendre_Poly_Value = Legendre_Poly(l, m, 1, T_Locations(td))


            Ylm_Table(m, l, td, pd, te, pe) = Norm_Storage                            &   ! Normalization Factor
                                            * Legendre_Poly_Value(1)                  &   ! Legendre Polynomial
                                            * CDEXP(CMPLX(0.0_idp,m*P_Locations(pd),idp))   ! exp(im phi)


!            PRINT*,l,m,P_Locations(pd),CDEXP(CMPLX(0.0_idp,m*P_Locations(pd),idp))

        END DO ! m Loop
        END DO ! l Loop
        END DO ! td Loop
        END DO ! pd Loop
    END DO ! te Loop
END DO ! pe Loop





DO pe = 0,Num_P_Elements-1

     !                                                          !
    !!    Map Phi Locations from [-1,1] space to real space.    !!
     !                                                          !
    P_Locations = Map_From_X_Space(plocs(pe), plocs(pe + 1), INT_P_LOCATIONS)

    DO te = 0,Num_T_Elements-1

         !                                                          !
        !!   Map Theta Locations from [-1,1] space to real space.   !!
         !                                                          !
        T_Locations = Map_From_X_Space(tlocs(te), tlocs(te + 1), INT_T_LOCATIONS)

        CSC_VAL(:) = 1.0_idp/DSIN(T_Locations(:))
        COT_VAL(:) = 1.0_idp/DTAN(T_Locations(:))


        DO pd = 1,NUM_P_QUAD_POINTS
        DO td = 1,NUM_T_QUAD_POINTS
        DO l = 0,L_LIMIT
        DO m = -M_VALUES(l),M_VALUES(l)

            REAL_L = REAL(l, idp)

            tpd_loc = (td-1)*NUM_P_QUAD_POINTS + pd
            lm_loc = Map_To_lm(l,m)
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


REAL(idp), DIMENSION(0:Ord)          ::  Lagrange_Poly_Values
REAL(idp), DIMENSION(0:Ord)          ::  Lagrange_DRV_Values

INTEGER                                     ::  Eval_Point
INTEGER                                     ::  rd, d, dp,dd






#ifdef POSEIDON_AMREX_FLAG

INTEGER                                     ::  lvl, elem
INTEGER                                     ::  iNLE, Cur_Table
REAL(idp)                                   ::  ra, rb, wl
REAL(idp), DIMENSION(1:Num_R_Quad_Points)   ::  Local_R



DO lvl = 0,AMReX_Num_Levels-1
iNLE = 2**lvl-1
DO elem = 0,iNLE

    Cur_Table = iNLE+elem
    wl = 2.0_idp/2.0_idp**lvl
    ra = REAL(     -1 + wl*elem, Kind = idp)
    rb = REAL( -1 + wl*(elem+1), Kind = idp)


    Local_R(:) = Map_From_X_Space(ra, rb, Int_R_Locations(:))


    DO Eval_Point = 1,Num_R_Quad_Points
        

        Lagrange_Poly_Values = Lagrange_Poly(Local_R(Eval_Point), Ord, FEM_Node_xlocs)
        Lagrange_DRV_Values  = Lagrange_Poly_Deriv(Local_R(Eval_Point), Ord, FEM_Node_xlocs)

        LagPoly_MultiLayer_Table( :, Eval_Point, 0, Cur_Table) = Lagrange_Poly_Values
        LagPoly_MultiLayer_Table( :, Eval_Point, 1, Cur_Table) = Lagrange_DRV_Values


    END DO

END DO
END DO

#endif


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














!+203+##########################################################################!
!                                                                               !
!             Initialize_Ylm_Norm_Tables_AMReX                                  !
!                                                                               !
!###############################################################################!
SUBROUTINE Initialize_Ylm_Norm_Tables_AMReX


INTEGER                                                         ::  l, m

! Ylm = Norm_Factor(l,m) * Legendre_Poly(l,m,theta) * exp(im, phi)
!
! Create tables for
!       Norm_Factor
!       Sqrt_Factor (Like Norm_Factor but used for derivatives )
!

DO l = 0,L_LIMIT
DO m = -l,l
    Ylm_Norm_Table( m,l ) = Norm_Factor(l,m)
    Ylm_Sqrt_Table( m,l ) = CMPLX( Sqrt_Factor(l,m), 0.0_idp, idp )
END DO ! m Loop
END DO ! l Loop



END SUBROUTINE Initialize_Ylm_Norm_Tables_AMReX






!+203+##########################################################################!
!                                                                               !
!          Initialize_Normed_Legendre_Tables_on_Level                           !
!                                                                               !
!###############################################################################!
SUBROUTINE Initialize_Normed_Legendre_Tables_on_Level( iEU, iEL, Level )

INTEGER,   DIMENSION(3), INTENT(IN)                 ::  iEU
INTEGER,   DIMENSION(3), INTENT(IN)                 ::  iEL
INTEGER,                 INTENT(IN)                 ::  Level


INTEGER                                             ::  te, l, m
INTEGER                                             ::  lm, teb

REAL(idp), DIMENSION(1:Num_T_Quad_Points)           ::  Cur_T_Locs

REAL(idp), DIMENSION(-L_Limit:L_Limit,              &
                     0:L_Limit,                     &
                     1:Num_T_Quad_Points,           &
                     0:AMReX_Max_Grid_Size(2)-1 )   ::  LegPoly_Table



! Ylm = Norm_Factor(l,m) * Legendre_Poly(l,m,theta) * exp(im,phi)
!     = rBT_NormedLegendre(l,m,theta) * exp(-im,phi)
!
! Create table for
!       rBT_NormedLegendre(l,m,theta)
!
!


!PRINT*,"A2",ALLOCATED(rBT_NormedLegendre),ALLOCATED(Int_T_Locations)
rBT_NormedLegendre = 0.0_idp
DO te = iEL(2), iEU(2)
    teb = te-iEL(2)
    Cur_T_locs = Level_dx(Level,2)/2.0_idp * (Int_T_Locations(:) + 1.0_idp + 2.0_idp*te )

    DO l = 0,L_LIMIT
    DO m = -l,l

        LegPoly_Table(m, l, :, teb) = Legendre_Poly(l,m,Num_T_Quad_Points,Cur_T_Locs)

    END DO ! m Loop
    END DO ! l Loop
END DO ! te Loop



!PRINT*,"B2"
DO te = iEL(2), iEU(2)
DO l = 0,L_LIMIT
DO m = -l,l

    lm = Map_To_lm(l,m)
    teb = te-iEL(2)
    rBT_NormedLegendre(m,l,:,teb) = Ylm_Norm_Table( m, l )      &
                                  * LegPoly_Table(m,l,:,teb)
!    PRINT*,te,l,m
!    PRINT*,Ylm_Norm_Table( m, l )
!    PRINT*,LegPoly_Table(m,l,:,teb)
!    PRINT*,"+++++++++++++++"
END DO ! m Loop
END DO ! l Loop
END DO ! te Loop
!PRINT*,"C2"



END SUBROUTINE Initialize_Normed_Legendre_Tables_on_Level










!+203+##########################################################################!
!                                                                               !
!          Initialize_Ylm_Tables_on_Elem                                  !
!                                                                               !
!###############################################################################!
SUBROUTINE Initialize_Ylm_Tables_on_Elem( iTE, iPE, iEL, Level )

INTEGER,    INTENT(IN)                              ::  iTE
INTEGER,    INTENT(IN)                              ::  iPE
INTEGER,    INTENT(IN)                              ::  iEL(3)
INTEGER,    INTENT(IN)                              ::  Level

INTEGER                                             ::  l, m, lm
INTEGER                                             ::  td, pd, tpd
INTEGER                                             ::  teb

REAL(idp), DIMENSION(1:Num_T_Quad_Points)           ::  Cot_Val
REAL(idp), DIMENSION(1:Num_T_Quad_Points)           ::  Csc_Val

REAL(idp), DIMENSION(1:Num_T_Quad_Points)           ::  tlocs
REAL(idp), DIMENSION(1:Num_P_Quad_Points)           ::  plocs

REAL(idp), DIMENSION(-L_LIMIT:L_LIMIT)       :: M_POWER_TABLE

!PRINT*,"In Init_Ylm_Tables",iTE,iPE

teb = iTE - iEL(2)

tlocs = Level_dx(Level,2)/2.0_idp * (Int_T_Locations(:) + 1.0_idp + 2.0_idp*iTE )
plocs = Level_dx(Level,3)/2.0_idp * (Int_P_Locations(:) + 1.0_idp + 2.0_idp*iPE )



Csc_Val = 1.0_idp/DSIN(tlocs(:))
Cot_Val = 1.0_idp/DTAN(tlocs(:))



M_POWER_TABLE(0) = 1.0_idp
IF ( L_LIMIT > 0 ) THEN
    DO m = 1, L_LIMIT

        M_POWER_TABLE(m) = -1.0_idp*M_POWER_TABLE(m-1)
        M_POWER_TABLE(-m) = -1.0_idp*M_POWER_TABLE(m-1)
    END DO
END IF



DO td = 1,Num_T_Quad_Points
DO pd = 1,Num_P_Quad_Points
DO l  = 0,L_Limit
DO m  = -l,l

    tpd = Map_To_tpd(td,pd)
    lm = Map_To_lm(l,m)


    Ylm_Elem_Values( lm, tpd )      = rBT_NormedLegendre(m,l,td,teb)           &
                                    * CDEXP(CMPLX(0.0_idp,m*plocs(pd),idp))


!    PRINT*,lm,tpd, rBT_NormedLegendre(m,l,td,teb),CDEXP(CMPLX(0.0_idp,m*plocs(pd),idp))

    Ylm_Elem_dt_Values( lm, tpd )   = REAL( l, idp ) * COT_VAL(td)              &
                                      * rBT_NormedLegendre(m,l,td,teb)         &
                                      * CDEXP(CMPLX(0.0_idp,m*plocs(pd),idp))   &
                                    - Ylm_Sqrt_Table(m,l) * CSC_VAL(td)         &
                                      * rBT_NormedLegendre(m,l-1,td,teb)       &
                                      * CDEXP(CMPLX(0.0_idp,m*plocs(pd),idp))

    Ylm_Elem_dp_Values( lm, tpd )   = CMPLX(0,m,idp)                            &
                                    * rBT_NormedLegendre(m,l,td,teb)           &
                                    * CDEXP(CMPLX(0.0_idp,m*plocs(pd),idp))

    
    Ylm_Elem_CC_Values( tpd, lm )   = M_Power_Table(m)                          &
                                    * rBT_NormedLegendre(-m,l,td,teb)          &
                                    * CDEXP(CMPLX(0.0_idp,-m*plocs(pd),idp))



END DO ! m Loop
END DO ! l Loop
END DO ! pd Loop
END DO ! td Loop

END SUBROUTINE Initialize_Ylm_Tables_on_Elem








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
