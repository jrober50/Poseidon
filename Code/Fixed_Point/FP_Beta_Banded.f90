  !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE FP_Beta_Banded                                                               !##!
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

USE Poseidon_Parameters, &
            ONLY :  DEGREE,                     &
                    L_LIMIT,                    &
                    Verbose_Flag

USE Poseidon_IO_Parameters, &
            ONLY :  Poseidon_LinSys_Dir

    
USE Variables_Quadrature, &
            ONLY :  Num_R_Quad_Points,          &
                    Num_T_Quad_Points,          &
                    Num_P_Quad_Points,          &
                    Num_TP_Quad_Points,         &
                    Int_R_Locations,            &
                    Int_T_Locations,            &
                    Int_P_Locations,            &
                    Int_R_Weights,              &
                    Int_T_Weights,              &
                    Int_P_Weights


USE Variables_Mesh, &
            ONLY :  Num_R_Elements,             &
                    Num_T_Elements,             &
                    Num_P_Elements,             &
                    rlocs,                      &
                    tlocs,                      &
                    plocs

USE Variables_Derived, &
            ONLY :  Num_R_Nodes,                &
                    ULM_Length,                 &
                    LM_Length,                  &
                    Var_Dim,                    &
                    Beta_Prob_Dim

USE Variables_Tables, &
            ONLY :  LPT_LPT,                    &
                    Ylm_Values,             &
                    Ylm_dt_Values,          &
                    Ylm_dp_Values,          &
                    Ylm_CC_Values,          &
                    Ylm_CC_dt_Values,       &
                    Ylm_CC_dp_Values

USE Variables_IO, &
            ONLY :  File_Suffix

USE Functions_Quadrature, &
            ONLY :  Initialize_LG_Quadrature

USE Variables_FP, &
            ONLY :  Matrix_Format,              &
                    Beta_Diagonals,             &
                    Beta_Bandwidth,             &
                    Beta_IPIV,                  &
                    Beta_MVL_Banded,            &
                    Beta_MVL_Diagonal,          &
                    CFA_EQ_Flags,               &
                    CFA_EQ_Map,                 &
                    CFA_Mat_Map,                &
                    First_Column_Beta_Storage,  &
                    Last_Column_Beta_Storage,   &
                    Beta_Factorized_Flag

USE Poseidon_IO_Module, &
            ONLY :  Open_New_File,              &
                    CLOCK_IN

USE IO_FP_Linear_System, &
            ONLY :  Output_Laplace,             &
                    Output_Laplace_Beta,        &
                    Output_Source_Beta

USE Variables_BC, &
            ONLY :  INNER_CFA_BC_VALUES,        &
                    OUTER_CFA_BC_VALUES,        &
                    INNER_CFA_BC_TYPE,          &
                    OUTER_CFA_BC_TYPE

USE FP_Functions_Mapping, &
            ONLY :  FP_Beta_Array_Map,          &
                    FP_FEM_Node_Map,            &
                    FP_LM_Map


USE MPI

IMPLICIT NONE



CONTAINS






!+101+###########################################################################!
!                                                                                !
!           Initialize_Laplace_Matrices_Full                                     !
!                                                                                !
!################################################################################!
SUBROUTINE Initialize_Beta_MVL_Banded()



INTEGER                                                 ::  l, m, lp, mp, lm, lm_loc, lpmp_loc
INTEGER                                                 ::  re, te, pe, rep
INTEGER                                                 ::  rd, td, pd, tpd
INTEGER                                                 ::  d, dp
INTEGER                                                 ::  i, j, ui, uj
INTEGER                                                 ::  row, col

REAL(KIND = idp)                                        ::  deltar, TODR
REAL(KIND = idp)                                        ::  L_Lp1

INTEGER                                                 ::  Mat_Loc

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             ::  CUR_R_LOCS
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             ::  CUR_T_LOCS
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             ::  R_SQUARE

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: SIN_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: COTAN_VAL


COMPLEX(idp), DIMENSION(0:DEGREE)                       ::  Reusable_Values
COMPLEX(idp), DIMENSION(0:DEGREE)                       ::  Reusable_Values_b



REAL(KIND = idp)                                        ::  deltar_overtwo,     &
                                                            deltat_overtwo,     &
                                                            deltap_overtwo

REAL(KIND = idp), ALLOCATABLE, DIMENSION( :,:,: )       :: RR_Factor
REAL(KIND = idp), ALLOCATABLE, DIMENSION( :,:,: )       :: DRR_Factor
REAL(KIND = idp), ALLOCATABLE, DIMENSION( :,:,: )       :: RDR_Factor
REAL(KIND = idp), ALLOCATABLE, DIMENSION( :,:,: )       :: DRDR_Factor


COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION( :,:,: )    :: TP_TP_Integrals

INTEGER                                                 :: INFO

REAL(idp), DIMENSION(1:3)                               :: Timer
REAL(idp), DIMENSION(1:3)                               :: Timer_Tots

Beta_MVL_Banded = 0.0_idp

ALLOCATE( CUR_R_LOCS(1:Num_R_Quad_Points) )
ALLOCATE( CUR_T_LOCS(1:Num_T_Quad_Points) )
ALLOCATE( R_SQUARE(1:Num_R_Quad_Points)   )

ALLOCATE( SIN_SQUARE( 1:NUM_TP_QUAD_POINTS ) )
ALLOCATE( COTAN_VAL( 1:NUM_TP_QUAD_POINTS ) )

ALLOCATE( RR_Factor( 1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE )    )
ALLOCATE( RDR_Factor( 1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE )    )
ALLOCATE( DRR_Factor( 1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE )    )
ALLOCATE( DRDR_Factor( 1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE )    )

ALLOCATE( TP_TP_Integrals( 1:LM_Length, 1:LM_Length, 1:16) )


IF ( Verbose_Flag ) THEN
    PRINT*,"-Initializing Modified Vector Laplacian Matrix.  Format: Banded. "
END IF


DO re = 0,NUM_R_ELEMENTS-1
DO d = 0,DEGREE
DO ui = 1,3
DO l = 0,L_LIMIT


    L_Lp1 = REAL( l*(l+1), idp )

    deltar        = rlocs(re+1) - rlocs(re)
    TODR          = 2.0_idp/deltar
    CUR_R_LOCS(:) = (deltar/2.0_idp) * (Int_R_Locations(:)+1.0_idp) + rlocs(re)
    R_SQUARE      = CUR_R_LOCS**2



    Reusable_Values = 0.0_idp
    DO dp = 0,DEGREE
        DO rd = 1,Num_R_Quad_Points

            Reusable_Values(dp) = Reusable_Values(dp) &
                                - R_SQUARE(rd) * LPT_LPT(rd,d,dp,1,1)   &
                                    * TODR * Int_R_Weights(rd)            &
                                + L_Lp1 * LPT_LPT(rd,d,dp,0,0)          &
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
END DO ! re Loop




! At this point Laplace Matrix_Beta should contain the discretized Laplace operator for the three
! Beta components.  To test set up three Poisson/Laplace problems and assign one to each
! Beta component.  There is no coupling yet, so this should solve the three independently.

!Call Output_Laplace_Beta(Beta_MVL_Banded,Beta_Prob_Dim, Beta_Prob_Dim, "Z")


DO re = 0,Num_R_Elements-1
    deltar = rlocs(re+1) - rlocs(re)
    
    TODR = 2.0_idp/deltar
    CUR_R_LOCS(:) = (deltar/2.0_idp) * (Int_R_Locations(:)+1.0_idp) + rlocs(re)
    R_SQUARE = CUR_R_LOCS**2

    CALL Calc_RR_Values( R_Square, TODR, RR_Factor, DRR_Factor, RDR_Factor, DRDR_Factor )




    DO te = 0,Num_T_Elements-1

        DeltaT_OverTwo = (tlocs(te + 1) - tlocs(te))/2.0_idp
        CUR_T_LOCS(:) = DeltaT_OverTwo * (INT_T_LOCATIONS(:)+1.0_idp) + tlocs(te)




        DO pe = 0,Num_P_Elements-1



            DeltaP_OverTwo = (plocs(pe + 1) - plocs(pe))/2.0_idp


            CALL Calc_TP_Values( DeltaT_OverTwo, DeltaP_OverTwo, Cur_T_Locs, te, pe,                &
                                 Sin_Square, Cotan_Val,                                             &
                                 TP_TP_Integrals            )



            DO d = 0,DEGREE
                DO l = 0,L_LIMIT
                DO m = -l,l

                    CALL Calc_Beta1_Terms( re, d, l, m,                                       &
                                            RR_Factor, dRR_Factor, dRdR_Factor,                &
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
        END DO ! pe Loop
    END DO ! te Loop
END DO ! re Loop


!    PRINT*,"Work_Mat"
!    DO ui = 1,3
!    re = Num_R_Elements-1
!    DO d = 0,Degree
!    DO l = 1,LM_Length
!    DO uj = 1,3
!    rep = Num_R_Elements-1
!    DO dp = 0,Degree
!    DO lp = 1,LM_Length
!
!        i = FP_Beta_Array_Map(re,d,ui,l)
!        j = FP_Beta_Array_Map(rep,dp,uj,lp)
!
!        PRINT*,i,j,Beta_MVL_Banded(Beta_Bandwidth+i-j,j)
!    END DO
!    END DO
!    END DO
!    END DO
!    END DO
!    END DO



!Call Output_Laplace_Beta(Beta_MVL_Banded,Beta_Prob_Dim, Beta_Prob_Dim)
Beta_Factorized_Flag = .FALSE.


END SUBROUTINE Initialize_Beta_MVL_Banded

















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

REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ), INTENT(IN)  :: RR_Factor
REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ), INTENT(IN)  :: DRR_Factor
REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ), INTENT(IN)  :: DRDR_Factor

COMPLEX(idp), DIMENSION( 1:LM_Length, 1:LM_Length, 1:16),       INTENT(IN)  :: TP_TP_Integrals


REAL(idp),    DIMENSION( 1:Num_R_Quad_Points ),                 INTENT(IN)  :: Cur_R_Locs
REAL(idp),    DIMENSION( 1:Num_R_Quad_Points ),                 INTENT(IN)  :: R_Square


INTEGER                                                     :: d, rd, ui, uj
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

REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ), INTENT(IN)  :: RR_Factor
REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ), INTENT(IN)  :: DRR_Factor
REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ), INTENT(IN)  :: RdR_Factor
REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ), INTENT(IN)  :: DRDR_Factor

COMPLEX(idp), DIMENSION( 1:LM_Length, 1:LM_Length, 1:16),       INTENT(IN)  :: TP_TP_Integrals


REAL(idp),    DIMENSION( 1:Num_R_Quad_Points ),                 INTENT(IN)  :: Cur_R_Locs
REAL(idp),    DIMENSION( 1:Num_R_Quad_Points ),                 INTENT(IN)  :: R_Square


INTEGER                                                     :: d, rd, ui, uj
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

REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ), INTENT(IN)    :: RR_Factor
REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ), INTENT(IN)    :: DRR_Factor
REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ), INTENT(IN)    :: RDR_Factor
REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ), INTENT(IN)    :: DRDR_Factor

COMPLEX(idp), DIMENSION( 1:LM_Length, 1:LM_Length, 1:16),       INTENT(IN)    :: TP_TP_Integrals

REAL(idp),    DIMENSION( 1:Num_R_Quad_Points ),                 INTENT(IN)    :: Cur_R_Locs
REAL(idp),    DIMENSION( 1:Num_R_Quad_Points ),                 INTENT(IN)    :: R_Square


INTEGER                                                     :: d, rd, ui, uj
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








!+403+###########################################################################!
!                                                                                !
!                   CALC_RR_Values                                               !
!                                                                                !
!################################################################################!
SUBROUTINE CALC_RR_Values( R_Square, TODR,          &
                           RR_Factor, DRR_Factor,   &
                           RDR_Factor, DRDR_Factor  )

REAL(idp), DIMENSION(1:Num_R_Quad_Points),                      INTENT(IN)      :: R_Square
REAL(idp),                                                      INTENT(IN)      :: TODR

REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ), INTENT(INOUT)   :: RR_Factor
REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ), INTENT(INOUT)   :: DRR_Factor
REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ), INTENT(INOUT)   :: RDR_Factor
REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ), INTENT(INOUT)   :: DRDR_Factor


INTEGER                                                                         :: d, dp
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)                                     :: R_Int_Weights



ALLOCATE( R_Int_Weights( 1:NUM_R_QUAD_POINTS ) )




R_Int_Weights(:) = R_SQUARE(:) * INT_R_WEIGHTS(:)/TODR
DO d = 0,DEGREE
    DO dp = 0,DEGREE

        RR_Factor(:, d, dp)      = R_Int_Weights(:) * LPT_LPT(:, d, dp, 0, 0)
        DRR_Factor(:, d, dp)     = R_Int_Weights(:) * LPT_LPT(:, d, dp, 1, 0) * TODR
        RDR_Factor(:, d, dp)     = R_Int_Weights(:) * LPT_LPT(:, d, dp, 0, 1) * TODR
        DRDR_Factor(:, d, dp)    = R_Int_Weights(:) * LPT_LPT(:, d, dp, 1, 1) * TODR * TODR

    END DO
END DO



END SUBROUTINE CALC_RR_Values






!+403+###########################################################################!
!                                                                                !
!                   CALC_TP_Values                                       !
!                                                                                !
!################################################################################!
SUBROUTINE CALC_TP_Values( DeltaT_OverTwo, DeltaP_OverTwo, Cur_T_Locs, te, pe,                  &
                            Sin_Square, Cotan_Val,                                              &
                            TP_TP_Integrals         )

REAL(idp),                                                                    INTENT(IN)    :: DeltaT_OverTwo
REAL(idp),                                                                    INTENT(IN)    :: DeltaP_OverTwo

REAL(idp),    DIMENSION( 1:Num_T_Quad_Points ),                               INTENT(IN)    :: Cur_T_Locs

INTEGER,                                                                      INTENT(IN)    ::  te, pe

REAL(idp),    DIMENSION( 1:Num_TP_Quad_Points ),                              INTENT(INOUT) :: Sin_Square
REAL(idp),    DIMENSION( 1:Num_TP_Quad_Points ),                              INTENT(INOUT) :: Cotan_Val


COMPLEX(idp), DIMENSION( 1:LM_Length, 1:LM_Length, 1:16), INTENT(INOUT)                 :: TP_TP_Integrals


INTEGER                                                                                     :: td, pd, tpd
INTEGER                                                                                     :: lm_loc, lpmp_loc
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)                                                 :: TP_Int_Weights

COMPLEX(idp)                                                                                :: Test

ALLOCATE( TP_Int_Weights( 1:NUM_TP_QUAD_POINTS) )



DO td = 1,NUM_T_QUAD_POINTS
    DO pd = 1,NUM_P_QUAD_POINTS
        tpd = (td-1)*NUM_P_QUAD_POINTS + pd

        Sin_Square(tpd) = DSIN( Cur_T_Locs(td) )*DSIN( Cur_T_Locs(td) )
        Cotan_Val(tpd)  = 1.0_idp/DTAN( CUR_T_LOCS(td) )
!        TP_Int_Weights( (td-1)*NUM_P_QUAD_POINTS + pd ) = DSIN( Cur_T_Locs(td) )                &
!                                                        * DeltaT_OverTwo * INT_T_WEIGHTS(td)    &
!                                                        * DeltaP_OverTwo * INT_P_WEIGHTS(pd)
 
        TP_Int_Weights( (td-1)*NUM_P_QUAD_POINTS + pd ) = DSIN( Cur_T_Locs(td) )                &
                                                        * DeltaT_OverTwo * INT_T_WEIGHTS(td)    &
                                                        * INT_P_WEIGHTS(pd)

END DO
END DO




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
    TP_TP_Integrals( lm_loc, lpmp_loc, 2 ) = SUM( Ylm_Values( lm_loc, :, te, pe )               &
                                                     * Ylm_CC_dp_Values( :, lpmp_loc, te, pe)   &
                                                     * TP_Int_Weights(:)                        )


    ! Ylm * Y^lpmp * Cotan
    TP_TP_Integrals( lm_loc, lpmp_loc, 2 ) = SUM( Ylm_Values( lm_loc, :, te, pe )               &
                                                     * Ylm_CC_Values( :, lpmp_loc, te, pe)      &
                                                     * TP_Int_Weights(:)                        &
                                                     * Cotan_Val(:)                             )

    ! d Ylm/dt * Y^lpmp     SUM( dTP_TP_Factor(:,lm_loc,lpmp_loc))
    TP_TP_Integrals( lm_loc, lpmp_loc, 3 ) = SUM( Ylm_dt_Values( lm_loc, :, te, pe )            &
                                                    * Ylm_CC_Values( :, lpmp_loc, te, pe)       &
                                                    * TP_Int_Weights(:)                         )

    ! Ylm * d Y^lpmp/dp
    TP_TP_Integrals( lm_loc, lpmp_loc, 4 ) = SUM( Ylm_Values( lm_loc, :, te, pe )               &
                                                    * Ylm_CC_dt_Values( :, lpmp_loc, te, pe)    &
                                                    * TP_Int_Weights(:)                         )

    ! d Ylm/dp * Y^lpmp
    TP_TP_Integrals( lm_loc, lpmp_loc, 5 ) = SUM( Ylm_dp_Values( lm_loc, :, te, pe )           &
                                                    * Ylm_CC_Values( :, lpmp_loc, te, pe)       &
                                                    * TP_Int_Weights(:) )

    ! SUM( dTP_dTP_Factor(:,lm_loc,lpmp_loc) )
    TP_TP_Integrals( lm_loc, lpmp_loc, 6 ) = SUM( Ylm_dt_Values( lm_loc, :, te, pe )             &
                                                    * Ylm_CC_dt_Values( :, lpmp_loc, te, pe)    &
                                                    * TP_Int_Weights(:) )

    ! SUM( dTP_TP_Factor(:,lm_loc,lpmp_loc) * Cotan_Val(:) )
    TP_TP_Integrals( lm_loc, lpmp_loc, 7 ) = SUM( Ylm_dt_Values( lm_loc, :, te, pe )           &
                                                    * Ylm_CC_Values( :, lpmp_loc, te, pe)       &
                                                    * TP_Int_Weights(:)                         &
                                                    * Cotan_Val(:) )

    ! SUM( TP_TP_Factor(:,lm_loc,lpmp_loc) / Sin_Square(:) )
    TP_TP_Integrals( lm_loc, lpmp_loc, 8 ) = SUM( Ylm_Values( lm_loc, :, te, pe )              &
                                                    * Ylm_CC_Values( :, lpmp_loc, te, pe)       &
                                                    * TP_Int_Weights(:)                         &
                                                    / Sin_Square(:) )

    ! SUM( TP_TP_Factor(:,lm_loc,lpmp_loc) * (1-Cotan_Val(:)**2) )
    TP_TP_Integrals( lm_loc, lpmp_loc, 9 ) = SUM( Ylm_Values( lm_loc, :, te, pe )              &
                                                    * Ylm_CC_Values( :, lpmp_loc, te, pe)       &
                                                    * TP_Int_Weights(:)                         &
                                                    * (1-Cotan_Val(:)**2) )

    ! SUM( dTP_TdP_Factor(:,lm_loc,lpmp_loc) )
    TP_TP_Integrals( lm_loc, lpmp_loc, 10 ) = SUM( Ylm_dt_Values( lm_loc, :, te, pe )             &
                                                    * Ylm_CC_dp_Values( :, lpmp_loc, te, pe)    &
                                                    * TP_Int_Weights(:) )

    ! SUM( TdP_TP_Factor(:,lm_loc,lpmp_loc) * Cotan_Val(:)   )
    TP_TP_Integrals( lm_loc, lpmp_loc, 11 ) = SUM( Ylm_dp_Values( lm_loc, :, te, pe )           &
                                                    * Ylm_CC_Values( :, lpmp_loc, te, pe)       &
                                                    * TP_Int_Weights(:)                         &
                                                    * Cotan_Val(:)   )

    ! SUM( TdP_TP_Factor(:,lm_loc,lpmp_loc) / Sin_Square(:)     )
    TP_TP_Integrals( lm_loc, lpmp_loc, 12 ) = SUM( Ylm_dp_Values( lm_loc, :, te, pe )           &
                                                    * Ylm_CC_Values( :, lpmp_loc, te, pe)       &
                                                    * TP_Int_Weights(:)                         &
                                                    / Sin_Square(:)     )

    ! SUM( TdP_dTP_Factor(:,lm_loc,lpmp_loc)/ Sin_Square(:) )
    TP_TP_Integrals( lm_loc, lpmp_loc, 13 ) = SUM( Ylm_dp_Values( lm_loc, :, te, pe )             &
                                                    * Ylm_CC_dt_Values( :, lpmp_loc, te, pe)    &
                                                    * TP_Int_Weights(:)                         &
                                                    / Sin_Square(:) )

    ! SUM( TdP_TP_Factor(:,lm_loc,lpmp_loc) * Cotan_Val(:) / Sin_Square(:)     )
    TP_TP_Integrals( lm_loc, lpmp_loc, 14 ) = SUM( Ylm_dp_Values( lm_loc, :, te, pe )           &
                                                    * Ylm_CC_Values( :, lpmp_loc, te, pe)       &
                                                    * TP_Int_Weights(:)                         &
                                                    * Cotan_Val(:)                              &
                                                    / Sin_Square(:)     )

    ! SUM( TdP_TdP_Factor(:,lm_loc,lpmp_loc)/Sin_Square(:) )
    TP_TP_Integrals( lm_loc, lpmp_loc, 15 ) = SUM( Ylm_dp_Values( lm_loc, :, te, pe )             &
                                                    * Ylm_CC_dp_Values( :, lpmp_loc, te, pe)    &
                                                    * TP_Int_Weights(:)                         &
                                                    / Sin_Square(:) )

    ! SUM( dTP_TP_Factor(:,lm_loc,lpmp_loc) * Cotan_Val(:) )
    TP_TP_Integrals( lm_loc, lpmp_loc, 16 ) = SUM( Ylm_dt_Values( lm_loc, :, te, pe )           &
                                                    * Ylm_CC_Values( :, lpmp_loc, te, pe)       &
                                                    * TP_Int_Weights(:)                         &
                                                    * Cotan_Val(:) )


END DO
END DO









END SUBROUTINE CALC_TP_Values


!+501+###########################################################################!
!                                                                                !
!                   Factorize_Beta_Banded                                       !
!                                                                                !
!################################################################################!
SUBROUTINE Factorize_Beta_Banded()


INTEGER                                                 :: Col, Row
INTEGER                                                 :: lm

INTEGER                                                 ::  i, j
INTEGER                                                 ::  re, d, l, ui
INTEGER                                                 ::  rep, dp, lp, uj


INTEGER                                                 :: INFO
COMPLEX(idp),   DIMENSION(2*Beta_Prob_Dim)              :: WORK
REAL(idp),      DIMENSION(Beta_Prob_Dim)                :: RWORK
REAL(idp)                                               :: RCOND_One
REAL(idp)                                               :: NORM

REAL(idp),      DIMENSION(2)                            :: Timer



IF ( Verbose_Flag ) THEN
    PRINT*,"--In Factorize_Beta_Banded."
END IF

timer(1) = MPI_Wtime()

!   Dirichlet BCs modify the stiffness matrix so we modify it now.
!   But to apply the BCs we will need values from the original matrix,
!   so those values are stored before we modify the matrix.
!
CALL DIRICHLET_BC_Beta_Banded_Mat()

!PRINT*,"Work_Mat"
!DO ui = 1,3
!re = Num_R_Elements-1
!DO d = 0,Degree
!DO l = 1,LM_Length
!DO uj = 1,3
!rep = Num_R_Elements-1
!DO dp = 0,Degree
!DO lp = 1,LM_Length
!
!    i = FP_Beta_Array_Map(re,d,ui,l)
!    j = FP_Beta_Array_Map(rep,dp,uj,lp)
!
!    PRINT*,i,j,Beta_MVL_Banded(Beta_Bandwidth+i-j,j)
!END DO
!END DO
!END DO
!END DO
!END DO
!END DO



CALL Jacobi_PC_MVL_Banded()




NORM = 0.0_idp
DO i = 1,Beta_Prob_Dim

    NORM = MAX( NORM, ABS(SUM(Beta_MVL_Banded(:,i) ) ) )

END DO



CALL ZGBTRF( Beta_Prob_Dim,             &
             Beta_Prob_Dim,             &
             Beta_Diagonals,            &
             Beta_Diagonals,            &
             Beta_MVL_Banded,           &
             3*Beta_Diagonals+1,        &
             Beta_IPIV,                 &
             INFO                       )

IF (INFO .NE. 0) THEN
    print*,"ZGBTRF has failed with INFO = ",INFO
ELSE

    Beta_Factorized_Flag = .TRUE.

END IF



CALL ZGBCON( '1',                   &
             Beta_Prob_Dim,         &
             Beta_Diagonals,        &
             Beta_Diagonals,        &
             Beta_MVL_Banded,       &
             3*Beta_Diagonals+1,    &
             Beta_IPIV,             &
             NORM,                  &
             RCOND_One,             &
             WORK,                  &
             RWORK,                 &
             INFO                   )
IF (INFO .NE. 0) THEN
    print*,"ZGBCON has failed with INFO = ",INFO
ELSE
    IF ( Verbose_Flag ) THEN
        PRINT*,"RCOND = ",RCOND_One," Norm = ",Norm
    END IF
END IF
!PRINT*,"STOPing in Factorize_Beta_Banded"
!STOP

timer(2) = MPI_Wtime()
CALL Clock_In(Timer(2)-Timer(1),19)

END SUBROUTINE Factorize_Beta_Banded







!+501+###########################################################################!
!                                                                                !
!                   Jacobi_PC_MVL_Banded                                       !
!                                                                                !
!################################################################################!
SUBROUTINE Jacobi_PC_MVL_Banded()

INTEGER                                         ::  i, j
INTEGER                                         ::  Row, Col
INTEGER                                         ::  RE, REp
INTEGER                                         ::  ui, d, l
INTEGER                                         ::  uj, dp, lp

! Store inverse diagonal
DO i = 1,Beta_Prob_Dim

    Row = Beta_Bandwidth + i
    Beta_MVL_Diagonal(i) = 1.0_idp/Beta_MVL_Banded(Beta_Bandwidth,i)

END DO






! Row : RE = 0, D = 0
DO ui = 1,3
DO l  = 1,LM_Length

    Row = Beta_Bandwidth + FP_Beta_Array_Map(0,0,ui,l)


    ! Column : RE = 0, D = 0
    DO uj = 1,3
    DO lp = 1,LM_Length
    
        Col = FP_Beta_Array_Map(0,0,uj,lp)

        Beta_MVL_Banded(Row-Col,Col) = Beta_MVL_Banded(Row-Col,Col)*Beta_MVL_Diagonal(Row-Beta_Bandwidth)

!        PRINT*,Row-Beta_Bandwidth,Col,Beta_MVL_Banded(Row-Col,Col)

    END DO  ! lp
    END DO  ! uj


END DO ! ui
END DO ! l







DO RE = 0,Num_R_Elements-1

    ! Row : RE , D = 0
    DO ui = 1,3
    DO l  = 1,LM_Length

        Row = Beta_Bandwidth + FP_Beta_Array_Map(re, 0, ui, l)


        ! Column : RE , D = 1,DEGREE
        DO dp = 1,DEGREE
        DO uj = 1,3
        DO lp = 1,LM_Length

            Col = FP_Beta_Array_Map(RE,dp,uj,lp)

            Beta_MVL_Banded(Row-Col,Col) = Beta_MVL_Banded(Row-Col,Col)*Beta_MVL_Diagonal(Row-Beta_Bandwidth)
!            PRINT*,Row-Beta_Bandwidth,Col,Beta_MVL_Banded(Row-Col,Col)

        END DO  ! lp
        END DO  ! uj
        END DO  ! dp

    END DO ! l
    END DO ! ui




    ! Row RE D = 1,Degree
    DO d  = 1,DEGREE
    DO ui = 1,3
    DO l  = 1,LM_Length

        Row = Beta_Bandwidth + FP_Beta_Array_Map(RE,d,ui,l)


        ! Column RE, D = 0
        DO uj = 1,3
        DO lp = 1,LM_Length

            Col = FP_Beta_Array_Map(RE,0,uj,lp)

            Beta_MVL_Banded(Row-Col,Col) = Beta_MVL_Banded(Row-Col,Col)*Beta_MVL_Diagonal(Row-Beta_Bandwidth)
!            PRINT*,Row-Beta_Bandwidth,Col,Beta_MVL_Banded(Row-Col,Col)


        END DO  ! lp
        END DO  ! uj





        ! Column RE, D = 1,DEGREE
        DO dp = 1,DEGREE
        DO uj = 1,3
        DO lp = 1,LM_Length

            Col = FP_Beta_Array_Map(RE,dp,uj,lp)

            Beta_MVL_Banded(Row-Col,Col) = Beta_MVL_Banded(Row-Col,Col)*Beta_MVL_Diagonal(Row-Beta_Bandwidth)
!            PRINT*,Row-Beta_Bandwidth,Col,Beta_MVL_Banded(Row-Col,Col)


        END DO  ! lp
        END DO  ! uj
        END DO  ! dp



    END DO  ! l
    END DO  ! ui
    END DO  ! d


END DO ! RE





END SUBROUTINE Jacobi_PC_MVL_Banded









!+501+###########################################################################!
!                                                                                !
!                   Jacobi_PC_MVL_Banded                                       !
!                                                                                !
!################################################################################!
SUBROUTINE Jacobi_PC_MVL_Banded_Vector( Work_Vec )

COMPLEX(idp),   DIMENSION(1:Beta_Prob_Dim),     INTENT(INOUT)       ::  Work_Vec



! Multiply diagonal and work vec

Work_Vec(:)  = Work_Vec(:)*Beta_MVL_Diagonal(:)


END SUBROUTINE Jacobi_PC_MVL_Banded_Vector








!+501+###########################################################################!
!                                                                                !
!                   DIRICHLET_BC_Beta_Banded_Mat                                       !
!                                                                                !
!################################################################################!
SUBROUTINE DIRICHLET_BC_Beta_Banded_Mat()


INTEGER                                                 :: Col, Row
INTEGER                                                 :: re, d, ui, lm
INTEGER                                                 :: rep, dp, uj, lpmp, lp

INTEGER                                                 ::  i, j



DO ui = 1,3

    !############################################!
    !#                                          #!
    !#          Inner Boundary Condition        #!
    !#                                          #!
    !############################################!
    IF (INNER_CFA_BC_TYPE(ui) == "D") THEN

        Col = FP_Beta_Array_Map(0,0,ui,0)
    

        DO d = 0,Degree
        DO lm = 1,LM_Length

            Row = Beta_Bandwidth                    &
                + FP_Beta_Array_Map(0,d,ui,lm)

            First_Column_Beta_Storage(lm,d,ui) = Beta_MVL_Banded(Row-Col,Col)
   
        END DO ! l Loop
        END DO ! d Loop



       DO d = 0,Degree
       DO lm = 1,LM_Length

           Row = Beta_Bandwidth                    &
               + FP_Beta_Array_Map(0,d,ui,lm)


           Beta_MVL_Banded(Row-Col,Col) = 0.0_idp

       END DO ! l Loop
       END DO ! d Loop




       DO d = 0,Degree
       DO lm = 1,LM_Length

           Row = Beta_Bandwidth + FP_Beta_Array_Map(0,0,ui,0)
           Col = FP_Beta_Array_Map(0,d,ui,lm)

           Beta_MVL_Banded(Row-Col,Col) = 0.0_idp

       END DO ! l Loop
       END DO ! d Loop

       

       DO lm = 1,LM_Length

           Row = FP_Beta_Array_Map(0,0,ui,lm)
           Col = FP_Beta_Array_Map(0,0,ui,lm)
 
           Beta_MVL_Banded(Row-Col,Col) = 1.0_idp

       END DO ! l Loop
        
    END IF






    !############################################!
    !#                                          #!
    !#          Outer Boundary Condition        #!
    !#                                          #!
    !############################################!
    IF (OUTER_CFA_BC_TYPE(ui) == "D") THEN

        !
        ! Save the needed values first
        !

        DO d = 0,Degree
        DO lm = 1,LM_Length

            Row = FP_Beta_Array_Map(Num_R_Elements-1,d,ui,lm)+Beta_Bandwidth
            Col = FP_Beta_Array_Map(Num_R_Elements-1,Degree,ui,1)
 
            Last_Column_Beta_Storage(lm,d,ui) = Beta_MVL_Banded(Row-Col,Col)
    
        END DO ! l Loop
        END DO ! d Loop



        ! Clear the Rows
        DO lm = 1,LM_Length

            Row = FP_Beta_Array_Map(Num_R_Elements-1,Degree,ui,lm) + Beta_Bandwidth

            DO uj = 1,3
            DO dp = 0,Degree
            DO lp = 1,LM_Length

                Col = FP_Beta_Array_Map(Num_R_Elements-1,dp,uj,lp)

                Beta_MVL_Banded(Row-Col,Col) = 0.0_idp

            END DO ! lp
            END DO ! dp
            END DO ! uj
        END DO ! l Loop



        ! Clear the Columns
        DO lm = 1,LM_Length

            Col = FP_Beta_Array_Map(Num_R_Elements-1,Degree,ui,lm)

            Beta_MVL_Banded(:,Col) = 0.0_idp

        END DO ! l Loop




        
        DO lm = 1,LM_Length

            Row = Beta_Bandwidth + FP_Beta_Array_Map(Num_R_Elements-1,Degree,ui,lm)
            Col = FP_Beta_Array_Map(Num_R_Elements-1,Degree,ui,lm)

            Beta_MVL_Banded(Row-Col,Col) = 1.0_idp

        END DO ! lm Loop


    END IF

END DO


END SUBROUTINE DIRICHLET_BC_Beta_Banded_Mat


END MODULE FP_Beta_Banded
