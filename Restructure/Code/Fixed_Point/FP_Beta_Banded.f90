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
            ONLY :  Open_New_File

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


IMPLICIT NONE

CONTAINS






!+101+###########################################################################!
!                                                                                !
!           Initialize_Laplace_Matrices_Full                                     !
!                                                                                !
!################################################################################!
SUBROUTINE Initialize_Beta_MVL_Banded()



INTEGER                                                 ::  l, m, lp, mp, lm, lm_loc, lpmp_loc
INTEGER                                                 ::  re, te, pe
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

REAL(KIND = idp)                                        ::  deltar_overtwo,     &
                                                            deltat_overtwo,     &
                                                            deltap_overtwo

REAL(KIND = idp), ALLOCATABLE, DIMENSION( :,:,: )       :: RR_Factor
REAL(KIND = idp), ALLOCATABLE, DIMENSION( :,:,: )       :: DRR_Factor
REAL(KIND = idp), ALLOCATABLE, DIMENSION( :,:,: )       :: RDR_Factor
REAL(KIND = idp), ALLOCATABLE, DIMENSION( :,:,: )       :: DRDR_Factor

COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION( :,:,: )    :: TP_TP_Factor
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION( :,:,: )    :: dTP_TP_Factor
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION( :,:,: )    :: TdP_TP_Factor
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION( :,:,: )    :: TP_dTP_Factor
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION( :,:,: )    :: TP_TdP_Factor
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION( :,:,: )    :: dTP_dTP_Factor
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION( :,:,: )    :: dTP_TdP_Factor
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION( :,:,: )    :: TdP_dTP_Factor
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION( :,:,: )    :: TdP_TdP_Factor

INTEGER                                                 :: INFO




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

ALLOCATE( TP_TP_Factor( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length) )
ALLOCATE( dTP_TP_Factor( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length) )
ALLOCATE( TdP_TP_Factor( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length) )
ALLOCATE( TP_dTP_Factor( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length) )
ALLOCATE( TP_TdP_Factor( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length) )
ALLOCATE( dTP_dTP_Factor( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length) )
ALLOCATE( dTP_TdP_Factor( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length) )
ALLOCATE( TdP_dTP_Factor( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length) )
ALLOCATE( TdP_TdP_Factor( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length) )



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
                                 TP_TP_Factor, dTP_TP_Factor, TdP_TP_Factor, TP_dTP_Factor,         &
                                 dTP_TdP_Factor, TP_TdP_Factor, TdP_TdP_Factor, dTP_dTP_Factor,     &
                                 TdP_dTP_Factor )



            DO d = 0,DEGREE
                DO l = 0,L_LIMIT
                DO m = -l,l

                    CALL Calc_Beta1_Terms( re, d, l, m,                                       &
                                            RR_Factor, dRR_Factor, dRdR_Factor,                &
                                            TP_TP_Factor, dTP_TP_Factor, TdP_TP_Factor,        &
                                            TP_dTP_Factor, TP_TdP_Factor,                      &
                                            Cur_R_Locs, R_Square, Cotan_Val )

                
                    CALL Calc_Beta2_Terms( re, d, l, m,                                         &
                                            RR_Factor, dRR_Factor, dRdR_Factor, RdR_Factor,     &
                                            TP_TP_Factor, dTP_TP_Factor, TdP_TP_Factor,         &
                                            TP_dTP_Factor, TP_TdP_Factor, dTP_dTP_Factor,       &
                                            dTP_TdP_Factor,                                     &
                                            Cur_R_Locs, R_Square, Sin_Square, Cotan_Val )



                    CALL Calc_Beta3_Terms( re, d, l, m,                                         &
                                            RR_Factor, dRR_Factor, dRdR_Factor, RdR_Factor,     &
                                            TP_TP_Factor, dTP_TP_Factor, TdP_TP_Factor,         &
                                            TP_dTP_Factor, TP_TdP_Factor, TdP_TdP_Factor,       &
                                            TdP_dTP_Factor,                                     &
                                            Cur_R_Locs, R_Square, Sin_Square, Cotan_Val )



                END DO  ! mp Loop
                END DO  ! lp Loop
            END DO  ! dp Loop
        END DO ! pe Loop
    END DO ! te Loop
END DO ! re Loop


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
                             TP_TP_Factor, dTP_TP_Factor, TdP_TP_Factor,        &
                             TP_dTP_Factor, TP_TdP_Factor,                      &
                             Cur_R_Locs, R_Square, Cotan_Val )

INTEGER,                                                                      INTENT(IN)    :: re, dp, lp, mp

REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ),               INTENT(IN)    :: RR_Factor
REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ),               INTENT(IN)    :: DRR_Factor
REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ),               INTENT(IN)    :: DRDR_Factor

COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length), INTENT(IN)    :: TP_TP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length), INTENT(IN)    :: dTP_TP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length), INTENT(IN)    :: TdP_TP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length), INTENT(IN)    :: TP_dTP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length), INTENT(IN)    :: TP_TdP_Factor


REAL(idp),    DIMENSION( 1:Num_R_Quad_Points ),                               INTENT(IN)    :: Cur_R_Locs
REAL(idp),    DIMENSION( 1:Num_R_Quad_Points ),                               INTENT(IN)    :: R_Square
REAL(idp),    DIMENSION( 1:Num_TP_Quad_Points ),                              INTENT(IN)    :: Cotan_Val

COMPLEX(idp)                                                                                :: tmpA, tmpB, tmpC


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


        DO rd = 1,Num_R_Quad_Points
            Beta_MVL_Banded(Row-Col, Col) = Beta_MVL_Banded(Row-Col, Col)   &
                                    - 1.0_idp/3.0_idp                       &   ! Term 1
                                        * SUM( TP_TP_Factor(:,lm_loc,lpmp_loc)  )      &
                                        * dRdR_Factor(rd, d, dp)            &
                                    - 8.0_idp/(3.0_idp * R_Square(rd))      &   ! Term 2
                                        * SUM( TP_TP_Factor(:,lm_loc,lpmp_loc)  )      &
                                        * RR_Factor(rd, d, dp)
        

        END DO ! rd Loop
    END DO ! lpmp_loc Loop







    uj = 2

    DO lm_loc = 1,LM_Length
        Col = FP_Beta_Array_Map(re,d,uj,lm_loc)


!        tmpA = 0.0_idp
!        tmpB = 0.0_idp
!        tmpC = 0.0_idp
        DO rd = 1,Num_R_Quad_Points
            Beta_MVL_Banded(Row-Col, Col) = Beta_MVL_Banded(Row-Col, Col)               &
                                    - 1.0_idp/3.0_idp                                   &   ! Term 1
                                        * SUM( TP_dTP_Factor(:,lm_loc,lpmp_loc)  )                 &
                                        * dRR_Factor(rd, d, dp)                         &
                                    - 1.0_idp/3.0_idp                                   &   ! Term 2
                                        * SUM( TP_TP_Factor(:,lm_loc,lpmp_loc) * Cotan_Val(:) )    &
                                        * dRR_Factor(rd, d, dp)                         &
                                    - 2.0_idp/CUR_R_LOCS(rd)                            &   ! Term 3
                                        * SUM( dTP_TP_Factor(:,lm_loc,lpmp_loc) )                  &
                                        * RR_Factor(rd, d, dp)                          &
                                    - 2.0_idp/CUR_R_LOCS(rd)                            &   ! Term 4
                                        * SUM( TP_TP_Factor(:,lm_loc,lpmp_loc) * Cotan_Val(:)  )   &
                                        * RR_Factor(rd, d, dp)

!            tmpA = tmpA- 1.0_idp/3.0_idp                                   &   ! Term 2
!                * SUM( TP_TP_Factor(:,lm_loc,lpmp_loc) * Cotan_Val(:) )    &
!                * dRR_Factor(rd, d, dp)
!
!            tmpB = tmpB - 2.0_idp/CUR_R_LOCS(rd)                            &   ! Term 3
!                * SUM( dTP_TP_Factor(:,lm_loc,lpmp_loc) )                  &
!                * RR_Factor(rd, d, dp)
!
!            tmpC = tmpC - 2.0_idp/CUR_R_LOCS(rd)                            &   ! Term 4
!                * SUM( TP_TP_Factor(:,lm_loc,lpmp_loc) * Cotan_Val(:)  )   &
!                * RR_Factor(rd, d, dp)

        END DO ! rd Loop

!        PRINT*,lm_loc,lpmp_loc,tmpA,tmpB,tmpC

!        PRINT*,re,d,lm_loc,lpmp_loc, Beta_MVL_Banded(Row-Col, Col)
!        PRINT*,lm_loc,lpmp_loc,                                 &
!                SUM( TP_dTP_Factor(:,lm_loc,lpmp_loc)  ),       &
!                SUM( TP_TP_Factor(:,lm_loc,lpmp_loc) * Cotan_Val(:) ),  &
!                SUM( dTP_TP_Factor(:,lm_loc,lpmp_loc) )

    END DO ! lpmp_loc Loop
    
!    PRINT*,"Stopping Here"
!    STOP








    uj = 3 ! beta^phi

    DO lm_loc = 1,LM_Length
        Col = FP_Beta_Array_Map(re,d,uj,lm_loc)

        DO rd = 1,Num_R_Quad_Points
            Beta_MVL_Banded(Row-Col, Col) = Beta_MVL_Banded(Row-Col, Col)   &
                                    - 1.0_idp/3.0_idp                       &   ! Term 1
                                        * SUM( TP_TdP_Factor(:,lm_loc,lpmp_loc) )      &
                                        * dRR_Factor(rd, d, dp)             &
                                    - 2.0_idp/CUR_R_LOCS(rd)                &   ! Term 2
                                        * SUM( TdP_TP_Factor(:,lm_loc,lpmp_loc) )      &
                                        * RR_Factor(rd, d, dp)

        END DO ! rd Loop
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
                             TP_TP_Factor, dTP_TP_Factor, TdP_TP_Factor,        &
                             TP_dTP_Factor, TP_TdP_Factor, dTP_dTP_Factor,      &
                             dTP_TdP_Factor,                                    &
                             Cur_R_Locs, R_Square, Sin_Square, Cotan_Val )

INTEGER,                                                                      INTENT(IN)    :: re, dp, lp, mp

REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ),               INTENT(IN)    :: RR_Factor
REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ),               INTENT(IN)    :: DRR_Factor
REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ),               INTENT(IN)    :: RdR_Factor
REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ),               INTENT(IN)    :: DRDR_Factor

COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length), INTENT(IN)    :: TP_TP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length), INTENT(IN)    :: dTP_TP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length), INTENT(IN)    :: TdP_TP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length), INTENT(IN)    :: TP_dTP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length), INTENT(IN)    :: TP_TdP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length), INTENT(IN)    :: dTP_dTP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length), INTENT(IN)    :: dTP_TdP_Factor

REAL(idp),    DIMENSION( 1:Num_R_Quad_Points ),                               INTENT(IN)    :: Cur_R_Locs
REAL(idp),    DIMENSION( 1:Num_R_Quad_Points ),                               INTENT(IN)    :: R_Square
REAL(idp),    DIMENSION( 1:Num_TP_Quad_Points ),                              INTENT(IN)    :: Sin_Square
REAL(idp),    DIMENSION( 1:Num_TP_Quad_Points ),                              INTENT(IN)    :: Cotan_Val


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

    DO rd = 1,Num_R_Quad_Points
        Beta_MVL_Banded(Row-Col, Col) = Beta_MVL_Banded(Row-Col, Col)               &
                                + 1.0_idp/(3.0_idp*R_Square(rd) )                   &  ! Term 1
                                    * SUM( dTP_TP_Factor(:,lm_loc,lpmp_loc))        &
                                    * RdR_Factor(rd, d, dp)                         &
                                - 8.0_idp/(3.0_idp*R_Square(rd)* Cur_R_Locs(rd))    &   ! Term 2
                                    * SUM( dTP_TP_Factor(:,lm_loc,lpmp_loc)  )      &
                                    * RR_Factor(rd, d, dp)

        

    END DO ! rd Loop
!    PRINT*,row, col, lm_loc, lpmp_loc, SUM( dTP_TP_Factor(:,lm_loc,lpmp_loc))
!    PRINT*,re,lm_loc, lpmp_loc,SUM( dTP_TP_Factor(:,lm_loc,lpmp_loc)),row, col, Beta_MVL_Banded(Row, Col)
END DO ! lpmp_loc Loop




uj = 2
DO lm_loc = 1,LM_Length
    Col = FP_Beta_Array_Map(re,d,uj,lm_loc)

    DO rd = 1,Num_R_Quad_Points
        Beta_MVL_Banded(Row-Col, Col) = Beta_MVL_Banded(Row-Col, Col)   &
                                - 1.0_idp/(3.0_idp*R_Square(rd) )       &   ! Term 1
                                    * SUM( dTP_dTP_Factor(:,lm_loc,lpmp_loc) )     &
                                    * RR_Factor(rd, d, dp)              &
                                - 1.0_idp/(3.0_idp*R_Square(rd) )       &   ! Term 2
                                    * SUM( dTP_TP_Factor(:,lm_loc,lpmp_loc) * Cotan_Val(:) )    &
                                    * RR_Factor(rd, d, dp)              &
                                + 2.0_idp/Cur_R_Locs(rd)                &   ! Term 3
                                    * SUM( TP_TP_Factor(:,lm_loc,lpmp_loc)  )      &
                                    * dRR_Factor(rd, d, dp)             &
                                - 1.0_idp/(3.0_idp*R_Square(rd))        &   ! Term 4
                                    * SUM( TP_TP_Factor(:,lm_loc,lpmp_loc) / Sin_Square(:) )     &
                                    * RR_Factor(rd, d, dp)              &
                                + 1.0_idp/R_Square(rd)                  &   ! Term 5
                                    * SUM( TP_TP_Factor(:,lm_loc,lpmp_loc) * (1-Cotan_Val(:)**2) )     &
                                    * RR_Factor(rd, d, dp)


    END DO ! rd Loop
END DO ! lpmp_loc Loop


uj = 3 ! beta^phi

DO lm_loc = 1,LM_Length
    Col = FP_Beta_Array_Map(re,d,uj,lm_loc)

    DO rd = 1,Num_R_Quad_Points
        Beta_MVL_Banded(Row-Col, Col) = Beta_MVL_Banded(Row-Col, Col)       &
                                - 1.0_idp/(3.0_idp*R_Square(rd) )       &   ! Term 1
                                    * SUM( dTP_TdP_Factor(:,lm_loc,lpmp_loc) ) &
                                    * RR_Factor(rd, d, dp)                  &
                                - 2.0_idp/R_Square(rd)                  &   ! Term 2
                                    * SUM( TdP_TP_Factor(:,lm_loc,lpmp_loc) * Cotan_Val(:)   )   &
                                    * RR_Factor(rd, d, dp)


    END DO ! rd Loop
END DO ! lpmp_loc Loop

END DO ! dp Loop



END SUBROUTINE Calc_Beta2_Terms













!+403+###########################################################################!
!                                                                                !
!                   Calc_Beta3_Terms                                       !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Beta3_Terms( re, dp, lp, mp,                                       &
                             RR_Factor, dRR_Factor, dRdR_Factor, RdR_Factor,    &
                             TP_TP_Factor, dTP_TP_Factor, TdP_TP_Factor,        &
                             TP_dTP_Factor, TP_TdP_Factor, TdP_TdP_Factor,      &
                             TdP_dTP_Factor,                                    &
                             Cur_R_Locs, R_Square, Sin_Square, Cotan_Val )

INTEGER,                                                                      INTENT(IN)    :: re, dp, lp, mp

REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ),               INTENT(IN)    :: RR_Factor
REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ),               INTENT(IN)    :: DRR_Factor
REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ),               INTENT(IN)    :: RDR_Factor
REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ),               INTENT(IN)    :: DRDR_Factor

COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length), INTENT(IN)    :: TP_TP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length), INTENT(IN)    :: dTP_TP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length), INTENT(IN)    :: TdP_TP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length), INTENT(IN)    :: TP_dTP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length), INTENT(IN)    :: TP_TdP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length), INTENT(IN)    :: TdP_TdP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length), INTENT(IN)    :: TdP_dTP_Factor

REAL(idp),    DIMENSION( 1:Num_R_Quad_Points ),                               INTENT(IN)    :: Cur_R_Locs
REAL(idp),    DIMENSION( 1:Num_R_Quad_Points ),                               INTENT(IN)    :: R_Square
REAL(idp),    DIMENSION( 1:Num_TP_Quad_Points ),                              INTENT(IN)    :: Cotan_Val
REAL(idp),    DIMENSION( 1:Num_TP_Quad_Points ),                              INTENT(IN)    :: Sin_Square

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

    DO rd = 1,Num_R_Quad_Points
        Beta_MVL_Banded(Row-Col, Col) = Beta_MVL_Banded(Row-Col, Col)               &
                                - 1.0_idp/(3.0_idp * R_Square(rd) )         &
                                    * SUM( TdP_TP_Factor(:, lm_loc, lpmp_loc) )        &
                                    * RdR_Factor(rd, d, dp)                 &
                                + 8.0_idp/(3.0_idp*Cur_R_Locs(rd) * R_Square(rd) )  & ! Term 1
                                    * SUM( TdP_TP_Factor(:,lm_loc,lpmp_loc) / Sin_Square(:)     )   &
                                    * RR_Factor(rd, d, dp)


    END DO ! rd Loop
END DO ! lpmp_loc Loop



uj = 2

DO lm_loc = 1,LM_Length
    Col = FP_Beta_Array_Map(re,d,uj,lm_loc)

    DO rd = 1,Num_R_Quad_Points
        Beta_MVL_Banded(Row-Col, Col) = Beta_MVL_Banded(Row-Col, Col)       &
                                - 1.0_idp/( 3.0_idp * R_Square(rd) )&
                                    * SUM( TdP_dTP_Factor(:,lm_loc,lpmp_loc)/ Sin_Square(:) )   &
                                    * RR_Factor(rd, d, dp)          &
                                + 8.0_idp/( 3.0_idp * R_Square(rd) )&   ! Term 2
                                    * SUM( TdP_TP_Factor(:,lm_loc,lpmp_loc) * Cotan_Val(:) / Sin_Square(:)     )      &
                                    * RR_Factor(rd, d, dp)



    END DO ! rd Loop
END DO ! lpmp_loc Loop



uj = 3 ! beta^phi
DO lm_loc = 1,LM_Length
    Col = FP_Beta_Array_Map(re,d,uj,lm_loc)

    DO rd = 1,Num_R_Quad_Points
        Beta_MVL_Banded(Row-Col, Col) = Beta_MVL_Banded(Row-Col, Col)       &
                                - 1.0_idp/(3.0_idp * R_Square(rd) ) &  ! Term 1
                                    * SUM( TdP_TdP_Factor(:,lm_loc,lpmp_loc)/Sin_Square(:) )   &
                                    * RR_Factor(rd, d, dp )         &
                                + 2.0_idp/(Cur_R_Locs(rd) )         & ! Term 2
                                    * SUM( TP_TP_Factor(:,lm_loc,lpmp_loc) )   &
                                    * dRR_Factor(rd, d, dp)         &
                                + 2.0_idp/( R_Square(rd) )          & ! Term 3
                                    * SUM( dTP_TP_Factor(:,lm_loc,lpmp_loc) * Cotan_Val(:) )    &
                                    * RR_Factor(rd, d, dp)


    END DO ! rd Loop
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
                            TP_TP_Factor, dTP_TP_Factor, TdP_TP_Factor, TP_dTP_Factor,          &
                            dTP_TdP_Factor, TP_TdP_Factor, TdP_TdP_Factor, dTP_dTP_Factor,      &
                            TdP_dTP_Factor )

REAL(idp),                                                                    INTENT(IN)    :: DeltaT_OverTwo
REAL(idp),                                                                    INTENT(IN)    :: DeltaP_OverTwo

REAL(idp),    DIMENSION( 1:Num_T_Quad_Points ),                               INTENT(IN)    :: Cur_T_Locs

INTEGER,                                                                      INTENT(IN)    ::  te, pe

REAL(idp),    DIMENSION( 1:Num_TP_Quad_Points ),                              INTENT(INOUT) :: Sin_Square
REAL(idp),    DIMENSION( 1:Num_TP_Quad_Points ),                              INTENT(INOUT) :: Cotan_Val

COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length), INTENT(INOUT) :: TP_TP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length), INTENT(INOUT) :: dTP_TP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length), INTENT(INOUT) :: TdP_TP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length), INTENT(INOUT) :: TP_dTP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length), INTENT(INOUT) :: dTP_TdP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length), INTENT(INOUT) :: TP_TdP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length), INTENT(INOUT) :: TdP_TdP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length), INTENT(INOUT) :: dTP_dTP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 1:LM_Length, 1:LM_Length), INTENT(INOUT) :: TdP_dTP_Factor


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
        TP_Int_Weights( (td-1)*NUM_P_QUAD_POINTS + pd ) = DSIN( Cur_T_Locs(td) )                &
                                                        * DeltaT_OverTwo * INT_T_WEIGHTS(td)    &
                                                        * DeltaP_OverTwo * INT_P_WEIGHTS(pd)
    END DO
END DO



DO lm_loc = 1,LM_Length
    DO lpmp_loc = 1,LM_Length

        TP_TP_Factor( :, lm_loc, lpmp_loc )  = Ylm_Values( lm_loc, :, te, pe )              &
                                                * Ylm_CC_Values( :, lpmp_loc, te, pe)       &
                                                * TP_Int_Weights(:)

        dTP_TP_Factor( :, lm_loc, lpmp_loc )  = Ylm_dt_Values( lm_loc, :, te, pe )           &
                                                * Ylm_CC_Values( :, lpmp_loc, te, pe)       &
                                                * TP_Int_Weights(:)
        

!        PRINT*,lm_loc,lpmp_loc
!        DO tpd = 1,Num_TP_Quad_Points
!            PRINT*,Ylm_dt_Values( lm_loc, tpd, te, pe ),dTP_TP_Factor(tpd, lm_loc, lpmp_loc )
!        END DO
!        PRINT*," "

        TdP_TP_Factor( :, lm_loc, lpmp_loc )  = Ylm_dp_Values( lm_loc, :, te, pe )           &
                                                * Ylm_CC_Values( :, lpmp_loc, te, pe)       &
                                                * TP_Int_Weights(:)

        TP_dTP_Factor(:, lm_loc, lpmp_loc ) = Ylm_Values( lm_loc, :, te, pe )                &
                                                * Ylm_CC_dt_Values( :, lpmp_loc, te, pe)    &
                                                * TP_Int_Weights(:)

        dTP_dTP_Factor(:, lm_loc, lpmp_loc ) = Ylm_dt_Values( lm_loc, :, te, pe )             &
                                                * Ylm_CC_dt_Values( :, lpmp_loc, te, pe)    &
                                                * TP_Int_Weights(:)

        TP_TdP_Factor(:, lm_loc, lpmp_loc ) = Ylm_Values( lm_loc, :, te, pe )                &
                                                * Ylm_CC_dp_Values( :, lpmp_loc, te, pe)    &
                                                * TP_Int_Weights(:)

        TdP_TdP_Factor(:, lm_loc, lpmp_loc ) = Ylm_dp_Values( lm_loc, :, te, pe )             &
                                                * Ylm_CC_dp_Values( :, lpmp_loc, te, pe)    &
                                                * TP_Int_Weights(:)

        dTP_TdP_Factor(:, lm_loc, lpmp_loc ) = Ylm_dt_Values( lm_loc, :, te, pe )             &
                                                * Ylm_CC_dp_Values( :, lpmp_loc, te, pe)    &
                                                * TP_Int_Weights(:)

        TdP_dTP_Factor(:, lm_loc, lpmp_loc ) = Ylm_dp_Values( lm_loc, :, te, pe )             &
                                                * Ylm_CC_dt_Values( :, lpmp_loc, te, pe)    &
                                                * TP_Int_Weights(:)


    END DO
END DO



END SUBROUTINE CALC_TP_Values


!+501+###########################################################################!
!                                                                                !
!                   Factorize_Beta_Banded                                       !
!                                                                                !
!################################################################################!
SUBROUTINE Factorize_Beta_Banded()


INTEGER                                                 :: ui, Col, Row
INTEGER                                                 :: d, lm

INTEGER                                                 ::  i, j, l, uj
INTEGER                                                 ::  re, dp, lp


INTEGER                                                 :: INFO
COMPLEX(idp),   DIMENSION(2*Beta_Prob_Dim)              :: WORK
REAL(idp),      DIMENSION(Beta_Prob_Dim)                :: RWORK
REAL(idp)                                               :: RCOND_One
REAL(idp)                                               :: NORM


!   Dirichlet BCs modify the stiffness matrix so we modify it now.
!   But to apply the BCs we will need values from the original matrix,
!   so those values are stored before we modify the matrix.
!
CALL DIRICHLET_BC_Beta_Banded_Mat()


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


INTEGER                                                 :: ui, Col, Row
INTEGER                                                 :: d, lm

INTEGER                                                 ::  i, j, l, uj
INTEGER                                                 ::  re, dp, lp


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
        DO d = 0,Degree
        DO lm = 1,LM_Length

            Row = FP_Beta_Array_Map(Num_R_Elements-1,Degree,ui,lm) + Beta_Bandwidth

            DO dp = 0,Degree
            DO lp = 1,LM_Length

                Col = FP_Beta_Array_Map(Num_R_Elements-1,dp,ui,lp)

                Beta_MVL_Banded(Row-Col,Col) = 0.0_idp

            END DO ! lp
            END DO ! dp
        END DO ! l Loop
        END DO ! d Loop



        ! Clear the Columns
        DO lm = 1,LM_Length

            Col = FP_Beta_Array_Map(Num_R_Elements-1,Degree,ui,lm)

            Beta_MVL_Banded(:,Col) = 0.0_idp

        END DO ! l Loop




        
        DO lm = 1,LM_Length

            Row = Beta_Bandwidth + FP_Beta_Array_Map(Num_R_Elements-1,Degree,ui,lm)
            Col = FP_Beta_Array_Map(Num_R_Elements-1,Degree,ui,lm)

            Beta_MVL_Banded(Row-Col,Col) = 1.0_idp

        END DO ! l Loop


    END IF

END DO


END SUBROUTINE DIRICHLET_BC_Beta_Banded_Mat


END MODULE FP_Beta_Banded
