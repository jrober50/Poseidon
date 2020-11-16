   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE FP_Functions_Laplace_Beta                                                         !##!
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
                    L_LIMIT

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
                    Laplace_Matrix_Beta,        &
                    CFA_EQ_Flags,               &
                    CFA_EQ_Map,                 &
                    CFA_Mat_Map

USE Poseidon_IO_Module, &
            ONLY :  Open_New_File


IMPLICIT NONE

CONTAINS






!+201+###########################################################################!
!                                                                                !
!           Initialize_Laplace_Matrices_Full                                     !
!                                                                                !
!################################################################################!
SUBROUTINE Initialize_Laplace_Matrices_Beta()

INTEGER                                                 ::  l, m, lp, mp, lm_loc, lpmp_loc
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



Laplace_Matrix_Beta = 0.0_idp

ALLOCATE( CUR_R_LOCS(1:Num_R_Quad_Points) )
ALLOCATE( CUR_T_LOCS(1:Num_T_Quad_Points) )
ALLOCATE( R_SQUARE(1:Num_R_Quad_Points)   )

ALLOCATE( SIN_SQUARE( 1:NUM_TP_QUAD_POINTS ) )
ALLOCATE( COTAN_VAL( 1:NUM_TP_QUAD_POINTS ) )

ALLOCATE( RR_Factor( 1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE )    )
ALLOCATE( RDR_Factor( 1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE )    )
ALLOCATE( DRR_Factor( 1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE )    )
ALLOCATE( DRDR_Factor( 1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE )    )

ALLOCATE( TP_TP_Factor( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1) )
ALLOCATE( dTP_TP_Factor( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1) )
ALLOCATE( TdP_TP_Factor( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1) )
ALLOCATE( TP_dTP_Factor( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1) )
ALLOCATE( TP_TdP_Factor( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1) )
ALLOCATE( dTP_dTP_Factor( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1) )
ALLOCATE( dTP_TdP_Factor( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1) )
ALLOCATE( TdP_dTP_Factor( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1) )
ALLOCATE( TdP_TdP_Factor( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1) )

PRINT*,"Initializing Beta Matrix..."




DO re = 0,NUM_R_ELEMENTS-1
    DO d = 0,DEGREE
        DO ui = 1,3
            DO l = 0,L_LIMIT
                L_Lp1 = REAL( l*(l+1), idp )

                deltar = rlocs(re+1) - rlocs(re)
                TODR = 2.0_idp/deltar
                CUR_R_LOCS(:) = (deltar/2.0_idp) * (Int_R_Locations(:)+1.0_idp) + rlocs(re)
                R_SQUARE = CUR_R_LOCS**2

    
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
                    j = (re*Degree + d)*LM_Length*3      &
                      + (ui - 1) * LM_Length            &
                      + (l*(l+1) + m ) + 1
                    DO dp = 0,Degree
                        i = (re*Degree + dp)*LM_Length*3      &
                          + (ui - 1) * LM_Length            &
                          + (l*(l+1) + m ) + 1

                        Laplace_Matrix_Beta(i,j) = Laplace_Matrix_Beta(i,j)     &
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

!Call Output_Laplace_Beta(Laplace_Matrix_Beta,Beta_Prob_Dim, Beta_Prob_Dim, "A")


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


!Call Output_Laplace_Beta(Laplace_Matrix_Beta,Beta_Prob_Dim, Beta_Prob_Dim,"B")

!PRINT*,"STOPPING at end of Initialize_Laplace_Matrices_Beta"
!STOP

END SUBROUTINE Initialize_Laplace_Matrices_Beta






!+403+###########################################################################!
!                                                                                !
!                   OUTPUT_JACOBIAN_MATRIX                                       !
!                                                                                !
!################################################################################!
SUBROUTINE OUTPUT_LAPLACE_Beta( Matrix, row_in, col_in, flag )

INTEGER,                        INTENT(IN)              ::  row_in, col_in
COMPLEX(idp), DIMENSION(1:Row_In,1:Col_in), INTENT(IN)                ::  Matrix


CHARACTER(LEN = 1 ), INTENT(IN)                          ::  flag

CHARACTER(LEN = 300)                                     ::  FILE_NAME
CHARACTER(LEN = 300)                                     ::  FILE_NAMEb
CHARACTER(LEN = 40)                                     ::  fmt


INTEGER                                                 ::  rows, cols

INTEGER                                                 ::  FILE_ID
INTEGEr                                                 ::  i,j

100 FORMAT (A,A,A,A,A,A)

fmt = '(ES24.16E3,SP,ES24.16E3,"i")'

cols = size( Laplace_Matrix_Beta, 1 )
rows = size( Laplace_Matrix_Beta, 2 )



WRITE(FILE_NAMEb,100) Poseidon_LinSys_Dir,"Beta_Laplace_Dim_",trim(File_Suffix),"_",flag,".out"
CALL OPEN_NEW_FILE( trim(FILE_NAMEb), FILE_ID, 300 )
WRITE(FILE_ID,*) Cols, Rows, NUM_R_Elements, DEGREE, L_LIMIT
CLOSE(FILE_ID)



WRITE(FILE_NAME,100) Poseidon_LinSys_Dir,"Beta_Laplace_",trim(File_Suffix),"_",flag,".out"
CALL OPEN_NEW_FILE( trim(FILE_NAME), FILE_ID, 300 )
DO i = 1,rows
    DO j = 1, cols
        WRITE(FILE_ID,fmt) Matrix(j,i)
    END DO
END DO
CLOSE( FILE_ID )



END SUBROUTINE OUTPUT_LAPLACE_Beta



!+403+###########################################################################!
!                                                                                !
!                   OUTPUT_JACOBIAN_MATRIX                                       !
!                                                                                !
!################################################################################!
SUBROUTINE OUTPUT_LAPLACE( Matrix, row_in, col_in, flag )

INTEGER,                        INTENT(IN)              ::  row_in, col_in
COMPLEX(idp), DIMENSION(0:Row_In-1,0:Col_in-1), INTENT(IN)                ::  Matrix


CHARACTER(LEN = 1 ), INTENT(IN)                          ::  flag

CHARACTER(LEN = 300)                                     ::  FILE_NAME
CHARACTER(LEN = 300)                                     ::  FILE_NAMEb
CHARACTER(LEN = 40)                                     ::  fmt


INTEGER                                                 ::  rows, cols

INTEGER                                                 ::  FILE_ID
INTEGEr                                                 ::  i,j

100 FORMAT (A,A,A,A,A,A)

fmt = '(ES24.16E3,SP,ES24.16E3,"i")'

cols = size( Laplace_Matrix_Beta, 1 )
rows = size( Laplace_Matrix_Beta, 2 )



WRITE(FILE_NAMEb,100) Poseidon_LinSys_Dir,"Beta_Laplace_Dim_",trim(File_Suffix),"_",flag,".out"
CALL OPEN_NEW_FILE( trim(FILE_NAMEb), FILE_ID, 300 )
WRITE(FILE_ID,*) Cols, Rows, NUM_R_Elements, DEGREE, L_LIMIT
CLOSE(FILE_ID)



WRITE(FILE_NAME,100) Poseidon_LinSys_Dir,"Beta_Laplace_",trim(File_Suffix),"_",flag,".out"
CALL OPEN_NEW_FILE( trim(FILE_NAME), FILE_ID, 300 )
DO i = 0,row_in-1
    DO j = 0, col_in-1
        WRITE(FILE_ID,fmt) Matrix(j,i)
    END DO
END DO
CLOSE( FILE_ID )



END SUBROUTINE OUTPUT_LAPLACE






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
!                   OUTPUT_JACOBIAN_MATRIX                                       !
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

COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1), INTENT(INOUT) :: TP_TP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1), INTENT(INOUT) :: dTP_TP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1), INTENT(INOUT) :: TdP_TP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1), INTENT(INOUT) :: TP_dTP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1), INTENT(INOUT) :: dTP_TdP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1), INTENT(INOUT) :: TP_TdP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1), INTENT(INOUT) :: TdP_TdP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1), INTENT(INOUT) :: dTP_dTP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1), INTENT(INOUT) :: TdP_dTP_Factor


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



DO lm_loc = 0,LM_LENGTH-1
    DO lpmp_loc = 0,LM_LENGTH-1

        TP_TP_Factor( :, lm_loc, lpmp_loc )  = Ylm_Values( lm_loc, :, te, pe )              &
                                                * Ylm_CC_Values( :, lpmp_loc, te, pe)       &
                                                * TP_Int_Weights(:)

        dTP_TP_Factor( :, lm_loc, lpmp_loc )  = Ylm_dt_Values( lm_loc, :, te, pe )           &
                                                * Ylm_CC_Values( :, lpmp_loc, te, pe)       &
                                                * TP_Int_Weights(:)

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







!+403+###########################################################################!
!                                                                                !
!                   Calc_Beta1_Terms                                       !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Beta1_Terms( re, d, l, m,                                       &
                             RR_Factor, dRR_Factor, dRdR_Factor,                &
                             TP_TP_Factor, dTP_TP_Factor, TdP_TP_Factor,        &
                             TP_dTP_Factor, TP_TdP_Factor,                      &
                             Cur_R_Locs, R_Square, Cotan_Val )

INTEGER,                                                                      INTENT(IN)    :: re, d, l, m

REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ),               INTENT(IN)    :: RR_Factor
REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ),               INTENT(IN)    :: DRR_Factor
REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ),               INTENT(IN)    :: DRDR_Factor

COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1), INTENT(IN)    :: TP_TP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1), INTENT(IN)    :: dTP_TP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1), INTENT(IN)    :: TdP_TP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1), INTENT(IN)    :: TP_dTP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1), INTENT(IN)    :: TP_TdP_Factor


REAL(idp),    DIMENSION( 1:Num_R_Quad_Points ),                               INTENT(IN)    :: Cur_R_Locs
REAL(idp),    DIMENSION( 1:Num_R_Quad_Points ),                               INTENT(IN)    :: R_Square
REAL(idp),    DIMENSION( 1:Num_TP_Quad_Points ),                              INTENT(IN)    :: Cotan_Val


INTEGER                                                     :: dp, lp, mp, rd, ui, uj
INTEGER                                                     :: row, col


ui = 1
Row = (re*Degree + d)* 3 * LM_Length       &
        + (ui - 1) * LM_Length             &
        + (l*(l+1) + m ) + 1


DO dp = 0,Degree
uj = 1
DO lp = 0,L_Limit
DO mp = -lp,lp
    Col = (re*Degree + dp) * 3 * LM_Length   &
        + (uj - 1) * LM_Length              &
        + (lp*(lp+1) + mp ) + 1

    DO rd = 1,Num_R_Quad_Points
        Laplace_Matrix_Beta(Row, Col) = Laplace_Matrix_Beta(Row, Col)           &
                                        - 1.0_idp/3.0_idp                       &   ! Term 1
                                            * SUM( TP_TP_Factor(:,l,lp)  )      &
                                            * dRdR_Factor(rd, d, dp)            &
                                        - 8.0_idp/(3.0_idp * R_Square(rd))      &   ! Term 2
                                            * SUM( TP_TP_Factor(:,l,lp)  )      &
                                            * RR_Factor(rd, d, dp)
    

    END DO ! rd Loop

END DO ! lp Loop
END DO ! mp Loop


uj = 2

DO lp = 0,L_Limit
DO mp = -lp,lp


    Col = (re*Degree + dp) * 3 * LM_Length  &
        + (uj - 1) * LM_Length          &
        + (lp*(lp+1) + mp ) + 1


    DO rd = 1,Num_R_Quad_Points
        Laplace_Matrix_Beta(Row, Col) = Laplace_Matrix_Beta(Row, Col)                       &
                                        - 1.0_idp/3.0_idp                                   &   ! Term 1
                                            * SUM( TP_dTP_Factor(:,l,lp)  )                 &
                                            * dRR_Factor(rd, d, dp)                         &
                                        - 1.0_idp/3.0_idp                                   &   ! Term 2
                                            * SUM( TP_TP_Factor(:,l,lp) * Cotan_Val(:) )    &
                                            * dRR_Factor(rd, d, dp)                         &
                                        - 2.0_idp/CUR_R_LOCS(rd)                            &   ! Term 3
                                            * SUM( dTP_TP_Factor(:,l,lp) )                  &
                                            * RR_Factor(rd, d, dp)                          &
                                        - 2.0_idp/CUR_R_LOCS(rd)                            &   ! Term 4
                                            * SUM( TP_TP_Factor(:,l,lp) * COTAN_VAL(:)  )   &
                                            * RR_Factor(rd, d, dp)

    END DO ! rd Loop

END DO ! lp Loop
END DO ! mp Loop



uj = 3 ! beta^phi

DO lp = 0,L_Limit
DO mp = -lp,lp
    Col = (re*Degree + dp) * 3 * LM_Length  &
        + (uj - 1) * LM_Length          &
        + (lp*(lp+1) + mp ) + 1


    DO rd = 1,Num_R_Quad_Points
        Laplace_Matrix_Beta(Row, Col) = Laplace_Matrix_Beta(Row, Col)           &
                                        - 1.0_idp/3.0_idp                       &   ! Term 1
                                            * SUM( TP_TdP_Factor(:,l,lp) )      &
                                            * dRR_Factor(rd, d, dp)             &
                                        - 2.0_idp/CUR_R_LOCS(rd)                &   ! Term 2
                                            * SUM( TdP_TP_Factor(:,l,lp) )      &
                                            * RR_Factor(rd, d, dp)

    END DO ! rd Loop
END DO ! lp Loop
END DO ! mp Loop

END DO ! dp Loop


END SUBROUTINE Calc_Beta1_Terms





!+403+###########################################################################!
!                                                                                !
!                   Calc_Beta2_Terms                                       !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Beta2_Terms( re, d, l, m,                                       &
                             RR_Factor, dRR_Factor, dRdR_Factor, RdR_Factor,    &
                             TP_TP_Factor, dTP_TP_Factor, TdP_TP_Factor,        &
                             TP_dTP_Factor, TP_TdP_Factor, dTP_dTP_Factor,      &
                             dTP_TdP_Factor,                                    &
                             Cur_R_Locs, R_Square, Sin_Square, Cotan_Val )

INTEGER,                                                                      INTENT(IN)    :: re, d, l, m

REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ),               INTENT(IN)    :: RR_Factor
REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ),               INTENT(IN)    :: DRR_Factor
REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ),               INTENT(IN)    :: RdR_Factor
REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ),               INTENT(IN)    :: DRDR_Factor

COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1), INTENT(IN)    :: TP_TP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1), INTENT(IN)    :: dTP_TP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1), INTENT(IN)    :: TdP_TP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1), INTENT(IN)    :: TP_dTP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1), INTENT(IN)    :: TP_TdP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1), INTENT(IN)    :: dTP_dTP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1), INTENT(IN)    :: dTP_TdP_Factor

REAL(idp),    DIMENSION( 1:Num_R_Quad_Points ),                               INTENT(IN)    :: Cur_R_Locs
REAL(idp),    DIMENSION( 1:Num_R_Quad_Points ),                               INTENT(IN)    :: R_Square
REAL(idp),    DIMENSION( 1:Num_TP_Quad_Points ),                              INTENT(IN)    :: Sin_Square
REAL(idp),    DIMENSION( 1:Num_TP_Quad_Points ),                              INTENT(IN)    :: Cotan_Val


INTEGER                                                     :: dp, lp, mp, rd, ui, uj
INTEGER                                                     :: row, col


                    ! ui = 1
ui = 2
Row = (re*Degree + d)* 3 * LM_Length       &
        + (ui - 1) * LM_Length             &
        + (l*(l+1) + m ) + 1


DO dp = 0,Degree

uj = 1

DO lp = 0,L_Limit
DO mp = -lp,lp
    Col = (re*Degree + dp) * 3 * LM_Length      &
        + (uj - 1) * LM_Length                  &
        + (lp*(lp+1) + mp ) + 1

    DO rd = 1,Num_R_Quad_Points
        Laplace_Matrix_Beta(Row, Col) = Laplace_Matrix_Beta(Row, Col)               &
                                        - 1.0_idp/(3.0_idp*R_Square(rd) )           &           ! Term 1
                                            * SUM( dTP_TP_Factor(:,l,lp))           &
                                            * RdR_Factor(rd, d, dp)                 &
                                        + 8.0_idp/(3.0_idp*R_Square(rd)* Cur_R_Locs(rd))    &   ! Term 2
                                            * SUM( dTP_TP_Factor(:,l,lp)  )         &
                                            * RR_Factor(rd, d, dp)



    END DO ! rd Loop
END DO ! lp Loop
END DO ! mp Loop



uj = 2

DO lp = 0,L_Limit
DO mp = -lp,lp
    Col = (re*Degree + dp) * 3 * LM_Length  &
        + (uj - 1) * LM_Length          &
        + (lp*(lp+1) + mp ) + 1

    DO rd = 1,Num_R_Quad_Points
        Laplace_Matrix_Beta(Row, Col) = Laplace_Matrix_Beta(Row, Col)           &
                                        - 1.0_idp/(R_Square(rd) )               &   ! Term 1
                                            * SUM( dTP_dTP_Factor(:,l,lp) )     &
                                            * RR_Factor(rd, d, dp)              &
                                        - 1.0_idp/(3.0_idp*R_Square(rd) )       &   ! Term 2
                                            * SUM( dTP_TP_Factor(:, l, lp) )    &
                                            * RR_Factor(rd, d, dp)              &
                                        + 2.0_idp/Cur_R_Locs(rd)                &   ! Term 3
                                            * SUM( TP_TP_Factor(:,l,lp)  )      &
                                            * dRR_Factor(rd, d, dp)             &
                                        - 1.0_idp/(3.0_idp*R_Square(rd))        &   ! Term 4
                                            * SUM( TP_TP_Factor(:,l,lp) / Sin_Square(:) )     &
                                            * RR_Factor(rd, d, dp)              &
                                        + 1.0_idp/R_Square(rd)                  &   ! Term 5
                                            * SUM( TP_TP_Factor(:,l,lp) * (1-Cotan_Val(:)**2) )     &
                                            * RR_Factor(rd, d, dp)


    END DO ! rd Loop
END DO ! lp Loop
END DO ! mp Loop



uj = 3 ! beta^phi

DO lp = 0,L_Limit
DO mp = -lp,lp
    Col = (re*Degree + dp) * 3 * LM_Length  &
        + (uj - 1) * LM_Length          &
        + (lp*(lp+1) + mp ) + 1

    DO rd = 1,Num_R_Quad_Points
        Laplace_Matrix_Beta(Row, Col) = Laplace_Matrix_Beta(Row, Col)       &
                                        - 1.0_idp/(3.0_idp*R_Square(rd) )   &   ! Term 1
                                            * SUM( dTP_TdP_Factor(:,l,lp) ) &
                                            * RR_Factor(rd, d, dp)          &
                                        - 2.0_idp/R_Square(rd)              &   ! Term 2
                                            * SUM( TdP_TP_Factor(:,l,lp) * Cotan_Val(:)   )   &
                                            * RR_Factor(rd, d, dp)


    END DO ! rd Loop
END DO ! lp Loop
END DO ! mp Loop

END DO ! dp Loop



END SUBROUTINE Calc_Beta2_Terms













!+403+###########################################################################!
!                                                                                !
!                   Calc_Beta3_Terms                                       !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Beta3_Terms( re, d, l, m,                                       &
                             RR_Factor, dRR_Factor, dRdR_Factor, RdR_Factor,    &
                             TP_TP_Factor, dTP_TP_Factor, TdP_TP_Factor,        &
                             TP_dTP_Factor, TP_TdP_Factor, TdP_TdP_Factor,      &
                             TdP_dTP_Factor,                                    &
                             Cur_R_Locs, R_Square, Sin_Square, Cotan_Val )

INTEGER,                                                                      INTENT(IN)    :: re, d, l, m

REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ),               INTENT(IN)    :: RR_Factor
REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ),               INTENT(IN)    :: DRR_Factor
REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ),               INTENT(IN)    :: RDR_Factor
REAL(idp), DIMENSION(1:NUM_R_QUAD_POINTS, 0:DEGREE, 0:DEGREE ),               INTENT(IN)    :: DRDR_Factor

COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1), INTENT(IN)    :: TP_TP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1), INTENT(IN)    :: dTP_TP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1), INTENT(IN)    :: TdP_TP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1), INTENT(IN)    :: TP_dTP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1), INTENT(IN)    :: TP_TdP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1), INTENT(IN)    :: TdP_TdP_Factor
COMPLEX(idp), DIMENSION( 1:NUM_TP_QUAD_POINTS, 0:LM_LENGTH-1, 0:LM_LENGTH-1), INTENT(IN)    :: TdP_dTP_Factor

REAL(idp),    DIMENSION( 1:Num_R_Quad_Points ),                               INTENT(IN)    :: Cur_R_Locs
REAL(idp),    DIMENSION( 1:Num_R_Quad_Points ),                               INTENT(IN)    :: R_Square
REAL(idp),    DIMENSION( 1:Num_TP_Quad_Points ),                              INTENT(IN)    :: Cotan_Val
REAL(idp),    DIMENSION( 1:Num_TP_Quad_Points ),                              INTENT(IN)    :: Sin_Square

INTEGER                                                     :: dp, lp, mp, rd, ui, uj
INTEGER                                                     :: row, col


ui = 3
Row = (re*Degree + d)* 3 * LM_Length       &
        + (ui - 1) * LM_Length             &
        + (l*(l+1) + m ) + 1



DO dp = 0,Degree

uj = 1

DO lp = 0,L_Limit
DO mp = -lp,lp
    Col = (re*Degree + dp) * 3 * LM_Length      &
        + (uj - 1) * LM_Length                  &
        + (lp*(lp+1) + mp ) + 1

    DO rd = 1,Num_R_Quad_Points
        Laplace_Matrix_Beta(Row, Col) = Laplace_Matrix_Beta(Row, Col)               &
                                        - 1.0_idp/(3.0_idp * R_Square(rd) )         &
                                            * SUM( TdP_TP_Factor(:, l, lp) )        &
                                            * RdR_Factor(rd, d, dp)                 &
                                        + 8.0_idp/(3.0_idp*Cur_R_Locs(rd) * R_Square(rd) )  & ! Term 1
                                            * SUM( TdP_TP_Factor(:,l,lp) / Sin_Square(:)     )   &
                                            * RR_Factor(rd, d, dp)


    END DO ! rd Loop
END DO ! lp Loop
END DO ! mp Loop



uj = 2

DO lp = 0,L_Limit
DO mp = -lp,lp
    Col = (re*Degree + dp) * 3 * LM_Length  &
        + (uj - 1) * LM_Length          &
        + (lp*(lp+1) + mp ) + 1

    DO rd = 1,Num_R_Quad_Points
        Laplace_Matrix_Beta(Row, Col) = Laplace_Matrix_Beta(Row, Col)       &
                                        - 1.0_idp/( 3.0_idp * R_Square(rd) )&
                                            * SUM( TdP_dTP_Factor(:,l,lp)/ Sin_Square(:) )   &
                                            * RR_Factor(rd, d, dp)          &
                                        + 8.0_idp/( 3.0_idp * R_Square(rd) )&   ! Term 2
                                            * SUM( TdP_TP_Factor(:,l,lp) * Cotan_Val(:) / Sin_Square(:)     )      &
                                            * RR_Factor(rd, d, dp)



    END DO ! rd Loop
END DO ! lp Loop
END DO ! mp Loop



uj = 3 ! beta^phi

DO lp = 0,L_Limit
DO mp = -lp,lp
    Col = (re*Degree + dp) * 3 * LM_Length  &
        + (uj - 1) * LM_Length          &
        + (lp*(lp+1) + mp ) + 1

    DO rd = 1,Num_R_Quad_Points
        Laplace_Matrix_Beta(Row, Col) = Laplace_Matrix_Beta(Row, Col)       &
                                        - 1.0_idp/(3.0_idp * R_Square(rd) ) &  ! Term 1
                                            * SUM( TdP_TdP_Factor(:,l,lp)/Sin_Square(:) )   &
                                            * RR_Factor(rd, d, dp )         &
                                        + 2.0_idp/(Cur_R_Locs(rd) )         & ! Term 2
                                            * SUM( TP_TP_Factor(:,l,lp) )   &
                                            * dRR_Factor(rd, d, dp)         &
                                        + 2.0_idp/( R_Square(rd) )          & ! Term 3
                                            * SUM( dTP_TP_Factor(:,l,lp) * Cotan_Val(:) )    &
                                            * RR_Factor(rd, d, dp)


    END DO ! rd Loop
END DO ! lp Loop
END DO ! mp Loop

END DO ! dp Loop








END SUBROUTINE Calc_Beta3_Terms








END MODULE FP_Functions_Laplace_Beta
