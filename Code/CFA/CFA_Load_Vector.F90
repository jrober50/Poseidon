   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE CFA_Load_Vector_Module                                                       !##!
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
            ONLY :  pi,                         &
                    TwoPi

USE Poseidon_Units_Module, &
            ONLY :  GR_Source_Scalar

USE Poseidon_Parameters, &
            ONLY :  Degree,                     &
                    L_Limit,                    &
                    Eq_Flags

USE Parameters_Variable_Indices, &
            ONLY :  iU_CF,                      &
                    iU_LF,                      &
                    iU_S1,                      &
                    iU_S2,                      &
                    iU_S3,                      &
                    iVB_S


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
                    Int_P_Weights,              &
                    Int_TP_Weights

USE Variables_Mesh, &
            ONLY :  Num_R_Elements,             &
                    Num_T_Elements,             &
                    Num_P_Elements,             &
                    rlocs,                      &
                    tlocs,                      &
                    plocs

USE Variables_Tables, &
            ONLY :  Ylm_Values,                 &
                    Ylm_dt_Values,              &
                    Ylm_dp_Values,              &
                    Ylm_CC_Values,              &
                    Lagrange_Poly_Table

USE Variables_Derived, &
            ONLY :  Num_R_Nodes,                &
                    LM_Length

USE Variables_Vectors, &
            ONLY :  cVA_Coeff_Vector,          &
                    cVB_Coeff_Vector,          &
                    cVA_Load_Vector,         &
                    cVB_Load_Vector

USE Functions_Jacobian, &
            ONLY :  JCBN_kappa_FUNCTION_3D_ALL,     &
                    JCBN_BIGK_FUNCTION

USE Maps_Fixed_Point, &
            ONLY :  FP_Array_Map_TypeB

USE Maps_Domain, &
            ONLY :  Map_To_FEM_Node

USE CFA_Source_Terms_Module, &
            ONLY :  Calc_Source_Terms


USE MPI




IMPLICIT NONE

REAL(KIND = idp)                                        :: Time_S

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CUR_R_LOCS
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CUR_T_LOCS
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CUR_P_LOCS

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: R_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: R_CUBED

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: SIN_VAL_B
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: SIN_VAL
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: SIN_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: TP_SIN_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CSC_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: COTAN_VAL

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)           :: RSIN_SQUARE


REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: R_Int_Weights
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: TP_Int_Weights

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)           :: CUR_VAL_PSI
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)           :: CUR_VAL_ALPHAPSI
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: CUR_VAL_BETA

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: CUR_DRV_PSI
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: CUR_DRV_ALPHAPSI
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:,:)       :: CUR_DRV_BETA

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)           :: Beta_DRV_Trace

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: SOURCE_TERMS

CONTAINS



!+101+###########################################################################!
!                                                                                !
!           Calc_Load_Vector                                                   !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_CFA_Load_Vector()


REAL(KIND = idp),DIMENSION(1:4)                             ::  Timer

INTEGER                                                     ::  re, te, pe,     &
                                                                rd, tpd, td, pd

REAL(KIND = idp)                                                ::  deltar_overtwo,     &
                                                                    deltat_overtwo,     &
                                                                    deltap_overtwo



Timer = 0.0_idp
cVA_Load_Vector = 0.0_idp
cVB_Load_Vector = 0.0_idp


!PRINT*,"**WARNING** Create_CFA_Load_Vector hacked, Lm loop limited."

DO re = 0,Num_R_Elements-1

    DELTAR_OVERTWO = 0.5_idp *(rlocs(re + 1) - rlocs(re))
    CUR_R_LOCS(:) = deltar_overtwo * (Int_R_Locations(:)+1.0_idp) + rlocs(re)


    R_SQUARE(:) = CUR_R_LOCS(:)*CUR_R_LOCS(:)
    R_CUBED(:) = R_SQUARE(:)*CUR_R_LOCS(:)



    DO pe = 0,Num_P_Elements-1
        deltap_overtwo = 0.5_idp * (plocs(pe + 1)-plocs(pe))
        CUR_P_LOCS(:) = deltap_overtwo * (Int_P_Locations+1.0_idp) + plocs(pe)


        DO te = 0,Num_T_Elements-1



            deltat_overtwo = 0.5_idp*(tlocs(te + 1) - tlocs(te))
            CUR_T_LOCS(:) = deltat_overtwo * (Int_T_Locations(:)+1.0_idp) + tlocs(te)
            



            COTAN_VAL(:) = 1.0_idp/DTAN(CUR_T_LOCS(:))
            SIN_VAL(:) = DSIN(CUR_T_LOCS(:))

            SIN_SQUARE(:) = SIN_VAL(:)*SIN_VAL(:)
            CSC_SQUARE(:) = 1.0_idp/SIN_SQUARE(:)
            DO td = 1,Num_T_Quad_Points
            DO pd = 1,Num_P_Quad_Points
                tpd = (td-1)*Num_P_Quad_Points + pd
                SIN_VAL_B(tpd) =   SIN_VAL(td)
                TP_SIN_SQUARE(tpd) = SIN_SQUARE(td)
            END DO
            END DO

            DO rd = 1,Num_R_Quad_Points

                RSIN_SQUARE(:,rd) = R_SQUARE(rd)*SIN_SQUARE(:)

            END DO



!            PRINT*,"re,te,pe",re,te,pe
!            PRINT*,"Before Calc_CFA_Current_Values"
            CALL Calc_CFA_Current_Values(  re, te, pe,     &
                                          DELTAR_OVERTWO, &
                                          DELTAT_OVERTWO, &
                                          DELTAP_OVERTWO  )

!            PRINT*,"Before Calc_CFA_Source_Terms"
            CALL Calc_CFA_Source_Terms(    re, te, pe )

!            PRINT*,"Before Create_CFA_Load_Vector"
            CALL Create_CFA_Load_Vector( re, te, pe,       &
                                        DELTAR_OVERTWO,   &
                                        SIN_VAL_B         )



        END DO ! te Loop
    END DO ! pe Loop


END DO ! re Loop


!PRINT*,"Source Vector"
!DO lm = 1,LM_Length
!    PRINT*,"Lm_loc = ",lm
!    PRINT*,CFA_Load_Vector(:,lm,1)
!END DO



END SUBROUTINE Calc_CFA_Load_Vector










!+202+###########################################################################!
!                                                                                !
!                  Calc_Current_Values          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_CFA_Current_Values( re, te, pe,                                  &
                                    DELTAR_OVERTWO,                             &
                                    DELTAT_OVERTWO,                             &
                                    DELTAP_OVERTWO                              )

INTEGER, INTENT(IN)                                             ::  re, te, pe


REAL(KIND = idp), INTENT(IN)                                    ::  DELTAR_OVERTWO,     &
                                                                    DELTAT_OVERTWO,     &
                                                                    DELTAP_OVERTWO



COMPLEX(KIND = idp), DIMENSION(1:5)                             ::  Tmp_U_Value,        &
                                                                    Tmp_U_R_DRV_Value,  &
                                                                    Tmp_U_T_DRV_Value,  &
                                                                    Tmp_U_P_DRV_Value



INTEGER                                                         ::  tpd, td, pd, rd
INTEGER                                                         ::  ui, d
INTEGER                                                         ::  Here, There


                          !                                                 !
                         !!                                                 !!
                        !!!          Initialize Local Quadratures           !!!
                         !!                                                 !!
                          !                                                 !
R_Int_Weights(:) = DELTAR_OVERTWO * R_SQUARE(:) * Int_R_Weights(:)

DO td = 1,Num_T_Quad_Points
DO pd = 1,Num_P_Quad_Points
!    TP_Int_Weights( (td-1)*Num_P_Quad_Points + pd ) = SIN_VAL(td)                           &
!                                                    * DELTAT_OVERTWO * Int_T_Weights(td)    &
!                                                    * DELTAP_OVERTWO * Int_P_Weights(pd)

    TP_Int_Weights( (td-1)*Num_P_Quad_Points + pd ) = SIN_VAL(td)                           &
                                                    * DELTAT_OVERTWO * Int_T_Weights(td)    &
                                                    * Int_P_Weights(pd)
END DO
END DO





DO rd = 1,Num_R_Quad_Points
DO tpd = 1,Num_TP_Quad_Points

    Tmp_U_Value = 0.0_idp
    Tmp_U_R_DRV_Value = 0.0_idp
    Tmp_U_T_DRV_Value = 0.0_idp
    Tmp_U_P_DRV_Value = 0.0_idp

    
    DO ui = iU_CF,iU_LF
    DO d = 0,DEGREE

        Here = Map_To_FEM_Node(re,d)

    


        TMP_U_Value(ui)         = TMP_U_Value(ui)                           &
                                + SUM( cVA_Coeff_Vector( Here, :, ui )     &
                                * Ylm_Values( :, tpd, te, pe )       )      &
                                * Lagrange_Poly_Table( d, rd, 0 )


        TMP_U_R_DRV_Value(ui)   = TMP_U_R_DRV_Value(ui)                     &
                                + SUM( cVA_Coeff_Vector( Here, :, ui )     &
                                * Ylm_Values( :, tpd, te, pe )       )      &
                                * Lagrange_Poly_Table( d, rd, 1 )           &
                                / DELTAR_OVERTWO


        TMP_U_T_DRV_Value(ui)   = TMP_U_T_DRV_Value(ui)                     &
                                + SUM( cVA_Coeff_Vector( Here, :, ui )     &
                                * Ylm_dt_Values( :, tpd, te, pe)     )      &
                                * Lagrange_Poly_Table( d, rd, 0)

        TMP_U_P_DRV_Value(ui)   = TMP_U_P_DRV_Value(ui)                     &
                                + SUM( cVA_Coeff_Vector( Here, :, ui )     &
                                * Ylm_dp_Values( :, tpd, te, pe)     )      &
                                * Lagrange_Poly_Table( d, rd, 0)




    END DO  ! d
    END DO  ! ui


    DO ui = iU_S1,iU_S3
    DO d = 0,DEGREE
        Here  = FP_Array_Map_TypeB(ui,iVB_S,re,d,1)
        There = FP_Array_Map_TypeB(ui,iVB_S,re,d,LM_Length)



        TMP_U_Value(ui)         = TMP_U_Value(ui)                           &
                                + SUM( cVB_Coeff_Vector(Here:There,iVB_S)  &
                                * Ylm_Values( :, tpd, te, pe )       )      &
                                * Lagrange_Poly_Table( d, rd, 0 )


        TMP_U_R_DRV_Value(ui)   = TMP_U_R_DRV_Value(ui)                     &
                                + SUM( cVB_Coeff_Vector(Here:There,iVB_S)  &
                                * Ylm_Values( :, tpd, te, pe )       )      &
                                * Lagrange_Poly_Table( d, rd, 1 )           &
                                / DELTAR_OVERTWO


        TMP_U_T_DRV_Value(ui)   = TMP_U_T_DRV_Value(ui)                     &
                                + SUM( cVB_Coeff_Vector(Here:There,iVB_S)  &
                                * Ylm_dt_Values( :, tpd, te, pe)     )      &
                                * Lagrange_Poly_Table( d, rd, 0)

        TMP_U_P_DRV_Value(ui)   = TMP_U_P_DRV_Value(ui)                     &
                                + SUM( cVB_Coeff_Vector(Here:There,iVB_S)  &
                                * Ylm_dp_Values( :, tpd, te, pe)     )      &
                                * Lagrange_Poly_Table( d, rd, 0)




    END DO  ! d
    END DO  ! ui

    CUR_VAL_PSI( tpd, rd )         = REAL(Tmp_U_Value(1), KIND = idp)
    CUR_DRV_PSI( tpd, rd, 1 )      = REAL(Tmp_U_R_DRV_Value(1), KIND = idp)
    CUR_DRV_PSI( tpd, rd, 2 )      = REAL(Tmp_U_T_DRV_Value(1), KIND = idp)
    CUR_DRV_PSI( tpd, rd, 3 )      = REAL(Tmp_U_P_DRV_Value(1), KIND = idp)


    CUR_VAL_ALPHAPSI( tpd, rd )    = REAL(Tmp_U_Value(2), KIND = idp)
    CUR_DRV_ALPHAPSI( tpd, rd, 1 ) = REAL(Tmp_U_R_DRV_Value(2), KIND = idp)
    CUR_DRV_ALPHAPSI( tpd, rd, 2 ) = REAL(Tmp_U_T_DRV_Value(2), KIND = idp)
    CUR_DRV_ALPHAPSI( tpd, rd, 3 ) = REAL(Tmp_U_P_DRV_Value(2), KIND = idp)



    CUR_VAL_BETA( tpd, rd, 1 )     = REAL(Tmp_U_Value(3), KIND = idp)
    CUR_VAL_BETA( tpd, rd, 2 )     = REAL(Tmp_U_Value(4), KIND = idp)
    CUR_VAL_BETA( tpd, rd, 3 )     = REAL(Tmp_U_Value(5), KIND = idp)


    CUR_DRV_BETA( tpd, rd, 1, 1 )  = REAL(Tmp_U_R_DRV_Value(3), KIND = idp)
    CUR_DRV_BETA( tpd, rd, 2, 1 )  = REAL(Tmp_U_R_DRV_Value(4), KIND = idp)
    CUR_DRV_BETA( tpd, rd, 3, 1 )  = REAL(Tmp_U_R_DRV_Value(5), KIND = idp)

    CUR_DRV_BETA( tpd, rd, 1, 2 )  = REAL(Tmp_U_T_DRV_Value(3), KIND = idp)
    CUR_DRV_BETA( tpd, rd, 2, 2 )  = REAL(Tmp_U_T_DRV_Value(4), KIND = idp)
    CUR_DRV_BETA( tpd, rd, 3, 2 )  = REAL(Tmp_U_T_DRV_Value(5), KIND = idp)

    CUR_DRV_BETA( tpd, rd, 1, 3 )  = REAL(Tmp_U_P_DRV_Value(3), KIND = idp)
    CUR_DRV_BETA( tpd, rd, 2, 3 )  = REAL(Tmp_U_P_DRV_Value(4), KIND = idp)
    CUR_DRV_BETA( tpd, rd, 3, 3 )  = REAL(Tmp_U_P_DRV_Value(5), KIND = idp)

    Beta_DRV_Trace( tpd, rd )      = CUR_DRV_BETA( tpd, rd, 1, 1 )              &
                                   + CUR_DRV_BETA( tpd, rd, 2, 2 )              &
                                   + CUR_DRV_BETA( tpd, rd, 3, 3 )


END DO ! tpd
END DO ! rd

END SUBROUTINE Calc_CFA_Current_Values





!+202+###########################################################################!
!                                                                                !
!                  Calc_CFA_Source_Terms          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_CFA_Source_Terms( re, te, pe )

INTEGER, INTENT(IN)                                                     ::  re, te, pe
INTEGER                                                                 ::  pd, td, rd,     &
                                                                            i, tpd


REAL(KIND = idp), DIMENSION(1:11)                                       ::  PSI_POWER
REAL(KIND = idp), DIMENSION(1:4)                                        ::  ALPHAPSI_POWER


REAL(KIND = idp), DIMENSION(1:3,1:3)                                    ::  Kappa_Array
REAL(KIND = idp), DIMENSION(1:3)                                        ::  n_Array

REAL(KIND = idp)                                                        ::  BigK_Value


Kappa_Array = 0.0_idp
n_Array     = 0.0_idp
BigK_Value  = 0.0_idp


DO rd = 1,Num_R_Quad_Points
DO td = 1,Num_T_Quad_Points
DO pd = 1,Num_P_Quad_Points

    tpd = (td-1)*Num_P_Quad_Points + pd


    PSI_POWER(1) = CUR_VAL_PSI( tpd, rd)
    DO i = 2,11
        PSI_POWER(i) = PSI_POWER(i-1)*PSI_POWER(1)
    END DO

    ALPHAPSI_POWER(1) = CUR_VAL_ALPHAPSI( tpd, rd)
    DO i = 2,4
        ALPHAPSI_POWER(i) = ALPHAPSI_POWER(i-1)*ALPHAPSI_POWER(1)
    END DO


    ! K_{ij}K^{ij} = Psi^{14}/AlphaPsi^{2} * BIGK
!    PRINT*,"Before BigK_Value"
    BigK_Value = JCBN_BIGK_FUNCTION( rd, tpd,                                                        &
                                     Num_R_Quad_Points, Num_TP_Quad_Points,                          &
                                     CUR_VAL_BETA, CUR_DRV_BETA,                                     &
                                     CUR_R_LOCS(rd), R_SQUARE(rd), SIN_SQUARE(td), CSC_SQUARE(td),   &
                                     RSIN_SQUARE(td, rd), COTAN_VAL(td)                              )




    Kappa_Array = JCBN_kappa_FUNCTION_3D_ALL(  rd, tpd,                                        &
                                               Num_R_Quad_Points, Num_TP_Quad_Points,          &
                                               CUR_R_LOCS(rd), R_SQUARE(rd), R_CUBED(rd),      &
                                               RSIN_SQUARE(td, rd),                            &
                                               COTAN_VAL(td),                                   &
                                               CUR_VAL_BETA, CUR_DRV_BETA                      )


!    PRINT*,"Before n_Array"
    n_Array(:) = CUR_DRV_ALPHAPSI( tpd, rd, : ) / ALPHAPSI_POWER(1)   &
                        - 7.0_idp * CUR_DRV_PSI( tpd, rd, : )/ PSI_POWER(1)



!    PRINT*,"Before Calc_Source_Terms"
    CALL Calc_Source_Terms( Source_Terms,                                      &
                            re, te, pe,                                        &
                            td, pd, tpd, rd,                                   &
                            PSI_POWER, ALPHAPSI_POWER,                         &
                            BigK_Value, n_Array, Kappa_Array    )


    

END DO  ! pd loop
END DO  ! td loop
END DO  ! rd loop




END SUBROUTINE Calc_CFA_Source_Terms









!+204+###########################################################################!
!                                                                                !
!                  CREATE_3D_RHS_VECTOR                                          !
!                                                                                !
!################################################################################!
SUBROUTINE Create_CFA_Load_Vector( re, te, pe, DELTAR_OVERTWO, SIN_VAL )



INTEGER, INTENT(IN)                                                     ::  re, te, pe

REAL(KIND = idp), INTENT(IN)                                            ::  DELTAR_OVERTWO
REAL(KIND = idp), DIMENSION(1:Num_TP_Quad_Points), INTENT(IN)           ::  SIN_VAL

INTEGER                                                                 ::  rd, d, lm_loc, ui

INTEGER                                                                 ::  Current_i_Location

COMPLEX(KIND = idp), DIMENSION(1:5)                                     ::  RHS_TMP



DO ui = iU_CF,iU_LF
    IF ( Eq_Flags(ui) == 1 ) THEN

        
        DO lm_loc = 1,LM_Length
        DO d = 0,DEGREE


            RHS_TMP = 0.0_idp
            DO rd = 1,Num_R_Quad_Points

                RHS_TMP(ui) =  RHS_TMP(ui)                                          &
                                 + SUM( Source_Terms( :, rd, ui )                    &
                                       * Ylm_CC_Values( :, lm_loc, te, pe)         &
                                       * TP_Int_Weights(:)                     )   &
                               * Lagrange_Poly_Table( d, rd, 0)                     &
                               * R_Int_Weights(rd)

            END DO  ! rd Loop
            

            Current_i_Location = Map_To_FEM_Node(re,d)
            cVA_Load_Vector(Current_i_Location,lm_loc,ui)                &
                = cVA_Load_Vector(Current_i_Location,lm_loc,ui)          &
                + RHS_TMP(ui)


        END DO  ! d Loop
        END DO  ! lm_loc Loop
    END IF


END DO





DO ui = iU_S1,iU_S3
    IF( Eq_Flags(ui) == 1 ) THEN

        DO d = 0,DEGREE
        DO lm_loc = 1,LM_Length
        
            RHS_TMP = 0.0_idp
            DO rd = 1,Num_R_Quad_Points

                RHS_TMP(ui) =  RHS_TMP(ui)                                          &
                                + SUM( Source_Terms( :, rd, ui )                    &
                                        * Ylm_CC_Values( :, lm_loc, te, pe)         &
                                        * TP_Int_Weights(:)                     )   &
                                * Lagrange_Poly_Table( d, rd, 0)                     &
                                * R_Int_Weights(rd)


            END DO  ! rd Loop

            Current_i_Location = FP_Array_Map_TypeB(ui,iVB_S,re,d,lm_loc)

            cVB_Load_Vector(Current_i_Location,iVB_S)                &
                = cVB_Load_Vector(Current_i_Location,iVB_S)          &
                + RHS_TMP(ui)


        END DO  ! lm_loc Loop
        END DO  ! d Loop

    END IF

END DO ! ui





END SUBROUTINE Create_CFA_Load_Vector













!+701+###########################################################################!
!                                                                                !
!           Allocate_Master_Build_Variables                                      !
!                                                                                !
!################################################################################!
SUBROUTINE Allocate_CFA_Source_Variables()



ALLOCATE( CUR_R_LOCS(1:Num_R_Quad_Points) )
ALLOCATE( CUR_T_LOCS(1:Num_T_Quad_Points) )
ALLOCATE( CUR_P_LOCS(1:Num_P_Quad_Points) )


ALLOCATE( R_SQUARE(1:Num_R_Quad_Points) )
ALLOCATE( R_CUBED(1:Num_R_Quad_Points) )

ALLOCATE( SIN_VAL_B( 1:Num_TP_Quad_Points ) )
ALLOCATE( SIN_VAL( 1:Num_T_Quad_Points ) )
ALLOCATE( SIN_SQUARE( 1:Num_T_Quad_Points ) )
ALLOCATE( TP_SIN_SQUARE( 1:Num_TP_Quad_Points ) )
ALLOCATE( CSC_SQUARE( 1:Num_T_Quad_Points ) )
ALLOCATE( COTAN_VAL( 1:Num_T_Quad_Points ) )

ALLOCATE( RSIN_SQUARE( 1:Num_T_Quad_Points, 1:Num_R_Quad_Points ) )


ALLOCATE( R_Int_Weights( 1:Num_R_Quad_Points ) )
ALLOCATE( TP_Int_Weights( 1:Num_TP_Quad_Points) )

ALLOCATE( CUR_VAL_PSI(1:Num_TP_Quad_Points, 1:Num_R_Quad_Points) )
ALLOCATE( CUR_VAL_ALPHAPSI(1:Num_TP_Quad_Points, 1:Num_R_Quad_Points) )
ALLOCATE( CUR_VAL_BETA(1:Num_TP_Quad_Points, 1:Num_R_Quad_Points, 1:3) )

ALLOCATE( CUR_DRV_PSI(1:Num_TP_Quad_Points, 1:Num_R_Quad_Points, 1:3) )
ALLOCATE( CUR_DRV_ALPHAPSI(1:Num_TP_Quad_Points, 1:Num_R_Quad_Points, 1:3) )
ALLOCATE( CUR_DRV_BETA(1:Num_TP_Quad_Points, 1:Num_R_Quad_Points, 1:3, 1:3) )

ALLOCATE( Beta_DRV_Trace(1:Num_TP_Quad_Points, 1:Num_R_Quad_Points) )

ALLOCATE( Source_Terms( 1:Num_TP_Quad_Points, 1:Num_R_Quad_Points, 1:5) )


END SUBROUTINE Allocate_CFA_Source_Variables



!+702+###########################################################################!
!                                                                                !
!           Deallocate_Master_Build_Variables                                    !
!                                                                                !
!################################################################################!
SUBROUTINE Deallocate_CFA_Source_Variables()

DEALLOCATE( CUR_R_LOCS )
DEALLOCATE( CUR_T_LOCS )
DEALLOCATE( CUR_P_LOCS )

DEALLOCATE( R_SQUARE )
DEALLOCATE( R_CUBED )

DEALLOCATE( SIN_VAL_B )
DEALLOCATE( SIN_VAL )
DEALLOCATE( SIN_SQUARE )
DEALLOCATE( TP_SIN_SQUARE )
DEALLOCATE( CSC_SQUARE )
DEALLOCATE( COTAN_VAL )

DEALLOCATE( RSIN_SQUARE )


DEALLOCATE( R_Int_Weights )
DEALLOCATE( TP_Int_Weights )

DEALLOCATE( CUR_VAL_PSI )
DEALLOCATE( CUR_VAL_ALPHAPSI )
DEALLOCATE( CUR_VAL_BETA )

DEALLOCATE( CUR_DRV_PSI )
DEALLOCATE( CUR_DRV_ALPHAPSI )
DEALLOCATE( CUR_DRV_BETA )

DEALLOCATE( Beta_DRV_Trace )

DEALLOCATE( Source_Terms )


END SUBROUTINE Deallocate_CFA_Source_Variables

END MODULE CFA_Load_Vector_Module
