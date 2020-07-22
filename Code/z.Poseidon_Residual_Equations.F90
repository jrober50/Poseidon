   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_Residual_Equations_Module                                           !##!
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

USE MPI

USE OMP_LIB

USE Units_Module, &
                                ONLY :  GR_Source_Scalar


USE Poseidon_Constants_Module, &
                                ONLY :  idp,                &
                                        pi,                 &
                                        TwoPi,              &
                                        OneThird,           &
                                        TwoThirds,          &
                                        FourThirds,         &
                                        OneThirtySecond

USE Poseidon_Parameters, &
                                ONLY :  DOMAIN_DIM,                 &
                                        DEGREE,                     &
                                        L_LIMIT,                    &
                                        NUM_SHELLS,                 &
                                        NUM_SUBSHELLS,              &
                                        NUM_SUBSHELLS_PER_SHELL,    &
                                        NUM_CFA_VARS,               &
                                        nPROCS_POSEIDON,            &
                                        NUM_R_QUAD_POINTS,          &
                                        NUM_T_QUAD_POINTS,          &
                                        NUM_P_QUAD_POINTS,          &
                                        NUM_QUAD_DOF,               &
                                        NUM_R_ELEMS_PER_BLOCK,      &
                                        NUM_T_ELEMS_PER_BLOCK,      &
                                        NUM_P_ELEMS_PER_BLOCK,      &
                                        NUM_R_ELEMS_PER_SUBSHELL,   &
                                        NUM_BLOCK_THETA_ROWS,       &
                                        NUM_BLOCK_PHI_COLUMNS,      &
                                        NUM_R_ELEMS_PER_SHELL,      &
                                        CUR_ITERATION

USE Poseidon_Variables_Module, &
                                ONLY :  NUM_R_ELEMENTS,             &
                                        NUM_T_ELEMENTS,             &
                                        NUM_P_ELEMENTS,             &
                                        NUM_TP_QUAD_POINTS,         &
                                        rlocs, tlocs, plocs,        &
                                        NUM_R_NODES,                &
                                        INT_R_LOCATIONS,            &
                                        INT_T_LOCATIONS,            &
                                        INT_P_LOCATIONS,            &
                                        INT_R_WEIGHTS,              &
                                        INT_T_WEIGHTS,              &
                                        INT_P_WEIGHTS,              &
                                        VAR_DIM,                    &
                                        Block_Source_E,             &
                                        Block_Source_S,             &
                                        Block_Source_Si,            &
                                        Ylm_Table_Block,            &
                                        Ylm_Values,                 &
                                        Ylm_dt_Values,              &
                                        Ylm_dp_Values,              &
                                        Ylm_CC_Values,              &
                                        Ylm_CC_dt_Values,           &
                                        Ylm_CC_dp_Values,           &
                                        Lagrange_Poly_Table,        &
                                        LPT_LPT,                    &
                                        Coefficient_Vector,         &
                                        nPROCS_SHELL,               &
                                        myID_Poseidon,              &
                                        myID_Shell,                 &
                                        myID_SubShell,              &
                                        myID_PETSc,                 &
                                        INNER_CFA_BC_VALUES,        &
                                        OUTER_CFA_BC_VALUES,        &
                                        INNER_CFA_BC_TYPE,          &
                                        OUTER_CFA_BC_TYPE,          &
                                        PROB_DIM,                   &
                                        Block_PROB_DIM,             &
                                        Coefficient_Vector,         &
                                        LM_LENGTH,                  &
                                        M_VALUES,                   &
                                        ULM_LENGTH,                 &
                                        Local_Length,               &
                                        POSEIDON_COMM_WORLD,        &
                                        POSEIDON_COMM_SHELL,        &
                                        POSEIDON_COMM_PETSC,        &
                                        Block_STF_MAT,              &
                                        ELEM_PROB_DIM,              &
                                        ELEM_PROB_DIM_SQR,          &
                                        SUBSHELL_PROB_DIM,          &
                                        BLOCK_ELEM_STF_MATVEC,      &
                                        myShell,                    &
                                        Block_RHS_Vector,           &
                                        Matrix_Location,            &
                                        LM_Location


USE DRIVER_Parameters,  &
                    ONLY :  DRIVER_FRAME

USE Poseidon_IO_Module, &
                    ONLY :  OPEN_NEW_FILE

USE Poseidon_Mapping_Functions_Module, &
                                ONLY :  Map_To_X_Space,             &
                                        CFA_ALL_Matrix_Map

USE Jacobian_Internal_Functions_Module, &
                                ONLY :  JCBN_kappa_FUNCTION_3D_ALL,    &
                                        JCBN_BIGK_FUNCTION

IMPLICIT NONE



REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CUR_R_LOCS
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CUR_T_LOCS
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CUR_P_LOCS


REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: R_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: R_CUBED

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: SIN_VAL
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: SIN_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: TP_SIN_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: COS_VAL
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: COS_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CSC_VAL
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: CSC_SQUARE
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)             :: COTAN_VAL

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)           :: RSIN_SQUARE

COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)          :: PHI_EXP
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)          :: PHI_TWOEXP


REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)           :: CUR_VAL_PSI
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)           :: CUR_VAL_ALPHAPSI
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: CUR_VAL_BETA

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: CUR_DRV_PSI
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)         :: CUR_DRV_ALPHAPSI
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:,:)       :: CUR_DRV_BETA
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:)      :: CUR_LAPLACIAN_VAL



CONTAINS

!+102+###########################################################################!
!                                                                                !
!           Output_Residual                                                   !
!                                                                                !
!################################################################################!
SUBROUTINE Output_Residual()

CHARACTER(LEN = 65), DIMENSION(:), ALLOCATABLE              ::  Filenames
INTEGER, DIMENSION(:), ALLOCATABLE                          ::  File_IDs
INTEGER                                                     ::  Num_Files
INTEGER                                                     ::  Num_CFA_Eqs
INTEGER                                                     ::  Num_Shift_Eqs


REAL(KIND = idp)                                            ::  Delta_Theta,        &
                                                                Delta_Phi

REAL(KIND = idp), DIMENSION(:,:,:), ALLOCATABLE             ::  Residual_Holder
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  R_Holder,           &
                                                                T_Holder,           &
                                                                P_Holder,           &
                                                                Output_rc,          &
                                                                Output_re,          &
                                                                Output_dr

INTEGER                                                     ::  re,           &
                                                                te,           &
                                                                pe,           &
                                                                Global_re,          &
                                                                Global_te,          &
                                                                Global_pe

INTEGER                                                     ::  Block_T_Begin,      &
                                                                Block_P_Begin

REAL(KIND = idp)                                            ::  TWOOVER_DELTAR,     &
                                                                deltar_overtwo,     &
                                                                deltat_overtwo,     &
                                                                deltap_overtwo

INTEGER                                                     ::  NUM_RADIAL_SAMPLES, &
                                                                NUM_THETA_RAYS,     &
                                                                NUM_PHI_RAYS

INTEGER                                                     ::  rd, td, pd, tpd, i


REAL(KIND = idp), DIMENSION(1:11)                           ::  PSI_POWER
REAL(KIND = idp), DIMENSION(1:4)                            ::  ALPHAPSI_POWER
REAL(KIND = idp)                                            ::  JCBN_BIGK_VALUE
REAL(KIND = idp), DIMENSION(1:3)                            ::  JCBN_n_ARRAY
REAL(KIND = idp)                                            ::  Beta_Source_Prefix

INTEGER                                                     ::  Num_Entries
REAL(KIND = idp), DIMENSION(1:5)                            ::  RMS_VALUE


CALL Allocate_Residual_Variables()


100 FORMAT (A,I2.2,A,I2.2,A)

! Define Parameters

Num_Shift_Eqs = 1
Num_CFA_Eqs = 2 + Num_Shift_Eqs
Num_Files = 4 + Num_CFA_EQs


Block_T_Begin = MOD(myID_Shell,NUM_BLOCK_THETA_ROWS)*NUM_T_ELEMS_PER_BLOCK
Block_P_Begin = (myID_Shell/NUM_BLOCK_THETA_ROWS)*NUM_P_ELEMS_PER_BLOCK


!
!! Allocate Holders for Variables
!

ALLOCATE( Filenames(1:Num_Files) )
ALLOCATE( File_IDs(1:Num_Files) )

ALLOCATE( Residual_Holder(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:Num_CFA_eqs) )


! Open Output Files
WRITE(Filenames(1),100)"OUTPUT/Poseidon_Objects/Residual/Residual_DIM_F",DRIVER_FRAME,"_I",CUR_ITERATION,".out"
WRITE(Filenames(2),100)"OUTPUT/Poseidon_Objects/Residual/R_Values_F",DRIVER_FRAME,"_I",CUR_ITERATION,".out"
WRITE(Filenames(3),100)"OUTPUT/Poseidon_Objects/Residual/T_Values_F",DRIVER_FRAME,"_I",CUR_ITERATION,".out"
WRITE(Filenames(4),100)"OUTPUT/Poseidon_Objects/Residual/P_Values_F",DRIVER_FRAME,"_I",CUR_ITERATION,".out"
WRITE(Filenames(5),100)"OUTPUT/Poseidon_Objects/Residual/EQ1_Residual_F",DRIVER_FRAME,"_I",CUR_ITERATION,".out"
WRITE(Filenames(6),100)"OUTPUT/Poseidon_Objects/Residual/EQ2_Residual_F",DRIVER_FRAME,"_I",CUR_ITERATION,".out"
WRITE(Filenames(7),100)"OUTPUT/Poseidon_Objects/Residual/EQ3_Residual_F",DRIVER_FRAME,"_I",CUR_ITERATION,".out"
IF( Num_CFA_EQs > 3 ) THEN
    WRITE(Filenames(8),100)"OUTPUT/Poseidon_Objects/Residual/EQ4_Residual_F",DRIVER_FRAME,"_I",CUR_ITERATION,".out"
    IF ( Num_CFA_EQs > 4 ) THEN
        WRITE(Filenames(9),100)"OUTPUT/Poseidon_Objects/Residual/EQ5_Residual_F",DRIVER_FRAME,"_I",CUR_ITERATION,".out"
    END IF
END IF


File_IDs = [(181 + i, i=1,Num_Files)]
DO i = 1,Num_Files
    CALL OPEN_NEW_FILE( Filenames(i), File_IDs(i) )
END DO











WRITE(FILE_IDs(1),*)NUM_R_ELEMS_PER_BLOCK,NUM_T_ELEMS_PER_BLOCK,NUM_P_ELEMS_PER_BLOCK
WRITE(FILE_IDs(1),*)NUM_R_QUAD_POINTS, NUM_T_QUAD_POINTS, NUM_P_QUAD_POINTS


DO pe = 0, NUM_P_ELEMS_PER_BLOCK-1

    Global_pe = Block_P_Begin + pe

    deltap_overtwo = (plocs(Global_pe + 1) - plocs(Global_pe))/2.0_idp
    CUR_P_LOCS(:) = deltap_overtwo * (INT_P_LOCATIONS(:)+1.0_idp) + plocs(Global_pe)

    PHI_EXP(:) = EXP( CMPLX(0, -CUR_P_LOCS(:), KIND = idp) )
    PHI_TWOEXP(:) = EXP( CMPLX(0, -2.0_idp*CUR_P_LOCS(:), KIND = idp) )


    DO te = 0, NUM_T_ELEMS_PER_BLOCK-1

        Global_te = Block_T_Begin + te

        deltat_overtwo = (tlocs(Global_te + 1) - tlocs(Global_te))/2.0_idp
        CUR_T_LOCS(:) = deltat_overtwo * (INT_T_LOCATIONS(:)+1.0_idp) + tlocs(Global_te)


        COTAN_VAL(:) = 1.0_idp/DTAN(CUR_T_LOCS(:))
        CSC_VAL(:) = 1.0_idp/DSIN(CUR_T_LOCS(:))
        SIN_VAL(:) = DSIN(CUR_T_LOCS(:))
        COS_VAL(:) = DCOS(CUR_T_LOCS(:))

        COS_SQUARE(:) = COS_VAL(:)*COS_VAL(:)
        SIN_SQUARE(:) = SIN_VAL(:)*SIN_VAL(:)
        CSC_SQUARE(:) = CSC_VAL(:)*CSC_VAL(:)

        DO td = 1,NUM_T_QUAD_POINTS
            DO pd = 1,NUM_P_QUAD_POINTS
                tpd = (td-1)*NUM_P_QUAD_POINTS + pd
                TP_SIN_SQUARE(tpd) = SIN_SQUARE(td)
            END DO
        END DO

        !
        !    Move through radius
        !
        DO re = 0,NUM_R_ELEMS_PER_BLOCK-1

            Global_re = myShell*NUM_R_ELEMS_PER_SHELL + re


            deltar_overtwo = (rlocs(Global_re + 1) - rlocs(Global_re))/2.0_idp
            TWOOVER_DELTAR = 1.0_idp/deltar_overtwo
            CUR_R_LOCS(:) = deltar_overtwo * (INT_R_LOCATIONS(:)+1.0_idp) + rlocs(Global_re)


            R_SQUARE(:) = CUR_R_LOCS(:)*CUR_R_LOCS(:)
            R_CUBED(:) = R_SQUARE(:)*CUR_R_LOCS(:)

            DO rd = 1,NUM_R_QUAD_POINTS
                RSIN_SQUARE(:,rd) = R_SQUARE(rd)*SIN_SQUARE(:)
            END DO




            CALL Calc_3D_Current_Values(Global_re , te  , pe,               &
                                        DELTAR_OVERTWO, DELTAT_OVERTWO, DELTAP_OVERTWO  )




            CALL CALC_LAPLACIAN_VALUES(re, te, pe, NUM_CFA_EQs, TWOOVER_DELTAR)




            DO rd = 1,NUM_R_QUAD_POINTS
                DO td = 1,NUM_T_QUAD_POINTS
                    DO pd = 1,NUM_P_QUAD_POINTS

                        tpd = (td-1)*NUM_P_QUAD_POINTS + pd


                        JCBN_BIGK_VALUE = JCBN_BIGK_FUNCTION( rd, tpd,                                &
                                      CUR_VAL_BETA, CUR_DRV_BETA,                                     &
                                      CUR_R_LOCS(rd), R_SQUARE(rd), SIN_SQUARE(td), CSC_SQUARE(td),   &
                                      RSIN_SQUARE(td, rd), COTAN_VAL(td)                              )


                        PSI_POWER(1) = CUR_VAL_PSI( tpd, rd)
                        DO i = 2,11
                            PSI_POWER(i) = PSI_POWER(i-1)*PSI_POWER(1)
                        END DO

                        ALPHAPSI_POWER(1) = CUR_VAL_ALPHAPSI( tpd, rd)
                        DO i = 2,4
                            ALPHAPSI_POWER(i) = ALPHAPSI_POWER(i-1)*ALPHAPSI_POWER(1)
                        END DO
                        Beta_Source_Prefix = 16.0_idp * pi * ALPHAPSI_POWER(1) * PSI_POWER(3) * GR_Source_Scalar

                        JCBN_n_ARRAY(:) = CUR_DRV_ALPHAPSI( tpd, rd, : ) / ALPHAPSI_POWER(1)   &
                                            - 7 * CUR_DRV_PSI( tpd, rd, : )/ PSI_POWER(1)

                        Residual_Holder(tpd, rd, 1) = CUR_LAPLACIAN_VAL(tpd,rd,1)                               &
                                                  + twopi                                                       &
                                                  * GR_Source_Scalar                                            &
                                                  * PSI_POWER(5)                                                &
                                                  * Block_Source_E(rd, td, pd, re, te, pe)    &
                                                  + PSI_POWER(7)/ (16.0_idp * ALPHAPSI_POWER(2)) * JCBN_BIGK_VALUE

                        Residual_Holder(tpd, rd, 2) = CUR_LAPLACIAN_VAL(tpd,rd,2)                                   &
                                                    - twopi                                                         &
                                                    * GR_Source_Scalar                                              &
                                                    * ALPHAPSI_POWER(1)                                             &
                                                    * PSI_POWER(4)                                                  &
                                                    * (Block_Source_E(rd, td, pd, re, te, pe)                       &
                                                        + 2.0_idp * Block_Source_S(rd, td, pd, re, te, pe)  )       &
                                                    - 7.0_idp*PSI_POWER(6)/ (16.0_idp * ALPHAPSI_POWER(1)) * JCBN_BIGK_VALUE

                        Residual_Holder(tpd, rd, 3) = CUR_LAPLACIAN_VAL(tpd,rd,2)                                               &
                                                    - Beta_Source_Prefix * Block_Source_Si(rd, td, pd, re, te, pe, 1)           &
                                                    + ( 8.0_idp/(3.0_idp * R_SQUARE(rd) )                                       &
                                                        - 4.0_idp* JCBN_n_ARRAY(1)/(3.0_idp * CUR_R_LOCS(rd)) )                 &
                                                    * CUR_VAL_BETA(tpd, rd, 1)                                                  &
                                                    + ( 2.0_idp * COTAN_VAL(td)/CUR_R_LOCS(rd)                                  &
                                                      - TwoThirds * JCBN_n_ARRAY(1)*COTAN_VAL(td)  )                            &
                                                    * CUR_VAL_BETA(tpd, rd, 2 )                                                 &
                                                    + ( FourThirds * JCBN_n_ARRAY(1) )                                          &
                                                    * CUR_DRV_BETA(tpd, rd, 1, 1 )                                              &
                                                    + ( JCBN_n_ARRAY(2)/R_SQUARE(rd) )                                          &
                                                    * CUR_DRV_BETA(tpd, rd, 2, 1 )                                              &
                                                    + ( JCBN_n_ARRAY(3)/RSIN_SQUARE(td,rd) )                                    &
                                                    * CUR_DRV_BETA(tpd, rd, 3, 1 )                                              &
                                                    + ( JCBN_n_ARRAY(2) - COTAN_VAL(td)/3.0_idp )                               &
                                                    * CUR_DRV_BETA(tpd, rd, 1, 2 )                                              &
                                                    + ( 2.0_idp * FourThirds / CUR_R_LOCS(rd) - TwoThirds * JCBN_n_ARRAY(1) )   &
                                                    * CUR_DRV_BETA(tpd, rd, 2, 2 )                                              &
                                                    + ( JCBN_n_ARRAY(3)  )                                                      &
                                                    * CUR_DRV_BETA(tpd, rd, 1, 3 )                                              &
                                                    + ( 2.0_idp * FourThirds / CUR_R_LOCS(rd) - TwoThirds * JCBN_n_Array(1) )   &
                                                    * CUR_DRV_BETA(tpd, rd, 3, 3 )


                        RMS_VALUE(1:3) = RMS_VALUE(1:3) + ABS(Residual_Holder(tpd, rd, 1:3))**2


!                        PRINT*,RE,rd,REAL(CUR_LAPLACIAN_VAL(tpd,rd,1),KIND = idp),                          &
!                              +twopi                                                                        &
!                              * GR_Source_Scalar                                                            &
!                              * PSI_POWER(5)                                                                &
!                              * Block_Source_E(rd, td, pd, re, te, pe)                    &
!                              + PSI_POWER(7)/ (16.0_idp * ALPHAPSI_POWER(2)) * JCBN_BIGK_VALUE,             &
!                             Residual_Holder(tpd,rd,1)

!                        PRINT*,REAL(CUR_LAPLACIAN_VAL(tpd,rd,2), KIND = idp),                   &
!                               - twopi                                                          &
!                                * GR_Source_Scalar                                              &
!                                * ALPHAPSI_POWER(1)                                             &
!                                * PSI_POWER(4)                                                  &
!                                * (Block_Source_E(rd, td, pd, re, te, pe)                       &
!                                    + 2.0_idp * Block_Source_S(rd, td, pd, re, te, pe)  )       &
!                                - 7.0_idp*PSI_POWER(6)/ (16.0_idp * ALPHAPSI_POWER(1)) * JCBN_BIGK_VALUE,   &
!                                Residual_Holder(tpd,rd,1)


!                        PRINT*,Residual_Holder(tpd,rd,:)

                    END DO ! pd Loop
                END DO ! td Loop
            END DO ! rd Loop

            
            WRITE(File_IDs(4),*) CUR_P_LOCS
            WRITE(File_IDs(5),*) Residual_Holder(:,:,1)
            WRITE(File_IDs(6),*) Residual_Holder(:,:,2)
            WRITE(File_IDs(7),*) Residual_Holder(:,:,3)
    

        END DO  ! pe Loop
        WRITE(File_IDs(3),*) CUR_T_LOCS
    END DO  ! te Loop
    WRITE(File_IDs(2),*) CUR_R_LOCS
END DO  ! re Loop


DO i = 1,Num_Files
    CLOSE( FILE_IDs(i) )
END DO


NUM_ENTRIES = NUM_R_ELEMS_PER_BLOCK         &
            * NUM_T_ELEMS_PER_BLOCK         &
            * NUM_P_ELEMS_PER_BLOCK         &
            * NUM_R_QUAD_POINTS             &
            * NUM_T_QUAD_POINTS             &
            * NUM_P_QUAD_POINTS

PRINT*,"Residual RMS ",SQRT( RMS_VALUE(1:3)/NUM_ENTRIES )," Euclidean ",SQRT( RMS_Value(1:3) )




CALL Deallocate_Residual_Variables()


END SUBROUTINE Output_Residual








!+102+###########################################################################!
!                                                                                !
!           Output_Residual                                                   !
!                                                                                !
!################################################################################!
SUBROUTINE Output_Final_Residual()

CHARACTER(LEN = 65), DIMENSION(:), ALLOCATABLE              ::  Filenames
INTEGER, DIMENSION(:), ALLOCATABLE                          ::  File_IDs
INTEGER                                                     ::  Num_Files
INTEGER                                                     ::  Num_CFA_Eqs
INTEGER                                                     ::  Num_Shift_Eqs


REAL(KIND = idp)                                            ::  Delta_Theta,        &
                                                                Delta_Phi

REAL(KIND = idp), DIMENSION(:,:,:), ALLOCATABLE             ::  Residual_Holder
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  R_Holder,           &
                                                                T_Holder,           &
                                                                P_Holder,           &
                                                                Output_rc,          &
                                                                Output_re,          &
                                                                Output_dr

INTEGER                                                     ::  re,           &
                                                                te,           &
                                                                pe,           &
                                                                Global_re,          &
                                                                Global_te,          &
                                                                Global_pe

INTEGER                                                     ::  Block_T_Begin,      &
                                                                Block_P_Begin

REAL(KIND = idp)                                            ::  TWOOVER_DELTAR,     &
                                                                deltar_overtwo,     &
                                                                deltat_overtwo,     &
                                                                deltap_overtwo

INTEGER                                                     ::  NUM_RADIAL_SAMPLES, &
                                                                NUM_THETA_RAYS,     &
                                                                NUM_PHI_RAYS

INTEGER                                                     ::  rd, td, pd, tpd, i


REAL(KIND = idp), DIMENSION(1:11)                           ::  PSI_POWER
REAL(KIND = idp), DIMENSION(1:4)                            ::  ALPHAPSI_POWER
REAL(KIND = idp)                                            ::  JCBN_BIGK_VALUE
REAL(KIND = idp), DIMENSION(1:3)                            ::  JCBN_n_ARRAY
REAL(KIND = idp)                                            ::  Beta_Source_Prefix

INTEGER                                                     ::  Num_Entries
REAL(KIND = idp), DIMENSION(1:5)                            ::  RMS_VALUE


CALL Allocate_Residual_Variables()


100 FORMAT (A,I4.4,A,I2.2,A)

! Define Parameters

Num_Shift_Eqs = 1
Num_CFA_Eqs = 2 + Num_Shift_Eqs
Num_Files = 4 + Num_CFA_EQs + 1


Block_T_Begin = MOD(myID_Shell,NUM_BLOCK_THETA_ROWS)*NUM_T_ELEMS_PER_BLOCK
Block_P_Begin = (myID_Shell/NUM_BLOCK_THETA_ROWS)*NUM_P_ELEMS_PER_BLOCK


!
!! Allocate Holders for Variables
!

ALLOCATE( Filenames(1:Num_Files) )
ALLOCATE( File_IDs(1:Num_Files) )

ALLOCATE( Residual_Holder(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:Num_CFA_eqs) )


! Open Output Files
WRITE(Filenames(1),100)"OUTPUT/Poseidon_Objects/Residual/Residual_DIM_RE",NUM_R_ELEMENTS,"_D",DEGREE,".out"
WRITE(Filenames(2),100)"OUTPUT/Poseidon_Objects/Residual/R_Values_RE",NUM_R_ELEMENTS,"_D",DEGREE,".out"
WRITE(Filenames(3),100)"OUTPUT/Poseidon_Objects/Residual/T_Values_RE",NUM_R_ELEMENTS,"_D",DEGREE,".out"
WRITE(Filenames(4),100)"OUTPUT/Poseidon_Objects/Residual/P_Values_RE",NUM_R_ELEMENTS,"_D",DEGREE,".out"
WRITE(Filenames(5),100)"OUTPUT/Poseidon_Objects/Residual/EQ1_Residual_RE",NUM_R_ELEMENTS,"_D",DEGREE,".out"
WRITE(Filenames(6),100)"OUTPUT/Poseidon_Objects/Residual/EQ2_Residual_RE",NUM_R_ELEMENTS,"_D",DEGREE,".out"
WRITE(Filenames(7),100)"OUTPUT/Poseidon_Objects/Residual/EQ3_Residual_RE",NUM_R_ELEMENTS,"_D",DEGREE,".out"
IF( Num_CFA_EQs > 3 ) THEN
    WRITE(Filenames(8),100)"OUTPUT/Poseidon_Objects/Residual/EQ4_Residual_RE",NUM_R_ELEMENTS,"_D",DEGREE,".out"
    IF ( Num_CFA_EQs > 4 ) THEN
        WRITE(Filenames(9),100)"OUTPUT/Poseidon_Objects/Residual/EQ5_Residual_RE",NUM_R_ELEMENTS,"_D",DEGREE,".out"
    END IF
END IF
WRITE(Filenames(Num_Files),100)"OUTPUT/Poseidon_Objects/Residual/EQ1_Res_Par_RE",NUM_R_ELEMENTS,"_D",DEGREE,".out"

File_IDs = [(181 + i, i=1,Num_Files)]
DO i = 1,Num_Files
    CALL OPEN_NEW_FILE( Filenames(i), File_IDs(i) )
END DO











WRITE(FILE_IDs(1),*)NUM_R_ELEMS_PER_BLOCK,NUM_T_ELEMS_PER_BLOCK,NUM_P_ELEMS_PER_BLOCK
WRITE(FILE_IDs(1),*)NUM_R_QUAD_POINTS, NUM_T_QUAD_POINTS, NUM_P_QUAD_POINTS


DO pe = 0, NUM_P_ELEMS_PER_BLOCK-1

    Global_pe = Block_P_Begin + pe

    deltap_overtwo = (plocs(Global_pe + 1) - plocs(Global_pe))/2.0_idp
    CUR_P_LOCS(:) = deltap_overtwo * (INT_P_LOCATIONS(:)+1.0_idp) + plocs(Global_pe)

    PHI_EXP(:) = EXP( CMPLX(0, -CUR_P_LOCS(:), KIND = idp) )
    PHI_TWOEXP(:) = EXP( CMPLX(0, -2.0_idp*CUR_P_LOCS(:), KIND = idp) )


    DO te = 0, NUM_T_ELEMS_PER_BLOCK-1

        Global_te = Block_T_Begin + te

        deltat_overtwo = (tlocs(Global_te + 1) - tlocs(Global_te))/2.0_idp
        CUR_T_LOCS(:) = deltat_overtwo * (INT_T_LOCATIONS(:)+1.0_idp) + tlocs(Global_te)


        COTAN_VAL(:) = 1.0_idp/DTAN(CUR_T_LOCS(:))
        CSC_VAL(:) = 1.0_idp/DSIN(CUR_T_LOCS(:))
        SIN_VAL(:) = DSIN(CUR_T_LOCS(:))
        COS_VAL(:) = DCOS(CUR_T_LOCS(:))

        COS_SQUARE(:) = COS_VAL(:)*COS_VAL(:)
        SIN_SQUARE(:) = SIN_VAL(:)*SIN_VAL(:)
        CSC_SQUARE(:) = CSC_VAL(:)*CSC_VAL(:)

        DO td = 1,NUM_T_QUAD_POINTS
            DO pd = 1,NUM_P_QUAD_POINTS
                tpd = (td-1)*NUM_P_QUAD_POINTS + pd
                TP_SIN_SQUARE(tpd) = SIN_SQUARE(td)
            END DO
        END DO

        !
        !    Move through radius
        !
        DO re = 0,NUM_R_ELEMS_PER_BLOCK-1

            Global_re = myShell*NUM_R_ELEMS_PER_SHELL + re


            deltar_overtwo = (rlocs(Global_re + 1) - rlocs(Global_re))/2.0_idp
            TWOOVER_DELTAR = 1.0_idp/deltar_overtwo
            CUR_R_LOCS(:) = deltar_overtwo * (INT_R_LOCATIONS(:)+1.0_idp) + rlocs(Global_re)


            R_SQUARE(:) = CUR_R_LOCS(:)*CUR_R_LOCS(:)
            R_CUBED(:) = R_SQUARE(:)*CUR_R_LOCS(:)

            DO rd = 1,NUM_R_QUAD_POINTS
                RSIN_SQUARE(:,rd) = R_SQUARE(rd)*SIN_SQUARE(:)
            END DO




            CALL Calc_3D_Current_Values(Global_re , te  , pe,               &
                                        DELTAR_OVERTWO, DELTAT_OVERTWO, DELTAP_OVERTWO  )




            CALL CALC_LAPLACIAN_VALUES(re, te, pe, NUM_CFA_EQs, TWOOVER_DELTAR)




            DO rd = 1,NUM_R_QUAD_POINTS
                DO td = 1,NUM_T_QUAD_POINTS
                    DO pd = 1,NUM_P_QUAD_POINTS

                        tpd = (td-1)*NUM_P_QUAD_POINTS + pd


                        JCBN_BIGK_VALUE = JCBN_BIGK_FUNCTION( rd, tpd,                                &
                                      CUR_VAL_BETA, CUR_DRV_BETA,                                     &
                                      CUR_R_LOCS(rd), R_SQUARE(rd), SIN_SQUARE(td), CSC_SQUARE(td),   &
                                      RSIN_SQUARE(td, rd), COTAN_VAL(td)                              )


                        PSI_POWER(1) = CUR_VAL_PSI( tpd, rd)
                        DO i = 2,11
                            PSI_POWER(i) = PSI_POWER(i-1)*PSI_POWER(1)
                        END DO

                        ALPHAPSI_POWER(1) = CUR_VAL_ALPHAPSI( tpd, rd)
                        DO i = 2,4
                            ALPHAPSI_POWER(i) = ALPHAPSI_POWER(i-1)*ALPHAPSI_POWER(1)
                        END DO
                        Beta_Source_Prefix = 16.0_idp * pi * ALPHAPSI_POWER(1) * PSI_POWER(3) * GR_Source_Scalar

                        JCBN_n_ARRAY(:) = CUR_DRV_ALPHAPSI( tpd, rd, : ) / ALPHAPSI_POWER(1)   &
                                            - 7 * CUR_DRV_PSI( tpd, rd, : )/ PSI_POWER(1)

                        Residual_Holder(tpd, rd, 1) = CUR_LAPLACIAN_VAL(tpd,rd,1)                               &
                                                  + twopi                                                       &
                                                  * GR_Source_Scalar                                            &
                                                  * PSI_POWER(5)                                                &
                                                  * Block_Source_E(rd, td, pd, re, te, pe)    &
                                                  + PSI_POWER(7)/ (16.0_idp * ALPHAPSI_POWER(2)) * JCBN_BIGK_VALUE

                        Residual_Holder(tpd, rd, 2) = CUR_LAPLACIAN_VAL(tpd,rd,2)                                   &
                                                    - twopi                                                         &
                                                    * GR_Source_Scalar                                              &
                                                    * ALPHAPSI_POWER(1)                                             &
                                                    * PSI_POWER(4)                                                  &
                                                    * (Block_Source_E(rd, td, pd, re, te, pe)                       &
                                                        + 2.0_idp * Block_Source_S(rd, td, pd, re, te, pe)  )       &
                                                    - 7.0_idp*PSI_POWER(6)/ (16.0_idp * ALPHAPSI_POWER(1)) * JCBN_BIGK_VALUE

                        Residual_Holder(tpd, rd, 3) = CUR_LAPLACIAN_VAL(tpd,rd,2)                                               &
                                                    - Beta_Source_Prefix * Block_Source_Si(rd, td, pd, re, te, pe, 1)           &
                                                    + ( 8.0_idp/(3.0_idp * R_SQUARE(rd) )                                       &
                                                        - 4.0_idp* JCBN_n_ARRAY(1)/(3.0_idp * CUR_R_LOCS(rd)) )                 &
                                                    * CUR_VAL_BETA(tpd, rd, 1)                                                  &
                                                    + ( 2.0_idp * COTAN_VAL(td)/CUR_R_LOCS(rd)                                  &
                                                      - TwoThirds * JCBN_n_ARRAY(1)*COTAN_VAL(td)  )                            &
                                                    * CUR_VAL_BETA(tpd, rd, 2 )                                                 &
                                                    + ( FourThirds * JCBN_n_ARRAY(1) )                                          &
                                                    * CUR_DRV_BETA(tpd, rd, 1, 1 )                                              &
                                                    + ( JCBN_n_ARRAY(2)/R_SQUARE(rd) )                                          &
                                                    * CUR_DRV_BETA(tpd, rd, 2, 1 )                                              &
                                                    + ( JCBN_n_ARRAY(3)/RSIN_SQUARE(td,rd) )                                    &
                                                    * CUR_DRV_BETA(tpd, rd, 3, 1 )                                              &
                                                    + ( JCBN_n_ARRAY(2) - COTAN_VAL(td)/3.0_idp )                               &
                                                    * CUR_DRV_BETA(tpd, rd, 1, 2 )                                              &
                                                    + ( 2.0_idp * FourThirds / CUR_R_LOCS(rd) - TwoThirds * JCBN_n_ARRAY(1) )   &
                                                    * CUR_DRV_BETA(tpd, rd, 2, 2 )                                              &
                                                    + ( JCBN_n_ARRAY(3)  )                                                      &
                                                    * CUR_DRV_BETA(tpd, rd, 1, 3 )                                              &
                                                    + ( 2.0_idp * FourThirds / CUR_R_LOCS(rd) - TwoThirds * JCBN_n_Array(1) )   &
                                                    * CUR_DRV_BETA(tpd, rd, 3, 3 )


                        RMS_VALUE(1:3) = RMS_VALUE(1:3) + ABS(Residual_Holder(tpd, rd, 1:3))**2




!                        PRINT*,REAL(CUR_LAPLACIAN_VAL(tpd,rd,2), KIND = idp),                   &
!                               - twopi                                                          &
!                                * GR_Source_Scalar                                              &
!                                * ALPHAPSI_POWER(1)                                             &
!                                * PSI_POWER(4)                                                  &
!                                * (Block_Source_E(rd, td, pd, re, te, pe)                       &
!                                    + 2.0_idp * Block_Source_S(rd, td, pd, re, te, pe)  )       &
!                                - 7.0_idp*PSI_POWER(6)/ (16.0_idp * ALPHAPSI_POWER(1)) * JCBN_BIGK_VALUE,   &
!                                Residual_Holder(tpd,rd,1)


!                        PRINT*,Residual_Holder(tpd,rd,:)

                    END DO ! pd Loop
                END DO ! td Loop

                WRITE(File_IDs(Num_Files),*)RE,rd,REAL(CUR_LAPLACIAN_VAL(1,rd,1),KIND = idp),                          &
                  twopi                                                                        &
                  * GR_Source_Scalar                                                            &
                  * PSI_POWER(5)                                                                &
                  * Block_Source_E(rd, 1, 1, re, te, pe)                    &
                  + PSI_POWER(7)/ (16.0_idp * ALPHAPSI_POWER(2)) * JCBN_BIGK_VALUE,             &
                 Residual_Holder(1,rd,1)

            END DO ! rd Loop



            
            WRITE(File_IDs(2),*) CUR_R_LOCS
            WRITE(File_IDs(5),*) Residual_Holder(:,:,1)
            WRITE(File_IDs(6),*) Residual_Holder(:,:,2)
            WRITE(File_IDs(7),*) Residual_Holder(:,:,3)
    

        END DO  ! pe Loop
        WRITE(File_IDs(3),*) CUR_T_LOCS
    END DO  ! te Loop
    WRITE(File_IDs(4),*) CUR_P_LOCS
END DO  ! pe Loop


DO i = 1,Num_Files
    CLOSE( FILE_IDs(i) )
END DO


NUM_ENTRIES = NUM_R_ELEMS_PER_BLOCK         &
            * NUM_T_ELEMS_PER_BLOCK         &
            * NUM_P_ELEMS_PER_BLOCK         &
            * NUM_R_QUAD_POINTS             &
            * NUM_T_QUAD_POINTS             &
            * NUM_P_QUAD_POINTS

PRINT*,"Residual RMS ",SQRT( RMS_VALUE(1:3)/NUM_ENTRIES )," Euclidean ",SQRT( RMS_Value(1:3) )




CALL Deallocate_Residual_Variables()


END SUBROUTINE Output_Final_Residual








!+202+###########################################################################!
!                                                                                !
!                  Calc_Current_Values          !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Laplacian_Values( re, te, pe,               &
                                    NUM_CFA_EQs,            &
                                    TWOOVER_DELTAR          )

INTEGER, INTENT(IN)                                         ::  re, te, pe
INTEGER, INTENT(IN)                                         ::  NUM_CFA_EQs

REAL(KIND = idp ), INTENT(IN)                               ::  TWOOVER_DELTAR

REAL(KIND = idp ), DIMENSION(1:NUM_R_QUAD_POINTS)           ::  TMP_VAL



INTEGER                                                     :: u, d, l, m, lm_loc
INTEGER                                                     :: current_location
INTEGER                                                     :: rd


CUR_LAPLACIAN_VAL = 0.0_idp
DO u = 1,NUM_CFA_EQs
    DO d = 0,DEGREE
        DO l = 0,L_Limit
            DO m = -M_VALUES(l),M_VALUES(l)

                lm_loc = l*(l+1) + m


                Current_Location =  Matrix_Location( u, l, m, re, d )



                TMP_VAL(:) = Lagrange_Poly_Table(d, :, 2 )          &
                                * TWOOVER_DELTAR                    &
                                * TWOOVER_DELTAR                    &
                            + 2.0_idp/CUR_R_LOCS(:)                 &
                                * Lagrange_Poly_Table(d, :, 1 )     &
                                * TWOOVER_DELTAR                    &
                            - REAL( l*(l+1), KIND = idp )           &
                                /R_SQUARE(:)                        &
                                * Lagrange_Poly_Table(d, :, 0 )
            
 


                DO rd = 1,NUM_R_QUAD_POINTS
                    CUR_LAPLACIAN_VAL(:,rd,u) = CUR_LAPLACIAN_VAL(:,rd,u)   &
                                    + Coefficient_Vector(Current_Location)  &
                                    * TMP_VAL(rd)                           &
                                    * Ylm_Values( lm_loc, :, te, pe )
                END DO ! rd Loop

            END DO ! m Loop
        END DO  ! l Loop
    END DO ! d Loop
END DO  ! u Loop








END SUBROUTINE Calc_Laplacian_Values








SUBROUTINE Calc_3D_Current_Values( re, te, pe,                                  &
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



INTEGER                                                         ::  l, m, d, dp,        &
                                                                    rd, td, pd, tpd,    &
                                                                    ui


INTEGER                                                         ::  Here, There
INTEGER                                                         ::  lm_loc, lpmp_loc

COMPLEX(KIND = idp), DIMENSION(1:5)                             ::  Local_Coefficients




!$OMP PARALLEL DEFAULT(none)                                            &
!$OMP PRIVATE( tpd, rd, l, m, d, ui,                                    &
!$OMP           TMP_U_Value, Tmp_U_R_DRV_Value, Tmp_U_T_DRV_Value,      &
!$OMP           Tmp_U_P_DRV_Value,                                      &
!$OMP           local_coefficients,                                     &
!$OMP           Here, There,                                            &
!$OMP           lm_loc                                              )   &
!$OMP SHARED( re, te, pe,                                               &
!$OMP           DEGREE,                                                 &
!$OMP           NUM_TP_QUAD_POINTS, NUM_R_QUAD_POINTS,                  &
!$OMP           CUR_VAL_PSI, CUR_VAL_ALPHAPSI, CUR_VAL_BETA,            &
!$OMP           CUR_DRV_PSI, CUR_DRV_ALPHAPSI, CUR_DRV_BETA,            &
!$OMP           CUR_R_LOCS, CUR_T_LOCS, CUR_P_LOCS,                     &
!$OMP           R_SQUARE,                                               &
!$OMP           SIN_VAL,                                                &
!$OMP           myID_Poseidon,                                          &
!$OMP           DELTAR_OVERTWO, DELTAT_OVERTWO, DELTAP_OVERTWO,         &
!$OMP           Lagrange_Poly_Table,                                    &
!$OMP           Ylm_Values, Ylm_dt_Values, Ylm_dp_Values,               &
!$OMP           Coefficient_Vector,                                     &
!$OMP           Matrix_Location,                                        &
!$OMP           LM_Location,                                            &
!$OMP           LM_Length, ULM_LENGTH                           )




!$OMP DO SCHEDULE(dynamic), COLLAPSE(2)
DO rd = 1,NUM_R_QUAD_POINTS

    DO tpd = 1,NUM_TP_QUAD_POINTS

         !                                                  !
        !!   Set/Reset temporary value holder to zero.      !!
         !                                                  !
        Tmp_U_Value = 0.0_idp
        Tmp_U_R_DRV_Value = 0.0_idp
        Tmp_U_T_DRV_Value = 0.0_idp
        Tmp_U_P_DRV_Value = 0.0_idp


        DO d = 0,DEGREE

            DO ui = 1,NUM_CFA_VARS

                Here = CFA_All_Matrix_Map(ui, 0, re, d)
                There = Here + LM_LENGTH - 1

!                PRINT*,Coefficient_Vector(Here:There)

                TMP_U_Value(ui)         = TMP_U_Value(ui)                           &
                                        + SUM( Coefficient_Vector( Here:There )     &
                                        * Ylm_Values( :, tpd, te, pe )       )      &
                                        * Lagrange_Poly_Table( d, rd, 0 )

                TMP_U_R_DRV_Value(ui)   = TMP_U_R_DRV_Value(ui)                     &
                                        + SUM( Coefficient_Vector( Here:There )     &
                                        * Ylm_Values( :, tpd, te, pe )       )      &
                                        * Lagrange_Poly_Table( d, rd, 1 )           &
                                        / DELTAR_OVERTWO

                TMP_U_T_DRV_Value(ui)   = TMP_U_T_DRV_Value(ui)                     &
                                        + SUM( Coefficient_Vector( Here:There )     &
                                        * Ylm_dt_Values( :, tpd, te, pe)     )       &
                                        * Lagrange_Poly_Table( d, rd, 0)

                TMP_U_P_DRV_Value(ui)   = TMP_U_P_DRV_Value(ui)                     &
                                        + SUM( Coefficient_Vector( Here:There )     &
                                        * Ylm_dp_Values( :, tpd, te, pe)     )      &
                                        * Lagrange_Poly_Table( d, rd, 0)

            END DO ! ui Loop

        END DO  !   d Loop


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
        CUR_DRV_BETA( tpd, rd, 2, 1 )  = REAL(Tmp_U_T_DRV_Value(3), KIND = idp)
        CUR_DRV_BETA( tpd, rd, 3, 1 )  = REAL(Tmp_U_P_DRV_Value(3), KIND = idp)
        CUR_DRV_BETA( tpd, rd, 1, 2 )  = REAL(Tmp_U_R_DRV_Value(4), KIND = idp)
        CUR_DRV_BETA( tpd, rd, 2, 2 )  = REAL(Tmp_U_T_DRV_Value(4), KIND = idp)
        CUR_DRV_BETA( tpd, rd, 3, 2 )  = REAL(Tmp_U_P_DRV_Value(4), KIND = idp)
        CUR_DRV_BETA( tpd, rd, 1, 3 )  = REAL(Tmp_U_R_DRV_Value(5), KIND = idp)
        CUR_DRV_BETA( tpd, rd, 2, 3 )  = REAL(Tmp_U_T_DRV_Value(5), KIND = idp)
        CUR_DRV_BETA( tpd, rd, 3, 3 )  = REAL(Tmp_U_P_DRV_Value(5), KIND = idp)

!        Beta_DRV_Trace( tpd, rd )      = REAL(Tmp_U_R_DRV_Value(3), KIND = idp)     &
!                                       + REAL(Tmp_U_T_DRV_Value(4), KIND = idp)     &
!                                       + REAL(Tmp_U_P_DRV_Value(5), KIND = idp)
        
    END DO  !   tpd Loop

END DO  !   rd Loop
!$OMP END DO

!$OMP END PARALLEL








END SUBROUTINE Calc_3D_Current_Values






!+601+###########################################################################!
!                                                                                !
!           Allocate_Residual_Variables                                      !
!                                                                                !
!################################################################################!
SUBROUTINE Allocate_Residual_Variables()


ALLOCATE( CUR_R_LOCS(1:NUM_R_QUAD_POINTS) )
ALLOCATE( CUR_T_LOCS(1:NUM_T_QUAD_POINTS) )
ALLOCATE( CUR_P_LOCS(1:NUM_P_QUAD_POINTS) )


ALLOCATE( R_SQUARE(1:NUM_R_QUAD_POINTS) )
ALLOCATE( R_CUBED(1:NUM_R_QUAD_POINTS) )

ALLOCATE( SIN_VAL( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( SIN_SQUARE( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( TP_SIN_SQUARE( 1:NUM_TP_QUAD_POINTS ) )
ALLOCATE( COS_VAL( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( COS_SQUARE( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( CSC_VAL( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( CSC_SQUARE( 1:NUM_T_QUAD_POINTS ) )
ALLOCATE( COTAN_VAL( 1:NUM_T_QUAD_POINTS ) )

ALLOCATE( RSIN_SQUARE( 1:NUM_T_QUAD_POINTS, 1:NUM_R_QUAD_POINTS ) )

ALLOCATE( PHI_EXP( 1:NUM_P_QUAD_POINTS ) )
ALLOCATE( PHI_TWOEXP( 1:NUM_P_QUAD_POINTS ) )

ALLOCATE( CUR_VAL_PSI(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS) )
ALLOCATE( CUR_VAL_ALPHAPSI(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS) )
ALLOCATE( CUR_VAL_BETA(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:3) )

ALLOCATE( CUR_DRV_PSI(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:3) )
ALLOCATE( CUR_DRV_ALPHAPSI(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:3) )
ALLOCATE( CUR_DRV_BETA(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:3, 1:3) )


ALLOCATE( CUR_LAPLACIAN_VAL(1:NUM_TP_QUAD_POINTS, 1:NUM_R_QUAD_POINTS, 1:5) )


END SUBROUTINE Allocate_Residual_Variables









!+602+###########################################################################!
!                                                                                !
!           Deallocate_Residual_Variables                                    !
!                                                                                !
!################################################################################!
SUBROUTINE Deallocate_Residual_Variables()


DEALLOCATE( CUR_R_LOCS )
DEALLOCATE( CUR_T_LOCS )
DEALLOCATE( CUR_P_LOCS )

DEALLOCATE( R_SQUARE )
DEALLOCATE( R_CUBED )

DEALLOCATE( SIN_VAL )
DEALLOCATE( SIN_SQUARE )
DEALLOCATE( TP_SIN_SQUARE )
DEALLOCATE( COS_VAL )
DEALLOCATE( COS_SQUARE )
DEALLOCATE( CSC_VAL )
DEALLOCATE( CSC_SQUARE )
DEALLOCATE( COTAN_VAL )

DEALLOCATE( RSIN_SQUARE )

DEALLOCATE( PHI_EXP )
DEALLOCATE( PHI_TWOEXP )

DEALLOCATE( CUR_VAL_PSI )
DEALLOCATE( CUR_VAL_ALPHAPSI )
DEALLOCATE( CUR_VAL_BETA )

DEALLOCATE( CUR_DRV_PSI )
DEALLOCATE( CUR_DRV_ALPHAPSI )
DEALLOCATE( CUR_DRV_BETA )

DEALLOCATE( CUR_LAPLACIAN_VAL )

END SUBROUTINE Deallocate_Residual_Variables









END MODULE Poseidon_Residual_Equations_Module
