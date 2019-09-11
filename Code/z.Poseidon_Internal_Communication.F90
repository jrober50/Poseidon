   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_Internal_Communication_Module                                       !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   Poseidon_CFA_Block_Share                                            !##!
!##!                                                                                !##!
!##!    +201+   Poseidon_Distribute_Solution                                        !##!
!##!                                                                                !##!
!##!    +301+   Poseidon_Newton_Block_Share                                         !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!



USE Poseidon_Constants_Module, &
            ONLY :  idp, pi,speed_of_light


USE Units_Module, &
            ONLY :  Set_Units, Grav_Constant_G


USE DRIVER_Parameters,  &
            ONLY :  Analytic_Solution,      &
                    Shift_Solution,         &
                    Enclosed_Mass,          &
                    DRIVER_R_LOCS,         &
                    nPROCS,                 &
                    myID,                   &
                    myID_Theta,             &
                    myID_Phi,               &
                    Ratio_T_BNDLperBLCK,    &
                    Ratio_P_BNDLperBLCK,    &
                    Ratio_BNDLperBLCK,      &
                    NUM_ENTRIES,            &
                    SELFSIM_R_VALS,         &
                    SELFSIM_T


USE Poseidon_Parameters, &
            ONLY :  DOMAIN_DIM,             &
                    DEGREE,                 &
                    L_LIMIT,                &
                    NUM_CFA_VARS,           &
                    DATA_DIST_MODE,         &
                    NUM_R_ELEMS_PER_SHELL,  &
                    NUM_SHELLS,             &
                    NUM_BLOCKS_PER_SHELL,   &
                    NUM_BLOCK_THETA_ROWS,   &
                    NUM_BLOCK_PHI_COLUMNS,  &
                    nPROCS_POSEIDON,        &
                    STF_MAPPING_FLAG,       &
                    NUM_R_ELEMS_PER_BLOCK,  &
                    NUM_T_ELEMS_PER_BLOCK,  &
                    NUM_P_ELEMS_PER_BLOCK,  &
                    NUM_R_QUAD_POINTS,      &
                    NUM_T_QUAD_POINTS,      &
                    NUM_P_QUAD_POINTS,      &
                    R_COARSEN_FACTOR,       &
                    T_COARSEN_FACTOR,       &
                    P_COARSEN_FACTOR,       &
                    SOL_DIST_SCHEME


USE Poseidon_Variables_Module, &
            ONLY :  ierr,                                               &
                    myID_Poseidon,                                      &
                    myID_Shell,                                         &
                    myID_PETSc,                                         &
                    myID_Dist,                                          &
                    myShell,                                            &
                    FirstCall_Flag,                                     &
                    Stiffness_Matrix_Initialized_Flag,                  &
                    Matrix_Cholesky_Factorized_Flag,                    &
                    rlocs,                                              &
                    NUM_R_ELEMENTS,                                     &
                    NUM_T_ELEMENTS,                                     &
                    NUM_P_ELEMENTS,                                     &
                    R_INNER, R_OUTER,                                   &
                    Block_Source_E,                                     &
                    Block_Source_S,                                     &
                    Block_Source_Si,                                    &
                    POSEIDON_COMM_WORLD,                                &
                    POSEIDON_COMM_SHELL,                                &
                    POSEIDON_COMM_DIST,                                 &
                    Prob_Dim,                                           &
                    Block_Prob_Dim,                                     &
                    LM_Length,                                          &
                    ULM_Length,                                         &
                    Local_Length,                                       &
                    Coefficient_Vector,                                 &
                    rlocs, tlocs, plocs



USE Additional_Functions_Module, &
            ONLY :  Lagrange_Poly,                                      &
                    Spherical_Harmonic,                                 &
                    Map_To_X_Space, Map_From_X_Space,                   &
                    Initialize_LGL_Quadrature,                          &
                    Initialize_LGL_Quadrature_Locations,                &
                    Initialize_LG_Quadrature_Locations,                 &
                    MVMULT_FULL








use mpi



IMPLICIT NONE

CONTAINS


!+101+##########################################################################!
!                                                                               !
!                           Poseidon_CFA_Block_Share                            !
!                                                                               !
!###############################################################################!
SUBROUTINE Poseidon_CFA_Block_Share( Proc_id, Proc_Theta_id, Proc_Phi_id,                    &
                                     My_Source_E, My_Source_S, My_Source_Si,                 &
                                     Local_RE_Dim, Local_TE_Dim, Local_PE_Dim,               &
                                     Local_RD_Dim, Local_TD_Dim, Local_PD_Dim,               &
                                     Input_R_Quad, Input_T_Quad, Input_P_Quad,               &
                                     Left_Limit, Right_Limit,                                &
                                     R_Elements_Input, T_Elements_Input, P_Elements_Input,   &
                                     Input_Delta_R, Input_Delta_T, Input_Delta_P,            &
                                     Block_E, Block_S, Block_Si                              )




INTEGER, INTENT(IN)                                 ::  Proc_id,        &
                                                        Proc_Theta_id,  &
                                                        Proc_Phi_id


REAL(KIND = idp), INTENT(IN), DIMENSION(    1:Local_RD_Dim*Local_TD_Dim*Local_PD_Dim,             &
                                            0:Local_RE_Dim-1,           &
                                            0:Local_TE_Dim-1,           &
                                            0:Local_PE_Dim-1  )         ::  My_Source_E,    &
                                                                            My_Source_S

REAL(KIND = idp), INTENT(IN), DIMENSION(    1:Local_RD_Dim*Local_TD_Dim*Local_PD_Dim,             &
                                            0:Local_RE_Dim-1,           &
                                            0:Local_TE_Dim-1,           &
                                            0:Local_PE_Dim-1,           &
                                            1:DOMAIN_DIM            )   :: My_Source_Si





INTEGER, INTENT(IN)                                                     ::  Local_RE_Dim,   &
                                                                            Local_TE_Dim,   &
                                                                            Local_PE_Dim,   &
                                                                            Local_RD_Dim,   &
                                                                            Local_TD_Dim,   &
                                                                            Local_PD_Dim


REAL(KIND = idp), DIMENSION(1:Local_RD_Dim), INTENT(IN)                 ::  Input_R_Quad
REAL(KIND = idp), DIMENSION(1:Local_TD_Dim), INTENT(IN)                 ::  Input_T_Quad
REAL(KIND = idp), DIMENSION(1:Local_PD_Dim), INTENT(IN)                 ::  Input_P_Quad

REAL(KIND = idp), INTENT(IN)                                            ::  Left_Limit,     &
                                                                            Right_Limit

INTEGER, INTENT(IN)                                                     ::  R_Elements_Input,       &
                                                                            T_Elements_Input,       &
                                                                            P_Elements_Input

REAL(KIND = idp), DIMENSION(1:R_Elements_Input), INTENT(IN)             ::  Input_Delta_R
REAL(KIND = idp), DIMENSION(1:T_Elements_Input), INTENT(IN)             ::  Input_Delta_T
REAL(KIND = idp), DIMENSION(1:P_Elements_Input), INTENT(IN)             ::  Input_Delta_P


REAL(KIND = idp), INTENT(INOUT), DIMENSION( 1:NUM_R_QUAD_POINTS,        &
                                            1:NUM_T_QUAD_POINTS,        &
                                            1:NUM_P_QUAD_POINTS,        &
                                            0:NUM_R_ELEMS_PER_BLOCK-1,  &
                                            0:NUM_T_ELEMS_PER_BLOCK-1,  &
                                            0:NUM_P_ELEMS_PER_BLOCK-1 ) ::  Block_E

REAL(KIND = idp), INTENT(INOUT), DIMENSION( 1:NUM_R_QUAD_POINTS,        &
                                            1:NUM_T_QUAD_POINTS,        &
                                            1:NUM_P_QUAD_POINTS,        &
                                            0:NUM_R_ELEMS_PER_BLOCK-1,  &
                                            0:NUM_T_ELEMS_PER_BLOCK-1,  &
                                            0:NUM_P_ELEMS_PER_BLOCK-1 ) ::  Block_S

REAL(KIND = idp), INTENT(INOUT), DIMENSION( 1:NUM_R_QUAD_POINTS,        &
                                            1:NUM_T_QUAD_POINTS,        &
                                            1:NUM_P_QUAD_POINTS,        &
                                            0:NUM_R_ELEMS_PER_BLOCK-1,  &
                                            0:NUM_T_ELEMS_PER_BLOCK-1,  &
                                            0:NUM_P_ELEMS_PER_BLOCK-1,  &
                                            1:DOMAIN_DIM            )   ::  Block_Si




INTEGER                                         ::  Block_Num
Integer                                         ::  Shell, Block_indx
INTEGER                                         ::  Input_Shell_R_Begin,      &
                                                    Input_Shell_R_End,        &
                                                    Input_TE_Begin,      &
                                                    Input_PE_Begin


INTEGER                                         ::  Input_Shell_R_Loc,        &
                                                    Input_Shell_T_Loc,        &
                                                    Input_Shell_P_Loc

INTEGER                                         ::  i, j, k, d, n
INTEGER                                         ::  ii, jj, kk
INTEGER                                         ::  rd, td, pd
INTEGER                                         ::  jd, kd
INTEGER                                         ::  jq, kq


INTEGER                                         ::  Block_R_Loc,        &
                                                    Block_T_Loc,        &
                                                    Block_P_Loc

INTEGER                                         ::  Input_RE_Start,     &
                                                    Input_TE_Start,     &
                                                    Input_PE_Start

INTEGER                                                     ::  re, te, pe,                 &
                                                                Local_R, Local_T, Local_P,  &
                                                                Input_R, Input_T, Input_P,  &
                                                                Local_Here, Input_Here

INTEGER                                                     ::  Num_Local_DOF,              &
                                                                Num_Input_DOF,              &
                                                                Num_Coarse_DOF





REAL(KIND = idp), DIMENSION(1:NUM_R_QUAD_POINTS)            ::  Local_R_Locations
REAL(KIND = idp), DIMENSION(1:NUM_T_QUAD_POINTS)            ::  Local_T_Locations
REAL(KIND = idp), DIMENSION(1:NUM_P_QUAD_POINTS)            ::  Local_P_Locations

REAL(KIND = idp), DIMENSION(1:Local_RD_Dim)                 ::  Input_R_Locations
REAL(KIND = idp), DIMENSION(1:Local_TD_Dim)                 ::  Input_T_Locations
REAL(KIND = idp), DIMENSION(1:Local_PD_Dim)                 ::  Input_P_Locations







REAL(KIND = idp), DIMENSION(:,:), ALLOCATABLE                 ::  R_Lag_Poly_Values
REAL(KIND = idp), DIMENSION(:,:), ALLOCATABLE                 ::  T_Lag_Poly_Values
REAL(KIND = idp), DIMENSION(:,:), ALLOCATABLE                 ::  P_Lag_Poly_Values

REAL(KIND = idp), DIMENSION(:,:), ALLOCATABLE               ::  Translation_Matrix

REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Local_E_Coeffs,         &
                                                                Local_S_Coeffs

REAL(KIND = idp), DIMENSION(:,:), ALLOCATABLE               ::  Local_Si_Coeffs


REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)                 ::  Input_CR_Locations
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)                 ::  Input_CT_Locations
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)                 ::  Input_CP_Locations

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:,:)           ::  TMP_E_STORAGE,         &
                                                                TMP_S_STORAGE

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:,:,:)         ::  TMP_Si_STORAGE


REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Send_Message

REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Recv_Message

REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Delta_FREoCRE
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Delta_FTEoCTE
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Delta_FPEoCPE

REAL(KIND = idp)                                            ::  Input_Delta_x
REAL(KIND = idp)                                            ::  Delta_x_Ratio

INTEGER                                                     ::  Size_of_Send
INTEGER                                                     ::  Size_of_Send_b



INTEGER                                                     ::  Elems_Per_Block

INTEGER                                                     ::  Cur_Elem_Num,           &
                                                                Next_Elem_Num

INTEGER                                                     ::  Start_Here, End_Here
INTEGER                                                     ::  Start_Here_B, End_Here_B
INTEGER                                                     ::  Start_Here_C, End_Here_C

INTEGER                                                     ::  block_j, block_k, Recv_Addr

INTEGER                                                     ::  TheirID, TheirID_theta, TheirID_phi

INTEGER                                                     ::  Block_RE_Begin,Block_TE_Begin, Block_PE_Begin

INTEGER                                                     ::  Block_RE_Loc,                &
                                                                Block_TE_Loc,                &
                                                                Block_PE_Loc

INTEGER                                                     ::  Block_Input_RE_Begin,        &
                                                                Block_Input_TE_Begin,        &
                                                                Block_Input_PE_Begin

INTEGER                                                     ::  Global_Input_R_Loc,          &
                                                                Global_Input_T_Loc,          &
                                                                Global_Input_P_Loc

INTEGER                                                     ::  Input_RAD_Elems_Per_Shell
INTEGER                                                     ::  Block_Index

INTEGER                                                     ::  print_indx
INTEGER                                                     ::  request
INTEGER, DIMENSION(MPI_STATUS_SIZE)                         ::  status

INTEGER                                                     ::  Coarsen_Factor

INTEGER                                                     ::  Coarse_RD_Dim,  &
                                                                Coarse_TD_Dim,  &
                                                                Coarse_PD_Dim

INTEGER                                                     ::  Shift_Factor_A, &
                                                                Shift_Factor_B, &
                                                                Shift_Factor_C, &
                                                                Shift_Factor_D, &
                                                                Shift_Factor_E

INTEGER                                                     ::  P_Shift_Loc,    &
                                                                T_Shift_Loc,    &
                                                                R_Shift_Loc

INTEGER             ::  Cur_Elem_Num_B,  &
                        Local_Here_B





Num_Input_DOF = Local_RD_Dim * Local_TD_Dim * Local_PD_Dim

Num_Local_DOF = NUM_R_QUAD_POINTS * NUM_T_QUAD_POINTS * NUM_P_QUAD_POINTS

INPUT_RAD_ELEMS_PER_SHELL = NUM_R_ELEMS_PER_SHELL*R_Coarsen_Factor

Elems_Per_Block = Local_TE_Dim*Local_PE_Dim*INPUT_RAD_ELEMS_PER_SHELL
Size_of_Send = Elems_Per_Block * Num_Input_DOF
Size_of_Send_b = Elems_Per_Block * Num_Input_DOF*NUM_CFA_VARS

Coarsen_Factor = R_Coarsen_Factor*T_Coarsen_Factor*P_Coarsen_Factor

Coarse_RD_Dim = R_Coarsen_Factor * Local_RD_Dim
Coarse_TD_Dim = T_Coarsen_Factor * Local_TD_Dim
Coarse_PD_Dim = P_Coarsen_Factor * Local_PD_Dim

Num_Coarse_DOF = Coarse_RD_Dim * Coarse_TD_Dim * Coarse_PD_Dim

Shift_Factor_A = R_Coarsen_Factor*T_Coarsen_Factor*Local_RD_Dim*Local_TD_Dim*Local_PD_Dim
Shift_Factor_B = R_Coarsen_Factor*Local_RD_Dim*Local_TD_Dim
Shift_Factor_C = Local_RD_Dim
Shift_Factor_D = R_Coarsen_Factor*T_Coarsen_Factor*Local_RD_Dim*Local_TD_Dim
Shift_Factor_E = R_Coarsen_Factor*Local_RD_Dim


Delta_x_Ratio = 2.0_idp/(Right_Limit - Left_Limit)


ALLOCATE(Translation_Matrix(1:Num_Coarse_DOF, 1:Num_Local_DOF))

ALLOCATE(Local_E_Coeffs(1:Num_Local_DOF))
ALLOCATE(Local_S_Coeffs(1:Num_Local_DOF))
ALLOCATE(Local_Si_Coeffs(1:Num_Local_DOF,1:DOMAIN_DIM))


ALLOCATE( TMP_E_STORAGE(0:Num_Coarse_DOF-1,             &
                        0:NUM_R_ELEMS_PER_BLOCK-1,      &
                        0:NUM_T_ELEMS_PER_BLOCK-1,      &
                        0:NUM_P_ELEMS_PER_BLOCK-1)      )

ALLOCATE( TMP_S_STORAGE(0:Num_Coarse_DOF-1,             &
                        0:NUM_R_ELEMS_PER_BLOCK-1,      &
                        0:NUM_T_ELEMS_PER_BLOCK-1,      &
                        0:NUM_P_ELEMS_PER_BLOCK-1)      )

ALLOCATE( TMP_Si_STORAGE(0:Num_Coarse_DOF-1,            &
                         0:NUM_R_ELEMS_PER_BLOCK-1,     &
                         0:NUM_T_ELEMS_PER_BLOCK-1,     &
                         0:NUM_P_ELEMS_PER_BLOCK-1,     &
                         1:DOMAIN_DIM)                  )


ALLOCATE( Input_CR_Locations(0:Coarse_RD_Dim - 1) )
ALLOCATE( Input_CT_Locations(0:Coarse_TD_Dim - 1) )
ALLOCATE( Input_CP_Locations(0:Coarse_PD_Dim - 1) )


ALLOCATE( Send_Message(0:Size_of_Send_b+2) )
ALLOCATE( Recv_Message(0:Size_of_Send_b+2) )


ALLOCATE( Delta_FREoCRE(0:R_Coarsen_Factor) )
ALLOCATE( Delta_FTEoCTE(0:T_Coarsen_Factor) )
ALLOCATE( Delta_FPEoCPE(0:P_Coarsen_Factor) )



ALLOCATE( R_Lag_Poly_Values(1:Coarse_RD_Dim,1:NUM_R_QUAD_POINTS) )
ALLOCATE( T_Lag_Poly_Values(1:Coarse_TD_Dim,1:NUM_T_QUAD_POINTS) )
ALLOCATE( P_Lag_Poly_Values(1:Coarse_PD_Dim,1:NUM_P_QUAD_POINTS) )


Local_R_Locations = Initialize_LG_Quadrature_Locations(NUM_R_QUAD_POINTS)
Local_T_Locations = Initialize_LG_Quadrature_Locations(NUM_T_QUAD_POINTS)
Local_P_Locations = Initialize_LG_Quadrature_Locations(NUM_P_QUAD_POINTS)



Input_R_Locations = Input_R_Quad - Left_Limit
Input_T_Locations = Input_T_Quad - Left_Limit
Input_P_Locations = Input_P_Quad - Left_Limit


!$OMP PARALLEL DEFAULT(none)                                                    &
!$OMP PRIVATE(   Shell, k, j, i, kd, jd, kq, jq, d,                             &
!$OMP            kk, jj, ii,                                                    &
!$OMP            Local_R, Local_T, Local_P, Input_R, Input_T, Input_P,          &
!$OMP            Block_Index,                                                   &
!$OMP            Input_Shell_R_Begin, Input_Shell_R_End,                        &
!$OMP            Input_Shell_R_Loc, Cur_Elem_Num, Start_Here, End_Here,         &
!$OMP            Start_Here_B, End_Here_B, Start_Here_C, End_Here_C,            &
!$OMP            Local_Here, Input_Here,                                        &
!$OMP            block_j, block_k, Recv_Addr,                                   &
!$OMP            TheirID, TheirID_theta, TheirID_phi,                           &
!$OMP            Input_TE_Begin, Input_PE_Begin,                                &
!$OMP            Block_RE_Begin, Block_TE_Begin, Block_PE_Begin,                &
!$OMP            Block_Input_TE_Begin, Block_Input_PE_Begin,                    &
!$OMP            Global_Input_R_Loc, Global_Input_T_Loc, Global_Input_P_Loc,    &
!$OMP            Block_RE_Loc, Block_TE_Loc, Block_PE_Loc,                      &
!$OMP            R_Shift_Loc, T_Shift_Loc, P_Shift_Loc,                         &
!$OMP            Delta_FREoCRE, Delta_FTEoCTE, Delta_FPEoCPE,                   &
!$OMP            Input_PE_Start, Input_TE_Start, Input_RE_Start,                &
!$OMP            Input_CP_Locations, Input_CT_Locations, Input_CR_Locations,    &
!$OMP            R_Lag_Poly_Values, T_Lag_Poly_Values, P_Lag_Poly_Values,       &
!$OMP            Local_E_Coeffs, Local_S_Coeffs, Local_Si_Coeffs,               &
!$OMP            Translation_Matrix                                        )    &
!$OMP SHARED(    Local_PE_Dim, Local_TE_Dim, INPUT_RAD_ELEMS_PER_SHELL,         &
!$OMP            Local_PD_Dim, Local_TD_Dim, Local_RD_Dim,                      &
!$OMP            NUM_R_ELEMS_PER_BLOCK, NUM_T_ELEMS_PER_BLOCK, NUM_P_ELEMS_PER_BLOCK, &
!$OMP            R_COARSEN_FACTOR, T_COARSEN_FACTOR, P_COARSEN_FACTOR,          &
!$OMP            DOMAIN_DIM, NUM_SHELLS, NUM_BLOCK_THETA_ROWS, NUM_BLOCKS_PER_SHELL,  &
!$OMP            NUM_R_QUAD_POINTS, NUM_T_QUAD_POINTS, NUM_P_QUAD_POINTS,       &
!$OMP            rlocs, tlocs, plocs,                                           &
!$OMP            Input_Delta_P, Input_Delta_T, Input_Delta_R,                   &
!$OMP            Ratio_T_BNDLperBLCK, Ratio_P_BNDLperBLCK,                      &
!$OMP            Ratio_BNDLperBLCK,                                             &
!$OMP            Send_Message, Recv_Message,                                    &
!$OMP            Size_of_Send, Size_of_Send_b,                                  &
!$OMP            My_Source_E, My_Source_S, My_Source_Si,                        &
!$OMP            Num_Input_DOF, Num_Local_DOF, Num_Coarse_DOF,                  &
!$OMP            Coarse_RD_Dim, Coarse_TD_Dim, Coarse_PD_Dim,                   &
!$OMP            Shift_Factor_A, Shift_Factor_B, Shift_Factor_C,                &
!$OMP            Shift_Factor_D, Shift_Factor_E,                                &
!$OMP            Delta_x_Ratio,                                                 &
!$OMP            Local_R_Locations, Local_T_Locations, Local_P_Locations,       &
!$OMP            Input_R_Locations, Input_T_Locations, Input_P_Locations,       &
!$OMP            myID, myID_Theta, myID_Phi, myShell, myID_Shell, ierr,         &
!$OMP            request, sTatus, MPI_STATUS_IGNORE,                            &
!$OMP            TMP_E_STORAGE, TMP_S_STORAGE, TMP_Si_STORAGE,                  &
!$OMP            Block_E, Block_S, Block_Si                                     )


Delta_FREoCRE(0) = 0.0_idp
Delta_FTEoCTE(0) = 0.0_idp
Delta_FPEoCPE(0) = 0.0_idp



DO Shell = 0, NUM_SHELLS - 1


    !
    !   Package and Send Data
    !

    Input_Shell_R_Begin = Shell * INPUT_RAD_ELEMS_PER_SHELL
    Input_Shell_R_End =  (Shell + 1) * INPUT_RAD_ELEMS_PER_SHELL - 1


    !$OMP DO SCHEDULE(dynamic), COLLAPSE(3)
    DO k = 0,Local_PE_Dim - 1
        DO j = 0,Local_TE_Dim - 1
            DO i = 0,INPUT_RAD_ELEMS_PER_SHELL - 1

                Input_Shell_R_Loc = Input_Shell_R_Begin + i

                Cur_Elem_Num = (k * Local_TE_Dim + j )          &
                                * INPUT_RAD_ELEMS_PER_SHELL     &
                             + i




                Start_Here = Cur_Elem_Num*Num_Input_DOF + 3
                End_Here = (Cur_Elem_Num+1)*Num_Input_DOF+2


                Send_Message(Start_Here:End_Here) = My_Source_E(1:NUM_Input_DOF, Input_Shell_R_Loc, j, k )



                Start_Here = Cur_Elem_Num*Num_Input_DOF + 3 + Size_of_Send
                End_Here = (Cur_Elem_Num+1)*Num_Input_DOF+2 + Size_of_Send


                Send_Message(Start_Here:End_Here) = My_Source_S(1:NUM_Input_DOF, Input_Shell_R_Loc, j, k )




                DO d = 1,DOMAIN_DIM

                    Start_Here = Cur_Elem_Num*Num_Input_DOF + 3 + (d+1)*Size_of_Send
                    End_Here = (Cur_Elem_Num+1)*Num_Input_DOF+ 2 + (d+1)*Size_of_Send

                    Send_Message(Start_Here:End_Here) = My_Source_Si(1:NUM_Input_DOF, Input_Shell_R_Loc, j, k, d)

                END DO ! d Loop




            END DO ! k Loop

        END DO ! j Loop

    END DO ! i Loop
    !$OMP END DO



    !$OMP MASTER
    Send_Message(0) = REAL( myID, idp )
    Send_Message(1) = REAL( myID_Theta,idp)
    Send_Message(2) = REAL( myID_Phi, idp)

    block_j = myID_Theta/Ratio_T_BNDLperBLCK
    block_k = myID_Phi/Ratio_P_BNDLperBLCK

    Recv_Addr = block_k * NUM_BLOCK_THETA_ROWS + block_j + Shell*NUM_BLOCKS_PER_SHELL

    CALL MPI_ISend( Send_Message, Size_of_Send_b+3, MPI_DOUBLE, Recv_Addr, 6, MPI_COMM_WORLD, request, ierr)






    !
    !   If Process is responsible for a block in the current shell
    !   then receive data.
    !
    IF ( myShell == Shell ) THEN


        DO Block_index = 1,Ratio_BNDLperBLCK


            !$OMP MASTER
            CALL MPI_Recv( Recv_Message, Size_of_Send_b+3, MPI_DOUBLE, MPI_ANY_SOURCE,      &
                                6, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr                  )
            !$OMP END MASTER
            !$OMP BARRIER


            TheirID = INT(Recv_Message(0))
            TheirID_theta = INT(Recv_Message(1))
            TheirID_phi = INT(Recv_Message(2))


!            Bundle_List(Block_index) = TheirID

            Input_TE_Begin = TheirID_theta*LOCAL_TE_DIM
            Input_PE_Begin = TheirID_phi*LOCAL_PE_DIM

            Block_RE_Begin = Shell*NUM_R_ELEMS_PER_BLOCK
            Block_TE_Begin = MOD(myID_Shell,NUM_BLOCK_THETA_ROWS)*NUM_T_ELEMS_PER_BLOCK
            Block_PE_Begin = (myID_Shell/NUM_BLOCK_THETA_ROWS)*NUM_P_ELEMS_PER_BLOCK

            Block_Input_TE_Begin = Block_TE_Begin * T_Coarsen_Factor
            Block_Input_PE_Begin = Block_PE_Begin * P_Coarsen_Factor



            !
            !  Unpack Data into Temporary Block Storage
            !

            !$OMP DO SCHEDULE(dynamic), COLLAPSE(3)
            DO k = 0,Local_PE_Dim - 1
                DO j = 0,Local_TE_Dim - 1
                    DO i = 0,INPUT_RAD_ELEMS_PER_SHELL - 1

                        Global_Input_R_Loc = Input_Shell_R_Begin + i
                        Block_RE_Loc = i/ R_COARSEN_FACTOR
                        R_Shift_Loc = MOD(i,R_Coarsen_Factor)

                        Global_Input_T_Loc = Input_TE_Begin + j
                        Block_TE_Loc = (Global_Input_T_Loc - Block_Input_TE_Begin)/T_Coarsen_Factor
                        T_Shift_Loc = MOD(j,T_Coarsen_Factor)


                        Global_Input_P_LOC = Input_PE_Begin + k
                        Block_PE_Loc = (Global_Input_P_Loc - Block_Input_PE_Begin)/ P_Coarsen_Factor
                        P_Shift_Loc = MOD(k,P_Coarsen_Factor)



                        Cur_Elem_Num = (k * Local_TE_Dim + j )          &
                                        * INPUT_RAD_ELEMS_PER_SHELL     &
                                     + i



                        DO kd = 0,Local_PD_DIM-1

                            kq = kd*SHIFT_FACTOR_D

                            DO jd = 0,Local_TD_Dim-1

                                jq = jd*SHIFT_FACTOR_E


                                ! Location in Block Quadrature arrangement
                                Start_Here = P_Shift_Loc*Shift_Factor_A      &
                                           + T_Shift_Loc*Shift_Factor_B      &
                                           + R_Shift_Loc*Shift_Factor_C      &
                                           + kq + jq

                                End_Here = Start_Here + Local_RD_DIM - 1


                                ! Location in Recieved Message
                                Start_Here_B = Cur_Elem_Num*Num_Input_DOF       &
                                             + kd*Local_RD_Dim*Local_TD_Dim     &
                                             + jd*Local_RD_Dim                  &
                                             + 3

                                End_Here_B = Start_Here_B + Local_RD_Dim - 1





                                TMP_E_STORAGE( Start_Here:End_Here, Block_RE_Loc, Block_TE_Loc, Block_PE_Loc )     &
                                        = Recv_Message( Start_Here_B:End_Here_B )



                                ! Location in Recieved Message
                                Start_Here_C = Start_Here_B + Size_of_Send
                                End_Here_C = Start_Here_C + Local_RD_Dim - 1

                                TMP_S_STORAGE( Start_Here:End_Here, Block_RE_Loc, Block_TE_Loc, Block_PE_Loc )     &
                                        = Recv_Message( Start_Here_C:End_Here_C )



                                DO d = 1, DOMAIN_DIM

                                    ! Location in Recieved Message
                                    Start_Here_C = Start_Here_B + (d+1)*Size_of_Send
                                    End_Here_C = Start_Here_C + Local_RD_Dim - 1

                                    TMP_Si_STORAGE( Start_Here:End_Here, Block_RE_Loc, Block_TE_Loc, Block_PE_Loc, d )  &
                                            = Recv_Message( Start_Here_C:End_Here_C )

                                END DO ! d Loop


                            END DO ! td Loop
                        END DO ! kd Loop





                    END DO ! i Loop
                END DO ! j Loop
            END DO ! k Loop
            !$OMP END DO

        END DO ! Block_Index Loop

        !
        !   Translate from Input Quadrature to Poseidon Quadrate
        !

        !$OMP DO SCHEDULE(dynamic), COLLAPSE(3)
        DO k = 0,NUM_P_ELEMS_PER_BLOCK - 1
            DO j = 0,NUM_T_ELEMS_PER_BLOCK - 1
                DO i = 0,NUM_R_ELEMS_PER_BLOCK - 1


                    Input_RE_Start = i * R_Coarsen_Factor
                    Delta_FREoCRE(1:R_Coarsen_Factor) = Input_Delta_R(Input_RE_Start+1:Input_RE_Start+R_Coarsen_Factor)     &
                                                      / (rlocs(Block_RE_Begin + i + 1) - rlocs(Block_RE_Begin + i))

                    Input_TE_Start = Block_Input_TE_Begin + j * T_Coarsen_Factor
                    Delta_FTEoCTE(1:T_Coarsen_Factor) = Input_Delta_T(Input_TE_Start+1:Input_TE_Start+T_Coarsen_Factor)     &
                                                      / (tlocs(Block_TE_Begin + j + 1) - tlocs(Block_TE_Begin + j))

                    Input_PE_Start = Block_Input_PE_Begin + k * P_Coarsen_Factor
                    Delta_FPEoCPE(1:P_Coarsen_Factor) = Input_Delta_P(Input_PE_Start+1:Input_PE_Start+P_Coarsen_Factor)     &
                                                      / (plocs(Block_PE_Begin + k + 1) - plocs(Block_PE_Begin + k))

                    DO kk = 0,P_Coarsen_Factor - 1

                        Start_Here = Local_PD_Dim*kk
                        End_Here = Start_Here + Local_PD_Dim - 1

                        Input_CP_Locations(Start_Here:End_Here) = Delta_x_Ratio                     &
                                                                    * Delta_FPEoCPE(kk+1)           &
                                                                    * Input_P_Locations(:)          &
                                                                + 2 * SUM(Delta_FPEoCPE(0:kk))     &
                                                                - 1

                    END DO




                    DO jj = 0,T_Coarsen_Factor - 1

                        Start_Here = Local_TD_Dim*jj
                        End_Here = Start_Here + Local_TD_Dim - 1

                        Input_CT_Locations(Start_Here:End_Here) = Delta_x_Ratio                     &
                                                                    * Delta_FTEoCTE(jj+1)           &
                                                                    * Input_T_Locations(:)          &
                                                                + 2 * SUM(Delta_FTEoCTE(0:jj))     &
                                                                - 1
                    END DO

                    DO ii = 0,R_Coarsen_Factor - 1

                        Start_Here = Local_RD_Dim*ii
                        End_Here = Start_Here + Local_RD_Dim - 1

                        Input_CR_Locations(Start_Here:End_Here) = Delta_x_Ratio                     &
                                                                    * Delta_FREoCRE(ii+1)           &
                                                                    * Input_R_Locations(:)          &
                                                                + 2 * SUM(Delta_FREoCRE(0:ii))     &
                                                                - 1
                    END DO







                    !
                    !   Create Translation Matrix
                    !

                    DO Local_R = 1,NUM_R_QUAD_POINTS
                        R_Lag_Poly_Values(:,Local_R) = Lagrange_Poly(Local_R_Locations(Local_R),  &
                                                                       Coarse_RD_Dim-1,             &
                                                                       Input_CR_Locations           )

                    END DO

                    DO Local_T = 1,NUM_T_QUAD_POINTS
                        T_Lag_Poly_Values(:,Local_T) = Lagrange_Poly(Local_T_Locations(Local_T),  &
                                                                       Coarse_TD_Dim-1,             &
                                                                       Input_CT_Locations           )

                    END DO

                    DO Local_P = 1,NUM_P_QUAD_POINTS
                        P_Lag_Poly_Values(:,Local_P) = Lagrange_Poly(Local_P_Locations(Local_P),  &
                                                                       Coarse_PD_Dim-1,             &
                                                                       Input_CP_Locations           )

                    END DO


                    DO Local_P = 1,NUM_P_QUAD_POINTS
                        DO Local_T = 1,NUM_T_QUAD_POINTS
                            DO Local_R = 1,NUM_R_QUAD_POINTS

                                Local_Here = ((Local_P-1) * NUM_T_QUAD_POINTS + Local_T-1 )   &
                                             * NUM_R_QUAD_POINTS                              &
                                           + Local_R

                                DO Input_P = 1,Coarse_PD_Dim
                                    DO Input_T = 1,Coarse_TD_Dim

                                            Start_Here = ((Input_P-1)*Coarse_TD_DIM + Input_T - 1)  &
                                                            * Coarse_RD_Dim

                                            End_Here = Start_Here + Coarse_RD_Dim

                                            Translation_Matrix(Start_Here+1:End_Here, Local_Here)  =                &
                                                                      R_Lag_Poly_Values(1:Coarse_RD_Dim,Local_R)    &
                                                                    * T_Lag_Poly_Values(Input_T,Local_T)            &
                                                                    * P_Lag_Poly_Values(Input_P,Local_P)

                                    END DO  !   Input_T Loop
                                END DO  !   Input_P Loop

                            END DO  !   Local_R Loop
                        END DO  !   Local_T Loop
                    END DO  !   Local_P looop

                    !
                    !   MV_Mult to Convert
                    !
                    DO Local_Here = 1,Num_Local_DOF

                        Local_E_Coeffs(Local_Here) = DOT_PRODUCT( Translation_Matrix(:,Local_Here), &
                                                                  TMP_E_Storage(:,i,j,k)            )

                        Local_S_Coeffs(Local_Here) = DOT_PRODUCT( Translation_Matrix(:,Local_Here), &
                                                                  TMP_S_Storage(:,i,j,k)            )

                        DO d = 1,DOMAIN_DIM
                            Local_Si_Coeffs(Local_Here,d) = DOT_PRODUCT( Translation_Matrix(:,Local_Here), &
                                                                         TMP_Si_Storage(:,i,j,k,d)         )
                        END DO

                    END DO






                    DO Local_P = 1,NUM_P_QUAD_POINTS
                        DO Local_T = 1,NUM_T_QUAD_POINTS

                             !                                                           !
                            !!   Input current elements source values into the global    !!
                            !!   source vector.                                          !!
                             !                                                           !

                            Local_Here = ((Local_P-1) * NUM_T_QUAD_POINTS + Local_T-1 )   &
                                             * NUM_R_QUAD_POINTS


                            Block_E(1:NUM_R_QUAD_POINTS, Local_T, Local_P, i, j, k )      &
                                = Local_E_Coeffs(Local_Here+1:Local_Here+NUM_R_QUAD_POINTS)

                            Block_S(1:NUM_R_QUAD_POINTS, Local_T, Local_P, i, j, k )      &
                                = Local_S_Coeffs(Local_Here+1:Local_Here+NUM_R_QUAD_POINTS)

                            Do d = 1,DOMAIN_DIM

                                Block_Si(1:NUM_R_QUAD_POINTS, Local_T, Local_P, i, j, k, d )   &
                                    = Local_Si_Coeffs(Local_Here+1:Local_Here+NUM_R_QUAD_POINTS,d)

                            END DO ! d Loop

                        END DO ! Local_T Loop
                    END DO ! Local_P Loop


                END DO ! i Loop
            END DO ! j Loop
        END DO ! k Loop
        !$OMP END DO

    END IF
    !$OMP MASTER
    CALL MPI_WAIT(request, status, ierr)
    !$OMP END MASTER

    !$OMP BARRIER


END DO ! Shell Loop
!$OMP END PARALLEL





END SUBROUTINE Poseidon_CFA_Block_Share











!+201+##########################################################################!
!                                                                               !
!                           Poseidon_Distribute_Solution                        !
!                                                                               !
!###############################################################################!
SUBROUTINE Poseidon_Distribute_Solution()




INTEGER                                                 ::  i, j
INTEGER                                                 ::  ierr

INTEGER                                                 ::  Length_A,       &
                                                            Length_B

INTEGER                                                 ::  Start_Here,     &
                                                            End_Here

COMPLEX(kind = idp), DIMENSION(:), ALLOCATABLE          ::  Buffer_A,       &
                                                            Buffer_B,       &
                                                            Buffer_C

COMPLEX(kind = idp),DIMENSION(:), ALLOCATABLE           ::  Buffer_2

INTEGER, DIMENSION(0:NUM_SHELLS-1)                      ::  ID_Array







Length_A = BLOCK_PROB_DIM - ULM_LENGTH
Length_B = BLOCK_PROB_DIM


ALLOCATE( Buffer_A(0:Length_A-1 ) )
ALLOCATE( Buffer_B(0:Length_B-1 ) )
ALLOCATE( Buffer_C(0:Local_Length-1) )








IF ( SOL_DIST_SCHEME == 1 ) THEN
!
!   POSEIDON_COMM_PETSC BCASTS TO ALL
!

    ID_Array = (/ (i*NUM_BLOCKS_PER_SHELL, i=0,NUM_SHELLS-1, 1) /)
    DO i = 0,NUM_SHELLS-2

        Start_Here = i*NUM_R_ELEMS_PER_BLOCK*DEGREE*ULM_LENGTH
        End_Here = Start_Here + Length_A


!        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

        IF ( myID_PETSC == i ) THEN
            Buffer_A = Coefficient_Vector(Start_Here:End_Here-1)
            CALL MPI_BCAST( Buffer_A, Length_A, MPI_DOUBLE_COMPLEX, ID_Array(i), MPI_COMM_WORLD, ierr )

        ELSE
            CALL MPI_BCAST( Buffer_A, Length_A, MPI_DOUBLE_COMPLEX, ID_Array(i), MPI_COMM_WORLD, ierr )
            Coefficient_Vector(Start_Here:End_Here-1) = Buffer_A

        END IF


    END DO



    i = NUM_SHELLS-1
    Start_Here = i*NUM_R_ELEMS_PER_BLOCK*DEGREE*ULM_LENGTH
    End_Here = Start_Here + Length_B


    IF ( myID_PETSC == i ) THEN

        Buffer_B = Coefficient_Vector(Start_Here:End_Here-1)

        CALL MPI_BCAST( Buffer_B, Length_B, MPI_DOUBLE_COMPLEX, ID_Array(i), MPI_COMM_WORLD, ierr )

    ELSE

        CALL MPI_BCAST( Buffer_B, Length_B, MPI_DOUBLE_COMPLEX, ID_Array(i), MPI_COMM_WORLD, ierr )
        Coefficient_Vector(Start_Here:End_Here-1) = Buffer_B

    END IF





ELSE IF ( SOL_DIST_SCHEME == 2 ) THEN
!
!   Use POSEIDON_COMM_DIST and call simultanious BCASTS
!

    ALLOCATE( Buffer_2(0:PROB_DIM-1) )



    IF ( myID_DIST == 0 ) THEN

        Buffer_2 = Coefficient_Vector
        CALL MPI_BCAST( Buffer_2, PROB_DIM, MPI_DOUBLE_COMPLEX, 0, POSEIDON_COMM_DIST, ierr )

    ELSE

        CALL MPI_BCAST( Buffer_2, PROB_DIM, MPI_DOUBLE_COMPLEX, 0, POSEIDON_COMM_DIST, ierr )
        Coefficient_Vector = Buffer_2

    END IF






ELSE IF ( SOL_DIST_SCHEME == 3 ) THEN
!
!   TWO-STEP BACK UP POSEIDON_COMM TREE
!


    !
    !   BCAST Local Coeffs to Shell
    !
    IF ( POSEIDON_COMM_SHELL .NE. MPI_COMM_NULL ) THEN

        Start_Here = myShell*NUM_R_ELEMS_PER_BLOCK*DEGREE*ULM_LENGTH
        End_Here = Start_Here + Local_Length

        IF ( myID_SHELL == 0 ) THEN

            Buffer_C = Coefficient_Vector(Start_Here:End_Here-1)
            CALL MPI_BCAST( Buffer_C, LOCAL_LENGTH, MPI_DOUBLE_COMPLEX, 0, POSEIDON_COMM_SHELL, ierr )

        ELSE

            CALL MPI_BCAST( Buffer_C, LOCAL_LENGTH, MPI_DOUBLE_COMPLEX, 0, POSEIDON_COMM_SHELL, ierr )
            Coefficient_Vector(Start_Here:End_Here-1) = Buffer_C

        END IF


    END IF



    !
    !   Send to all rays intersecting block
    !




END IF








END SUBROUTINE Poseidon_Distribute_Solution
























!+301+##########################################################################!
!                                                                               !
!                           Poseidon_Newton_Block_Share                         !
!                                                                               !
!###############################################################################!
SUBROUTINE Poseidon_Newton_Block_Share( Proc_id, Proc_Theta_id, Proc_Phi_id,                    &
                                        My_Source_E,                                            &
                                        Local_RE_Dim, Local_TE_Dim, Local_PE_Dim,               &
                                        Local_RD_Dim, Local_TD_Dim, Local_PD_Dim,               &
                                        Input_R_Quad, Input_T_Quad, Input_P_Quad,               &
                                        Left_Limit, Right_Limit,                                &
                                        R_Elements_Input, T_Elements_Input, P_Elements_Input,   &
                                        Input_Delta_R, Input_Delta_T, Input_Delta_P,            &
                                        Block_E                                                 )




INTEGER, INTENT(IN)                                 ::  Proc_id,        &
                                                        Proc_Theta_id,  &
                                                        Proc_Phi_id


REAL(KIND = idp), INTENT(IN), DIMENSION(    1:Local_RD_Dim*Local_TD_Dim*Local_PD_Dim,             &
                                            0:Local_RE_Dim-1,           &
                                            0:Local_TE_Dim-1,           &
                                            0:Local_PE_Dim-1  )         ::  My_Source_E




INTEGER, INTENT(IN)                                                     ::  Local_RE_Dim,   &
                                                                            Local_TE_Dim,   &
                                                                            Local_PE_Dim,   &
                                                                            Local_RD_Dim,   &
                                                                            Local_TD_Dim,   &
                                                                            Local_PD_Dim


REAL(KIND = idp), DIMENSION(1:Local_RD_Dim), INTENT(IN)                 ::  Input_R_Quad
REAL(KIND = idp), DIMENSION(1:Local_TD_Dim), INTENT(IN)                 ::  Input_T_Quad
REAL(KIND = idp), DIMENSION(1:Local_PD_Dim), INTENT(IN)                 ::  Input_P_Quad

REAL(KIND = idp), INTENT(IN)                                            ::  Left_Limit,     &
                                                                            Right_Limit

INTEGER, INTENT(IN)                                                     ::  R_Elements_Input,       &
                                                                            T_Elements_Input,       &
                                                                            P_Elements_Input

REAL(KIND = idp), DIMENSION(1:R_Elements_Input), INTENT(IN)             ::  Input_Delta_R
REAL(KIND = idp), DIMENSION(1:T_Elements_Input), INTENT(IN)             ::  Input_Delta_T
REAL(KIND = idp), DIMENSION(1:P_Elements_Input), INTENT(IN)             ::  Input_Delta_P


REAL(KIND = idp), INTENT(INOUT), DIMENSION( 1:NUM_R_QUAD_POINTS,        &
                                            1:NUM_T_QUAD_POINTS,        &
                                            1:NUM_P_QUAD_POINTS,        &
                                            0:NUM_R_ELEMS_PER_BLOCK-1,  &
                                            0:NUM_T_ELEMS_PER_BLOCK-1,  &
                                            0:NUM_P_ELEMS_PER_BLOCK-1 ) ::  Block_E



INTEGER                                         ::  Block_Num
Integer                                         ::  Shell, Block_indx
INTEGER                                         ::  Input_Shell_R_Begin,      &
                                                    Input_Shell_R_End,        &
                                                    Input_TE_Begin,      &
                                                    Input_PE_Begin


INTEGER                                         ::  Input_Shell_R_Loc,        &
                                                    Input_Shell_T_Loc,        &
                                                    Input_Shell_P_Loc

INTEGER                                         ::  i, j, k, d, n
INTEGER                                         ::  ii, jj, kk
INTEGER                                         ::  rd, td, pd
INTEGER                                         ::  jd, kd
INTEGER                                         ::  jq, kq


INTEGER                                         ::  Block_R_Loc,        &
                                                    Block_T_Loc,        &
                                                    Block_P_Loc

INTEGER                                         ::  Input_RE_Start,     &
                                                    Input_TE_Start,     &
                                                    Input_PE_Start

INTEGER                                                     ::  re, te, pe,                 &
                                                                Local_R, Local_T, Local_P,  &
                                                                Input_R, Input_T, Input_P,  &
                                                                Local_Here, Input_Here

INTEGER                                                     ::  Num_Local_DOF,              &
                                                                Num_Input_DOF,              &
                                                                Num_Coarse_DOF





REAL(KIND = idp), DIMENSION(1:NUM_R_QUAD_POINTS)            ::  Local_R_Locations
REAL(KIND = idp), DIMENSION(1:NUM_T_QUAD_POINTS)            ::  Local_T_Locations
REAL(KIND = idp), DIMENSION(1:NUM_P_QUAD_POINTS)            ::  Local_P_Locations

REAL(KIND = idp), DIMENSION(1:Local_RD_Dim)                 ::  Input_R_Locations
REAL(KIND = idp), DIMENSION(1:Local_TD_Dim)                 ::  Input_T_Locations
REAL(KIND = idp), DIMENSION(1:Local_PD_Dim)                 ::  Input_P_Locations







REAL(KIND = idp), DIMENSION(:,:), ALLOCATABLE                 ::  R_Lag_Poly_Values
REAL(KIND = idp), DIMENSION(:,:), ALLOCATABLE                 ::  T_Lag_Poly_Values
REAL(KIND = idp), DIMENSION(:,:), ALLOCATABLE                 ::  P_Lag_Poly_Values

REAL(KIND = idp), DIMENSION(:,:), ALLOCATABLE               ::  Translation_Matrix

REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Local_E_Coeffs,         &
                                                                Local_S_Coeffs

REAL(KIND = idp), DIMENSION(:,:), ALLOCATABLE               ::  Local_Si_Coeffs


REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)                 ::  Input_CR_Locations
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)                 ::  Input_CT_Locations
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)                 ::  Input_CP_Locations

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:,:)           ::  TMP_E_STORAGE,         &
                                                                TMP_S_STORAGE

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:,:,:,:)         ::  TMP_Si_STORAGE


REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Send_Message

REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Recv_Message

REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Delta_FREoCRE
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Delta_FTEoCTE
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Delta_FPEoCPE

REAL(KIND = idp)                                            ::  Input_Delta_x
REAL(KIND = idp)                                            ::  Delta_x_Ratio

INTEGER                                                     ::  Size_of_Send
INTEGER                                                     ::  Size_of_Send_b



INTEGER                                                     ::  Elems_Per_Block

INTEGER                                                     ::  Cur_Elem_Num,           &
                                                                Next_Elem_Num

INTEGER                                                     ::  Start_Here, End_Here
INTEGER                                                     ::  Start_Here_B, End_Here_B
INTEGER                                                     ::  Start_Here_C, End_Here_C

INTEGER                                                     ::  block_j, block_k, Recv_Addr

INTEGER                                                     ::  TheirID, TheirID_theta, TheirID_phi

INTEGER                                                     ::  Block_RE_Begin,Block_TE_Begin, Block_PE_Begin

INTEGER                                                     ::  Block_RE_Loc,                &
                                                                Block_TE_Loc,                &
                                                                Block_PE_Loc

INTEGER                                                     ::  Block_Input_RE_Begin,        &
                                                                Block_Input_TE_Begin,        &
                                                                Block_Input_PE_Begin

INTEGER                                                     ::  Global_Input_R_Loc,          &
                                                                Global_Input_T_Loc,          &
                                                                Global_Input_P_Loc

INTEGER                                                     ::  Input_RAD_Elems_Per_Shell
INTEGER                                                     ::  Block_Index

INTEGER                                                     ::  print_indx
INTEGER                                                     ::  request
INTEGER, DIMENSION(MPI_STATUS_SIZE)                         ::  status

INTEGER                                                     ::  Coarsen_Factor

INTEGER                                                     ::  Coarse_RD_Dim,  &
                                                                Coarse_TD_Dim,  &
                                                                Coarse_PD_Dim

INTEGER                                                     ::  Shift_Factor_A, &
                                                                Shift_Factor_B, &
                                                                Shift_Factor_C, &
                                                                Shift_Factor_D, &
                                                                Shift_Factor_E

INTEGER                                                     ::  P_Shift_Loc,    &
                                                                T_Shift_Loc,    &
                                                                R_Shift_Loc

INTEGER             ::  Cur_Elem_Num_B,  &
                        Local_Here_B


!PRINT*,"In Poseidon_Coarse_Block_Share A"



Num_Input_DOF = Local_RD_Dim * Local_TD_Dim * Local_PD_Dim

Num_Local_DOF = NUM_R_QUAD_POINTS * NUM_T_QUAD_POINTS * NUM_P_QUAD_POINTS

INPUT_RAD_ELEMS_PER_SHELL = NUM_R_ELEMS_PER_SHELL*R_Coarsen_Factor

Elems_Per_Block = Local_TE_Dim*Local_PE_Dim*INPUT_RAD_ELEMS_PER_SHELL
Size_of_Send = Elems_Per_Block * Num_Input_DOF
Size_of_Send_b = Elems_Per_Block * Num_Input_DOF

Coarsen_Factor = R_Coarsen_Factor*T_Coarsen_Factor*P_Coarsen_Factor

Coarse_RD_Dim = R_Coarsen_Factor * Local_RD_Dim
Coarse_TD_Dim = T_Coarsen_Factor * Local_TD_Dim
Coarse_PD_Dim = P_Coarsen_Factor * Local_PD_Dim

Num_Coarse_DOF = Coarse_RD_Dim * Coarse_TD_Dim * Coarse_PD_Dim

Shift_Factor_A = R_Coarsen_Factor*T_Coarsen_Factor*Local_RD_Dim*Local_TD_Dim*Local_PD_Dim
Shift_Factor_B = R_Coarsen_Factor*Local_RD_Dim*Local_TD_Dim
Shift_Factor_C = Local_RD_Dim
Shift_Factor_D = R_Coarsen_Factor*T_Coarsen_Factor*Local_RD_Dim*Local_TD_Dim
Shift_Factor_E = R_Coarsen_Factor*Local_RD_Dim

Delta_x_Ratio = 2.0_idp/(Right_Limit - Left_Limit)


ALLOCATE(Translation_Matrix(1:Num_Coarse_DOF, 1:Num_Local_DOF))

ALLOCATE(Local_E_Coeffs(1:Num_Local_DOF))
ALLOCATE(Local_S_Coeffs(1:Num_Local_DOF))
ALLOCATE(Local_Si_Coeffs(1:Num_Local_DOF,1:DOMAIN_DIM))


ALLOCATE( TMP_E_STORAGE(0:Num_Coarse_DOF-1,             &
                        0:NUM_R_ELEMS_PER_BLOCK-1,      &
                        0:NUM_T_ELEMS_PER_BLOCK-1,      &
                        0:NUM_P_ELEMS_PER_BLOCK-1)      )

ALLOCATE( TMP_S_STORAGE(0:Num_Coarse_DOF-1,             &
                        0:NUM_R_ELEMS_PER_BLOCK-1,      &
                        0:NUM_T_ELEMS_PER_BLOCK-1,      &
                        0:NUM_P_ELEMS_PER_BLOCK-1)      )

ALLOCATE( TMP_Si_STORAGE(0:Num_Coarse_DOF-1,            &
                         0:NUM_R_ELEMS_PER_BLOCK-1,     &
                         0:NUM_T_ELEMS_PER_BLOCK-1,     &
                         0:NUM_P_ELEMS_PER_BLOCK-1,     &
                         1:DOMAIN_DIM)                  )


ALLOCATE( Input_CR_Locations(0:Coarse_RD_Dim - 1) )
ALLOCATE( Input_CT_Locations(0:Coarse_TD_Dim - 1) )
ALLOCATE( Input_CP_Locations(0:Coarse_PD_Dim - 1) )


ALLOCATE( Send_Message(0:Size_of_Send_b+2) )
ALLOCATE( Recv_Message(0:Size_of_Send_b+2) )


ALLOCATE( Delta_FREoCRE(0:R_Coarsen_Factor) )
ALLOCATE( Delta_FTEoCTE(0:T_Coarsen_Factor) )
ALLOCATE( Delta_FPEoCPE(0:P_Coarsen_Factor) )



ALLOCATE( R_Lag_Poly_Values(1:Coarse_RD_Dim,1:NUM_R_QUAD_POINTS) )
ALLOCATE( T_Lag_Poly_Values(1:Coarse_TD_Dim,1:NUM_T_QUAD_POINTS) )
ALLOCATE( P_Lag_Poly_Values(1:Coarse_PD_Dim,1:NUM_P_QUAD_POINTS) )


Local_R_Locations = Initialize_LG_Quadrature_Locations(NUM_R_QUAD_POINTS)
Local_T_Locations = Initialize_LG_Quadrature_Locations(NUM_T_QUAD_POINTS)
Local_P_Locations = Initialize_LG_Quadrature_Locations(NUM_P_QUAD_POINTS)



Input_R_Locations = Input_R_Quad - Left_Limit
Input_T_Locations = Input_T_Quad - Left_Limit
Input_P_Locations = Input_P_Quad - Left_Limit


!$OMP PARALLEL DEFAULT(none)                                                    &
!$OMP PRIVATE(   Shell, k, j, i, kd, jd, kq, jq, d,                             &
!$OMP            kk, jj, ii,                                                    &
!$OMP            Local_R, Local_T, Local_P, Input_R, Input_T, Input_P,          &
!$OMP            Block_Index,                                                   &
!$OMP            Input_Shell_R_Begin, Input_Shell_R_End,                        &
!$OMP            Input_Shell_R_Loc, Cur_Elem_Num, Start_Here, End_Here,         &
!$OMP            Start_Here_B, End_Here_B, Start_Here_C, End_Here_C,            &
!$OMP            Local_Here, Input_Here,                                        &
!$OMP            block_j, block_k, Recv_Addr,                                   &
!$OMP            TheirID, TheirID_theta, TheirID_phi,                           &
!$OMP            Input_TE_Begin, Input_PE_Begin,                                &
!$OMP            Block_RE_Begin, Block_TE_Begin, Block_PE_Begin,                &
!$OMP            Block_Input_TE_Begin, Block_Input_PE_Begin,                    &
!$OMP            Global_Input_R_Loc, Global_Input_T_Loc, Global_Input_P_Loc,    &
!$OMP            Block_RE_Loc, Block_TE_Loc, Block_PE_Loc,                      &
!$OMP            R_Shift_Loc, T_Shift_Loc, P_Shift_Loc,                         &
!$OMP            Delta_FREoCRE, Delta_FTEoCTE, Delta_FPEoCPE,                   &
!$OMP            Input_PE_Start, Input_TE_Start, Input_RE_Start,                &
!$OMP            Input_CP_Locations, Input_CT_Locations, Input_CR_Locations,    &
!$OMP            R_Lag_Poly_Values, T_Lag_Poly_Values, P_Lag_Poly_Values,       &
!$OMP            Local_E_Coeffs, Local_S_Coeffs, Local_Si_Coeffs,               &
!$OMP            Translation_Matrix                                        )    &
!$OMP SHARED(    Local_PE_Dim, Local_TE_Dim, INPUT_RAD_ELEMS_PER_SHELL,         &
!$OMP            Local_PD_Dim, Local_TD_Dim, Local_RD_Dim,                      &
!$OMP            NUM_R_ELEMS_PER_BLOCK, NUM_T_ELEMS_PER_BLOCK, NUM_P_ELEMS_PER_BLOCK, &
!$OMP            R_COARSEN_FACTOR, T_COARSEN_FACTOR, P_COARSEN_FACTOR,          &
!$OMP            DOMAIN_DIM, NUM_SHELLS, NUM_BLOCK_THETA_ROWS, NUM_BLOCKS_PER_SHELL,  &
!$OMP            NUM_R_QUAD_POINTS, NUM_T_QUAD_POINTS, NUM_P_QUAD_POINTS,       &
!$OMP            rlocs, tlocs, plocs,                                           &
!$OMP            Input_Delta_P, Input_Delta_T, Input_Delta_R,                   &
!$OMP            Ratio_T_BNDLperBLCK, Ratio_P_BNDLperBLCK,                      &
!$OMP            Ratio_BNDLperBLCK,                                             &
!$OMP            Send_Message, Recv_Message,                                    &
!$OMP            Size_of_Send, Size_of_Send_b,                                  &
!$OMP            My_Source_E,                                                   &
!$OMP            Num_Input_DOF, Num_Local_DOF, Num_Coarse_DOF,                  &
!$OMP            Coarse_RD_Dim, Coarse_TD_Dim, Coarse_PD_Dim,                   &
!$OMP            Shift_Factor_A, Shift_Factor_B, Shift_Factor_C,                &
!$OMP            Shift_Factor_D, Shift_Factor_E,                                &
!$OMP            Delta_x_Ratio,                                                 &
!$OMP            Local_R_Locations, Local_T_Locations, Local_P_Locations,       &
!$OMP            Input_R_Locations, Input_T_Locations, Input_P_Locations,       &
!$OMP            myID, myID_Theta, myID_Phi, myShell, myID_Shell, ierr,         &
!$OMP            request, sTatus, MPI_STATUS_IGNORE,                            &
!$OMP            TMP_E_STORAGE, TMP_S_STORAGE, TMP_Si_STORAGE,                  &
!$OMP            Block_E                                                        )


Delta_FREoCRE(0) = 0.0_idp
Delta_FTEoCTE(0) = 0.0_idp
Delta_FPEoCPE(0) = 0.0_idp

!Input_TE_Begin = Proc_Theta_id*NUM_LOC_T_ELEMENTS
!Input_PE_Begin = Proc_Phi_id*NUM_LOC_P_ELEMENTS


!PRINT*,"IN Poseidon_Coarse_Block_Share B"


DO Shell = 0, NUM_SHELLS - 1


    !
    !   Package and Send Data
    !

    Input_Shell_R_Begin = Shell * INPUT_RAD_ELEMS_PER_SHELL
    Input_Shell_R_End =  (Shell + 1) * INPUT_RAD_ELEMS_PER_SHELL - 1


    !$OMP DO SCHEDULE(dynamic), COLLAPSE(3)
    DO k = 0,Local_PE_Dim - 1
        DO j = 0,Local_TE_Dim - 1
            DO i = 0,INPUT_RAD_ELEMS_PER_SHELL - 1

                Input_Shell_R_Loc = Input_Shell_R_Begin + i

                Cur_Elem_Num = (k * Local_TE_Dim + j )          &
                                * INPUT_RAD_ELEMS_PER_SHELL     &
                             + i




                Start_Here = Cur_Elem_Num*Num_Input_DOF + 3
                End_Here = (Cur_Elem_Num+1)*Num_Input_DOF+2


                Send_Message(Start_Here:End_Here) = My_Source_E(1:NUM_Input_DOF, Input_Shell_R_Loc, j, k )








            END DO ! k Loop

        END DO ! j Loop

    END DO ! i Loop
    !$OMP END DO



    !$OMP MASTER
    Send_Message(0) = REAL( myID, idp )
    Send_Message(1) = REAL( myID_Theta,idp)
    Send_Message(2) = REAL( myID_Phi, idp)

    block_j = myID_Theta/Ratio_T_BNDLperBLCK
    block_k = myID_Phi/Ratio_P_BNDLperBLCK

    Recv_Addr = block_k * NUM_BLOCK_THETA_ROWS + block_j + Shell*NUM_BLOCKS_PER_SHELL

    CALL MPI_ISend( Send_Message, Size_of_Send_b+3, MPI_DOUBLE, Recv_Addr, 6, MPI_COMM_WORLD, request, ierr)





    !
    !   If Process is responsible for a block in the current shell
    !   then receive data.
    !
    IF ( myShell == Shell ) THEN


        DO Block_index = 1,Ratio_BNDLperBLCK


            !$OMP MASTER
            CALL MPI_Recv( Recv_Message, Size_of_Send_b+3, MPI_DOUBLE, MPI_ANY_SOURCE,      &
                                6, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr                  )

            !$OMP END MASTER
            !$OMP BARRIER



            TheirID = INT(Recv_Message(0))
            TheirID_theta = INT(Recv_Message(1))
            TheirID_phi = INT(Recv_Message(2))


!            Bundle_List(Block_index) = TheirID

            Input_TE_Begin = TheirID_theta*LOCAL_TE_DIM
            Input_PE_Begin = TheirID_phi*LOCAL_PE_DIM

            Block_RE_Begin = Shell*NUM_R_ELEMS_PER_BLOCK
            Block_TE_Begin = MOD(myID_Shell,NUM_BLOCK_THETA_ROWS)*NUM_T_ELEMS_PER_BLOCK
            Block_PE_Begin = (myID_Shell/NUM_BLOCK_THETA_ROWS)*NUM_P_ELEMS_PER_BLOCK

            Block_Input_TE_Begin = Block_TE_Begin * T_Coarsen_Factor
            Block_Input_PE_Begin = Block_PE_Begin * P_Coarsen_Factor



            !
            !  Unpack Data into Temporary Block Storage
            !

            !$OMP DO SCHEDULE(dynamic), COLLAPSE(3)
            DO k = 0,Local_PE_Dim - 1
                DO j = 0,Local_TE_Dim - 1
                    DO i = 0,INPUT_RAD_ELEMS_PER_SHELL - 1

                        Global_Input_R_Loc = Input_Shell_R_Begin + i
                        Block_RE_Loc = i/ R_COARSEN_FACTOR
                        R_Shift_Loc = MOD(i,R_Coarsen_Factor)

                        Global_Input_T_Loc = Input_TE_Begin + j
                        Block_TE_Loc = (Global_Input_T_Loc - Block_Input_TE_Begin)/T_Coarsen_Factor
                        T_Shift_Loc = MOD(j,T_Coarsen_Factor)


                        Global_Input_P_LOC = Input_PE_Begin + k
                        Block_PE_Loc = (Global_Input_P_Loc - Block_Input_PE_Begin)/ P_Coarsen_Factor
                        P_Shift_Loc = MOD(k,P_Coarsen_Factor)



                        Cur_Elem_Num = (k * Local_TE_Dim + j )          &
                                        * INPUT_RAD_ELEMS_PER_SHELL     &
                                     + i



                        DO kd = 0,Local_PD_DIM-1

                            kq = kd*SHIFT_FACTOR_D

                            DO jd = 0,Local_TD_Dim-1

                                jq = jd*SHIFT_FACTOR_E


                                ! Location in Block Quadrature arrangement
                                Start_Here = P_Shift_Loc*Shift_Factor_A      &
                                           + T_Shift_Loc*Shift_Factor_B      &
                                           + R_Shift_Loc*Shift_Factor_C      &
                                           + kq + jq

                                End_Here = Start_Here + Local_RD_DIM - 1


                                ! Location in Recieved Message
                                Start_Here_B = Cur_Elem_Num*Num_Input_DOF       &
                                             + kd*Local_RD_Dim*Local_TD_Dim     &
                                             + jd*Local_RD_Dim                  &
                                             + 3

                                End_Here_B = Start_Here_B + Local_RD_Dim - 1





                                TMP_E_STORAGE( Start_Here:End_Here, Block_RE_Loc, Block_TE_Loc, Block_PE_Loc )     &
                                        = Recv_Message( Start_Here_B:End_Here_B )



                                ! Location in Recieved Message
                                Start_Here_C = Start_Here_B + Size_of_Send
                                End_Here_C = Start_Here_C + Local_RD_Dim - 1

                                TMP_S_STORAGE( Start_Here:End_Here, Block_RE_Loc, Block_TE_Loc, Block_PE_Loc )     &
                                        = Recv_Message( Start_Here_C:End_Here_C )



                                DO d = 1, DOMAIN_DIM

                                    ! Location in Recieved Message
                                    Start_Here_C = Start_Here_B + (d+1)*Size_of_Send
                                    End_Here_C = Start_Here_C + Local_RD_Dim - 1

                                    TMP_Si_STORAGE( Start_Here:End_Here, Block_RE_Loc, Block_TE_Loc, Block_PE_Loc, d )  &
                                            = Recv_Message( Start_Here_C:End_Here_C )

                                END DO ! d Loop


                            END DO ! td Loop
                        END DO ! kd Loop





                    END DO ! i Loop
                END DO ! j Loop
            END DO ! k Loop
            !$OMP END DO


        END DO ! Block_Index Loop




        !
        !   Translate from Input Quadrature to Poseidon Quadrate
        !


        !$OMP DO SCHEDULE(dynamic), COLLAPSE(3)
        DO k = 0,NUM_P_ELEMS_PER_BLOCK - 1
            DO j = 0,NUM_T_ELEMS_PER_BLOCK - 1
                DO i = 0,NUM_R_ELEMS_PER_BLOCK - 1

                    Input_RE_Start = i * R_Coarsen_Factor
                    Delta_FREoCRE(1:R_Coarsen_Factor) = Input_Delta_R(Input_RE_Start+1:Input_RE_Start+R_Coarsen_Factor)     &
                                                      / (rlocs(Block_RE_Begin + i + 1) - rlocs(Block_RE_Begin + i))

                    Input_TE_Start = Block_Input_TE_Begin + j * T_Coarsen_Factor
                    Delta_FTEoCTE(1:T_Coarsen_Factor) = Input_Delta_T(Input_TE_Start+1:Input_TE_Start+T_Coarsen_Factor)     &
                                                      / (tlocs(Block_TE_Begin + j + 1) - tlocs(Block_TE_Begin + j))

                    Input_PE_Start = Block_Input_PE_Begin + k * P_Coarsen_Factor
                    Delta_FPEoCPE(1:P_Coarsen_Factor) = Input_Delta_P(Input_PE_Start+1:Input_PE_Start+P_Coarsen_Factor)     &
                                                      / (plocs(Block_PE_Begin + k + 1) - plocs(Block_PE_Begin + k))


                    DO kk = 0,P_Coarsen_Factor - 1

                        Start_Here = Local_PD_Dim*kk
                        End_Here = Start_Here + Local_PD_Dim - 1

                        Input_CP_Locations(Start_Here:End_Here) = Delta_x_Ratio                     &
                                                                    * Delta_FPEoCPE(kk+1)           &
                                                                    * Input_P_Locations(:)          &
                                                                + 2 * SUM(Delta_FPEoCPE(0:kk))     &
                                                                - 1

                    END DO




                    DO jj = 0,T_Coarsen_Factor - 1

                        Start_Here = Local_TD_Dim*jj
                        End_Here = Start_Here + Local_TD_Dim - 1

                        Input_CT_Locations(Start_Here:End_Here) = Delta_x_Ratio                     &
                                                                    * Delta_FTEoCTE(jj+1)           &
                                                                    * Input_T_Locations(:)          &
                                                                + 2 * SUM(Delta_FTEoCTE(0:jj))     &
                                                                - 1
                    END DO

                    DO ii = 0,R_Coarsen_Factor - 1

                        Start_Here = Local_RD_Dim*ii
                        End_Here = Start_Here + Local_RD_Dim - 1

                        Input_CR_Locations(Start_Here:End_Here) = Delta_x_Ratio                     &
                                                                    * Delta_FREoCRE(ii+1)           &
                                                                    * Input_R_Locations(:)          &
                                                                + 2 * SUM(Delta_FREoCRE(0:ii))     &
                                                                - 1
                    END DO








                    !
                    !   Create Translation Matrix
                    !

                    DO Local_R = 1,NUM_R_QUAD_POINTS
                        R_Lag_Poly_Values(:,Local_R) = Lagrange_Poly(Local_R_Locations(Local_R),  &
                                                                       Coarse_RD_Dim-1,             &
                                                                       Input_CR_Locations           )

                    END DO


                    DO Local_T = 1,NUM_T_QUAD_POINTS
                        T_Lag_Poly_Values(:,Local_T) = Lagrange_Poly(Local_T_Locations(Local_T),  &
                                                                       Coarse_TD_Dim-1,             &
                                                                       Input_CT_Locations           )

                    END DO

                    DO Local_P = 1,NUM_P_QUAD_POINTS
                        P_Lag_Poly_Values(:,Local_P) = Lagrange_Poly(Local_P_Locations(Local_P),  &
                                                                       Coarse_PD_Dim-1,             &
                                                                       Input_CP_Locations           )

                    END DO



                    DO Local_P = 1,NUM_P_QUAD_POINTS
                        DO Local_T = 1,NUM_T_QUAD_POINTS
                            DO Local_R = 1,NUM_R_QUAD_POINTS

                                Local_Here = ((Local_P-1) * NUM_T_QUAD_POINTS + Local_T-1 )   &
                                             * NUM_R_QUAD_POINTS                              &
                                           + Local_R

                                DO Input_P = 1,Coarse_PD_Dim
                                    DO Input_T = 1,Coarse_TD_Dim

                                            Start_Here = ((Input_P-1)*Coarse_TD_DIM + Input_T - 1)  &
                                                            * Coarse_RD_Dim

                                            End_Here = Start_Here + Coarse_RD_Dim

                                            Translation_Matrix(Start_Here+1:End_Here, Local_Here)  =                &
                                                                      R_Lag_Poly_Values(1:Coarse_RD_Dim,Local_R)    &
                                                                    * T_Lag_Poly_Values(Input_T,Local_T)            &
                                                                    * P_Lag_Poly_Values(Input_P,Local_P)

                                    END DO  !   Input_T Loop
                                END DO  !   Input_P Loop

                            END DO  !   Local_R Loop
                        END DO  !   Local_T Loop
                    END DO  !   Local_P looop


                    !
                    !   MV_Mult to Convert
                    !

                    DO Local_Here = 1,Num_Local_DOF

                            Local_E_Coeffs(Local_Here) = DOT_PRODUCT( Translation_Matrix(:,Local_Here), &
                                                                      TMP_E_Storage(:,i,j,k)            )


                    END DO







                    DO Local_P = 1,NUM_P_QUAD_POINTS
                        DO Local_T = 1,NUM_T_QUAD_POINTS

                             !                                                           !
                            !!   Input current elements source values into the global    !!
                            !!   source vector.                                          !!
                             !                                                           !

                            Local_Here = ((Local_P-1) * NUM_T_QUAD_POINTS + Local_T-1 )   &
                                             * NUM_R_QUAD_POINTS


                            Block_E(1:NUM_R_QUAD_POINTS, Local_T, Local_P, i, j, k )      &
                                = Local_E_Coeffs(Local_Here+1:Local_Here+NUM_R_QUAD_POINTS)


                        END DO ! Local_T Loop
                    END DO ! Local_P Loop





                END DO ! i Loop
            END DO ! j Loop
        END DO ! k Loop
        !$OMP END DO

    END IF

    !$OMP MASTER
    CALL MPI_WAIT(request, status, ierr)
    !$OMP END MASTER

    !$OMP BARRIER

!    CALL MPI_BARRIER(POSEIDON_COMM_WORLD,ierr)

END DO ! Shell Loop
!$OMP END PARALLEL


!IF (myID == 0 ) THEN
!    PRINT*,"Block_E"
!    PRINT*, Block_E
!    PRINT*," "
!END IF



END SUBROUTINE Poseidon_Newton_Block_Share



















END MODULE Poseidon_Internal_Communication_Module
