   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Initialization_Poseidon                                                      !##!
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

USE Units_Module, &
                ONLY :  Set_Units

USE Poseidon_Parameters, &
                ONLY :  Domain_Dim,             &
                        Degree,                 &
                        L_Limit,                &
                        Solver_Type,            &
                        Verbose_Flag,           &
                        Convergence_Criteria,   &
                        Num_CFA_Vars,           &
                        Max_Iterations

USE Variables_IO, &
                ONLY :  RUN_REPORT_FILE_ID,     &
                        File_Suffix

USE Variables_Functions, &
                ONLY :  LM_Location

USE Variables_Quadrature, &
                ONLY :  Num_R_Quad_Points,      &
                        Num_T_Quad_Points,      &
                        Num_P_Quad_Points,      &
                        Num_TP_Quad_Points,     &
                        Num_Quad_DOF

USE Variables_Mesh, &
                ONLY :  Num_R_Elements,         &
                        Num_T_Elements,         &
                        Num_P_Elements,         &
                        Num_Loc_R_Elements,     &
                        Num_Loc_T_Elements,     &
                        Num_Loc_P_Elements,     &
                        rlocs,                  &
                        tlocs,                  &
                        plocs,                  &
                        drlocs,                 &
                        dtlocs,                 &
                        dplocs,                 &
                        R_Inner,                &
                        R_Outer,                &
                        R_Coarsen_Factor,       &
                        T_Coarsen_Factor,       &
                        P_Coarsen_Factor,       &
                        locs_set,               &
                        dlocs_set

USE Variables_Derived, &
                ONLY :  Prob_Dim,               &
                        Block_Prob_Dim,         &
                        SubShell_Prob_Dim,      &
                        Elem_Prob_Dim_Sqr,      &
                        Var_Dim,                &
                        Block_Var_Dim,          &
                        Num_Off_Diagonals,      &
                        ULM_Length,             &
                        Num_R_Nodes
                
USE Variables_MPI, &
                ONLY :  myID_Poseidon,          &
                        nProcs_Poseidon,        &
                        Num_Blocks_Per_Shell,   &
                        Num_R_Elems_Per_Block,  &
                        Num_T_Elems_Per_Block,  &
                        Num_P_Elems_Per_Block,  &
                        Num_R_Elems_Per_Shell,  &
                        Num_Shells,             &
                        Num_SubShells,          &
                        Num_SubShells_Per_Shell

USE Allocation_Core, &
                ONLY :  Allocate_Poseidon_CFA_Variables

USE Allocation_Mesh, &
                ONLY :  Allocate_Mesh

USE Initialization_Mesh, &
                ONLY :  Initialize_Mesh

USE Initialization_Quadrature, &
                ONLY :  Initialize_Quadrature

USE Initialization_MPI, &
                ONLY :  Initialize_MPI

USE Initialization_Tables, &
                ONLY :  Initialize_Tables

USE Initialization_Derived, &
                ONLY :  Initialize_Derived

USE Initialization_FP, &
                ONLY :  Initialize_FP

USE Initialization_NR, &
                ONLY :  Initialize_NR

USE Functions_Mapping, &
                ONLY :  CFA_3D_LM_Map

IMPLICIT NONE


PUBLIC :: Initialize_Poseidon


CONTAINS

 !+101+####################################################################################!
!                                                                                           !
!       Initialize_Poseidon                                                                 !
!                                                                                           !
 !#########################################################################################!
SUBROUTINE Initialize_Poseidon( Dimensions_Option,                      &
                                FEM_Degree_Option,                      &
                                L_Limit_Option,                         &
                                Units_Option,                           &
                                Domain_Edge_Option,                     &
                                NE_Option,                              &
                                NQ_Option,                              &
                                Coarsen_Option,                         &
                                r_Option, t_Option, p_Option,           &
                                dr_Option, dt_Option, dp_Option,        &
                                Solver_Type_Option,                     &
                                CFA_Eq_Flags_Option,                    &
                                nProcs_Option,                          &
                                Suffix_Flag_Option,                     &
                                Frame_Option,                           &
                                Verbose_Option                       )


CHARACTER(LEN=1),        INTENT(IN), OPTIONAL               ::  Units_Option

INTEGER,                 INTENT(IN), OPTIONAL               ::  FEM_Degree_Option
INTEGER,                 INTENT(IN), OPTIONAL               ::  L_Limit_Option

INTEGER,   DIMENSION(3), INTENT(IN), OPTIONAL               ::  NQ_Option
INTEGER,   DIMENSION(3), INTENT(IN), OPTIONAL               ::  NE_Option
REAL(idp), DIMENSION(2), INTENT(IN), OPTIONAL               ::  Domain_Edge_Option

INTEGER,                 INTENT(IN), OPTIONAL               ::  Solver_Type_Option
LOGICAL,                 INTENT(IN), OPTIONAL               ::  Verbose_Option

INTEGER,                 INTENT(IN), OPTIONAL               ::  Dimensions_Option

INTEGER,   DIMENSION(3), INTENT(IN), OPTIONAL               ::  Coarsen_Option

REAL(idp), DIMENSION(:), INTENT(IN), OPTIONAL               ::  r_Option
REAL(idp), DIMENSION(:), INTENT(IN), OPTIONAL               ::  t_Option
REAL(idp), DIMENSION(:), INTENT(IN), OPTIONAL               ::  p_Option
REAL(idp), DIMENSION(:), INTENT(IN), OPTIONAL               ::  dr_Option
REAL(idp), DIMENSION(:), INTENT(IN), OPTIONAL               ::  dt_Option
REAL(idp), DIMENSION(:), INTENT(IN), OPTIONAL               ::  dp_Option

INTEGER,   DIMENSION(5), INTENT(IN), OPTIONAL               ::  CFA_EQ_Flags_Option
INTEGER,                 INTENT(IN), OPTIONAL               ::  nProcs_Option

CHARACTER(LEN=10),       INTENT(IN), OPTIONAL               ::  Suffix_Flag_Option
INTEGER,                 INTENT(IN), OPTIONAL               ::  Frame_Option

IF ( PRESENT( Verbose_Option ) ) THEN
    Verbose_Flag = Verbose_Option
ELSE
    Verbose_Flag = .TRUE.
END IF

IF ( Verbose_Flag .EQV. .TRUE. ) THEN
    PRINT*,"Initializing Poseidon..."
END IF






IF ( PRESENT( Units_Option ) ) THEN
    CALL Set_Units(Units_Option)
ELSE
    CALL Set_Units("G")
END IF

IF ( PRESENT( FEM_Degree_Option ) ) THEN
    Degree = FEM_Degree_Option
ELSE
    Degree = 1
END IF


IF ( PRESENT( L_Limit_Option ) ) THEN
    L_Limit = L_Limit_Option
ELSE
    L_Limit = 0
END IF


IF ( PRESENT(nProcs_Option) ) THEN
    nProcs_Poseidon = nProcs_Option
ELSE
    nProcs_Poseidon = 1
END IF

IF ( PRESENT( NQ_Option ) ) THEN
    Num_R_Quad_Points = NQ_Option(1)
    Num_T_Quad_Points = NQ_Option(2)
    Num_P_Quad_Points = NQ_Option(3)
ELSE
    Num_R_Quad_Points = 1
    Num_T_Quad_Points = 1
    Num_P_Quad_Points = 1
END IF
Num_TP_Quad_Points = Num_T_Quad_Points*Num_P_Quad_Points
Num_Quad_DOF       = Num_R_Quad_Points*Num_TP_Quad_Points




!-------------------------------------------!
!                                           !
!                   Mesh                    !
!                                           !
!-------------------------------------------!
IF ( PRESENT( NE_Option ) ) THEN
    Num_R_Elements = NE_Option(1)
    Num_T_Elements = NE_Option(2)
    Num_P_Elements = NE_Option(3)
ELSE
    Num_R_Elements = 1
    Num_T_Elements = 1
    Num_P_Elements = 1
END IF



IF ( PRESENT( Coarsen_Option) ) THEN
    R_Coarsen_Factor = Coarsen_Option(1)
    T_Coarsen_Factor = Coarsen_Option(2)
    P_Coarsen_Factor = Coarsen_Option(3)
ELSE
    R_Coarsen_Factor = 1
    T_Coarsen_Factor = 1
    P_Coarsen_Factor = 1
END IF
Num_Loc_R_Elements = Num_R_Elements
Num_Loc_T_Elements = Num_R_Elements
Num_Loc_P_Elements = Num_R_Elements





CALL Allocate_Mesh()

IF ( PRESENT( Domain_Edge_Option ) ) THEN
    R_Inner = Domain_Edge_Option(1)
    R_Outer = Domain_Edge_Option(2)
ELSE
    R_Inner = 0.0_idp
    R_Outer = 1.0_idp
END IF

IF ( PRESENT( r_Option ) ) THEN
    rlocs = r_Option
    locs_set(1) = .TRUE.
END IF
IF ( PRESENT( dr_Option ) ) THEN
    drlocs = dr_Option
    dlocs_set(1) = .TRUE.
END IF

IF ( PRESENT( t_Option ) ) THEN
    tlocs = t_Option
    locs_set(2) = .TRUE.
END IF
IF ( PRESENT( dt_Option ) ) THEN
    dtlocs = dt_Option
    dlocs_set(2) = .TRUE.
END IF

IF ( PRESENT( p_Option ) ) THEN
    plocs = p_Option
    locs_set(3) = .TRUE.
END IF
IF ( PRESENT( dp_Option ) ) THEN
    dplocs = dp_Option
    dlocs_set(3) = .TRUE.
END IF


CALL Initialize_Mesh( )





IF ( PRESENT( Solver_Type_Option ) ) THEN
    Solver_Type = Solver_Type_Option
ELSE
    Solver_Type = 2
END IF



IF ( PRESENT(Suffix_Flag_Option) ) THEN


    IF ( Suffix_Flag_Option == "Params") THEN

        WRITE(File_Suffix,'(A,I4.4,A,I2.2,A,I2.2)')"RE",Num_R_Elements,"_D",Degree,"_L",L_Limit

    ELSEIF ( SUffix_Flag_Option == "Frame") THEN
        
        IF ( PRESENT(Frame_Option) ) THEN
            WRITE(File_Suffix,'(I5.5)') Frame_Option
        ELSE
            WRITE(File_Suffix,'(I5.5)') 1
        END IF
    END IF
ELSE
    WRITE(File_Suffix,'(I5.5)') 1
END IF



IF ( PRESENT( Dimensions_Option ) ) THEN
    Domain_Dim = Dimensions_Option
ELSE
    Domain_Dim = 3
END IF

LM_Location => CFA_3D_LM_Map

CALL Allocate_Poseidon_CFA_Variables()
CALL Initialize_Derived()
CALL Initialize_Quadrature()
CALL Initialize_MPI()
CALL Initialize_Tables()


IF ( Solver_Type == 1 ) THEN
    CALL Initialize_NR()
ELSE IF ( Solver_Type == 2 ) THEN
    CALL Initialize_FP(CFA_EQ_Flags_Option)
END IF


CALL Output_Setup_Table()



IF ( Verbose_Flag ) THEN
    PRINT*,"Poseidon Initialization Complete"
END IF

END SUBROUTINE Initialize_Poseidon







 !+201+####################################################################################!
!                                                                                           !
!       Initialize_Poseidon_From_File                                                       !
!                                                                                           !
 !#########################################################################################!
SUBROUTINE Initialize_Poseidon_From_File()




END SUBROUTINE Initialize_Poseidon_From_File






!+301+######################################################################################!
!                                                                                           !
!       OUTPUT_SETUP_TABLE                                                                  !
!                                                                                           !
!###########################################################################################!
SUBROUTINE OUTPUT_SETUP_TABLE( )

INTEGER                                     ::  NUM_JACOBIAN_ELEMENTS
INTEGER                                     ::  NONZEROS
INTEGER                                     ::  BLOCK_NONZEROS
REAL(KIND = idp)                            ::  SPARSITY

LOGICAL                                     ::  EX

1401 FORMAT('------------- POSEIDON PARAMETERS --------------'/)
1402 FORMAT('         ACTIVE DIMENSIONS = ',I12.1)
1403 FORMAT('                    DEGREE = ',I12.1)
1404 FORMAT('                   L_LIMIT = ',I12.1)
1405 FORMAT('              NUM_CFA_VARS = ',I12.1)
1406 FORMAT('                    nPROCS = ',I12.1)
1411 FORMAT('  POSEIDON RADIAL ELEMENTS = ',I12.1)
1412 FORMAT('   POSEIDON THETA ELEMENTS = ',I12.1)
1413 FORMAT('     POSEIDON PHI ELEMENTS = ',I12.1)
1414 FORMAT('          R_COARSEN_FACTOR = ',I12.1)
1415 FORMAT('          T_COARSEN_FACTOR = ',I12.1)
1416 FORMAT('          P_COARSEN_FACTOR = ',I12.1)
1417 FORMAT('               NUM_R_NODES = ',I12.1)
1418 FORMAT('          NUMBER OF SHELLS = ',I12.1)
1419 FORMAT('       NUMBER OF SUBSHELLS = ',I12.1)
1420 FORMAT('       SUBSHELLS PER SHELL = ',I12.1)
1421 FORMAT('      NUM BLOCKS PER SHELL = ',I12.1)
1422 FORMAT('     NUM_R_ELEMS_PER_SHELL = ',I12.1)
1423 FORMAT('     NUM_T_ELEMS_PER_BLOCK = ',I12.1)
1424 FORMAT('     NUM_P_ELEMS_PER_BLOCK = ',I12.1)
1425 FORMAT('                   VAR_DIM = ',I12.1)
1426 FORMAT('             BLOCK_VAR_DIM = ',I12.1)
1427 FORMAT('                  PROB_DIM = ',I12.1)
1428 FORMAT('            BLOCK_PROB_DIM = ',I12.1)
1429 FORMAT('         SUBSHELL_PROB_DIM = ',I12.1)
1430 FORMAT('         NUM_OFF_DIAGONALS = ',I12.1)
1431 FORMAT('     NUM JACOBIAN ELEMENTS = ',I12.1)
1432 FORMAT(' NUM NNZ JACOBIAN ELEMENTS = ',I12.1)
1433 FORMAT('        JACOBIAN SPARSITY  = ',F12.6)
1434 FORMAT('         NUM NNZ PER BLOCK = ',I12.1)
1435 FORMAT('              Inner Radius = ',ES20.12E3)
1436 FORMAT('              Outer Radius = ',ES20.12E3)
1437 FORMAT('      # Radial Quad Points = ',I12.1)
1438 FORMAT('       # Theta Quad Points = ',I12.1)
1439 FORMAT('         # Phi Quad Points = ',I12.1)
1440 FORMAT('        Maximum Iterations = ',I12.1)
1441 FORMAT('      Convergence Criteria = ',ES20.12E3)


IF ( Verbose_Flag ) THEN

    NUM_JACOBIAN_ELEMENTS = PROB_DIM*PROB_DIM

    NONZEROS = NUM_R_ELEMENTS*ELEM_PROB_DIM_SQR          &
             - (NUM_R_ELEMENTS-1)*ULM_LENGTH*ULM_LENGTH

    BLOCK_NONZEROS = NUM_R_ELEMS_PER_BLOCK*ELEM_PROB_DIM_SQR          &
                   - (NUM_R_ELEMS_PER_BLOCK-1)*ULM_LENGTH*ULM_LENGTH

    SPARSITY = REAL(NONZEROS,KIND =idp)/REAL(NUM_JACOBIAN_ELEMENTS, KIND = idp)


    IF  ( myID_Poseidon == 0 )  THEN
        WRITE(*,1401)
        WRITE(*,1402)DOMAIN_DIM
        WRITE(*,1403)DEGREE
        WRITE(*,1404)L_LIMIT
        WRITE(*,1405)NUM_CFA_VARS
        WRITE(*,1406)nPROCS_POSEIDON
        WRITE(*,1411)NUM_R_ELEMENTS
        WRITE(*,1412)NUM_T_ELEMENTS
        WRITE(*,1413)NUM_P_ELEMENTS
        WRITE(*,1414)R_COARSEN_FACTOR
        WRITE(*,1415)T_COARSEN_FACTOR
        WRITE(*,1416)P_COARSEN_FACTOR
        WRITE(*,1417)NUM_R_NODES
        WRITE(*,1418)NUM_SHELLS
        WRITE(*,1419)NUM_SUBSHELLS
        WRITE(*,1420)NUM_SUBSHELLS_PER_SHELL
        WRITE(*,1421)NUM_BLOCKS_PER_SHELL
        WRITE(*,1422)NUM_R_ELEMS_PER_SHELL
        WRITE(*,1423)NUM_T_ELEMS_PER_BLOCK
        WRITE(*,1424)NUM_P_ELEMS_PER_BLOCK
        WRITE(*,1425)VAR_DIM
        WRITE(*,1426)BLOCK_VAR_DIM
        WRITE(*,1427)PROB_DIM
        WRITE(*,1428)BLOCK_PROB_DIM
        WRITE(*,1429)SUBSHELL_PROB_DIM
        WRITE(*,1430)NUM_OFF_DIAGONALS
        WRITE(*,1431)NUM_JACOBIAN_ELEMENTS
        WRITE(*,1432)NONZEROS
        WRITE(*,1433)SPARSITY
        WRITE(*,1434)BLOCK_NONZEROS
        WRITE(*,1435)R_INNER
        WRITE(*,1436)R_OUTER
        WRITE(*,1437)NUM_R_QUAD_POINTS
        WRITE(*,1438)NUM_T_QUAD_POINTS
        WRITE(*,1439)NUM_P_QUAD_POINTS
        WRITE(*,1440)MAX_ITERATIONS
        WRITE(*,1441)CONVERGENCE_CRITERIA

    END IF

    IF ( RUN_REPORT_FILE_ID .NE. -1 ) THEN
        WRITE(RUN_REPORT_FILE_ID,1401)
        WRITE(RUN_REPORT_FILE_ID,1402)DOMAIN_DIM
        WRITE(RUN_REPORT_FILE_ID,1403)DEGREE
        WRITE(RUN_REPORT_FILE_ID,1404)L_LIMIT
        WRITE(RUN_REPORT_FILE_ID,1405)NUM_CFA_VARS
        WRITE(RUN_REPORT_FILE_ID,1406)nPROCS_POSEIDON
        WRITE(RUN_REPORT_FILE_ID,1411)NUM_R_ELEMENTS
        WRITE(RUN_REPORT_FILE_ID,1412)NUM_T_ELEMENTS
        WRITE(RUN_REPORT_FILE_ID,1413)NUM_P_ELEMENTS
        WRITE(RUN_REPORT_FILE_ID,1414)R_COARSEN_FACTOR
        WRITE(RUN_REPORT_FILE_ID,1415)T_COARSEN_FACTOR
        WRITE(RUN_REPORT_FILE_ID,1416)P_COARSEN_FACTOR
        WRITE(RUN_REPORT_FILE_ID,1417)NUM_R_NODES
        WRITE(RUN_REPORT_FILE_ID,1418)NUM_SHELLS
        WRITE(RUN_REPORT_FILE_ID,1419)NUM_SUBSHELLS
        WRITE(RUN_REPORT_FILE_ID,1420)NUM_SUBSHELLS_PER_SHELL
        WRITE(RUN_REPORT_FILE_ID,1421)NUM_BLOCKS_PER_SHELL
        WRITE(RUN_REPORT_FILE_ID,1422)NUM_R_ELEMS_PER_SHELL
        WRITE(RUN_REPORT_FILE_ID,1423)NUM_T_ELEMS_PER_BLOCK
        WRITE(RUN_REPORT_FILE_ID,1424)NUM_P_ELEMS_PER_BLOCK
        WRITE(RUN_REPORT_FILE_ID,1425)VAR_DIM
        WRITE(RUN_REPORT_FILE_ID,1426)BLOCK_VAR_DIM
        WRITE(RUN_REPORT_FILE_ID,1427)PROB_DIM
        WRITE(RUN_REPORT_FILE_ID,1428)BLOCK_PROB_DIM
        WRITE(RUN_REPORT_FILE_ID,1429)SUBSHELL_PROB_DIM
        WRITE(RUN_REPORT_FILE_ID,1430)NUM_OFF_DIAGONALS
        WRITE(RUN_REPORT_FILE_ID,1431)NUM_JACOBIAN_ELEMENTS
        WRITE(RUN_REPORT_FILE_ID,1432)NONZEROS
        WRITE(RUN_REPORT_FILE_ID,1433)SPARSITY
        WRITE(RUN_REPORT_FILE_ID,1434)BLOCK_NONZEROS
        WRITE(RUN_REPORT_FILE_ID,1435)R_INNER
        WRITE(RUN_REPORT_FILE_ID,1436)R_OUTER
        WRITE(RUN_REPORT_FILE_ID,1437)NUM_R_QUAD_POINTS
        WRITE(RUN_REPORT_FILE_ID,1438)NUM_T_QUAD_POINTS
        WRITE(RUN_REPORT_FILE_ID,1439)NUM_P_QUAD_POINTS
        WRITE(RUN_REPORT_FILE_ID,1440)MAX_ITERATIONS
        WRITE(RUN_REPORT_FILE_ID,1441)CONVERGENCE_CRITERIA
        WRITE(RUN_REPORT_FILE_ID,'(/ / / / / / / /)')


    END IF
END IF


END SUBROUTINE OUTPUT_SETUP_TABLE




END MODULE Initialization_Poseidon
