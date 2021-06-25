   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
PROGRAM Beta_Mapping                                                                !##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!

USE Poseidon_Kinds_Module, &
                ONLY :  idp

USE Units_Module, &
                ONLY :  C_Square,       &
                        Set_Units,      &
                        Shift_Units,    &
                        Centimeter,     &
                        Second,         &
                        Gram

USE Initialization_Poseidon, &
                ONLY :  Initialize_Poseidon


USE Poseidon_Parameters, &
                ONLY :  DEGREE,                     &
                        L_LIMIT

USE Variables_IO, &
                ONLY :  Write_Results_Flag,     &
                        Write_Results_R_Samps,  &
                        Write_Results_T_Samps

USE Variables_MPI, &
                ONLY :  ierr

USE Variables_FP, &
                ONLY :  FP_Coeff_Vector,        &
                        FP_Coeff_Vector_Beta,   &
                        FP_Source_Vector_Beta,  &
                        Laplace_Matrix_Beta,    &
                        Beta_Diagonals,         &
                        Beta_MVL_Banded,        &
                        Beta_Bandwidth,         &
                        Matrix_Format,          &
                        Beta_Factorized_Flag

USE FP_Beta_Banded, &
                ONLY :  Factorize_Beta_Banded,  &
                        Jacobi_PC_MVL_Banded_Vector, &
                        Jacobi_PC_MVL_Banded,       &
                        DIRICHLET_BC_Beta_Banded_Mat


USE Variables_Derived,  &
                ONLY :  Beta_Prob_Dim,          &
                        LM_Length

USE Variables_Functions, &
                ONLY :  Potential_Solution

USE Functions_FP, &
                ONLY :  Laplace_Test_Sol

USE Functions_Mesh, &
                ONLY :  Create_3D_Mesh

USE Functions_Matrix, &
                ONLY :  MVMULT_FULL

USE Poseidon_IO_Module, &
                ONLY :  Open_Run_Report_File,       &
                        Output_Run_Report,        &
                        Output_Final_Results

USE Poseidon_Main_Module, &
                ONLY :  Poseidon_CFA_Set_Uniform_Boundary_Conditions,   &
                        Poseidon_Close

USE FP_System_Solvers_Module , &
                ONLY :  Solve_FP_System,        &
                        Solve_FP_System_Beta

USE IO_Print_Results, &
                ONLY :  Print_Results

USE FP_Coeffs_ReadWrite, &
                ONLY :  Output_FP_Coeffs,       &
                        ReadIn_FP_Coeffs


USE FP_Functions_BC,  &
                ONLY :  Dirichlet_BC_Beta,              &
                        DIRICHLET_BC_Beta_Banded

USE Linear_Solvers_And_Preconditioners, &
                ONLY :  PRECOND_CONJ_GRAD_CCS,          &
                        JACOBI_CONDITIONING,            &
                        Jacobi_Conditioning_Beta

USE FP_Functions_Mapping, &
                ONLY :  FP_Beta_Array_Map


USE MPI


IMPLICIT NONE

!                                       !
!   Poseidon Initialization Variables   !
!                                       !


INTEGER                                                 ::  FEM_Degree_Input
INTEGER                                                 ::  L_Limit_Input

INTEGER                                                 ::  Dimension_Input

INTEGER                                                 ::  Mesh_Type
INTEGER                                                 ::  Solver_Type

CHARACTER(LEN = 1)                                      ::  Units_Input

INTEGER, DIMENSION(3)                                   ::  NE
INTEGER, DIMENSION(3)                                   ::  NQ
REAL(idp), DIMENSION(2)                                 ::  Domain_Edge

REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  dx_c, dy_c, dz_c
REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  x_e, y_e, z_e
REAL(idp), DIMENSION(:), ALLOCATABLE                    ::  x_c, y_c, z_c

LOGICAL                                                 ::  Verbose

INTEGER,   DIMENSION(5)                                 ::  CFA_EQs
REAL(idp)                                               ::  Shift_Vector_BC
CHARACTER(LEN=1), DIMENSION(1:5)                        ::  INNER_BC_TYPES, OUTER_BC_TYPES
REAL(idp), DIMENSION(1:5)                               ::  INNER_BC_VALUES, OUTER_BC_VALUES

CHARACTER(LEN=10)                                       ::  Suffix_Input
CHARACTER(LEN=1)                                        ::  Suffix_Tail

INTEGER                                                 ::  RE_Index, RE_Index_Min, RE_Index_Max
INTEGER                                                 ::  Degree_Min, Degree_Max
INTEGER                                                 ::  L_Limit_Min, L_Limit_Max

INTEGER                                                 ::  Num_RE
INTEGER, DIMENSION(1:9)                                 ::  RE_Table
CHARACTER(LEN=1), DIMENSION(1:7)                        ::  Letter_Table


COMPLEX(idp), DIMENSION(:), ALLOCATABLE                 ::  Output
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:)                              ::  WORK_VEC
COMPLEX(KIND = idp), ALLOCATABLE, DIMENSION(:,:)                            ::  WORK_MAT



INTEGER                                                 ::  RE, D, ui, LM, here, ldim
INTEGER                                                 ::  Col, dp, l, lp, row, uj


CALL MPI_INIT(ierr)

!!                                       !!
!!   Poseidon Initialization Variables   !!
!!                                       !!

Units_Input         = "U"
Solver_Type         = 2
Suffix_Input        = "Params"

Letter_Table        = (/ "A","B","C","D","E","F","G" /)
RE_Table            = (/ 3, 80, 160, 240, 320, 400, 600, 800, 1000 /)

Dimension_Input     = 3

Mesh_Type           = 1
Domain_Edge(1)      = 1.0_idp   ! Inner Radius
Domain_Edge(2)      = 1E9_idp   ! Outer Radius

RE_Index_Min        = 3
RE_Index_Max        = 3

Degree_Min          = 1
Degree_Max          = 1

L_Limit_Min         = 1
L_Limit_Max         = 1

!Verbose             = .FALSE.
Verbose             = .TRUE.
CFA_Eqs             = [0, 0, 1, 0, 0]

NQ(1) = 10        ! Number of Radial Quadrature Points
NQ(2) = 15        ! Number of Theta Quadrature Points
NQ(3) = 15       ! Number of Phi Quadrature Points


Write_Results_R_Samps = 256
Write_Results_T_Samps = 256

CALL Set_Units(Units_Input)
Domain_Edge = Domain_Edge*Centimeter

DO RE_Index = RE_Index_Min, RE_Index_Max
    DO FEM_Degree_Input = Degree_Min, Degree_Max
        DO L_Limit_Input = L_Limit_Min, L_Limit_Max


            

            PRINT*,"Iteration Starting"
            WRITE(*,'(I4.4,A,I2.2,A,I2.2,A)')RE_Table(RE_Index)," (",RE_Index," of ",RE_Index_Max,")"
            WRITE(*,'(A,I2.2,A,I2.2)')"Degree: ",FEM_Degree_Input," of ",Degree_Max
            WRITE(*,'(A,I2.2,A,I2.2)')"L_Limit: ",L_Limit_Input," of ",L_Limit_Max

            NE(1) = RE_Table(RE_Index)      ! Number of Radial Elements
            NE(2) = 1                       ! Number of Theta Elements
            NE(3) = 1                       ! Number of Phi Elements

            Suffix_Tail = Letter_Table(1)


            ALLOCATE( x_e(0:NE(1)), y_e(0:NE(2)), z_e(0:NE(3)) )
            ALLOCATE( x_c(1:NE(1)), y_c(1:NE(2)), z_c(1:NE(3)) )
            ALLOCATE( dx_c(1:NE(1)), dy_c(1:NE(2)), dz_c(1:NE(3)) )



            CALL Open_Run_Report_File()

            CALL Create_3D_Mesh( Mesh_Type,         &
                                 Domain_Edge(1),    &
                                 Domain_Edge(2),    &
                                 NE(1),             &
                                 NE(2),             &
                                 NE(3),             &
                                 x_e, x_c, dx_c,    &
                                 y_e, y_c, dy_c,    &
                                 z_e, z_c, dz_c,    &
                                 Zoom = 1.032034864238313_idp )



            CALL Initialize_Poseidon &
                (   Dimensions_Option       = Dimension_Input,  &
                    FEM_Degree_Option       = FEM_Degree_Input, &
                    L_Limit_Option          = L_Limit_Input,    &
                    Units_Option            = Units_Input,      &
                    Domain_Edge_Option      = Domain_Edge,      &
                    NE_Option               = NE,               &
                    NQ_Option               = NQ,               &
                    r_Option                = x_e,              &
                    t_Option                = y_e,              &
                    p_Option                = z_e,              &
            !        dr_Option               = dx_c,             &
            !        dt_Option               = dy_c,             &
            !        dp_Option               = dz_c              &
                    Suffix_Flag_Option       = Suffix_Input,    &
                    Suffix_Tail_Option       = Suffix_Tail,     &
                    Solver_Type_Option       = Solver_Type,     &
                    CFA_Eq_Flags_Option      = CFA_Eqs,         &
                    Verbose_Option           = Verbose          )




            INNER_BC_TYPES = (/"N", "N","N","N","N"/)
            OUTER_BC_TYPES = (/"D", "D","D","D","D"/)

            Shift_Vector_BC = -1.0E2_idp
            OUTER_BC_VALUES = (/0.0_idp, 0.0_idp, Shift_Vector_BC, 0.0_idp, 0.0_idp /)

            CALL Poseidon_CFA_Set_Uniform_Boundary_Conditions("I", INNER_BC_TYPES, INNER_BC_VALUES)
            CALL Poseidon_CFA_Set_Uniform_Boundary_Conditions("O", OUTER_BC_TYPES, OUTER_BC_VALUES)

            


            ALLOCATE( Output(1:Beta_Prob_Dim) )
            Output = 0.0_idp
    
            CALL ReadIn_FP_Coeffs( NE(1), FEM_Degree_Input, 1, Suffix_Tail )



            ! Multiply Shift Matrix and Coeff Vector
            IF ( Matrix_Format == 'CCS' ) THEN

                ALLOCATE( WORK_VEC( 1:Beta_Prob_Dim ) )
                Work_Vec = 0.0_idp


!                IF ( .NOT. Beta_Factorized_Flag ) THEN
!
!                    CALL Factorize_Beta_Banded()
!
!                END IF

                CALL DIRICHLET_BC_Beta_Banded_Mat()


                CALL Jacobi_PC_MVL_Banded()



 
!                PRINT*,"Work_Vec A"
!                PRINT*,Work_Vec
                CALL DIRICHLET_BC_Beta_Banded(Beta_Prob_Dim, Work_Vec )


!                PRINT*,"Work_Vec B"
!                PRINT*,Work_Vec

                CALL Jacobi_PC_MVL_Banded_Vector( Work_Vec )



     
!                PRINT*,"Work_Vec C"
!                PRINT*,Work_Vec

                CALL ZGBMV('N',                     &
                            Beta_Prob_Dim,          &
                            Beta_Prob_Dim,          &
                            Beta_Diagonals,       &
                            2*Beta_Diagonals,         &
                            1.0_idp,                &
                            Beta_MVL_Banded,        &
                            3*Beta_Diagonals+1,     &
                            FP_Coeff_Vector_Beta,   &
                            1,                      &
                            0.0_idp,                &
                            Output,                 &
                            1                       )

!                PRINT*,Beta_MVL_Banded
!                PRINT*," "
!                PRINT*,FP_Coeff_Vector_Beta
!                PRINT*," "
!                PRINT*,FP_Source_Vector_Beta(:)
!                PRINT*," "


!                PRINT*,"ZGBMV"
!                PRINT*,Output - Work_Vec
!                PRINT*," "




            ELSE IF ( Matrix_Format == 'Full' ) THEN



                ALLOCATE( WORK_VEC( 1:Beta_Prob_Dim ) )
                ALLOCATE( WORK_MAT( 1:Beta_Prob_Dim, 1:Beta_Prob_Dim ) )

                WORK_MAT(:,:) = Laplace_Matrix_Beta(:,:)
                WORK_VEC(:)   = 0.0_idp

!                DO re = 1,Beta_Prob_Dim
!                DO ui = 1,Beta_Prob_Dim
!                    PRINT*,ui,re,Work_Mat(ui,re)
!                END DO
!                END DO



!                PRINT*,"Work_Vec A"
!                PRINT*,Work_Vec

                CALL DIRICHLET_BC_Beta(WORK_MAT, WORK_VEC)




!                PRINT*,"Work_Vec B"
!                PRINT*,Work_Vec

                CALL JACOBI_CONDITIONING_Beta(WORK_MAT, WORK_VEC, Beta_Prob_Dim, Beta_Prob_Dim)





!                PRINT*,"Work_Vec C"
!                PRINT*,Work_Vec

!                CALL ZGEMV('N',                     &
!                            Beta_Prob_Dim,          &
!                            Beta_Prob_Dim,          &
!                            1.0_idp,                &
!                            Laplace_Matrix_Beta,    &
!                            Beta_Prob_Dim,          &
!                            FP_Coeff_Vector_Beta,   &
!                            1,                      &
!                            0.0_idp,                &
!                            Output,                 &
!                            1                       )
!
!                PRINT*,"ZGEMV"
!                PRINT*,Output
!                PRINT*," "

                
                
                Output = MVMULT_FULL( Work_Mat,            &
                                      FP_Coeff_Vector_Beta(:),        &
                                      Beta_Prob_Dim, Beta_Prob_Dim    )
!                PRINT*,"MVMULT_FULL"
!                PRINT*,Output - Work_Vec
!                PRINT*," "
                

            END IF




            ! Analyze output
            PRINT*," "


!            
!            DO ui = 1,3
!            DO RE = 0,NE(1)-1
!            DO d  = 0,FEM_Degree_Input
!            DO lm = 0,LM_Length-1
!
!                Here = (RE*FEM_Degree_Input + d)*3*LM_Length     &
!                     +  (ui - 1)*LM_Length                       &
!                     +  lm + 1
!
!!                PRINT*,ui, RE+1,d,lm,Here,Output(Here)-Work_Vec(Here)
!!                PRINT*,Here,Output(Here),Work_Vec(Here),Output(Here)-Work_Vec(Here)
!                PRINT*,ui,RE,d,lm,FP_Coeff_Vector_Beta(Here)
!
!            END DO ! lm
!            END DO ! ui
!            END DO ! d
!            END DO ! RE
            PRINT*," "
            PRINT*,"Maximum value of the output vector is ",MAXVAL(abs(Output-Work_Vec)),        &
                   " at ",maxloc(abs(Output))



            CALL Poseidon_Close()



            DEALLOCATE( x_e, y_e, z_e )
            DEALLOCATE( x_c, y_c, z_c )
            DEALLOCATE( dx_c, dy_c, dz_c )
            DEALLOCATE( Output )

        END DO ! L_Limit
    END DO ! Degree
END DO ! RE_Index


WRITE(*,'(//A18//)')"DING! You're Done!"

END PROGRAM Beta_Mapping

