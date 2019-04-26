   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE IO_Functions_Module                                                          !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains the subroutines used to facilitate the output or input of parts    !##!
!##!        of the linear system to/from a file.                                    !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   OUTPUT_COEFFS_TO_FILE                                               !##!
!##!    +102+   INPUT_COEFFS_FROM_FILE                                              !##!
!##!    +103+   OUTPUT_SRC_TO_FILE                                                  !##!
!##!    +104+   OUTPUT_STF_MATRIX_TO_FILE                                           !##!
!##!                                                                                !##!
!##!    +201+   MATRIX_TO_AIJ                                                       !##!
!##!    +202+   MATRIX_TO_CRS                                                       !##!
!##!    +203+   MATRIX_TO_CCS                                                       !##!
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



USE Poseidon_Constants_Module, &
                    ONLY : idp, pi

USE Poseidon_Parameters, &
            ONLY :  DOMAIN_DIM,             &
                    DEGREE,                 &
                    L_LIMIT,                &
                    NUM_R_QUAD_POINTS,      &
                    NUM_T_QUAD_POINTS,      &
                    NUM_P_QUAD_POINTS,      &
                    NUM_R_ELEMS_PER_BLOCK

USE CHIMERA_Parameters,  &
            ONLY :  Analytic_Solution

USE Global_Variables_And_Parameters, &
                    ONLY :  NUM_R_NODES,                                    &
                            NUM_R_ELEMENTS,                                 &
                            RHS_Vector,                                     &
                            Coefficient_Vector,                             &
                            rlocs, tlocs, plocs,                            &
                            R_INNER, R_OUTER,                               &
                            INT_R_LOCATIONS,                                &
                            INT_T_LOCATIONS,                                &
                            INT_P_LOCATIONS,                                &
                            VAR_DIM,                                        &
                            NUM_OFF_DIAGONALS,                              &
                            Num_Timer_Calls,                                &
                            Time_Table
 

USE Additional_Functions_Module, &
                    ONLY :  Map_From_x_Space,                               &
                            Initialize_LGL_Quadrature_Locations,            &
                            Lagrange_Poly,                                  &
                            Lagrange_Poly_Deriv,                            &
                            Lagrange_Second_Deriv,                          &
                            CFA_Matrix_Map



!USE CFA_3D_Master_Build_Module, &
!                    ONLY :  Calc_3D_Values_At_Location




IMPLICIT NONE


!*F&S*==========================================!
!                                               !
!           Functions & Subroutines             !
!                                               !
!===============================================!
CONTAINS












!########################################################################################!
!
!
!########################################################################################!

SUBROUTINE WRITE_NR_MATRIX( Matrix, I_DIM, J_DIM, Cur_Iteration )


COMPLEX(KIND = idp), DIMENSION(0:I_DIM, 0:J_DIM )           :: Matrix
INTEGER, INTENT(IN)                                         :: I_DIM, J_DIM

INTEGER, OPTIONAL                                           :: Cur_Iteration




INTEGER                                                     :: i, j



CHARACTER(LEN = 29)                                          ::  FILENAME
CHARACTER(LEN = 37)                                          ::  FILENAMEA
CHARACTER(LEN = 23)                                          ::  FILE_DIR="OUTPUT/PETSC/NEW_M_ITER"
CHARACTER(LEN = 4)                                           ::  FILE_EXT=".out"






IF ( PRESENT(Cur_Iteration) )  THEN

    WRITE(FILENAMEA,'(A,A,I3.3,A)') FILE_DIR, "_", CUR_ITERATION, FILE_EXT
    OPEN(unit = 1, file = FILENAMEA)


ELSE



    WRITE(FILENAME,'(A,A)') FILE_DIR,FILE_EXT
    OPEN(unit = 1, file = FILENAME)



END IF





DO i = 0, J_DIM

    WRITE(1,*) REAL( Matrix(:, i), KIND = 4 )

END DO

CLOSE( unit = 1)






END SUBROUTINE WRITE_NR_MATRIX












!########################################################################################!
!
!
!########################################################################################!

SUBROUTINE WRITE_NR_VECTOR( Vector, I_DIM, Cur_Iteration )


COMPLEX(KIND = idp), DIMENSION(0:I_DIM )                    :: Vector
INTEGER, INTENT(IN)                                         :: I_DIM

INTEGER, OPTIONAL                                           :: Cur_Iteration




INTEGER                                                     :: i, j



CHARACTER(LEN = 27)                                          ::  FILENAME
CHARACTER(LEN = 38)                                          ::  FILENAMEA
CHARACTER(LEN = 29)                                          ::  FILE_DIR="OUTPUT/CFA/COMPARE/NEW_V_ITER"
CHARACTER(LEN = 4)                                           ::  FILE_EXT=".out"






IF ( PRESENT(Cur_Iteration) )  THEN

    WRITE(FILENAMEA,'(A,A,I3.3,A)') FILE_DIR, "_", CUR_ITERATION, FILE_EXT
    OPEN(unit = 1, file = FILENAMEA)


ELSE



    WRITE(FILENAME,'(A,A)') FILE_DIR,FILE_EXT
    OPEN(unit = 1, file = FILENAME)



END IF





DO i = 0, I_DIM

    WRITE(1,*) REAL( Vector(i), KIND = idp )

END DO

CLOSE( unit = 1)






END SUBROUTINE WRITE_NR_VECTOR















SUBROUTINE OUTPUT_COMPARE_MAP


                                                                            ! MAP AXIS NUMBERS
INTEGER, PARAMETER                                      :: MAP_AXIS_A = 1   ! MUST BE DIFFERENT
INTEGER, PARAMETER                                      :: MAP_AXIS_B = 3   ! 1 = X
                                                                            ! 2 = Y
                                                                            ! 3 = Z

INTEGER, PARAMETER                                      ::  NUM_X_SAMPLES = 256
INTEGER, PARAMETER                                      ::  NUM_Y_SAMPLES = 256

INTEGER                                                 ::  i, j

REAL(KIND = idp)                                        ::  deltar, r, theta, phi, &
                                                            deltax, deltay, map_x, map_y

REAL(KIND = idp)                                        ::  Local_Psi,          &
                                                            Local_AlphaPsi,     &
                                                            Local_Beta1,        &
                                                            Local_Beta2,        &
                                                            Local_Beta3,        &
                                                            Local_Solution



REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)              :: TMP_PSI,        &
                                                            TMP_ALPHAPSI,   &
                                                            TMP_BETA1,      &
                                                            TMP_BETA2,      &
                                                            TMP_BETA3,      &
                                                            TMP_SOLUTION










IF (MAP_AXIS_A .EQ. MAP_AXIS_B) THEN

    PRINT*,"Mapping axes must be different"

ELSE


    ALLOCATE(TMP_PSI(0:NUM_X_SAMPLES), TMP_ALPHAPSI(0:NUM_X_SAMPLES), TMP_BETA1(0:NUM_X_SAMPLES),   &
             TMP_BETA2(0:NUM_X_SAMPLES), TMP_BETA3(0:NUM_X_SAMPLES), TMP_SOLUTION(0:NUM_X_SAMPLES) )


    deltax = (1.40_idp*R_OUTER)/REAL(NUM_X_SAMPLES-1, KIND = idp)
    deltay = (1.40_idp*R_OUTER)/REAL(NUM_Y_SAMPLES-1, KIND = idp)



    OPEN(unit = 13, file = 'OUTPUT/COMPARE_MAP/COMPARE_MAP.PROP')

    WRITE(13,*) MAP_AXIS_A, MAP_AXIS_B, deltax, deltay

    CLOSE(unit = 13)



    OPEN(unit = 41, file = 'OUTPUT/COMPARE_MAP/COMPARE_MAP_MESH.DATA')
    OPEN(unit = 42, file = 'OUTPUT/COMPARE_MAP/COMPARE_MAP_PSI.DATA')
    OPEN(unit = 43, file = 'OUTPUT/COMPARE_MAP/COMPARE_MAP_ALPHAPSI.DATA')
    OPEN(unit = 44, file = 'OUTPUT/COMPARE_MAP/COMPARE_MAP_BETA1.DATA')
    OPEN(unit = 45, file = 'OUTPUT/COMPARE_MAP/COMPARE_MAP_BETA2.DATA')
    OPEN(unit = 46, file = 'OUTPUT/COMPARE_MAP/COMPARE_MAP_BETA3.DATA')
    OPEN(unit = 47, file = 'OUTPUT/COMPARE_MAP/COMPARE_MAP_POTENTIAL.DATA')



    DO i = 0,NUM_X_SAMPLES
        DO j = 0,NUM_Y_SAMPLES


            map_x = -.70_idp*R_OUTER + (i) * deltax
            map_y = -.70_idp*R_OUTER + (j) * deltay

            r = sqrt(map_x*map_x + map_y*map_y)


            IF (MIN(MAP_AXIS_A, MAP_AXIS_B) .EQ. 1) THEN

                IF( MAX(MAP_AXIS_A,MAP_AXIS_B) .EQ. 2) THEN


                    theta = pi/2.0_idp
                    phi = atan2(map_y,map_x)


                ELSE IF (MAX(MAP_AXIS_A,MAP_AXIS_B) .EQ. 3) THEN

                    theta = acos(ABS(map_x)/r)

                    IF (map_x .LE. 0) THEN

                        phi = pi

                    ELSE

                        phi = 0.0_idp

                    END IF

                END IF

            ELSE IF (MIN(MAP_AXIS_A,MAP_AXIS_B) .EQ. 2) THEN

                theta = acos(abs(map_x)/r)

                IF (map_x .LE. 0) THEN

                    phi = 1.50_idp * pi

                ELSE

                    phi = pi/2.0_idp

                END IF



            END IF


!            CALL Calc_3D_Values_At_Location( r, theta, phi, Local_Psi,         &
!                                            Local_AlphaPsi, Local_Beta1,      &
!                                            Local_Beta2, Local_Beta3          )



 !           TMP_PSI(j) = Local_Psi
 !           TMP_ALPHAPSI(j) = Local_AlphaPsi
 !           TMP_Beta1(j) = Local_Beta1
 !           TMP_BETA2(j) = Local_Beta2
 !           TMP_BETA3(j) = Local_Beta3

            TMP_SOLUTION(j) = Analytic_Solution(r,theta,phi)




            !PRINT*,r,theta,phi,REALPART(potential), ABS(TMP-REALPART(potential))
            !print*,r,REALPART(potential),TMP, ABS(TMP-REALPART(potential))

        END DO


  !      WRITE(41,*) map_x
  !      WRITE(42,*) TMP_PSI
  !      WRITE(43,*) TMP_ALPHAPSI
  !      WRITE(44,*) TMP_BETA1
  !      WRITE(45,*) TMP_BETA2
  !      WRITE(46,*) TMP_BETA3
  !      WRITE(47,*) TMP_SOLUTION


    END DO

    CLOSE(unit = 42)


END IF



END SUBROUTINE OUTPUT_COMPARE_MAP



!########################################################################################!
!
!
!########################################################################################!

SUBROUTINE WRITE_BDDtoFULL_MATRIX( Matrix, J_DIM, I_DIM, Cur_Iteration )


COMPLEX(KIND = idp), DIMENSION(0:J_DIM, 0:I_DIM )           :: Matrix
INTEGER, INTENT(IN)                                         :: I_DIM, J_DIM

INTEGER, OPTIONAL                                           :: Cur_Iteration




INTEGER                                                     ::  j,              &
                                                                cur_i_loc,      &
                                                                cur_j_loc,      &
                                                                cur_jp_loc,     &
                                                                re, te, pe,     &
                                                                l, m, d,        &
                                                                lp, mp, dp,     &
                                                                u, up


INTEGER                                                     ::  L_SHIFT,        &
                                                                Sparse_SHIFT


COMPLEX(KIND = idp), DIMENSION(0:I_DIM)                     ::  TMP_VECTOR


CHARACTER(LEN = 29)                                         ::  FILENAME
CHARACTER(LEN = 37)                                         ::  FILENAMEA
CHARACTER(LEN = 23)                                         ::  FILE_DIR="OUTPUT/PETSC/NEW_M_ITER"
CHARACTER(LEN = 4)                                          ::  FILE_EXT=".out"


REAL(KIND = 4), DIMENSION(0:I_DIM, 0:I_DIM)                 ::  TMP_MATRIX



IF ( PRESENT(Cur_Iteration) )  THEN

    WRITE(FILENAMEA,'(A,A,I3.3,A)') FILE_DIR, "_", CUR_ITERATION, FILE_EXT
    OPEN(unit = 1, file = FILENAMEA)

ELSE

    WRITE(FILENAME,'(A,A)') FILE_DIR,FILE_EXT
    OPEN(unit = 1, file = FILENAME)

END IF



L_SHIFT = 5*(L_LIMIT+1)*(L_LIMIT+1)


TMP_MATRIX = 0.0_idp
re = 0





Print*,"Output Matrix Size "








DO re = 0, NUM_R_ELEMENTS - 1

    do d = 0, DEGREE

        do l = 0,L_LIMIT

            do m = -l,l

                DO u = 0,4

                    Do dp = 0, DEGREE

                        DO lp = 0,L_LIMIT

                            DO mp = -lp,lp

                                DO up = 0,4

                                    cur_i_loc = (re*DEGREE + dp)*L_SHIFT + 5*(lp*(lp+1) + mp) + up

                                    cur_j_loc = 2*NUM_OFF_DIAGONALS                                     &
                                                + ( re*DEGREE + d ) * L_SHIFT + 5*(l*(l+1) + m) + u     &
                                                - cur_i_loc

                                    cur_jp_loc = (re*DEGREE + d ) * L_SHIFT + 5*(l*(l+1) + m) + u


                                    TMP_MATRIX(cur_jp_loc, cur_i_loc) = matrix(cur_j_loc, cur_i_loc)



                                END DO ! up Loop

                            END DO ! u Loop

                        END DO ! mp loop

                    END DO ! lp Loop

                END DO ! dp Loop


            END DO ! m Loop

        END DO ! l loop

    END DO ! d loop

END DO ! re




DO j = 0, I_DIM

    WRITE(1,*) REAL( TMP_MATRIX(j, :), KIND = 4 )

END DO

CLOSE( unit = 1)





CLOSE( unit = 1)






END SUBROUTINE WRITE_BDDtoFULL_MATRIX





!########################################################################################!
!
!
!########################################################################################!

SUBROUTINE WRITE_BDDtoFULL_MATRIX_Block( Matrix, J_DIM, I_DIM, Cur_Iteration )


COMPLEX(KIND = idp), DIMENSION(0:J_DIM, 0:I_DIM )           :: Matrix
INTEGER, INTENT(IN)                                         :: I_DIM, J_DIM

INTEGER, OPTIONAL                                           :: Cur_Iteration




INTEGER                                                     ::  j,              &
                                                                cur_i_loc,      &
                                                                cur_j_loc,      &
                                                                cur_jp_loc,     &
                                                                re, te, pe,     &
                                                                l, m, d,        &
                                                                lp, mp, dp,     &
                                                                u, up


INTEGER                                                     ::  L_SHIFT,        &
                                                                Sparse_SHIFT


COMPLEX(KIND = idp), DIMENSION(0:I_DIM)                     ::  TMP_VECTOR


CHARACTER(LEN = 29)                                         ::  FILENAME
CHARACTER(LEN = 37)                                         ::  FILENAMEA
CHARACTER(LEN = 23)                                         ::  FILE_DIR="OUTPUT/PETSC/NEW_M_ITER"
CHARACTER(LEN = 4)                                          ::  FILE_EXT=".out"


REAL(KIND = 4), DIMENSION(0:I_DIM, 0:I_DIM)                 ::  TMP_MATRIX



IF ( PRESENT(Cur_Iteration) )  THEN

    WRITE(FILENAMEA,'(A,A,I3.3,A)') FILE_DIR, "_", CUR_ITERATION, FILE_EXT
    OPEN(unit = 1, file = FILENAMEA)

ELSE

    WRITE(FILENAME,'(A,A)') FILE_DIR,FILE_EXT
    OPEN(unit = 1, file = FILENAME)

END IF



L_SHIFT = 5*(L_LIMIT+1)*(L_LIMIT+1)


TMP_MATRIX = 0.0_idp
re = 0














DO re = 0, NUM_R_ELEMS_PER_BLOCK - 1

    do d = 0, DEGREE

        do l = 0,L_LIMIT

            do m = -l,l

                DO u = 0,4

                    Do dp = 0, DEGREE

                        DO lp = 0,L_LIMIT

                            DO mp = -lp,lp

                                DO up = 0,4

                                    cur_i_loc = (re*DEGREE + dp)*L_SHIFT + 5*(lp*(lp+1) + mp) + up

                                    cur_j_loc = NUM_OFF_DIAGONALS                                       &
                                                + ( re*DEGREE + d ) * L_SHIFT + 5*(l*(l+1) + m) + u     &
                                                - cur_i_loc

                                    cur_jp_loc = (re*DEGREE + d ) * L_SHIFT + 5*(l*(l+1) + m) + u


                                    TMP_MATRIX(cur_jp_loc, cur_i_loc) = matrix(cur_j_loc, cur_i_loc)



                                END DO ! up Loop

                            END DO ! u Loop

                        END DO ! mp loop

                    END DO ! lp Loop

                END DO ! dp Loop


            END DO ! m Loop

        END DO ! l loop

    END DO ! d loop

END DO ! re




DO j = 0, I_DIM

    WRITE(1,*) REAL( TMP_MATRIX(j, :), KIND = 4 )

END DO

CLOSE( unit = 1)





CLOSE( unit = 1)






END SUBROUTINE WRITE_BDDtoFULL_MATRIX_Block

























!########################################################################################!
!
!
!########################################################################################!
SUBROUTINE CLOCK_IN( Time, Ident )


REAL(KIND = idp), INTENT(IN)                 ::   Time
INTEGER, INTENT(IN)                          ::   Ident


Time_Table(Ident) = Time

END SUBROUTINE CLOCK_IN









!########################################################################################!
!
!
!########################################################################################!
SUBROUTINE PRINT_TIMETABLE(Ident)


INTEGER,INTENT(IN)           :: Ident

PRINT*,"=========================================================================="
PRINT*," "
PRINT*,"                            Time Table for Process",Ident
PRINT*," "
PRINT*,"=========================================================================="
PRINT*,"                    Initialize Time : ",Time_Table(1)
PRINT*," Input/Communicate Source Data Time : ",Time_Table(2)
PRINT*,"     Input Boundary Conditions Time : ",Time_Table(3)
PRINT*,"        CFA_3D_Apply_BCs_Part1 Time : ",Time_Table(4)
PRINT*,"----------------------------------------------------------"
PRINT*," ||     Calc_3D_Current_Values Time : ",Time_Table(5)
PRINT*," ||    CREATE_3D_SubJcbn_Terms Time : ",Time_Table(6)
PRINT*," ||       CREATE_3D_RHS_VECTOR Time : ",Time_Table(7)
PRINT*,"\  /     CREATE_3D_JCBN_MATRIX Time : ",Time_Table(8)
PRINT*,"-\/ ------------------------------------------------------"
PRINT*,"CREATE_3D_NONLAPLACIAN_STF_MAT Time : ",Time_Table(9)
PRINT*,"REDUCE_3D_NONLAPLACIAN_STF_MAT Time : ",Time_Table(10)
PRINT*,"FINISH_3D_NONLAPLACIAN_STF_MAT Time : ",Time_Table(11)
PRINT*,"          FINISH_3D_RHS_VECTOR Time : ",Time_Table(12)
PRINT*,"        CFA_3D_Apply_BCs_Part2 Time : ",Time_Table(13)
PRINT*,"                    CFA_Solver Time : ",Time_Table(14)
PRINT*,"        CFA_Coefficient_Update Time : ",Time_Table(15)
PRINT*,"   CFA_Coefficient_Share_PETSc Time : ",Time_Table(16)
PRINT*,"         CFA_Convergence_Check Time : ",Time_Table(17)
PRINT*,"               Total Iteration Time : ",Time_Table(18)
PRINT*,"             Poseidon_Dist_Sol Time : ",Time_Table(19)
PRINT*,"=========================================================================="






END SUBROUTINE PRINT_TIMETABLE
































END MODULE IO_Functions_Module
