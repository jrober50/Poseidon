   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_BC_Module                                                           !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   CFA_3D_Apply_BCs_Part1                                              !##!
!##!    +102+   CFA_3D_Apply_BCs_Part2                                              !##!
!##!    +103+   CFA_3D_Dirichlet_BCs_Part1                                          !##!
!##!    +104+   CFA_3D_Dirichlet_BCs_Part2                                          !##!
!##!    +105+   CFA_3D_Neumann_BCs                                                  !##!
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
            ONLY :  pi

USE Units_Module, &
            ONLY :  C_Square,           &
                    GR_Source_Scalar

USE Variables_Functions, &
            ONLY :  Potential_Solution,     &
                    Matrix_Location


USE Poseidon_Parameters, &
            ONLY :  DEGREE,                     &
                    L_LIMIT,                    &
                    NUM_CFA_VARS

USE Variables_MPI,  &
            ONLY :  Num_SubShells,              &
                    Num_R_Elems_Per_Block,      &
                    myID_PETSc

USE Variables_Mesh, &
            ONLY :  Num_R_Elements

USE Variables_NR, &
            ONLY :  Coefficient_Vector,         &
                    BLOCK_ELEM_STF_MATVEC,      &
                    Block_RHS_Vector

USE Variables_BC, &
            ONLY :  INNER_CFA_BC_VALUES,        &
                    OUTER_CFA_BC_VALUES,        &
                    INNER_CFA_BC_TYPE,          &
                    OUTER_CFA_BC_TYPE
            
USE Variables_Derived, &
            ONLY :  Prob_Dim,                   &
                    Elem_Prob_Dim
                    
USE Variables_Tables, &
            ONLY :  M_Values

USE Functions_Math, &
            ONLY :  Lagrange_Poly

USE Functions_Mapping, &
            ONLY :  Map_To_X_Space

USE Functions_Quadrature, &
            ONLY :  Initialize_LG_Quadrature,       &
                    Initialize_LGL_Quadrature


IMPLICIT NONE


CONTAINS


!+101+###########################################################################!
!                                                                                !
!                  CFA_3D_Apply_BCs                                                  !
!                                                                                !
!################################################################################!
SUBROUTINE CFA_3D_Apply_BCs_Part1(  )

CALL CFA_3D_Dirichlet_BCs_Part1()
CALL CFA_3D_Neumann_BCs()

END SUBROUTINE CFA_3D_Apply_BCs_Part1



!+102+###########################################################################!
!                                                                                !
!                  CFA_3D_Apply_BCs                                                  !
!                                                                                !
!################################################################################!
SUBROUTINE CFA_3D_Apply_BCs_Part2(  )

CALL CFA_3D_Dirichlet_BCs_Part2()
CALL CFA_3D_Neumann_BCs()

END SUBROUTINE CFA_3D_Apply_BCs_Part2


!+103+###########################################################################!
!                                                                                !
!                  CFA_3D_Dirichlet_BCs                                                  !
!                                                                                !
!################################################################################!
SUBROUTINE CFA_3D_Dirichlet_BCs_Part1( )

INTEGER                                             ::  i, l, m, ui, d

INTEGER                                             ::  Value_Location

!*!
!*!     Steps 1
!*!
DO ui = 1,NUM_CFA_VARS


    !*!
    !*!     Innner BCs
    !*!
    IF ( INNER_CFA_BC_TYPE(ui) == 'D' ) THEN

        DO l = 0,L_LIMIT
            DO m = -M_VALUES(l),M_VALUES(l)

                Value_Location =  Matrix_Location( ui, l, m, 0, 0 )
                Coefficient_Vector( Value_Location ) = 0.0_idp

            END DO
        END DO

        Value_Location =  Matrix_Location( ui, 0, 0, 0, 0 )
        Coefficient_Vector( Value_Location ) = 2.0_idp * sqrt(pi) * INNER_CFA_BC_VALUES(ui)
    END IF





    !*!
    !*!     Outer BCs
    !*!
    IF ( OUTER_CFA_BC_TYPE(ui) == 'D' ) THEN

        DO l = 1,L_LIMIT
            DO m = -M_VALUES(l),M_VALUES(l)

                Value_Location =  Matrix_Location( ui, l, m, NUM_R_ELEMENTS-1, DEGREE )
                Coefficient_Vector( Value_Location ) = 0.0_idp

            END DO
        END DO

        Value_Location =  Matrix_Location( ui, 0, 0, NUM_R_ELEMENTS-1, DEGREE  )
        Coefficient_Vector( Value_Location ) = 2.0_idp * sqrt(pi) * OUTER_CFA_BC_VALUES(ui)

    END IF

END DO


END SUBROUTINE CFA_3D_Dirichlet_BCs_Part1






!+104+###########################################################################!
!                                                                                !
!                  CFA_Dirichlet_BCs                                                  !
!                                                                                !
!################################################################################!
SUBROUTINE CFA_3D_Dirichlet_BCs_Part2( )



INTEGER                                             ::  i, l, m, ui, d

INTEGER                                             ::  Value_Location



!*!
!*!     Step 3
!*!
DO ui = 1,NUM_CFA_VARS



    !*!
    !*!     Innner BCs
    !*!

    IF ( INNER_CFA_BC_TYPE(ui) == 'D' .AND. myID_PETSc == 0 ) THEN

        DO l = 0,L_LIMIT
            DO m = -M_VALUES(l),M_VALUES(l)

                Value_Location =  MAtrix_Location( ui, l, m, 0, 0 )
                Block_RHS_VECTOR(Value_Location ) = 0.0_idp

                Value_Location =  MAtrix_Location( ui, l, m, 0, 0 )
                DO i = 0,ELEM_PROB_DIM-1
                    ! Clear the Column !
                    BLOCK_ELEM_STF_MATVEC( i*ELEM_PROB_DIM+value_location, 0)=0.0_idp

                END DO

                ! Clear the Row !
                BLOCK_ELEM_STF_MATVEC(Value_Location*ELEM_PROB_DIM:(Value_Location+1)*ELEM_PROB_DIM-1, 0) = 0.0_idp
                BLOCK_ELEM_STF_MATVEC(Value_Location*ELEM_PROB_DIM+Value_Location, 0) = 1.0_idp


            END DO
        END DO
            !*!
            !*!     Modify the Stiffness Matrix !
            !*!




    END IF






    !*!
    !*!     Outer BCs
    !*!

    IF ( OUTER_CFA_BC_TYPE(ui) == 'D' .AND. myID_PETSC == NUM_SUBSHELLS-1 ) THEN

        DO l = 0,L_LIMIT
            DO m = -M_VALUES(l),M_VALUES(l)

                Value_Location =  Matrix_Location( ui, l, m, NUM_R_ELEMS_PER_BLOCK-1, DEGREE )
                Block_RHS_VECTOR(Value_Location ) = 0.0_idp

                    !*!
                    !*!     Modify the Stiffness Matrix !
                    !*!
                Value_Location =  Matrix_Location( ui, l, m, 0, DEGREE )

                DO i = 0,ELEM_PROB_DIM-1

                    ! Clear the Column !
                    BLOCK_ELEM_STF_MATVEC( i*ELEM_PROB_DIM + Value_Location, NUM_R_ELEMS_PER_BLOCK-1) = 0.0_idp

                END DO
                ! Clear the Row !
                BLOCK_ELEM_STF_MATVEC( Value_Location*ELEM_PROB_DIM:(Value_Location+1)*ELEM_PROB_DIM-1,       &
                                     NUM_R_ELEMS_PER_BLOCK-1) = 0.0_idp

                BLOCK_ELEM_STF_MATVEC( Value_Location*ELEM_PROB_DIM+Value_Location, NUM_R_ELEMS_PER_BLOCK-1) = 1.0_idp

             END DO
         END DO

    END IF

END DO ! ui Loop


END SUBROUTINE CFA_3D_Dirichlet_BCs_Part2







!+105+###########################################################################!
!                                                                                !
!                  CFA_Neumann_BCs                                               !
!                                                                                !
!################################################################################!
SUBROUTINE CFA_3D_Neumann_BCs( )


! Nothing needs to be done. Score!


END SUBROUTINE CFA_3D_Neumann_BCs




!+401+###########################################################################!
!                                                                                !
!              Calc_Shift_BC_1D                                                  !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Shift_BC_1D( Shift_Vector_BC,                                    &
                          NUM_R_QUAD, NUM_T_QUAD, NUM_P_QUAD,                 &
                          NUM_R_ELEM, NUM_T_ELEM, NUM_P_ELEM, Dimn,           &
                          Sr, r_locs                                          )


REAL(KIND = idp),      INTENT( OUT )                               ::  Shift_Vector_BC

INTEGER,               INTENT( IN )                                ::  NUM_R_QUAD,  &
                                                                       NUM_T_QUAD,  &
                                                                       NUM_P_QUAD,  &
                                                                       NUM_R_ELEM,  &
                                                                       NUM_T_ELEM,  &
                                                                       NUM_P_ELEM,  &
                                                                       Dimn

REAL(KIND = idp), DIMENSION( 1:NUM_R_QUAD*NUM_T_QUAD*NUM_P_QUAD,  &
                             0:NUM_R_ELEM-1,                      &
                             0:NUM_T_ELEM-1,                      &
                             0:NUM_P_ELEM-1,                      &
                             1:Dimn ), INTENT(IN)                  ::  Sr

REAL(KIND = idp), DIMENSION( 0:NUM_R_ELEM ),  INTENT(IN)          ::  r_locs


COMPLEX(KIND = idp)                                               ::  Basis_Funcs,  &
                                                                      Tmp_Psi,      &
                                                                      Tmp_AlphaPsi

INTEGER                                                           ::  Ord,          &
                                                                      Current_Location


INTEGER                                                           ::  i, j, l, m, d, re, reb

REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                       :: x_locs,     &
                                                                     ri_locs,    &
                                                                     wi

REAL(KIND = idp), DIMENSION(:,:), ALLOCATABLE                     :: rij_locs,    &
                                                                     PSI_10,     &
                                                                     Sr_New

REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                       :: AlphaPsi
REAL(KIND = idp)                                                  :: Psi
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                       :: Beta_Tmp
REAL(KIND = idp)                                                  :: Inner_Int, Outer_Int

REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                       :: xlocP, yi, LagP
REAL(KIND = idp)                                                  :: x_tmp
Ord = 6




ALLOCATE( Beta_Tmp(1:NUM_R_ELEM + 1) )

ALLOCATE( LagP(0:DEGREE) )
ALLOCATE( xlocP(0:DEGREE) )
ALLOCATE( yi(0:DEGREE) )

ALLOCATE( x_locs(1:Ord) )
ALLOCATE( ri_locs(1:Ord) )
ALLOCATE( wi(1:Ord)  )
ALLOCATE( AlphaPsi(1:Ord) )

ALLOCATE( rij_locs(1:Ord,1:Ord) )
ALLOCATE( PSI_10(1:Ord,1:Ord)  )
ALLOCATE( Sr_New(1:Ord,1:Ord)  )

CALL Initialize_LG_Quadrature( Ord, x_locs, wi )
CALL Initialize_LGL_Quadrature( DEGREE, xlocP, yi)



Beta_tmp = 0.0_idp

DO re = 0,NUM_R_ELEM-1

   ! Calculate the r locations for the Outer Integral's Quadrature Points !
   ri_locs(:) = (r_locs(re+1)-r_locs(re))/2.0_idp * ( x_locs(:) + 1.0_idp ) + r_locs(re)

   ! Calculate the Alpha Psi values at each of the Outer Integral's Quadrature Points !
   DO i = 1,Ord
      AlphaPsi(i) =  1.0_idp + 0.5_idp*Potential_Solution(ri_locs(i),0.0_idp,0.0_idp)/C_Square
   END DO


   DO i = 1,Ord

      ! Calculate the Quadrature Points for each of the Inner Integrals
      rij_locs(:,i) = (ri_locs(i) - 0.0_idp)/2.0_idp *(x_locs(:) + 1.0_idp) + 0.0_idp


      ! Calculate Psi^10 values at each of the Inner Quadrature Points
      DO j = 1,Ord

           Psi = 1.0_idp - 0.5_idp*Potential_Solution(rij_locs(j,i),0.0_idp,0.0_idp)/C_Square
           Psi_10(j,i) = Psi**10
            
      END DO
      ! Calculate Sr values at each quadrature Point.
      DO j = 1,Ord
         DO reb = 0,NUM_R_ELEM - 1


            IF ( (rij_Locs(j,i) > r_locs(reb)) .AND. (rij_Locs(j,i) .LE. r_locs(reb+1)) ) THEN

!               Sr_New(j,i) = 1.0_idp/(r_locs(reb+1) - r_locs(reb))               &
!                           * ( Sr(1,reb,0,0,1)*(r_locs(reb+1) - rij_locs(j,i))   &
!                              +Sr(1,reb+1,0,0,1)*(rij_locs(j,i) - r_locs(reb))   )

                Sr_New(j,i) = Sr(1,reb,0,0,1)

               exit
            END IF
         END DO ! reb Loop

      END DO

   END DO ! i Loop


   ! Do the Outer Integral
   Outer_Int = 0.0_idp
   DO i = 1,Ord


     ! Do the Inner Integrals
     Inner_Int = 0.0_idp
     DO j = 1,Ord

        Inner_Int = Inner_Int                                 &
                  + rij_locs(j,i)*rij_locs(j,i)*rij_locs(j,i) &
                  * PSI_10(j,i)                               &
                  * Sr_New(j,i)                               &
                  * wi(j)

!                  PRINT*,"PSI_10",PSI_10(j,i),"Sr_New",Sr_New(j,i),"ri_locs(i)",ri_locs(i)

      END DO ! j Loop

!      PRINT*,"Inner Int",Inner_Int*( ri_locs(i) - 0.0_idp)/2.0_idp
      Outer_Int = Outer_Int                                       &
                + AlphaPsi(i)                                     &
                / (ri_locs(i)*ri_locs(i)*ri_locs(i)*ri_locs(i))   &    ! *** Changed
                * Inner_Int                                       &
                * ( ri_locs(i) - 0.0_idp)/2.0_idp                 &
                * wi(i)

   END DO ! i Loop





   IF ( re == 0 ) THEN

      Beta_Tmp(re+1) = (3.0_idp/2.0_idp)                         &
                               * 8.0_idp*pi* GR_Source_Scalar  &
                               * r_locs(re+1)                             &
                               * ( r_locs(re+1) - R_locs(re))/2.0_idp     &
                               * Outer_Int

   ELSE

      Beta_Tmp(Re+1) = r_locs(re+1)/r_locs(re)                   &
                               * Beta_Tmp(re)                    &
                               + (3.0_idp/2.0_idp)                        &
                               * 8.0_idp*pi* GR_Source_Scalar  &
                               * r_locs(re+1)                             &
                               * ( r_locs(re+1) - r_locs(re))/2.0_idp     &
                               * Outer_Int

   END IF

END DO ! re loop

Shift_Vector_BC = Beta_Tmp(NUM_R_ELEM)




END SUBROUTINE Calc_Shift_BC_1D


END MODULE Poseidon_BC_Module
