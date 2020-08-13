   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_Guess_Module                                                        !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   Initialize_Guess_Values                                             !##!
!##!    +102+   Initialize_Special_Guess_Values                                     !##!
!##!                                                                                !##!
!##!    +201+   Calc_Shift_1D                                                       !##!
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
                        ONLY :  idp, pi,                 &
                                TwoThirds,               &
                                FourThirds

USE Units_Module, &
                        ONLY :  Speed_of_Light,          &
                                C_Square


USE Poseidon_Parameters, &
                        ONLY :  DOMAIN_DIM,             &
                                DEGREE,                 &
                                L_LIMIT,                &
                                NUM_BLOCK_THETA_ROWS,   &
                                NUM_T_ELEMS_PER_BLOCK,  &
                                NUM_P_ELEMS_PER_BLOCK,  &
                                NUM_R_QUAD_POINTS,      &
                                NUM_T_QUAD_POINTS,      &
                                NUM_P_QUAD_POINTS

USE DRIVER_Parameters, &
                        ONLY :  myID,                   &
                                Potential_Solution,      &
                                Shift_Solution


USE Poseidon_Variables_Module, &
                        ONLY :  NUM_R_ELEMENTS,         &
                                NUM_T_ELEMENTS,         &
                                NUM_P_ELEMENTS,         &
                                NUM_TP_QUAD_POINTS,       &
                                rlocs, tlocs, plocs,    &
                                NUM_R_NODES,            &
                                INT_R_LOCATIONS,        &
                                INT_T_LOCATIONS,        &
                                INT_P_LOCATIONS,        &
                                VAR_DIM,                &
                                PROB_DIM,               &
                                ULM_LENGTH,             &
                                Ylm_Table_Block,        &
                                Ylm_Values,             &
                                Ylm_dt_Values,          &
                                Ylm_dp_Values,          &
                                Ylm_CC_Values,          &
                                Ylm_CC_dt_Values,       &
                                Ylm_CC_dp_Values,       &
                                LM_LENGTH,              &
                                myID_Shell,             &
                                Lagrange_Poly_Table,    &
                                LPT_LPT,                &
                                Coefficient_Vector,     &
                                LM_Location,            &
                                Matrix_Location,        &
                                M_VALUES,               &
                                R_Inner,                &
                                R_Outer


USE Poseidon_Mapping_Functions_Module, &
        ONLY :  Map_To_X_Space,         &
                Map_From_X_Space

USE Poseidon_Quadrature_Module, &
        ONLY :  Initialize_LG_Quadrature,               &
                Initialize_LGL_Quadrature,              &
                Initialize_LGL_Quadrature_Locations


USE Poseidon_Math_Functions_Module, &
        ONLY :  Lagrange_Poly,      &
                Spherical_Harmonic

IMPLICIT NONE


CONTAINS



!+101+###########################################################################!
!                                                                                !
!                  Initialize_Guess_Values          !
!                                                                                !
!################################################################################!
SUBROUTINE Initialize_Guess_Values()


INTEGER                                     ::  beta_i, j

REAL(KIND = idp)                            ::  Beta_Start,     &
                                                delta_Beta



INTEGER                                     ::  CUR_PSI_LOC,        &
                                                CUR_ALPHPSI_LOC,    &
                                                CUR_BETA_LOC


INTEGER                                     :: re, rd



!
!   Empty Space Initial Guess
!
Coefficient_Vector(:) = 0.0_idp


Beta_Start = 0.0_idp
!Beta_Start = 1.0E-12
!Beta_Start = 0.05000_idp*2.0_idp*sqrt(pi)


delta_Beta = 0.0_idp
!delta_Beta = 0.0001_idp
!delta_Beta = Beta_Start/VAR_DIM







DO re = 0,NUM_R_ELEMENTS - 1


    DO rd = 0,DEGREE


        ! 2 sqrt(pi) is Ylm normalization factor


        CUR_PSI_LOC = Matrix_Location( 1, 0, 0, re, rd )
        Coefficient_Vector(CUR_PSI_LOC) = 0.9_idp * 2.0_idp * sqrt(pi)




        CUR_ALPHPSI_LOC = Matrix_Location( 2, 0, 0, re, rd )
        Coefficient_Vector(CUR_ALPHPSI_LOC) = 0.9_idp * 2.0_idp * sqrt(pi)


        DO beta_i = 1,DOMAIN_DIM-1

            CUR_BETA_LOC = Matrix_Location( 2+beta_i, 0, 0, re, rd )
            Coefficient_Vector(CUR_BETA_LOC) = Beta_Start - delta_Beta*(re*DEGREE+rd)

        END DO


    END DO


END DO




END SUBROUTINE Initialize_Guess_Values








!+102+###########################################################################!
!                                                                                !
!                  Initialize_Guess_Values          !
!                                                                                !
!################################################################################!
SUBROUTINE Initialize_Special_Guess_Values()


INTEGER                                 :: re, d, Here


REAL(KIND = idp), DIMENSION(0:DEGREE)                           :: R_Values


REAL(KIND = idp), DIMENSION(0:DEGREE)                           ::  Local_Locations


INTEGER                                                         ::  CUR_PSI_LOC,    &
                                                                    CUR_ALPHPSI_LOC,&
                                                                    CUR_SHIFT_LOC

INTEGER          :: Matrix_Map


Local_Locations = Initialize_LGL_Quadrature_Locations(DEGREE)


!
!   Empty Space Initial Guess
!
Coefficient_Vector = 0.0_idp




!PRINT*,"ALPHA_PSI initial guess set to one"


DO re = 0,NUM_R_ELEMENTS - 1

    R_Values = Map_From_X_Space(rlocs(re), rlocs(re + 1), Local_Locations)

    DO d = 0,DEGREE
 
        ! 2 sqrt(pi) is Ylm normalization factor


        CUR_PSI_LOC = Matrix_Location( 1, 0, 0, re, d )

        Coefficient_Vector(CUR_PSI_LOC) = 2.0_idp * sqrt(pi)                                        &
                                        * ( 1.0_idp - 0.5_idp                                       &
                                            * Potential_Solution(R_Values(d),0.0_idp,0.0_idp)/C_Square  )



        CUR_ALPHPSI_LOC = Matrix_Location( 2, 0, 0, re, d )
!        Coefficient_Vector(CUR_ALPHPSI_LOC) = 2.0_idp * sqrt(pi)
        Coefficient_Vector(CUR_ALPHPSI_LOC) = 2.0_idp * sqrt(pi)                                    &
                                        * ( 1.0_idp + 0.5_idp                                       &
                                            * Potential_Solution(R_Values(d),0.0_idp,0.0_idp)/C_Square  )



        CUR_SHIFT_LOC = Matrix_Location( 3, 0, 0, re, d )
        Coefficient_Vector(Cur_Shift_Loc) = 2.0_idp*sqrt(pi)*Shift_Solution(R_Values(d),rlocs,NUM_R_ELEMENTS)

    END DO
END DO



END SUBROUTINE Initialize_Special_Guess_Values





!+201+###########################################################################!
!                                                                                !
!              Calc_Shift_1D                                                  !
!                                                                                !
!################################################################################!
SUBROUTINE Calc_Shift_1D( Shift_Vector,                                       &
                             NUM_R_QUAD, NUM_T_QUAD, NUM_P_QUAD,                 &
                             NUM_R_ELEM, NUM_T_ELEM, NUM_P_ELEM, Dimn,           &
                             Sr, r_locs                                          )


REAL(KIND = idp), DIMENSION( 0:NUM_R_ELEM ),INTENT( OUT )          ::  Shift_Vector

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
REAL(KIND = idp)                                                  :: Beta_Tmp
REAL(KIND = idp)                                                  :: Inner_Int

REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                       :: xlocP, yi, LagP
REAL(KIND = idp)                                                  :: x_tmp
Ord = 6






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

Shift_Vector(0) = 0.0_idp
DO re = 0,NUM_R_ELEM-1

   ! Calculate the r locations for the Outer Integral's Quadrature Points !
   ri_locs(:) = (r_locs(re+1)-r_locs(re))/2.0_idp * ( x_locs(:) + 1.0_idp ) + r_locs(re)




   ! Calculate the Alpha Psi values at each of the Outer Integral's Quadrature Points !
   DO i = 1,Ord

      x_tmp = Map_To_X_Space(r_locs(re), r_locs(re+1), ri_locs(i))
      LagP = Lagrange_Poly(x_tmp, DEGREE, xlocP)
      Tmp_AlphaPsi = 0.0_idp

      DO l = 0,L_LIMIT
         DO m = -l,l
            DO d = 0,DEGREE

               Current_Location = Matrix_Location( 2, l, m, re, d )
               Basis_Funcs = Spherical_Harmonic(l,m,0.0_idp,0.0_idp)*LagP(d)

               Tmp_AlphaPsi = Tmp_AlphaPsi                                          &
                            + Coefficient_Vector(Current_Location)                  &
                            * Basis_Funcs


             END DO ! d Loop
          END DO ! m Loop
       END DO ! l Loop

    

       AlphaPsi(i) = REAL( Tmp_AlphaPsi, KIND = idp )

   END DO ! i Loop





   DO i = 1,Ord

      ! Calculate the Quadrature Points for each of the Inner Integrals
      rij_locs(:,i) = (ri_locs(i) - R_Inner)/2.0_idp *(x_locs(:) + 1.0_idp) + R_Inner


      ! Calculate Psi^10 values at each of the Inner Quadrature Points
      DO j = 1,Ord
         DO reb = 0,NUM_R_ELEM - 1
            IF ( (rij_Locs(j,i) > r_locs(reb)) .AND. (rij_Locs(j,i) < rlocs(reb+1)) ) THEN

               x_tmp = Map_To_X_Space(r_locs(reb), r_locs(reb+1), rij_locs(j,i))
               LagP = Lagrange_Poly(x_tmp, DEGREE, xlocP)
               Tmp_Psi = 0.0_idp

               DO l = 0,L_LIMIT
                  DO m = -l,l
                     DO d = 0,DEGREE

                        Current_Location = Matrix_Location( 1, l, m, reb, d )
                        Basis_Funcs = Spherical_Harmonic(l,m,0.0_idp,0.0_idp)*LagP(d)

                        Tmp_Psi = Tmp_Psi                                                    &
                                + Coefficient_Vector(Current_Location+0)                     &
                                * Basis_Funcs


                     END DO ! d Loop
                  END DO ! m Loop
               END DO ! l Loop

               ! Interpolate Sr Value to Quadrature Point !
               Sr_New(j,i) = Sr(1,reb,0,0,1)

            END IF
         END DO ! reb Loop

         Psi = REAL( Tmp_Psi, KIND = idp)
         Psi_10(j,i) = Psi**10

      END DO

   END DO ! i Loop



   ! Do the Outer Integral
   Beta_Tmp = 0.0_idp
   DO i = 1,Ord


     ! Do the Inner Integrals
     Inner_Int = 0.0_idp
     DO j = 1,Ord

        Inner_Int = Inner_Int                                 &
                  + rij_locs(j,i)*rij_locs(j,i)*rij_locs(j,i) &
                  * PSI_10(j,i)                               &
                  * Sr_New(j,i)                               &
                  * wi(j)


!         WRITE(*,'(ES22.15,1X,ES22.15,1X,ES22.5)')rij_locs(j,i)*rij_locs(j,i)*rij_locs(j,i),PSI_10(j,i),Sr_New(j,i)
      END DO ! j Loop

!      WRITE(*,'(I3,1X,I3,1X,ES22.15)') re,i,Inner_Int*(ri_locs(i) - R_Inner)/2.0_idp


      Beta_Tmp = Beta_Tmp                                        &
               + AlphaPsi(i)                                     &
               / (ri_locs(i)*ri_locs(i))                         &
               * Inner_Int                                       &
               * ( ri_locs(i) - R_Inner)/2.0_idp                 &
               * wi(i)

   

   END DO ! i Loop






   IF ( re == 0 ) THEN

      Shift_Vector(re+1) = (3.0_idp/2.0_idp)                        &
                         * r_locs(re+1)                             &
                         * ( r_locs(re+1) - R_locs(re))/2.0_idp        &
                         * Beta_Tmp

   ELSE

      Shift_Vector(re+1) = r_locs(re+1)/r_locs(re)                  &
                         * Shift_Vector(re)                         &
                         + (3.0_idp/2.0_idp)                        &
                         * r_locs(re+1)                             &
                         * ( r_locs(re+1) - r_locs(re))/2.0_idp        &
                         * Beta_Tmp

   END IF



END DO ! re loop







END SUBROUTINE Calc_Shift_1D





END MODULE Poseidon_Guess_Module
