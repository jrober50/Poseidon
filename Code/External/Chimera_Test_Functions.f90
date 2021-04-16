   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE CHIMERA_TEST_FUNCS_Module                                                    !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
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
                        ONLY :  pi, eps

USE Units_Module, &
                        ONLY :  Grav_Constant_G,        &
                                C_Square

USE DRIVER_Parameters, &
                        ONLY :  DRIVER_R_ELEMS,        &
                                DRIVER_T_ELEMS,        &
                                DRIVER_P_ELEMS,        &
                                DRIVER_R_LOCS,         &
                                DRIVER_T_LOCS,         &
                                DRIVER_P_LOCS,         &
                                DRIVER_Delta_R,        &
                                DRIVER_Delta_T,        &
                                DRIVER_Delta_P,        &
                                Enclosed_Mass,          &
                                DRIVER_Potential,      &
                                DRIVER_SHIFT_VAL,      &
                                DRIVER_E,              &
                                DRIVER_S,              &
                                DRIVER_Si,             &
                                DRIVER_T_LOCS_LOCAL,   &
                                DRIVER_P_LOCS_LOCAL,   &
                                myID,                   &
                                myID_theta,             &
                                myID_phi,               &
                                ij_ray_dim,             &
                                ik_ray_dim,             &
                                MPI_COMM_GRID,          &
                                DRIVER_PROCS,          &
                                DRIVER_y_PROCS,        &
                                DRIVER_z_PROCS


IMPLICIT NONE


CONTAINS



!+101+###########################################################################!
!                                                                                !
!               INIT_CHIMERA_POTENTIAL                                           !
!                                                                                !
!################################################################################!
SUBROUTINE INIT_CHIMERA_POTENTIAL(  k_l, k_u, j_l, j_u, i_l, i_u, Density )

INTEGER, INTENT(IN)                                                 ::  k_l, k_u, &
                                                                        j_l, j_u, &
                                                                        i_l, i_u


REAL(KIND = idp), DIMENSION(k_l:k_u,j_l:j_u,i_l:i_u), INTENT(IN)    ::  Density

INTEGER                                                             ::  i

REAL(KIND = idp), DIMENSION(0:DRIVER_R_ELEMS)                       ::  CHIMERA_Potential


Enclosed_Mass(0) = (4.0_idp/3.0_idp)*pi*DENSITY(k_l,j_l,1)     &
                                       *DRIVER_R_LOCS(0)**3


DO i = 1,DRIVER_R_ELEMS

    Enclosed_Mass(i) = Enclosed_Mass(i-1)                                 &
                       + (4.0_idp/3.0_idp)*pi*DENSITY(k_l,j_l,i)  &
                                          *( DRIVER_R_LOCS(i)**3       &
                                           - DRIVER_R_LOCS(i-1)**3)

END DO



CHIMERA_Potential(DRIVER_R_ELEMS) = -Grav_Constant_G * Enclosed_Mass(DRIVER_R_ELEMS)       &
                                                      /DRIVER_R_LOCS(DRIVER_R_ELEMS)
DO i = DRIVER_R_ELEMS-1,2, -1

     CHIMERA_Potential(i) = CHIMERA_Potential(i+1) - Grav_Constant_G                         &
                                                   * Enclosed_Mass(i)                        &
                                                   / (DRIVER_R_LOCS(i)*DRIVER_R_LOCS(i))   &
                                                   * DRIVER_Delta_R(i)

END DO
CHIMERA_Potential(1) = CHIMERA_Potential(2) - 3*Grav_Constant_G*Enclosed_Mass(2)     &
                                              /(2*DRIVER_R_LOCS(2))

CHIMERA_Potential(0) = CHIMERA_Potential(1)



Driver_Potential = CHIMERA_Potential

END SUBROUTINE INIT_CHIMERA_POTENTIAL












!+101+###########################################################################!
!                                                                                !
!               INIT_CHIMERA_SHIFT_VAL                                           !
!                                                                                !
!################################################################################!
SUBROUTINE INIT_CHIMERA_SHIFT_VAL(  Num_Nodes,                                  &
                                    NUM_R_ELEM, NUM_T_ELEM, NUM_P_ELEM,         &
                                    Sr, r_locs                                  )


INTEGER, DIMENSION(1:3), INTENT( IN )                              ::  Num_Nodes
INTEGER,                 INTENT( IN )                              ::  NUM_R_ELEM,  &
                                                                       NUM_T_ELEM,  &
                                                                       NUM_P_ELEM

REAL(KIND = idp), DIMENSION( 1:NUM_Nodes(1)*Num_Nodes(2)*Num_Nodes(3),  &
                             0:NUM_R_ELEM-1,                            &
                             0:NUM_T_ELEM-1,                            &
                             0:NUM_P_ELEM-1,                            &
                             1:3 ), INTENT(IN)                    ::  Sr

REAL(KIND = idp), DIMENSION( 0:NUM_R_ELEM ),  INTENT(IN)          ::  r_locs


COMPLEX(KIND = idp)                                               ::  Basis_Funcs,  &
                                                                      Tmp_Psi,      &
                                                                      Tmp_AlphaPsi

INTEGER                                                           ::  Ord,          &
                                                                      Current_Location


INTEGER                                                           ::  i, j, l, m, d, re, reb


REAL(KIND = idp)                                                  :: csqr

REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                       :: x_locs,     &
                                                                     ri_locs,    &
                                                                     wi

REAL(KIND = idp), DIMENSION(:,:), ALLOCATABLE                     :: rij_locs,    &
                                                                     PSI_10,     &
                                                                     Sr_New

REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                       :: AlphaPsi
REAL(KIND = idp)                                                  :: Psi
REAL(KIND = idp)                                                  :: Outer_Int
REAL(KIND = idp)                                                  :: Inner_Int

REAL(KIND = idp)                                                  :: x_tmp
Ord = 6


csqr = C_SQUARE

ALLOCATE( x_locs(1:Ord) )
ALLOCATE( ri_locs(1:Ord) )
ALLOCATE( wi(1:Ord)  )
ALLOCATE( AlphaPsi(1:Ord) )

ALLOCATE( rij_locs(1:Ord,1:Ord) )
ALLOCATE( PSI_10(1:Ord,1:Ord)  )
ALLOCATE( Sr_New(1:Ord,1:Ord)  )

CALL Initialize_LG_Quadrature( Ord, x_locs, wi )


DRIVER_SHIFT_VAL(0) = 0.0_idp
DO re = 0,NUM_R_ELEM-1

   ! Calculate the r locations for the Outer Integral's Quadrature Points !
   ri_locs(:) = (r_locs(re+1)-r_locs(re))/2.0_idp * ( x_locs(:) + 1.0_idp ) + r_locs(re)

   ! Calculate the Alpha Psi values at each of the Outer Integral's Quadrature Points !
   DO i = 1,Ord
      AlphaPsi(i) =  1.0_idp + 0.5_idp*Test_Chimera_Simulated_Potential(ri_locs(i),0.0_idp,0.0_idp)/csqr
   END DO


   DO i = 1,Ord

      ! Calculate the Quadrature Points for each of the Inner Integrals
      rij_locs(:,i) = (ri_locs(i) - 0.0_idp)/2.0_idp *(x_locs(:) + 1.0_idp) + 0.0_idp


      ! Calculate Psi^10 values at each of the Inner Quadrature Points
      DO j = 1,Ord

           Psi = 1.0_idp - 0.5_idp*Test_Chimera_Simulated_Potential(rij_locs(j,i),0.0_idp,0.0_idp)/csqr
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

      DRIVER_SHIFT_VAL(re+1) = (3.0_idp/2.0_idp)                         &
                               * 8.0_idp*pi* Grav_Constant_G/(csqr*csqr)  &
                               * r_locs(re+1)                             &
                               * ( r_locs(re+1) - R_locs(re))/2.0_idp     &
                               * Outer_Int

   ELSE

      DRIVER_SHIFT_VAL(re+1) = r_locs(re+1)/r_locs(re)                   &
                               * DRIVER_SHIFT_VAL(re)                    &
                               + (3.0_idp/2.0_idp)                        &
                               * 8.0_idp*pi* Grav_Constant_G/(csqr*csqr)  &
                               * r_locs(re+1)                             &
                               * ( r_locs(re+1) - r_locs(re))/2.0_idp     &
                               * Outer_Int

   END IF

!    PRINT*,"Driver_Shift_Val",Driver_Shift_Val(re+1)
END DO ! re loop


END SUBROUTINE INIT_CHIMERA_SHIFT_VAL







!+301+##################################################################
!
!   Test_Chimera_Simulated_Potential
!
!#######################################################################
PURE FUNCTION Test_Chimera_Simulated_Potential(r, theta, phi)

REAL(KIND = idp),INTENT(IN)         :: r, theta, phi
REAL(KIND = idp)                    :: Test_Chimera_Simulated_Potential
REAL(KIND = idp)                    :: Sol
INTEGER                             :: i


IF ( r >= DRIVER_R_LOCS(DRIVER_R_ELEMS) ) THEN

     Sol = DRIVER_POTENTIAL(DRIVER_R_ELEMS)

ELSE

   DO i = 0,DRIVER_R_ELEMS-1

      if (( r >= DRIVER_R_LOCS(i) ) .AND. ( r < DRIVER_R_LOCS(i+1))) THEN

!          Sol = ( (r-DRIVER_R_LOCS(i))*DRIVER_POTENTIAL(i+1)        &
!                 +(DRIVER_R_LOCS(i+1)-r)*DRIVER_POTENTIAL(i) )      &
!                / (DRIVER_R_LOCS(i+1) - DRIVER_R_LOCS(i))
           Sol = DRIVER_POTENTIAL(i+1)

          EXIT
      END IF
   END DO
END IF

TEST_Chimera_Simulated_Potential = Sol


END FUNCTION Test_Chimera_Simulated_Potential





!+201+###########################################################################!
!                                                                                !
!                         SELFSIM_SHIFT_SOL                                       !
!                                                                                !
!################################################################################!
FUNCTION Test_Chimera_Simulated_Shift( r, r_locs,  NUM_R_ELEM )

REAL(KIND = idp), INTENT(IN)                              :: r
INTEGER, INTENT(IN)                                       :: NUM_R_ELEM

REAL(KIND = idp), DIMENSION(0:NUM_R_ELEM), INTENT(IN)     :: r_locs


REAL(KIND = idp)                 :: Test_Chimera_Simulated_Shift

INTEGER                          :: cur_entry
INTEGER                          :: i

REAL(KIND = idp)                 :: deltar

cur_entry = 0
DO i = 0,NUM_R_ELEM-1

   IF ( r == 0 ) THEN

       cur_entry = 0


   ELSE IF (( r > r_locs(i) ) .AND. ( r <= r_locs(i+1) ) ) THEN

       cur_entry = i

   END IF

END DO


deltar = r_locs(cur_entry+1) - r_locs(cur_entry)

Test_Chimera_Simulated_Shift = (1.0_idp/deltar)                                     &
                 *( DRIVER_SHIFT_VAL(cur_entry)*(r_locs(cur_entry+1) - r)         &
                   +DRIVER_SHIFT_VAL(cur_entry+1)*(r - r_locs(cur_entry))         )



END FUNCTION Test_Chimera_Simulated_Shift








!+503+##################################################################!
!                                                                       !
!   Initialize_LG_Quadrature - Calculate the Legendre-Gauss Quadrature  !
!                              node locations and weights.              !
!                                                                       !
!#######################################################################!
SUBROUTINE Initialize_LG_Quadrature(Ord, xloc, weights)

INTEGER, INTENT(IN)                                     ::  Ord
REAL(KIND = idp), INTENT(INOUT), DIMENSION(1:Ord)       ::  xloc, weights

INTEGER                                                 :: i, j, m
REAL(KIND = idp)                                        :: p1, p2, p3, pp, z, z1


m = (Ord + 1)/2

DO i = 1,m

    z = cos(pi * (i-0.25_idp)/(Ord + 0.5_idp))
    z1 = 42.0_idp

    DO WHILE ( ABS(z - z1) .GT. eps)
        p1 = 1.0_idp
        p2 = 0.0_idp

        DO j = 1, Ord

            p3 = p2
            p2 = p1
            p1 = ((2.0_idp*j - 1.0_idp)*z*p2 - (j - 1.0_idp)*p3)/j

        END DO


        pp = Ord*(z*p1 - p2)/(z*z-1.0_idp)
        z1 = z
        z = z1 - p1/pp


    END DO

    xloc(i) = -z
    xloc(Ord-i+1) = +z

    weights(i) = 2.0_idp/((1.0_idp - z*z)*pp*pp)
    weights(Ord-i+1) = weights(i)

END DO

END SUBROUTINE Initialize_LG_Quadrature




END MODULE CHIMERA_TEST_FUNCS_Module
