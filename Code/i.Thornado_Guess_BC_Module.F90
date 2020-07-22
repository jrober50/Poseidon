   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Thornado_Guess_BC_Module                                                     !##!
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

USE Poseidon_Constants_Module, &
                            ONLY :  idp,        &
                                    pi,         &
                                    eps,        &
                                    Speed_of_Light, &
                                    Grav_Constant_G

USE Poseidon_Variables_Module,  &
                            ONLY:   NEWTON_SHIFT_VEC,   &
                                    NEWTON_POT_VEC,     &
                                    Block_Source_Si,    &
                                    rlocs,              &
                                    NUM_R_ELEMENTS,     &
                                    Matrix_Location,    &
                                    Coefficient_Vector


USE Poseidon_Parameters, &
                            ONLY :  NUM_R_QUAD_POINTS,  &
                                    DEGREE,             &
                                    DOMAIN_DIM

IMPLICIT NONE

LOGICAL                                 :: Potential_Set_Flag = .FALSE.
LOGICAL                                 :: Shift_Set_Flag = .FALSE.



CONTAINS


!+101+###########################################################################!
!                                                                                !
!                  CALC_LAPSE_BC_VALUE                                           !
!                                                                                !
!################################################################################!
SUBROUTINE CALC_BC_VALUES(Num_R_Input_Nodes, Left_Limit, Right_Limit,       &
                            Node_Locations, r_locs, NUM_R_ELEM,             &
                            Density, Momentum_Density,                      &
                            Lapse_BC, ConFact_BC, Shift_BC                  )

INTEGER, INTENT(IN)                                             :: Num_R_Input_Nodes
REAL(KIND = idp), INTENT(IN)                                    :: Left_Limit
REAL(KIND = idp), INTENT(IN)                                    :: Right_Limit
REAL(KIND = idp), DIMENSION(1:Num_R_Input_Nodes), INTENT(IN)    :: Node_Locations
REAL(KIND = idp), DIMENSION(0:NUM_R_ELEM), INTENT(IN)           :: r_locs
INTEGER, INTENT(IN)                                             :: NUM_R_ELEM

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:Num_R_Input_Nodes,   &
                                         1:Num_R_Elem,           &
                                         1:1,                   &
                                         1:1  )                 :: Density

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:Num_R_Input_Nodes,   &
                                         1:Num_R_Elem,          &
                                         1:1,                   &
                                         1:1  )                 :: Momentum_Density

REAL(KIND = idp), INTENT(INOUT)                                 ::  Lapse_BC
REAL(KIND = idp), INTENT(INOUT)                                 ::  ConFact_BC
REAL(KIND = idp), INTENT(INOUT)                                 ::  Shift_BC



CALL CALC_LAPSE_BC_VALUE( Num_R_Input_Nodes, Left_Limit, Right_Limit,   &
                            Node_Locations, r_locs, NUM_R_ELEM,           &
                            Density, Lapse_BC                           )

CALL CALC_CONFACT_BC_VALUE( Num_R_Input_Nodes, Left_Limit, Right_Limit,   &
                                  Node_Locations, r_locs, NUM_R_ELEM,           &
                                  Density, Confact_BC                           )

CALL CALC_SHIFT_BC_VALUE( Num_R_Input_Nodes, Left_Limit, Right_Limit,       &
                                Node_Locations, r_locs, NUM_R_ELEM,               &
                                Density, Momentum_Density, Shift_BC               )




END SUBROUTINE CALC_BC_VALUES





!+101+###########################################################################!
!                                                                                !
!                  CALC_LAPSE_BC_VALUE                                           !
!                                                                                !
!################################################################################!
SUBROUTINE CALC_LAPSE_BC_VALUE( Num_R_Input_Nodes, Left_Limit, Right_Limit,   &
                                Node_Locations, r_locs, NUM_R_ELEM,           &
                                Density, Lapse_BC                           )

INTEGER, INTENT(IN)                                             :: Num_R_Input_Nodes
REAL(KIND = idp), INTENT(IN)                                    :: Left_Limit
REAL(KIND = idp), INTENT(IN)                                    :: Right_Limit
REAL(KIND = idp), DIMENSION(1:Num_R_Input_Nodes), INTENT(IN)    :: Node_Locations
REAL(KIND = idp), DIMENSION(0:NUM_R_ELEM), INTENT(IN)           :: r_locs
INTEGER, INTENT(IN)                                             :: NUM_R_ELEM

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:Num_R_Input_Nodes,   &
                                         1:Num_R_Elem,          &
                                         1:1,                   &
                                         1:1                    )   :: Density


REAL(KIND = idp), INTENT(INOUT)                                 :: Lapse_BC




IF ( Potential_Set_Flag .eqv. .FALSE.) THEN

    CALL Set_Potential( Num_R_Input_Nodes, Left_Limit, Right_Limit,       &
                        Node_Locations, r_locs, NUM_R_ELEM,               &
                        Density                                           )

END IF


Lapse_BC = 1.0_idp + 0.5_idp                                                          &
                    * NEWTON_POTENTIAL_1D_VAL(r_locs(Num_R_Elem), r_locs, NUM_R_ELEM)   &
                    / (Speed_of_Light*Speed_of_Light)

END SUBROUTINE CALC_LAPSE_BC_VALUE




!+102+###########################################################################!
!                                                                                !
!                  CALC_CONFACT_BC_VALUE                                           !
!                                                                                !
!################################################################################!
SUBROUTINE CALC_CONFACT_BC_VALUE( Num_R_Input_Nodes, Left_Limit, Right_Limit,   &
                                  Node_Locations, r_locs, NUM_R_ELEM,           &
                                  Density, Confact_BC                           )

INTEGER, INTENT(IN)                                             :: Num_R_Input_Nodes
REAL(KIND = idp), INTENT(IN)                                    :: Left_Limit
REAL(KIND = idp), INTENT(IN)                                    :: Right_Limit
REAL(KIND = idp), DIMENSION(1:Num_R_Input_Nodes), INTENT(IN)    :: Node_Locations
REAL(KIND = idp), DIMENSION(0:NUM_R_ELEM), INTENT(IN)           :: r_locs
INTEGER, INTENT(IN)                                             :: NUM_R_ELEM

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:Num_R_Input_Nodes,   &
                                         1:Num_R_Elem,          &
                                         1:1,                   &
                                         1:1                    )   :: Density


REAL(KIND = idp), INTENT(INOUT)                                 :: Confact_BC




IF ( Potential_Set_Flag .eqv. .FALSE.) THEN

    CALL Set_Potential( Num_R_Input_Nodes, Left_Limit, Right_Limit,       &
                        Node_Locations, r_locs, NUM_R_ELEM,               &
                        Density                                           )

END IF


Confact_BC = 1.0_idp - 0.5_idp                                                          &
                    * NEWTON_POTENTIAL_1D_VAL(r_locs(Num_R_Elem), r_locs, NUM_R_ELEM)   &
                    / (Speed_of_Light*Speed_of_Light)


END SUBROUTINE CALC_CONFACT_BC_VALUE




!+103+###########################################################################!
!                                                                                !
!                  CALC_SHIFT_BC_VALUE                                           !
!                                                                                !
!################################################################################!
SUBROUTINE CALC_SHIFT_BC_VALUE( Num_R_Input_Nodes, Left_Limit, Right_Limit,       &
                                Node_Locations, r_locs, NUM_R_ELEM,               &
                                Density, Momentum_Density, Shift_BC               )


INTEGER, INTENT(IN)                                             :: Num_R_Input_Nodes
REAL(KIND = idp), INTENT(IN)                                    :: Left_Limit
REAL(KIND = idp), INTENT(IN)                                    :: Right_Limit
REAL(KIND = idp), DIMENSION(1:Num_R_Input_Nodes), INTENT(IN)    :: Node_Locations
REAL(KIND = idp), DIMENSION(0:NUM_R_ELEM), INTENT(IN)           :: r_locs
INTEGER, INTENT(IN)                                             :: NUM_R_ELEM

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:Num_R_Input_Nodes,   &
                                         1:Num_R_Elem,           &
                                         1:1,                   &
                                         1:1  )                 :: Density

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:Num_R_Input_Nodes,   &
                                         1:Num_R_Elem,          &
                                         1:1,                   &
                                         1:1  )                 :: Momentum_Density

REAL(KIND = idp), INTENT(INOUT)                                 :: Shift_BC



IF( Shift_Set_Flag .eqv. .FALSE. ) THEN

    CALL Set_Shift( Num_R_Input_Nodes, Left_Limit, Right_Limit,       &
                    Node_Locations, r_locs, NUM_R_ELEM,               &
                    Density, Momentum_Density                         )
END IF


Shift_BC = NEWTON_SHIFT_1D_VAL( r_locs(Num_R_Elem), r_locs,  NUM_R_ELEM )





END SUBROUTINE CALC_SHIFT_BC_VALUE








!+201+###########################################################################!
!                                                                                !
!                  SET_POTENTIAL                                                 !
!                                                                                !
!################################################################################!
SUBROUTINE SET_POTENTIAL( Num_R_Input_Nodes, Left_Limit, Right_Limit,       &
                          Node_Locations, r_locs, NUM_R_ELEM,               &
                          Density                                           )


INTEGER, INTENT(IN)                                             :: Num_R_Input_Nodes
REAL(KIND = idp), INTENT(IN)                                    :: Left_Limit
REAL(KIND = idp), INTENT(IN)                                    :: Right_Limit
REAL(KIND = idp), DIMENSION(1:Num_R_Input_Nodes), INTENT(IN)    :: Node_Locations
REAL(KIND = idp), DIMENSION(0:NUM_R_ELEM), INTENT(IN)           :: r_locs
INTEGER, INTENT(IN)                                             :: NUM_R_ELEM

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:Num_R_Input_Nodes,   &
                                         1:Num_R_Elem,          &
                                         1:1,                   &
                                         1:1                    )   :: Density



REAL(KIND = idp), DIMENSION( 1:Num_R_Input_Nodes, 1:Num_R_Elem )        :: Enclosed_Mass
REAL(KIND = idp)                                                        :: deltax
REAL(KIND = idp), DIMENSION(1:Num_R_Input_Nodes)                        :: Quad_Locations
REAL(KIND = idp)                                                        :: delta_ratio
REAL(KIND = idp), DIMENSION( 1:Num_R_Input_Nodes, 1:Num_R_Elem )        :: cur_r_loc
REAL(KIND = idp)                                                        :: last_r_loc
REAL(KIND = idp)                                                        :: last_mass
REAL(KIND = idp)                                                        :: Shell_Width
INTEGER     :: i, re, q



deltax = Right_Limit - Left_Limit
Quad_Locations = Node_Locations - Left_Limit

Enclosed_Mass = 0.0_idp
Last_R_Loc = 0.0_idp
Last_Mass = 0.0_idp
DO re = 1,Num_R_Elem
    DO q = 1,Num_R_Input_Nodes
        delta_ratio = (r_locs(re+1) - r_locs(re))/deltax
        cur_r_loc(q,re) = r_locs(i) + Quad_Locations(q)*delta_ratio
        Shell_Width = cur_r_loc(q,re)**3 - Last_R_Loc**3

        Enclosed_Mass(q,re) = Last_Mass + (4.0_idp/3.0_idp)* pi*Density(q,re,1,1)*Shell_Width


        Last_Mass = Enclosed_Mass(q,re)
        Last_R_Loc = cur_r_loc(q,re)
    END DO
END DO



re = Num_R_Elem
q = Num_R_Input_Nodes
NEWTON_POT_VEC(q,re) = -GRAV_Constant_G*Enclosed_Mass(q,re)/cur_r_loc(q,re)
Last_R_Loc = cur_r_loc(q,re)
DO q = 1,Num_R_Input_Nodes-2,-1

    NEWTON_POT_VEC(q,re) = NEWTON_POT_VEC(q+1,re)                               &
                         - Grav_Constant_G*Enclosed_Mass(q,re)                        &
                         * ( cur_r_loc(q,re) - cur_r_loc(q+1,re) )                           &
                         / (cur_r_loc(q,re)*cur_r_loc(q,re) )


END DO

DO re = 1,Num_R_ELEM-1,-1
    q = Num_R_Input_Nodes
    NEWTON_POT_VEC(q,re) = NEWTON_POT_VEC(1,re+1)                               &
                         - Grav_Constant_G*Enclosed_Mass(q,re)                        &
                         * ( cur_r_loc(1,re+1) - cur_r_loc(q,re))                &
                         / (cur_r_loc(q,re)*cur_r_loc(q,re) )

    DO q = 1,Num_R_Input_Nodes-1,-1

        NEWTON_POT_VEC(q,re) = NEWTON_POT_VEC(q+1,re)                               &
                             - Grav_Constant_G*Enclosed_Mass(q,re)                        &
                             * ( cur_r_loc(q,re) - cur_r_loc(q+1,re) )                           &
                             / (cur_r_loc(q,re)*cur_r_loc(q,re) )


    END DO
END DO




END SUBROUTINE SET_POTENTIAL





!+201+###########################################################################!
!                                                                                !
!                  SET_SHIFT                                                     !
!                                                                                !
!################################################################################!
SUBROUTINE SET_SHIFT( Num_R_Input_Nodes, Left_Limit, Right_Limit,       &
                      Node_Locations, r_locs, NUM_R_ELEM,               &
                      Density, Momentum_Density                         )


INTEGER, INTENT(IN)                                             :: Num_R_Input_Nodes
REAL(KIND = idp), INTENT(IN)                                    :: Left_Limit
REAL(KIND = idp), INTENT(IN)                                    :: Right_Limit
REAL(KIND = idp), DIMENSION(1:Num_R_Input_Nodes), INTENT(IN)    :: Node_Locations
REAL(KIND = idp), DIMENSION(0:NUM_R_ELEM), INTENT(IN)           :: r_locs
INTEGER, INTENT(IN)                                             :: NUM_R_ELEM

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:Num_R_Input_Nodes,   &
                                         1:Num_R_Elem,      &
                                         1:1,                   &
                                         1:1  )                 :: Density

REAL(KIND = idp), INTENT(IN), DIMENSION( 1:Num_R_Input_Nodes,   &
                                         1:Num_R_Elem,      &
                                         1:1,                   &
                                         1:1  )                 :: Momentum_Density

REAL(KIND = idp)                                                        :: deltax
REAL(KIND = idp), DIMENSION(1:Num_R_Input_Nodes)                        :: Quad_Locations
REAL(KIND = idp)                                                        :: delta_ratio
REAL(KIND = idp), DIMENSION( 1:Num_R_Input_Nodes, 1:Num_R_Elem )    :: cur_r_loc
REAL(KIND = idp)                                                        :: last_R_loc

INTEGER                                                           ::  i, j, l, m, d, re, reb,q
INTEGER                                                           ::  Ord

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
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                       :: Inner_Int

REAL(KIND = idp)                                                  :: x_tmp
Ord = 6


csqr = Speed_of_Light*Speed_of_Light


ALLOCATE( x_locs(1:Ord) )
ALLOCATE( ri_locs(1:Ord) )
ALLOCATE( wi(1:Ord)  )
ALLOCATE( AlphaPsi(1:Ord) )

ALLOCATE( Inner_Int(1:Ord) )

ALLOCATE( rij_locs(1:Ord,1:Ord) )
ALLOCATE( PSI_10(1:Ord,1:Ord)  )
ALLOCATE( Sr_New(1:Ord,1:Ord)  )

CALL Initialize_LG_Quadrature( Ord, x_locs, wi )


IF ( Potential_Set_Flag .eqv. .FALSE.) THEN

    CALL Set_Potential( Num_R_Input_Nodes, Left_Limit, Right_Limit,       &
                        Node_Locations, r_locs, NUM_R_ELEM,               &
                        Density                                           )

END IF


deltax = Right_Limit - Left_Limit
Quad_Locations = Node_Locations - Left_Limit

Last_R_Loc = 0.0_idp
DO re = 1,Num_R_Elem
    DO q = 1,Num_R_Input_Nodes

        delta_ratio = (r_locs(re+1) - r_locs(re))/deltax
        cur_r_loc(q,re) = r_locs(i) + Quad_Locations(q)*delta_ratio

        Last_R_Loc = cur_r_loc(q,re)
    END DO
END DO



DO re = 1,Num_R_ELEM

    ! Calculate the r locations for the Outer Integral's Quadrature Points !
    ri_locs(:) = (r_locs(re+1)-r_locs(re))/2.0_idp * ( x_locs(:) + 1.0_idp ) + r_locs(re)


    DO i = 1,NUM_R_INPUT_NODES
        AlphaPsi(i) =  1.0_idp + 0.5_idp*NEWTON_POTENTIAL_1D_VAL( ri_locs(i), r_locs, NUM_R_ELEM )/csqr
    END DO

   DO i = 1,Ord

      ! Calculate the Quadrature Points for each of the Inner Integrals
      rij_locs(:,i) = (ri_locs(i) - 0.0_idp)/2.0_idp *(x_locs(:) + 1.0_idp) + 0.0_idp


      ! Calculate Psi^10 values at each of the Inner Quadrature Points
      DO j = 1,Ord

           Psi = 1.0_idp - 0.5_idp*NEWTON_POTENTIAL_1D_VAL(rij_locs(j,i),r_locs, NUM_R_ELEM)/csqr
           Psi_10(j,i) = Psi**10

      END DO

      ! Calculate Sr values at each quadrature Point.
      DO j = 1,Ord
    
            Sr_New(j,i) = 1.0_idp/(r_locs(reb+1) - r_locs(reb))                                     &
                           * ( Block_Source_Si(1,1,1,re-1,0,0,1)*(r_locs(reb+1) - rij_locs(j,i))    &
                              +Block_Source_Si(Num_R_Quad_Points,1,1,re-1,0,0,1)*(rij_locs(j,i) - r_locs(reb))   )



      END DO

   END DO ! i Loop




    ! Do the Inner Integrals
    Inner_Int = 0.0_idp
    DO i = 1,Ord
        DO j = 1,Ord

            Inner_Int(i) = Inner_Int(i)                           &
                      + rij_locs(j,i)*rij_locs(j,i)*rij_locs(j,i) &
                      * PSI_10(j,i)                               &
                      * Sr_New(j,i)                               &
                      * wi(j)

        !                  PRINT*,"PSI_10",PSI_10(j,i),"Sr_New",Sr_New(j,i),"ri_locs(i)",ri_locs(i)

        END DO ! j Loop
    END DO ! i Loop




   ! Do the Outer Integral
   Outer_Int = 0.0_idp
   DO i = 1,Ord
!      PRINT*,"Inner Int",Inner_Int*( ri_locs(i) - 0.0_idp)/2.0_idp

      Outer_Int = Outer_Int                                       &
                + AlphaPsi(i)                                     &
                / (ri_locs(i)*ri_locs(i)*ri_locs(i)*ri_locs(i))   &    ! *** Changed
                * Inner_Int(i)                                    &
                * ( ri_locs(i) - 0.0_idp)/2.0_idp                 &
                * wi(i)

   END DO ! i Loop



   IF ( re == 0 ) THEN

      NEWTON_SHIFT_VEC(re+1) = (3.0_idp/2.0_idp)                        &
                               * 8.0_idp*pi* Grav_Constant_G/(csqr*csqr)  &
                               * r_locs(re+1)                             &
                               * ( r_locs(re+1) - R_locs(re))/2.0_idp     &
                               * Outer_Int

   ELSE

      NEWTON_SHIFT_VEC(re+1) = r_locs(re+1)/r_locs(re)                     &
                               * NEWTON_SHIFT_VEC(re)                      &
                               + (3.0_idp/2.0_idp)                        &
                               * 8.0_idp*pi* Grav_Constant_G/(csqr*csqr)  &
                               * r_locs(re+1)                             &
                               * ( r_locs(re+1) - r_locs(re))/2.0_idp     &
                               * Outer_Int

   END IF




END DO ! re







END SUBROUTINE SET_SHIFT






!+301+###########################################################################!
!                                                                                !
!                         NEWTON_POTENTIAL_VAL                                       !
!                                                                                !
!################################################################################!
FUNCTION NEWTON_POTENTIAL_1D_VAL( r, r_locs, NUM_R_ELEM )

REAL(KIND = idp), INTENT(IN)                              :: r
REAL(KIND = idp), DIMENSION(0:NUM_R_ELEM), INTENT(IN)     :: r_locs
INTEGER, INTENT(IN)                                       :: NUM_R_ELEM

REAL(KIND = idp)                 :: NEWTON_POTENTIAL_1D_VAL

INTEGER                          :: cur_entry
INTEGER                          :: i

REAL(KIND = idp)                 :: deltar


DO i = 0,NUM_R_ELEM-1

   IF ( r == 0 ) THEN

       cur_entry = 0

   ELSE IF (( r > r_locs(i) ) .AND. ( r <= r_locs(i+1) ) ) THEN

       cur_entry = i

   END IF
 
END DO


deltar = r_locs(cur_entry+1) - r_locs(cur_entry)

NEWTON_POTENTIAL_1D_VAL = (1.0_idp/deltar)    &
                 *( NEWTON_POT_VEC(1,cur_entry)*(r_locs(cur_entry+1) - r)         &
                   +NEWTON_POT_VEC(NUM_R_QUAD_POINTS-1,cur_entry)*(r - r_locs(cur_entry))         )



END FUNCTION NEWTON_POTENTIAL_1D_VAL



!+302+###########################################################################!
!                                                                                !
!                         NEWTON_SHIFT_VAL                                       !
!                                                                                !
!################################################################################!
FUNCTION NEWTON_SHIFT_1D_VAL( r, r_locs,  NUM_R_ELEM )

REAL(KIND = idp), INTENT(IN)                              :: r
REAL(KIND = idp), DIMENSION(0:NUM_R_ELEM), INTENT(IN)     :: r_locs
INTEGER, INTENT(IN)                                       :: NUM_R_ELEM


REAL(KIND = idp)                 :: NEWTON_SHIFT_1D_VAL

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

NEWTON_SHIFT_1D_VAL = (1.0_idp/deltar)    &
                 *( NEWTON_SHIFT_VEC(cur_entry)*(r_locs(cur_entry+1) - r)         &
                   +NEWTON_SHIFT_VEC(cur_entry+1)*(r - r_locs(cur_entry))         )



END FUNCTION NEWTON_SHIFT_1D_VAL






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




!+101+###########################################################################!
!                                                                                !
!                  Initialize_Flat_Space_Guess_Values          !
!                                                                                !
!################################################################################!
SUBROUTINE Initialize_Flat_Space_Guess_Values()


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




END SUBROUTINE Initialize_Flat_Space_Guess_Values


END MODULE Thornado_Guess_BC_Module
