   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE SelfSimilar_Module                                                            !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
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
USE Poseidon_Constants_Module, &
            ONLY :  idp,             &
                    eps,             &
                    pi,              &
                    Grav_Constant_G, &
                    Speed_of_Light

USE CHIMERA_PARAMETERS,  &
            ONLY :  NUM_ENTRIES,       &
                    SELFSIM_R_VALS,    &
                    SELFSIM_POT_VALS,  &
                    SELFSIM_SHIFT_VALs,&
                    Analytic_Solution


IMPLICIT NONE

CONTAINS


!+101+###########################################################################!
!                                                                                !
!                  UNPACK_SELF_SIMILAR                                           !
!                                                                                !
!################################################################################!
SUBROUTINE UNPACK_SELF_SIMILAR( t_in, kappa, gamma,                     &
                                Num_Nodes, INPUT_R_QUAD,                &
                                NUM_R_ELEM, NUM_T_ELEM, NUM_P_ELEM,     &
                                Delta_r, r_locs,                        &
                                Input_E, Input_S, Input_Si )


REAL(KIND = idp),               INTENT(IN)                                  ::  t_in, kappa, gamma
INTEGER,    DIMENSION(1:3),     INTENT(IN)                                  ::  Num_Nodes
REAL(KIND = idp), DIMENSION(1:NUM_NODES(1)),INTENT(IN)                      ::  INPUT_R_QUAD
INTEGER,                        INTENT(IN)                                  ::  NUM_R_ELEM
INTEGER,                        INTENT(IN)                                  ::  NUM_T_ELEM
INTEGER,                        INTENT(IN)                                  ::  NUM_P_ELEM

REAL(KIND = idp), DIMENSION(1:NUM_R_ELEM), INTENT(INOUT)                         ::  Delta_R
REAL(KIND = idp), DIMENSION(0:NUM_R_ELEM), INTENT(INOUT)                         ::  r_locs

REAL(KIND = idp), DIMENSION(1:Num_Nodes(1)*Num_Nodes(2)*Num_Nodes(3),               &
                            0:NUM_R_ELEM-1, 0:NUM_T_ELEM-1, 0:NUM_P_ELEM-1    ),    &
                            INTENT(INOUT)                                               ::  Input_E

REAL(KIND = idp), DIMENSION(1:Num_Nodes(1)*Num_Nodes(2)*Num_Nodes(3),               &
                            0:NUM_R_ELEM-1, 0:NUM_T_ELEM-1, 0:NUM_P_ELEM-1    ),    &
                            INTENT(INOUT)                                               ::  Input_S

REAL(KIND = idp), DIMENSION(1:Num_Nodes(1)*Num_Nodes(2)*Num_Nodes(3),               &
                            0:NUM_R_ELEM-1, 0:NUM_T_ELEM-1, 0:NUM_P_ELEM-1,         &
                            1:3 ),  INTENT(INOUT)                                       ::  Input_Si


REAL(KIND = idp)                                                                        ::  t


REAL(KIND = idp), DIMENSION(:),ALLOCATABLE                                              ::  Input_X,    &
                                                                                            Input_D,    &
                                                                                            Input_V,    &
                                                                                            Input_M

REAL(KIND = idp)                                                                        ::  V_Factor,R_Factor

REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                                             ::  Enclosed_Mass
REAL(KIND = idp), DIMENSION(:),ALLOCATABLE                                              ::  Input_R


CHARACTER(LEN=128)                  :: line


INTEGER                                     :: NUM_LINES
INTEGER                                     :: CUR_LINE


INTEGER                             :: nread
INTEGER                             :: iskipp
INTEGER                             :: istat

101 FORMAT (a128)
111 FORMAT (1x,e16.10)
121 FORMAT (21x,e12.10)
131 FORMAT (37x,e13.10)
141 FORMAT (55x,e13.10)

nread = 42

OPEN(UNIT=nread, FILE='Input/YahilHomologousCollapse_Gm_130.dat', STATUS='OLD', IOSTAT=istat)
IF ( istat .NE. 0 ) THEN

    PRINT*,"Could not open 'Input/YahilHomologousCollapse_Gm_130.dat'. "

END IF
REWIND(nread)
NUM_LINES = 0
DO

    READ(nread, 101, IOSTAT=istat) line
    IF ( istat .NE. 0 ) EXIT
    NUM_LINES = NUM_LINES + 1

END DO



NUM_LINES = NUM_LINES -1

ALLOCATE( Input_X(1:NUM_LINES), Input_D(1:NUM_LINES), Input_V(1:NUM_LINES) )
ALLOCATE( Input_R(1:NUM_LINES) )
ALLOCATE( Input_M(1:NUM_LINES) )
ALLOCATE( Enclosed_Mass(1:NUM_LINES)  )

REWIND(nread)
CUR_LINE = 1
DO

    READ(nread, 101, IOSTAT=istat) line
    IF ( istat .NE. 0 ) EXIT

    IF ( line(1:3) == '#  ' ) CYCLE

    READ( line, 111) Input_X(CUR_LINE)
    READ( line, 121) Input_D(CUR_LINE)
    READ( line, 131) Input_V(CUR_LINE)
    READ( line, 141) Input_M(CUR_LINE)

!    PRINT*,"Input_X",CUR_LINE,Input_X(Cur_Line)
    CUR_LINE = CUR_LINE + 1

END DO

CLOSE(UNIT=nread,STATUS='keep',IOSTAT=istat)

t = t_in/1000.0_idp

R_Factor = SQRT(kappa)                                &
        *(Grav_Constant_G**((1.0_idp-gamma)/2.0_idp))       &
        *((t)**(2.0_idp-gamma))

Input_R = R_Factor*Input_X


Enclosed_Mass = kappa**(1.50_idp)                                   &
              * Grav_Constant_G**((1.0_idp-3.0_idp*gamma)/2.0_idp)  &
              * (t**(4.0_idp - 3.0_idp*gamma))                      &
              * Input_M



IF ( .TRUE. ) THEN

    PRINT*,"Using Yahil self-similar profile with following parameters."
    PRINT*,"Time = ",t
    PRINT*,"Kappa = ",Kappa
    PRINT*,"Gamma = ",Gamma

END IF

IF (.FALSE.) THEN

    PRINT*,"sqrt(kappa)",SQRT(kappa),"kappa**(1.50_idp)",kappa**(1.50_idp)
    PRINT*,"(Grav_Constant_G**((1.0_idp-gamma)/2.0_idp))",(Grav_Constant_G**((1.0_idp-gamma)/2.0_idp))
    PRINT*,"Grav_Constant_G**((1.0_idp-3.0_idp*gamma)/2.0_idp)",Grav_Constant_G**((1.0_idp-3.0_idp*gamma)/2.0_idp)
    PRINT*,"((t)**(2.0_idp-gamma))",((t)**(2.0_idp-gamma))
    PRINT*,"(t**(4.0_idp - 3.0_idp*gamma))",(t**(4.0_idp - 3.0_idp*gamma))
END IF

!PRINT*,"Input_R"
!PRINT*,Input_R
!PRINT*," "



CALL CONVERT_SELF_SIMILAR(  t, kappa, gamma,                        &
                            Num_Nodes, NUM_LINES, INPUT_R_QUAD,     &
                            NUM_R_ELEM, NUM_T_ELEM, NUM_P_ELEM,     &
                            Delta_R, r_locs,                        &
                            Input_R, Input_D, Input_V,              &
                            Input_E, Input_S, Input_Si              )


!PRINT*,"Input_E"
!PRINT*,Input_E
!PRINT*," "
!PRINT*,"Input_M"
!PRINT*,Input_M
!PRINT*," "
!PRINT*,"Enclosed_Mass"
!PRINT*,Enclosed_mass
!PRINT*," "



CALL CREATE_SELFSIM_NEWT_SOL( NUM_LINES, Input_R, Enclosed_Mass )
CALL CREATE_SELFSIM_SHIFT_SOL( Num_Nodes, NUM_R_ELEM, NUM_T_ELEM, NUM_P_ELEM, Input_Si, r_locs )



 5000 RETURN
END SUBROUTINE UNPACK_SELF_SIMILAR










!+201+###########################################################################!
!                                                                                !
!                  CONVERT_SELF_SIMILAR                                          !
!                                                                                !
!################################################################################!
SUBROUTINE CONVERT_SELF_SIMILAR( t, kappa, gamma,                        &
                                 Num_Nodes, NUM_LINES, INPUT_R_QUAD,     &
                                 NUM_R_ELEM, NUM_T_ELEM, NUM_P_ELEM,     &
                                 Delta_R, r_locs,                        &
                                 Input_R, Input_D, Input_V,              &
                                 Input_E, Input_S, Input_Si              )


REAL(KIND = idp),               INTENT(IN)                                  ::  t, kappa, gamma
INTEGER,    DIMENSION(1:3),     INTENT(IN)                                  ::  Num_Nodes
INTEGER,    INTENT(IN)                                                      ::  Num_LINES
REAL(KIND = idp), DIMENSION(1:NUM_NODES(1)),INTENT(IN)                      ::  INPUT_R_QUAD
INTEGER,                        INTENT(IN)                                  ::  NUM_R_ELEM
INTEGER,                        INTENT(IN)                                  ::  NUM_T_ELEM
INTEGER,                        INTENT(IN)                                  ::  NUM_P_ELEM

REAL(KIND = idp), DIMENSION(NUM_R_ELEM), INTENT(IN)                         ::  Delta_R
REAL(KIND = idp), DIMENSION(0:NUM_R_ELEM), INTENT(IN)                       ::  r_locs

REAL(KIND = idp), DIMENSION(1:NUM_LINES), INTENT(IN)                                    ::  Input_R,    &
                                                                                            Input_D,    &
                                                                                            Input_V

REAL(KIND = idp), DIMENSION(1:Num_Nodes(1)*Num_Nodes(2)*Num_Nodes(3),               &
                            0:NUM_R_ELEM-1, 0:NUM_T_ELEM-1, 0:NUM_P_ELEM-1    ),    &
                            INTENT(INOUT)                                               ::  Input_E

REAL(KIND = idp), DIMENSION(1:Num_Nodes(1)*Num_Nodes(2)*Num_Nodes(3),               &
                            0:NUM_R_ELEM-1, 0:NUM_T_ELEM-1, 0:NUM_P_ELEM-1    ),    &
                            INTENT(INOUT)                                               ::  Input_S

REAL(KIND = idp), DIMENSION(1:Num_Nodes(1)*Num_Nodes(2)*Num_Nodes(3),               &
                            0:NUM_R_ELEM-1, 0:NUM_T_ELEM-1, 0:NUM_P_ELEM-1,         &
                            1:3 ),  INTENT(INOUT)                                       ::  Input_Si



INTEGER                                                                     ::  re, te, pe, &
                                                                                rd, td, pd



INTEGER                                                                     ::  nd, line, line_min
REAL(KIND = idp), DIMENSION(0:1)                                            ::  xlocs
REAL(KIND = idp)                                                            ::  x
REAL(KIND = idp)                                                            ::  Density, Velocity

REAL(KIND = idp)                                                            ::  Pressure, Energy
REAL(KIND = idp)                                                            ::  vsqr, LF_sqr
REAL(KIND = idp)                                                            ::  E, Si, S

REAL(KIND = idp), DIMENSION(0:1)                                            :: LagPoly_Vals
REAL(KIND = idp), DIMENSION(0:NUM_NODES(1)-1)                               :: CUR_R_LOCS
REAL(KIND = idp)                                                            :: deltar_overtwo

REAL(KIND = idp)                                                            :: D_FACTOR, V_FACTOR
REAL(KIND = idp)                                                            :: r_Cur


REAL(KIND = idp)                                                            :: csqr
REAL(KIND = idp)                                                            :: Specific_Enthalpy

110 FORMAT (11X,A1,24X,A2,22X,A2,20X,A7,14X,A17,11X,A8,16X,A6)
111 FORMAT (11X,A1,24X,A1,22X,A2,22X,A1,22X,A3)
112 FORMAT (11X,A1,23X,A3,17X,A8)

113 FORMAT (E22.15,3X,E22.15,3X,E22.15,3X,E22.15,3X,E22.15,3X,E22.15,3X,E22.15)
114 FORMAT (E22.15,3x,ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15)
115 FORMAT (E22.15,3X,ES22.15,3X,ES22.15)

xlocs(0) = -1.0_idp
xlocs(1) = 1.0_idp

csqr = Speed_Of_Light*Speed_Of_Light

D_FACTOR = 1.0_idp /(Grav_Constant_G*t*t )

V_FACTOR = SQRT(kappa)                          &
         * Grav_Constant_G**((1.0_idp-gamma)/2)  &
         * t**(1.0_idp - gamma)

!PRINT*,"D_FACTOR",D_FACTOR,"V_FACTOR",V_FACTOR

r_cur = 0.0_idp
line_min = 1


OPEN( UNIT = 42, file = 'OUTPUT/Sources.out')
OPEN( UNIT = 43, file = 'OUTPUT/Base_Sources.out')


!WRITE(*,110)"r","Si","S","Density","Specific Enthalpy","Velocity","% of c%"
WRITE(42,111)"r","E","Si","S","v/c"
WRITE(43,112)"r","rho","velocity"



DO re = 0,NUM_R_ELEM-1

    deltar_overtwo = Delta_R(Re+1)/2.0_idp
    CUR_R_LOCS(:) = deltar_overtwo * (INPUT_R_QUAD(:)+1.0_idp) + r_cur
    r_cur = r_cur + Delta_R(re+1)

    DO rd = 0, NUM_NODES(1)-1
        DO line = line_min,NUM_LINES-1


            IF ( ( CUR_R_LOCS(rd) > Input_R(Line) ) .AND. ( CUR_R_LOCS(rd) <= Input_R(Line + 1) ) ) THEN


                line_min = line
                x = MAP_TO_X_SPACE(Input_R(Line),Input_R(Line+1),CUR_R_LOCS(rd))

                LagPoly_Vals = Lagrange_Poly(x, 1, xlocs)

                ! Interpolate Self-Similar Values to Input locations
                Density = (INPUT_D(line)*LagPoly_Vals(0) + INPUT_D(line+1)*LagPoly_Vals(1))*D_FACTOR
!                Velocity = (Input_V(line)*LagPoly_Vals(0) + INPUT_V(line+1)*LagPoly_Vals(1))*V_FACTOR
                Velocity = 0.0_idp

!                PRINT*,"**********************  Velocity Zeroed  **************************"
!                Velocity = 0.0_idp

                ! Calculate Usable Quantities
                Pressure = kappa * Density**gamma
                Energy = Pressure/(gamma - 1.0_idp)
                Specific_Enthalpy = csqr + (Energy + Pressure)/Density

                vsqr = Velocity*Velocity
                LF_sqr = 1.0_idp/(1.0_idp - vsqr/csqr)



                !  Calculate CFA Input Values
                E  = Density*Specific_Enthalpy*LF_sqr - Pressure
                Si = Density*Specific_Enthalpy*LF_sqr*Velocity/csqr
                S  = Density*Specific_Enthalpy*LF_sqr*vsqr/(csqr) + 3.0_idp * Pressure



!                Si = 0.0_idp
!                S  = 0.0_idp  
!                PRINT*,"*****************  S,Si Zeroed  ********************"

!                WRITE(*,113)CUR_R_LOCS(0),Si, S, Density, Specific_Enthalpy,Velocity,Velocity/Speed_of_Light
                WRITE(42,114)CUR_R_LOCS(rd),E,Si, S,Velocity/Speed_of_Light
                WRITE(43,115)CUR_R_LOCS(rd),Density,Velocity

                DO te = 0,NUM_T_ELEM-1
                    DO pe = 0,NUM_P_ELEM-1

                        DO td = 0,NUM_NODES(2)-1
                            DO pd = 0,NUM_NODES(3)-1

                                nd = pd*NUM_NODES(2)*NUM_NODES(1)   &
                                   + td*NUM_NODES(1)                &
                                   + rd + 1

                                Input_E(nd, re,te,pe) = E
                                Input_Si(nd, re, te, pe, 1) = Si
                                Input_S(nd, re, te, pe) = S

                            END DO ! pd
                        END DO ! td

                    END DO ! pe
                END DO ! te


            END IF
        END DO ! Line
    END DO ! rd
END DO ! re


!PRINT*,"Input_D"
!PRINT*,Input_D
!PRINT*," "
!PRINT*," "
!PRINT*,"Input_E"
!PRINT*,Input_E
!PRINT*," "
!PRINT*," "
!PRINT*,"Input_Si"
!PRINT*,Input_Si(:,:,:,:,1)



END SUBROUTINE CONVERT_SELF_SIMILAR











!+201+###########################################################################!
!                                                                                !
!                  CREATE_SELFSIM_NEWT_SOL                                       !
!                                                                                !
!################################################################################!
SUBROUTINE CREATE_SELFSIM_NEWT_SOL( Line_Count, R_Values, Enclosed_Mass )

INTEGER, INTENT(IN)                                             :: Line_Count
REAL(KIND = idp), DIMENSION(1:Line_Count), INTENT(IN)           :: R_Values
REAL(KIND = idp), DIMENSION(1:Line_Count), INTENT(IN)           :: Enclosed_Mass

INTEGER     :: i

PRINT*,"Line_Count",Line_Count

NUM_ENTRIES = Line_Count
ALLOCATE( SELFSIM_R_VALS(0:NUM_ENTRIES) )
ALLOCATE( SELFSIM_POT_VALS(0:NUM_ENTRIES) )


SELFSIM_R_VALS(0) = 0.0_idp
SELFSIM_R_VALS(1:NUM_ENTRIES) = R_Values(1:NUM_ENTRIES)
!SELFSIM_POT_VALS = - Grav_Constant_G*Enclosed_Mass(:)/R_Values(:)

SELFSIM_POT_VALS(NUM_ENTRIES) = -GRAV_Constant_G*Enclosed_Mass(Num_Entries)/R_Values(Num_Entries)

DO i = NUM_ENTRIES-1,1,-1

     SELFSIM_POT_VALS(i) = SELFSIM_POT_VALS(i+1) - Grav_Constant_G*Enclosed_Mass(i)             &
                                                 * (SELFSIM_r_Vals(i+1)-SELFSIM_R_Vals(i))                  &
                                                 / (SELFSIM_R_Vals(i)*SELFSIM_R_Vals(i))

END DO
!SELFSIM_POT_VALS(0) = SELFSIM_POT_VALS(1)/2.0_idp

SELFSIM_POT_VALS(0) = SELFSIM_POT_VALS(1) - 3*Grav_Constant_G*Enclosed_Mass(1)           &
                                            /(2*SELFSIM_R_VALS(1))

!PRINT*,"SELFSIM_R_VALS"
!PRINT*,SELFSIM_R_VALS
!PRINT*," "
!PRINT*,"Enclosed_Mass"
!PRINT*,Enclosed_Mass
!PRINT*," "
!PRINT*,"SELFSIM_POT_VALS"
!PRINT*,SELFSIM_POT_VALS
!PRINT*," "

!CALL SELFSIM_NEWT_SUB( SELFSIM_R_VALS(100) + (SELFSIM_R_VALS(101)-SELFSIM_R_VALS(100))/2.0_idp)


END SUBROUTINE CREATE_SELFSIM_NEWT_SOL



!+201+###########################################################################!
!                                                                                !
!                         SELFSIM_NEWT_SOL                                       !
!                                                                                !
!################################################################################!
FUNCTION SELFSIM_NEWT_SOL( r, theta, phi )

REAL(KIND = idp), INTENT(IN)     :: r, theta, phi
REAL(KIND = idp)                 :: SELFSIM_NEWT_SOL

INTEGER                          :: cur_entry
INTEGER                          :: i

REAL(KIND = idp)                 :: deltar


DO i = 0,NUM_ENTRIES-1

   IF ( r == 0 ) THEN

       cur_entry = 0

   
   ELSE IF (( r > SELFSIM_R_VALS(i) ) .AND. ( r <= SELFSIM_R_VALS(i+1) ) ) THEN

       cur_entry = i

   END IF
 
END DO


deltar = SELFSIM_R_VALS(cur_entry+1) - SELFSIM_R_VALS(cur_entry)

SELFSIM_NEWT_SOL = (1.0_idp/deltar)    &
                 *( SELFSIM_POT_VALS(cur_entry)*(SELFSIM_R_VALS(cur_entry+1) - r)         &
                   +SELFSIM_POT_VALS(cur_entry+1)*(r - SELFSIM_R_VALS(cur_entry))         )



END FUNCTION SELFSIM_NEWT_SOL




!+401+###########################################################################!
!                                                                                !
!              CREATE_SELFSIM_SHIFT_SOL                                          !
!                                                                                !
!################################################################################!
SUBROUTINE CREATE_SELFSIM_SHIFT_SOL( Num_Nodes,                                  &
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


csqr = Speed_of_Light*Speed_of_Light


ALLOCATE( SELFSIM_SHIFT_VALS(0:NUM_R_ELEM) ) 



ALLOCATE( x_locs(1:Ord) )
ALLOCATE( ri_locs(1:Ord) )
ALLOCATE( wi(1:Ord)  )
ALLOCATE( AlphaPsi(1:Ord) )

ALLOCATE( rij_locs(1:Ord,1:Ord) )
ALLOCATE( PSI_10(1:Ord,1:Ord)  )
ALLOCATE( Sr_New(1:Ord,1:Ord)  )

CALL Initialize_LG_Quadrature( Ord, x_locs, wi )


SELFSIM_SHIFT_VALS(0) = 0.0_idp
DO re = 0,NUM_R_ELEM-1

   ! Calculate the r locations for the Outer Integral's Quadrature Points !
   ri_locs(:) = (r_locs(re+1)-r_locs(re))/2.0_idp * ( x_locs(:) + 1.0_idp ) + r_locs(re)

   ! Calculate the Alpha Psi values at each of the Outer Integral's Quadrature Points !
   DO i = 1,Ord
      AlphaPsi(i) =  1.0_idp + 0.5_idp*SELFSIM_NEWT_SOL(ri_locs(i),0.0_idp,0.0_idp)/csqr
   END DO 


   DO i = 1,Ord

      ! Calculate the Quadrature Points for each of the Inner Integrals
      rij_locs(:,i) = (ri_locs(i) - 0.0_idp)/2.0_idp *(x_locs(:) + 1.0_idp) + 0.0_idp


      ! Calculate Psi^10 values at each of the Inner Quadrature Points
      DO j = 1,Ord

           Psi = 1.0_idp - 0.5_idp*SELFSIM_NEWT_SOL(rij_locs(j,i),0.0_idp,0.0_idp)/csqr
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

      SELFSIM_SHIFT_VALS(re+1) = (3.0_idp/2.0_idp)                        &
                               * 8.0_idp*pi* Grav_Constant_G/(csqr*csqr)  &
                               * r_locs(re+1)                             &
                               * ( r_locs(re+1) - R_locs(re))/2.0_idp     &
                               * Outer_Int

   ELSE

      SELFSIM_SHIFT_VALS(re+1) = r_locs(re+1)/r_locs(re)                  &
                               * SELFSIM_SHIFT_VALS(re)                   &
                               + (3.0_idp/2.0_idp)                        &
                               * 8.0_idp*pi* Grav_Constant_G/(csqr*csqr)  &
                               * r_locs(re+1)                             &
                               * ( r_locs(re+1) - r_locs(re))/2.0_idp     &
                               * Outer_Int

   END IF


END DO ! re loop



END SUBROUTINE CREATE_SELFSIM_SHIFT_SOL




!+201+###########################################################################!
!                                                                                !
!                         SELFSIM_SHIFT_SOL                                       !
!                                                                                !
!################################################################################!
FUNCTION SELFSIM_SHIFT_SOL( r, r_locs,  NUM_R_ELEM )

REAL(KIND = idp), INTENT(IN)                              :: r
INTEGER, INTENT(IN)                                       :: NUM_R_ELEM

REAL(KIND = idp), DIMENSION(0:NUM_R_ELEM), INTENT(IN)     :: r_locs


REAL(KIND = idp)                 :: SELFSIM_SHIFT_SOL

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

SELFSIM_SHIFT_SOL = (1.0_idp/deltar)    &
                 *( SELFSIM_SHIFT_VALS(cur_entry)*(r_locs(cur_entry+1) - r)         &
                   +SELFSIM_SHIFT_VALS(cur_entry+1)*(r - r_locs(cur_entry))         )



END FUNCTION SELFSIM_SHIFT_SOL


!+201+###########################################################################!
!                                                                                !
!                  CREATE_SELFSIM_NEWT_SOL                                       !
!                                                                                !
!################################################################################!
SUBROUTINE SELFSIM_NEWT_SUB( r )

REAL(KIND = idp), INTENT(IN)     :: r
REAL(KIND = idp)                 :: SELFSIM_NEWT_SOL

INTEGER                          :: cur_entry
INTEGER                          :: i

REAL(KIND = idp)                 :: deltar


DO i = 1,NUM_ENTRIES-1

   IF (( r > SELFSIM_R_VALS(i) ) .AND. ( r <= SELFSIM_R_VALS(i+1) ) ) THEN

       cur_entry = i

   END IF

END DO


SELFSIM_NEWT_SOL = (1.0_idp/(SELFSIM_R_VALS(cur_entry+1) - SELFSIM_R_VALS(cur_entry)))    &
                 *( SELFSIM_POT_VALS(cur_entry)*(SELFSIM_R_VALS(cur_entry+1) - r)         &
                   +SELFSIM_POT_VALS(cur_entry+1)*(r - SELFSIM_R_VALS(cur_entry))         )


END SUBROUTINE SELFSIM_NEWT_SUB





 !+101+################################################################!
!                                                                       !
!   Lagrange_Poly - Calculates the value of the Lagrange Polynomial     !
!                                                                       !
!                                                                       !
 !#####################################################################!
PURE FUNCTION Lagrange_Poly(x, Ord, xlocs)


INTEGER, INTENT(IN)                                 :: Ord
REAL(KIND = idp), INTENT(IN)                        :: x
REAL(KIND = idp), INTENT(IN), DIMENSION(0:Ord)      :: xlocs


INTEGER                                             :: i,j
REAL(KIND = idp), DIMENSION(0:Ord)                  :: tmp
REAL(KIND = idp), DIMENSION(0:Ord)                  :: Lagrange_Poly



tmp(:) = 1.0_idp

DO j = 0,Ord
    DO i = 0,Ord
        IF (i .NE. j) THEN
            tmp(i) = tmp(i) * (x - xlocs(j))/(xlocs(i) - xlocs(j))
        END IF
    END DO
END DO



Lagrange_Poly = tmp


END FUNCTION Lagrange_Poly




!+301+##########################################################!
!                                                               !
!      Map_To_X_Space - maps r value between ra, and rb to x    !
!                   space such that x in [-1,1].                !
!                                                               !
!###############################################################!
PURE ELEMENTAL FUNCTION Map_To_X_Space(ra, rb, r)

REAL(KIND = idp)                            ::  Map_To_X_Space
REAL(KIND = idp), intent(in)                ::  ra, rb
REAL(KIND = idp), intent(in)                ::  r

Map_To_X_Space = (2.0_idp*(r - ra))/(rb - ra) - 1.0_idp


END FUNCTION Map_To_X_Space







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









END MODULE SelfSimilar_Module

