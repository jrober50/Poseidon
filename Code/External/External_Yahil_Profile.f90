   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE External_Yahil_Profile_Module                                                !##!
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
USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Poseidon_Numbers_Module, &
            ONLY :  pi, eps

USE Poseidon_Message_Routines_Module, &
            ONLY :  Driver_Init_Message,    &
                    Warning_Message

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag

USE Poseidon_Units_Module, &
            ONLY :  Grav_Constant_G,    &
                    C_Square,           &
                    GR_Source_Scalar,   &
                    Centimeter,         &
                    Second,             &
                    Millisecond,         &
                    Erg,                &
                    Gram

USE Variables_MPI, &
            ONLY :  myID_Poseidon,      &
                    MasterID_Poseidon

USE Variables_External, &
            ONLY :  NUM_ENTRIES,        &
                    SELFSIM_R_VALS,     &
                    SELFSIM_POT_VALS,   &
                    SELFSIM_SHIFT_VALs, &
                    SELFSIM_V_SWITCH,   &
                    OUTPUT_PRIMATIVES_FLAG

USE Variables_Functions, &
            ONLY :  Potential_Solution,             &
                    Shift_Solution

USE Allocation_Yahil_Profile, &
            ONLY :  Allocate_Yahil_Profile,         &
                    Deallocate_Yahil_Profile

USE IO_Output_Sources_Module, &
            ONLY :  Output_Primatives

USE Maps_Quadrature, &
            ONLY :  Quad_Map

USE Timer_Routines_Module, &
            ONLY :  TimerStart,                     &
                    TimerStop

USE Timer_Variables_Module, &
            ONLY :  Timer_Core_Init_Test_Problem

USE Flags_IO_Module, &
            ONLY :  lPF_IO_Flags,           &
                    iPF_IO_Print_Setup


IMPLICIT NONE

CONTAINS


!+101+###########################################################################!
!                                                                                !
!                  Initialize_Yahil_Sources                                           !
!                                                                                !
!################################################################################!
SUBROUTINE Initialize_Yahil_Sources( t_in, kappa, gamma, ecc,                &
                                    Num_Nodes, INPUT_R_QUAD, INPUT_T_QUAD,  &
                                    NUM_R_ELEM, NUM_T_ELEM, NUM_P_ELEM,     &
                                    Delta_r, r_locs, t_locs,                &
                                    Input_E, Input_S, Input_Si )


REAL(idp),                          INTENT(IN)                      ::  t_in, kappa, gamma, ecc
INTEGER,  DIMENSION(1:3),           INTENT(IN)                      ::  Num_Nodes
REAL(idp),DIMENSION(1:NUM_NODES(1)),INTENT(IN)                      ::  INPUT_R_QUAD
REAL(idp),DIMENSION(1:NUM_NODES(2)),INTENT(IN)                      ::  INPUT_T_QUAD

INTEGER,                            INTENT(IN)                      ::  NUM_R_ELEM
INTEGER,                            INTENT(IN)                      ::  NUM_T_ELEM
INTEGER,                            INTENT(IN)                      ::  NUM_P_ELEM

REAL(idp), DIMENSION(1:NUM_R_ELEM), INTENT(IN)                      ::  Delta_R
REAL(idp), DIMENSION(0:NUM_R_ELEM), INTENT(IN)                      ::  r_locs
REAL(idp), DIMENSION(0:NUM_T_ELEM), INTENT(IN)                      ::  t_locs


REAL(idp), DIMENSION( 1:Num_Nodes(1)*Num_Nodes(2)*Num_Nodes(3),     &
                      0:NUM_R_ELEM-1,                               &
                      0:NUM_T_ELEM-1,                               &
                      0:NUM_P_ELEM-1    ), INTENT(INOUT)            ::  Input_E

REAL(idp), DIMENSION( 1:Num_Nodes(1)*Num_Nodes(2)*Num_Nodes(3),     &
                      0:NUM_R_ELEM-1,                               &
                      0:NUM_T_ELEM-1,                               &
                      0:NUM_P_ELEM-1    ), INTENT(INOUT)            ::  Input_S

REAL(idp), DIMENSION( 1:Num_Nodes(1)*Num_Nodes(2)*Num_Nodes(3),     &
                      0:NUM_R_ELEM-1,                               &
                      0:NUM_T_ELEM-1,                               &
                      0:NUM_P_ELEM-1,1:3), INTENT(INOUT)            ::  Input_Si


REAL(idp)                                                           ::  t


REAL(idp), DIMENSION(:),ALLOCATABLE                                 ::  Input_X,    &
                                                                        Input_D,    &
                                                                        Input_V,    &
                                                                        Input_M

REAL(idp)                                                           ::  R_Factor
REAL(idp)                                                           ::  Kappa_wUnits


REAL(idp), DIMENSION(:), ALLOCATABLE                                ::  Enclosed_Mass
REAL(idp), DIMENSION(:),ALLOCATABLE                                 ::  Input_R


CHARACTER(LEN=128)                          :: line

INTEGER                                     :: NUM_LINES
INTEGER                                     :: CUR_LINE


INTEGER                                     :: nread
INTEGER                                     :: istat

IF ( Verbose_Flag ) CALL Driver_Init_Message('Initializing Source Profile : Yahil-Lattimer Self-Similar Collapse.')
CALL TimerStart( Timer_Core_Init_Test_Problem )




101 FORMAT (a128)


nread = 42
OPEN(UNIT=nread, FILE='../../Input/YahilHomologousCollapse_Gm_130.dat', STATUS='OLD', IOSTAT=istat)
IF ( istat .NE. 0 ) THEN
    PRINT*,"Could not open 'Input/YahilHomologousCollapse_Gm_130.dat'. "
END IF
REWIND(nread)
NUM_LINES = 0
DO
    READ(nread, *, IOSTAT=istat) line
    IF ( istat .NE. 0 ) EXIT
    NUM_LINES = NUM_LINES + 1
END DO
NUM_LINES = NUM_LINES -1



ALLOCATE( Input_X(1:NUM_LINES), Input_D(1:NUM_LINES), Input_V(1:NUM_LINES) )
ALLOCATE( Input_R(1:NUM_LINES) )
ALLOCATE( Input_M(1:NUM_LINES) )
ALLOCATE( Enclosed_Mass(1:NUM_LINES)  )

NUM_ENTRIES = NUM_LINES

CALL Allocate_Yahil_Profile(Num_Entries)




REWIND(nread)
CUR_LINE = 1
READ(nread,*)
DO

    READ(nread, 101, IOSTAT=istat) line
    IF ( istat .NE. 0 ) EXIT


    READ(line,*) Input_X(CUR_LINE), Input_D(CUR_LINE), Input_V(CUR_LINE), Input_M(CUR_LINE)
    CUR_LINE = CUR_LINE + 1
    
END DO

CLOSE(UNIT=nread,STATUS='keep',IOSTAT=istat)








IF ( .TRUE. ) THEN
    Kappa_wUnits = Kappa*((Erg/Centimeter**3)/(Gram/Centimeter**3)**Gamma)
ELSE
    Kappa_wUnits = 18.394097187796024_idp
    PRINT*,"Kappa_wUnits over written.",Kappa_wUnits
END IF



t = t_in*Millisecond


R_Factor = SQRT(Kappa_wUnits)                                &
        *(Grav_Constant_G**((1.0_idp-gamma)/2.0_idp))       &
        *((t)**(2.0_idp-gamma))


Input_R = R_Factor*Input_X


Enclosed_Mass = Kappa_wUnits**(1.50_idp)                                   &
              * Grav_Constant_G**((1.0_idp-3.0_idp*gamma)/2.0_idp)  &
              * (t**(4.0_idp- 3.0_idp*gamma))                      &
              * Input_M


IF ( myID_Poseidon == MasterID_Poseidon ) THEN
IF ( lPF_IO_Flags(iPF_IO_Print_Setup)   ) THEN
    WRITE(*,'(A)')'------------- Test Parameters ----------------'
    WRITE(*,'(A)')' Source Configuration : Yahil Self-Similar Collapse Profile'
    WRITE(*,'(A,ES12.5,A)') ' - Yahil Time      : ', t_in,' ms'
    WRITE(*,'(A,ES12.5)')   ' - Kappa           : ', Kappa_wUnits
    WRITE(*,'(A,ES12.5)')   ' - Gamma           : ', Gamma
    WRITE(*,'(A,ES12.5,A)')   ' - Central Density : ', Input_D(1)/(Grav_Constant_G*t*t )/(Gram/Centimeter**3),' g/cm^3'
    WRITE(*,'(/)')
END IF
END IF



IF ( SELFSIM_V_SWITCH == 1 ) THEN
    CALL Warning_Message("The Yahil profile has had its velocity zeroed.")
    Input_V = 0.0_idp
END IF



CALL CONVERT_SELF_SIMILAR_3D(  t, Kappa_wUnits, gamma, ecc,                   &
                                Num_Nodes, NUM_LINES,                   &
                                INPUT_R_QUAD, INPUT_T_QUAD,             &
                                NUM_R_ELEM, NUM_T_ELEM, NUM_P_ELEM,     &
                                Delta_R, r_locs, t_locs,                &
                                Input_D, Input_V, Input_X,             &
                                Input_E, Input_S, Input_Si              )





CALL CREATE_SELFSIM_NEWT_SOL( NUM_LINES, Input_R, Enclosed_Mass )
CALL CREATE_SELFSIM_SHIFT_SOL( Num_Nodes, NUM_R_ELEM, NUM_T_ELEM, NUM_P_ELEM, Input_Si, r_locs )

Potential_Solution => SELFSIM_NEWT_SOL
Shift_Solution => SELFSIM_SHIFT_SOL




!CALL Deallocate_Yahil_Profile()

CALL TimerStop( Timer_Core_Init_Test_Problem )

END SUBROUTINE Initialize_Yahil_Sources









!+101+###########################################################################!
!                                                                                !
!                  Initialize_Yahil_Sources                                           !
!                                                                                !
!################################################################################!
SUBROUTINE Initialize_Yahil_Density( t_in, kappa, gamma, ecc,               &
                                    Num_Nodes, INPUT_R_QUAD,                &
                                    NUM_R_ELEM, NUM_T_ELEM, NUM_P_ELEM,     &
                                    Delta_r, r_locs,                        &
                                    Density )


REAL(idp),               INTENT(IN)                              ::  t_in, kappa, gamma, ecc
INTEGER,    DIMENSION(1:3),     INTENT(IN)                        ::  Num_Nodes
REAL(idp), DIMENSION(1:NUM_NODES(1)),INTENT(IN)                  ::  INPUT_R_QUAD

INTEGER,                        INTENT(IN)                              ::  NUM_R_ELEM
INTEGER,                        INTENT(IN)                              ::  NUM_T_ELEM
INTEGER,                        INTENT(IN)                              ::  NUM_P_ELEM

REAL(idp), DIMENSION(1:NUM_R_ELEM), INTENT(IN)                   ::  Delta_R
REAL(idp), DIMENSION(0:NUM_R_ELEM), INTENT(IN)                   ::  r_locs


REAL(idp), DIMENSION(1:Num_Nodes(1)*Num_Nodes(2)*Num_Nodes(3),   &
                            0:NUM_R_ELEM-1,                             &
                            0:NUM_T_ELEM-1,                             &
                            0:NUM_P_ELEM-1    ),  INTENT(INOUT)         ::  Density


REAL(idp)                                                                        ::  t


REAL(idp), DIMENSION(:),ALLOCATABLE                                              ::  Input_X,    &
                                                                                            Input_D,    &
                                                                                            Input_V,    &
                                                                                            Input_M

REAL(idp)                                                                        ::  R_Factor
REAL(idp)                                                                        ::  Kappa_wUnits


REAL(idp), DIMENSION(:), ALLOCATABLE                                             ::  Enclosed_Mass
REAL(idp), DIMENSION(:),ALLOCATABLE                                              ::  Input_R


CHARACTER(LEN=128)                  :: line


INTEGER                                     :: NUM_LINES
INTEGER                                     :: CUR_LINE


INTEGER                             :: nread
INTEGER                             :: istat


IF ( Verbose_Flag ) CALL Driver_Init_Message('Initializing Density Profile : Yahil-Lattimer Self-Similar Collapse.')
CALL TimerStart( Timer_Core_Init_Test_Problem )




101 FORMAT (a128)


nread = 42
OPEN(UNIT=nread, FILE='../../Input/YahilHomologousCollapse_Gm_130.dat', STATUS='OLD', IOSTAT=istat)
IF ( istat .NE. 0 ) THEN
    PRINT*,"Could not open 'Input/YahilHomologousCollapse_Gm_130.dat'. "
END IF
REWIND(nread)
NUM_LINES = 0
DO
    READ(nread, *, IOSTAT=istat) line
    IF ( istat .NE. 0 ) EXIT
    NUM_LINES = NUM_LINES + 1
END DO
NUM_LINES = NUM_LINES -1



ALLOCATE( Input_X(1:NUM_LINES), Input_D(1:NUM_LINES), Input_V(1:NUM_LINES) )
ALLOCATE( Input_R(1:NUM_LINES) )
ALLOCATE( Input_M(1:NUM_LINES) )
ALLOCATE( Enclosed_Mass(1:NUM_LINES)  )

NUM_ENTRIES = NUM_LINES

CALL Allocate_Yahil_Profile(Num_Entries)




REWIND(nread)
CUR_LINE = 1
READ(nread,*)
DO

    READ(nread, 101, IOSTAT=istat) line
    IF ( istat .NE. 0 ) EXIT


    READ(line,*) Input_X(CUR_LINE), Input_D(CUR_LINE), Input_V(CUR_LINE), Input_M(CUR_LINE)
    CUR_LINE = CUR_LINE + 1
    
END DO

CLOSE(UNIT=nread,STATUS='keep',IOSTAT=istat)




IF ( .TRUE. ) THEN
    Kappa_wUnits = Kappa*((Erg/Centimeter**3)/(Gram/Centimeter**3)**Gamma)
ELSE
    Kappa_wUnits = 18.394097187796024_idp
    PRINT*,"Kappa_wUnits over written.",Kappa_wUnits
END IF



t = t_in*Millisecond

R_Factor = SQRT(Kappa_wUnits)                                &
        *(Grav_Constant_G**((1.0_idp-gamma)/2.0_idp))       &
        *((t)**(2.0_idp-gamma))


Input_R = R_Factor*Input_X


Enclosed_Mass = Kappa_wUnits**(1.50_idp)                                   &
              * Grav_Constant_G**((1.0_idp-3.0_idp*gamma)/2.0_idp)  &
              * (t**(4.0_idp- 3.0_idp*gamma))                      &
              * Input_M


IF ( myID_Poseidon == MasterID_Poseidon ) THEN
IF ( lPF_IO_Flags(iPF_IO_Print_Setup)   ) THEN
    WRITE(*,'(A)')'------------- Test Parameters ----------------'
    WRITE(*,'(A)')' Source Configuration : Yahil Self-Similar Collapse Profile'
    WRITE(*,'(A,ES12.5,A)') ' - Yahil Time      : ', t_in,' ms'
    WRITE(*,'(A,ES12.5)')   ' - Kappa           : ', Kappa_wUnits
    WRITE(*,'(A,ES12.5)')   ' - Gamma           : ', Gamma
    WRITE(*,'(A,ES12.5,A)')   ' - Central Density : ', Input_D(1)/(Grav_Constant_G*t*t )/(Gram/Centimeter**3),' g/cm^3'
    WRITE(*,'(/)')
END IF
END IF




CALL Calculate_Yahil_Density(   t, kappa, gamma,                        &
                                Num_Nodes, NUM_LINES,                   &
                                INPUT_R_QUAD,                           &
                                NUM_R_ELEM, NUM_T_ELEM, NUM_P_ELEM,     &
                                Delta_R, r_locs,                        &
                                Input_D,                                &
                                Input_X,                                &
                                Density                                )




CALL CREATE_SELFSIM_NEWT_SOL( NUM_LINES, Input_R, Enclosed_Mass )

Potential_Solution => SELFSIM_NEWT_SOL


CALL TimerStop( Timer_Core_Init_Test_Problem )

END SUBROUTINE Initialize_Yahil_Density








!+201+###########################################################################!
!                                                                                !
!                  CONVERT_SELF_SIMILAR_3D                                       !
!                                                                                !
!################################################################################!
SUBROUTINE CONVERT_SELF_SIMILAR_3D( t, kappa, gamma, ecc,                  &
                                    Num_Nodes, NUM_LINES,                   &
                                    INPUT_R_QUAD, INPUT_T_QUAD,             &
                                    NUM_R_ELEM, NUM_T_ELEM, NUM_P_ELEM,     &
                                    Delta_R, r_locs, t_locs,                &
                                    Input_D, Input_V, Input_X,              &
                                    Input_E, Input_S, Input_Si              )


REAL(idp),               INTENT(IN)                                  ::  t, kappa, gamma, ecc
INTEGER,    DIMENSION(1:3),     INTENT(IN)                                  ::  Num_Nodes
INTEGER,    INTENT(IN)                                                      ::  Num_LINES
REAL(idp), DIMENSION(1:NUM_NODES(1)),INTENT(IN)                      ::  INPUT_R_QUAD
REAL(idp), DIMENSION(1:NUM_NODES(2)),INTENT(IN)                      ::  INPUT_T_QUAD

INTEGER,                        INTENT(IN)                                  ::  NUM_R_ELEM
INTEGER,                        INTENT(IN)                                  ::  NUM_T_ELEM
INTEGER,                        INTENT(IN)                                  ::  NUM_P_ELEM

REAL(idp), DIMENSION(NUM_R_ELEM),   INTENT(IN)                       ::  Delta_R
REAL(idp), DIMENSION(0:NUM_R_ELEM), INTENT(IN)                       ::  r_locs
REAL(idp), DIMENSION(0:NUM_T_ELEM), INTENT(IN)                       ::  t_locs


REAL(idp), DIMENSION(1:NUM_LINES), INTENT(IN)                        ::  Input_D,    &
                                                                                Input_V,    &
                                                                                Input_X

REAL(idp), DIMENSION(1:Num_Nodes(1)*Num_Nodes(2)*Num_Nodes(3),       &
                            0:NUM_R_ELEM-1,                                 &
                            0:NUM_T_ELEM-1,                                 &
                            0:NUM_P_ELEM-1    ),    INTENT(INOUT)           ::  Input_E

REAL(idp), DIMENSION(1:Num_Nodes(1)*Num_Nodes(2)*Num_Nodes(3),       &
                            0:NUM_R_ELEM-1,                                 &
                            0:NUM_T_ELEM-1,                                 &
                            0:NUM_P_ELEM-1    ),    INTENT(INOUT)           ::  Input_S

REAL(idp), DIMENSION(1:Num_Nodes(1)*Num_Nodes(2)*Num_Nodes(3),       &
                            0:NUM_R_ELEM-1,                                 &
                            0:NUM_T_ELEM-1,                                 &
                            0:NUM_P_ELEM-1, 1:3 ),  INTENT(INOUT)           ::  Input_Si


INTEGER                                                                     ::  re, te, pe, &
                                                                                rd, td, pd

INTEGER                                                                     ::  nd, line, line_min
REAL(idp), DIMENSION(0:1)                                            ::  xlocs
REAL(idp)                                                            ::  x
REAL(idp)                                                            ::  Density, Velocity

REAL(idp)                                                            ::  Pressure, Energy
REAL(idp)                                                            ::  vsqr, LF_sqr
REAL(idp)                                                            ::  E, Si, S

REAL(idp), DIMENSION(0:1)                                            ::  LagPoly_Vals
REAL(idp), DIMENSION(1:NUM_NODES(1))                                 ::  CUR_R_LOCS

REAL(idp)                                                            ::  D_FACTOR,               &
                                                                                V_FACTOR,               &
                                                                                X_Factor,               &
                                                                                M_FACTOR

REAL(idp)                                                            ::  ecc_sqr, ooomes
REAL(idp)                                                            ::  Specific_Enthalpy

INTEGER                                                                     ::  Num_Radial_Points
INTEGER                                                                     ::  Num_Theta_Points

REAL(idp), DIMENSION(:), ALLOCATABLE                                 ::  Density_Holder,         &
                                                                                Velocity_Holder

REAL(idp)                                                            ::  xloc
REAL(idp)                                                            ::  E_Units

REAL(idp), DIMENSION(:), ALLOCATABLE                                 ::  DX_Holder
REAL(idp), DIMENSION(:), ALLOCATABLE                                 ::  VX_Holder

REAL(idp), DIMENSION(:), ALLOCATABLE                                 ::  r_Holder


Num_Radial_Points = NUM_R_ELEM*NUM_NODES(1)
Num_Theta_Points = NUM_T_ELEM*NUM_NODES(2)
ALLOCATE( Density_Holder(1:Num_Radial_Points) )
ALLOCATE( Velocity_Holder(1:Num_Radial_Points) )
ALLOCATE( DX_Holder(1:Num_Radial_Points) )
ALLOCATE( VX_Holder(1:Num_Radial_Points) )
ALLOCATE( r_Holder(1:Num_Radial_Points) )

E_Units = Erg/Centimeter**3

xlocs(0) = -1.0_idp
xlocs(1) = 1.0_idp

ecc_sqr = ecc*ecc
ooomes = 1.0_idp/(1.0_idp- ecc_sqr)


D_FACTOR = 1.0_idp/(Grav_Constant_G*t*t )

V_FACTOR = SQRT(kappa)                          &
         * Grav_Constant_G**((1.0_idp-gamma)/2.0_idp)  &
         * t**(1.0_idp- gamma)

X_Factor = kappa**(-0.5_idp)                               &
        *(Grav_Constant_G**((gamma-1.0_idp)/2.0_idp))       &
        *((t)**(gamma-2.0_idp))

M_FACTOR = kappa**(3.0_idp/2.0_idp)                             &
         * Grav_Constant_G**((1.0_idp-3.0_idp*gamma)/2.0_idp)   &
         * t**(4.0_idp- 3.0_idp* gamma )




DO pe = 0,NUM_P_ELEM-1
DO te = 0,NUM_T_ELEM-1

    line_min = 1
    DO re = 0,NUM_R_ELEM-1
        CUR_R_LOCS(:) = Delta_R(Re+1) * (INPUT_R_QUAD(:)+0.5_idp) + r_locs(re)
    
        DO rd = 1,NUM_NODES(1)
            xloc = CUR_R_LOCS(rd)*X_Factor
!            CALL Find_Line_SUB(xloc, Input_X, NUM_LINES)
            line = Find_Line(xloc, Input_X, NUM_LINES)


            x = MAP_TO_X_SPACE(Input_X(Line),Input_X(Line+1),xloc)
            IF ( x > 1 ) THEN
                x = 1
            END IF
            LagPoly_Vals = Lagrange_Poly(x, 1, xlocs)


            ! Interpolate Self-Similar Values to Input locations
            Density  = (INPUT_D(line)*LagPoly_Vals(0) + INPUT_D(line+1)*LagPoly_Vals(1))*D_FACTOR
            Velocity = (INPUT_V(line)*LagPoly_Vals(0) + INPUT_V(line+1)*LagPoly_Vals(1))*V_FACTOR

         
            DX_Holder(re*Num_Nodes(1)+rd) = (INPUT_D(line)*LagPoly_Vals(0) + INPUT_D(line+1)*LagPoly_Vals(1))
            VX_Holder(re*Num_Nodes(1)+rd) = (INPUT_V(line)*LagPoly_Vals(0) + INPUT_V(line+1)*LagPoly_Vals(1))

            Density_Holder(re*NUM_NODES(1)+rd)  = Density
            Velocity_Holder(re*NUM_NODES(1)+rd) = Velocity
            r_Holder(re*NUM_NODES(1)+rd) = CUR_R_LOCS(rd)

            ! Calculate Usable Quantities
            Pressure = kappa * Density**gamma
            Energy = Pressure/(gamma - 1.0_idp)
            Specific_Enthalpy = C_Square + (Energy + Pressure)/Density

            
            vsqr = Velocity*Velocity
            LF_sqr = 1.0_idp/(1.0_idp- vsqr/C_Square)
            

            !  Calculate CFA Input Values
            E  = Density*Specific_Enthalpy*LF_sqr - Pressure
            Si = Density*Specific_Enthalpy*LF_sqr*Velocity/C_Square
            S  = Density*Specific_Enthalpy*LF_sqr*vsqr/C_Square + 3.0_idp* Pressure

!            PRINT*,re,rd,cur_r_locs(rd),E,S,Si
!            PRINT*,re,Velocity
!            PRINT*,Density/(gram/centimeter**3),Velocity,Specific_Enthalpy,LF_Sqr
!            PRINT*,Density/(gram/centimeter**3),Density

            
!            PRINT*,E,Specific_Enthalpy,LF_Sqr,Density/(gram/centimeter**3),Pressure

            DO pd = 1,NUM_NODES(3)
            DO td = 1,NUM_NODES(2)

                nd = Quad_Map(rd,td,pd)

                Input_E(nd, re,te,pe) = E
                Input_Si(nd, re, te, pe, 1) = Si
                Input_Si(nd, re, te, pe, 2) = 0.0_idp
                Input_Si(nd, re, te, pe, 3) = 0.0_idp
                Input_S(nd, re, te, pe) = S

            END DO ! td
            END DO ! pd
        END DO ! rd
    END DO ! re
END DO ! te
END DO ! pe

IF ( OUTPUT_PRIMATIVES_FLAG == 1 ) THEN

    CALL Output_Primatives( Density_Holder, Velocity_Holder, r_Holder, Num_Radial_Points)

END IF




DEALLOCATE( Density_Holder  )
DEALLOCATE( Velocity_Holder )
DEALLOCATE( DX_Holder       )
DEALLOCATE( VX_Holder       )


END SUBROUTINE CONVERT_SELF_SIMILAR_3D








!+201+###########################################################################!
!                                                                                !
!                  Calculate_Yahil_Density                                       !
!                                                                                !
!################################################################################!
SUBROUTINE Calculate_Yahil_Density( t, kappa, gamma,                        &
                                    Num_Nodes, NUM_LINES,                   &
                                    INPUT_R_QUAD,                           &
                                    NUM_R_ELEM, NUM_T_ELEM, NUM_P_ELEM,     &
                                    Delta_R, r_locs,                        &
                                    Input_D,                                &
                                    Input_X,                                &
                                    Output_D                                )


REAL(idp),                              INTENT(IN)                  ::  t, kappa, gamma
INTEGER,    DIMENSION(1:3),             INTENT(IN)                  ::  Num_Nodes
INTEGER,                                INTENT(IN)                  ::  Num_LINES
REAL(idp),  DIMENSION(1:NUM_NODES(1)),  INTENT(IN)                  ::  INPUT_R_QUAD

INTEGER,                                INTENT(IN)                  ::  NUM_R_ELEM
INTEGER,                                INTENT(IN)                  ::  NUM_T_ELEM
INTEGER,                                INTENT(IN)                  ::  NUM_P_ELEM

REAL(idp),  DIMENSION(NUM_R_ELEM),      INTENT(IN)                  ::  Delta_R
REAL(idp),  DIMENSION(0:NUM_R_ELEM),    INTENT(IN)                  ::  r_locs

REAL(idp),  DIMENSION(1:NUM_LINES),     INTENT(IN)                  ::  Input_D
REAL(idp),  DIMENSION(1:NUM_LINES),     INTENT(IN)                  ::  Input_X

REAL(idp),  DIMENSION(1:Num_Nodes(1)*Num_Nodes(2)*Num_Nodes(3),     &
                            0:NUM_R_ELEM-1,                         &
                            0:NUM_T_ELEM-1,                         &
                            0:NUM_P_ELEM-1    ),    INTENT(INOUT)   ::  Output_D


INTEGER                                                             ::  re
INTEGER                                                             ::  rd, td, pd

INTEGER                                                             ::  nd, line, line_min
REAL(idp)                                                           ::  xloc
REAL(idp), DIMENSION(0:1)                                           ::  xlocs
REAL(idp)                                                           ::  x

REAL(idp), DIMENSION(0:1)                                           ::  LagPoly_Vals
REAL(idp), DIMENSION(1:NUM_NODES(1))                                ::  CUR_R_LOCS

REAL(idp)                                                           ::  D_FACTOR
REAL(idp)                                                           ::  X_FACTOR

INTEGER                                                             ::  Num_Radial_Points

REAL(idp), DIMENSION(:), ALLOCATABLE                                ::  D_Holder
REAL(idp), DIMENSION(:), ALLOCATABLE                                ::  V_Holder
REAL(idp), DIMENSION(:), ALLOCATABLE                                ::  R_Holder



Num_Radial_Points = NUM_R_ELEM*NUM_NODES(1)

ALLOCATE( D_Holder(1:Num_Radial_Points) )
ALLOCATE( V_Holder(1:Num_Radial_Points) )
ALLOCATE( R_Holder(1:Num_Radial_Points) )

V_Holder = 0.0_idp


xlocs(0) = -1.0_idp
xlocs(1) = 1.0_idp

D_FACTOR = 1.0_idp/(Grav_Constant_G*t*t )

X_Factor = kappa**(-0.5_idp)                               &
        *(Grav_Constant_G**((gamma-1.0_idp)/2.0_idp))       &
        *((t)**(gamma-2.0_idp))



line_min = 1
DO re = 0,NUM_R_ELEM-1
    CUR_R_LOCS(:) = Delta_R(Re+1) * (INPUT_R_QUAD(:)+0.5_idp) + r_locs(re)



    DO rd = 1,NUM_NODES(1)
        xloc = CUR_R_LOCS(rd)*X_Factor
    !            CALL Find_Line_SUB(xloc, Input_X, NUM_LINES)
        line = Find_Line(xloc, Input_X, NUM_LINES)

        x = MAP_TO_X_SPACE(Input_X(Line),Input_X(Line+1),xloc)
        IF ( x > 1 ) THEN
            x = 1
        END IF
        LagPoly_Vals = Lagrange_Poly(x, 1, xlocs)


        D_Holder(re*NUM_NODES(1)+rd) = D_FACTOR                              &
                                     * ( INPUT_D(line)*LagPoly_Vals(0)       &
                                       + INPUT_D(line+1)*LagPoly_Vals(1)     )

        r_Holder(re*NUM_NODES(1)+rd) = CUR_R_LOCS(rd)

    END DO ! rd

END DO ! RE




DO re = 0,NUM_R_ELEM-1
DO pd = 1,NUM_NODES(3)
DO td = 1,NUM_NODES(2)
DO rd = 1,Num_Nodes(1)

    nd = Quad_Map(rd,td,pd)

    Output_D(nd, re,:,:) = D_Holder(re*Num_Nodes(1)+rd)

END DO ! rd
END DO ! td
END DO ! pd
END DO ! re





IF ( OUTPUT_PRIMATIVES_FLAG == 1 ) THEN

    CALL Output_Primatives( D_Holder, V_Holder, R_Holder, Num_Radial_Points)

END IF




DEALLOCATE( D_Holder )
DEALLOCATE( V_Holder )
DEALLOCATE( R_Holder )



END SUBROUTINE Calculate_Yahil_Density









!+201+###########################################################################!
!                                                                                !
!                  CREATE_SELFSIM_NEWT_SOL                                       !
!                                                                                !
!################################################################################!
SUBROUTINE CREATE_SELFSIM_NEWT_SOL( Line_Count, R_Values, Enclosed_Mass )

INTEGER, INTENT(IN)                                             :: Line_Count
REAL(idp), DIMENSION(1:Line_Count), INTENT(IN)           :: R_Values
REAL(idp), DIMENSION(1:Line_Count), INTENT(IN)           :: Enclosed_Mass

INTEGER     :: i



SELFSIM_R_VALS(0) = 0.0_idp
SELFSIM_R_VALS(1:NUM_ENTRIES) = R_Values(1:NUM_ENTRIES)

SELFSIM_POT_VALS(NUM_ENTRIES) = -GRAV_Constant_G*Enclosed_Mass(Num_Entries)/R_Values(Num_Entries)

DO i = NUM_ENTRIES-1,1,-1

    SELFSIM_POT_VALS(i) = SELFSIM_POT_VALS(i+1)                         &
                        - Grav_Constant_G*Enclosed_Mass(i)              &
                        * (SELFSIM_r_Vals(i+1)-SELFSIM_R_Vals(i))       &
                        / (SELFSIM_R_Vals(i)*SELFSIM_R_Vals(i))


END DO

SELFSIM_POT_VALS(0) = SELFSIM_POT_VALS(1)                           &
                    - 3*Grav_Constant_G*Enclosed_Mass(1)            &
                    /(2*SELFSIM_R_VALS(1))



END SUBROUTINE CREATE_SELFSIM_NEWT_SOL








!+201+###########################################################################!
!                                                                                !
!                         SELFSIM_NEWT_SOL                                       !
!                                                                                !
!################################################################################!
FUNCTION SELFSIM_NEWT_SOL( r, theta, phi )

REAL(idp), INTENT(IN)     :: r, theta, phi
REAL(idp)                 :: SELFSIM_NEWT_SOL

INTEGER                          :: cur_entry
INTEGER                          :: i

REAL(idp)                 :: deltar

DO i = 0,NUM_ENTRIES-1

    IF ( r == 0 ) THEN
        cur_entry = 0

    ELSE IF (( r > SELFSIM_R_VALS(i) ) .AND. ( r <= SELFSIM_R_VALS(i+1) ) ) THEN
        cur_entry = i

    ELSE IF ( r > SELFSIM_R_VALS(Num_Entries) ) THEN
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
!                                                        !
!                                                                                !
!################################################################################!
SUBROUTINE CREATE_SELFSIM_SHIFT_SOL( Num_Nodes,                                  &
                                     NUM_R_ELEM, NUM_T_ELEM, NUM_P_ELEM,         &
                                     Sr, r_locs                                  )


INTEGER, DIMENSION(1:3), INTENT( IN )                              ::  Num_Nodes
INTEGER,                 INTENT( IN )                              ::  NUM_R_ELEM,  &
                                                                       NUM_T_ELEM,  &
                                                                       NUM_P_ELEM

REAL(idp), DIMENSION( 1:NUM_Nodes(1)*Num_Nodes(2)*Num_Nodes(3),  &
                             0:NUM_R_ELEM-1,                            &
                             0:NUM_T_ELEM-1,                            &
                             0:NUM_P_ELEM-1,                            &
                             1:3 ), INTENT(IN)                    ::  Sr

REAL(idp), DIMENSION( 0:NUM_R_ELEM ),  INTENT(IN)          ::  r_locs


INTEGER                                                           ::  Ord

INTEGER                                                           ::  i, j, re, reb

REAL(idp), DIMENSION(:), ALLOCATABLE                       :: x_locs,     &
                                                                     ri_locs,    &
                                                                     wi

REAL(idp), DIMENSION(:,:), ALLOCATABLE                     :: rij_locs,    &
                                                                     PSI_10,     &
                                                                     Sr_New

REAL(idp), DIMENSION(:), ALLOCATABLE                       :: AlphaPsi
REAL(idp)                                                  :: Psi
REAL(idp)                                                  :: Outer_Int
REAL(idp)                                                  :: Inner_Int

Ord = 6



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
   ri_locs(:) = (r_locs(re+1)-r_locs(re))/2.0_idp* ( x_locs(:) + 1.0_idp) + r_locs(re)

   ! Calculate the Alpha Psi values at each of the Outer Integral's Quadrature Points !
   DO i = 1,Ord
      AlphaPsi(i) =  1.0_idp+ 0.5_idp*SELFSIM_NEWT_SOL(ri_locs(i),0.0_idp,0.0_idp)/C_Square
   END DO


   DO i = 1,Ord

      ! Calculate the Quadrature Points for each of the Inner Integrals
      rij_locs(:,i) = (ri_locs(i) - 0.0_idp)/2.0_idp*(x_locs(:) + 1.0_idp) + 0.0_idp


      ! Calculate Psi^10 values at each of the Inner Quadrature Points
      DO j = 1,Ord

           Psi = 1.0_idp- 0.5_idp*SELFSIM_NEWT_SOL(rij_locs(j,i),0.0_idp,0.0_idp)/C_Square
           Psi_10(j,i) = Psi**10

      END DO

      ! Calculate Sr values at each quadrature Point.
      DO j = 1,Ord
         DO reb = 0,NUM_R_ELEM - 1


            IF ( (rij_Locs(j,i) > r_locs(reb)) .AND. (rij_Locs(j,i) .LE. r_locs(reb+1)) ) THEN


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


      END DO ! j Loop

      Outer_Int = Outer_Int                                       &
                + AlphaPsi(i)                                     &
                / (ri_locs(i)*ri_locs(i)*ri_locs(i)*ri_locs(i))   &
                * Inner_Int                                       &
                * ( ri_locs(i) - 0.0_idp)/2.0_idp                &
                * wi(i)

   END DO ! i Loop



   IF ( re == 0 ) THEN

      SELFSIM_SHIFT_VALS(re+1) = (3.0_idp/2.0_idp)                      &
                               * 8.0_idp*pi* GR_Source_Scalar           &
                               * r_locs(re+1)                           &
                               * ( r_locs(re+1) - r_locs(re))/2.0_idp  &
                               * Outer_Int

   ELSE
        
      SELFSIM_SHIFT_VALS(re+1) = r_locs(re+1)/r_locs(re)                &
                               * SELFSIM_SHIFT_VALS(re)                 &
                               + (3.0_idp/2.0_idp)                      &
                               * 8.0_idp*pi*GR_Source_Scalar            &
                               * r_locs(re+1)                           &
                               * ( r_locs(re+1) - r_locs(re))/2.0_idp  &
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

REAL(idp), INTENT(IN)                              :: r
INTEGER, INTENT(IN)                                       :: NUM_R_ELEM

REAL(idp), DIMENSION(0:NUM_R_ELEM), INTENT(IN)     :: r_locs


REAL(idp)                 :: SELFSIM_SHIFT_SOL

INTEGER                          :: cur_entry
INTEGER                          :: i

REAL(idp)                 :: deltar

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

REAL(idp), INTENT(IN)     :: r
REAL(idp)                 :: SELFSIM_NEWT_SOL

INTEGER                          :: cur_entry
INTEGER                          :: i


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
REAL(idp), INTENT(IN)                        :: x
REAL(idp), INTENT(IN), DIMENSION(0:Ord)      :: xlocs


INTEGER                                             :: i,j
REAL(idp), DIMENSION(0:Ord)                  :: tmp
REAL(idp), DIMENSION(0:Ord)                  :: Lagrange_Poly



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

REAL(idp)                            ::  Map_To_X_Space
REAL(idp), intent(in)                ::  ra, rb
REAL(idp), intent(in)                ::  r

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
REAL(idp), INTENT(INOUT), DIMENSION(1:Ord)       ::  xloc, weights

INTEGER                                                 :: i, j, m
REAL(idp)                                        :: p1, p2, p3, pp, z, z1


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

    weights(i) = 2.0_idp/((1.0_idp- z*z)*pp*pp)
    weights(Ord-i+1) = weights(i)

END DO

END SUBROUTINE Initialize_LG_Quadrature






!+503+##################################################################!
!                                                                       !
!                                                                       !
!#######################################################################!
PURE INTEGER FUNCTION Find_Line( x, x_list, list_len )

REAL(idp), INTENT(IN)            :: x, x_list(list_len)
INTEGER, INTENT(IN)                     :: list_len

INTEGER                                 :: up, down, mid

up = list_len+1
down = 1
DO WHILE (up - down > 1)
    mid = (up + down)/2
    IF ( (x_list(list_len)>=x_list(1)).eqv.(x>=x_list(mid)) ) THEN
        down = mid
    ELSE
        up = mid
    END IF
END DO

IF ( x == x_list(1) ) THEN
    Find_Line = 1
ELSEIF ( x == x_list(list_len) ) THEN
    Find_Line = list_len - 1
ELSEIF ( down == List_len ) THEN
    Find_Line = List_len - 1
ELSE
    Find_Line = down
END IF

END FUNCTION Find_Line









!+101+###########################################################################!
!                                                                                !
!                  Initialize_Yahil_Sources                                           !
!                                                                                !
!################################################################################!
FUNCTION Calc_Yahil_Central_E( t_in, kappa, gamma )


REAL(idp),  INTENT(IN)                      ::  t_in
REAL(idp),  INTENT(IN)                      ::  kappa
REAL(idp),  INTENT(IN)                      ::  gamma

REAL(idp)                                   ::  Calc_Yahil_Central_E


REAL(idp)                                   ::  t


REAL(idp), DIMENSION(:),ALLOCATABLE         ::  Input_X,    &
                                                Input_D,    &
                                                Input_V,    &
                                                Input_M

REAL(idp)                                   ::  Density, Velocity

REAL(idp)                                   ::  Pressure, Energy
REAL(idp)                                   ::  vsqr, LF_sqr


REAL(idp)                                   ::  D_FACTOR,               &
                                                V_FACTOR

REAL(idp)                                   ::  Specific_Enthalpy

REAL(idp)                                   ::  Kappa_WUnits

CHARACTER(LEN=128)                          :: line

INTEGER                                     :: NUM_LINES
INTEGER                                     :: CUR_LINE


INTEGER                                     :: nread
INTEGER                                     :: istat





101 FORMAT (a128)


nread = 42
OPEN(UNIT=nread, FILE='../../Input/YahilHomologousCollapse_Gm_130.dat', STATUS='OLD', IOSTAT=istat)
IF ( istat .NE. 0 ) THEN
    PRINT*,"Could not open 'Input/YahilHomologousCollapse_Gm_130.dat'. "
END IF
REWIND(nread)
NUM_LINES = 0
DO
    READ(nread, *, IOSTAT=istat) line
    IF ( istat .NE. 0 ) EXIT
    NUM_LINES = NUM_LINES + 1
END DO
NUM_LINES = NUM_LINES -1



ALLOCATE( Input_X(1:NUM_LINES), Input_D(1:NUM_LINES), Input_V(1:NUM_LINES) )
ALLOCATE( Input_M(1:NUM_LINES) )

NUM_ENTRIES = NUM_LINES





REWIND(nread)
CUR_LINE = 1
READ(nread,*)
DO

    READ(nread, 101, IOSTAT=istat) line
    IF ( istat .NE. 0 ) EXIT


    READ(line,*) Input_X(CUR_LINE), Input_D(CUR_LINE), Input_V(CUR_LINE), Input_M(CUR_LINE)
    CUR_LINE = CUR_LINE + 1
    
END DO

CLOSE(UNIT=nread,STATUS='keep',IOSTAT=istat)


t = t_in*Millisecond

IF ( .TRUE. ) THEN
    Kappa_wUnits = Kappa*((Erg/Centimeter**3)/(Gram/Centimeter**3)**Gamma)
ELSE
    Kappa_wUnits = 18.394097187796024_idp
    PRINT*,"Kappa_wUnits over written.",Kappa_wUnits
END IF


D_FACTOR = 1.0_idp/(Grav_Constant_G*t*t )
V_FACTOR = SQRT(kappa)                          &
         * Grav_Constant_G**((1.0_idp-gamma)/2.0_idp)  &
         * t**(1.0_idp- gamma)


! Interpolate Self-Similar Values to Input locations
Density  = INPUT_D(1)*D_FACTOR
Velocity = INPUT_V(1)*V_FACTOR

! Calculate Usable Quantities
Pressure = kappa * Density**gamma
Energy = Pressure/(gamma - 1.0_idp)

Specific_Enthalpy = C_Square + (Energy + Pressure)/Density


vsqr = Velocity*Velocity
LF_sqr = 1.0_idp/(1.0_idp- vsqr/C_Square)
            

!  Calculate CFA Input Values
Calc_Yahil_Central_E  = Density*Specific_Enthalpy*LF_sqr - Pressure








END FUNCTION  Calc_Yahil_Central_E













END MODULE External_Yahil_Profile_Module

