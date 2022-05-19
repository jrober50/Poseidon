   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Driver_Init_Data_On_Level_Module                                    !##!
!##!                                                                         !##!
!##!_________________________________________________________________________!##!
!##!                                                                         !##!
!##!                                                                         !##!
!##!=========================================================================!##!
!##!                                                                         !##!
!##!    Contains:                                                            !##!
!##!                                                                         !##!
!##!                                                                         !##!
!###############################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !##########################################################################!


!*D*================================!
!                                   !
!           Dependencies            !
!                                   !
!===================================!
USE ISO_C_BINDING

#ifdef POSEIDON_AMREX_FLAG
use amrex_base_module

USE amrex_box_module,       ONLY: &
  amrex_box
USE amrex_boxarray_module,  ONLY: &
  amrex_boxarray,         &
  amrex_boxarray_build,   &
  amrex_boxarray_destroy
USE amrex_distromap_module, ONLY: &
  amrex_distromap,       &
  amrex_distromap_build, &
  amrex_distromap_destroy
USE amrex_multifab_module,  ONLY: &
  amrex_multifab, &
  amrex_multifab_build
#endif


USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag

USE Poseidon_Units_Module, &
            ONLY :  Grav_Constant_G,    &
                    Speed_of_Light,     &
                    C_Square,           &
                    GR_Source_Scalar,   &
                    Centimeter,         &
                    Second,             &
                    Millisecond,         &
                    Erg,                &
                    Gram


USE Variables_Source, &
            ONLY :  Block_Source_E,             &
                    Block_Source_S,             &
                    Block_Source_Si

USE Variables_Quadrature, &
            ONLY :  INT_R_LOCATIONS,            &
                    NUM_R_QUAD_POINTS,          &
                    NUM_T_QUAD_POINTS,          &
                    NUM_P_QUAD_POINTS


USE Variables_External, &
            ONLY :  SelfSim_T,              &
                    SelfSim_Kappa,          &
                    SelfSim_Gamma

USE Poseidon_AMReX_Multilayer_Utilities_Module, &
            ONLY :  Find_Coarsest_Parent

USE Variables_Tables, &
            ONLY :  Level_dx

USE Functions_Math, &
            ONLY :  Lagrange_Poly

USE Timer_Routines_Module, &
            ONLY :  TimerStart,                     &
                    TimerStop

USE Timer_VAriables_Module, &
            ONLY :  Timer_Core_Init_Test_Problem

USE Variables_Functions, &
            ONLY :  Potential_Solution

USE Maps_Quadrature, &
            ONLY :  Quad_Map


USE Variables_Interface, &
            ONLY :  Caller_Set,                     &
                    Caller_nLevels,                 &
                    Caller_NQ,                      &
                    Caller_Quad_DOF,                &
                    Caller_xL,                      &
                    Caller_RQ_xlocs,                &
                    Caller_TQ_xlocs,                &
                    Caller_PQ_xlocs

IMPLICIT NONE


CONTAINS




 !+101+################################################################!
!                                                                       !
!          Poseidon_Init_Data_On_Level                                  !
!                                                                       !
 !#####################################################################!
SUBROUTINE Poseidon_Init_Data_On_Level(   Level,          &
                                            BLo, BHi,       &
                                            SLo, SHi,       &
                                            nComps,         &
                                            Src             )

INTEGER,                INTENT(IN)          ::  Level
INTEGER,                INTENT(IN)          ::  BLo(3), BHi(3)
INTEGER,                INTENT(IN)          ::  SLo(3), SHi(3)
INTEGER,                INTENT(IN)          ::  nComps
REAL(idp),              INTENT(INOUT)       ::  Src(SLo(1):SHi(1),  &
                                                    SLo(2):SHi(2),  &
                                                    SLo(3):SHi(3),  &
                                                    nComps          )

REAL(idp)                                   ::  M_Factor
REAL(idp)                                   ::  R_Factor
REAL(idp)                                   ::  V_Factor
REAL(idp)                                   ::  X_Factor
REAL(idp)                                   ::  D_Factor

REAL(idp), DIMENSION(:),ALLOCATABLE         ::  Input_X,    &
                                                Input_D,    &
                                                Input_V,    &
                                                Input_M,    &
                                                Input_R


CHARACTER(LEN=128)                          :: line_chr

INTEGER                                     :: Num_Lines
INTEGER                                     :: Cur_Line
INTEGER                                     :: Line_Min
INTEGER                                     :: nread
INTEGER                                     :: istat


INTEGER                                     ::  re, te, pe
INTEGER                                     ::  rd, td, pd
INTEGER                                     ::  Num_DOF
INTEGER                                     ::  Here
INTEGER                                     ::  i

REAL(idp)                                   ::  t
REAL(idp)                                   ::  Kappa_WUnits
REAL(idp)                                   ::  E_Units

REAL(idp)                                   ::  x
REAL(idp)                                   ::  xwidth
REAL(idp)                                   ::  xloc
REAL(idp), DIMENSION(0:1)                   ::  xlocs
REAL(idp), DIMENSION(1:Caller_NQ(1))        ::  cur_r_locs
REAL(idp)                                   ::  DROT
REAL(idp), DIMENSION(0:1)                   ::  LagPoly_Vals

REAL(idp)                                   ::  velocity
REAL(idp)                                   ::  vsqr
REAL(idp)                                   ::  LF_Sqr

REAL(idp)                                   ::  Psi_Holder

REAL(idp)                                   ::  Density
REAL(idp)                                   ::  Pressure
REAL(idp)                                   ::  Potential

REAL(idp), DIMENSION(:), ALLOCATABLE        ::  Enclosed_Mass
REAL(idp), DIMENSION(:), ALLOCATABLE        ::  Newtonian_Potential
REAL(idp)                                   ::  Energy
REAL(idp)                                   ::  Specific_Enthalpy

REAL(idp)                                   ::  S
REAL(idp)                                   ::  Si
REAL(idp)                                   ::  E

101 FORMAT (a128)


CALL TimerStart( Timer_Core_Init_Test_Problem )


nread = 42
OPEN(UNIT=nread, FILE='../../Input/YahilHomologousCollapse_Gm_130.dat', STATUS='OLD', IOSTAT=istat)
IF ( istat .NE. 0 ) THEN
    PRINT*,"Could not open 'Input/YahilHomologousCollapse_Gm_130.dat'. "
END IF
REWIND(nread)
NUM_LINES = 0
DO
    READ(nread, *, IOSTAT=istat) line_chr
    IF ( istat .NE. 0 ) EXIT
    NUM_LINES = NUM_LINES + 1
END DO
NUM_LINES = NUM_LINES -1



ALLOCATE( Input_X(1:NUM_LINES), Input_D(1:NUM_LINES), Input_V(1:NUM_LINES) )
ALLOCATE( Input_R(1:NUM_LINES) )
ALLOCATE( Input_M(1:NUM_LINES) )
ALLOCATE( Enclosed_Mass(1:NUM_LINES)  )
ALLOCATE( Newtonian_Potential(0:NUM_LINES)  )






REWIND(nread)
CUR_LINE = 1
READ(nread,*)
DO

    READ(nread, 101, IOSTAT=istat) line_chr
    IF ( istat .NE. 0 ) EXIT


    READ(line_chr,*) Input_X(CUR_LINE), Input_D(CUR_LINE), Input_V(CUR_LINE), Input_M(CUR_LINE)
    CUR_LINE = CUR_LINE + 1
    
END DO

CLOSE(UNIT=nread,STATUS='keep',IOSTAT=istat)




IF ( .TRUE. ) THEN
    Kappa_wUnits = SelfSim_Kappa*((Erg/Centimeter**3)/(Gram/Centimeter**3)**SelfSim_Gamma)
ELSE
    Kappa_wUnits = 18.394097187796024_idp
    PRINT*,"Kappa_wUnits over written.",Kappa_wUnits
END IF



t = SelfSim_T*Millisecond

R_Factor = SQRT(Kappa_wUnits)                                &
        *(Grav_Constant_G**((1.0_idp-SelfSim_Gamma)/2.0_idp))       &
        *((t)**(2.0_idp-SelfSim_Gamma))


Input_R = R_Factor*Input_X


!Input_V = 0.0_idp



E_Units = Erg/Centimeter**3




D_Factor = 1.0_idp/(Grav_Constant_G*t*t )

V_Factor = SQRT(kappa_WUnits)                                   &
         * Grav_Constant_G**((1.0_idp-SelfSim_Gamma)/2.0_idp)   &
         * t**(1.0_idp- SelfSim_Gamma)

X_Factor = kappa_WUnits**(-0.5_idp)                                 &
        *(Grav_Constant_G**((SelfSim_Gamma-1.0_idp)/2.0_idp))       &
        *((t)**(SelfSim_Gamma-2.0_idp))

M_Factor = kappa_WUnits**(3.0_idp/2.0_idp)                              &
         * Grav_Constant_G**((1.0_idp-3.0_idp*SelfSim_Gamma)/2.0_idp)   &
         * t**(4.0_idp- 3.0_idp* SelfSim_Gamma )


!PRINT*,D_Factor,V_Factor,X_Factor,M_Factor

Enclosed_Mass = M_Factor * Input_M
Newtonian_Potential(Num_Lines) = -GRAV_Constant_G*Enclosed_Mass(Num_Lines)/Input_R(Num_Lines)
DO i = Num_Lines-1,1,-1

     Newtonian_Potential(i) = Newtonian_Potential(i+1)              &
                            - Grav_Constant_G*Enclosed_Mass(i)      &
                            * (Input_R(i+1)-Input_R(i))             &
                            / (Input_R(i)*Input_R(i))

END DO
Newtonian_Potential(0) = Newtonian_Potential(1)                 &
                       - 3*Grav_Constant_G*Enclosed_Mass(1)     &
                       /(2*Input_R(1))




Num_DOF = nComps/5


!PRINT*,Int_R_Locations/2.0_idp



xlocs(0) = -1.0_idp
xlocs(1) = +1.0_idp
xwidth  = Caller_xL(2)-Caller_xL(1)

DROT = Level_dx(Level,1)/xwidth


DO pe = BLo(3),BHi(3)
DO te = BLo(2),BHi(2)

line_min = 1
DO re = BLo(1),BHi(1)

    Cur_R_Locs(:) = DROT * (Caller_RQ_xlocs(:) + 1.0_idp + re*xwidth)
    
!    PRINT*,re,DROT*2.0_idp
!    PRINT*,Cur_R_Locs


    DO rd = 1,Num_R_Quad_Points

        xloc = Cur_R_Locs(rd)*X_Factor
        cur_line = Find_Line(xloc, Input_X, NUM_LINES)

        x = MAP_TO_X_SPACE(Input_X(Cur_Line),Input_X(Cur_Line+1),xloc)
        IF ( x > 1 ) THEN
            x = 1
        END IF
        LagPoly_Vals = Lagrange_Poly(x, 1, xlocs)


        ! Interpolate Self-Similar Values to Input locations
        Density  = (INPUT_D(Cur_Line)*LagPoly_Vals(0) + INPUT_D(Cur_Line+1)*LagPoly_Vals(1))*D_FACTOR
        Velocity = (INPUT_V(Cur_Line)*LagPoly_Vals(0) + INPUT_V(Cur_Line+1)*LagPoly_Vals(1))*V_FACTOR
        Potential= Newtonian_Potential(Cur_Line)  *LagPoly_Vals(0)                                      &
                 + Newtonian_Potential(Cur_Line+1)*LagPoly_Vals(1)

        ! Calculate Usable Quantities
        Pressure = kappa_WUnits * Density**SelfSim_Gamma
        Energy = Pressure/(SelfSim_Gamma - 1.0_idp)
        Specific_Enthalpy = C_Square + (Energy + Pressure)/Density

        
        vsqr = Velocity*Velocity
        LF_sqr = 1.0_idp/(1.0_idp- vsqr/C_Square)
        

        !  Calculate CFA Input Values
        E  = Density*Specific_Enthalpy*LF_sqr - Pressure
        S  = Density*Specific_Enthalpy*LF_sqr*vsqr/C_Square + 3.0_idp* Pressure
        Si = Density*Specific_Enthalpy*LF_sqr*Velocity/C_Square
        

        Psi_Holder = (1.0_idp - 0.5_idp*Potential/C_Square)**6

        E  = E*Psi_Holder
        S  = S*Psi_Holder
        Si = Si*Psi_Holder
        


!        PRINT*,Level,re,rd,Cur_R_Locs(rd),E,S,Si

        DO pd = 1,Num_P_Quad_Points
        DO td = 1,Num_T_Quad_Points


            here = Quad_Map(rd,td,pd)

!            PRINT*,re,te,pe,rd,td,pd,here,2*Num_DOF+Here,Si

            Src(re,te,pe,0*Num_DOF+Here) = E
            Src(re,te,pe,1*Num_DOF+Here) = S
            Src(re,te,pe,2*Num_DOF+Here) = Si
            Src(re,te,pe,3*Num_DOF+Here) = 0.0_idp
            Src(re,te,pe,4*Num_DOF+Here) = 0.0_idp




        END DO ! td
        END DO ! pd
    END DO ! rd

END DO ! re
END DO ! te
END DO ! pe



CALL TimerStop( Timer_Core_Init_Test_Problem )


END SUBROUTINE Poseidon_Init_Data_On_Level











 





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








END MODULE Driver_Init_Data_On_Level_Module
