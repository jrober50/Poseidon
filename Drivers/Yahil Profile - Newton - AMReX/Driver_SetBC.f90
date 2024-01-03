   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Driver_SetBC_Module                                                   !##!
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
USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Poseidon_Parameters, &
            ONLY :  Verbose_Flag

USE Poseidon_Message_Routines_Module, &
            ONLY :  Driver_Init_Message


USE Poseidon_Units_Module, &
            ONLY :  Grav_Constant_G,    &
                    GravPot_Units,      &
                    Speed_of_Light,     &
                    C_Square,           &
                    GR_Source_Scalar,   &
                    Centimeter,         &
                    Second,             &
                    Millisecond,         &
                    Erg,                &
                    Gram

USE Variables_Functions, &
            ONLY :  Potential_Solution

USE Variables_External, &
            ONLY :  SelfSim_T,              &
                    SelfSim_Kappa,          &
                    SelfSim_Gamma,          &
                    NUM_ENTRIES,        &
                    SELFSIM_R_VALS,     &
                    SELFSIM_POT_VALS

USE Allocation_Yahil_Profile, &
            ONLY :  Allocate_Yahil_Profile,               &
                    Deallocate_Yahil_Profile

USE Variables_Mesh, &
            ONLY :  R_Outer

USE Poseidon_Interface_Boundary_Conditions, &
            ONLY :  Poseidon_Set_Uniform_Boundary_Conditions


USE External_Yahil_Profile_Module, &
            ONLY :  SELFSIM_NEWT_SOL,                       &
                    Create_Yahil_Newtonian_Solution_Coeffs

IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!   Driver_SetBC                                                                !
!                                                                               !
!###############################################################################!
SUBROUTINE Driver_SetBC( )


CHARACTER(LEN=1)                                    ::  Inner_BC_Types
CHARACTER(LEN=1)                                    ::  Outer_BC_Types
REAL(idp)                                           ::  Inner_BC_Values
REAL(idp)                                           ::  Outer_BC_Values


REAL(idp)                                           ::  R_Factor
REAL(idp)                                           ::  Kappa_wUnits


REAL(idp), DIMENSION(:), ALLOCATABLE                ::  Enclosed_Mass
REAL(idp), DIMENSION(:),ALLOCATABLE                 ::  Input_R
REAL(idp), DIMENSION(:),ALLOCATABLE                 ::  Input_X,    &
                                                        Input_D,    &
                                                        Input_V,    &
                                                        Input_M
REAL(idp)                                           ::  t

CHARACTER(LEN=128)                                  :: line

INTEGER                                             :: NUM_LINES
INTEGER                                             :: CUR_LINE


INTEGER                                             :: nread
INTEGER                                             :: istat



IF ( Verbose_Flag ) CALL Driver_Init_Message('Calculating boundary conditions.')



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




Kappa_wUnits = SelfSim_Kappa*((Erg/Centimeter**3)/(Gram/Centimeter**3)**SelfSim_Gamma)
t = SelfSim_T*Millisecond

R_Factor = SQRT(Kappa_wUnits)                                &
            *(Grav_Constant_G**((1.0_idp-SelfSim_Gamma)/2.0_idp))       &
            *((t)**(2.0_idp-SelfSim_Gamma))

Input_R = R_Factor*Input_X


Enclosed_Mass = Kappa_wUnits**(1.50_idp)                                   &
              * Grav_Constant_G**((1.0_idp-3.0_idp*SelfSim_Gamma)/2.0_idp)  &
              * (t**(4.0_idp- 3.0_idp*SelfSim_Gamma))                      &
              * Input_M

CALL Create_Yahil_Newtonian_Solution_Coeffs( NUM_LINES, Input_R, Enclosed_Mass )
Potential_Solution => SELFSIM_NEWT_SOL




Inner_BC_Types = "N"
Outer_BC_Types = "D"


Inner_BC_Values = 0.0_idp
Outer_BC_Values = Potential_Solution(R_Outer, 0.0_idp, 0.0_idp)

IF ( Verbose_Flag ) CALL Driver_Init_Message('Setting boundary conditions.')



CALL Poseidon_Set_Uniform_Boundary_Conditions("I", Inner_BC_Types, Inner_BC_Values)
CALL Poseidon_Set_Uniform_Boundary_Conditions("O", Outer_BC_Types, Outer_BC_Values)



END SUBROUTINE Driver_SetBC



END MODULE Driver_SetBC_Module

