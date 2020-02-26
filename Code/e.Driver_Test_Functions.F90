    !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Driver_Test_Functions_Module                                                 !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!   Contains the functions and subroutines used in the test codes provided.      !##!
!##!    These functions initialize specific sources that produce known potentials   !##!
!##!    which allow for comparison of the code's solution to the problem and the    !##!
!##!    known analytic solution.                                                    !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   Poseidon_Initialize_Test_Problem                                    !##!
!##!                                                                                !##!
!##!    +201+   Test_Source_Spherical_Symmetry_No_Surface                           !##!
!##!                                                                                !##!
!##!    +301+   Test_Spherical_Symmetry_No_Surface                                  !##!
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
                ONLY :  idp, pi, Speed_of_Light,Grav_Constant_G


USE DRIVER_Parameters, &
                ONLY :  DRIVER_R_ELEMS,                        &
                        DRIVER_T_ELEMS,                        &
                        DRIVER_P_ELEMS,                        &
                        DRIVER_R_LOCS,                         &
                        DRIVER_R_INPUT_NODES,                  &
                        DRIVER_T_INPUT_NODES,                  &
                        DRIVER_P_INPUT_NODES,                  &
                        DRIVER_INNER_RADIUS,                   &
                        DRIVER_OUTER_RADIUS,                   &
                        DRIVER_POTENTIAL,                      &
                        DRIVER_SHIFT_VAL,                      &
                        DRIVER_E,                              &
                        DRIVER_S,                              &
                        DRIVER_Si,                             &
                        myID,                                   &
                        myID_Theta,                             &
                        myID_Phi,                               &
                        POWER_A,                                &
                        RHO_O,                                  &
                        Analytic_Solution,                      &
                        Shift_Solution,                         &
                        SELFSIM_T,                              &
                        SELFSIM_KAPPA,                          &
                        SELFSIM_GAMMA,                          &
                        SELFSIM_ECC

USE CHIMERA_TEST_FUNCS_Module, &
                ONLY :  Test_Chimera_Simulated_Potential,       &
                        Test_Chimera_Simulated_Shift


USE Driver_Additional_Functions_Module, &
                ONLY :  Map_From_X_Space,                       &
                        Initialize_LG_Quadrature_Locations,     &
                        Generate_Defined_Mesh

USE SelfSimilar_Module, &
            ONLY :  UNPACK_SELF_SIMILAR,                        &
                    SELFSIM_NEWT_SOL,                           &
                    SELFSIM_SHIFT_SOL
 
USE mpi




IMPLICIT NONE











!*F&S*==========================================!
!                                               !
!           Functions & Subroutines             !
!                                               !
!===============================================!
CONTAINS

 !+102+############################################################################!
!                                                                                   !
!       Poseidon_Initialize_Test_Problem                                            !
!                                                                                   !
!      This code initializes the source for the selected problem as well as points  !
!       the code to the analytic solution associated with the given source.         !
!                                                                                   !
!-----------------------------------------------------------------------------------!
!                                                                                   !
!   Input                                                                           !
!                                                                                   !
!           Test_Number     -   Integer, Selects problem to be solved.              !
!                                                                                   !
!                                   1)  Spherically Symmetric Source with radial    !
!                                           dependence.                             !
!                                                                                   !
!                                   2) Spherically Symmetric source with radial     !
!                                       dependence and include surface discontinuity!
!                                                                                   !
!                                   3) MacLaurin Ellipsoid                          !
!                                                                                   !
 !#################################################################################!
SUBROUTINE Poseidon_Initialize_CFA_Test_Problem_CHIMERA(    Test_Number,                                &
                                                            Num_Global_R, Num_Global_T, Num_Global_P,   &
                                                            Num_Local_R, Num_Local_T, Num_Local_P,      &
                                                            NUM_LOCS, R_QUAD, T_QUAD,                   &
                                                            deltar, r_locs, t_locs,                     &
                                                            Input_E, Input_S, Input_Si )




INTEGER,                        INTENT(IN)                                  ::  Test_Number

INTEGER,                        INTENT(IN)                                  ::  Num_Global_R,    &
                                                                                Num_Global_T,    &
                                                                                Num_Global_P

INTEGER,                        INTENT(IN)                                  ::  Num_Local_R,    &
                                                                                Num_Local_T,    &
                                                                                Num_Local_P

INTEGER, DIMENSION(1:3),        INTENT(IN)                                  ::  NUM_LOCS
 
REAL(KIND = idp), DIMENSION(1:NUM_LOCS(1)), INTENT(IN)                      ::  R_QUAD
REAL(KIND = idp), DIMENSION(1:NUM_LOCS(2)), INTENT(IN)                      ::  T_QUAD


REAL(KIND = idp), DIMENSION(1:NUM_GLOBAL_R), INTENT(INOUT)                     ::  deltar
REAL(KIND = idp), DIMENSION(0:NUM_GLOBAL_R), INTENT(INOUT)                     ::  r_locs
REAL(KIND = idp), DIMENSION(0:NUM_GLOBAL_T), INTENT(INOUT)                     ::  t_locs


REAL(KIND = idp), DIMENSION(    1:NUM_LOCS(1)*NUM_LOCS(2)*NUM_LOCS(3), &
                                0:NUM_Local_R-1,                      &
                                0:NUM_Local_T-1,                      &
                                0:NUM_Local_P-1 ), INTENT(INOUT)              ::  Input_E


REAL(KIND = idp), DIMENSION(    1:NUM_LOCS(1)*NUM_LOCS(2)*NUM_LOCS(3), &
                                0:NUM_Local_R-1,                      &
                                0:NUM_Local_T-1,                      &
                                0:NUM_Local_P-1    ), INTENT(INOUT)           ::  Input_S


REAL(KIND = idp), DIMENSION(    1:NUM_LOCS(1)*NUM_LOCS(2)*NUM_LOCS(3), &
                                0:NUM_Local_R-1,                      &
                                0:NUM_Local_T-1,                      &
                                0:NUM_Local_P-1,                      &
                                1:3),                 INTENT(INOUT)         ::  Input_Si



INTEGER, DIMENSION(1:3)                                                     :: QUAD_NUMBERS




INTEGER                                                                     ::  i, j, k,    &
                                                                                Num_DOF,    &
                                                                                nPROC,      &
                                                                                ierr,       &
                                                                                my_i_start, &
                                                                                my_j_start



INTEGER                                                                     ::  re, te, pe, &
                                                                                rd, td, pd



REAL(KIND = idp), DIMENSION(1:NUM_LOCS(1)*NUM_LOCS(2)*NUM_LOCS(3),     &
                            1:Num_Global_R,                         &
                            1:Num_Global_T,                         &
                            1:Num_Global_P)                                 ::  Test_Source_Input




REAL(KIND = idp), DIMENSION(1:NUM_LOCS(1)*NUM_LOCS(2)*NUM_LOCS(3),  &
                            0:NUM_Local_R-1, 0:NUM_Local_T-1, 0:NUM_Local_P-1    )       ::  rho




INTEGER                                                                                 :: Output_Here





REAL(KIND = idp)                                            ::  rho_c, rho_here, gamma,     &
                                                                h, epsilon, Pressure,           &
                                                                r_here, t_here, p_here,         &
                                                                R_Surface, Velocity, Lorentz,   &
                                                                D, Polytropic_K, Vel_Shift



REAL(KIND = idp)    :: csqr




csqr = Speed_of_Light*Speed_of_Light


QUAD_NUMBERS(1) = NUM_LOCS(1)
QUAD_NUMBERS(2) = NUM_LOCS(2)
QUAD_NUMBERS(3) = NUM_LOCS(3)


IF (Test_Number .EQ. 1) THEN


    CALL Test_Source_Spherical_Symmetry(QUAD_NUMBERS, Test_Source_Input,  &
                                        Num_Global_R, Num_Global_T, Num_Global_P)

    my_i_start = myid_theta*Num_Local_T
    my_j_start = myid_phi*Num_Local_P

    DO pe = 0, Num_Local_P-1
        DO te = 0,Num_Local_T-1
            DO re = 0,NUM_Local_R-1

                rho(:, re, te, pe) = Test_Source_Input(:,re+1,              &
                                                        my_i_start + te+1,  &
                                                        my_j_start + pe+1   )
            END DO
        END DO
    END DO

    Input_E = (Grav_Constant_G/csqr) * rho
    Input_S = 0.0_idp
    Input_Si = 0.0_idp


    Analytic_Solution => Test_Spherical_Symmetry_No_Surface







ELSE IF (Test_Number == 2) THEN




    Analytic_Solution => Test_Chimera_Simulated_Potential
    Shift_Solution => Test_Chimera_Simulated_Shift

    Input_E = DRIVER_E
    Input_S = DRIVER_S
    Input_Si = DRIVER_Si




ELSE IF ( Test_Number == 3) THEN


    CALL UNPACK_SELF_SIMILAR( SELFSIM_T, SELFSIM_KAPPA, SELFSIM_GAMMA, SELFSIM_ECC,  &
                              Num_Locs, R_Quad, T_QUAD,                 &
                              NUM_LOCAL_R, NUM_LOCAL_T, NUM_LOCAL_P,    &
                              deltar, r_locs, t_locs,                   &
                              Input_E, Input_S, Input_Si                )




    Analytic_Solution => SELFSIM_NEWT_SOL
    Shift_Solution => SELFSIM_SHIFT_SOL


END IF






END SUBROUTINE Poseidon_Initialize_CFA_Test_Problem_CHIMERA











!+201+######################
!
!   Spherically Symmetric
!
!###########################
SUBROUTINE Test_Source_Spherical_Symmetry(Num_Nodes, Test_Source_Input, R_Elements, T_Elements, P_Elements)



INTEGER,                        INTENT(IN)                                  ::  R_Elements,     &
                                                                                T_Elements,     &
                                                                                P_Elements


INTEGER,    DIMENSION(1:3),     INTENT(IN)                                  ::  Num_Nodes


REAL(KIND = idp), DIMENSION(1:Num_Nodes(1)*Num_Nodes(2)*Num_Nodes(3),   &
                            1:R_ELEMENTS,                           &
                            1:T_ELEMENTS,                           &
                            1:P_ELEMENTS),                          &
                            INTENT(INOUT)                                   ::  Test_Source_Input

                         !                        !
                        !!  Subroutine Variables  !!
                         !                        !
INTEGER                                                     ::  re, te, pe, i,                  &
                                                                Output_R, Output_T, Output_P,   &
                                                                l, m, d,                        &
                                                                Output_Here,                    &
                                                                Num_Output_DOF


REAL(KIND = idp), DIMENSION(1:Num_Nodes(1))                 ::  Output_R_X_Locations


REAL(KIND = idp), DIMENSION(1:Num_Nodes(1))                 ::  Output_R_Locations

REAL(KIND = idp)                                            ::  TMP

REAL(KIND = idp)                                            ::  Star_Surface



                             !                    !
                            !!  Global Variables  !!
                             !                    !

!INTEGER,                                                   ::  DRIVER_R_ELEMS,             &
!                                                               DRIVER_T_ELEMS,             &
!                                                               DRIVER_P_ELEMS


REAL(KIND = idp), DIMENSION(1:R_Elements) :: Delta_R_Vector
REAL(KIND = idp), DIMENSION(0:R_Elements) :: rlocations

rlocations = 0.0_idp

Delta_R_Vector = (DRIVER_OUTER_RADIUS - 0.0_idp) / REAL(R_Elements )
Star_Surface = DRIVER_OUTER_RADIUS + 1.0_idp


CALL Generate_Defined_Mesh(R_Elements, 0.0_idp, Delta_R_Vector, rlocations)



                          !                                                 !
                         !!                                                 !!
                        !!!               Calculate Local DOF               !!!
                         !!                                                 !!
                          !                                                 !

Num_Output_DOF = Num_Nodes(1) * Num_Nodes(2) * Num_Nodes(3)



                          !                                                 !
                         !!                                                 !!
                        !!!     Map Output Locations to [ -1, 1 ] Space      !!!
                         !!                                                 !!
                          !                                                 !

Output_R_X_Locations = Initialize_LG_Quadrature_Locations(Num_Nodes(1))




DO re = 0,R_ELEMENTS-1


     !                                                          !
    !!   Map Radial Locations from [-1,1] space to real space.   !!
     !                                                          !


    Output_R_Locations = Map_From_X_Space(rlocations(re), rlocations(re + 1), Output_R_X_Locations)


    DO te = 0,T_ELEMENTS-1



        DO pe = 0,P_ELEMENTS-1


             !                                          !
            !!   Set/Reset Output Vector Location       !!
             !                                          !
            Output_Here = 1





            DO Output_P = 1,Num_Nodes(3)

                DO Output_T = 1,Num_Nodes(2)

                    DO Output_R = 1,Num_Nodes(1)


                        TMP = 1.0_idp







                        IF (Output_R_Locations(Output_R) .LE. Star_Surface) THEN


                            IF (POWER_A >= 0) THEN

                                DO i = 1,POWER_A

                                    TMP = TMP*Output_R_Locations(Output_R)

                                END DO
                            END IF


                            IF (POWER_A < 0) THEN
                                DO i = 1,ABS(POWER_A)

                                    TMP = TMP/Output_R_Locations(Output_R)

                                END DO
                            END IF

                            Test_Source_Input(Output_Here, re+1, te+1, pe+1) = RHO_O*TMP


                        ELSE

                            Test_Source_Input(Output_Here, re+1, te+1, pe+1) = 0.0_idp


                        END IF



                         !                                      !
                        !!   Increment output vector location   !!
                         !                                      !
                        Output_Here = Output_Here + 1





                    END DO  !   Output_R Loop

                END DO  !   Output_T Loop

            END DO  !   Output_P Loop







        END DO  ! pe Loop

    END DO  ! te Loop

END DO  !   re Loop




!PRINT*,"In TEST_SOURCE : ", Test_Source_Input


END SUBROUTINE Test_Source_Spherical_Symmetry
































!+301+##################################################################
!
!   Test_Spherical_Symmetry_No_Surface
!
!#######################################################################
PURE FUNCTION Test_Spherical_Symmetry_No_Surface(r, theta, phi)


REAL(KIND = idp),INTENT(IN)         :: r, theta, phi
REAL(KIND = idp)                    :: Test_Spherical_Symmetry_No_Surface

INTEGER                             :: i
REAL(KIND = idp)                    :: TMP


REAL(KIND = idp)                    ::  r_out_sqr,          &
                                        r_out_cube,         &
                                        enclosed_mass




r_out_sqr = DRIVER_OUTER_RADIUS*DRIVER_OUTER_RADIUS
r_out_cube = r_out_sqr * DRIVER_OUTER_RADIUS

enclosed_mass = (4.0_idp/(3.0_idp+POWER_A))* pi *RHO_O * DRIVER_OUTER_RADIUS**(3+POWER_A)

TMP = 1.0_idp
IF (POWER_A > -2) THEN



!    Test_Spherical_Symmetry_No_Surface = - Grav_Constant_G * enclosed_mass *(3*r_out_sqr - r*r)/(2*r_out_cube)
     Test_Spherical_Symmetry_No_Surface = - Grav_Constant_G * 4.0_idp*pi*RHO_O                     &
                                        * ( DRIVER_OUTER_RADIUS**(2+POWER_A)/(2+POWER_A)          &
                                          - r**(2+POWER_A)/(6+5*POWER_A+POWER_A*POWER_A)   )





ELSE IF (POWER_A  .EQ. -3) THEN



    Test_Spherical_Symmetry_No_Surface = Grav_Constant_G * 4.0_idp * pi *RHO_O*(log(r)+1)/r



ELSE IF (POWER_A  < -3) THEN


     Test_Spherical_Symmetry_No_Surface = - Grav_Constant_G * 4.0_idp*pi*RHO_O                     &
                                        * ( DRIVER_OUTER_RADIUS**(2+POWER_A)/(2+POWER_A)          &
                                          - r**(2+POWER_A)/(6+5*POWER_A+POWER_A*POWER_A)   )

END IF





END FUNCTION Test_Spherical_Symmetry_No_Surface









END MODULE Driver_Test_Functions_Module
