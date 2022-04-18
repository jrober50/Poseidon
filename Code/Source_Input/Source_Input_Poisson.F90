   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE  Source_Intput_Poisson                                                !##!
!##!                                                                         !##!
!##!_________________________________________________________________________!##!
!##!                                                                         !##!
!##!    +i01+   Poseidon_Newtonian_Source_Input                              !##!
!##!                                                                         !##!
!##!    +201+   Poseidon_Input_Sources_Poisson_1D                           !##!
!##!    +202+   Poseidon_Input_Sources_Poisson_2D                           !##!
!##!    +203+   Poseidon_Input_Sources_Poisson_3D                           !##!
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

USE Variables_Poisson, &
            ONLY :  Source_Terms

USE Variables_Mesh, &
            ONLY :  NUM_R_ELEMENTS,     &
                    NUM_T_ELEMENTS,     &
                    NUM_P_ELEMENTS

USE Variables_Quadrature, &
            ONLY :  NUM_R_QUAD_POINTS,          &
                    NUM_T_QUAD_POINTS,          &
                    NUM_P_QUAD_POINTS,          &
                    NUM_TP_QUAD_POINTS,         &
                    INT_R_LOCATIONS,            &
                    INT_T_LOCATIONS,            &
                    INT_P_LOCATIONS

USE Functions_Quadrature, &
            ONLY :  Initialize_LG_Quadrature,           &
                    Initialize_LG_Quadrature_Locations

USE Maps_X_Space, &
            ONLY :  Map_To_X_Space

USE Functions_Math, &
            ONLY :  Lagrange_Poly

USE Timer_Routines_Module, &
            ONLY :  TimerStart,                         &
                    TimerStop

USE Timer_Variables_Module, &
            ONLY :  Timer_Poisson_SourceInput,          &
                    Timer_Poisson_SourceInput_PartA,    &
                    Timer_Poisson_SourceInput_PartB


IMPLICIT NONE


 !+i01+############################################################################!
!                                                                                   !
!   Poseidon_Newtonian_Source_Input -   Interface allowing for simplified source    !
!                                       input depending on dimensionality of source.!
!                                                                                   !
!-----------------------------------------------------------------------------------!
!                                                                                   !
!   Acceptable Call Forms (See Associated Routine for description of variables)     !
!                                                                                   !
!   1 Dimension :   Poseidon_Input_Sources_Poisson_1D                               !
!                                                                                   !
!           Left_Limit      -       Integer                                         !
!           Right_Limit     -       Integer                                         !
!           Num_Nodes       -       Integer                                         !
!           Input_R_Nodes   -       Real Vector, Dimension(1:Num_Nodes)             !
!           Rho             -       Real Vector, Dimension( 1:Num_Nodes,        &   !
!                                                           1:Num_R_Elements,   &   !
!                                                           1,                  &   !
!                                                           1                   )   !
!                                                                                   !
!   2 Dimensions :  Poseidon_Input_Sources_Poisson_2D                               !
!                                                                                   !
!           Left_Limit      -       Integer                                         !
!           Right_Limit     -       Integer                                         !
!           Num_Nodes       -       Integer Vector, Dimension(1:2)                  !
!           Input_R_Nodes   -       Real Vector, Dimension(1:Num_Nodes(1))          !
!           Input_T_Nodes   -       Real Vector, Dimension(1:Num_Nodes(2))          !
!           Rho             -       Real Vector, Dimension( 1:Num_Nodes,        &   !
!                                                           1:Num_R_Elements,   &   !
!                                                           1:Num_T_Elements,   &   !
!                                                           1                   )   !
!                                                                                   !
!   3 Dimensions :  Poseidon_Input_Sources_Poisson_3D                               !
!                                                                                   !
!           Left_Limit      -       Integer                                         !
!           Right_Limit     -       Integer                                         !
!           Num_Nodes       -       Integer Vector, Dimension(1:3)                  !
!           Input_R_Nodes   -       Real Vector, Dimension(1:Num_Nodes(1))          !
!           Input_T_Nodes   -       Real Vector, Dimension(1:Num_Nodes(2))          !
!           Input_P_Nodes   -       Real Vector, Dimension(1:Num_Nodes(3))          !
!           Rho             -       Real Vector, Dimension( 1:Num_Nodes,        &   !
!                                                           1:Num_R_Elements,   &   !
!                                                           1:Num_T_Elements,   &   !
!                                                           1:Num_P_Elements    )   !
!                                                                                   !
 !#################################################################################!
INTERFACE Poseidon_Newtonian_Source_Input

    PROCEDURE   Poseidon_Input_Sources_Poisson_1D,     &
                Poseidon_Input_Sources_Poisson_2D,     &
                Poseidon_Input_Sources_Poisson_3D

END INTERFACE



CONTAINS



 !+201+####################################################################################!
!                                                                                           !
!                               Poseidon_Input_Sources_Poisson_1D                          !
!                                                                                           !
!-------------------------------------------------------------------------------------------!
!                                                                                           !
!       This subroutine provides a simplified input method for source values and locations  !
!   for 1-Dimensional systems. The input process interpolates the source values using       !
!   Lagrange Interpolating Polynomials to the locations needed for the source integrations  !
!   required to construct the source vector for the linear systems created by Poseidon      !
!                                                                                           !
!       The 1-Dimensional version of this subroutine takes the input, restructures it, and  !
!   sends it along with filler variables to the subroutine,                                 !
!   Poseidon_Input_Sources_Poisson_3D (+203+). The filler variables instruct that           !
!   subroutine to act only in one dimension.                                                !
!                                                                                           !
!-------------------------------------------------------------------------------------------!
!                                                                                           !
!   Input Variables                                                                         !
!                                                                                           !
!   Left_Limit, Right_Limit -   Integers, Left and Right limits of the input quadrature     !
!                               space.                                                      !
!                                                                                           !
!   Num_Nodes               -   Integer, Number of radial nodes                             !
!                                                                                           !
!   Input_R_Quad            -   Real Vectors, Locations of radial nodes in quadrature space !
!                                                                                           !
!   Rho                     -   Real Vector, Source values at each location in each element.!
!                                                                                           !
 !#########################################################################################!
SUBROUTINE Poseidon_Input_Sources_Poisson_1D(   Left_Limit,             &
                                                Right_Limit,            &
                                                Num_Nodes,              &
                                                Input_R_Quad,           &
                                                Rho                         )





                             !                  !
                            !!  Input Variables !!
                             !                  !
REAL(KIND = idp),                               INTENT(IN)  ::  Left_Limit,             &
                                                                Right_Limit


INTEGER,                                        INTENT(IN)  ::  Num_Nodes


REAL(KIND = idp), DIMENSION(1:Num_Nodes),       INTENT(IN)  ::  Input_R_Quad


REAL(KIND = idp), DIMENSION(1:Num_Nodes,                    &
                            1:NUM_R_ELEMENTS,               &
                            1:NUM_T_ELEMENTS,               &
                            1:NUM_P_ELEMENTS),  INTENT(IN)  ::  Rho





                         !                        !
                        !!  Subroutine Variables  !!
                         !                        !

INTEGER, DIMENSION(1:3)                                     ::  Num_Nodes_Vectorized

REAL(KIND = idp), DIMENSION(1:1)                            ::  Input_T_Quad
REAL(KIND = idp), DIMENSION(1:1)                            ::  Input_P_Quad





 !                                                                   !
!!   The 3D subroutine requires a length 3 vector for Num_Nodes.     !!
 !                                                                   !
Num_Nodes_Vectorized(1) = Num_Nodes
Num_Nodes_Vectorized(2) = 1
Num_Nodes_Vectorized(3) = 1




   !                                                                     !
  !!   Sets the location of the single node in the Theta and Phi         !!
 !!!   dimensions to 1.  This choice is arbitrary as this dimension is   !!!
!!!!   being dismissed, which is achieved by using the lowert order      !!!!
 !!!   Lagrange Polynomial, L(x) = 1, which is independent of input      !!!
  !!   location.                                                         !!
   !                                                                     !
Input_T_Quad = 1.0_idp
Input_P_Quad = 1.0_idp









CALL Poseidon_Input_Sources_Poisson_3D( Left_Limit,             &
                                        Right_Limit,            &
                                        Num_Nodes_Vectorized,   &
                                        Input_R_Quad,           &
                                        Input_T_Quad,           &
                                        Input_P_Quad,           &
                                        Rho                     )








END SUBROUTINE Poseidon_Input_Sources_Poisson_1D









 !+202+####################################################################################!
!                                                                                           !
!                               Poseidon_Input_Sources_Poisson_2D                           !
!                                                                                           !
!-------------------------------------------------------------------------------------------!
!                                                                                           !
!       This subroutine provides a simplified input method for source values and locations  !
!   for 2-Dimensional systems. The input process interpolates the source values using       !
!   Lagrange Interpolating Polynomials to the locations needed for the source integrations  !
!   required to construct the source vector for the linear systems created by Poseidon      !
!                                                                                           !
!       The 2-Dimensional version of this subroutine takes the input, restructures it, and  !
!   sends it along with filler variables to the subroutine,                                 !
!   Poseidon_Input_Sources_Poisson_3D (+203+). The filler variables instruct that           !
!   subroutine to act only in two dimensions.                                               !
!                                                                                           !
!-------------------------------------------------------------------------------------------!
!                                                                                           !
!   Input Variables     :   * = R, T                                                        !
!                                                                                           !
!   Left_Limit, Right_Limit -   Integers, Left and Right limits of the input quadrature     !
!                               space.                                                      !
!                                                                                           !
!   Num_Nodes               -   Integer Vector, Number of nodes per dimension               !
!                                                                                           !
!   Input_*_Quad            -   Real Vectors, Locations of nodes in each dimension in       !
!                               quadrature space                                            !
!                                                                                           !
!   Rho                     -   Real Vector, Source values at each location in each element.!
!                                                                                           !
 !#########################################################################################!
SUBROUTINE Poseidon_Input_Sources_Poisson_2D(   Left_Limit,             &
                                                Right_Limit,            &
                                                Num_Nodes,              &
                                                Input_R_Quad,           &
                                                Input_T_Quad,           &
                                                Rho                     )


                             !                  !
                            !!  Input Variables !!
                             !                  !
REAL(KIND = idp),                               INTENT(IN)      ::  Left_Limit,             &
                                                                    Right_Limit


INTEGER, DIMENSION(1:2),                        INTENT(IN)      ::  Num_Nodes


REAL(KIND = idp), DIMENSION(1:Num_Nodes(1)),    INTENT(IN)      ::  Input_R_Quad
REAL(KIND = idp), DIMENSION(1:Num_Nodes(2)),    INTENT(IN)      ::  Input_T_Quad  !??? Right Dimension ???!


REAL(KIND = idp), DIMENSION(1:Num_Nodes(1)*Num_Nodes(2),        &
                            1:NUM_R_ELEMENTS,                   &
                            1:NUM_T_ELEMENTS,                   &
                            1:NUM_P_ELEMENTS),  INTENT(IN)      ::  Rho







                         !                        !
                        !!  Subroutine Variables  !!
                         !                        !

INTEGER, DIMENSION(1:3)                                     ::  Num_Nodes_Vectorized


REAL(KIND = idp), DIMENSION(1:1)                            ::  Input_P_Quad






 !                                                                   !
!!   The 3D subroutine requires a length 3 vector for Num_Nodes.     !!
 !                                                                   !
Num_Nodes_Vectorized(1) = Num_Nodes(1)
Num_Nodes_Vectorized(2) = Num_Nodes(2)
Num_Nodes_Vectorized(3) = 1



  !                                                                     !
 !!   Sets the location of the single node in the Phi dimension to      !!
!!!   1.  This choice is arbitrary as this dimension is being           !!!
!!!   dismissed, which is achieved by using the lowert order Lagrange   !!!
 !!   Polynomial, L(x) = 1, which is independent of input location.     !!
  !                                                                     !
Input_P_Quad = 1.0_idp





CALL Poseidon_Input_Sources_Poisson_3D( Left_Limit,             &
                                        Right_Limit,            &
                                        Num_Nodes_Vectorized,   &
                                        Input_R_Quad,           &
                                        Input_T_Quad,           &
                                        Input_P_Quad,           &
                                        Rho                     )







END SUBROUTINE Poseidon_Input_Sources_Poisson_2D
















 !+203+####################################################################################!
!                                                                                           !
!                               Poseidon_Input_Sources_Poisson_3D                           !
!                                                                                           !
!-------------------------------------------------------------------------------------------!
!                                                                                           !
!       This subroutine provides an input method for source values and locations for        !
!    3-Dimensional systems. The input process interpolates the source values using Lagrange !
!   Interpolating Polynomials to the locations needed for the source integrations required  !
!   to construct the source vectors for the linear systems created by Poseidon.             !
!                                                                                           !
!       This subroutine is capable of handling 1- and 2-Dimensional sources if the input    !
!   is properly established. To simplify the input of lower dimensional systems, the        !
!   subroutines, Poseidon_Input_Sources_Poisson_1D (+201+), and                             !
!   Poseidon_Input_Sources_Poisson_2D (+202+) have been created. These subroutines require  !
!   less input, and generate filler variables that limit this subroutine actions to the     !
!   appropriote dimensionality.                                                             !
!                                                                                           !
!-------------------------------------------------------------------------------------------!
!                                                                                           !
!   Input Variables     :   * = R, T, P                                                     !
!                                                                                           !
!   Left_Limit, Right_Limit -   Integers, Left and Right limits of the input quadrature     !
!                               space.                                                      !
!                                                                                           !
!   Num_Nodes               -   Integer Vector, Number of nodes per dimension               !
!                                                                                           !
!   Input_*_Quad            -   Real Vectors, Locations of nodes in each dimension in       !
!                               quadrature space                                            !
!                                                                                           !
!   Rho                     -   Real Vector, Source values at each location in each element.!
!                                                                                           !
 !#########################################################################################!
SUBROUTINE Poseidon_Input_Sources_Poisson_3D(   Left_Limit,                 &
                                                Right_Limit,                &
                                                Num_Nodes,                  &
                                                Input_R_Quad,               &
                                                Input_T_Quad,               &
                                                Input_P_Quad,               &
                                                Rho                           )


                             !                  !
                            !!  Input Variables !!
                             !                  !
REAL(KIND = idp), INTENT(IN)                                            ::  Left_Limit,     &
                                                                            Right_Limit


INTEGER, DIMENSION(1:3), INTENT(IN)                                     ::  Num_Nodes


REAL(KIND = idp), DIMENSION(1:Num_Nodes(1)),    INTENT(IN)              ::  Input_R_Quad
REAL(KIND = idp), DIMENSION(1:Num_Nodes(2)),    INTENT(IN)              ::  Input_T_Quad
REAL(KIND = idp), DIMENSION(1:Num_Nodes(3)),    INTENT(IN)              ::  Input_P_Quad



REAL(KIND = idp), DIMENSION(1:Num_Nodes(1)*Num_Nodes(2)*Num_Nodes(3),   &
                            1:NUM_R_ELEMENTS,                           &
                            1:NUM_T_ELEMENTS,                           &
                            1:NUM_P_ELEMENTS),  INTENT(IN)              ::  Rho










                         !                        !
                        !!  Subroutine Variables  !!
                         !                        !
INTEGER                                                     ::  re, te, pe, i,                &
                                                                Local_R, Local_T, Local_P,  &
                                                                Input_R, Input_T, Input_P,  &
                                                                Local_Here, Input_Here,     &
                                                                Num_Local_DOF,              &
                                                                Num_Input_DOF


REAL(KIND = idp), DIMENSION(1:Num_R_Quad_Points)            ::  Local_R_Locations
REAL(KIND = idp), DIMENSION(1:Num_T_Quad_Points)            ::  Local_T_Locations
REAL(KIND = idp), DIMENSION(1:Num_P_Quad_Points)            ::  Local_P_Locations

REAL(KIND = idp), DIMENSION(1:Num_Nodes(1))                 ::  Input_R_Locations
REAL(KIND = idp), DIMENSION(1:Num_Nodes(2))                 ::  Input_T_Locations
REAL(KIND = idp), DIMENSION(1:Num_Nodes(3))                 ::  Input_P_Locations



REAL(KIND = idp), DIMENSIOn(1:Num_Nodes(1))                 ::  R_Lag_Poly_Values
REAL(KIND = idp), DIMENSIOn(1:Num_Nodes(2))                 ::  T_Lag_Poly_Values
REAL(KIND = idp), DIMENSIOn(1:Num_Nodes(3))                 ::  P_Lag_Poly_Values

REAL(KIND = idp), DIMENSION(:,:), ALLOCATABLE               ::  Translation_Matrix

REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Local_Coefficients


CALL TimerStart( Timer_Poisson_SourceInput )


                          !                                                 !
                         !!                                                 !!
                        !!!     Calculate Local DOF and Allocate Space      !!!
                         !!                                                 !!
                          !                                                 !

Num_Input_DOF = Num_Nodes(1) * Num_Nodes(2) * Num_Nodes(3)

Num_Local_DOF = Num_R_Quad_Points * Num_T_Quad_Points * Num_P_Quad_Points



ALLOCATE(Translation_Matrix(1:Num_Local_DOF, 1:Num_Input_DOF))

ALLOCATE(Local_Coefficients(1:Num_Local_DOF))




                          !                                                 !
                         !!                                                 !!
                        !!!          Initialize Local Quadratures           !!!
                         !!                                                 !!
                          !                                                 !

Local_R_Locations = Int_R_Locations
Local_T_Locations = Int_T_Locations
Local_P_Locations = Int_P_Locations






                          !                                                 !
                         !!                                                 !!
                        !!!     Map Input Locations to [ -1, 1 ] Space      !!!
                         !!                                                 !!
                          !                                                 !

Input_R_Locations = Map_To_X_Space(Left_Limit, Right_Limit, Input_R_Quad)
Input_T_Locations = Map_To_X_Space(Left_Limit, Right_Limit, Input_T_Quad)
Input_P_Locations = Map_To_X_Space(Left_Limit, Right_Limit, Input_P_Quad)





                          !                                                 !
                         !!                                                 !!
                        !!!             Create Translation Matrix           !!!
                         !!                                                 !!
                          !                                                 !

CALL TimerStart( Timer_Poisson_SourceInput_PartA )
 !                               !
!!   Set Input vector location   !!
 !                               !
Local_Here = 1




DO Local_P = 1,Num_P_Quad_Points


     !                                                                               !
    !!   Calculate Lagrange Polynomial values for the current local phi location     !!
     !                                                                               !
    P_Lag_Poly_Values = Lagrange_Poly(Local_P_Locations(Local_P), Num_Nodes(3)-1, Input_P_Locations)



    DO Local_T = 1, Num_T_Quad_Points


         !                                                                               !
        !!   Calculate Lagrange Polynomial values for the current local theta location   !!
         !                                                                               !
        T_Lag_Poly_Values = Lagrange_Poly(Local_T_Locations(Local_T), Num_Nodes(2)-1, Input_T_Locations)



        DO Local_R = 1,Num_R_Quad_Points


             !                                                                               !
            !!   Calculate Lagrange Polynomial values for the current local radial location  !!
             !                                                                               !
            R_Lag_Poly_Values = Lagrange_Poly(Local_R_Locations(Local_R), Num_Nodes(1)-1, Input_R_Locations)


             !                                       !
            !!   Set/Reset Input vector location     !!
             !                                       !
            Input_Here = 1



            DO Input_P = 1,Num_Nodes(3)
            DO Input_T = 1,Num_Nodes(2)
            DO Input_R = 1,Num_Nodes(1)


              !                                                                                         !
             !!   Calculate Translation_Matrix element. Each element is the product of the Lagrange     !!
            !!!   interpolating polynomial for each dimension. When a dimension is being omited the     !!!
            !!!   lowest order Lagrange Polynomial, L(x) = 1, is used allowing for any dimensionality   !!!
             !!   to be selected.                                                                       !!
              !                                                                                         !
                Translation_Matrix(Local_Here, Input_Here)  = R_Lag_Poly_Values(Input_R)        &
                                                            * T_Lag_Poly_Values(Input_T)        &
                                                            * P_Lag_Poly_Values(Input_P)




                 !                               !
                !! Update Input vector location  !!
                 !                               !
                Input_Here = Input_Here + 1


            END DO  !   Input_R Loop
            END DO  !   Input_T Loop
            END DO  !   Input_P Loop




             !                               !
            !! Update Local vector location  !!
             !                               !
            Local_Here = Local_Here + 1


        END DO  !   LocaL_R Loop

    END DO  !   Local_T Loop

END DO  !   Local_P looop



CALL TimerStop(  Timer_Poisson_SourceInput_PartA )
CALL TimerStart( Timer_Poisson_SourceInput_PartB )






                              !                                                     !
                             !!                                                     !!
                            !!!     Multiply Input Values and Translation Matrix    !!!
                             !!                                                     !!
                              !                                                     !


Source_Terms = 0.0_idp
#if defined(POSEIDON_OPENMP_OL_FLAG)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(6) &
    !$OMP PRIVATE(  ) &
    !$OMP REDUCTION( MIN: TimeStep )
#elif defined(POSEIDON_OPENACC_FLAG)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(6) &
    !$ACC PRIVATE(  ) &
    !$ACC PRESENT( ) &
    !$ACC REDUCTION( MIN: TimeStep )
#elif defined(POSEIDON_OPENMP_FLAG)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)  &
    !$OMP PRIVATE(  re,te,pe, i   )
#endif

DO pe = 0, NUM_P_ELEMENTS-1
DO te = 0, NUM_T_ELEMENTS-1
DO re = 0, NUM_R_ELEMENTS-1
DO i = 1,Num_Local_DOF

    Source_Terms(i,re,te,pe) = SUM( Translation_Matrix(i,:)     &
                                    * Rho(:,re+1,te+1,pe+1)     )


END DO
END DO  ! RE
END DO  ! TE
END DO  ! PE


#if defined(POSEIDON_OPENMP_OL_FLAG)
    !$OMP END PARALLEL DO SIMD
#elif defined(POSEIDON_OPENACC_FLAG)
    !$ACC END PARALLEL LOOP
#elif defined(POSEIDON_OPENMP_FLAG)
    !$OMP END PARALLEL DO SIMD
#endif



CALL TimerStop(Timer_Poisson_SourceInput_PartB)
CALL TimerStop(Timer_Poisson_SourceInput)

 !                                                              !
!!  Deallocate Translation Matrix, and Local Coefficient Vector !!
 !                                                              !

DEALLOCATE(Translation_Matrix, Local_Coefficients)




END SUBROUTINE Poseidon_Input_Sources_Poisson_3D





END MODULE
