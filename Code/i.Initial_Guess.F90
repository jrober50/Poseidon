   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Initial_Guess_Module                                                         !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   Initialize_Flat_Space_Guess_Values                                  !##!
!##!    +102+   Initialize_Special_Guess_Values                                     !##!
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
                        ONLY :  idp, pi,                &
                                TwoThirds,              &
                                FourThirds

USE Units_Module, &
                        ONLY :  C_Square

USE Poseidon_Parameters, &
                        ONLY :  DOMAIN_DIM,             &
                                DEGREE,                 &
                                NUM_CFA_VARS

USE DRIVER_Parameters, &
                        ONLY :  myID,                   &
                                Potential_Solution,     &
                                Shift_Solution,         &
                                Potential_Sol_Flag,     &
                                Shift_Sol_Flag,         &
                                Solver_Mode

USE Poseidon_Variables_Module, &
                        ONLY :  NUM_R_ELEMENTS,         &
                                rlocs,                  &
                                Matrix_Location,        &
                                Coefficient_Vector

USE Driver_Additional_Functions_Module, &
                        ONLY :  Map_From_X_Space,                       &
                                Initialize_LGL_Quadrature_Locations

USE Poseidon_IO_Module, &
                        ONLY :  OPEN_EXISTING_FILE

USE Poseidon_FP_Mapping_Functions_Module,   &
                        ONLY :  FP_Vector_Map

USE Poseidon_FP_Variables_Module, &
                        ONLY :  FP_Coeff_Vector


IMPLICIT NONE

CONTAINS

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


INTEGER                                     :: re, rd, Here



PRINT*,"Initializing Flat Space Guess"
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



IF ( SOLVER_MODE == 2 ) THEN

    FP_Coeff_Vector = 0.0_idp
    DO re = 0,NUM_R_ELEMENTS - 1

        DO rd = 0, Degree
            ! 2 sqrt(pi) is Ylm normalization factor

            Here = FP_Vector_Map(re, rd)

            FP_Coeff_Vector(Here,0,1) = 1.0_idp * 2.0_idp * sqrt(pi)
            FP_Coeff_Vector(Here,0,2) = 1.0_idp * 2.0_idp * sqrt(pi)
            FP_Coeff_Vector(Here,0,3) = Beta_Start - delta_Beta*(re*DEGREE+rd)
    !        Coefficient_Vector(Cur_Shift_Loc) = 0.0_idp
            
        END DO
    END DO

    

ELSE



    DO re = 0,NUM_R_ELEMENTS - 1


        DO rd = 0,DEGREE


            ! 2 sqrt(pi) is Ylm normalization factor


            CUR_PSI_LOC = Matrix_Location( 1, 0, 0, re, rd )
            Coefficient_Vector(CUR_PSI_LOC) = 1.0_idp * 2.0_idp * sqrt(pi)




            CUR_ALPHPSI_LOC = Matrix_Location( 2, 0, 0, re, rd )
            Coefficient_Vector(CUR_ALPHPSI_LOC) = 1.0_idp * 2.0_idp * sqrt(pi)


            DO beta_i = 1,DOMAIN_DIM-1

                CUR_BETA_LOC = Matrix_Location( 2+beta_i, 0, 0, re, rd )
                Coefficient_Vector(CUR_BETA_LOC) = Beta_Start - delta_Beta*(re*DEGREE+rd)

            END DO


        END DO


    END DO

END IF


END SUBROUTINE Initialize_Flat_Space_Guess_Values








!+102+###########################################################################!
!                                                                                !
!                  Initialize_Calculated_Guess_Values                            !
!                                                                                !
!################################################################################!
SUBROUTINE Initialize_Calculated_Guess_Values()


INTEGER                                                         :: re, d, Here


REAL(KIND = idp), DIMENSION(0:DEGREE)                           :: R_Values


REAL(KIND = idp), DIMENSION(0:DEGREE)                           ::  Local_Locations


INTEGER                                                         ::  CUR_PSI_LOC,    &
                                                                    CUR_ALPHPSI_LOC,&
                                                                    CUR_SHIFT_LOC



Print*,"Initialize_Calculated_Guess_Values"


Local_Locations = Initialize_LGL_Quadrature_Locations(DEGREE)


!
!   Empty Space Initial Guess
!
Coefficient_Vector = 0.0_idp
IF ( SOLVER_MODE == 2 ) THEN

    FP_Coeff_Vector = 0.0_idp

    DO re = 0,NUM_R_ELEMENTS - 1

        R_Values = Map_From_X_Space(rlocs(re), rlocs(re + 1), Local_Locations)

        DO d = 0,DEGREE
     
            ! 2 sqrt(pi) is Ylm normalization factor

            Here = FP_Vector_Map(re, d)

            FP_Coeff_Vector(Here,0,1) = 2.0_idp * sqrt(pi)                                          &
                                            * ( 1.0_idp - 0.5_idp                                       &
                                                * Potential_Solution(R_Values(d),0.0_idp,0.0_idp)/C_Square  )

            FP_Coeff_Vector(Here,0,2) = 2.0_idp * sqrt(pi)                                          &
                                            * ( 1.0_idp + 0.5_idp                                       &
                                                * Potential_Solution(R_Values(d),0.0_idp,0.0_idp)/C_Square  )

            FP_Coeff_Vector(Here,0,3) = 2.0_idp*sqrt(pi)*Shift_Solution(R_Values(d),rlocs,NUM_R_ELEMENTS)
    !        Coefficient_Vector(Cur_Shift_Loc) = 0.0_idp
            
        END DO
    END DO



ELSE
    DO re = 0,NUM_R_ELEMENTS - 1

        R_Values = Map_From_X_Space(rlocs(re), rlocs(re + 1), Local_Locations)

        DO d = 0,DEGREE
     
            ! 2 sqrt(pi) is Ylm normalization factor

            CUR_PSI_LOC = Matrix_Location( 1, 0, 0, re, d )
            Coefficient_Vector(CUR_PSI_LOC) = 2.0_idp * sqrt(pi)                                        &
                                            * ( 1.0_idp - 0.5_idp                                       &
                                                * Potential_Solution(R_Values(d),0.0_idp,0.0_idp)/C_Square  )



            CUR_ALPHPSI_LOC = Matrix_Location( 2, 0, 0, re, d )
            
            Coefficient_Vector(CUR_ALPHPSI_LOC) = 2.0_idp * sqrt(pi)                                    &
                                            * ( 1.0_idp + 0.5_idp                                       &
                                                * Potential_Solution(R_Values(d),0.0_idp,0.0_idp)/C_Square  )


            CUR_SHIFT_LOC = Matrix_Location( 3, 0, 0, re, d )
            Coefficient_Vector(Cur_Shift_Loc) = 2.0_idp*sqrt(pi)*Shift_Solution(R_Values(d),rlocs,NUM_R_ELEMENTS)
    !        Coefficient_Vector(Cur_Shift_Loc) = 0.0_idp
            
        END DO
    END DO
END IF


END SUBROUTINE Initialize_Calculated_Guess_Values









!+102+###########################################################################!
!                                                                                !
!                  Load_Initial_Guess_From_File                                  !
!                                                                                !
!################################################################################!
SUBROUTINE Load_Initial_Guess_From_File()


INTEGER                                                 ::  Frame_Num, Iter_Num

CHARACTER(LEN = 57)                                     ::  FILE_NAME
CHARACTER(LEN = 40)                                     ::  fmt

INTEGER                                                 ::  FILE_ID
INTEGER                                                 ::  i

INTEGER                                                 ::  istat


100 FORMAT (A,I2.2,A,I2.2,A)
101 FORMAT (I5.5," ",I2.2," ",I2.2)

fmt = '(ES24.16E3,SP,ES24.16E3,"i")'
!fmt = '(F16.10,SP,F16.10,"i")'



Frame_Num = 1
Iter_Num = 2


WRITE(FILE_NAME,100)"Data3/Poseidon_Coeffs/COEFF_VEC_F",Frame_Num,"_I",Iter_Num,".out"
CALL OPEN_EXISTING_FILE( FILE_NAME, FILE_ID, istat )

READ(FILE_ID,*)Coefficient_Vector


CLOSE(FILE_ID)

END SUBROUTINE Load_Initial_Guess_From_File










!+401+###########################################################################!
!                                                                                !
!                  CALC_SHIFT_BC_VALUE                                           !
!                                                                                !
!################################################################################!
SUBROUTINE CALC_SHIFT_BC_VALUE()



IF ( Shift_Sol_Flag .eqv. .FALSE. ) THEN
    IF ( Potential_Sol_Flag .eqv. .FALSE. ) THEN

        


    END IF
END IF

END SUBROUTINE CALC_SHIFT_BC_VALUE









!!+401+###########################################################################!
!!                                                                                !
!!                  Inital_Guess_Shift_Vector_1D                               !
!!                                                                                !
!!################################################################################!
!SUBROUTINE Inital_Guess_Shift_Vector_1D( Num_RE, Num_TE, Num_PE,   &
!                                            Num_Nodes,                &
!                                            Mesh_Width, Mesh_Center, Node_Locs,    &
!                                            Psi, Alpha, Sr,           &
!                                            Shift                     )
!
!INTEGER, INTENT(IN)                                                                 ::  Num_RE, &
!                                                                                        Num_TE, &
!                                                                                        Num_PE
!INTEGER, INTENT(IN)                                                                 ::  Num_Nodes
!
!REAL(KIND = idp), DIMENSION(1:Num_RE), INTENT(IN)                                   ::  Mesh_Width
!REAL(KIND = idp), DIMENSION(1:Num_RE), INTENT(IN)                                   ::  Mesh_Center
!REAL(KIND = idp), DIMESNION(1:Num_Nodes), INTENT(IN)                                ::  Node_Locs
!
!REAL(KIND = idp), DIMENSION(1:Num_Nodes, 1:Num_RE, 1:Num_TE, 1:Num_PE), INTENT(IN)  ::  Psi
!REAL(KIND = idp), DIMENSION(1:Num_Nodes, 1:Num_RE, 1:Num_TE, 1:Num_PE), INTENT(IN)  ::  Alpha
!REAL(KIND = idp), DIMENSION(1:Num_Nodes, 1:Num_RE, 1:Num_TE, 1:Num_PE), INTENT(IN)  ::  Sr
!
!REAL(KIND = idp), DIMENSION(1:Num_Nodes, 1:Num_RE, 1:Num_TE, 1:Num_PE), INTENT(OUT) ::  Shift
!
!
!
!
!INTEGER                                                           ::  Ord
!INTEGER                                                           ::  i, j, l, m, d, re, reb
!
!
!
!REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                       :: x_locs,     &
!                                                                     ri_locs,    &
!                                                                     wi
!
!REAL(KIND = idp), DIMENSION(:,:), ALLOCATABLE                     :: rij_locs,    &
!                                                                     PSI_10,     &
!                                                                     Sr_New
!
!REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                       :: AlphaPsi
!REAL(KIND = idp)                                                  :: Psi
!REAL(KIND = idp)                                                  :: Outer_Int
!REAL(KIND = idp)                                                  :: Inner_Int
!
!Ord = 6
!
!ALLOCATE( x_locs(1:Ord) )
!ALLOCATE( ri_locs(1:Ord) )
!ALLOCATE( wi(1:Ord)  )
!ALLOCATE( AlphaPsi(1:Ord) )
!
!ALLOCATE( rij_locs(1:Ord,1:Ord) )
!ALLOCATE( PSI_10(1:Ord,1:Ord)  )
!ALLOCATE( Sr_New(1:Ord,1:Ord)  )
!
!
!CALL Initialize_LG_Quadrature( Ord, x_locs, wi )
!
!
!
!DO re = 1,Num_RE
!
!    ! Calculate the r locations for the Outer Integral's Quadrature Points !
!    ri_locs(:) = Mesh_Width(re)*(x_locs(:)/2.0_idp)  + Mesh_Center(re)
!
!
!
!
!END DO
!
!
!
!
!
!
!
!END SUBROUTINE Inital_Guess_Shift_Vector_1D









END MODULE Initial_Guess_Module
