   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Constraint_Equations_Module                                                  !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!        Contains the functions used to calculate the components of the          !##!
!##!    extended Jacobian matrix as well as the derivative coefficients. These      !##!
!##!    are used to construct the stiffness matrix.                                 !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   CFA_3D_Master_Build                                                 !##!
!##!                                                                                !##!
!##!    +201+   CREATE_3D_NONLAPLACIAN_SOE                                          !##!
!##!    +202+   Calc_3D_Current_Values                                              !##!
!##!    +203+   Calc_3D_SubJcbn_Terms                                               !##!
!##!    +204+   CREATE_3D_RHS_VECTOR                                                !##!
!##!    +205+   CREATE_3D_JCBN_MATRIX                                               !##!
!##!                                                                                !##!
!##!    +301+   REDUCE_3D_NONLAPLACIAN_SOE                                          !##!
!##!                                                                                !##!
!##!    +401+   FINISH_3D_JACOBIAN_MATRIX                                           !##!
!##!                                                                                !##!
!##!    +501+   FINISH_3D_RHS_VECTOR                                                !##!
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



USE MPI

USE OMP_LIB

USE Units_Module, &
                                ONLY :  GR_Source_Scalar

USE Poseidon_Constants_Module, &
                                ONLY :  idp,                &
                                        pi,                 &
                                        TwoPi,              &
                                        OneThird,           &
                                        TwoThirds,          &
                                        FourThirds,         &
                                        OneThirtySecond

USE IO_Functions_Module, &
                                ONLY :  OPEN_NEW_FILE






 !+101+############################################################################!
!                                                                                   !
!                     OUTPUT_CONSTRAINT_EQUATION_RESULTS                            !
!                                                                                   !
 !#################################################################################!
SUBROUTINE OUTPUT_CONSTRAINT_EQUATION_RESULTS()


INTEGER                                                     ::  NUM_SAMPLES


CHARACTER(LEN = 50), DIMENSION(:), ALLOCATABLE              ::  Filenames
INTEGER, DIMENSION(:), ALLOCATABLE                          ::  File_IDs
INTEGER                                                     ::  Num_Files

INTEGER                                                     ::  NUM_THETA_RAYS,     &
                                                                NUM_PHI_RAYS,       &
                                                                NUM_RADIAL_SAMPLES

REAL(KIND = idp), DIMENSION(:,:,:), ALLOCATABLE             ::  Lapse_Holder,       &
                                                                ConForm_Holder
REAL(KIND = idp), DIMENSION(:,:,:,:), ALLOCATABLE           ::  Shift_Holder
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  R_Holder,           &
                                                                T_Holder,           &
                                                                P_Holder





Num_Files = 7

ALLOCATE( Filenames(1:Num_Files) )
ALLOCATE( File_IDs(1:Num_Files) )

WRITE(Filenames(1),'(A34,I5.5,A4)')"OUTPUT/CONSTRAINT_EQS/Hamiltonian_",DRIVER_FRAME,".out"
WRITE(Filenames(2),'(A33,I5.5,A4)')"OUTPUT/CONSTRAINT_EQS/Momentum_1_",DRIVER_FRAME,".out"
WRITE(Filenames(3),'(A33,I5.5,A4)')"OUTPUT/CONSTRAINT_EQS/Momentum_1_",DRIVER_FRAME,".out"
WRITE(Filenames(4),'(A33,I5.5,A4)')"OUTPUT/CONSTRAINT_EQS/Momentum_1_",DRIVER_FRAME,".out"
WRITE(Filenames(5),'(A32,I5.5,A4)')"OUTPUT/CONSTRAINT_EQS/R_VALUES_",DRIVER_FRAME,".out"
WRITE(Filenames(6),'(A32,I5.5,A4)')"OUTPUT/CONSTRAINT_EQS/T_VALUES_",DRIVER_FRAME,".out"
WRITE(Filenames(7),'(A32,I5.5,A4)')"OUTPUT/CONSTRAINT_EQS/P_VALUES_",DRIVER_FRAME,".out"


DO i = 1,Num_Files
    CALL OPEN_NEW_FILE( Filenames(i), File_IDs(i) )
END DO




! Create Output Spacing
NUM_PHI_RAYS = 1
IF ( NUM_PHI_RAYS == 1 ) THEN
    DELTA_PHI = pi/2.0_idp
!    DELTA_PHI = 0.0_idp
ELSE
    DELTA_PHI = 2.0_idp*pi/(NUM_PHI_RAYS-1)
END IF

NUM_THETA_RAYS = 20
DELTA_THETA = pi/(NUM_THETA_RAYS+1)

NUM_RADIAL_SAMPLES = 1000
CALL Create_Logarithmic_1D_Mesh( R_INNER, R_OUTER, NUM_RADIAL_SAMPLES,  &
                                 output_re, output_rc, output_dr        )

ALLOCATE( Lapse_Holder(1:NUM_PHI_RAYS, 1:NUM_THETA_RAYS, 1:NUM_RADIAL_SAMPLES) )
ALLOCATE( ConForm_Holder(1:NUM_PHI_RAYS, 1:NUM_THETA_RAYS, 1:NUM_RADIAL_SAMPLES) )
ALLOCATE( Shift_Holder(1:3,1:NUM_PHI_RAYS, 1:NUM_THETA_RAYS, 1:NUM_RADIAL_SAMPLES) )
ALLOCATE( R_Holder(1:NUM_RADIAL_SAMPLES) )
ALLOCATE( T_Holder(1:NUM_THETA_RAYS) )
ALLOCATE( P_Holder(1:NUM_PHI_RAYS) )









END SUBROUTINE OUTPUT_CONSTRAINT_EQUATION_RESULTS














END MODULE Constraint_Equations_Module
