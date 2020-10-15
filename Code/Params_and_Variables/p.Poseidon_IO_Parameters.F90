   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_IO_Parameters                                                       !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains parameters used in Poseidon's input and output routines.           !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!



USE Poseidon_Constants_Module, &
            ONLY :  idp


IMPLICIT NONE


CHARACTER(LEN = 24), PARAMETER        :: Poseidon_Reports_Dir   = "Poseidon_Output/Reports/"
CHARACTER(LEN = 42), PARAMETER        :: Poseidon_IterReports_Dir = "Poseidon_Output/Reports/Iteration_Reports/"


CHARACTER(LEN = 24), PARAMETER        :: Poseidon_Results_Dir   = "Poseidon_Output/Results/"


CHARACTER(LEN = 24), PARAMETER        :: Poseidon_Objects_Dir   = "Poseidon_Output/Objects/"
CHARACTER(LEN = 32), PARAMETER        :: Poseidon_Sources_Dir   = "Poseidon_Output/Objects/Sources/"
CHARACTER(LEN = 33), PARAMETER        :: Poseidon_Residual_Dir  = "Poseidon_Output/Objects/Residual/"
CHARACTER(LEN = 38), PARAMETER        :: Poseidon_LinSys_Dir    = "Poseidon_Output/Objects/Linear_System/"



CHARACTER(LEN = 32),DIMENSION(1:5), PARAMETER   :: &
CFA_Var_Names = ['Conformal Factor                  ', &
                 'Lapse Function                    ', &
                 'Radial Shift Component            ', &
                 'Theta Shift Component             ', &
                 'Phi Shift Component               ' ]

CHARACTER(LEN = 32),DIMENSION(1:3), PARAMETER   :: &
Solver_Names = ['Newton-Raphson                  ', &
                'Fixed Point Iteration           ', &
                'Jacobian Free GMRES             ' ]

END MODULE Poseidon_IO_Parameters
