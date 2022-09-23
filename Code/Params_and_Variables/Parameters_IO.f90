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



USE Poseidon_Kinds_Module, &
                    ONLY : idp


IMPLICIT NONE


CHARACTER(LEN = 24), PARAMETER        :: Poseidon_Reports_Dir     = "Poseidon_Output/Reports/"
CHARACTER(LEN = 42), PARAMETER        :: Poseidon_IterReports_Dir = "Poseidon_Output/Reports/Iteration_Reports/"


CHARACTER(LEN = 24), PARAMETER        :: Poseidon_Results_Dir   = "Poseidon_Output/Results/"


CHARACTER(LEN = 24), PARAMETER        :: Poseidon_Objects_Dir   = "Poseidon_Output/Objects/"
CHARACTER(LEN = 29), PARAMETER        :: Poseidon_Mesh_Dir      = "Poseidon_Output/Objects/Mesh/"
CHARACTER(LEN = 32), PARAMETER        :: Poseidon_Sources_Dir   = "Poseidon_Output/Objects/Sources/"
CHARACTER(LEN = 33), PARAMETER        :: Poseidon_Residual_Dir  = "Poseidon_Output/Objects/Residual/"
CHARACTER(LEN = 38), PARAMETER        :: Poseidon_LinSys_Dir    = "Poseidon_Output/Objects/Linear_System/"
CHARACTER(LEN = 37), PARAMETER        :: Poseidon_Coeffs_Dir    = "Poseidon_Output/Objects/Coefficients/"


INTEGER, PARAMETER                                      :: N_Variables = 8
CHARACTER(LEN = 34),DIMENSION(N_Variables), PARAMETER   :: &
CFA_Var_Names = ['Conformal Factor                  ', &
                 'Lapse Function                    ', &
                 'Radial Shift Component            ', &
                 'Theta Shift Component             ', &
                 'Phi Shift Component               ', &
                 'Radial X Component                ', &
                 'Theta X Component                 ', &
                 'Phi X Component                   ' ]


INTEGER, PARAMETER                                      :: N_VecVars = 2
CHARACTER(LEN = 15),DIMENSION(N_VecVars), PARAMETER   :: &
CFA_VecVar_Names = ['Shift Component', &
                    'X Component    ']


INTEGER, PARAMETER                                      :: N_ShortVars = 8
CHARACTER(LEN = 10),DIMENSION(N_ShortVars), PARAMETER   :: &
CFA_ShortVars = ['ConFactor ', &
                 'Lapse     ', &
                 'Beta1     ', &
                 'Beta2     ', &
                 'Beta3     ', &
                 'X1        ', &
                 'X2        ', &
                 'X3        ' ]


INTEGER, PARAMETER                                          :: N_Kij_ShortVars = 6
CHARACTER(LEN = 2),DIMENSION(N_Kij_ShortVars), PARAMETER    :: &
Kij_ShortVars = ['K11', &
                 'K12', &
                 'K13', &
                 'K21', &
                 'K22', &
                 'K33' ]


INTEGER, PARAMETER                                      :: N_Methods = 3
CHARACTER(LEN = 32),DIMENSION(N_Methods), PARAMETER     :: &
Method_Names = ['Newton-Raphson                  ', &
                'Fixed Point Iteration           ', &
                'XCFC with Fixed Point Iteration ' ]

INTEGER, PARAMETER                                      ::  N_Meshes = 5
CHARACTER(LEN = 35),DIMENSION(N_Meshes), PARAMETER      ::  &
Mesh_Names = [  'Uniform                            ',      &
                'Logarithmic                        ',      &
                'Uniform Core & Logarithmic Exterior',      &
                'Zoom                               ',      &
                'MacLaurin Shell                    '       ]


END MODULE Poseidon_IO_Parameters
