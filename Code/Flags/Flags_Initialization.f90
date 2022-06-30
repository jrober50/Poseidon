   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Flags_Initialization_Module                                                  !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains parameters used to define the running of Poseidon.                 !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!

USE Poseidon_Kinds_Module, &
            ONLY :  idp

USE Poseidon_Message_Routines_Module,   &
            ONLY :  Warning_Message

USE Flags_Core_Module, &
            ONLY :  lPF_Core_Flags,         &
                    iPF_Core_Poisson_Mode

IMPLICIT NONE

!  Initialization Flags
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Num_Flags          = 13

!INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Units_Set          = 1
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Expansion_Params   = 1
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Matrices           = 2
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_IO_Params          = 3
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_FP_Params          = 4
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Mesh               = 5
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Quadrature         = 6
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Tables             = 7
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_MPI                = 8
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_MTGV               = 9
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Method_Vars        = 10
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Alloc_LinSys       = 11
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Caller_Vars        = 12
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Alloc_Source       = 13

LOGICAL,    PUBLIC, DIMENSION(1:iPF_Init_Num_Flags)     ::  lPF_Init_Flags





!  Initialization Flags - Matrices
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Matrices_Num_Flags         = 4

INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Matrices_Type_A            = 1
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Matrices_Type_A_Cholesky   = 2
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Matrices_Type_B            = 3
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Matrices_Type_B_LU         = 4

LOGICAL,    PUBLIC, DIMENSION(1:iPF_Init_Matrices_Num_Flags)    ::  lPF_Init_Matrices_Flags






!  Initialization Flags - Quadrature
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Quad_Num_Flags = 3

INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Quad_Params    = 1
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Quad_Alloc     = 2
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Quad_Init      = 3

LOGICAL,    PUBLIC, DIMENSION(1:iPF_Init_Quad_Num_Flags)        ::  lPF_Init_Quad_Flags



!  Initialization Flags - Mesh
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Mesh_Num_Flags = 3

INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Mesh_Params    = 1
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Mesh_Alloc     = 2
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Mesh_Init      = 3

LOGICAL,    PUBLIC, DIMENSION(1:iPF_Init_Mesh_Num_Flags)        ::  lPF_Init_Mesh_Flags



!  Initialization Flags - Tables
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Tables_Num_Flags = 2

INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Tables_Alloc     = 1
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_Tables_Init      = 2

LOGICAL,    PUBLIC, DIMENSION(1:iPF_Init_Tables_Num_Flags)        ::  lPF_Init_Tables_Flags



!  Initialization Flags - Maps, Tables and General Variables (MTGV)
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_MTGV_Num_Flags = 2

INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_MTGV_Derived   = 1
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_MTGV_TransMat  = 2

LOGICAL,    PUBLIC, DIMENSION(1:iPF_Init_MTGV_Num_Flags)        ::  lPF_Init_MTGV_Flags



!  Initialization Flags - AMReX
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_AMReX_Num_Flags = 2

INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_AMReX_Params    = 1
INTEGER,    PUBLIC, PARAMETER       ::  iPF_Init_AMReX_Maps      = 2

LOGICAL,    PUBLIC, DIMENSION(1:iPF_Init_AMReX_Num_Flags)     ::  lPF_Init_AMReX_Flags





CONTAINS



 !+101+####################################################!
!                                                           !
!          Poseidon_Initialization_Check                    !
!                                                           !
 !#########################################################!
LOGICAL FUNCTION Poseidon_Initialization_Check()


lPF_Init_Flags(iPF_Init_Matrices) = Poseidon_Init_Matrices_Check()
lPF_Init_Flags(iPF_Init_Quadrature) = Poseidon_Init_Quadrature_Check()
lPF_Init_Flags(iPF_Init_MTGV) = Poseidon_Init_MTGV_Check()
lPF_Init_Flags(iPF_Init_Mesh) = Poseidon_Init_Mesh_Check()
lPF_Init_Flags(iPF_Init_Tables) = Poseidon_Init_Tables_Check()



IF ( ALL(lPF_Init_Flags) ) THEN
    Poseidon_Initialization_Check = .TRUE.
ELSE
!    IF ( .NOT. lPF_Init_Flags(iPF_Init_Units_Set))         &
!        CALL Warning_Message('Poseidon Initialization Check Failed. Units Not Set.')

    IF ( .NOT. lPF_Init_Flags(iPF_Init_Expansion_Params))         &
        CALL Warning_Message('Poseidon Initialization Check Failed. Expansion parameters not set.')

    IF ( .NOT. lPF_Init_Flags(iPF_Init_Matrices)) THEN
        CALL Warning_Message('Poseidon Initialization Check Failed. Matrices not initialized.')
    END IF

    IF ( .NOT. lPF_Init_Flags(iPF_Init_IO_Params))         &
        CALL Warning_Message('Poseidon Initialization Check Failed. IO parameters not set.')

    IF ( .NOT. lPF_Init_Flags(iPF_Init_FP_Params))         &
        CALL Warning_Message('Poseidon Initialization Check Failed. Fixed Point parameters not set.')

    IF ( .NOT. lPF_Init_Flags(iPF_Init_Mesh))         &
        CALL Warning_Message('Poseidon Initialization Check Failed. Mesh not initialized.')

    IF ( .NOT. lPF_Init_Flags(iPF_Init_Quadrature))         &
        CALL Warning_Message('Poseidon Initialization Check Failed. Quadrature not initialized.')

    IF ( .NOT. lPF_Init_Flags(iPF_Init_Tables))         &
        CALL Warning_Message('Poseidon Initialization Check Failed. Tables not initialized.')

    IF ( .NOT. lPF_Init_Flags(iPF_Init_MPI))         &
        CALL Warning_Message('Poseidon Initialization Check Failed. MPI variables not initialized.')

    IF ( .NOT. lPF_Init_Flags(iPF_Init_MTGV))       &
        CALL Warning_Message('Poseidon Initialization Check Failed. Maps, tables and general variables not initialized.')

    IF ( .NOT. lPF_Init_Flags(iPF_Init_Method_Vars))         &
        CALL Warning_Message('Poseidon Initialization Check Failed. Method variables not initalized.')

    IF ( .NOT. lPF_Init_Flags(iPF_Init_Alloc_LinSys))         &
        CALL Warning_Message('Poseidon Initialization Check Failed. Linear system variables not allocated.')

    IF ( .NOT. lPF_Init_Flags(iPF_Init_Caller_Vars))         &
        CALL Warning_Message('Poseidon Initialization Check Failed. Caller variables not initalized.')

    IF ( .NOT. lPF_Init_Flags(iPF_Init_Alloc_Source))         &
        CALL Warning_Message('Poseidon Initialization Check Failed. Source variables not allocated.')

    Poseidon_Initialization_Check = .FALSE.
END IF

END FUNCTION Poseidon_Initialization_Check



 !+201+####################################################!
!                                                           !
!          Poseidon_Init_Matrices_Check                     !
!                                                           !
 !#########################################################!
LOGICAL FUNCTION Poseidon_Init_Matrices_Check()

IF ( lPF_Core_Flags(iPF_Core_Poisson_Mode) ) THEN
    Poseidon_Init_Matrices_Check = Poseidon_Init_Matrices_Check_Poisson()
ELSE
    Poseidon_Init_Matrices_Check = Poseidon_Init_Matrices_Check_XCFC()
END IF ! lPF_Flags_Core(iPF_Core_Poisson_Mode)

END FUNCTION Poseidon_Init_Matrices_Check


 !+202+####################################################!
!                                                           !
!          Poseidon_Init_Matrices_Check_XCFC                !
!                                                           !
 !#########################################################!
LOGICAL FUNCTION Poseidon_Init_Matrices_Check_XCFC()

IF ( ALL( [ lPF_Init_Matrices_Flags(iPF_Init_Matrices_Type_A), &
            lPF_Init_Matrices_Flags(iPF_Init_Matrices_Type_B)] ) ) THEN
    Poseidon_Init_Matrices_Check_XCFC = .TRUE.
ELSE

    IF ( .NOT. lPF_Init_Matrices_Flags(iPF_Init_Matrices_Type_A))         &
        CALL Warning_Message('Poseidon Initialization Check Failed. Scalar laplace matrix not built.')

!    IF ( .NOT. lPF_Init_Matrices_Flags(iPF_Init_Matrices_Type_A_Cholesky))         &
!        CALL Warning_Message('Poseidon Initialization Check Failed. Scalar laplace matrix not factorized.')

    IF ( .NOT. lPF_Init_Matrices_Flags(iPF_Init_Matrices_Type_B))         &
        CALL Warning_Message('Poseidon Initialization Check Failed. Vector laplace matrix not built.')

!    IF ( .NOT. lPF_Init_Matrices_Flags(iPF_Init_Matrices_Type_B_LU))         &
!        CALL Warning_Message('Poseidon Initialization Check Failed. Vector laplace matrix not factorized.')

    Poseidon_Init_Matrices_Check_XCFC = .FALSE.
END IF


END FUNCTION Poseidon_Init_Matrices_Check_XCFC


 !+202+####################################################!
!                                                           !
!          Poseidon_Init_Matrices_Check_Poisson             !
!                                                           !
 !#########################################################!
LOGICAL FUNCTION Poseidon_Init_Matrices_Check_Poisson()

IF ( lPF_Init_Matrices_Flags(iPF_Init_Matrices_Type_A) ) THEN
    Poseidon_Init_Matrices_Check_Poisson = .TRUE.
ELSE

    IF ( .NOT. lPF_Init_Matrices_Flags(iPF_Init_Matrices_Type_A))         &
        CALL Warning_Message('Poseidon Initialization Check Failed. Scalar laplace matrix not built.')

!    IF ( .NOT. lPF_Init_Matrices_Flags(iPF_Init_Matrices_Type_A_Cholesky))         &
!        CALL Warning_Message('Poseidon Initialization Check Failed. Scalar laplace matrix not factorized.')

    Poseidon_Init_Matrices_Check_Poisson = .FALSE.
END IF

END FUNCTION Poseidon_Init_Matrices_Check_Poisson

 !+301+####################################################!
!                                                           !
!          Poseidon_Init_Quadrature_Check                   !
!                                                           !
 !#########################################################!
LOGICAL FUNCTION Poseidon_Init_Quadrature_Check()

IF ( ALL(lPF_Init_Quad_Flags) ) THEN
    Poseidon_Init_Quadrature_Check = .TRUE.
ELSE
    
    IF ( .NOT. lPF_Init_Quad_Flags(iPF_Init_Quad_Params))         &
        CALL Warning_Message('Poseidon Initialization Check Failed. Quadrature parameters not set.')

    IF ( .NOT. lPF_Init_Quad_Flags(iPF_Init_Quad_Alloc))         &
        CALL Warning_Message('Poseidon Initialization Check Failed. Quadrature variables not allocated.')

    IF ( .NOT. lPF_Init_Quad_Flags(iPF_Init_Quad_Init))         &
        CALL Warning_Message('Poseidon Initialization Check Failed. Quadrature variables not initialized.')


    Poseidon_Init_Quadrature_Check = .FALSE.
END IF

END FUNCTION Poseidon_Init_Quadrature_Check




 !+401+####################################################!
!                                                           !
!          Poseidon_Init_Mesh_Check                         !
!                                                           !
 !#########################################################!
LOGICAL FUNCTION Poseidon_Init_Mesh_Check()

IF ( ALL(lPF_Init_Mesh_Flags) ) THEN
    Poseidon_Init_Mesh_Check = .TRUE.
ELSE
    
    IF ( .NOT. lPF_Init_Mesh_Flags(iPF_Init_Mesh_Params))         &
        CALL Warning_Message('Poseidon Initialization Check Failed. Mesh parameters not set.')

    IF ( .NOT. lPF_Init_Mesh_Flags(iPF_Init_Mesh_Alloc))         &
        CALL Warning_Message('Poseidon Initialization Check Failed. Mesh variables not allocated.')

    IF ( .NOT. lPF_Init_Mesh_Flags(iPF_Init_Mesh_Init))         &
        CALL Warning_Message('Poseidon Initialization Check Failed. Mesh variables not initialized.')


    Poseidon_Init_Mesh_Check = .FALSE.
END IF

END FUNCTION Poseidon_Init_Mesh_Check



 !+501+####################################################!
!                                                           !
!          Poseidon_Init_Mesh_Check                         !
!                                                           !
 !#########################################################!
LOGICAL FUNCTION Poseidon_Init_Tables_Check()

IF ( ALL(lPF_Init_Tables_Flags) ) THEN
    Poseidon_Init_Tables_Check = .TRUE.
ELSE
    
    IF ( .NOT. lPF_Init_Mesh_Flags(iPF_Init_Tables_Alloc))         &
        CALL Warning_Message('Poseidon Initialization Check Failed. Mesh variables not allocated.')

    IF ( .NOT. lPF_Init_Mesh_Flags(iPF_Init_Tables_Init))         &
        CALL Warning_Message('Poseidon Initialization Check Failed. Mesh variables not initialized.')


    Poseidon_Init_Tables_Check = .FALSE.
END IF

END FUNCTION Poseidon_Init_Tables_Check





 !+601+####################################################!
!                                                           !
!          Poseidon_Init_MTGV_Check                   !
!                                                           !
 !#########################################################!
LOGICAL FUNCTION Poseidon_Init_MTGV_Check()

IF ( ALL(lPF_Init_MTGV_Flags) ) THEN
     Poseidon_Init_MTGV_Check = .TRUE.
ELSE

    IF ( .NOT. lPF_Init_MTGV_Flags(iPF_Init_MTGV_Derived))         &
        CALL Warning_Message('Poseidon Initialization Check Failed. Derived variables not set.')

    IF ( .NOT. lPF_Init_MTGV_Flags(iPF_Init_MTGV_TransMat))         &
        CALL Warning_Message('Poseidon Initialization Check Failed. Translation matrix not built.')


     Poseidon_Init_MTGV_Check = .FALSE.
END IF


END FUNCTION Poseidon_Init_MTGV_Check














END MODULE Flags_Initialization_Module


