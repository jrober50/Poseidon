DIM                3                   DIM        ! Dimensions
DEGREE             1                   DEGREE     ! FEM Degree
LLIMIT             0                   LLIMIT     ! Spherical Harmonic Expansion Limit

PPROC              1                   PPROC      ! Number of Processes

NSHELL             1                   NSHELL     ! Number of Shells
NSSHEL             1                   NSSHEL     ! Number of Subshells

NBPSHL             1                   NBPSHL     ! Number of Blocks Per Shell, NBPSHL = NBTROW*NBPCOL
NBTROW             1                   NBTROW     ! Number of Theta Blocks Per Shell
NBPCOL             1                   NBPCOL     ! Number of Phi Blocks Per Shell

NREPS             50                   NREPS      ! Number of Radial Elements Per Shell
NREPSS            50                   NREPSS     ! Number of Radial Elements Per Subshell

NTEPB              1                   NTEPB      ! Number of Theta Elements Per Block
NPEPB              1                   NPEPB      ! Number of Phi Elements Per Block


PRQ                5                   PRQ        ! Number of Radial Quadrature Points Per Element
PTQ                3                   PTQ        ! Number of Theta Quadrature Points Per Element
PPQ                1                   PPQ        ! Number of Phi Quadrature Points Per Elements

RCF                1                   RCF        ! Radial Coarsening Factor
TCF                1                   TCF        ! Theta Coarsening Factor
PCF                1                   PCF        ! Phi Coarsening Factor

MI                 5                   MI         ! Maximum Newton-Raphson Iterations
CC          1.00E-12                   CC         ! Convergence Criteria
