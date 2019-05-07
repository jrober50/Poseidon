DIM                3                   DIM        ! Dimensions
DEGREE             1                   DEGREE     ! FEM Degree
LLIMIT             0                   LLIMIT     ! Spherical Harmonic Expansion Limit

PPROC              1                   PPROC      ! Number of Processes

NSHELL             1                   NSHELL     ! Number of Shells
NSSHEL             1                   NSSHEL     ! Number of Subshells

NBPSHL             1                   NBPSHL     ! Number of Blocks Per Shell, NBPSHL = NBTROW*NBPCOL
NBTROW             1                   NBTROW     ! Number of Theta Blocks Per Shell
NBPCOL             1                   NBPCOL     ! Number of Phi Blocks Per Shell

NREPS            700                   NREPS      ! Number of Radial Elements Per Shell
NREPSS           700                   NREPSS     ! Number of Radial Elements Per Subshell

NTEPB              1                   NTEPB      ! Number of Theta Elements Per Block
NPEPB              1                   NPEPB      ! Number of Phi Elements Per Block

PRQ                5                   PRQ        ! Number of Radial Quadrature Points Per Element
PTQ                3                   PTQ        ! Number of Theta Quadrature Points Per Element
PPQ                1                   PPQ        ! Number of Phi Quadrature Points Per Elements

MI                 1                   MI         ! Maximum Newton-Raphson Iterations
CC          1.00E-10                   CC         ! Convergence Criteria


WRTTT              3                   WRTTT      ! Write Timetable Flag        : 0=Off, 1=To Screen, 2=To File, 3=Both
WRTIR              3                   WRTIT      ! Write Iteration Report Flag : 0=Off, 1=To Screen, 2=To File, 3=Both
IRNS              20                   IRNS       ! Number of Samples in each Iteration Report
WRTRS              1                   WRTRS      ! Write Results to File Flag  : 0=Off, 1=To File
