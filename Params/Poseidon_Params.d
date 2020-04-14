DIM                3                   DIM        ! Dimensions
DEGREE             1                   DEGREE     ! FEM Degree
LLIMIT             1                   LLIMIT     ! Spherical Harmonic Expansion Limit

NREPS             10                   NREPS      ! Number of Radial Elements Per Shell
NREPSS            10                   NREPSS     ! Number of Radial Elements Per SubShell

NTEPB              1                   NTEPB      ! Number of Theta Elements Per Block
NPEPB              1                   NPEPB      ! Number of Phi Elements Per Block

PPROC              1                   PPROC      ! Number of Processes

NSHELL             1                   NSHELL     ! Number of Shells
NSSHEL             1                   NSSHEL     ! Number of Subshells

NBPSHL             1                   NBPSHL     ! Number of Blocks Per Shell, NBPSHL = NBTROW*NBPCOL
NBTROW             1                   NBTROW     ! Number of Theta Blocks Per Shell
NBPCOL             1                   NBPCOL     ! Number of Phi Blocks Per Shell

PRQ                8                   PRQ        ! Number of Radial Quadrature Points Per Element
PTQ                8                   PTQ        ! Number of Theta Quadrature Points Per Element
PPQ                8                   PPQ        ! Number of Phi Quadrature Points Per Elements

MI                10                   MI         ! Maximum Newton-Raphson Iterations
CC           1.00E-8                   CC         ! Convergence Criteria

OSTF               0                   OSTF       ! Output Setup Table Flag    :  0=Off(Default), 1=On
OMF                1                   OMF        ! Write Jacobian Matrix to File
ORF                1                   ORF        ! Write RHS Vector to File
WRTTT              0                   WRTTT      ! Write Timetable Flag        : 0=Off, 1=To Screen, 2=To File, 3=Both
WRTIR              1                   WRTIT      ! Write Iteration Report Flag : 0=Off, 1=To Screen, 2=To File, 3=Both
IRNS              20                   IRNS       ! Number of Samples in each Iteration Report
WRTRS              1                   WRTRS      ! Write Results to File Flag  : 0=Off, 1=To File

RSMPS           1000                   RSMPS      ! Number of Radial Samples for Results Output
TSMPS              1                   TSMPS      ! Number of Theta Samples for Results Output
PSMPS              1                   PSMPS      ! Number of Phi Samples for Results Output


NPS                0                   NPS        ! Use the New PETSc SNES routines

