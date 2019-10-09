RE              2000                   RE	   ! Number of Radial Elements
CE               400                   CE          ! Number of Radial Elements in Core
TE                 1                   TE          ! Number of Theta Elements
PE                 1                   PE          ! Number of Phi Elements

RQ                 5                   RQ          ! Number of Radial Source Points Per Element
TQ                 1                   TQ          ! Number of Theta Source Points Per Element
PQ                 1                   PQ          ! Number of Phi Source Points Per Element

DIM                3                   DIM         ! Number of Dimensions

PROC               1                   PROC        ! Total Number of Processes, PROC = yPROC*zPROC
yPROC              1                   yPROC       ! Number of Processes in Theta Direction
zPROC              1                   zPROC       ! Number of Processes in Phi Direction

LL         -0.50E+00                   LL          ! Left Limit for Source x-Space
RL          0.50E+00                   RL          ! Right Limit for Source x-Space

DTN                3                   DTN         ! Test Number, 1 = Sphere Sym, 2 = CHIMERA, 3 = Yahil

SOF                1                   SOF         ! Source Output Flag,  1 = On, 0 = Off
ROF                1                   ROF         ! Results Output Flag, 1 = On, 0 = Off
RRF                1                   RRF         ! Run Report Flag,     1 = On, 0 = Off
FRF                1                   FRF         ! Frame Report Flag,   1 = On, 0 = Off


MT                 3                   MT          ! Mesh Type, 1 = Uniform, 2 = Logarithmic, 3 = Split 
IR          0.00E+00                   IR          ! Inner Boundary Distance
CR          1.00E+06                   CR          ! Core Boundary Distance 
OR          1.00E+09                   OR          ! Outer Boundary Distance



! Used if DTN == 2 !

CSF                1                  CSF          ! CHIMERA Start Frame
CEF                1                  CEF          ! CHIMERA End Frame


! Used if DTN == 3 !

YST          5.0E-01                   SST         ! Self-Similar Start Time Parameter
YET          5.0E-01                   YET         ! Self-Similar Final Time Parameter
YNF                1                   YNF         ! Self-Similar Number of Frames
                                                   ! If one frame, uses YST

SSK        9.539E+14                   SSK         ! Self-Similar Kappa Parameter
SSG          1.3E+00                   SSG         ! Self-Similar Gamma Parameter
SSE          0.00000                   SSE         ! Self-Similar Eccentricity Parameter


! Used if DTN == 1 !

DEN          1.0E+10                   DEN         ! Base Density for Spherical Test Profile
PWA               -1                   PWA         ! Density's Radial Dependence, rho = DEN*r^PWA
