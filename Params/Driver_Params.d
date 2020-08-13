RE               256                   RE	   ! Number of Radial Elements
CE               250                   CE          ! Number of Radial Elements in Core
TE                 1                   TE          ! Number of Theta Elements
PE                 1                   PE          ! Number of Phi Elements

RQ                 2                   RQ          ! Number of Radial Source Points Per Element
TQ                 1                   TQ          ! Number of Theta Source Points Per Element
PQ                 1                   PQ          ! Number of Phi Source Points Per Element

DIM                3                   DIM         ! Number of Dimensions

PROC               1                   PROC        ! Total Number of Processes, PROC = yPROC*zPROC
yPROC              1                   yPROC       ! Number of Processes in Theta Direction
zPROC              1                   zPROC       ! Number of Processes in Phi Direction

LL         -0.50E+00                   LL          ! Left Limit for Source x-Space
RL          0.50E+00                   RL          ! Right Limit for Source x-Space

DTN                3                   DTN         ! Test Number, 1 = Sphere Sym, 2 = CHIMERA, 3 = Yahil

SOF                0                   SOF         ! Source Output Flag,  1 = On, 0 = Off
ROF                1                   ROF         ! Results Output Flag, 1 = On, 0 = Off
RRF                1                   RRF         ! Run Report Flag,     1 = On, 0 = Off
FRF                1                   FRF         ! Frame Report Flag,   1 = On, 0 = Off


MT                 4                   MT          ! Mesh Type, 1 = Uniform, 2 = Logarithmic, 3 = Split(1&2), 4 = Geometric 
IR          0.00E+00                   IR          ! Inner Boundary Distance
CR          1.00E+05                   CR          ! Core Boundary Distance 
OR          1.00E+10                   OR          ! Outer Boundary Distance
GZ    1.032034864238313                GZ 

FGF                1                   FGF         ! First Guess Flag, 1=FlatSpace, 2=Newtonian Approx. (Test only), 3=Load From File
SGF                4                   SGF         ! Subsequent Guess Flag, 1,2,3 same as above, 4 = Use Previous Frames Solution


! Used if DTN == 2 !

CSF               61                   CSF          ! CHIMERA Start Frame
CEF               70                   CEF          ! CHIMERA End Frame


! Used if DTN == 3 !

YST          1.5E+02                   SST         ! Self-Similar Start Time Parameter, milliseconds
YET          5.0E-01                   YET         ! Self-Similar Final Time Parameter, milliseconds
YNF                1                   YNF         ! Self-Similar Number of Frames
                                                   ! If one frame, uses YST

SSK        9.539E+14                   SSK         ! Self-Similar Kappa Parameter
SSG          1.3E+00                   SSG         ! Self-Similar Gamma Parameter
SSE          0.00000                   SSE         ! Self-Similar Eccentricity Parameter

SSV                1                   SSV         ! Self-Similar Velocity on/off (On == 1)

! Used if DTN == 1 !

DEN          1.0E+10                   DEN         ! Base Density for Spherical Test Profile
PWA               -1                   PWA         ! Density's Radial Dependence, rho = DEN*r^PWA
