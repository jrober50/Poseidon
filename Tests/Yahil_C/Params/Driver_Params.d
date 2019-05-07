RE               100                   RE	   ! Number of Radial Elements
CE                 1                   CE          ! Number of Radial Elements in Core
TE                 1                   TE          ! Number of Theta Elements
PE                 1                   PE          ! Number of Phi Elements

RQ                 1                   RQ          ! Number of Radial Source Points Per Element
TQ                 1                   TQ          ! Number of Theta Source Points Per Element
PQ                 1                   PQ          ! Number of Phi Source Points Per Element

DIM                3                   DIM         ! Number of Dimensions

PROC               1                   PROC        ! Total Number of Processes, PROC = yPROC*zPROC
yPROC              1                   yPROC       ! Number of Processes in Theta Direction
zPROC              1                   zPROC       ! Number of Processes in Phi Direction

MT                 3                   MT          ! Mesh Type, 1 = Uniform, 2 = Logarithmic, 3 = Split 
IR          0.00E+00                   IR          ! Inner Boundary Distance
CR          3.00E+06                   CR          ! Core Boundary Distance 
OR          1.00E+09                   OR          ! Outer Boundary Distance

SST          5.0E+01                   SST         ! Self-Similar Time Parameter
SSK        9.539E+14                   SSK         ! Self-Similar Kappa Parameter
SSG          1.3E+00                   SSG         ! Self-Similar Gamma Parameter
