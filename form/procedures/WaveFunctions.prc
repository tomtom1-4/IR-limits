#procedure WaveFunctions
********************************************************************************
* Transversality for polarisation vectors
* Equations-of-motion for spinors
* Additional Dirac algebra

****************************************

****************************************
* polarisation vectors
* assumed transverse to q


id p5.E5 = 0;
id p4.E4 = 0;
id p3.E3 = 0;
id UBar(p3) * g(p3) = 0;
id g(p4) * V(p4) = 0;


.sort

#call Kinematics

#endprocedure