#procedure FeynmanRules
********************************************************************************
* Notation for external parton polarization vectors and bispinors
* Standard color indices for external particles
* Light-cone-gauge Feynman rules for gluons
* Traces for closed fermion loops


****************************************
* gluons

id EpsStar(mu?,q,cOli1?) = E(mu);

b DV;

.sort
Keep Brackets;
#do i=1,1
  id once DV(v1?,v2?,k?,0) = DS(k)*(-d_(v1,v2));*+(k(v1)*n(v2)+k(v2)*n(v1))*Eik(k, n));
  if (count(DV,1) != 0) redefine i "0";
  b DV;
  .sort:DV;
  Keep Brackets;
#enddo

b VVV;
.sort
Keep Brackets;

#do i=1,1
  id once VVV(k1?,v1?,cOli1?,k2?,v2?,cOli2?,k3?,v3?,cOli3?) =
    ((k1(v3)-k2(v3))*d_(v1,v2)+(k2(v1)-k3(v1))*d_(v2,v3)+(k3(v2)-k1(v2))*d_(v1,v3))
    *(-i_*cOlf(cOli1, cOli2, cOli3));
  if (count(VVV,1) != 0) redefine i "0";
  b VVV;
  .sort:VVV;
  Keep Brackets;
#enddo
.sort

#do i=1,1
  id once VVVV(k1?,v1?,cOli1?,k2?,v2?,cOli2?,k3?,v3?,cOli3?, k4?, v4?, cOli4?) =
     (cOlf(cOli1, cOli2, cOli10)*cOlf(cOli3, cOli4, cOli10)*(d_(v1, v3)*d_(v2, v4) - d_(v1, v4)*d_(v2, v3))
     +cOlf(cOli1, cOli3, cOli10)*cOlf(cOli2, cOli4, cOli10)*(d_(v1, v2)*d_(v3, v4) - d_(v1, v4)*d_(v2, v3))
     +cOlf(cOli1, cOli4, cOli10)*cOlf(cOli2, cOli3, cOli10)*(d_(v1, v2)*d_(v3, v4) - d_(v1, v3)*d_(v2, v4)));
  sum cOli10;
  if (count(VVV,1) != 0) redefine i "0";
  b VVV;
  .sort:VVV;
  Keep Brackets;
#enddo

.sort

****************************************
* quarks
id DS(p?, 0) = DS(p);


.sort
id EikV(p?, mu?) = -2 * p(mu);
.sort


mul nl^`nl';
mul 1/`sym';

.sort
#endprocedure