#procedure FeynmanRules
********************************************************************************
* Notation for external parton polarization vectors and bispinors
* Standard color indices for external particles
* Light-cone-gauge Feynman rules for gluons
* Traces for closed fermion loops

****************************************
* gluons

id EpsStar(mu?,p3,cOli1?) = E3(mu)*d_(cOli1,c3);
id EpsStar(mu?,p4,cOli1?) = E4(mu) * d_(cOli1, c4);
id EpsStar(mu?,p5,cOli1?) = E5(mu) * d_(cOli1, c5);


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
  id once VVV(k1?,v1?,k2?,v2?,k3?,v3?) =
    (k1(v3)-k2(v3))*d_(v1,v2)+(k2(v1)-k3(v1))*d_(v2,v3)+(k3(v2)-k1(v2))*d_(v1,v3);
  id once VVV(0,v1?,k2?,v2?,k3?,v3?) =
    (-k2(v3))*d_(v1,v2)+(k2(v1)-k3(v1))*d_(v2,v3)+(k3(v2))*d_(v1,v3);
  id once VVV(k1?,v1?,0,v2?,k3?,v3?) =
    (k1(v3))*d_(v1,v2)+(-k3(v1))*d_(v2,v3)+(k3(v2)-k1(v2))*d_(v1,v3);
  id once VVV(k1?,v1?,k2?,v2?,0,v3?) =
    (k1(v3)-k2(v3))*d_(v1,v2)+(k2(v1))*d_(v2,v3)+(-k1(v2))*d_(v1,v3);
  if (count(VVV,1) != 0) redefine i "0";
  b VVV;
  .sort:VVV;
  Keep Brackets;
#enddo

.sort

****************************************
* quarks
id UfBar(p3,0,cOli3?) = UBar(p3)*d_(cOli3, c3);
id Vf(p4, 0, cOli4?) = V(p4)*d_(cOli4, c4);

id DS(p?, 0) = DS(p);

id SF(v?,k?,0) = g_(v,k)*DS(k);
id G(v1?,v2?) = g_(v1,v2);
.sort


****************************************
* Wilson lines

id EikV(p?, mu?) = -2 * p(mu);
id VWilson(p1, cOli1?) = d_(cOli1,c1);
id VWilson(p2, cOli1?) = d_(cOli1,c2);
.sort

id Quad = 1;


.sort
#endprocedure