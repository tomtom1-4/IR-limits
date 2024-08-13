#-
#: IncDir procedures
#: SmallExtension 200M
#: WorkSpace      4G
#: MaxTermSize 4M
Off Statistics;

#include declarations.h
Format 160;

auto i mui;
auto v pi;
v E1, E2, E3;
f J1, Gamma2, Gamma3;
f W1, W2, W3Dip, W3Quad;
cf S2, S3, S4, gamma2, gamma3, j1;

l J2 = 1/2*(J1(pi1, c1, q1, i1)*J1(pi2, c2, q2, i2) + J1(pi2, c2, q2, i2)*J1(pi1, c1, q1, i1))
  + Gamma2(pi1, c1, c2, q1, q2, i1);

l J3 = (1/6*(J1(pi1, c1, q1, i1)*J1(pi2, c2, q2, i2)*J1(pi3, c3, q3, i3)
           +J1(pi1, c1, q1, i1)*J1(pi3, c3, q3, i3)*J1(pi2, c2, q2, i2)
           +J1(pi2, c2, q2, i2)*J1(pi1, c1, q1, i1)*J1(pi3, c3, q3, i3)
           +J1(pi2, c2, q2, i2)*J1(pi3, c3, q3, i3)*J1(pi1, c1, q1, i1)
           +J1(pi3, c3, q3, i3)*J1(pi1, c1, q1, i1)*J1(pi2, c2, q2, i2)
           +J1(pi3, c3, q3, i3)*J1(pi2, c2, q2, i2)*J1(pi1, c1, q1, i1))
      +1/2*(J1(pi1, c1, q1, i1)*Gamma2(pi2, c2, c3, q2, q3, i2)
           +Gamma2(pi2, c2, c3, q2, q3, i2)*J1(pi1, c1, q1, i1)
           +J1(pi1, c2, q2, i1)*Gamma2(pi2, c1, c3, q1, q3, i2)
           +Gamma2(pi2, c1, c3, q1, q3, i2)*J1(pi1, c2, q2, i1)
           +J1(pi1, c3, q3, i1)*Gamma2(pi2, c1, c2, q1, q2, i2)
           +Gamma2(pi2, c1, c2, q1, q2, i2)*J1(pi1, c3, q3, i1))
      +Gamma3(pi1, c1, c2, c3, q1, q2, q3, i1));

#do i=1,3
#do j=1,3
  id J1(pi`i', c`j', q`j', i`i') = j1(pi`i', q`j', mu`j')*E`j'(mu`j')*T(i`i', c`j')*summe(i`i');
  sum cOli1, cOli2, mu`j';
  #do k=`j'+1,3
    id once Gamma2(pi2?, c`j', c`k', q`j', q`k', i2?) = summe(i2)*i_*cOlf(c`j', c`k', cOli1)*T(i2, cOli1)*gamma2(pi2, q`j', q`k', mu`j', mu`k')*E`j'(mu`j')*E`k'(mu`k');
    sum cOli1, mu`j', mu`k';
  #enddo
#enddo
#enddo

id Gamma3(pi1, c1, c2, c3, q1, q2, q3, i1) = summe(i1)*cOlf(c1, c2, cOli1)*cOlf(cOli1, c3, cOli2)*T(i1, cOli2)*gamma3(pi1, q1, q2, q3, mu1, mu2, mu3)*E1(mu1)*E2(mu2)*E3(mu3)
  *(1 + replace_(q3, q1, q1, q3, c3, c1, c1, c3, mu3, mu1, mu1, mu3, E3, E1, E1, E3) + replace_(q3, q2, q2, q3, c3, c2, c2, c3, mu3, mu2, mu2, mu3, E3, E2, E2, E3));
sum mu1, mu2, mu3;
*id Gamma3(pi1, c1, c2, c3, q1, q2, q3, i1) = summe(i1)*cOlf(c1, c2, cOli1)*cOlf(cOli1, c3, cOli2)*T(i1, cOli2)
*  *den(pi1.q1 + pi1.q2 + pi1.q3)*(1/12*pi1.E1*pi1.E2*pi1.E3*(3*pi1.q3 - pi1.q1 - pi1.q2)*den(pi1.q2*pi1.q3*(pi1.q1 + pi1.q2))
*    + pi1.E3*(pi.q3 - pi1.q1 - pi1.q2)*den(pi1.q3*(pi1.q1 + pi1.q2)*2*q1.q2)*(1/2*E1.E2*pi1.q1 + pi1.E2*q2.E1)
*    + den(2*(q1.q2 + q1.q3 + q2.q3)*2*q1.q2)*(2*q1.q2*pi1.E1*E2.E3 + 2*q2.E1*E2.E3*(pi1.q3 - pi1.q1 - pi1.q2) + 4*q3.E1*q1.E2*pi1.E3
*    + 4*q2.E1*pi1.E2*(q1.E3 + q2.E3) + E1.E2*(2*q2.q3*pi1.E3 + q1.E3*(pi1.q1 + pi1.q3 - 3*pi1.q2))))*(1 - replace_(E1, E2, E2, E1, q1, q2, q2, q1))
*    * (1 + replace_(E1, E3, E3, E1, q1, q3, q3, q1, c1, c3, c3, c1) + replace_(E2, E3, E3, E2, q2, q3, q3, q2, c2, c3, c3, c2));
sum cOli1, cOli2;
.sort

* Define complex conjugate of J3
* All expressions are symmetrized, i.e. the color operators don't need to be explicitly reversed (?)
l J2Conj = J2*replace_(i_, -i_, i1, i3, i2, i4, pi1, pi3, pi2, pi4);
l J3Conj = J3*replace_(i_, -i_, i1, i4, pi1, pi4, i2, i5, pi2, pi5, i3, i6, pi3, pi6);

* rename color indices
InExpression J3Conj, J2Conj;
  mul replace_(N1_?, cOli1, N2_?, cOli2, N3_?, cOli3, N4_?, cOli4, N5_?, cOli5, N6_?, cOli6, N7_?, cOli7, N8_?, cOli8, N9_?, cOli9, N10_?, cOli10);
EndInExpression;
.sort

l square = J3Conj*J3;
l square2 = J2Conj*J2;
sum cOli1,...,cOli8;
.sort
l Catani = (1/6*(W1(pi1, pi2, q1, i1, i2)*W1(pi3, pi4, q2, i3, i4)*W1(pi5, pi6, q3, i5, i6)
                +W1(pi1, pi2, q1, i1, i2)*W1(pi5, pi6, q3, i5, i6)*W1(pi3, pi4, q2, i3, i4)
                +W1(pi3, pi4, q2, i3, i4)*W1(pi1, pi2, q1, i1, i2)*W1(pi5, pi6, q3, i5, i6)
                +W1(pi3, pi4, q2, i3, i4)*W1(pi5, pi6, q3, i5, i6)*W1(pi1, pi2, q1, i1, i2)
                +W1(pi5, pi6, q3, i5, i6)*W1(pi1, pi2, q1, i1, i2)*W1(pi3, pi4, q2, i3, i4)
                +W1(pi5, pi6, q3, i5, i6)*W1(pi3, pi4, q2, i3, i4)*W1(pi1, pi2, q1, i1, i2))
           +1/2*(W1(pi1, pi2, q1, i1, i2)*W2(pi3, pi4, q2, q3, i3, i4)
                +W2(pi3, pi4, q2, q3, i3, i4)*W1(pi1, pi2, q1, i1, i2))
               *(1 + replace_(q1, q2, q2, q1) + replace_(q1, q3, q3, q1))
          +W3Dip(pi1, pi2, q1, q2, q3, i1, i2) + W3Quad(pi1, pi2, pi3, pi4, q1, q2, q3, i1, i2, i3, i4));

l Catani2 = 1/2*(W1(pi1, pi2, q1, i1, i2)*W1(pi3, pi4, q2, i3, i4) + W1(pi1, pi2, q2, i1, i2)*W1(pi3, pi4, q1, i3, i4))
  + W2(pi1, pi2, q1, q2, i1, i2);

.sort

repeat;
  id once W1(pi1?, pi2?, q1?, i1?, i2?) = -summe(i1)*summe(i2)*T(i1, cOli1)*T(i2, cOli1)*j1(pi1, q1, mu1)*j1(pi2, q1, mu1);
  sum cOli1, mu1;
  id once W2(pi1?, pi2?, q1?, q2?, i1?, i2?) = -CA*summe(i1)*summe(i2)*T(i1, cOli1)*T(i2, cOli1)*1/2
    *(S2(pi1, pi2, q1, q2) + S2(pi2, pi1, q1, q2) - 0*S2(pi1, pi1, q1, q2) - 0*S2(pi2, pi2, q1, q2));
  sum cOli1;
  id once W3Dip(pi1, pi2, q1, q2, q3, i1, i2) = -CA^2*summe(i1)*summe(i2)*T(i1, cOli1)*T(i2, cOli1)*S3(pi1, pi2, q1, q2, q3);
  sum cOli1;
  id once W3Quad(pi1, pi2, pi3, pi4, q1, q2, q3, i1, i2, i3, i4) = 1/2*cOlf(cOli1, cOli2, cOli5)*cOlf(cOli5, cOli3, cOli4)*summe(i1)*summe(i2)*summe(i3)*summe(i4)
    *(T(i3, cOli1)*(T(i1, cOli3)*T(i2, cOli4) + T(i2, cOli4)*T(i1, cOli3))*T(i4, cOli2) + T(i4, cOli2)*(T(i1, cOli3)*T(i2, cOli4) + T(i2, cOli4)*T(i1, cOli3))*T(i3, cOli1))
    *S4(pi1, pi2, pi3, pi4, q1, q2, q3);
  sum cOli1,...,cOli5;
endrepeat;
id S2(pi1?, pi2?, q1?, q2?) = -(+ gamma2(pi1,q1,q2,N2_?,N3_?)*gamma2(pi2,q1,q2,N2_?,N3_?)
          + gamma2(pi1,q1,q2,N2_?,N3_?)*j1(pi1,q1,N2_?)*j1(pi2,q2,N3_?)
          - gamma2(pi2,q1,q2,N2_?,N3_?)*j1(pi1,q1,N2_?)*j1(pi2,q2,N3_?)
          - 1/2*j1(pi1,q1,N2_?)^2*j1(pi1,q2,N3_?)*j1(pi2,q2,N3_?)
          + 3/4*j1(pi1,q1,N2_?)*j1(pi1,q2,N3_?)*j1(pi2,q1,N2_?)*j1(pi2,q2,N3_?)
          - 1/2*j1(pi1,q1,N2_?)*j1(pi2,q1,N2_?)*j1(pi2,q2,N3_?)^2);

*id S2(pi1?, pi2?, q1?, q2?) = (1 - (4 - d)/2)*den(q1.q2^2)*(pi1.q1*pi2.q2 + pi1.q2*pi2.q1)*den(pi1.q1 + pi1.q2)*den(pi2.q1 + pi2.q2)
*  - pi1.pi2^2/2*den(pi1.q1*pi2.q2*pi1.q2*pi2.q1)*(2 - (pi1.q1*pi2.q2 + pi1.q2*pi2.q1)*den(pi1.q1 + pi1.q2)*den(pi2.q1 + pi2.q2))
*  + pi1.pi2/2*den(q1.q2)*(2*den(pi1.q1*pi2.q2) + 2*den(pi2.q1*pi1.q2) - den(pi1.q1 + pi1.q2)*den(pi2.q1 + pi2.q2)*(4 + (pi1.q1*pi2.q2 + pi1.q2*pi2.q1)^2*den(pi1.q1*pi2.q2*pi1.q2*pi2.q1)));

id S3(pi1?, pi2?, q1?, q2?, q3?) = (1/2*gamma3(pi2, q1, q2, q3, mu1, mu2, mu3)*(gamma3(pi1, q1, q2, q3, mu1, mu2, mu3) + gamma3(pi1, q1, q3, q2, mu1, mu3, mu2) + gamma2(pi1, q1, q2, mu1, mu2)*j1(pi2, q3, mu3) - gamma2(pi2, q1, q2, mu1, mu2)*j1(pi1, q3, mu3)
  + gamma2(pi1, q1, q3, mu1, mu3)*j1(pi2, q2, mu2) - gamma2(pi2, q1, q3, mu1, mu3)*j1(pi1, q2, mu2) + 1/2*j1(pi2, q1, mu1)*j1(pi1, q2, mu2)*(j1(pi1, q3, mu3) + j1(pi2, q3, mu3)))
  + 1/2*gamma2(pi1, q1, q2, mu1, mu2)*(gamma2(pi2, q1, q2, mu1, mu2)*(3/4*j1(pi1, q3, mu3)*j1(pi2, q3, mu3) - 1/2*j1(pi1, q3, mu3)*j1(pi1, q3, mu3)) - 1/2*gamma2(pi1, q1, q2, mu1, mu2)*j1(pi1, q3, mu3)*j1(pi2, q3, mu3)
  + 1/4*gamma2(pi2, q1, q3, mu1, mu3)*(j1(pi2, q2, mu2)*j1(pi1, q3, mu3) + 2*j1(pi1, q2, mu2)*j1(pi2, q3, mu3) - 2*j1(pi1, q2, mu2)*j1(pi1, q3, mu3)) - 1/2*gamma2(pi1, q1, q3, mu1, mu3)*j1(pi2, q2, mu2)*j1(pi1, q3, mu3)
  + j1(pi1, q1, mu1)*j1(pi2, q2, mu2)*(7/4*j1(pi1, q3, mu3)*j1(pi2, q3, mu3) - 3/4*j1(pi1, q3, mu3)*j1(pi1, q3, mu3) - 1/2*j1(pi2, q3, mu3)*j1(pi2, q3, mu3)))
  + 1/2*j1(pi1, q1, mu1)*j1(pi1, q2, mu2)*j1(pi1, q3, mu3)*j1(pi2, q3, mu3)*(1/3*j1(pi1, q1, mu1)*j1(pi1, q2, mu2) - 5/6*j1(pi1, q1, mu1)*j1(pi2, q2, mu2) + 31/72*j1(pj, q1, mu1)*j1(pj, q2, mu2))
  + 1/16*j1(pi1, q1, mu1)^2*j1(pi1, q2, mu2)*j1(pi2, q2, mu2)*j1(pi2, q3, mu3)^2)
  *(1 + replace_(q1, q2, q2, q1, mu1, mu2, mu2, mu1)
      + replace_(q1, q3, q3, q1, mu1, mu3, mu3, mu1)
      + replace_(q2, q3, q3, q2, mu2, mu3, mu3, mu2)
      + replace_(q1, q2, q2, q3, q3, q1, mu1, mu2, mu2, mu3, mu3, mu1)
      + replace_(q1, q3, q3, q2, q2, q1, mu1, mu3, mu3, mu2, mu2, mu1));

id S4(pi1, pi2, pi3, pi4, q1, q2, q3) = ((1/2*gamma3(pi1, q1, q2, q3, mu1, mu2, mu3)*j1(pi3, q1, mu1) - 7/24*j1(pi1, q1, mu1)*j1(pi3, q2, mu2)*j1(pi2, q1, mu1)*j1(pi4, q3, mu3))*j1(pi4, q2, mu2)*j1(pi2, q3, mu3)
  - gamma2(pi1, q1, q2, mu1, mu2)*(1/2*j1(pi2, q1, mu1)*j1(pi3, q2, mu2)*j1(pi3, q3, mu3)*j1(pi4, q3, mu3) + 1/4*j1(pi4, q1, mu1)*j1(pi3, q2, mu2)*j1(pi1, q3, mu3)*j1(pi2, q3, mu3)
  + 1/2*gamma2(pi3, q1, q3, mu1, mu3)*j1(pi2, q2, mu2)*j1(pi4, q3, mu3)))
  *(1 + replace_(q1, q2, q2, q1, mu1, mu2, mu2, mu1)
      + replace_(q1, q3, q3, q1, mu1, mu3, mu3, mu1)
      + replace_(q2, q3, q3, q2, mu2, mu3, mu3, mu2)
      + replace_(q1, q2, q2, q3, q3, q1, mu1, mu2, mu2, mu3, mu3, mu1)
      + replace_(q1, q3, q3, q2, q2, q1, mu1, mu3, mu3, mu2, mu2, mu1));

*id gamma2(pi1?, q1?, q2?, mu1?, mu2?) = den(pi1.q1 + pi1.q2)*(pi1(mu1)*pi1(mu2)*(den(2*pi1.q1) - den(2*pi1.q2)) + den(q1.q2)*(pi1(mu1)*q1(mu2) - pi1(mu2)*q2(mu1) + 1/2*d_(mu1, mu2)*(pi1.q2 - pi1.q1)));

*id gamma3(pi1?, q1?, q2?, q3?, mu1?, mu2?, mu3?) = den(pi1.q1 + pi1.q2 + pi1.q3)*(1/12*pi1(mu1)*pi1(mu2)*pi1(mu3)*(3*pi1.q3 - pi1.q1 - pi1.q2)*den(pi1.q2*pi1.q3*(pi1.q1 + pi1.q2))
*  + pi1(mu3)*(pi1.q3 - pi1.q1 - pi1.q2)*den(pi1.q3*(pi1.q1 + pi1.q2)*2*q1.q2)*(1/2*d_(mu1, mu2)*pi1.q1 + pi1(mu2)*q2(mu1))
*  + den(2*(q1.q2 + q1.q3 + q2.q3)*2*(q1.q2))*(2*q1.q2*pi1(mu1)*d_(mu2, mu3) + 2*q2(mu1)*d_(mu2, mu3)*(pi1.q3 - pi1.q1 - pi1.q2) + 4*q3(mu1)*q1(mu2)*pi1(mu3)
*  + 4*q2(mu1)*pi1(mu2)*(q1(mu3) + q2(mu3)) + d_(mu1, mu2)*(2*q2.q3*pi1(mu3) + q1(mu3)*(pi1.q1 + pi1.q3 - 3*pi1.q2))))
*  *(1 - replace_(q1, q2, q2, q1, mu1, mu2, mu2, mu1));
*id j1(pi1?, q1?, mu1?) = pi1(mu1)*den(pi1.q1);
.sort
l difference = square - Catani;
l difference2 = square2 - Catani2;
.sort

* Use spin sums to replace Polarization vectors
ct EE1, EE2, EE3;
drop J3, J3Conj, J2, J2Conj;
toTensor functions E1, EE1;
toTensor functions E2, EE2;
toTensor functions E3, EE3;
id EE1(mu1?, mu2?) = -d_(mu1, mu2);
id EE2(mu1?, mu2?) = -d_(mu1, mu2);
id EE3(mu1?, mu2?) = -d_(mu1, mu2);
sum c1, c2, c3;
.sort
*if(count(d, 1) != 1) discard;
renumber 1;
.sort
sum cOli1,...,cOli8;
mul replace_(N1_?, cOli1, N2_?, cOli2, N3_?, cOli3, N4_?, cOli4, N5_?, cOli5, N6_?, cOli6, N7_?, cOli7, N8_?, cOli8, N9_?, cOli9, N10_?, cOli10);
* order color operators
#do l=1,4
* jacobi
  id cOlf(cOli1,cOli4,cOli5)*cOlf(cOli2,cOli3,cOli5) = - cOlf(cOli1,cOli2,cOli5)*cOlf(cOli3,cOli4,cOli5) - cOlf(cOli1,cOli3,cOli5)*cOlf(cOli4,cOli2,cOli5);
  repeat;
    #do i=1,6
      #do j=1,6
        id once T(i`i', c2?!{cOli1})*T(i`j', cOli1)*summe(i`i')*summe(i`j') = T(i`j', cOli1)*T(i`i', c2)*summe(i`i')*summe(i`j') + i_*cOlf(c2, cOli1, cOli12)*T(i`i', cOli12)*summe(i`i')*replace_(i`j',i`i',pi`j',pi`i');
        sum cOli12;
        id once T(i`i', c2?!{cOli1,cOli2})*T(i`j', cOli2)*summe(i`i')*summe(i`j') = T(i`j', cOli2)*T(i`i', c2)*summe(i`i')*summe(i`j') + i_*cOlf(c2, cOli2, cOli13)*T(i`i', cOli13)*summe(i`i')*replace_(i`j',i`i',pi`j',pi`i');
        sum cOli13;
        id once T(i`i', c2?!{cOli1,cOli2,cOli3})*T(i`j', cOli3)*summe(i`i')*summe(i`j') = T(i`j', cOli3)*T(i`i', c2)*summe(i`i')*summe(i`j') + i_*cOlf(c2, cOli3, cOli14)*T(i`i', cOli14)*summe(i`i')*replace_(i`j',i`i',pi`j',pi`i');
        sum cOli14;
      #enddo
    #enddo
  endrepeat;
  repeat;
    id cOlf(cOli1?, cOli2?, cOli3?)*cOlf(cOli1?, cOli2?, cOli4?) = CA*delta(cOli3, cOli4);
    id cOlf(cOli1?, cOli3?, cOli4?)*cOlf(cOli2?, cOli5?, cOli6?)*cOlf(cOli3?, cOli5?, cOli7?)*cOlf(cOli4?, cOli6?, cOli7?) = CA^2/2*delta(cOli1, cOli2);
    id cOlf(cOli1?,cOli4?,cOli5?)*cOlf(cOli2?,cOli4?,cOli6?)*cOlf(cOli3?,cOli5?,cOli6?) = - cOlTr(cOli1,cOli2,cOli3)*i_*NF + cOlTr(cOli1,cOli3,cOli2)*i_*NF;
    id T(i1?, cOli1?)*delta(cOli1?, cOli2?) = T(i1, cOli2);
    id cOlf(cOli1?, cOli2?, cOli3?)*delta(cOli3?, cOli4?) = cOlf(cOli1, cOli2, cOli4);
    id cOlTr(cOli1, cOli3, cOli2) = cOlTr(cOli1, cOli2, cOli3) + i_*cOlf(cOli3, cOli2, cOli1)*TF;
    id TF*NF = CA;
  endrepeat;
  .sort
  sum cOli1,...,cOli8;
  mul replace_(N1_?, cOli1, N2_?, cOli2, N3_?, cOli3, N4_?, cOli4, N5_?, cOli5, N6_?, cOli6, N7_?, cOli7, N8_?, cOli8, N9_?, cOli9, N10_?, cOli10);
* Find products of color operators
  #if `l'>=2
    #do i=2,4
      if(match(T(i2?, cOli`i')*T(i3?, cOli`i')) && (match(T(i4?, cOli1)*T(i5?, cOli1))==0)) mul replace_(cOli1, cOli`i', cOli`i', cOli1);
    #enddo
  #endif
#enddo
sum cOli1,...,cOli8;
.sort
* rename
#do l=1,2
#do i=1,5
  #do j=`i'+1,6
    if(match(summe(i`i'))==0);
      if(match(summe(i`j')));
        mul replace_(i`j', i`i', pi`j', pi`i');
        redefine j "10";
        redefine i "1";
      endif;
    endif;
    .sort
  #enddo
#enddo
#enddo
.sort
repeat;
#do i=1,5
  #do j=`i'+1, 6
    id once T(i`j', cOli1?)*T(i`i', cOli2?) = T(i`j', cOli1)*T(i`i', cOli2)*replace_(i`j',i`i',i`i',i`j',pi`j',pi`i',pi`i',pi`j');
  #enddo
#enddo
endrepeat;
.sort
mul (1 + replace_(q1, q2, q2, q1) + replace_(q1, q3, q3, q1) + replace_(q2, q3, q3, q2))/4;
id once pi1.pi3 = pi1.pi3*replace_(i3,i2,i2,i3,pi3,pi2,pi2,pi3);
id once pi2.pi3 = pi2.pi3*replace_(i2,i1,i1,i3,i3,i2,pi2,pi1,pi1,pi3,pi3,pi2);
id once pi1.pi2*pi1.pi3 = pi1.pi2*pi1.pi3*(1 + replace_(i3, i2, i2, i3, pi3, pi2, pi2, pi3))/2;
id once pi1.pi2 = pi1.pi2*(1 + replace_(i1, i2, i2, i1, pi1, pi2, pi2, pi1))/2;
.sort

sum cOli1,...,cOli8;
mul replace_(N1_?, cOli1, N2_?, cOli2, N3_?, cOli3, N4_?, cOli4, N5_?, cOli5, N6_?, cOli6, N7_?, cOli7, N8_?, cOli8, N9_?, cOli9, N10_?, cOli10);
#do i=1,3
repeat;
  id cOlf(cOli1,cOli4,cOli5)*cOlf(cOli2,cOli3,cOli5) = - cOlf(cOli1, cOli2, cOli5)*cOlf(cOli3, cOli4, cOli5) - cOlf(cOli1, cOli3, cOli5)*cOlf(cOli4, cOli2, cOli5);
  #do k=1,2
    #do j=`k'+1,6
      id once T(i`j', cOli1?)*T(i`k', cOli2?)*summe(i`j')*summe(i`k') = T(i`k', cOli2)*T(i`j', cOli1)*summe(i`k')*summe(i`j') + i_*cOlf(cOli1, cOli2, cOli3)*T(i`k', cOli3)*summe(i`k')*replace_(i`j', i`k', pi`j', pi`k');
      sum cOli3;
    #enddo
  #enddo
  repeat;
    id cOlf(cOli1?, cOli2?, cOli3?)*cOlf(cOli1?, cOli2?, cOli4?) = CA*delta(cOli3, cOli4);
    id cOlf(cOli1?, cOli3?, cOli4?)*cOlf(cOli2?, cOli5?, cOli6?)*cOlf(cOli3?, cOli5?, cOli7?)*cOlf(cOli4?, cOli6?, cOli7?) = CA^2/2*delta(cOli1, cOli2);
    id cOlf(cOli1?,cOli4?,cOli5?)*cOlf(cOli2?,cOli4?,cOli6?)*cOlf(cOli3?,cOli5?,cOli6?) = - cOlTr(cOli1,cOli2,cOli3)*i_*NF + cOlTr(cOli1,cOli3,cOli2)*i_*NF;
    id T(i1?, cOli1?)*delta(cOli1?, cOli2?) = T(i1, cOli2);
    id cOlf(cOli1?, cOli2?, cOli3?)*delta(cOli3?, cOli4?) = cOlf(cOli1, cOli2, cOli4);
    id cOlTr(cOli1, cOli3, cOli2) = cOlTr(cOli1, cOli2, cOli3) + i_*cOlf(cOli3, cOli2, cOli1)*TF;
    id TF*NF = CA;
  endrepeat;
endrepeat;
.sort
  id once gamma2(pi1,q1,q2,mu1?,mu2?)*j1(pi2,q1,mu1?)*j1(pi3,q2,mu2?)*summe(i1)*summe(i2)*summe(i3) = gamma2(pi1,q1,q2,mu1,mu2)*j1(pi2,q1,mu1)*j1(pi3,q2,mu2)
                                                                                                  *summe(i1)*summe(i2)*summe(i3)*replace_(i2, i1, i1, i3, i3, i2, pi2, pi1, pi1, pi3, pi3, pi2);
  id once gamma2(pi1,q1,q2,mu1?,mu2?)*j1(pi2,q2,mu2?)*j1(pi3,q1,mu1?) = gamma2(pi1,q1,q2,mu1,mu2)*j1(pi2,q2,mu2)*j1(pi3,q1,mu1)*replace_(i3, i1, i1, i3, pi3, pi1, pi1, pi3);
  id once j1(pi1,q2,mu1?)*j1(pi2,q1,mu2?) = j1(pi1,q2,mu1)*j1(pi2,q1,mu2)*replace_(i1, i2, i2, i1, pi1, pi2, pi2, pi1);
.sort
  sum cOli1,...,cOli8;
  mul replace_(N1_?, cOli1, N2_?, cOli2, N3_?, cOli3, N4_?, cOli4, N5_?, cOli5, N6_?, cOli6, N7_?, cOli7, N8_?, cOli8, N9_?, cOli9, N10_?, cOli10);
#enddo
sum cOli1,...,cOli8;

#do l=1,2
#do i=1,5
  #do j=`i'+1,6
    if(match(summe(i`i'))==0);
      if(match(summe(i`j')));
        mul replace_(i`j', i`i', pi`j', pi`i');
        redefine j "10";
        redefine i "1";
      endif;
    endif;
    .sort
  #enddo
#enddo
#enddo
.sort
repeat;
#do i=1,5
  #do j=`i'+1, 6
    id once T(i`j', cOli1?)*T(i`i', cOli2?) = T(i`j', cOli1)*T(i`i', cOli2)*replace_(i`j',i`i',i`i',i`j',pi`j',pi`i',pi`i',pi`j');
  #enddo
#enddo
endrepeat;

.sort
#do i=1,6
  id pi`i'.q1 = si`i'1/2;
  id pi`i'.q2 = si`i'2/2;
  id pi`i'.q3 = si`i'3/2;
  id pi`i'.pi`i' = 0;
  argument den;
    id pi`i'.q1 = si`i'1/2;
    id pi`i'.q2 = si`i'2/2;
    id pi`i'.q3 = si`i'3/2;
    id pi`i'.pi`i' = 0;
  endargument;
#enddo
#do i=1,5
  #do j=`i'+1,6
    id pi`i'.pi`j' = si`i'i`j'/2;
    argument den;
      id pi`i'.pi`j' = si`i'i`j'/2;
    endargument;
  #enddo
#enddo
#do i=1,3
  id q`i'.q`i' = 0;
  argument den;
    id q`i'.q`i' = 0;
  endargument;
  #do j=`i'+1,3
    id q`i'.q`j' = s`i'`j'/2;
    argument den;
      id q`i'.q`j' = s`i'`j'/2;
    endargument;
  #enddo
#enddo
.sort
if((count(T, 1) != 5)) discard;
*if(count(cOlf, 1) != 0) discard;
#call FullSimplify
*if(count(d, 1) != 1) discard;
#call Simplify

b cOlf, T, delta, d, CA, summe, cOlTr, si1i2, si1i3,si2i3;
print+s square, Catani, difference;
.end

