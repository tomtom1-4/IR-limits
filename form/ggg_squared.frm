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
v q4, q5; * just a dummy vector
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
id once S2(pi1?, pi2?, q1?, q2?) = -(+ gamma2(pi1,q1,q2,mu1,mu2)*gamma2(pi2,q1,q2,mu1,mu2)
          + gamma2(pi1,q1,q2,mu1,mu2)*j1(pi1,q1,mu1)*j1(pi2,q2,mu2)
          - gamma2(pi2,q1,q2,mu1,mu2)*j1(pi1,q1,mu1)*j1(pi2,q2,mu2)
          - 1/2*j1(pi1,q1,mu1)^2*j1(pi1,q2,mu2)*j1(pi2,q2,mu2)
          + 3/4*j1(pi1,q1,mu1)*j1(pi1,q2,mu2)*j1(pi2,q1,mu1)*j1(pi2,q2,mu2)
          - 1/2*j1(pi1,q1,mu1)*j1(pi2,q1,mu1)*j1(pi2,q2,mu2)^2);
sum mu1, mu2;

id S3(pi1?, pi2?, q1?, q2?, q3?) = (1/2*gamma3(pi2, q1, q2, q3, mu1, mu2, mu3)*(gamma3(pi1, q1, q2, q3, mu1, mu2, mu3) + gamma3(pi1, q1, q3, q2, mu1, mu3, mu2) + gamma2(pi1, q1, q2, mu1, mu2)*j1(pi2, q3, mu3) - gamma2(pi2, q1, q2, mu1, mu2)*j1(pi1, q3, mu3)
  + gamma2(pi1, q1, q3, mu1, mu3)*j1(pi2, q2, mu2) - gamma2(pi2, q1, q3, mu1, mu3)*j1(pi1, q2, mu2) + 1/2*j1(pi2, q1, mu1)*j1(pi1, q2, mu2)*(j1(pi1, q3, mu3) + j1(pi2, q3, mu3)))
  + 1/2*gamma2(pi1, q1, q2, mu1, mu2)*(gamma2(pi2, q1, q2, mu1, mu2)*(3/4*j1(pi1, q3, mu3)*j1(pi2, q3, mu3) - 1/2*j1(pi1, q3, mu3)*j1(pi1, q3, mu3)) - 1/2*gamma2(pi1, q1, q2, mu1, mu2)*j1(pi1, q3, mu3)*j1(pi2, q3, mu3)
  + 1/4*gamma2(pi2, q1, q3, mu1, mu3)*(j1(pi2, q2, mu2)*j1(pi1, q3, mu3) + 2*j1(pi1, q2, mu2)*j1(pi2, q3, mu3) - 2*j1(pi1, q2, mu2)*j1(pi1, q3, mu3)) - 1/2*gamma2(pi1, q1, q3, mu1, mu3)*j1(pi2, q2, mu2)*j1(pi1, q3, mu3)
  + j1(pi1, q1, mu1)*j1(pi2, q2, mu2)*(7/4*j1(pi1, q3, mu3)*j1(pi2, q3, mu3) - 3/4*j1(pi1, q3, mu3)*j1(pi1, q3, mu3) - 1/2*j1(pi2, q3, mu3)*j1(pi2, q3, mu3)))
  + 1/2*j1(pi1, q1, mu1)*j1(pi1, q2, mu2)*j1(pi1, q3, mu3)*j1(pi2, q3, mu3)*(1/3*j1(pi1, q1, mu1)*j1(pi1, q2, mu2) - 5/6*j1(pi1, q1, mu1)*j1(pi2, q2, mu2) + 31/72*j1(pi2, q1, mu1)*j1(pi2, q2, mu2))
  + 1/16*j1(pi1, q1, mu1)^2*j1(pi1, q2, mu2)*j1(pi2, q2, mu2)*j1(pi2, q3, mu3)^2)
  *(1 + replace_(q1, q2, q2, q1, mu1, mu2, mu2, mu1)
      + replace_(q1, q3, q3, q1, mu1, mu3, mu3, mu1)
      + replace_(q2, q3, q3, q2, mu2, mu3, mu3, mu2)
      + replace_(q1, q2, q2, q3, q3, q1, mu1, mu2, mu2, mu3, mu3, mu1)
      + replace_(q1, q3, q3, q2, q2, q1, mu1, mu3, mu3, mu2, mu2, mu1));
sum mu1,...,mu3;

id S4(pi1, pi2, pi3, pi4, q1, q2, q3) = ((1/2*gamma3(pi1, q1, q2, q3, mu1, mu2, mu3)*j1(pi3, q1, mu1) - 7/24*j1(pi1, q1, mu1)*j1(pi3, q2, mu2)*j1(pi2, q1, mu1)*j1(pi4, q3, mu3))*j1(pi4, q2, mu2)*j1(pi2, q3, mu3)
  - gamma2(pi1, q1, q2, mu1, mu2)*(1/2*j1(pi2, q1, mu1)*j1(pi3, q2, mu2)*j1(pi3, q3, mu3)*j1(pi4, q3, mu3) + 1/4*j1(pi4, q1, mu1)*j1(pi3, q2, mu2)*j1(pi1, q3, mu3)*j1(pi2, q3, mu3)
  + 1/2*gamma2(pi3, q1, q3, mu1, mu3)*j1(pi2, q2, mu2)*j1(pi4, q3, mu3)))
  *(1 + replace_(q1, q2, q2, q1, mu1, mu2, mu2, mu1)
      + replace_(q1, q3, q3, q1, mu1, mu3, mu3, mu1)
      + replace_(q2, q3, q3, q2, mu2, mu3, mu3, mu2)
      + replace_(q1, q2, q2, q3, q3, q1, mu1, mu2, mu2, mu3, mu3, mu1)
      + replace_(q1, q3, q3, q2, q2, q1, mu1, mu3, mu3, mu2, mu2, mu1));
sum mu1,...,mu3;
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
renumber 1;
.sort
sum cOli1,...,cOli8;
sum mu1,...,mu4;
*
#do i=1,2
  #do j=`i',3
    id gamma2(pi1?, q`j', q`i', mu1?, mu2?) = -gamma2(pi1, q`i', q`j', mu2, mu1);
    id gamma3(pi1?, q`j', q`i', q4?, mu1?, mu2?, mu3?) = -gamma3(pi1, q`i', q`j', q4, mu2, mu1, mu3);
  #enddo
#enddo

#do l=1,5
  #do i=1,5
    #do j=`i'+1,6
      id once j1(pi`i', q2, mu1?)*j1(pi`j', q1, mu2?) = j1(pi`i', q2, mu1)*j1(pi`j', q1, mu2)*replace_(pi`i', pi`j', pi`j', pi`i', i`i', i`j', i`j', i`i');
      id once j1(pi`i', q3, mu1?)*j1(pi`j', q1, mu2?) = j1(pi`i', q3, mu1)*j1(pi`j', q1, mu2)*replace_(pi`i', pi`j', pi`j', pi`i', i`i', i`j', i`j', i`i');
      id once j1(pi`i', q3, mu1?)*j1(pi`j', q2, mu2?) = j1(pi`i', q3, mu1)*j1(pi`j', q2, mu2)*replace_(pi`i', pi`j', pi`j', pi`i', i`i', i`j', i`j', i`i');
    #enddo
  #enddo
  #do i=1,5
    #do j=`i'+1,6
      id once gamma2(pi`i', q1?, q2?, mu1?, mu2?)*j1(pi`j', q3?, mu3?) = gamma2(pi`i', q1, q2, mu1, mu2)*j1(pi`j', q3, mu3)*replace_(i`i', i`j', i`j', i`i', pi`i', pi`j', pi`j', pi`i');
      id once gamma3(pi`i', q1?, q2?, q3?, mu1?, mu2?, mu3?)*j1(pi`j', q4?, mu4?) = gamma3(pi`i', q1, q2, q3, mu1, mu2, mu3)*j1(pi`j', q4, mu4)*replace_(i`i', i`j', i`j', i`i', pi`i', pi`j', pi`j', pi`i');
      id once gamma3(pi`i', q1?, q2?, q3?, mu1?, mu2?, mu3?)*gamma2(pi`j', q4?, q5?, mu4?, mu5?) = gamma3(pi`i', q1, q2, q3, mu1, mu2, mu3)*gamma2(pi`j', q4, q5, mu4, mu5)*replace_(i`i', i`j', i`j', i`i', pi`i', pi`j', pi`j', pi`i');
    #enddo
  #enddo
  #do i=1,5
    #do j=`i'+1,6
      id once gamma2(pi`i', q2, q3, mu1?, mu2?)*gamma2(pi`j', q1, q2, mu3?, mu4?) = gamma2(pi`i', q2, q3, mu1, mu2)*gamma2(pi`j', q1, q2, mu3, mu4)*replace_(i`i', i`j', i`j', i`i', pi`i', pi`j', pi`j', pi`i');
      id once gamma2(pi`i', q1, q3, mu1?, mu2?)*gamma2(pi`j', q1, q2, mu3?, mu4?) = gamma2(pi`i', q1, q3, mu1, mu2)*gamma2(pi`j', q1, q2, mu3, mu4)*replace_(i`i', i`j', i`j', i`i', pi`i', pi`j', pi`j', pi`i');
      id once gamma2(pi`i', q2, q3, mu1?, mu2?)*gamma2(pi`j', q1, q3, mu3?, mu4?) = gamma2(pi`i', q2, q3, mu1, mu2)*gamma2(pi`j', q1, q3, mu3, mu4)*replace_(i`i', i`j', i`j', i`i', pi`i', pi`j', pi`j', pi`i');
    #enddo
  #enddo
* jacobi
  id cOlf(N1_?,N4_?,N5_?)*cOlf(N2_?,N3_?,N5_?) = - cOlf(N1_?,N2_?,N5_?)*cOlf(N3_?,N4_?,N5_?) - cOlf(N1_?,N3_?,N5_?)*cOlf(N4_?,N2_?,N5_?);
  repeat;
    #do k=1,5
      #do j=`k'+1,6
        id once T(i`j', cOli1?)*T(i`k', cOli2?)*summe(i`j')*summe(i`k') = T(i`k', cOli2)*T(i`j', cOli1)*summe(i`k')*summe(i`j') + i_*cOlf(cOli1, cOli2, cOli3)*T(i`k', cOli3)*summe(i`k')*replace_(i`j', i`k', pi`j', pi`k');
        sum cOli3;
      #enddo
  #enddo
  endrepeat;
  repeat;
    id cOlf(cOli1?, cOli2?, cOli3?)*cOlf(cOli1?, cOli2?, cOli4?) = CA*delta(cOli3, cOli4);
    id cOlf(cOli1?, cOli3?, cOli4?)*cOlf(cOli2?, cOli5?, cOli6?)*cOlf(cOli3?, cOli5?, cOli7?)*cOlf(cOli4?, cOli6?, cOli7?) = CA^2/2*delta(cOli1, cOli2);
    id cOlf(cOli1?,cOli4?,cOli5?)*cOlf(cOli2?,cOli4?,cOli6?)*cOlf(cOli3?,cOli5?,cOli6?) = CA/2*cOlf(cOli1, cOli2, cOli3);
    id T(i1?, cOli1?)*delta(cOli1?, cOli2?) = T(i1, cOli2);
    id cOlf(cOli1?, cOli2?, cOli3?)*delta(cOli3?, cOli4?) = cOlf(cOli1, cOli2, cOli4);
  endrepeat;
  .sort
* rename
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
sum cOli1,...,cOli8;

.sort
*if(count(summe, 1)<2) discard;
renumber 1;
.sort
b cOlf, T, delta, d, CA, summe, cOlTr, si1i2, si1i3,si2i3;
*b summe;
print+s difference, difference2;
.end

