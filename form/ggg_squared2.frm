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
s marker1, marker2;
s m1sq, m2sq, m3sq, m4sq, m5sq, m6sq;

l J2 = 1/2*(J1(pi1, c1, q1, i1)*J1(pi2, c2, q2, i2) + J1(pi2, c2, q2, i2)*J1(pi1, c1, q1, i1))
  + Gamma2(pi1, c1, c2, q1, q2, i1);

#do i=1,3
#do j=1,3
  id J1(pi`i', c`j', q`j', i`i') = j1(pi`i', q`j', mu`j')*E`j'(mu`j')*T(i`i', c`j')*summe(i`i');
  sum cOli1, cOli2, mu`j';
  #do k=`j'+1,3
    id Gamma2(pi2?, c`j', c`k', q`j', q`k', i2?) = summe(i2)*i_*cOlf(c`j', c`k', cOli1)*T(i2, cOli1)*gamma2(pi2, q`j', q`k', mu`j', mu`k')*E`j'(mu`j')*E`k'(mu`k');
    sum cOli1, mu`j', mu`k';
  #enddo
#enddo
#enddo

sum cOli1, cOli2;
.sort

* Define complex conjugate of J3
* All expressions are symmetrized, i.e. the color operators don't need to be explicitly reversed (?)
l J2Conj = J2*replace_(i_, -i_, i1, i3, i2, i4, pi1, pi3, pi2, pi4);
.sort
* rename color indices
InExpression J2Conj;
  id T(i1?, cOli11?)*cOlf(c1?, c2?, cOli11?) = T(i1, cOli1)*cOlf(c1, c2, cOli1);
  id cOlf(c1?,c2?,cOli12?)*cOlf(c3?,cOli1?,cOli12?) = cOlf(c1,c2,cOli2)*cOlf(c3,cOli1,cOli2);
EndInExpression;
.sort
b summe, cOlf, T;
print+s J2, J2Conj;
.sort

l square2 = J2Conj*J2;
sum cOli1, cOli2;
.sort
b summe cOlf, T;
print+s square2;
.sort

l Catani2 = 1/2*(W1(pi1, pi2, q1, i1, i2)*W1(pi3, pi4, q2, i3, i4) + W1(pi1, pi2, q2, i1, i2)*W1(pi3, pi4, q1, i3, i4))
  + W2(pi1, pi2, q1, q2, i1, i2);
.sort

repeat;
  id once W1(pi1?, pi2?, q1?, i1?, i2?) = -summe(i1)*summe(i2)*T(i1, cOli1)*T(i2, cOli1)*j1(pi1, q1, mu1)*j1(pi2, q1, mu1);
  sum cOli1, mu1;
  id once W2(pi1?, pi2?, q1?, q2?, i1?, i2?) = -CA*summe(i1)*summe(i2)*T(i1, cOli1)*T(i2, cOli1)*1/2
    *(S2(pi1, pi2, q1, q2)*2 + 0*S2(pi2, pi1, q1, q2) - 0*S2(pi1, pi1, q1, q2) - 0*S2(pi2, pi2, q1, q2));
  sum cOli1;
endrepeat;

id S2(pi1?, pi2?, q1?, q2?) = -(+ gamma2(pi1,q1,q2,N2_?,N3_?)*gamma2(pi2,q1,q2,N2_?,N3_?)
          + gamma2(pi1,q1,q2,N2_?,N3_?)*j1(pi1,q1,N2_?)*j1(pi2,q2,N3_?)
          - gamma2(pi2,q1,q2,N2_?,N3_?)*j1(pi1,q1,N2_?)*j1(pi2,q2,N3_?)
          + 3/4*j1(pi1,q1,N2_?)*j1(pi1,q2,N3_?)*j1(pi2,q1,N2_?)*j1(pi2,q2,N3_?)
          - 1/2*j1(pi1,q1,N2_?)^2*j1(pi1,q2,N3_?)*j1(pi2,q2,N3_?)
          - 1/2*j1(pi1,q1,N2_?)*j1(pi2,q1,N2_?)*j1(pi2,q2,N3_?)^2);

*id gamma2(pi1?, q1?, q2?, mu1?, mu2?) = den(pi1.q1 + pi1.q2)*(pi1(mu1)*pi1(mu2)*(den(2*pi1.q1) - den(2*pi1.q2)) + den(q1.q2)*(pi1(mu1)*q1(mu2) - pi1(mu2)*q2(mu1) + 1/2*d_(mu1, mu2)*(pi1.q2 - pi1.q1)));

*id j1(pi1?, q1?, mu1?) = pi1(mu1)*den(pi1.q1);

.sort
l difference2 = square2 - Catani2;
.sort


* Use spin sums to replace Polarization vectors
ct EE1, EE2, EE3;
drop J2, J2Conj;
toTensor functions E1, EE1;
toTensor functions E2, EE2;
b summe cOlf, T;
print+s square2;
.sort
id EE1(mu1?, mu2?) = -d_(mu1, mu2);
id EE2(mu1?, mu2?) = -d_(mu1, mu2);
sum c1, c2, c3;
.sort
renumber 1;
.sort
sum cOli1,...,cOli8;
mul replace_(N1_?, cOli1, N2_?, cOli2, N3_?, cOli3, N4_?, cOli4, N5_?, cOli5, N6_?, cOli6, N7_?, cOli7, N8_?, cOli8, N9_?, cOli9, N10_?, cOli10);
* order color operators
#do l=1,4
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
    id T(i1?, cOli1?)*delta(cOli1?, cOli2?) = T(i1, cOli2);
    id cOlf(cOli1?, cOli2?, cOli3?)*delta(cOli3?, cOli4?) = cOlf(cOli1, cOli2, cOli4);
  endrepeat;
  .sort
  sum cOli1,...,cOli8;
  mul replace_(N1_?, cOli1, N2_?, cOli2, N3_?, cOli3, N4_?, cOli4, N5_?, cOli5, N6_?, cOli6, N7_?, cOli7, N8_?, cOli8, N9_?, cOli9, N10_?, cOli10);
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
*mul (1 + replace_(q1, q2, q2, q1))/2;
*id once pi1.pi3 = pi1.pi3*replace_(i3,i2,i2,i3,pi3,pi2,pi2,pi3);
*id once pi2.pi3 = pi2.pi3*replace_(i2,i1,i1,i3,i3,i2,pi2,pi1,pi1,pi3,pi3,pi2);
*id once pi1.pi2*pi1.pi3 = pi1.pi2*pi1.pi3*(1 + replace_(i3, i2, i2, i3, pi3, pi2, pi2, pi3))/2;
*id once pi1.pi2 = pi1.pi2*(1 + replace_(i1, i2, i2, i1, pi1, pi2, pi2, pi1))/2;
.sort

#do i=1,2
repeat;
  #do k=1,2
    #do j=`k'+1,6
      id once T(i`j', cOli1?)*T(i`k', cOli2?)*summe(i`j')*summe(i`k') = T(i`k', cOli2)*T(i`j', cOli1)*summe(i`k')*summe(i`j') + i_*cOlf(cOli1, cOli2, cOli3)*T(i`k', cOli3)*summe(i`k')*replace_(i`j', i`k', pi`j', pi`k');
      sum cOli3;
    #enddo
  #enddo
  id cOlf(cOli1?, cOli2?, cOli3?)*cOlf(cOli1?, cOli2?, cOli4?) = CA*delta(cOli3, cOli4);
  id T(i1?, cOli1?)*delta(cOli1?, cOli2?) = T(i1, cOli2);
  id cOlf(cOli1?, cOli2?, cOli3?)*delta(cOli3?, cOli4?) = cOlf(cOli1, cOli2, cOli4);
endrepeat;
.sort
  id gamma2(pi1,q1,q2,mu1?,mu2?)*j1(pi2,q1,mu1?)*j1(pi3,q2,mu2?)*summe(i1)*summe(i2)*summe(i3) = gamma2(pi1,q1,q2,mu1,mu2)*j1(pi2,q1,mu1)*j1(pi3,q2,mu2)
                                                                                                  *summe(i1)*summe(i2)*summe(i3)*replace_(i2, i1, i1, i3, i3, i2, pi2, pi1, pi1, pi3, pi3, pi2);
  id gamma2(pi1,q1,q2,mu1?,mu2?)*j1(pi2,q2,mu2?)*j1(pi3,q1,mu1?) = gamma2(pi1,q1,q2,mu1,mu2)*j1(pi2,q2,mu2)*j1(pi3,q1,mu1)*replace_(i3, i1, i1, i3, pi3, pi1, pi1, pi3);
  id j1(pi1,q2,mu1?)*j1(pi2,q1,mu2?) = j1(pi1,q2,mu1)*j1(pi2,q1,mu2)*replace_(i1, i2, i2, i1, pi1, pi2, pi2, pi1);
.sort
#enddo

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

.sort
#do i=1,6
  id pi`i'.q1 = si`i'1/2;
  id pi`i'.q2 = si`i'2/2;
  id pi`i'.q3 = si`i'3/2;
  id pi`i'.pi`i' = m`i'sq;
  argument den;
    id pi`i'.q1 = si`i'1/2;
    id pi`i'.q2 = si`i'2/2;
    id pi`i'.q3 = si`i'3/2;
    id pi`i'.pi`i' = m`i'sq;
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

#call FullSimplify
#call Simplify


b cOlf, T, delta, d, CA, summe, cOlTr, si1i2, si1i3,si2i3;
print+s ;
.end

