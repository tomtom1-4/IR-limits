#-
#: IncDir procedures
#: SmallExtension 200M
#: WorkSpace      4G
#: MaxTermSize 4M
Off Statistics;

#include declarations.h
auto i mui;
auto v pi;
cf gamma11, gamma20, gamma21, gamma21Conj, j1;
cf phase(symmetric), phaseConj(symmetric);
cf Tsym;
f J10, J11, Gamma20, Gamma21, Gamma21Conj;
cf C, summeP(symmetric);
Format 160;

l J20 = 1/2*(J10(pi1, c1, q1, i1)*J10(pi2, c2, q2, i2) + J10(pi2, c2, q2, i2)*J10(pi1, c1, q1, i1))
  + Gamma20(pi1, c1, c2, q1, q2, i1);

l J21 = J11(pi1, pi2, c1, q1, i1, i2)*J10(pi3, c2, q2, i3) + J11(pi1, pi2, c2, q2, i1, i2)*J10(pi3, c1, q1, i3)
  + Gamma21(pi1, pi2, c1, c2, q1, q2, i1, i2);

l J20Conj = 1/2*(J10(pi2, c2, q2, i2)*J10(pi1, c1, q1, i1) + J10(pi1, c1, q1, i1)*J10(pi2, c2, q2, i2))
  - Gamma20(pi1, c1, c2, q1, q2, i1);

l J21Conj = -J10(pi3, c2, q2, i3)*J11(pi1, pi2, c1, q1, i1, i2) - J10(pi3, c1, q1, i3)*J11(pi1, pi2, c2, q2, i1, i2)
  + Gamma21Conj(pi1, pi2, c1, c2, q1, q2, i1, i2);

sum mu1, cOli1,cOli2;
.sort
#do i=1,3
  id J10(pi1?, c`i', q`i', i1?) = j1(pi1, q`i', mu`i')*E`i'(mu`i')*T(i1, c`i')*summe(i1);
  sum mu`i';
  id J11(pi1?, pi2?, c`i', q`i', i1?, i2?) = summeP(i1)*summeP(i2)*i_*cOlf(c`i', cOli1, cOli2)*T(i1, cOli1)*T(i2, cOli2)*gamma11(pi1, pi2, q`i', mu`i')*E`i'(mu`i')*phase(i1, i2);
  sum cOli1, cOli2, mu`i';
  #do k=`i'+1,3
    id once Gamma20(pi2?, c`i', c`k', q`i', q`k', i2?) = summe(i2)*i_*cOlf(c`i', c`k', cOli1)*T(i2, cOli1)*gamma20(pi2, q`i', q`k', mu`i', mu`k')*E`i'(mu`i')*E`k'(mu`k');
    sum cOli1, mu`i', mu`k';
    id once Gamma21(pi1?, pi2?, c`i', c`k', q`i', q`k', i1?, i2?) = summeP(i1)*summeP(i2)*cOlf(c`i', cOli1, cOli2)*cOlf(c`k', cOli3, cOli2)*T(i1, cOli1)*T(i2, cOli3)
        *gamma21(pi1, pi2, q`i', q`k', mu`i', mu`k')*E`i'(mu`i')*E`k'(mu`k');
    sum cOli1,...,cOli3, mu`i', mu`k';
    id once Gamma21Conj(pi1?, pi2?, c`i', c`k', q`i', q`k', i1?, i2?) = summeP(i1)*summeP(i2)*cOlf(c`i', cOli1, cOli2)*cOlf(c`k', cOli3, cOli2)*T(i1, cOli1)*T(i2, cOli3)
        *gamma21Conj(pi1, pi2, q`i', q`k', mu`i', mu`k')*E`i'(mu`i')*E`k'(mu`k');
    sum cOli1,...,cOli3, mu`i', mu`k';
  #enddo
#enddo

.sort
* rename color indices
InExpression J20Conj, J21Conj;
  mul replace_(i1, i4, i2, i5, i3, i6, pi1, pi4, pi2, pi5, pi3, pi6, phase, phaseConj);
  mul replace_(N1_?, cOli1, N2_?, cOli2, N3_?, cOli3, N4_?, cOli4, N5_?, cOli5, N6_?, cOli6, N7_?, cOli7, N8_?, cOli8, N9_?, cOli9, N10_?, cOli10);
EndInExpression;
.sort
l test1 = (T(i1, cOli1)*(T(i2, cOli2)*T(i3, cOli3) + T(i3, cOli3)*T(i2, cOli2))*T(i4, cOli4)
         + T(i4, cOli4)*(T(i3, cOli3)*T(i2, cOli2) + T(i2, cOli2)*T(i3, cOli3))*T(i1, cOli1))
        *cOlf(cOli1, cOli4, cOli5)*cOlf(cOli2, cOli3, cOli5)*(1 + replace_(i_, -i_, phase, phaseConj, gamma21, gamma21Conj))
        *(summe(i1)*summe(i2)*summe(i3)*summe(i4)*(-gamma11(pi3,pi4,q1,N3_?)*j1(pi1,q1,N3_?)*j1(pi2,q2,N4_?)*j1(pi1,q2,N4_?)*phase(i3,i4)/4
                                                   -gamma11(pi3,pi4,q2,N3_?)*j1(pi1,q2,N3_?)*j1(pi2,q1,N4_?)*j1(pi1,q1,N4_?)*phase(i3,i4)/4
                                                   -gamma11(pi2,pi3,q1,N6_?)*gamma20(pi1,q1,q2,N6_?,N7_?)*j1(pi4,q2,N7_?)*phase(i2,i3)/4
                                                   -gamma11(pi2,pi3,q2,N6_?)*gamma20(pi1,q2,q1,N6_?,N7_?)*j1(pi4,q1,N7_?)*phase(i2,i3)/4)
          +summeP(i1)*summeP(i2)*summe(i3)*summe(i4)*(+gamma21(pi1,pi2,q1,q2,N6_?,N7_?)*j1(pi3,q2,N7_?)*j1(pi4,q1,N6_?)/8
                                                      +gamma21(pi1,pi2,q2,q1,N6_?,N7_?)*j1(pi3,q1,N7_?)*j1(pi4,q2,N6_?)/8));

l test2 = i_*cOlf(cOli1, cOli2, cOli3)*T(i1, cOli1)*T(i2, cOli2)*T(i3, cOli3)*(1 + replace_(q1,q2,q2,q1))*(1 + replace_(i_, -i_, phase, phaseConj, gamma21, gamma21Conj))
  *(summeP(i1)*summeP(i2)*summeP(i3)*(-gamma11(pi2,pi3,q1,N4_?)*gamma20(pi1,q1,q2,N4_?,N5_?)*j1(pi1,q2,N5_?)*phase(i2,i3)*CA/2
                                   +gamma11(pi2,pi3,q1,N4_?)*gamma20(pi1,q1,q2,N4_?,N5_?)*j1(pi2,q2,N5_?)*phase(i2,i3)*CA/2
                                   +gamma11(pi2,pi3,q1,N4_?)*gamma20(pi2,q1,q2,N4_?,N5_?)*j1(pi1,q2,N5_?)*phase(i2,i3)*CA/2
                                   -gamma20(pi1,q1,q2,N4_?,N5_?)*gamma21(pi2,pi3,q1,q2,N4_?,N5_?)*CA/4
                                   +gamma11(pi2,pi3,q1,N4_?)*j1(pi1,q1,N4_?)*j1(pi1,q2,N5_?)*j1(pi2,q2,N5_?)*phase(i2,i3)*CA*3/4
                                   +gamma11(pi2,pi3,q1,N4_?)*j1(pi1,q1,N4_?)*j1(pi2,q2,N5_?)*j1(pi3,q2,N5_?)*phase(i2,i3)*CA*1/2
                                   -gamma11(pi2,pi3,q1,N4_?)*j1(pi1,q2,N5_?)*j1(pi2,q1,N4_?)*j1(pi2,q2,N5_?)*phase(i2,i3)*CA*1/4
                                   -gamma21(pi1,pi2,q1,q2,N4_?,N5_?)*j1(pi1,q1,N4_?)*j1(pi3,q2,N5_?)*CA/2
                                   +gamma21(pi1,pi2,q1,q2,N4_?,N5_?)*j1(pi2,q1,N4_?)*j1(pi3,q2,N5_?)*CA/4));

l test3= summeP(i1)*summeP(i2)*T(i1, cOli1)*T(i2, cOli1)*(1 + replace_(q1, q2, q2, q1))*(1 + replace_(i_, -i_, phase, phaseConj, gamma21, gamma21Conj))
  *(gamma21(pi1,pi2,q1,q2,N2_?,N3_?)*j1(pi1,q1,N2_?)*j1(pi1,q2,N3_?)*CA^2/8
   -gamma21(pi1,pi2,q1,q2,N2_?,N3_?)*j1(pi1,q1,N2_?)*j1(pi2,q2,N3_?)*CA^2/8
   -gamma11(pi1,pi2,q1,N2_?)*j1(pi1,q1,N2_?)*j1(pi1,q2,N3_?)*j1(pi2,q2,N3_?)*phase(i1,i2)*CA^2
   +gamma20(pi1,q1,q2,N2_?,N3_?)*gamma21(pi1,pi2,q1,q2,N2_?,N3_?)*CA^2/4);

l TreeG = -(T(i4,cOli1)*T(i5,cOli1)*j1(pi4,q1,N1_?)*j1(pi5,q1,N1_?))*summe(i4)*summe(i5);
l InterferenceG = -summeP(i1)*summeP(i2)*T(i1,cOli2)*T(i2,cOli2)*CA*(-gamma11(pi1,pi2,q2,N2_?)*j1(pi1,q2,N2_?))*(phase(i1,i2) + phaseConj(i1,i2))
  -summeP(i1)*summeP(i2)*summeP(i3)*T(i1,cOli2)*T(i2,cOli3)*T(i3,cOli4)*i_*cOlf(cOli2,cOli3,cOli4)
    *gamma11(pi2,pi3,q2,N2_?)*j1(pi1,q2,N2_?)*(phase(i2,i3) - phaseConj(i2,i3));
l test = (TreeG*InterferenceG + InterferenceG*TreeG)/2*(1 + replace_(q1,q2,q2,q1));
l interference = J21Conj*J20 + J20Conj*J21 - test - test1 - test2 - test3;
sum cOli1,...,cOli10;
.sort

* Use spin sums to replace Polarization vectors
ct EE1, EE2;
drop J20, J21, J20Conj, J21Conj;
*drop interference;
toTensor functions E1, EE1;
toTensor functions E2, EE2;
id EE1(mu1?, mu2?) = -d_(mu1, mu2);
id EE2(mu1?, mu2?) = -d_(mu1, mu2);
sum c1, c2, c3;
.sort
renumber 1;
.sort

#do l=1,2
* rename
#do i=1,5
  #do j=`i'+1,6
    if((match(summe(i`i'))==0) && (match(summeP(i`i'))==0));
      if(match(summe(i`j'))||match(summeP(i`j')));
        mul replace_(i`j', i`i', pi`j', pi`i');
        redefine j "10";
        redefine i "1";
      endif;
    endif;
    .sort
  #enddo
#enddo


if(match(summeP(i3))) mul replace_(i3, i1, i1, i3, pi3, pi1, pi1, pi3);
if(match(summeP(i4))) mul replace_(i4, i2, i2, i4, pi4, pi2, pi2, pi4);

.sort
if(count(summeP, 1)==0) id once summe(i1) = summeP(i1);
.sort
#do i=2,5
  if(count(summeP, 1)==1)
    id once summe(i`i')*summeP(i1) = summeP(i1)*summeP(i`i') + summeP(i1)*replace_(i`i', i1, pi`i', pi1);
#enddo
.sort
#do i=3,5
if(count(summeP, 1)==2)
  id summe(i`i')*summeP(i1)*summeP(i2) = summeP(i1)*summeP(i2)*summeP(i`i') + summeP(i1)*summeP(i2)*(replace_(i`i', i1, pi`i', pi1) + replace_(i`i', i2, pi`i', pi2));
#enddo
.sort
#do i=4,5
if(count(summeP, 1)==3)
  id summe(i`i')*summeP(i1)*summeP(i2)*summeP(i3) = summeP(i1)*summeP(i2)*summeP(i3)*summeP(i`i')
                                                  + summeP(i1)*summeP(i2)*summeP(i3)*(replace_(i`i', i1, pi`i', pi1)
                                                                                    + replace_(i`i', i2, pi`i', pi2)
                                                                                    + replace_(i`i', i3, pi`i', pi3));
#enddo
.sort
if(count(summeP, 1)==3)
  id summe(i5)*summeP(i1)*summeP(i2)*summeP(i4) = summeP(i1)*summeP(i2)*summeP(i4)*summeP(i5)
                                                  + summeP(i1)*summeP(i2)*summeP(i4)*(replace_(i5, i1, pi5, pi1)
                                                                                    + replace_(i5, i2, pi5, pi2)
                                                                                    + replace_(i5, i4, pi5, pi4));
if(count(summeP, 1)==4)
  id summe(i5)*summeP(i1)*summeP(i2)*summeP(i3)*summeP(i4) = summeP(i1)*summeP(i2)*summeP(i3)*summeP(i4)*summeP(i5)
                                                  + summeP(i1)*summeP(i2)*summeP(i3)*summeP(i4)*(replace_(i5, i1, pi5, pi1)
                                                                                               + replace_(i5, i2, pi5, pi2)
                                                                                               + replace_(i5, i3, pi5, pi3)
                                                                                               + replace_(i5, i4, pi5, pi4));
#enddo
b cOlf, T, Tc, summe, summeP;
print+s ;
.sort
*cf Tc;

#do i=1,6
  id once T(i1?, cOli1?) = Tc(i1, cOli1);
  repeat;
    id once Tc(i1?, ?args)*T(i1?, cOli1?) = Tc(i1, ?args, cOli1);
  endrepeat;
#enddo

.sort

*mul (1 - replace_(q1, q2, q2, q1));

#do l=1,5
#do i=1,5
  #do j=`i'+1,6
    id once j1(pi`j', q1, mu1?)*j1(pi`i', q2, mu2?) = j1(pi`j', q1, mu1)*j1(pi`i', q2, mu2)*replace_(i`i', i`j', i`j', i`i', pi`i', pi`j', pi`j', pi`i');
  #enddo
#enddo
#do i=1,5
  #do j=`i'+1,6
    id once gamma11(pi`i', pi10?, q1?, mu1?)*j1(pi`j', q2?, mu2?) = gamma11(pi`i', pi10, q1, mu1)*j1(pi`j', q2, mu2)*replace_(i`i', i`j', i`j', i`i', pi`i', pi`j', pi`j', pi`i');
    id once gamma11(pi10?, pi`i', q1?, mu1?)*j1(pi`j', q2?, mu2?) = gamma11(pi10, pi`i', q1, mu1)*j1(pi`j', q2, mu2)*replace_(i`i', i`j', i`j', i`i', pi`i', pi`j', pi`j', pi`i');
    id once gamma21(pi1,pi2,q1,q2,mu1?,mu2?)*j1(pi3,q1,mu1?)*j1(pi3,q2,mu2?) = gamma21(pi1,pi2,q1,q2,mu1,mu2)*j1(pi3,q1,mu1)*j1(pi3,q2,mu2)*replace_(i3, i1, i1, i2, i2, i3, pi3, pi1, pi1, pi2, pi2, pi3);
    id once gamma21Conj(pi1,pi2,q1,q2,mu1?,mu2?)*j1(pi3,q1,mu1?)*j1(pi3,q2,mu2?) = gamma21Conj(pi1,pi2,q1,q2,mu1,mu2)*j1(pi3,q1,mu1)*j1(pi3,q2,mu2)*replace_(i3, i1, i1, i2, i2, i3, pi3, pi1, pi1, pi2, pi2, pi3);
    id once gamma20(pi`i', q1, q2, mu1?, mu2?)*j1(pi`j', q3?, mu3?) = gamma20(pi`i', q1, q2, mu1, mu2)*j1(pi`j', q3, mu3)*replace_(i`i', i`j', i`j', i`i', pi`i', pi`j', pi`j', pi`i');
  #enddo
#enddo
#do i=1,5
  #do j=`i'+1,6
    id gamma11(pi`j', pi`i', q1?, mu1?) = -gamma11(pi`i', pi`j', q1, mu1);
    id gamma21(pi`j', pi`i', q1, q2, mu1?, mu2?) = gamma21(pi`j', pi`i', q1, q2, mu1, mu2)*replace_(i`i', i`j', i`j', i`i', pi`i', pi`j', pi`j', pi`i');
    id gamma21Conj(pi`j', pi`i', q1, q2, mu1?, mu2?) = gamma21Conj(pi`j', pi`i', q1, q2, mu1, mu2)*replace_(i`i', i`j', i`j', i`i', pi`i', pi`j', pi`j', pi`i');
  #enddo
#enddo
#do i=1,5
  #do j=`i'+1,6
    id once gamma11(pi`i', pi10?, q3?, mu1?)*gamma20(pi`j', q1, q2, mu2?, mu3?) = gamma11(pi`i', pi10, q3, mu1)*gamma20(pi`j', q1, q2, mu2, mu3)*replace_(i`i', i`j', i`j', i`i', pi`i', pi`j', pi`j', pi`i');
  #enddo
#enddo
#do i=1,5
  #do j=`i'+1,6
    id gamma11(pi`j', pi`i', q1?, mu1?) = -gamma11(pi`i', pi`j', q1, mu1);
  #enddo
#enddo

id gamma11(pi1?, pi1?, q1?, mu1?) = 0;
id j1(pi?, q1?, mu1?)*j1(pi?, q2?, mu1?) = 0;
id gamma20(pi1?, q2, q1, mu2?, mu1?) = -gamma20(pi1, q1, q2, mu1, mu2);
id gamma21(pi1?, pi2?, q2, q1, mu2?, mu1?) = gamma21(pi2, pi1, q1, q2, mu1, mu2);
id gamma21Conj(pi1?, pi2?, q2, q1, mu2?, mu1?) = gamma21Conj(pi2, pi1, q1, q2, mu1, mu2);
#enddo

.sort
#do l=1,4
  renumber 1;
  .sort
* jacobi
*  id cOlf(N1_?,N4_?,N5_?)*cOlf(N2_?,N3_?,N5_?) = - cOlf(N1_?,N2_?,N5_?)*cOlf(N3_?,N4_?,N5_?) - cOlf(N1_?,N3_?,N5_?)*cOlf(N4_?,N2_?,N5_?);
*  id cOlf(N1_?,N3_?,N5_?)*cOlf(N4_?,N2_?,N5_?) = - cOlf(N1_?,N2_?,N5_?)*cOlf(N3_?,N4_?,N5_?) - cOlf(N1_?,N4_?,N5_?)*cOlf(N2_?,N3_?,N5_?);
  id cOlf(N1_?,N2_?,N5_?)*cOlf(N3_?,N4_?,N5_?) = - cOlf(N1_?,N4_?,N5_?)*cOlf(N2_?,N3_?,N5_?) - cOlf(N1_?,N3_?,N5_?)*cOlf(N4_?,N2_?,N5_?);
  repeat;
    id once cOlf(cOli1?, cOli2?, cOli3?)*cOlf(cOli1?, cOli2?, cOli4?) = CA*delta(cOli3, cOli4);
    id once cOlf(cOli1?, cOli3?, cOli4?)*cOlf(cOli2?, cOli5?, cOli6?)*cOlf(cOli3?, cOli5?, cOli7?)*cOlf(cOli4?, cOli6?, cOli7?) = CA^2/2*delta(cOli1, cOli2);
    id once cOlf(cOli1?,cOli4?,cOli5?)*cOlf(cOli2?,cOli4?,cOli6?)*cOlf(cOli3?,cOli5?,cOli6?) = CA/2*cOlf(cOli1, cOli2, cOli3);
    id once Tc(i1?, ?arg1, cOli1?, ?arg2)*delta(cOli1?, cOli2?) = Tc(i1, ?arg1, cOli2, ?arg2);
    id once cOlf(cOli1?, cOli2?, cOli3?)*delta(cOli3?, cOli4?) = cOlf(cOli1, cOli2, cOli4);
    id once Tsym(?args1, cOli1?, ?args2)*delta(cOli1?, cOli3?) = Tsym(?args1, cOli3, ?args2);
    id once delta(cOli1?, cOli2?)*delta(cOli2?, cOli3?) = delta(cOli1, cOli3);
    id once Tc(i1?, ?arg1, cOli1?, cOli2?, ?arg2)*cOlf(cOli1?, cOli2?, cOli3?) = 1/2*i_*cOlf(cOli1, cOli2, cOli4)*Tc(i1, ?arg1, cOli4, ?arg2)*cOlf(cOli1, cOli2, cOli3);
    sum cOli4;
    id once Tc(i1?, ?arg1, cOli1?, cOli1?, ?arg2) = C(i1)*Tc(i1, ?arg1, ?arg2);
    id once Tc(i1, cOli1?, cOli2?)*Tc(i2, cOli2?) = Tc(i1, cOli2, cOli1)*Tc(i2, cOli2) + i_*cOlf(cOli1, cOli2, cOli11)*Tc(i1, cOli11)*Tc(i2, cOli2);
    id once Tc(i1?,cOli1?,cOli2?,cOli3?)*cOlf(cOli1?,cOli3?,cOli4?) = 1/2*Tc(i1,cOli2,cOli4)*i_*CA + Tc(i1,cOli11,cOli3)*cOlf(cOli1,cOli2,cOli11)*cOlf(cOli1,cOli3,cOli4)*i_;
    sum cOli11;
    id once Tc(i1?, cOli1?, cOli2?) = Tsym(i1, cOli1, cOli2)/2 + i_*cOlf(cOli1, cOli2, cOli12)*Tc(i1, cOli12)/2;
    sum cOli12;
    id once Tsym(i1?, cOli1?, cOli2?)*cOlf(cOli1?, cOli2?, cOli3?) = 0;
  endrepeat;
  .sort

  mul replace_(N1_?, cOli1, N2_?, cOli2, N3_?, cOli3, N4_?, cOli4, N5_?, cOli5, N6_?, cOli6, N7_?, cOli7, N8_?, cOli8, N9_?, cOli9, N10_?, cOli10);
  #do i=1,5
    #do j=`i'+1,6
      repeat;
        id once Tc(i1?, ?args1, cOli`j', cOli`i', ?args2) = Tc(i1, ?args1, cOli`i', cOli`j', ?args2) + i_*cOlf(cOli`j', cOli`i', cOli11)*Tc(i1, ?args1, cOli11, ?args2);
        sum cOli11;
      endrepeat;
      id Tsym(i1?, cOli`j', cOli`i') = Tsym(i1, cOli`i', cOli`j');
    #enddo
  #enddo
  sum cOli1,...,cOli10;
*  #if `l'==2
*    #do i=1,5
*      #do j=`i'+1,6
*        #do l=1,5
*          #do m=`l'+1,6
*            id Tc(i`i', cOli`m')*Tc(i`j', cOli`l') = Tc(i`i', cOli`m')*Tc(i`j', cOli`l')*replace_(i`i', i`j', i`j', i`i', pi`i', pi`j', pi`j', pi`i');
*          #enddo
*        #enddo
*      #enddo
*    #enddo
*  #endif
*  sum cOli1,...,cOli10;
#enddo
id phase(i1?, i2?)*phaseConj(i1?, i2?) = 1;
.sort
*mul (1 - replace_(i_, -i_, phaseConj, phase, phase, phaseConj, gamma21Conj, gamma21, gamma21, gamma21Conj))/2;
*mul replace_(i_, -i_, phaseConj, phase)/2;
*if(occurs(phaseConj) || occurs(gamma21Conj)); discard;
*else if((occurs(phase)==0) && (occurs(phaseConj)==0)); mul 1/2;
*endif;

if(count(Tsym, 1)==1);
  #do i=2,3
    id Tsym(i`i', cOli1?, cOli2?) = Tsym(i`i', cOli1, cOli2)*replace_(i`i', i1, i1, i`i', pi`i', pi1, pi1, pi`i');
  #enddo
endif;
*id Tc(i1,N1_?)*Tc(i2,N2_?)*Tc(i3,N3_?)*Tc(i4,N4_?)*cOlf(N1_?,N3_?,N5_?)*cOlf(N2_?,N4_?,N5_?)
*  = Tc(i1,N1_?)*Tc(i2,N2_?)*Tc(i3,N3_?)*Tc(i4,N4_?)*cOlf(N1_?,N3_?,N5_?)*cOlf(N2_?,N4_?,N5_?)*replace_(i2, i3, i3, i2, pi2, pi3, pi3, pi2);
*id Tc(i1,N1_?)*Tc(i2,N2_?)*Tc(i3,N3_?)*Tc(i4,N4_?)*cOlf(N1_?,N4_?,N5_?)*cOlf(N2_?,N3_?,N5_?)
*  = Tc(i1,N1_?)*Tc(i2,N2_?)*Tc(i3,N3_?)*Tc(i4,N4_?)*cOlf(N1_?,N4_?,N5_?)*cOlf(N2_?,N3_?,N5_?)*replace_(i4, i2, i2, i4, pi4, pi2, pi2, pi4);

id Tc(i1,N1_?)*Tc(i2,N2_?)*Tc(i3,N3_?)*Tc(i4,N4_?)*cOlf(N1_?,N4_?,N5_?)*cOlf(N2_?,N3_?,N5_?)
  = Tc(i1,N1_?)*Tc(i2,N2_?)*Tc(i3,N3_?)*Tc(i4,N4_?)*cOlf(N1_?,N4_?,N5_?)*cOlf(N2_?,N3_?,N5_?)*replace_(i4, i3, i3, i4, pi4, pi3, pi3, pi4);
id gamma11(pi1,pi2,q1?,mu1?)*j1(pi1,q2?,mu2?)*j1(pi2,q1?,mu1?)
  = gamma11(pi1,pi2,q1,mu1)*j1(pi1,q2,mu2)*j1(pi2,q1,mu1)*replace_(i1,i2,i2,i1,pi1,pi2,pi2,pi1);
renumber 1;
id Tsym(i1,N5_?,N4_?) = Tsym(i1, N4_?, N5_?);
id Tsym(i1,N2_?,N1_?) = Tsym(i1, N1_?, N2_?);
id gamma11(pi4,pi3,q1?,cOli1?) = - gamma11(pi3,pi4,q1,cOli1);
id gamma11(pi2,pi1,q1?,cOli1?) = - gamma11(pi1,pi2,q1,cOli1);
id gamma11(pi1,pi2,q1?,mu1?)*j1(pi1,q2?,mu2?)*j1(pi2,q1?,mu1?)
  = -gamma11(pi2,pi1,q1,mu1)*j1(pi1,q2,mu2)*j1(pi2,q1,mu1)*replace_(i1,i2,i2,i1,pi1,pi2,pi2,pi1);
id Tsym(i2,N1_?,N2_?) = Tsym(i2,N1_?,N2_?)*replace_(i2,i1,i1,i2,pi2,pi1,pi1,pi2);
id gamma11(pi2,pi1,q1?,cOli1?) = - gamma11(pi1,pi2,q1,cOli1);
id gamma11(pi1,pi2,q1?,N3_?)*j1(pi1,q2?,N4_?)*j1(pi2,q1?,N3_?)*j1(pi3,q2?,N4_?)
  = gamma11(pi1,pi2,q1,N3_?)*j1(pi1,q2,N4_?)*j1(pi2,q1,N3_?)*j1(pi3,q2,N4_?)*replace_(i2,i3,i3,i2,pi2,pi3,pi3,pi2);
renumber 1;
.sort
id Tsym(i1,N2_?,N1_?) = Tsym(i1, N1_?, N2_?);
id gamma11(pi1,pi2,q1,N6_?)*gamma20(pi3,q1,q2,N6_?,N7_?)*j1(pi1,q2,N7_?)
  = gamma11(pi1,pi2,q1,N6_?)*gamma20(pi3,q1,q2,N6_?,N7_?)*j1(pi1,q2,N7_?)*replace_(i2,i3,i3,i2,pi2,pi3,pi3,pi2);
.sort
id cOlf(N1_?,N3_?,N5_?)*cOlf(N2_?,N3_?,N4_?) = - cOlf(N1_?,N3_?,N2_?)*cOlf(N4_?,N3_?,N5_?) - cOlf(N1_?,N3_?,N4_?)*cOlf(N5_?,N3_?,N2_?);
.sort
id Tsym(i1,N5_?,N4_?) = Tsym(i1, N4_?, N5_?);
id gamma21(pi1,pi2,q1,q2,N2_?,N3_?)*j1(pi2,q1,N2_?)*j1(pi2,q2,N3_?)
  = gamma21(pi1,pi2,q1,q2,N2_?,N3_?)*j1(pi2,q1,N2_?)*j1(pi2,q2,N3_?)*replace_(i1,i2,i2,i1,pi1,pi2,pi2,pi1);
id gamma21(pi2,pi4,q1,q2,N6_?,N7_?)*j1(pi1,q1,N6_?)*j1(pi3,q2,N7_?)
  = gamma21(pi2,pi4,q1,q2,N6_?,N7_?)*j1(pi1,q1,N6_?)*j1(pi3,q2,N7_?)*replace_(i2,i3,i3,i2,pi2,pi3,pi3,pi2);
id gamma21Conj(pi1,pi2,q1,q2,N2_?,N3_?)*j1(pi2,q1,N2_?)*j1(pi2,q2,N3_?)
  = gamma21Conj(pi1,pi2,q1,q2,N2_?,N3_?)*j1(pi2,q1,N2_?)*j1(pi2,q2,N3_?)*replace_(i1,i2,i2,i1,pi1,pi2,pi2,pi1);
id gamma21Conj(pi2,pi4,q1,q2,N6_?,N7_?)*j1(pi1,q1,N6_?)*j1(pi3,q2,N7_?)
  = gamma21Conj(pi2,pi4,q1,q2,N6_?,N7_?)*j1(pi1,q1,N6_?)*j1(pi3,q2,N7_?)*replace_(i2,i3,i3,i2,pi2,pi3,pi3,pi2);
.sort
id cOlf(N1_?,N3_?,N5_?)*cOlf(N2_?,N3_?,N4_?) = - cOlf(N1_?,N3_?,N2_?)*cOlf(N4_?,N3_?,N5_?) - cOlf(N1_?,N3_?,N4_?)*cOlf(N5_?,N3_?,N2_?);
.sort
id gamma11(pi1,pi3,q2,N6_?)*gamma20(pi2,q1,q2,N7_?,N6_?)*j1(pi1,q1,N7_?) = gamma11(pi1,pi3,q2,N6_?)*gamma20(pi2,q1,q2,N7_?,N6_?)*j1(pi1,q1,N7_?)*replace_(i3,i2,i2,i3,pi3,pi2,pi2,pi3);
id gamma11(pi2,pi3,q2,N4_?)*gamma20(pi1,q1,q2,N5_?,N4_?)*j1(pi3,q1,N5_?) = gamma11(pi2,pi3,q2,N4_?)*gamma20(pi1,q1,q2,N5_?,N4_?)*j1(pi3,q1,N5_?)*replace_(i3,i2,i2,i3,pi3,pi2,pi2,pi3);
.sort
id cOlf(N1_?,N3_?,N4_?)*cOlf(N2_?,N3_?,N5_?) = - cOlf(N1_?,N3_?,N2_?)*cOlf(N5_?,N3_?,N4_?) - cOlf(N1_?,N3_?,N5_?)*cOlf(N4_?,N3_?,N2_?);
id cOlf(N3_?,N4_?,N5_?)*Tsym(i1,N5_?,N4_?) = 0;
.sort
id Tsym(i1?, N5_?, N4_?) = Tsym(i1, N4_?, N5_?);
renumber 1;
.sort
id cOlf(N3_?,N4_?,N5_?)*Tsym(i1,N4_?,N5_?) = 0;
.sort
id gamma11(pi2,pi3,q1,N4_?)*gamma20(pi1,q1,q2,N4_?,N5_?)*j1(pi3,q2,N5_?) = gamma11(pi2,pi3,q1,N4_?)*gamma20(pi1,q1,q2,N4_?,N5_?)*j1(pi3,q2,N5_?)*replace_(i3, i2, i2, i3, pi2, pi3, pi3, pi2);
id gamma11(pi3,pi2,q1?,mu1?) = - gamma11(pi2,pi3,q1,mu1);
id gamma20(pi3,q1,q2,N4_?,N5_?)*gamma21(pi1,pi2,q1,q2,N4_?,N5_?) = gamma20(pi3,q1,q2,N4_?,N5_?)*gamma21(pi1,pi2,q1,q2,N4_?,N5_?)*replace_(i3, i1, i1, i2, i2, i3, pi3, pi1, pi1, pi2, pi2, pi3);
*if(match(summeP(i5))==0)discard;
*if(match(Tc(i1?, cOli1?, cOli2?))==0)discard;
*mul (1 + replace_(q1, q2, q2, q1))/2;
*if(count(gamma21, 1)==1)discard;
*if(count(j1, 1)!=3)discard;
*if(count(Tc, 1)!=4)discard;

b T, cOlf, summe, summeP, Tc, Tsym;
*b gamma11, j1, gamma21, gamma20;
print+s interference, test, test1;

.end