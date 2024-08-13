#-
#: IncDir procedures
#: SmallExtension 200M
#: WorkSpace      4G
#: MaxTermSize 4M
Off Statistics;

#include declarations.h
auto i mui;
auto v pi;
cf gamma1, j1;
cf phase(symmetric), phaseConj(symmetric);
f Tsym1, Tsym2;
cf C, summeP;
Format 160;

l J0 = -summe(i1)*T(i1, c1)*j1(pi1, q1, mu1)*E1(mu1);

l J1 = summe(i1)*summe(i2)*i_*cOlf(c1, cOli1, cOli2)*T(i1, cOli1)*T(i2, cOli2)*gamma1(pi1, pi2, q1, mu1)*E1(mu1)*phase(i1, i2);

sum mu1, cOli1,cOli2;
.sort
* The color factors are symmetrized so an explicit reversal is not necessary
l J0Conj = J0*replace_(i_, -i_, i1, i3, i2, i4, pi1, pi3, pi2, pi4, phase, phaseConj);
l J1Conj = J1*replace_(i_, -i_, i1, i3, i2, i4, pi1, pi3, pi2, pi4, phase, phaseConj);
.sort
* rename color indices
InExpression J0Conj, J1Conj;
  mul replace_(N1_?, cOli1, N2_?, cOli2, N3_?, cOli3, N4_?, cOli4, N5_?, cOli5, N6_?, cOli6, N7_?, cOli7, N8_?, cOli8, N9_?, cOli9, N10_?, cOli10);
EndInExpression;
.sort

l interference = J1Conj*J0 + J0Conj*J1;
l test = T(i1, cOli1)*T(i2, cOli2)*T(i3, cOli3)*cOlf(cOli1, cOli2, cOli3)*summe(i1)*summe(i2)*summe(i3);
l square = J1Conj*J1;
sum cOli1,...,cOli10;
.sort

*(1 + replace_(i1, i3, i3, i1, i2, i4, i4, i2, pi1, pi3, pi3, pi1, pi2, pi4, pi4, pi2))/2;

.sort
* Use spin sums to replace Polarization vectors
ct EE1;
drop J0, J1, J0Conj, J1Conj;
toTensor functions E1, EE1;
id EE1(mu1?, mu2?) = -d_(mu1, mu2);
sum c1, c2, c3;
.sort
renumber 1;
.sort
id gamma1(pi1, pi2, q1, mu1?)*j1(pi3, q1, mu2?) = gamma1(pi1, pi2, q1, mu1)*j1(pi3, q1, mu2)*replace_(i3, i1, i1, i2, i2, i3, pi3, pi1, pi1, pi2, pi2, pi3);
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
b T, cOlf, summe, summeP;
print+s interference, square;
.sort
if(count(summe, 1)==3);
  id summe(i1)*summe(i2)*summe(i3) = summeP(i1)*summeP(i2)*summeP(i3) + summeP(i1)*summeP(i2)*(replace_(i3, i2, pi3, pi2) + replace_(i3, i1, pi3, pi1)) + summeP(i1)*summeP(i3)*replace_(i2, i1, pi2, pi1) + summeP(i1)*replace_(i3, i1, i2, i1, pi3, pi1, pi2, pi1);
else if(count(summe, 1)==4);
  id summe(i1)*summe(i2)*summe(i3)*summe(i4) = summeP(i1)*summeP(i2)*summeP(i3)*summeP(i4)
        + summeP(i1)*summeP(i2)*summeP(i3)*(replace_(i4, i3, pi4, pi3) + replace_(i4, i2, pi4, pi2) + replace_(i4, i1, pi4, pi1))
        + summeP(i1)*summeP(i2)*summeP(i4)*(replace_(i3, i2, pi3, pi2) + replace_(i3, i1, pi3, pi1))
        + summeP(i1)*summeP(i3)*summe(i4)*replace_(i2, i1, pi2, pi1)
        + summeP(i1)*summeP(i2)*(replace_(i4, i1, pi4, pi1, i3, i2, pi3, pi2) + replace_(i4, i2, pi4, pi2, i3, i1, pi3, pi1) + replace_(i4, i1, pi4, pi1, i3, i1, pi3, pi1) + replace_(i4, i2, pi4, pi2, i3, i2, pi3, pi2))
        + summeP(i1)*summeP(i3)*(replace_(i2, i1, pi2, pi1, i4, i3, pi4, pi3) + replace_(i4, i1, pi4, pi1, i2, i1, pi2, pi1))
        + summeP(i1)*summeP(i4)*(replace_(i3, i1, pi3, pi1, i2, i1, pi2, pi1))
        + summeP(i1)*replace_(i4, i1, pi4, pi1, i3, i1, pi3, pi1, i2, i1, pi2, pi1);

endif;

#do l=1,4
  #if `l'==1
*    id T(i1, cOli1?)*T(i2, cOli2?)*T(i3, cOli3?)*T(i4, cOli4?)*cOlf(cOli1?, cOli2?, cOli5?)*cOlf(cOli3?, cOli4?, cOli5?)
*        = Tsym2(i1, i2, i3, i4) - cOlf(cOli1, cOli2, cOli5)*cOlf(cOli3, cOli4, cOli5)*(T(i3, cOli3)*T(i4, cOli4)*T(i1, cOli1)*T(i2, cOli2));
*    id T(i1, cOli1?)*T(i2, cOli2?)*T(i3, cOli3?)*cOlf(cOli1?, cOli2?, cOli3?) = Tsym1(i1, i2, i3) - cOlf(cOli1, cOli2, cOli3)*(T(i1, cOli1)*T(i3, cOli3)*T(i2, cOli2)
*                                                                                                                             + T(i2, cOli2)*T(i3, cOli3)*T(i1, cOli1)
*                                                                                                                             + T(i3, cOli3)*T(i2, cOli2)*T(i1, cOli1));
  #endif
* jacobi
  id cOlf(N1_?,N4_?,N5_?)*cOlf(N2_?,N3_?,N5_?) = - cOlf(N1_?,N2_?,N5_?)*cOlf(N3_?,N4_?,N5_?) - cOlf(N1_?,N3_?,N5_?)*cOlf(N4_?,N2_?,N5_?);
  repeat;
    #do k=1,5
      #do j=`k'+1,6
        id once T(i`j', cOli1?)*T(i`k', cOli2?)*summe(i`j')*summe(i`k') = T(i`k', cOli2)*T(i`j', cOli1)*summe(i`k')*summe(i`j') + i_*cOlf(cOli1, cOli2, cOli3)*T(i`k', cOli3)*summe(i`k')*replace_(i`j', i`k', pi`j', pi`k');
        sum cOli3;
        id once T(i`j', cOli1?)*T(i`k', cOli2?)*summeP(i`j')*summeP(i`k') = T(i`k', cOli2)*T(i`j', cOli1)*summeP(i`k')*summeP(i`j');
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
    id T(i1?, cOli1?)*T(i1?, cOli2?)*cOlf(cOli1?, cOli2?, cOli3?) = 1/2*i_*cOlf(cOli1, cOli2, cOli4)*T(i1, cOli4)*cOlf(cOli1, cOli2, cOli3);
    id Tsym1(i2, i3, i1) = -Tsym1(i1, i3, i2)*replace_(i2, i3, i3, i2, pi2, pi3, pi3, pi2);
    sum cOli4;
  endrepeat;
  .sort
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
  #do i=1,3
    #do j=`i'+1,4
      id gamma1(pi`j', pi`i', q1, mu1?) = -gamma1(pi`i', pi`j', q1, mu1);
    #enddo
  #enddo
  #if `l'==2
    id T(i1, cOli1?)*T(i2, cOli2?)*T(i2, cOli3?) = T(i1, cOli1)*T(i2, cOli2)*T(i2, cOli3)*replace_(i1, i2, i2, i1, pi2, pi1, pi1, pi2);
  #endif
  id gamma1(pi1?, pi1?, q1, mu1?) = 0;
*  id gamma1(pi1,pi2,q1,mu1?)*gamma1(pi2,pi3,q1,mu2?) = gamma1(pi1,pi2,q1,mu1)*gamma1(pi2,pi3,q1,mu2)*replace_(i2, i1, i1, i2, pi2, pi1, pi1, pi2);
*  id phase(i3,i4)*phaseConj(i1,i2) = phase(i3,i4)*phaseConj(i1,i2)*replace_(i1, i3, i3, i1, i2, i4, i4, i2, pi1, pi3, pi3, pi1, pi2, pi4, pi4, pi2);
#enddo
id phase(i1?, i2?)*phaseConj(i1?, i2?) = 1;

b T, cOlf, summe, summeP, Tsym1, Tsym2;
print+s interference, square;

.end