#-
#: IncDir procedures
#: SmallExtension 200M
#: WorkSpace      4G
#: MaxTermSize 4M
Off Statistics;

#include declarations.h
auto i mui;
auto v pi;
cf summeP;
cf Jqq, JqqConj, JqqReal, JqqImag;


l J1Squared = -summe(i1)*summeP(i2)*summeP(i3)*den(q1.q2)^2/2*den(pi1.q1 + pi1.q2)
  *(cOlTr(cOli1, cOli2, cOli3)*T(i1, cOli1)*T(i2, cOli2)*T(i3, cOli3) + cOlTr(cOli3, cOli2, cOli1)*T(i3, cOli3)*T(i2, cOli2)*T(i1, cOli1)*replace_(Jqq, JqqConj))
  *(g_(0, q2, pi1, q1, pi2)*Jqq(q1, q2, pi2, pi3)*den(pi2.q1 + pi2.q2) - g_(0, q2, pi1, q1, pi3)*Jqq(q2, q1, pi3, pi2)*den(pi3.q1 + pi3.q2));
tracen 0;

sum cOli1,...,cOli3;
id summe(i1)*summeP(i2)*summeP(i3) = summeP(i1)*summeP(i2)*summeP(i3) + summeP(i2)*summeP(i3)*(replace_(i1, i3, pi1, pi3) + replace_(i1, i2, pi1, pi2));
if(occurs(i1) == 0) mul replace_(i2, i1, pi2, pi1, i3, i2, pi3, pi2);

#do i=1,6
  id once T(i1?, cOli1?) = Tc(i1, cOli1);
  repeat;
    id once Tc(i1?, ?args)*T(i1?, cOli1?) = Tc(i1, ?args, cOli1);
  endrepeat;
#enddo
.sort
id Tc(i1,cOli1?,cOli2?) = Tc(i1,cOli1,cOli2)*replace_(i1, i2, i2, i1, pi1, pi2, pi2, pi1);
*id Tc(i1?, cOli1?, cOli2?)*cOlf(cOli1?, cOli2?, cOli3?) = i_/2*CA*Tc(i1, cOli3);
id cOlTr(cOli1?, cOli2?, cOli3?) = TF/2*(dsym(cOli1, cOli2, cOli3) + i_*cOlf(cOli1, cOli2, cOli3));
.sort
id delta(cOli1?, cOli2?)*dsym(cOli1?, cOli3?, cOli4?) = dsym(cOli2, cOli3, cOli4);

repeat;
  id once dsym(cOli1?, cOli2?, cOli3?) = 1/TF*(cOlT(cOli10, cOli11, cOli1)*cOlT(cOli11, cOli12, cOli2)*cOlT(cOli12, cOli10, cOli3)
                                        + cOlT(cOli10, cOli11, cOli2)*cOlT(cOli11, cOli12, cOli1)*cOlT(cOli12, cOli10, cOli3));
  sum cOli10,...,cOli12;
endrepeat;
#call Cvitanovic
.sort
id delta(cOli1?, cOli2?)*Tc(i2?, cOli2?) = Tc(i2, cOli1);
id cOlTr(cOli1?, cOli2?, cOli3?) = TF/2*(dsym(cOli1, cOli2, cOli3) + i_*cOlf(cOli1, cOli2, cOli3));
id Tc(i1?, cOli1?, cOli2?)*cOlf(cOli1?, cOli2?, cOli3?) = i_/2*CA*Tc(i1, cOli3);
.sort
id CA = 2*TF*NF;
id Jqq(q2,q1,pi3,pi2) = Jqq(q2,q1,pi3,pi2)*replace_(i3, i2, i2, i3, pi3, pi2, pi2, pi3);
id JqqConj(q2,q1,pi3,pi2) = JqqConj(q2,q1,pi3,pi2)*replace_(i3, i2, i2, i3, pi3, pi2, pi2, pi3);
*if(occurs(cOlf)==0);
*  id summeP(i1)*summeP(i2)*summeP(i3) = summeP(i1)*summeP(i2)*(summe(i3) - replace_(i3, i1, pi3, pi1) - replace_(i3, i2, pi3, pi2));
*  id summeP(i1)*summeP(i2) = summe(i1)*(summe(i2) - replace_(i2, i1, pi2, pi1));
*endif;
*if(occurs(i3) && (occurs(i2) == 0)) mul replace_(i3, i2, pi3, pi2);
*repeat;
*  id once Tc(i1?, cOli1?, ?args)*Tc(i1?, cOli2?) = Tc(i1, cOli1, ?args, cOli2);
*endrepeat;
.sort
s  kinfactor;
id q1.q2*pi1.pi2 = kinfactor + q1.pi1*q2.pi2 + q1.pi2*q2.pi1;
id pi1?.pi1? = 0;

if(occurs(i_)==0) id JqqConj(?args) = 2*JqqReal(?args) - Jqq(?args);
if(occurs(i_)) id JqqConj(?args) = -2*JqqImag(?args)*i_ + Jqq(?args);
id JqqReal(q1,q2,pi2,pi1) = JqqReal(q1,q2,pi2,pi1)*replace_(i1, i2, i2, i1, pi1, pi2, pi2, pi1);
*id Tc(i1,N1_?)*Tc(i2,N2_?,N3_?) = -Tc(i1,N1_?)*Tc(i2,N2_?)*Tc(i3, N3_?)*summeP(i3) - Tc(i1, N1_?, N3_?)*Tc(i2, N2_?);


b cOlTr, T, Tc, summeP, summe, dsym, cOlf;*, Jqq, JqqConj, JqqReal, JqqImag;
print+s ;
.end