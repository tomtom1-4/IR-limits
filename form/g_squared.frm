#-
#: IncDir procedures
#: SmallExtension 200M
#: WorkSpace      4G
#: MaxTermSize 4M
Off Statistics;

#include declarations.h
auto i mui;
auto v pi;
cf gamma1, j1, gamma2Dipole, gamma2Tripole, gamma2DipoleConj, gamma2TripoleConj;
cf phase(symmetric), phaseConj(symmetric);
f Tsym1, Tsym2;
cf C, summeP;
Format 160;

l J0 = summe(i1)*T(i1, c1)*j1(pi1, q1, mu1)*E1(mu1);

l J1 = summe(i1)*summe(i2)*i_*cOlf(c1, cOli1, cOli2)*T(i1, cOli1)*T(i2, cOli2)*gamma1(pi1, pi2, q1, mu1)*E1(mu1)*phase(i1, i2);

l J2 = summeP(i1)*summeP(i2)*i_*cOlf(c1, cOli1, cOli2)*T(i1, cOli1)*T(i2, cOli2)*gamma2Dipole(pi1, pi2, q1, mu1)*E1(mu1)*CA
  + summeP(i1)*summeP(i2)*summeP(i3)*cOlf(c1, cOli1, cOli2)*cOlf(cOli2, cOli3, cOli4)*T(i1, cOli3)*T(i2, cOli4)*T(i3, cOli1)*gamma2Tripole(pi1, pi2, pi3, q1, mu1)*E1(mu1);

sum mu1, cOli1,cOli2;
.sort
* The color factors are symmetrized so an explicit reversal is not necessary
l J0Conj = J0*replace_(i_, -i_, i1, i4, i2, i5, i3, i6, pi1, pi4, pi2, pi5, pi3, pi6, gamma2Dipole, gamma2DipoleConj, gamma2Tripole, gamma2TripoleConj, phase, phaseConj);
l J1Conj = J1*replace_(i_, -i_, i1, i4, i2, i5, i3, i6, pi1, pi4, pi2, pi5, pi3, pi6, gamma2Dipole, gamma2DipoleConj, gamma2Tripole, gamma2TripoleConj, phase, phaseConj);
l J2Conj = J2*replace_(i_, -i_, i1, i4, i2, i5, i3, i6, pi1, pi4, pi2, pi5, pi3, pi6, gamma2Dipole, gamma2DipoleConj, gamma2Tripole, gamma2TripoleConj, phase, phaseConj);
.sort
* rename color indices
InExpression J0Conj, J1Conj, J2Conj;
  mul replace_(N1_?, cOli1, N2_?, cOli2, N3_?, cOli3, N4_?, cOli4, N5_?, cOli5, N6_?, cOli6, N7_?, cOli7, N8_?, cOli8, N9_?, cOli9, N10_?, cOli10);
EndInExpression;
.sort

l OneLoop = J1Conj*J0 + J0Conj*J1;
l TwoLoop = J1Conj*J1 + J2Conj*J0 + J0Conj*J2;

l TwoLoopTest = (1 + replace_(i_, -i_, phase, phaseConj, phaseConj, phase, gamma2Dipole, gamma2DipoleConj, gamma2Tripole, gamma2TripoleConj))*
              ((T(i1, cOli1)*(T(i2, cOli2)*T(i3, cOli3) + T(i3, cOli3)*T(i2, cOli2))*T(i4, cOli4)
               + T(i4, cOli4)*(T(i3, cOli3)*T(i2, cOli2) + T(i2, cOli2)*T(i3, cOli3))*T(i1, cOli1))
                *cOlf(cOli1, cOli4, cOli5)*cOlf(cOli2, cOli3, cOli5)
        *(summe(i1)*summeP(i2)*summeP(i3)*summeP(i4)*(-j1(pi1, q1, mu1)*gamma2Tripole(pi2, pi3, pi4, q1, mu1)/4)
         +summe(i1)*summe(i2)*summe(i3)*summe(i4)*(-gamma1(pi1, pi4, q1, mu1)*gamma1(pi2, pi3, q1, mu1)*phase(i1, i4)*phaseConj(i2, i3)/8))
    +CA*cOlf(cOli1, cOli2, cOli3)*Tc(i1, cOli1)*Tc(i2, cOli2)*Tc(i3, cOli3)*summeP(i1)*summeP(i2)*summeP(i3)
        *(gamma1(pi1, pi2, q1, mu1)*gamma1(pi1, pi3, q1, mu1)*phase(i1, i2)*phaseConj(i1, i3)*i_/2
          -j1(pi1,q1,N4_?)*gamma2Tripole(pi1,pi2,pi3,q1,N4_?)*i_/4
          -j1(pi1,q1,N4_?)*gamma2Tripole(pi2,pi1,pi3,q1,N4_?)*i_/4
          -j1(pi1,q1,N4_?)*gamma2Dipole(pi2,pi3,q1,N4_?)*i_)
    +Tc(i1, cOli1)*Tc(i2, cOli1)*summeP(i1)*summeP(i2)*CA^2
        *(gamma1(pi1, pi2, q1, N2_?)^2/4
         +j1(pi1,q1,N2_?)*gamma2Dipole(pi1,pi2,q1,N2_?)) );
l TwoLoopControl = TwoLoop - TwoLoopTest;
sum mu1, cOli1,...,cOli10;
.sort
* Use spin sums to replace Polarization vectors
ct EE1;
drop J0, J1, J0Conj, J1Conj, J2, J2Conj;
toTensor functions E1, EE1;
id EE1(mu1?, mu2?) = -d_(mu1, mu2);
sum c1, c2, c3;
.sort
renumber 1;
.sort
id gamma1(pi1, pi2, q1, mu1?)*j1(pi3, q1, mu2?) = gamma1(pi1, pi2, q1, mu1)*j1(pi3, q1, mu2)*replace_(i3, i1, i1, i2, i2, i3, pi3, pi1, pi1, pi2, pi2, pi3);
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

#do i=1,4
  #do j=`i'+1,5
    if(match(summe(i`i')*summeP(i`j')));
      mul replace_(i`j', i`i', i`i', i`j', pi`j', pi`i', pi`i', pi`j');
      redefine j "10";
    endif;
    .sort
  #enddo
#enddo

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
if(count(summeP, 1)==4)
  id summe(i5)*summeP(i1)*summeP(i2)*summeP(i3)*summeP(i4) = summeP(i1)*summeP(i2)*summeP(i3)*summeP(i4)*summeP(i5)
                                                  + summeP(i1)*summeP(i2)*summeP(i3)*summeP(i4)*(replace_(i5, i1, pi5, pi1)
                                                                                               + replace_(i5, i2, pi5, pi2)
                                                                                               + replace_(i5, i3, pi5, pi3)
                                                                                               + replace_(i5, i4, pi5, pi4));
#enddo
.sort
#do i=1,6
  id once T(i1?, cOli1?) = Tc(i1, cOli1);
  repeat;
    id once Tc(i1?, ?args)*T(i1?, cOli1?) = Tc(i1, ?args, cOli1);
  endrepeat;
#enddo

#do l=1,4
  #if `l'==1
*    id T(i1, cOli1?)*T(i2, cOli2?)*T(i3, cOli3?)*T(i4, cOli4?)*cOlf(cOli1?, cOli2?, cOli5?)*cOlf(cOli3?, cOli4?, cOli5?)
*        = Tsym2(i1, i2, i3, i4) - cOlf(cOli1, cOli2, cOli5)*cOlf(cOli3, cOli4, cOli5)*(T(i3, cOli3)*T(i4, cOli4)*T(i1, cOli1)*T(i2, cOli2));
*    id T(i1, cOli1?)*T(i2, cOli2?)*T(i3, cOli3?)*cOlf(cOli1?, cOli2?, cOli3?) = Tsym1(i1, i2, i3) - cOlf(cOli1, cOli2, cOli3)*(T(i1, cOli1)*T(i3, cOli3)*T(i2, cOli2)
*                                                                                                                             + T(i2, cOli2)*T(i3, cOli3)*T(i1, cOli1)
*                                                                                                                             + T(i3, cOli3)*T(i2, cOli2)*T(i1, cOli1));
  #endif
* jacobi
  if(match(Tc(i2, cOli1?, cOli2?)));
    id cOlf(N1_?,N2_?,N5_?)*cOlf(N3_?,N4_?,N5_?) = - cOlf(N1_?,N4_?,N5_?)*cOlf(N2_?,N3_?,N5_?) - cOlf(N1_?,N3_?,N5_?)*cOlf(N4_?,N2_?,N5_?);
  else;
    id cOlf(N1_?,N4_?,N5_?)*cOlf(N2_?,N3_?,N5_?) = - cOlf(N1_?,N2_?,N5_?)*cOlf(N3_?,N4_?,N5_?) - cOlf(N1_?,N3_?,N5_?)*cOlf(N4_?,N2_?,N5_?);
  endif;
  repeat;
    id cOlf(cOli1?, cOli2?, cOli3?)*cOlf(cOli1?, cOli2?, cOli4?) = CA*delta(cOli3, cOli4);
    id cOlf(cOli1?, cOli3?, cOli4?)*cOlf(cOli2?, cOli5?, cOli6?)*cOlf(cOli3?, cOli5?, cOli7?)*cOlf(cOli4?, cOli6?, cOli7?) = CA^2/2*delta(cOli1, cOli2);
    id cOlf(cOli1?,cOli4?,cOli5?)*cOlf(cOli2?,cOli4?,cOli6?)*cOlf(cOli3?,cOli5?,cOli6?) = CA/2*cOlf(cOli1, cOli2, cOli3);
    id Tc(i1?, cOli1?)*delta(cOli1?, cOli2?) = Tc(i1, cOli2);
    id cOlf(cOli1?, cOli2?, cOli3?)*delta(cOli3?, cOli4?) = cOlf(cOli1, cOli2, cOli4);
    id Tc(i1?, cOli1?, cOli2?)*cOlf(cOli1?, cOli2?, cOli3?) = 1/2*i_*cOlf(cOli1, cOli2, cOli4)*Tc(i1, cOli4)*cOlf(cOli1, cOli2, cOli3);
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
  #do i=2, 4
    id j1(pi`i', q1, mu1?) = j1(pi`i', q1, mu1)*replace_(pi`i', pi1, pi1, pi`i', i`i', i1, i1, i`i');
  #enddo
  id gamma1(pi1,pi4,q1,mu1?)*gamma1(pi2,pi3,q1,mu1?) = gamma1(pi1,pi4,q1,mu1)*gamma1(pi2,pi3,q1,mu1)*replace_(i4, i2, i2, i3, i3, i4, pi4, pi2, pi2, pi3, pi3, pi4);
  id phase(i1,i3)*phaseConj(i1,i2) = phase(i1,i3)*phaseConj(i1,i2)*replace_(i2, i3, i3, i2, pi2, pi3, pi3, pi2);
  id phase(i2,i3)*phaseConj(i1,i2) = phase(i2,i3)*phaseConj(i1,i2)*replace_(i1, i3, i3, i1, pi1, pi3, pi3, pi1);
  id phase(i3,i4)*phaseConj(i1,i2) = phase(i3,i4)*phaseConj(i1,i2)*replace_(i1, i3, i3, i1, i2, i4, i4, i2, pi1, pi3, pi3, pi1, pi2, pi4, pi4, pi2);
  id phase(i1,i2)*phaseConj(i2,i3) = phase(i1,i2)*phaseConj(i2,i3)*replace_(i2, i1, i1, i2, pi2, pi1, pi1, pi2);
  id phase(i2,i3)*phaseConj(i1,i3) = phase(i2,i3)*phaseConj(i1,i3)*replace_(i3, i1, i1, i2, i2, i3, pi3, pi1, pi1, pi2, pi2,pi3);
  id phaseConj(i2,i3)*phase(i1,i3) = phaseConj(i2,i3)*phase(i1,i3)*replace_(i3, i1, i1, i2, i2, i3, pi3, pi1, pi1, pi2, pi2,pi3);
  #do i=1,3
    #do j=`i'+1,4
      id gamma1(pi`j', pi`i', q1, mu1?) = -gamma1(pi`i', pi`j', q1, mu1);
      id gamma2Dipole(pi`j', pi`i', q1, mu1?) = - gamma2Dipole(pi`i', pi`j', q1, mu1);
      id gamma2DipoleConj(pi`j', pi`i', q1, mu1?) = - gamma2DipoleConj(pi`i', pi`j', q1, mu1);
      id gamma2Tripole(pi`j', pi10?, pi`i', q1, mu1?) = - gamma2Tripole(pi`i', pi10, pi`j', q1, mu1);
      id gamma2TripoleConj(pi`j', pi10?, pi`i', q1, mu1?) = - gamma2TripoleConj(pi`i', pi10, pi`j', q1, mu1);
    #enddo
  #enddo

  id gamma1(pi1?, pi1?, q1, mu1?) = 0;
*  id gamma1(pi1,pi2,q1,mu1?)*gamma1(pi2,pi3,q1,mu2?) = gamma1(pi1,pi2,q1,mu1)*gamma1(pi2,pi3,q1,mu2)*replace_(i2, i1, i1, i2, pi2, pi1, pi1, pi2);
*  id phase(i3,i4)*phaseConj(i1,i2) = phase(i3,i4)*phaseConj(i1,i2)*replace_(i1, i3, i3, i1, i2, i4, i4, i2, pi1, pi3, pi3, pi1, pi2, pi4, pi4, pi2);
#enddo
id phase(i1?, i2?)*phaseConj(i1?, i2?) = 1;

b Tc, cOlf, summe, summeP, Tsym1, Tsym2;
print+s ;

.end