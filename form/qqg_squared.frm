#-
#: IncDir procedures
#: SmallExtension 100M
#: WorkSpace      2G
*#: MaxTermSize 2M
Off Statistics;

#include declarations.h
s colFac; * = TF/2*(CF - TF/dF) = 7/24
v q13, q12, q23, q123;
s marker;

l Sij1 = (pi.q2/4*den(q2.q3*q2.q3*pi.q1*(pi.q2 + pi.q3))*(pi.pj*den(pj.q1)*(3*pj.q3*den(pj.q2 + pj.q3) - 2*pi.q3*den(pi.q2 + pi.q3)) - 2*pi.pi*pj.q3*den(pi.q1*(pj.q2 + pj.q3)))
  + pi.pj/8*den(q2.q3*pi.q1*(pi.q2 + pi.q3))*(-3*pi.pj*den(pj.q1*(pj.q2 + pj.q3)) + 2*pi.pi*(den(pi.q1*(pj.q2 + pj.q3)) + den(pj.q1*(pi.q2 + pi.q3)))))*(1 + replace_(q2, q3, q3, q2));

l Sija = den(pi.q1*(pj.q2 + pj.q3)*(pa.q1 + pa.q2 + pa.q3)*2*(q1.q2 + q1.q3 + q2.q3))*(den(q1.q2*q2.q3)*(pi.q1*(pj.q2*pa.q3 + pj.q3*pa.q2)
  + pi.q2*(pj.q1*pa.q3 + pj.q3*pa.q1 + 2*pj.q2*pa.q3 + 2*pj.q3*pa.q2)
  + pi.q3*(-pj.q1*pa.q2 + pj.q2*pa.q1) + q1.q3*(pi.pj*pa.q2 - pi.pa*pj.q2 - pj.pa*pi.q2))
  + den(q1.q2)*(pi.pa*pj.q1 - pi.pj*pa.q1 - pj.pa*(pi.q1 + 2*pi.q2))
  + den(q2.q3)*(pj.pa*pi.q3 - pi.pj*pa.q3 - pi.pa*pj.q3))*(1 - 1*replace_(q2, q3, q3, q2));

l Sij2 = den(2*q2.q3*(pi.q1 + pi.q2 + pi.q3))*(2*den(2*(q1.q2 + q1.q3 + q2.q3)*q2.q3)*(pj.q2*(2*mi^2*q1.q3 - (pi.q1 + pi.q2 + pi.q3)*pi.q3)*den(pi.q1*(pj.q2 + pj.q3))
    - pi.q2*(2*pi.pj*q1.q3 - (pi.q1 - pi.q2 + 3*pi.q3)*pj.q3)*den(pj.q1*(pi.q2 + pi.q3)))
  + den(2*(q1.q2 + q1.q3 + q2.q3)*q1.q2)*((2*(pi.q1 + pi.q2)*(pi.pj*q2.q3 - pi.q2*pj.q3 - pi.q3*pj.q2) + mi^2*(pj.q2*q1.q3 - pj.q1*q2.q3))*den(pi.q1*(pj.q2 + pj.q3))
    - (mi^2*((pj.q1 + 2*pj.q2)*q2.q3 + pj.q2*q1.q3) - 2*(pi.q2*pj.q1 + (pi.q1 + 2*pi.q2)*pj.q2)*pi.q3)*den(pj.q1*(pi.q2 + pi.q3)))
  + den(2*(q1.q2 + q1.q3 + q2.q3))*((pi.pj*(pi.q1 + pi.q2 + pi.q3) + mi^2*(pj.q2 - 2*pj.q1))*den(pi.q1*(pj.q2 + pj.q3)) - (3*mi^2*pj.q2 - pi.pj*pi.q1)*den(pj.q1*(pi.q2 + pi.q3)))
  + den(q2.q3)*(pi.q2*(den(pi.q1) - den(pi.q2 + pi.q3))*(mi^2*pj.q3*den(pi.q1*(pj.q2 + pj.q3)) - pi.pj*pi.q3*den(pj.q1*(pi.q2 + pi.q3))))
  + 1/2*mi^2*pi.pj*(den(pi.q2 + pi.q3) - den(pi.q1))*(den(pi.q1*(pj.q2 + pj.q3)) - den(pj.q1*(pi.q2 + pi.q3))))
  *(1 + 1*replace_(q2, q3, q3, q2));

l Sij3NAB = den((pi.q1 + pi.q2 + pi.q3)*(pj.q1 + pj.q2 + pj.q3))*(den((2*(q1.q2 + q1.q3 + q2.q3))^2)*(8*pi.pj*q1.q2*q1.q3*den(q2.q3^2)
  + (2*((d - 6)*(pi.q1 - pi.q3) - (2 + d)*pi.q2)*pj.q2*q1.q3)*den(q1.q2)*den(q2.q3) + q2.q3*den(q1.q2*q1.q3)*(-2*pi.pj*q2.q3
  + pi.q1*((4 - d)*pj.q1 + 2*pj.q2) + 2*pi.q2*(pj.q1 + 2*pj.q3)) + den(q2.q3)*(8*pi.pj*q1.q3 + pi.q1*((2 - d)*pj.q1
  + 2*(d - 6)*pj.q2) + 8*pi.q2*(pj.q3 - pj.q2)) + den(q1.q2)*(pi.pj*((d - 2)*q1.q3 - 4*q2.q3) + pi.q1*(2*(6 - d)*pj.q1
  + 2*(d - 4)*pj.q2 + (10 - d)*pj.q3) + 2*pi.q2*((d + 2)*pj.q1 - 2*pj.q2 + 6*pj.q3)
  + pi.q3*((d - 2)*pj.q1 + 8*pj.q2))
  + (4 - d)*pi.pj) + den(2*(q1.q2 + q1.q3 + q2.q3))*(den(q2.q3^2)*((den(pi.q1) - den(pi.q2 + pi.q3))*2*pi.q2*(2*pi.pj*q1.q3 - 2*pi.q2*pj.q3 - pi.q3*pj.q1))
  + den(q1.q2*q2.q3)*(den(pi.q1) - den(pi.q2 + pi.q3))*((2*pi.pj*pi.q2 - mi^2*pj.q2)*q1.q3 - 2*pi.q2*((pi.q1 + pi.q2)*pj.q3 + pi.q3*(pj.q1 + pj.q2)))
  + den(q2.q3)*((den(pi.q1) - den(pi.q2 + pi.q3))*(mi^2*(pj.q1 - pj.q2) + 2*pi.pj*(2*pi.q2 - pi.q1)))
  + den(q1.q2)*((den(pi.q1) - den(pi.q2 + pi.q3))*(2*pi.pj*pi.q2 + mi^2*pj.q1) - 4*pi.pj))
  + pi.pj*pi.q2*pj.q3*den(2*q2.q3^2)*(den(pi.q1) - den(pi.q2 + pi.q3))*(den(pj.q1) - den(pj.q2 + pj.q3))
  - pi.pj^2*den(4*q2.q3)*(den(pi.q1) - den(pi.q2 + pi.q3))*(den(pj.q1) - den(pj.q2 + pj.q3))
  )*(1 + replace_(q2, q3, q3, q2))*(1 + replace_(pi, pj, pj, pi));

l Sij3AB = den(4*(q1.q2 + q1.q3 + q2.q3)^2)*den((pi.q1 + pi.q2 + pi.q3)*(pj.q1 + pj.q2 + pj.q3))*(q2.q3*den(q1.q2*q1.q3)*(2*pi.pj*q2.q3 + (d - 4)*pi.q1*pj.q1
  - 2*pi.q1*(pj.q2 + pj.q3) - 4*pi.q2*pj.q3) + den(q1.q2)*(pi.pj*((d - 2)*q1.q3 + 4*q2.q3) + 2*(4 - d)*pi.q1*pj.q2
  + 2*(2 - d)*pi.q1*pj.q3 + 4*pi.q2*(pj.q2 - pj.q3)) + (d - 4)*pi.pj)*(1 + replace_(q2, q3, q3, q2))*(1 + replace_(pi, pj, pj, pi));


*#call cppFormat
*brackets den;
*.sort
*keep Brackets;
*format Mathematica;
*#write <results/gqq.m> "den[x_]:= 1/x;";
*#write <results/gqq.m> "Sgqqij1:= (%+E);", Sij1;
*#write <results/gqq.m> "Sgqqij2 = (%+E);", Sij2;
*#write <results/gqq.m> "Sgqqij3AB = (%+E);", Sij3AB;
*#write <results/gqq.m> "Sgqqij3NAB = (%+E);", Sij3NAB;
*#write <results/gqq.m> "Sgqqija = (%+E);", Sija;
*#write <results/gqq.m> "Sgqq = Sgij Sqqab (Ti.Tj * Ta.Tb + Ta.Tb * Ti.Tj) - TF*Ti.Tj*(CA*(Sgqqij1 + Sgqqij2 + 1/4*(Sgqqij3NAB - Sgqqij3AB)) + CF*Sgqqij3AB) + dsym*Ti*Tj*Ta*Sgqqija;"

.sort

l reducible = TF*summe(i)*summe(j)*summe(a)*summe(b)*(-pi.pb*den(pi.q1*pb.q1)*g_(0,q2, pj, q3, pa)*den(pj.q2 + pj.q3)*den(pa.q2 + pa.q3)*den(4*(q2.q3)^2))*1/4*(T(i,cOli1)*T(j,cOli2)*T(a,cOli2)*T(b,cOli1) + T(j,cOli2)*T(i,cOli1)*T(b,cOli1)*T(a,cOli2) + T(i,cOli1)*T(j,cOli2)*T(b,cOli1)*T(a,cOli2) + T(j,cOli2)*T(i,cOli1)*T(a,cOli2)*T(b,cOli1));

l reducibleCatani = TF*summe(i)*summe(j)*summe(a)*summe(b)*(-pi.pb*den(pi.q1*pb.q1)*(pj.q2*pa.q3 + pj.q3*pa.q2 - pj.pa*q2.q3)*den(pj.q2 + pj.q3)*den(pa.q2 + pa.q3)*den(q2.q3^2))*1/2*(T(i, cOli1)*T(b, cOli1)*T(j, cOli2)*T(a, cOli2) + T(j, cOli2)*T(a, cOli2)*T(i, cOli1)*T(b, cOli1))
  - TF*CA*summe(i)*summe(j)*T(i, cOli1)*T(j, cOli1)*Sij1;

l reducibleDifference = reducible - reducibleCatani;

l interference = 1/2*summe(i)*summe(j)*summe(a)*1/2*TF*dsym(cOli1, cOli2, cOli3)*((T(i, cOli1)*T(j, cOli2) + T(j, cOli2)*T(i, cOli1))*T(a, cOli3) + T(a, cOli3)*(T(i, cOli1)*T(j, cOli2) + T(j, cOli2)*T(i, cOli1)))
  *den(2*(q1.q2 + q1.q3 + q2.q3)*(pa.q1 + pa.q2 + pa.q3))*(g_(0, q2, pa, q13, pi, q3, pj)*den(2*q1.q3) - g_(0, q2, pi, q12, pa, q3, pj)*den(2*q1.q2))
  /s23*den(pj.q2 + pj.q3)*den(pi.q1)*replace_(q13, q1 + q3, q12, q1 + q2, q23, q2 + q3)
  - 1/2*summe(i)*summe(j)*summe(a)*TF*1/2*i_*cOlf(cOli1, cOli2, cOli3)*((T(i, cOli1)*T(j, cOli2) + T(j, cOli2)*T(i, cOli1))*T(a, cOli3) - T(a, cOli3)*(T(i, cOli1)*T(j, cOli2) + T(j, cOli2)*T(i, cOli1)))
  *den(pa.q1 + pa.q2 + pa.q3)*den(2*q2.q3)*den(pi.q1)*den(pj.q2 + pj.q3)*(g_(0, q2, pa, q3, pj)*den(2*q2.q3)*pa.pi*(den(pa.q1) - den(pa.q2 + pa.q3))
    + den(2*(q1.q2 + q1.q3 + q2.q3))*(den(2*q2.q3)*(2*(pa.q2 + pa.q3 - pa.q1)*g_(0, q2, pi, q3, pj) - 4*(pi.q2 + pi.q3)*g_(0, q2, pa, q3, pj) + 4*pa.pi*g_(0, q2, q1, q3, pj))
      - den(2*q1.q2)*g_(0, q2, pi, q12, pa, q3, pj) - den(2*q1.q3)*g_(0, q2, pa, q13, pi, q3, pj)))*replace_(q12, q1 + q2, q13, q1 + q3, q23, q2 + q3);

l interferenceCatani = - summe(i)*summe(j)*summe(a)*TF*dsym(cOli1, cOli2, cOli3)*T(i, cOli1)*T(j, cOli2)*T(a, cOli3)*Sija
                       - summe(i)*summe(j)*TF*CA*T(i, cOli1)*T(j, cOli1)*Sij2;

l interferenceDifference = interference - interferenceCatani;

l irreducibleNAB = summe(i)*summe(j)*1/2*i_*cOlf(cOli1, cOli2, cOli3)*T(i, cOli2)*1/2*(-i_)*cOlf(cOli1, cOli4, cOli3)*T(j, cOli4)*TF
  *g_(0, q2)*(den(pi.q1 + pi.q2 + pi.q3)*(pi(mu)*den(2*q2.q3)*(den(pi.q1) - den(pi.q2 + pi.q3))*g_(0, pi)
    + den(2*q1.q2 + 2*q1.q3 + 2*q2.q3)*(den(2*q2.q3)*(2*(pi.q2 + pi.q3 - pi.q1)*g_(0, mu) - 4*(q2(mu) + q3(mu))*g_(0, pi) + 4*pi(mu)*g_(0, q1)) - g_(0, mu, q12, pi)*den(2*q1.q2) - g_(0, pi, q13, mu)*den(2*q1.q3))))
  *g_(0, q3)*(den(pj.q1 + pj.q2 + pj.q3)*(pj(mu)*den(2*q2.q3)*(den(pj.q1) - den(pj.q2 + pj.q3))*g_(0, pj)
    + den(2*q1.q2 + 2*q1.q3 + 2*q2.q3)*(den(2*q2.q3)*(2*(pj.q2 + pj.q3 - pj.q1)*g_(0, mu) - 4*(q2(mu) + q3(mu))*g_(0, pj) + 4*pj(mu)*g_(0, q1)) - g_(0, pj, q12, mu)*den(2*q1.q2) - g_(0, mu, q13, pj)*den(2*q1.q3))))
  *(-1)*replace_(q12, q1 + q2, q13, q1 + q3);

l irreducibleNABCatani = -summe(i)*summe(j)*TF*CA/4*T(i, cOli1)*T(j, cOli1)*Sij3NAB;

l irreducibleNABDifference = irreducibleNAB - irreducibleNABCatani;

l irreducibleAB = replace_(q13, q1 + q3)*replace_(q12, q1 + q2)*summe(i)*summe(j)*colFac*T(i, cOli1)*T(j, cOli1)*(-1)*den(4*(q1.q2 + q1.q3 + q2.q3)^2*(pi.q1 + pi.q2 + pi.q3)*(pj.q1 + pj.q2 + pj.q3))
  *g_(0, q2)*(g_(0, pi, q13, mu)*den(2*q1.q3) - g_(0, mu, q12, pi)*den(2*q1.q2))*g_(0, q3)*(g_(0, mu, q13, pj)*den(2*q1.q3) - g_(0, pj, q12, mu)*den(2*q1.q2));

l irreducibleABCatani = summe(i)*summe(j)*T(i, cOli1)*T(j, cOli1)*(-colFac)*Sij3AB;

l irreducibleABDifference = irreducibleAB - irreducibleABCatani;
.sort
drop Sij1, Sija, Sij2, Sij3AB, Sij3NAB;
drop reducible, reducibleCatani; *works
drop interference, interferenceCatani; *works
drop irreducibleNAB, irreducibleNABCatani; *works
drop irreducibleAB, irreducibleABCatani; *works

tracen 0;
sum cOli1,...,cOli4, mu;
.sort
*drop Sij1;

#do i=1,5
repeat;
  id T(j, cOli2?)*T(i, cOli1?)*summe(i)*summe(j) =  T(i, cOli1)*T(j, cOli2)*summe(i)*summe(j) + i_*cOlf(cOli2, cOli1, cOli10)*T(i, cOli10)*summe(i)*replace_(j, i, pj, pi, sj1, si1, sj2, si2, sj3, si3, sja, sia, sjb, sib, sij, 4*mi^2, mj, mi);
  sum cOli10;
  id T(a, cOli2?)*T(i, cOli1?)*summe(i)*summe(a) =  T(i, cOli1)*T(a, cOli2)*summe(i)*summe(a) + i_*cOlf(cOli2, cOli1, cOli11)*T(i, cOli11)*summe(i)*replace_(a, i, pa, pi, sa1, si1, sa2, si2, sa3, si3, sia, 4*mi^2, sja, sij, sab, sib, ma, mi);
  sum cOli11;
endrepeat;
.sort
repeat;
  id T(j, cOli2?)*T(b, cOli1?)*summe(b)*summe(j) =  T(b, cOli1)*T(j, cOli2)*summe(b)*summe(j) + i_*cOlf(cOli2, cOli1, cOli10)*T(j, cOli10)*summe(j)*replace_(b, j, pb, pj, sb1, sj1, sb2, sj2, sb3, sj3, sjb, 4*mj^2, sab, sja, sib, sij, mb, mj);
  sum cOli10;
  id T(a, cOli2?)*T(b, cOli1?)*summe(b)*summe(a) =  T(b, cOli1)*T(a, cOli2)*summe(b)*summe(a) + i_*cOlf(cOli2, cOli1, cOli11)*T(a, cOli11)*summe(a)*replace_(b, a, pb, pa, sb1, sa1, sb2, sa2, sb3, sa3, sab, 4*ma^2, sib, sia, sjb, sja, mb, ma);
  sum cOli11;
endrepeat;
repeat;
  id T(a, cOli2?)*T(j, cOli1?)*summe(a)*summe(j) = T(j, cOli1)*T(a, cOli2)*summe(a)*summe(j) + i_*cOlf(cOli2, cOli1, cOli10)*T(j, cOli10)*summe(j)*replace_(a, j, pa, pj, sa1, sj1, sa2, sj2, sa3, sj3, sja, 4*mj^2, sia, sij, sab, sjb, ma, mj);
  sum cOli10;
endrepeat;
.sort
id T(a,cOli2?)*T(j,cOli3?)*dsym(cOli1?,cOli2?,cOli3?) = T(j,cOli3)*T(a,cOli2)*dsym(cOli1,cOli2,cOli3);
id cOlf(cOli1?, cOli2?, cOli3?)*cOlf(cOli1?, cOli2?, cOli4?) = CA*delta(cOli3, cOli4);
id delta(cOli1?, cOli2?)*T(i?, cOli2?) = T(i, cOli1);
id delta(cOli1?, cOli2?)*cOlf(cOli1?, cOli3?, cOli4?) = cOlf(cOli2, cOli3, cOli4);
id cOlf(cOli1?, cOli2?, cOli3?)*dsym(cOli4?, cOli2?, cOli3?) = 0;

* rename
if(match(summe(b)) && (match(summe(j))==0));
  mul replace_(b, j, pb, pj, sb1, sj1, sb2, sj2, sb3, sj3, sab, sja, sib, sij, mb, mj);
else if(match(summe(a)) && (match(summe(j))==0));
  mul replace_(a, j, pa, pj, sa1, sj1, sa2, sj2, sa3, sj3, sia, sij, sab, sjb, ma, mj);
else if(match(summe(b)) && (match(summe(a))==0));
  mul replace_(b, a, pb, pa, sb1, sa1, sb2, sa2, sb3, sa3, sib, sia, sjb, sja, mb, ma);
endif;

.sort
#call Kinematics
#call FullSimplify
#call Simplify
.sort
repeat;
  id sa3*den(sa1 + sa2 + sa3) = 1 - sa2*den(sa1 + sa2 + sa3) - sa1*den(sa1 + sa2 + sa3);
endrepeat;
repeat;
  id sa3*den(sa2 + sa3) = 1 - sa2*den(sa2 + sa3);
endrepeat;
if(occurs(sa1, sa2, sa3, sia, sja, ma) == 0) id summe(a) = 0;
.sort
repeat;
  id sj3*den(sj1 + sj2 + sj3) = 1 - sj2*den(sj1 + sj2 + sj3) - sj1*den(sj1 + sj2 + sj3);
endrepeat;
repeat;
  id sj3*den(sj2 + sj3) = 1 - sj2*den(sj2 + sj3);
endrepeat;
if(occurs(sj1, sj2, sj3, sij, sja, mj) == 0) id summe(j) = 0;
.sort
repeat;
  id si3*den(si1 + si2 + si3) = 1 - si2*den(si1 + si2 + si3) - si1*den(si1 + si2 + si3);
endrepeat;
repeat;
  id si3*den(si2 + si3) = 1 - si2*den(si2 + si3);
endrepeat;
if(occurs(si1, si2, si3, sij, sia, mi) == 0) id summe(i) = 0;
*symmetrize result
#if((`i'==2)||(`i'==4))
if((match(summe(a)) == 0) && (match(summe(b)) == 0));
  mul 1/2*(1 + replace_(i, j, j, i, pi, pj, pj, pi, si1, sj1, sj1, si1, si2, sj2, sj2, si2, si3, sj3, sj3, si3, mi, mj, mj, mi));
endif;
.sort
*#call Kinematics
*#call FullSimplify
*#call Simplify
#else if(`i'==3)
if(match(summe(i)) && match(summe(j)) && match(summe(a)) && (match(summe(b)) == 0));
  mul 1/3*(1 + replace_(i, j, j, a, a, i, pi, pj, pj, pa, pa, pi, si1, sj1, sj1, sa1, sa1, si1, si2, sj2, sj2, sa2, sa2, si2, si3, sj3, sj3, sa3, sa3, si3, sij, sja, sia, sij, sja, sia, mi, mj, mj, ma, ma, mi)
             + replace_(i, a, a, j, j, i, pi, pa, pa, pj, pj, pi, si1, sa1, sa1, sj1, sj1, si1, si2, sa2, sa2, sj2, sj2, si2, si3, sa3, sa3, sj3, sj3, si3, sia, sja, sij, sia, sja, sij, mi, ma, ma, mj, mj, mi) );
endif;
*.sort
*#call Kinematics
*#call FullSimplify
*#call Simplify
*.sort
#endif
#enddo
.sort

*#call FullSimplify
*#call Simplify
*mul replace_(si1, piq1*2);
*mul replace_(si2, piq2*2);
*mul replace_(si3, piq3*2);
*mul replace_(sj1, pjq1*2);
*mul replace_(sj2, pjq2*2);
*mul replace_(sj3, pjq3*2);
*mul replace_(s12, q1q2*2);
*mul replace_(s13, q1q3*2);
*mul replace_(s23, q2q3*2);
*mul replace_(sij, pipj*2);
*
*id T(i?, cOli1?) = 1;
*id summe(i?) = 1;
*id TF = 1;
*id CF = 1;
*id CA = 1;
*format C;

format Mathematica;
b den, cOlT, cOlf, T, summe, CA, TF, d, dsym, marker, sij, sia, sij, sja;
print+s ;
.sort

.end