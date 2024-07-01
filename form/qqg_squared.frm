#-
#: IncDir procedures
#: SmallExtension 100M
#: WorkSpace      2G
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

.sort

* rewrite to c++ format
s piq1, piq2, piq3, pjq1, pjq2, pjq3, q1q2, q1q3, q2q3, pipj, mi2, pkq1, pkq2, pkq3, pipk, pjpk;
*id pi.q1 = piq1;
*id pi.q2 = piq2;
*id pi.q3 = piq3;
*id pj.q1 = pjq1;
*id pj.q2 = pjq2;
*id pj.q3 = pjq3;
*id pa.q1 = pkq1;
*id pa.q2 = pkq2;
*id pa.q3 = pkq3;
*id q1.q2 = q1q2;
*id q1.q3 = q1q3;
*id q2.q3 = q2q3;
*id pi.pj = pipj;
*id pi.pa = pipk;
*id pj.pa = pjpk;
*id pi.pi = mi2;
*argument den;
*  id pi.q1 = piq1;
*  id pi.q2 = piq2;
*  id pi.q3 = piq3;
*  id pj.q1 = pjq1;
*  id pj.q2 = pjq2;
*  id pj.q3 = pjq3;
*  id pa.q1 = pkq1;
*  id pa.q2 = pkq2;
*  id pa.q3 = pkq3;
*  id q1.q2 = q1q2;
*  id q1.q3 = q1q3;
*  id q2.q3 = q2q3;
*  id pi.pj = pipj;
*  id pi.pa = pipk;
*  id pj.pa = pjpk;
*  id pi.pi = mi2;
*endargument;
*id sum(i?) = 1;
*id T(i?, cOli1?) = 1;
*id mi^2 = mi2;
*format C;
*b den;
*print reducible;
*.end

l reducible = TF*sum(i)*sum(j)*sum(a)*sum(b)*(-pi.pb*den(pi.q1*pb.q1)*g_(0,q2, pj, q3, pa)*den(pj.q2 + pj.q3)*den(pa.q2 + pa.q3)/s23^2)*1/4*(T(i,cOli1)*T(j,cOli2)*T(a,cOli2)*T(b,cOli1) + T(j,cOli2)*T(i,cOli1)*T(b,cOli1)*T(a,cOli2) + T(i,cOli1)*T(j,cOli2)*T(b,cOli1)*T(a,cOli2) + T(j,cOli2)*T(i,cOli1)*T(a,cOli2)*T(b,cOli1));

l reducibleCatani = TF*sum(i)*sum(j)*sum(a)*sum(b)*(-pi.pb*den(pi.q1*pb.q1)*(pj.q2*pa.q3 + pj.q3*pa.q2 - pj.pa*q2.q3)*den(pj.q2 + pj.q3)*den(pa.q2 + pa.q3)*den(q2.q3^2))*1/2*(T(i, cOli1)*T(b, cOli1)*T(j, cOli2)*T(a, cOli2) + T(j, cOli2)*T(a, cOli2)*T(i, cOli1)*T(b, cOli1))
  - TF*CA*sum(i)*sum(j)*T(i, cOli1)*T(j, cOli1)*Sij1;

l reducibleDifference = reducible - reducibleCatani;

l interference = -1/2*sum(i)*sum(j)*sum(a)*1/2*TF*dsym(cOli1, cOli2, cOli3)*((T(i, cOli1)*T(j, cOli2) + T(j, cOli2)*T(i, cOli1))*T(a, cOli3) + T(a, cOli3)*(T(i, cOli1)*T(j, cOli2) + T(j, cOli2)*T(i, cOli1)))
  *den(2*(q1.q2 + q1.q3 + q2.q3)*(pa.q1 + pa.q2 + pa.q3))*(g_(0, q2, pa, q13, pi, q3, pj)*den(2*q1.q3) - g_(0, q2, pi, q12, pa, q3, pj)*den(2*q1.q2))
  /s23*den(pj.q2 + pj.q3)*den(pi.q1)*replace_(q13, q1 + q3, q12, q1 + q2, q23, q2 + q3)
  - 1/2*sum(i)*sum(j)*sum(a)*TF*i_*cOlf(cOli1, cOli2, cOli3)*((T(i, cOli1)*T(j, cOli2) + T(j, cOli2)*T(i, cOli1))*T(a, cOli3) + T(a, cOli3)*(T(i, cOli1)*T(j, cOli2) + T(j, cOli2)*T(i, cOli1)))
  *den(pa.q1 + pa.q2 + pa.q3)*den(2*q2.q3)*den(pi.q1)*den(pj.q2 + pj.q3)*(g_(0, q2, pa, q3, pj)*den(2*q2.q3)*pa.pi*(den(pa.q1) - den(pa.q2 + pa.q3))
    + den(2*(q1.q2 + q1.q3 + q2.q3))*(den(2*q2.q3)*(2*(pa.q2 + pa.q3 - pa.q1)*g_(0, q2, pi, q3, pj) - 4*(pi.q2 + pi.q3)*g_(0, q2, pa, q3, pj) + 4*pa.pi*g_(0, q2, q1, q3, pj))
      - den(2*q1.q2)*g_(0, q2, pi, q12, pa, q3, pj) - den(2*q1.q3)*g_(0, q2, pa, q13, pi, q3, pj)))*replace_(q12, q1 + q2, q13, q1 + q3, q23, q2 + q3);

l interferenceCatani = - sum(i)*sum(j)*sum(a)*TF*dsym(cOli1, cOli2, cOli3)*T(i, cOli1)*T(j, cOli2)*T(a, cOli3)*Sija
                       - sum(i)*sum(j)*TF*CA*T(i, cOli1)*T(j, cOli1)*Sij2;

l interferenceDifference = interference + interferenceCatani;

l irreducibleNAB = sum(i)*sum(j)*1/2*i_*cOlf(cOli1, cOli2, cOli3)*T(i, cOli2)*1/2*(-i_)*cOlf(cOli1, cOli4, cOli3)*T(j, cOli4)*TF
  *g_(0, q2)*(den(pi.q1 + pi.q2 + pi.q3)*(pi(mu)*den(2*q2.q3)*(den(pi.q1) - den(pi.q2 + pi.q3))*g_(0, pi)
    + den(2*q1.q2 + 2*q1.q3 + 2*q2.q3)*(den(2*q2.q3)*(2*(pi.q2 + pi.q3 - pi.q1)*g_(0, mu) - 4*(q2(mu) + q3(mu))*g_(0, pi) + 4*pi(mu)*g_(0, q1)) - g_(0, mu, q12, pi)*den(2*q1.q2) - g_(0, pi, q13, mu)*den(2*q1.q3))))
  *g_(0, q3)*(den(pj.q1 + pj.q2 + pj.q3)*(pj(mu)*den(2*q2.q3)*(den(pj.q1) - den(pj.q2 + pj.q3))*g_(0, pj)
    + den(2*q1.q2 + 2*q1.q3 + 2*q2.q3)*(den(2*q2.q3)*(2*(pj.q2 + pj.q3 - pj.q1)*g_(0, mu) - 4*(q2(mu) + q3(mu))*g_(0, pj) + 4*pj(mu)*g_(0, q1)) - g_(0, pj, q12, mu)*den(2*q1.q2) - g_(0, mu, q13, pj)*den(2*q1.q3))))
  *(-1)*replace_(q12, q1 + q2, q13, q1 + q3);

l irreducibleNABCatani = -sum(i)*sum(j)*TF*CA/4*T(i, cOli1)*T(j, cOli1)*Sij3NAB;

l irreducibleNABDifference = irreducibleNAB - irreducibleNABCatani;

l irreducibleAB = replace_(q13, q1 + q3)*replace_(q12, q1 + q2)*sum(i)*sum(j)*colFac*T(i, cOli1)*T(j, cOli1)*(-1)*den(4*(q1.q2 + q1.q3 + q2.q3)^2*(pi.q1 + pi.q2 + pi.q3)*(pj.q1 + pj.q2 + pj.q3))
  *g_(0, q2)*(g_(0, pi, q13, mu)*den(2*q1.q3) - g_(0, mu, q12, pi)*den(2*q1.q2))*g_(0, q3)*(g_(0, mu, q13, pj)*den(2*q1.q3) - g_(0, pj, q12, mu)*den(2*q1.q2));

l irreducibleABCatani = sum(i)*sum(j)*T(i, cOli1)*T(j, cOli1)*(-colFac)*Sij3AB;

l irreducibleABDifference = irreducibleAB - irreducibleABCatani;

tracen 0;
sum cOli1,...,cOli3, mu;
.sort
*drop Sij1;

#do i=1,2
repeat;
  id T(j, cOli2?)*T(i, cOli1?)*sum(i)*sum(j) =  T(i, cOli1)*T(j, cOli2)*sum(i)*sum(j) + i_*cOlf(cOli2, cOli1, cOli10)*T(i, cOli10)*sum(i)*replace_(j, i, pj, pi, sj1, si1, sj2, si2, sj3, si3, sja, sia, sjb, sib, sij, 4*mi^2, mj, mi);
  sum cOli10;
  id T(a, cOli2?)*T(i, cOli1?)*sum(i)*sum(a) =  T(i, cOli1)*T(a, cOli2)*sum(i)*sum(a) + i_*cOlf(cOli2, cOli1, cOli11)*T(i, cOli11)*sum(i)*replace_(a, i, pa, pi, sa1, si1, sa2, si2, sa3, si3, sia, 4*mi^2, sja, sij, sab, sib, ma, mi);
  sum cOli11;
endrepeat;
.sort
repeat;
  id T(j, cOli2?)*T(b, cOli1?)*sum(b)*sum(j) =  T(b, cOli1)*T(j, cOli2)*sum(b)*sum(j) + i_*cOlf(cOli2, cOli1, cOli10)*T(j, cOli10)*sum(j)*replace_(b, j, pb, pj, sb1, sj1, sb2, sj2, sb3, sj3, sjb, 4*mj^2, sab, sja, sib, sij, mb, mj);
  sum cOli10;
  id T(a, cOli2?)*T(b, cOli1?)*sum(b)*sum(a) =  T(b, cOli1)*T(a, cOli2)*sum(b)*sum(a) + i_*cOlf(cOli2, cOli1, cOli11)*T(a, cOli11)*sum(a)*replace_(b, a, pb, pa, sb1, sa1, sb2, sa2, sb3, sa3, sab, 4*ma^2, sib, sia, sjb, sja, mb, ma);
  sum cOli11;
endrepeat;
.sort

id T(j, cOli2?)*T(a, cOli3?)*cOlf(cOli1?, cOli2?, cOli3?)*sum(j)*sum(a) = sum(j)*1/2*i_*cOlf(cOli2, cOli3, cOli4)*T(j, cOli4)*replace_(a, j, pa, pj, sa1, sj1, sa2, sj2, sa3, sj3, sja, 4*mj^2, sia, sij, sab, sjb, ma, mj)*cOlf(cOli1, cOli2, cOli3);
sum cOli4;
id T(a, cOli2?)*T(j, cOli3?)*cOlf(cOli1?, cOli2?, cOli3?)*sum(j)*sum(a) = sum(j)*1/2*i_*cOlf(cOli2, cOli3, cOli4)*T(j, cOli4)*replace_(a, j, pa, pj, sa1, sj1, sa2, sj2, sa3, sj3, sja, 4*mj^2, sia, sij, sab, sjb, ma, mj)*cOlf(cOli1, cOli2, cOli3);
sum cOli4;
id T(b, cOli2?)*T(j, cOli3?)*cOlf(cOli1?, cOli2?, cOli3?)*sum(b)*sum(j) = sum(b)*1/2*i_*cOlf(cOli2, cOli3, cOli4)*T(b, cOli4)*replace_(j, b, pj, pb, sj1, sb1, sj2, sb2, sj3, sb3, sja, sab, sjb, 4*mb^2, mj, mb, sij, sib)*cOlf(cOli1, cOli2, cOli3);
sum cOli4;
id T(b, cOli2?)*T(a, cOli3?)*cOlf(cOli1?, cOli2?, cOli3?)*sum(b)*sum(a) = sum(b)*1/2*i_*cOlf(cOli2, cOli3, cOli4)*T(b, cOli4)*replace_(a, b, pa, pb, sa1, sb1, sa2, sb2, sa3, sb3, sab, 4*mb^2, sia, sib, sja, sjb, ma, mb)*cOlf(cOli1, cOli2, cOli3);
sum cOli4;
id T(a,cOli2?)*T(j,cOli3?)*dsym(cOli1?,cOli2?,cOli3?) = T(j,cOli3)*T(a,cOli2)*dsym(cOli1,cOli2,cOli3);
id cOlf(cOli1?, cOli2?, cOli3?)*cOlf(cOli1?, cOli2?, cOli4?) = CA*delta(cOli3, cOli4);
id delta(cOli1?, cOli2?)*T(i?, cOli2?) = T(i, cOli1);
id delta(cOli1?, cOli2?)*cOlf(cOli1?, cOli3?, cOli4?) = cOlf(cOli2, cOli3, cOli4);
id cOlf(cOli1?, cOli2?, cOli3?)*dsym(cOli4?, cOli2?, cOli3?) = 0;

if(match(sum(b)) && (match(sum(j))==0));
  mul replace_(b, j, pb, pj, sb1, sj1, sb2, sj2, sb3, sj3, sab, sja, sib, sij, mb, mj);
else if(match(sum(a)) && (match(sum(j))==0));
  mul replace_(a, j, pa, pj, sa1, sj1, sa2, sj2, sa3, sj3, sia, sij, sab, sjb, ma, mj);
endif;
.sort
#call Kinematics
#call FullSimplify
#call Simplify
.sort
#if `i'==1
if((match(sum(a)) == 0) && (match(sum(b)) == 0));
  mul 1/2*(1 + replace_(i, j, j, i, pi, pj, pj, pi, si1, sj1, sj1, si1, si2, sj2, sj2, si2, si3, sj3, sj3, si3, mi, mj, mj, mi));
endif;
.sort
#endif
#enddo
.sort
#do i={si1, si2, si3, sij, sia, sib}
  mul replace_(`i', marker*`i');
#enddo

#call FullSimplify
#call Simplify

mul replace_(si1, piq1/2);
mul replace_(si2, piq2/2);
mul replace_(si3, piq3/2);
mul replace_(sj1, pjq1/2);
mul replace_(sj2, pjq2/2);
mul replace_(sj3, pjq3/2);
mul replace_(s12, q1q2/2);
mul replace_(s13, q1q3/2);
mul replace_(s23, q2q3/2);
mul replace_(sij, pipj/2);





format C;
b T, cOlf, TF, sum, i_, pb, den, mi, mj, dsym, marker, d, CA, colFac;
print+s reducible, reducibleCatani, reducibleDifference;
.sort

.end