#-
#: IncDir procedures
#: SmallExtension 100M
#: WorkSpace      2G
Off Statistics;

#include declarations.h
v q13, q12, q23, q123;
s marker;

l Sij1 = (pi.q2/4*den(q2.q3*q2.q3*pi.q1*(pi.q2 + pi.q3))*(pi.pj*den(pj.q1)*(3*pj.q3*den(pj.q2 + pj.q3) - 2*pi.q3*den(pi.q2 + pi.q3)) - 2*pi.pi*pj.q3*den(pi.q1*(pj.q2 + pj.q3)))
  + pi.pj/8*den(q2.q3*pi.q1*(pi.q2 + pi.q3))*(-3*pi.pj*den(pj.q1*(pj.q2 + pj.q3)) + 2*pi.pi*(den(pi.q1*(pj.q2 + pj.q3)) + den(pj.q1*(pi.q2 + pi.q3)))))*(1 + replace_(q2, q3, q3, q2));

l Sija = den(pi.q1*(pj.q2 + pj.q3)*(pa.q1 + pa.q2 + pa.q3)*2*(q1.q2 + q1.q3 + q2.q3))*(den(q1.q2*q2.q3)*(pi.q1*(pj.q2*pa.q3 + pj.q3*pa.q2)
  + pi.q2*(pj.q1*pa.q3 + pj.q3*pa.q1 + 2*pj.q2*pa.q3 + 2*pj.q3*pa.q2)
  + pi.q3*(-pj.q1*pa.q2 + pj.q2*pa.q1) + q1.q3*(pi.pj*pa.q2 - pi.pa*pj.q2 - pj.pa*pi.q2))
  + den(q1.q2)*(pi.pa*pj.q1 - pi.pj*pa.q1 - pj.pa*(pi.q1 + 2*pi.q2))
  + den(q2.q3)*(pj.pa*pi.q3 - pi.pj*pa.q3 - pi.pa*pj.q3))*(1 - replace_(q2, q3, q3, q2));

l Sij2 = den(2*q2.q3*(pi.q1 + pi.q2 + pi.q3))*(2*den(2*(q1.q2 + q1.q3 + q2.q3)*q2.q3)*(pj.q2*(2*mi^2*q1.q3 - (pi.q1 + pi.q2 + pi.q3)*pi.q3)*den(pi.q1*(pj.q2 + pj.q3))
    - pi.q2*(2*pi.pj*q1.q3 - (pi.q1 - pi.q2 + 3*pi.q3)*pj.q3)*den(pj.q1*(pi.q2 + pi.q3)))
  + den(2*(q1.q2 + q1.q3 + q2.q3)*q1.q2)*((2*(pi.q1 + pi.q2)*(pi.pj*q2.q3 - pi.q2*pj.q3 - pi.q3*pj.q2) + mi^2*(pj.q2*q1.q3 - pj.q1*q2.q3))*den(pi.q1*(pj.q2 + pj.q3))
    - (mi^2*((pj.q1 + 2*pj.q2)*q2.q3 + pj.q2*q1.q3) - 2*(pi.q2*pj.q1 + (pi.q1 + 2*pi.q2)*pj.q2)*pi.q3)*den(pj.q1*(pi.q2 + pi.q3)))
  + den(2*(q1.q2 + q1.q3 + q2.q3))*((pi.pj*(pi.q1 + pi.q2 + pi.q3) + mi^2*(pj.q2 - 2*pj.q1))*den(pi.q1*(pj.q2 + pj.q3)) - (3*mi^2*pj.q2 - pi.pj*pi.q1)*den(pj.q1*(pi.q2 + pi.q3)))
  + den(q2.q3)*(pi.q2*(den(pi.q1) - den(pi.q2 + pi.q3))*(mi^2*pj.q3*den(pi.q1*(pj.q2 + pj.q3)) - pi.pj*pi.q3*den(pj.q1*(pi.q2 + pi.q3))))
  + 1/2*mi^2*pi.pj*(den(pi.q2 + pi.q3) - den(pi.q1))*(den(pi.q1*(pj.q2 + pj.q3)) - den(pj.q1*(pi.q2 + pi.q3))))*(1 + replace_(q2, q3, q3, q2));

.sort


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

tracen 0;
sum cOli1,...,cOli3;
.sort
*drop Sij1;

#do i=1,2
repeat;
  id T(j, cOli2?)*T(i, cOli1?)*sum(i)*sum(j) =  T(i, cOli1)*T(j, cOli2)*sum(i)*sum(j) + i_*cOlf(cOli2, cOli1, cOli10)*T(i, cOli10)*sum(i)*replace_(j, i, pj, pi, sj1, si1, sj2, si2, sj3, si3, sja, sia, sjb, sib, sij, 4*mi^2, mj, mi);
  sum cOli10;
  id T(a, cOli2?)*T(i, cOli1?)*sum(i)*sum(a) =  T(i, cOli1)*T(a, cOli2)*sum(i)*sum(a) + i_*cOlf(cOli2, cOli1, cOli10)*T(i, cOli10)*sum(i)*replace_(a, i, pa, pi, sa1, si1, sa2, si2, sa3, si3, sia, 4*mi^2, sja, sij, sab, sib, ma, mi);
  sum cOli10;
endrepeat;
repeat;
  id T(j, cOli2?)*T(b, cOli1?)*sum(b)*sum(j) =  T(b, cOli1)*T(j, cOli2)*sum(b)*sum(j) + i_*cOlf(cOli2, cOli1, cOli10)*T(j, cOli10)*sum(j)*replace_(b, j, pb, pj, sb1, sj1, sb2, sj2, sb3, sj3, sjb, 4*mj^2, sab, sja, sib, sij, mb, mj);
  sum cOli10;
  id T(a, cOli2?)*T(b, cOli1?)*sum(b)*sum(a) =  T(b, cOli1)*T(a, cOli2)*sum(b)*sum(a) + i_*cOlf(cOli2, cOli1, cOli10)*T(a, cOli10)*sum(a)*replace_(b, a, pb, pa, sb1, sa1, sb2, sa2, sb3, sa3, sab, 4*ma^2, sib, sia, sjb, sja, mb, ma);
  sum cOli10;
endrepeat;

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

#call Kinematics
#call FullSimplify
#call Simplify
.sort
#if `i'==1
if((match(sum(a)) == 0) && (match(sum(b)) == 0));
  mul 1/2*(1 + replace_(i, j, j, i, pi, pj, pj, pi, si1, sj1, sj1, si1, si2, sj2, sj2, si2, si3, sj3, sj3, si3, mi, mj, mj, mi));
endif;
#endif
#enddo
.sort

#do x = {si1, sj1, sa1, sb1, s12, s13}
  mul replace_(`x', marker*`x');
#enddo

#call FullSimplify
#call Simplify

argument den;
  id marker = 0;
endargument;

#call FullSimplify
#call Simplify

if(count(marker, 1) > -2) discard;

format mathematica;
b T, cOlf, TF, sum, i_, pb, den, mi, mj, dsym, marker;
print[] ;
.sort

.end