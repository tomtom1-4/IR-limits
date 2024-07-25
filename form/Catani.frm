* Compare results with Catani et al
* arXiv: 2210.09397

#-
#: IncDir procedures
#: SmallExtension 100M
*#: MaxTermSize    10M
*#: WorkSpace      1G
Off Statistics;

#include declarations.h
s marker;

load results/qqg.sav;
load results/ggg.sav;
g softQQG = i_*sum;
g softGGG = i_*ggg;

l gammaGG = 2*den(s14 + s15)*(p1.E4*p1.E5/s14 + 2/s45*(p1.E4*p4.E5 + 1/2*E4.E5*s15/2))*(1 - replace_(p4, p5, p5, p4, s14, s15, s15, s14, E4, E5, E5, E4));

l gammaGGG = 2*den(s13 + s14 + s15)*(1/12*p1.E3*p1.E4*p1.E5*(3*s15 - (s13 + s14))*4/s14/s15*den(s13 + s14)
    + p1.E5*(s15 - s13 - s14)/s15*den(s13 + s14)/s34*2*(1/2*E3.E4*s13/2 + p1.E4*p4.E3)
    + den(s34 + s35 + s45)/s34*(s34*p1.E3*E4.E5 + 2*p4.E3*E4.E5*(s15 - s13 - s14)/2 + 4*p5.E3*p3.E4*p1.E5
    + 4*p4.E3*p1.E4*(p3.E5 + p4.E5) + E3.E4*(s45*p1.E5 + p3.E5*(s13 + s15 - 3*s14)/2)))
  *(1 - replace_(p3, p4, p4, p3, s13, s14, s14, s13, s35, s45, s45, s35, E3, E4, E4, E3));

.sort

l CataniQQG = -1/2*p1.E5/s15*2*TTsym(c5,cOli1)*cOlT(c3,c4,cOli1)*UBar(p3)*g(p1)*V(p4)/s34*den(s13+s14)*2
          +1/2*T(cOli1)*(cOlT(c3, cOli2, c5)*cOlT(cOli2, c4, cOli1) + cOlT(c3, cOli2, cOli1)*cOlT(cOli2, c4, c5))
            *E5(mu)*den(s34+s35+s45)*den(s13+s14+s15)*2*UBar(p3)*(g(p1)*(g(p5)+g(p4))*g(mu)/s45 - g(mu)*(g(p5)+g(p3))*g(p1)/s35)*V(p4)
          +1/2*T(cOli1)*i_*cOlf(c5,cOli1,cOli2)*cOlT(c3,c4,cOli2)*E5(mu)*den(s13+s14+s15)*2*UBar(p3)*(p1(mu)/s34*(2/s15 - 2*den(s13 + s14))*g(p1)
            +den(s34 + s35 + s45)*(1/s34*((s13+s14-s15)*g(mu) - 4*(p3(mu)+p4(mu))*g(p1) + 4*p1(mu)*g(p5)) - g(mu)*(g(p5)+g(p3))*g(p1)/s35 - g(p1)*(g(p5)+g(p4))*g(mu)/s45))*V(p4);

l CataniGGG = p1.E3*p1.E4*p1.E5*8/s13/s14/s15*1/6*(T(c3)*T(c4)*T(c5) + T(c3)*T(c5)*T(c4) + T(c4)*T(c3)*T(c5) + T(c4)*T(c5)*T(c3) + T(c5)*T(c3)*T(c4) + T(c5)*T(c4)*T(c3))
    + p1.E3*2/s13*gammaGG*i_*cOlf(c4, c5, cOli1)*1/2*(T(c3)*T(cOli1) + T(cOli1)*T(c3))
      *(1 + replace_(c3, c4, c4, c3, s13, s14, s14, s13, p3, p4, p4, p3, s45, s35, s35, s45, E3, E4, E4, E3)
          + replace_(c3, c5, c5, c3, s13, s15, s15, s13, p3, p5, p5, p3, s45, s34, s34, s45, E3, E5, E5, E3))
    + cOlf(c3, c4, cOli1)*cOlf(cOli1, c5, cOli2)*T(cOli2)*gammaGGG
      *(1 + replace_(E5, E3, E3, E5, p5, p3, p3, p5, s15, s13, s13, s15, s45, s34, s34, s45, c5, c3, c3, c5)
          + replace_(E5, E4, E4, E5, p5, p4, p4, p5, s15, s14, s14, s15, s35, s34, s34, s35, c5, c4, c4, c5));


sum cOli1,...,cOli7;
.sort
l differenceQQG = softQQG - CataniQQG;
l differenceGGG = softGGG - CataniGGG;

.sort
#call FullSimplify
#call Simplify
*b cOlf, T, marker;
.sort

id cOlT(c3,cOli1?,c5)*cOlT(cOli1?,c4,cOli2?) = cOlT(c3,cOli1,cOli2)*cOlT(cOli1,c4,c5) + i_*cOlf(c5, cOli2, cOli11)*cOlT(c3, c4, cOli11);
sum cOli11;
#call Dirac
#call WaveFunctions
.sort
repeat;
  id T(cOli1?)*T(c3) = T(c3)*T(cOli1) + i_*cOlf(cOli1, c3, cOli15)*T(cOli15);
  sum cOli15;
endrepeat;
repeat;
  id T(cOli1?!{c3})*T(c4) = T(c4)*T(cOli1) + i_*cOlf(cOli1, c4, cOli14)*T(cOli14);
  sum cOli14;
endrepeat;
repeat;
  id T(cOli1?!{c3,c4})*T(c5) = T(c5)*T(cOli1) + i_*cOlf(cOli1, c5, cOli13)*T(cOli13);
  sum cOli13;
endrepeat;
.sort
* Jacobi
id cOlf(c3, cOli1?, cOli2?)*cOlf(c4, c5, cOli2?) = -cOlf(c4, cOli1, cOli2)*cOlf(c5, c3, cOli2) - cOlf(c5, cOli1, cOli2)*cOlf(c3, c4, cOli2);

#call FullSimplify
#call Simplify

*if(count(T, 1) <= 1) discard;
*bracket cOlf, cOlT, T, TTsym;
ab cOlf;

format mathematica;
*print+s softQQG, CataniQQG, differenceQQG;
print+s CataniGGG, softGGG, differenceGGG;
.end