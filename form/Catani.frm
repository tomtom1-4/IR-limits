* Compare results with Catani et al
* arXiv: 2210.09397

#-
#: IncDir procedures
#: SmallExtension 100M
*#: MaxTermSize    10M
#: WorkSpace      1G
Off Statistics;

#include declarations.h

load results/qqg.sav;
g qqg = i_*sum;
l Catani =
          -1/2*p1.E5/s15*2*TTsym(c5,cOli1)*cOlT(c3,c4,cOli1)*UBar(p3)*g(p1)*V(p4)/s34*den(s13+s14)*2
          +1/2*T(cOli1)*(cOlT(c3, cOli2, c5)*cOlT(cOli2, c4, cOli1) + cOlT(c3, cOli2, cOli1)*cOlT(cOli2, c4, c5))
            *E5(mu)*den(s34+s35+s45)*den(s13+s14+s15)*2*UBar(p3)*(g(p1)*(g(p5)+g(p4))*g(mu)/s45 - g(mu)*(g(p5)+g(p3))*g(p1)/s35)*V(p4)
          +1/2*T(cOli1)*i_*cOlf(c5,cOli1,cOli2)*cOlT(c3,c4,cOli2)*E5(mu)*den(s13+s14+s15)*2*UBar(p3)*(p1(mu)/s34*(2/s15 - 2*den(s13 + s14))*g(p1)
            +den(s34 + s35 + s45)*(1/s34*((s13+s14-s15)*g(mu) - 4*(p3(mu)+p4(mu))*g(p1) + 4*p1(mu)*g(p5)) - g(mu)*(g(p5)+g(p3))*g(p1)/s35 - g(p1)*(g(p5)+g(p4))*g(mu)/s45))*V(p4);
sum cOli1,...,cOli7;
.sort
s marker;
l difference = qqg - Catani;
.sort
*id cOlf(cOli1?, cOli2?, cOli3?) = 0;
id cOlT(c3,cOli1?,c5)*cOlT(cOli1?,c4,cOli2?) = cOlT(c3,cOli1,cOli2)*cOlT(cOli1,c4,c5) + i_*cOlf(c5, cOli2, cOli11)*cOlT(c3, c4, cOli11);
sum cOli11;

#call Dirac
#call WaveFunctions

#call FullSimplify
#call Simplify

bracket cOlf, cOlT, T, TTsym;
format mathematica;
print+s qqg, Catani, difference;

.end