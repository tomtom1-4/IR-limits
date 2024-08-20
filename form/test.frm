#-
#: IncDir procedures
#: SmallExtension 100M
#: WorkSpace      2G
#: MaxTermSize 4M
Off Statistics;

#include declarations.h

l test = cOlf(cOli1, cOli2, cOli3)*cOlf(cOli4, cOli5, cOli3)*(cOlT(cOli6, cOli7, cOli1)*(cOlT(cOli7, cOli8, cOli4)*cOlT(cOli9, cOli10, cOli5) + cOlT(cOli9, cOli10, cOli5)*cOlT(cOli7, cOli8, cOli4))*cOlT(cOli10, cOli11, cOli2)
                                                             +cOlT(cOli9, cOli10, cOli2)*(cOlT(cOli6, cOli7, cOli4)*cOlT(cOli10, cOli11, cOli5) + cOlT(cOli10, cOli11, cOli5)*cOlT(cOli6, cOli7, cOli4))*cOlT(cOli7, cOli8, cOli1));
#call Cvitanovic
*
*id Tc(i1, cOli1, cOli2, cOli3) = Tc(i1, cOli2, cOli1, cOli3) + i_*cOlf(cOli1, cOli2, cOli5)*Tc(i1, cOli5, cOli3);
*sum cOli5;
*id Tc(i1?, ?arg1, cOli1?, cOli2?, ?arg2)*cOlf(cOli1?, cOli2?, cOli3?) = 1/2*i_*cOlf(cOli1, cOli2, cOli14)*Tc(i1, ?arg1, cOli14, ?arg2)*cOlf(cOli1, cOli2, cOli3);
*sum cOli14;
*id cOlf(cOli1?, cOli2?, cOli3?)*cOlf(cOli1?, cOli2?, cOli4?) = CA*delta(cOli3, cOli4);
*id delta(cOli4?, cOli4?) = NF;
*id delta(cOli4?, cOli5?)*Tc(i1?, cOli2?, cOli4?) = Tc(i1, cOli2, cOli5);


*b Tc, cOlf;
print test;
.end

