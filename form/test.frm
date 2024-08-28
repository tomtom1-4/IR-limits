#-
#: IncDir procedures
#: SmallExtension 100M
#: WorkSpace      2G
#: MaxTermSize 4M
Off Statistics;

#include declarations.h
i c11, c2, c12;

l test = (dsym(cOli1, cOli2, cOli3)*(cOlT(c2, c12, cOli4)*cOlT(c12, cOli5, cOli1)*cOlT(cOli5, cOli6, cOli2)*cOlT(cOli6, c11, cOli3)*cOlT(c11, c2, cOli4)))/(NF*CF);

repeat;
  id once dsym(cOli1?, cOli2?, cOli3?) = 1/TF*(cOlT(cOli10, cOli11, cOli1)*cOlT(cOli11, cOli12, cOli2)*cOlT(cOli12, cOli10, cOli3)
                                        + cOlT(cOli10, cOli11, cOli2)*cOlT(cOli11, cOli12, cOli1)*cOlT(cOli12, cOli10, cOli3));
  sum cOli10,...,cOli12;
endrepeat;
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

