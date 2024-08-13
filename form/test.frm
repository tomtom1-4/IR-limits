#-
#: IncDir procedures
#: SmallExtension 100M
#: WorkSpace      2G
#: MaxTermSize 4M
Off Statistics;

#include declarations.h

*l test = Tc(i1,cOli1,cOli2,cOli3)*cOlf(cOli1,cOli3,cOli4);
*
*id Tc(i1, cOli1, cOli2, cOli3) = Tc(i1, cOli2, cOli1, cOli3) + i_*cOlf(cOli1, cOli2, cOli5)*Tc(i1, cOli5, cOli3);
*sum cOli5;
*id Tc(i1?, ?arg1, cOli1?, cOli2?, ?arg2)*cOlf(cOli1?, cOli2?, cOli3?) = 1/2*i_*cOlf(cOli1, cOli2, cOli14)*Tc(i1, ?arg1, cOli14, ?arg2)*cOlf(cOli1, cOli2, cOli3);
*sum cOli14;
*id cOlf(cOli1?, cOli2?, cOli3?)*cOlf(cOli1?, cOli2?, cOli4?) = CA*delta(cOli3, cOli4);
*id delta(cOli4?, cOli4?) = NF;
*id delta(cOli4?, cOli5?)*Tc(i1?, cOli2?, cOli4?) = Tc(i1, cOli2, cOli5);
cf Int;
s s1,s2,s3,s4;
l test = Int(0, 1)*Int(-1,1);
id Int(s1?, s2?)*Int(s3?, s4?) = Int(s1+s3,s2+s4);


*b Tc, cOlf;
print test;
.end

