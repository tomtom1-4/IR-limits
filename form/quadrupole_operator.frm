#-
#: IncDir procedures
#: SmallExtension 100M
#: WorkSpace      2G
*#: MaxTermSize 2M
Off Statistics;

* flavor and color

s CA, CF, TF, NF, NA, dF;    * S(NF) with fundamental representation trace TF
auto i cOli; * other color indices
t delta(symmetric), cOlf(antisymmetric), kronecker(symmetric);
i i, j, k, l, a;
f T;

l Qijkl = 1/2*cOlf(cOli1, cOli2, cOli5)*cOlf(cOli5, cOli3, cOli4)*(T(k, cOli1)*(T(i, cOli3)*T(j, cOli4) + T(j, cOli4)*T(i, cOli3))*T(l, cOli2)
                                                                 + T(l, cOli2)*(T(i, cOli3)*T(j, cOli4) + T(j, cOli4)*T(i, cOli3))*T(k, cOli1));

sum cOli1,...,cOli5;
.sort
#do b=1,4
#do x={k,j,i}
  repeat;
  id T(a?, cOli1?)*T(`x', cOli2?) = T(`x', cOli2)*T(a, cOli1) + i_*cOlf(cOli1, cOli2, cOli3)*T(`x', cOli3)*kronecker(`x', a);
  sum cOli3;
  endrepeat;
  .sort
  id cOlf(cOli1?, cOli2?, cOli3?)*cOlf(cOli1?, cOli2?, cOli4?) = CA*delta(cOli3, cOli4);
  id delta(cOli1?, cOli2?)*T(i?, cOli2?) = T(i, cOli1);
  id delta(cOli1?, cOli2?)*cOlf(cOli1?, cOli3?, cOli4?) = cOlf(cOli2, cOli3, cOli4);
  if(count(kronecker, 1) == 2) id kronecker(i?, k?)*kronecker(j?, k?) = kronecker(i,k)*kronecker(j,k)*kronecker(i,j);
  #do a={j,k,l}
    id kronecker(i, `a')*T(`a', cOli1?) = T(i, cOli1)*kronecker(i, `a');
  #enddo
  #do a={k,l}
    id kronecker(j, `a')*T(`a', cOli1?) = T(j, cOli1)*kronecker(j, `a');
  #enddo
  #do a={l}
    id kronecker(k, `a')*T(`a', cOli1?) = T(k, cOli1)*kronecker(k, `a');
  #enddo
  .sort
#enddo
#enddo
renumber 1;
.sort


*b T;
b cOlf;
print+s;
.end
