#procedure Cvitanovic
********************************************************************************
* Color structure by Cvitanovic's algorithm

repeat;
  id once cOlf(cOlj1?,cOlj2?,cOlj3?) =
    -i_/TF*(+cOlT(cOlk1,cOlk2,cOlj1)*cOlT(cOlk2,cOlk3,cOlj2)*cOlT(cOlk3,cOlk1,cOlj3)
            -cOlT(cOlk1,cOlk2,cOlj3)*cOlT(cOlk2,cOlk3,cOlj2)*cOlT(cOlk3,cOlk1,cOlj1));
  sum cOlk1, cOlk2, cOlk3; * warning: these have dimension d !!!
endrepeat;

****************************************
* delta's are used not have problems with
* the dimension of the indices

id cOlT(cOli1?,cOli2?,cOlj1?)*cOlT(cOli3?,cOli4?,cOlj1?) =
  TF*(delta(cOli1,cOli4)*delta(cOli2,cOli3)-1/NF*delta(cOli1,cOli2)*delta(cOli3,cOli4));

repeat;
  id delta(cOli1?,cOli2?)*cOlT(cOli2?,cOli3?,?a) = cOlT(cOli1,cOli3,?a);
  id cOlT(cOli3?,cOli2?,?a)*delta(cOli2?,cOli1?) = cOlT(cOli3,cOli1,?a);
  id delta(cOli1?,cOli2?)*delta(cOli2?,cOli3?) = delta(cOli1,cOli3);
  id delta(cOli1?,cOli1?) = NF;
  id delta(cOli2?, c2?)*M(cOli2?, c1?) = M(c2, c1);
endrepeat;

****************************************
* notation of the color package, which
* we, however, do not use here

repeat id cOlT(cOli1?,cOli2?,?a)*cOlT(cOli2?,cOli3?,?b) = cOlT(cOli1,cOli3,?a,?b);
id cOlT(cOli1?,cOli1?,?a) = cOlTr(?a);
id cOlTr(cOli1?) = 0;
id cOlTr(cOli1?, cOli2?) = TF*delta(cOli1, cOli2);

****************************************
* QCD values

*id TF^m? = 1/2^n;
*id NF^m? = 3^n;

.sort

#endprocedure