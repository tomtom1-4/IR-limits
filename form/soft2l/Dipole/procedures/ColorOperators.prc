#procedure ColorOperators

id cOlT(c1, cOli4?, cOli5?) = T1(c1, cOli4, cOli5);
id cOlT(cOli2?, c2, c3?) = T2(cOli2, c2, c3);

repeat;
  id T1(c1?, cOli4?, cOli5?) * cOlT(cOli4?, cOli2?, cOli7?) = T1(c1, cOli4, cOli5) * T1(cOli4, cOli2, cOli7);
  id cOlT(cOli2?, cOli4?, cOli7?) * T2(cOli4?, c2?, cOli5?) =  T2(cOli2, cOli4, cOli7) * T2(cOli4, c2, cOli5);
endrepeat;

#do i = 1,2
  id T`i'(c1?, c2?, c3?) = T`i'(c3);
#enddo

repeat;
  id T2(c1?) * T1(c2?) = T1(c2) * T2(c1);
endrepeat;

if (count(T2, 1) == 1);
  id T2(c3?) = -T2(c3);
elseif(count(T2, 1) == 2);
  id T2(c1?) * T2(c2?) = T2(c2) * T2(c1);
elseif(count(T2, 1) == 3);
  id T2(c1?) * T2(c2?) * T2(c3?) = -T2(c3) * T2(c2) * T2(c1);
elseif(count(T2, 1) == 4);
  id T2(c1?) * T2(c2?) * T2(c3?) * T2(c4?) = T2(c4) * T2(c3) * T2(c2) * T2(c1);
endif;

* The remaining cOlT must be traces
id cOlT(cOli1?, cOli2?, cOli3?) * cOlT(cOli2?, cOli1?, cOli4?) = TF * delta(cOli3, cOli4);
id cOlT(cOli1?, cOli2?, cOli6?) * cOlT(cOli2?, cOli4?, c3) * cOlT(cOli4?, cOli1?, c4) = cOlT(cOli1, cOli2, c3) * cOlT(cOli2, cOli4, c4) * cOlT(cOli4, cOli1, cOli6);
id cOlT(cOli1?, cOli2?, cOli6?) * cOlT(cOli2?, cOli4?, c4) * cOlT(cOli4?, cOli1?, c3) = cOlT(cOli1, cOli2, c4) * cOlT(cOli2, cOli4, c3) * cOlT(cOli4, cOli1, cOli6);
id cOlT(cOli1?, cOli2?, c4) * cOlT(cOli2?, cOli4?, c3) * cOlT(cOli4?, cOli1?, cOli6?)
   = i_ * TF * cOlf(c4, c3, cOli6) + cOlT(cOli1, cOli2, c3) * cOlT(cOli2, cOli4, c4) * cOlT(cOli4, cOli1, cOli6);

if (occurs(T1) && (occurs(T2) == 0));
    mul sumi;
elseif((occurs(T1) == 0) && occurs(T2));
    mul sumj;
elseif (occurs(T1) && occurs(T2));
    mul sumij;
else;
*    exit "Symmetry has failed";
endif;

.sort

#endprocedure

