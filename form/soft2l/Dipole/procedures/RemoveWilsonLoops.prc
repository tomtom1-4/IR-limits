#procedure RemoveWilsonLoops

* Remove diagrams where the External Wilson lines do not couple to the blob/photon
repeat;
  id once cOlT(cOli1?, cOli2?, ?args)*cOlT(cOli2?, cOli3?, cOli11?) = cOlT(cOli1, cOli3, ?args, cOli11);
endrepeat;
.sort
id VWilson(p1, cOli0?)*cOlT(cOli0?, cOli1?, ?args)*VWilson(p2, cOli1?) = 0;
repeat;
  id once cOlT(cOli1?, cOli2?, cOli3?, cOli4?, ?args) = cOlT(cOli1, cOli10, cOli3)*cOlT(cOli10, cOli2, cOli4, ?args);
  sum cOli10;
endrepeat;

* Remove remaining diagrams with closed Wilson loops
splitarg (p1) EikV;
splitarg (p2) EikV;
id EikV(p3?, p1?, cOli1?) = EikV(p1, cOli1); * remove soft momenta
.sort
argument EikV, 1;
  if((occurs(p1) == 0) && (occurs(p2) == 0));
    discard;
  endif;
endargument;
id EikV(0, cOli1?) = 0;

#endprocedure