#procedure Simplify
********************************************************************************
* Simplification of rational functions

ab DS, Eik, E3, E4, E5, nl, cOlf, cOlT, cOlTr, T1, T2, TTsym, CA, CF, dsym, TF,  d, UBar, V, i_, sij, sia, sib, sja, sjb, sab, summe, T, cOlf;
*b s12, s13, s14, s15, s23, s24, s25, s34, s35, s45, m1, m2, si1, si2, si3, sj1, sj2, sj3, sa1, sa2, sa3, sb1, sb2, sb3, mi, mj, ma, mb, marker, den, rat, T, summe, TF, CA, colFac, cOlf, dsym;
*print+s;
.sort
PolyRatFun rat;
*Keep Brackets;
factarg den;
#do x={s12, s13, s14, s15, s23, s24, s25, s34, s35, s45, m1, m2, si1, si2, si3, sj1, sj2, sj3, sa1, sa2, sa3, sb1, sb2, sb3, mi, mj, ma, mb, marker, si11, si12, si13, si21, si22, si23, si31, si32, si33, si41, si42, si43, si51, si52, si53, si61, si62, si63}
  id `x'^m? = rat(`x'^m, 1);
#enddo

#do i=1,1
  id once den(m?, ?args1) * rat(s12?, s13?) = rat(s12, m * s13) * den(?args1);
  id once rat(s12?, s13?) = rat(s12, 1) * den2(s13);
  id once den = 1;
  if (count(den, 1) != 0);
    redefine i "0";
  endif;
  .sort
#enddo

.sort
PolyRatFun;
.sort
id den2(m?) = den(m);
factarg den;
repeat;
  id once den(m?, ?args1) = den2(m) * den(?args1);
endrepeat;
id den2(m?) = den(m);
id den = 1;
.sort
splitarg den;
id den(m?) = 1/m;
repeat;
  id once den(s12?, s13?, ?args1) = den(s12 + s13, ?args1);
endrepeat;
.sort
id rat(s12?, 1) = s12;
.sort

#endprocedure