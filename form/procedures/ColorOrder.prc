#procedure ColorOrder

if(count(summe, 1) >= 1);
  if(occurs(summeP) || (count(summe, 1) > 1));
    print "%t";
    exit "summe occoured in ColorOrder";
  endif;
endif;

#do i=1,10
  #do j=`i'+1,11
    repeat;
      id once summeP(i`j')*summeP(i`i')*T(i`j', c1?)*T(i`i', c2?) = summeP(i`j')*summeP(i`i')*T(i`i', c2)*T(i`j', c1);
*      id once cOld(i`j', c1?, c2?)*T(i`i', c3?) = T(i`i', c3)*cOld(i`j', c1, c2);
    endrepeat;
  #enddo
#enddo
repeat;
  id once disorder T(i1?, c1?)*T(i1?, c2?) = T(i1, c2)*T(i1, c1) + i_*cOlf(c1, c2, c10)*T(i1, c10);
  sum c10;
endrepeat;
*id disorder cOld(i1?, c2?, c1?) = cOld(i1, c1, c2);

repeat;
  id once T(i1?, c1?)*T(i1?, c1?) = C(i1);
  id once T(i1?, c1?)*T(i1?, c2?)*cOlf(c1?, c2?, c3?) =1/2*i_*Cg*T(i1, c3);
  id once cOlf(c1?, c2?, c3?)*cOlf(c1?, c2?, c4?) = Cg*d_(c3, c4);
*id cOlf(c1?, c2?, c3?)*cOld(i1?, c2?, c3?) = 0;
  id once T(i1?,c1?)*T(i1?,c2?)*T(i1?,c3?)*cOlf(b?,c1?,c3?) = T(i1, c2)*i_/2*Cg*T(i1, b) + i_*cOlf(c1, c2, c10)*T(i1, c10)*T(i1, c3)*cOlf(b, c1, c3);
  sum c10;
  id once T(i2,c2?)*T(i2,c3?)*cOlf(b,c1?,c2?) = (T(i2,c3)*T(i2,c2) + i_*cOlf(c2, c3, c11)*T(i2, c11))*cOlf(b, c1, c2);
  sum c11;
endrepeat;
#endprocedure