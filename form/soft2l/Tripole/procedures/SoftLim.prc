#procedure SoftLim

#call Kinematics

***   mul replace_ is used when the variables to replace are also inside the argument of functions    ***

***   Every expression that contains the p3 (soft momenta) gets multiplied by lam   ***
***    s13 and s23 receive a contribution, being p3.pi      ***

#do x={s1q,s2q, s3q}
  mul replace_(`x', lam * `x');
#enddo

*** By the routing the k1 and k2 are alwais the momenta of emitted gluon from the wilson line, soft by definition    ***

#do x={q, k1, k2}
  mul replace_(`x', lam * `x');
#enddo

#do x={q,k1,k2}
  splitarg(lam * `x');
#enddo

*** Split the momenta to obtain DS(hard,soft) ***

repeat;
  id DS(p1?, ?args, k1?, k2?) = DS(p1, ?args, k1 + k2);
  id DS(0, ?args, k1?, k2?) = DS(0, ?args, k1 + k2);
endrepeat;

*** EXtract the lam factor ***
if (match(DS(p1?)));
  Print "%t";
  exit "found non-physical hard momentum";
endif;

.sort;

id DS(p1?,q?)=Eik(q,p1)/2/lam;
id DS(0,p1?)=DS(p1)/lam^2;

argument DS, Eik;
  id lam = 1;
endargument;

argument DS;
  if (match(lam));
    exit "Soft limit failed for DS";
  endif;
endargument;

argument Eik;
  if (match(lam));
    print "%t";
    exit "Soft limit failed for Eikonal";
  endif;
endargument;

mul lam^8; * measure
mul lam;

if (count(lam,1)!=0);
  exit "Sub-leading term in lamba detected";
endif;
.sort


#endprocedure