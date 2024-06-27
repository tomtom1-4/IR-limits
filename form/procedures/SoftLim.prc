#procedure SoftLim
id DS(p1?) = den(p1.p1);
id DS(p1?, m1) = den(p1.p1 - m1^2);
#call Kinematics
print+s d`i';
.sort

mul replace_(p3, lam*p3);
mul replace_(p4, lam*p4);
mul replace_(p5, lam*p5);
id UBar(p3*lam) = lam*UBar(p3);
id V(p4*lam) = V(p4);
*id g(lam*p3?) = lam*g(p3);

#do x={s13, s23, s14, s24, s15, s25}
  mul replace_(`x', lam * `x');
#enddo

#do x={s34, s35, s45}
  mul replace_(`x', lam^2 * `x');
#enddo

.sort
keep brackets;
#do i=1,1
  splitarg ((lam)) den;
  if(match(den(0,?args))) redefine i "0";
  repeat id den(m?,n1?,n2?,?args) = den(m,n1+n2,?args);
  id den(0,m?) = 1/lam * den(m/lam);
  b den;
  .sort
  keep brackets;
#enddo
.sort
mul lam^3;

#$max = 0;
if (count(lam,-1) > $max) $max = count_(lam,-1);

b den;
.sort
keep brackets;
#if `$max' > 0
  id den(n1?,n2?) = den(n1)*sum_(m,0,`$max',(-n2*den(n1))^m);
#else
  id den(n1?,n2?) = den(n1);
#endif
.sort

if(count(lam,1)>0) discard;

factarg den;
repeat id den(n1?,n2?,?args) = den(n1)*den(n2)*den(?args);
id den = 1;

splitarg den;
id den(m?) = 1/m;
repeat id den(n1?,n2?,?args) = den(n1+n2,?args);

#endprocedure
