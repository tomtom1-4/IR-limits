#procedure Scaleless
* Remove scaleles integrals
* For an integral to be non-zero, the integral must depend on all three scales
* p1.p2, p1.p3, p2.p3

splitarg (p3) Eik;
id Eik(p3, p1?) = scales(p3, p1)*Eik(p3, p1);
id Eik(k1?, p3, p1?) = scales(p3, p1)*Eik(k1 + p3, p1);

splitarg (k1) Eik;
id Eik(k1, p1?) = scales1(p1)*Eik(k1, p1);
id Eik(p3?, k1, p1?) = scales1(p1)*Eik(p3 + k1, p1);
id k1.p1? = scales1(p1)*k.p1;

splitarg (k2) Eik;
id Eik(k2, p1?) = scales2(p1)*Eik(k2, p1);
id Eik(p3?, k2, p1?) = scales2(p1)*Eik(p3 + k2, p1);
id k2.p1? = scales2(p1)*k2.p1;

id DS(p3 + k2) = scales2(p3)*DS(p3 + k2);
id DS(p3 + k1) = scales1(p3)*DS(p3 + k1);
id DS(p3 + k1 + k2) = scales1(p3)*scales2(p3)*scales12*DS(p3 + k1 + k2);
id DS(k1 + k2) = DS(k1 + k2)*scales12;
id DS(k1 - k2) = DS(k1 - k2)*scales12;
id k1.k2 = k1.k2*scales12;

repeat;
  id scales1(p2?)*scales1(p1?) = scales1(p2, p1);
  id scales1(p1?, p1?) = scales(p1);
  id scales1(p1?, p2?)*scales1(p3?) = scales1(p1, p2, p3);
  id scales1(p1?, p1?, p2?) = scales(p1, p2);
  id scales2(p2?)*scales2(p1?) = scales2(p2, p1);
  id scales2(p1?, p1?) = scales2(p1);
  id scales2(p1?, p2?)*scales2(p3?) = scales2(p1, p2, p3);
  id scales2(p1?, p2?, p2?) = scales2(p1, p2);
endrepeat;

if(occurs(scales12));
  mul replace_(scales1, scales);
  mul replace_(scales2, scales);
endif;

id scales1(p1?) = 1;
id scales2(p1?) = 1;
id scales1(p1?, p2?) = scales(p1, p2);
id scales2(p1?, p2?) = scales(p1, p2);
id scales1(p1, p2, p3) = scales(p1, p2, p3);
id scales2(p1, p2, p3) = scales(p1, p2, p3);

repeat;
  id once scales(p1?, p2?)^2 = scales(p1, p2);
  id once scales(p1?, p2?, p3?)^2 = scales(p1, p2, p3);
endrepeat;



*if(count(Eik, 1) >= 5)discard;


#endprocedure