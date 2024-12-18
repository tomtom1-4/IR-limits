#procedure Color

repeat;
  id T2(c1?) * T1(c2?) = T1(c2) * T2(c1);
  id once disorder T1(c1?) * T1(c2?) = T1(c2) * T1(c1) + i_ * cOlf(c1, c2, cOli104) * T1(cOli104);
  id once disorder T2(c1?) * T2(c2?) = T2(c2) * T2(c1) + i_ * cOlf(c1, c2, cOli105) * T2(cOli105);
  sum cOli104, cOli105;
  renumber 1;
endrepeat;

sum cOli0,...,cOli20;
#do i=1,2
  id T`i'(c1?) * T`i'(c2?) * cOlf(c3?, c1?, c2?) = i_/2 * cOlf(c1,c2,cOli100) * T`i'(cOli100) * cOlf(c3,c1,c2);
  sum cOli100;
#enddo
id M(c1?, c2?) = 1;

* Implement casimir operators
id T1(c1?) * T1(c1?) = C1;
id T2(c1?) * T2(c1?) = C2;
id cOlf(c1?, c2?, c3?) * cOlf(c1?, c3?, c4?) = -CA * delta(c2, c4);
id delta(c1?, c2?) * T1(c2?) = T1(c1);
id delta(c1?, c2?) * T2(c2?) = T2(c1);
id delta(c1?, c2?) * cOlf(c2?, cOli1?, cOli2?) = cOlf(c1, cOli1, cOli2);
id cOlf(c1?,cOli1?, cOli2?)*cOlf(c2?,cOli2?,cOli3?)*cOlf(c3?,cOli3?,cOli1?) = CA/2 * cOlf(c1,c2,c3);
id cOlf(c3?,cOli2?, cOli3?)*cOlf(c4?,cOli2?,cOli4?)*cOlf(cOli1?,cOli3?,cOli4?) =  CA/2 * cOlf(cOli1,c3,c4);

*Jacobi Identity
id cOlf(c3,N3_?,N4_?)*cOlf(N1_?,N2_?,N4_?) = -cOlf(c3, N1_?,N4_?)*cOlf(N2_?,N3_?,N4_?) - cOlf(c3,N2_?,N4_?)*cOlf(N3_?, N1_?, N4_?);
sum cOlj8;

repeat;
  id T2(c1?) * T1(c2?) = T1(c2) * T2(c1);
  id once disorder T1(c1?) * T1(c2?) = T1(c2) * T1(c1) + i_ * cOlf(c1, c2, cOli101) * T1(cOli101);
  id once disorder T2(c1?) * T2(c2?) = T2(c2) * T2(c1) + i_ * cOlf(c1, c2, cOli102) * T2(cOli102);
  sum cOli101, cOli102;
  renumber 1;
endrepeat;

#endprocedure