#procedure Color

#do i=1,3
*  id once T`i'(cOli1?) = T`i'(cOli1)*d_(c`i', cOli1);
#enddo

* Implement casimir operators
id cOlf(c1?, c2?, c3?) * cOlf(c1?, c3?, c4?) = -CA * delta(c2, c4);
id delta(c1?, c2?) * T1(c2?) = T1(c1);
id delta(c1?, c2?) * T2(c2?) = T2(c1);
id delta(c1?, c2?) * cOlf(c2?, cOli1?, cOli2?) = cOlf(c1, cOli1, cOli2);
id cOlf(c1?,cOli1?, cOli2?)*cOlf(c2?,cOli2?,cOli3?)*cOlf(c3?,cOli3?,cOli1?) = CA/2 * cOlf(c1,c2,c3);
id cOlf(c3?,cOli2?, cOli3?)*cOlf(c4?,cOli2?,cOli4?)*cOlf(cOli1?,cOli3?,cOli4?) =  CA/2 * cOlf(cOli1,c3,c4);

*Jacobi Identity
id cOlf(c3,N3_?,N4_?)*cOlf(N1_?,N2_?,N4_?) = -cOlf(c3, N1_?,N4_?)*cOlf(N2_?,N3_?,N4_?) - cOlf(c3,N2_?,N4_?)*cOlf(N3_?, N1_?, N4_?);
sum cOlj8;


#endprocedure