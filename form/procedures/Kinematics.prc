#procedure Kinematics


id p12345 = p1+p2+p3+p4+p5;

argument DS,Eik,EikV;
  id p12345 = p1+p2+p3+p4+p5;
endargument;

#do i = 1,5
  #if ((`i' == 1) || (`i' == 2))
    id p`i'.p`i'=m`i'^2;
  #else
    id p`i'.p`i' = 0;
  #endif
  argument den;
    #if ((`i' == 1) || (`i' == 2))
      id p`i'.p`i'=(m`i')^2;
    #else
      id p`i'.p`i' = 0;
    #endif
  endargument;
  #do j = `i'+1,5
    id (p`i'.p`j')^m? = (s`i'`j'/2)^m;
    id (p`i'.p`j') = (s`i'`j'/2);
    argument den;
      id (p`i'.p`j')^m? = (s`i'`j'/2)^m;
      id (p`i'.p`j') = (s`i'`j'/2);
    endargument;
  #enddo
#enddo

id pi.pj = (sij - mi^2 - mj^2)/2;
id pi.pa = (sia - mi^2 - ma^2)/2;
id pi.pb = (sib - mi^2 - mb^2)/2;
id pj.pa = (sja - mj^2 - ma^2)/2;
id pj.pb = (sjb - mj^2 - mb^2)/2;
id pa.pb = (sab - ma^2 - mb^2)/2;
argument den;
  id pi.pj = (sij - mi^2 - mj^2)/2;
  id pi.pa = (sia - mi^2 - ma^2)/2;
  id pi.pb = (sib - mi^2 - mb^2)/2;
  id pj.pa = (sja - mj^2 - ma^2)/2;
  id pj.pb = (sjb - mj^2 - mb^2)/2;
  id pa.pb = (sab - ma^2 - mb^2)/2;
endargument;

#do i = 1, 3
  id pi.q`i' = (si`i' - mi^2)/2;
  id pj.q`i' = (sj`i' - mj^2)/2;
  id pa.q`i' = (sa`i' - ma^2)/2;
  id pb.q`i' = (sb`i' - mb^2)/2;
  id pi.pi = mi^2;
  id pj.pj = mj^2;
  id pa.pa = ma^2;
  id pb.pb = mb^2;
  id q`i'.q`i' = 0;
  argument den;
    id pi.q`i' = (si`i' - mi^2)/2;
    id pj.q`i' = (sj`i' - mj^2)/2;
    id pa.q`i' = (sa`i' - ma^2)/2;
    id pb.q`i' = (sb`i' - mb^2)/2;
    id pi.pi = mi^2;
    id pj.pj = mj^2;
    id pa.pa = ma^2;
    id pb.pb = mb^2;
    id q`i'.q`i' = 0;
  endargument;
  #do j = `i'+1,3
    id q`i'.q`j' = s`i'`j'/2;
    argument den;
      id q`i'.q`j' = s`i'`j'/2;
    endargument;
  #enddo
#enddo


mul replace_(mi, 0, mj, 0, ma, 0, mb, 0);


.sort

#endprocedure