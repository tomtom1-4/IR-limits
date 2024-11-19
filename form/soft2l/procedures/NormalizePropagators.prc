#procedure NormalizePropagators

********************************************************************************
* Simple properties and notation for ordinary and eikonal propagators

****************************************
* ordinary massless scalar propagators
id once DS(-p3 + k2) = DS(-p3 + k2)*replace_(k2, -k2);
id once DS(-p3 + k1) = DS(-p3 + k1)*replace_(k1, -k1);
id once Eik(k1 + k2,p1?)*Eik(p3 + k1 + k2,p2?) = Eik(k1 + k2,p1)*Eik(p3 + k1 + k2,p2)*replace_(k1, -k1 - k2);
id once DS(p3 + k2)^2 = DS(p3 + k2)^2*replace_(k2, k2 - p3);

splitarg (k1) DS;
splitarg (k2) DS;
id DS(-k1) = DS(k1);
id DS(-k2) = DS(k2);
id DS(p?,-k1) = DS(-p,k1);
id DS(p?, k2?, -k1) = DS(-p, -k2, k1);
id DS(p?, -k2) = DS(-p, k2);

id k1.k1*DS(k1) = 1;
id k2.k2*DS(k2) = 1;

id DS(p?,k?) = DS(p+k);
id DS(p?, k2?, k1?) = DS(p + k2 + k1);

****************************************
* eikonal propagators

#do i=1,2
  #do k={k1,k2};
    splitarg (`k') Eik;
    #do P={p1,p2}
      id Eik(?args, -`P') = -Eik(?args, `P');
      id Eik(-`k', `P') = -Eik(`k', `P');
      id Eik(p?,-`k', `P') = -Eik(-p,`k', `P');

      id `P'.`k'*Eik(`k', `P') = 1;
      id `P'.`k'*Eik(p?,`k', `P') = 1-`P'.p*Eik(p,`k',`P');

      id Eik(`k', `P')*Eik(p?{p3, k2, k1, -p3, -k1, -k2},`k', `P') = Eik(p, `P')*(Eik(`k',`P')-Eik(p,`k',`P'));
      id Eik(p3?,`k', `P')*Eik(p4?,`k', `P') = Eik(p3-p4,`P')*(Eik(p4,`k',`P')-Eik(p3,`k',`P'));

      id Eik(p?,k1?, `P') = Eik(k1+p,`P');
      id Eik(p3, p1) = 2/s13;
      id Eik(-p3, p1) = -2/s13;
      id Eik(p3, p2) = 2/s23;
      id Eik(-p3, p2) = -2/s23;
    #enddo
  #enddo
#enddo


#endprocedure