#procedure SoftLim

#call Kinematics

#do x={s13,s23}
  mul replace_(`x', lam * `x');
#enddo

#do x={p3, k1, k2}
  mul replace_(`x', lam * `x');
#enddo

splitarg (lam * k1) DS;
id DS(p?, -k1 * lam) =  DS(-p, k1 * lam);
id DS(0, lam * k1) = DS(k1) * lam^-2;
id DS(0, -lam * k1) = DS(k1) * lam^-2;


#do x = {p1, p2, -p1, -p2}
  #do y = {p3, k2, -p3, -k2, p3 + k2, -p3 - k2, p3 - k2, -p3 + k2}
    id DS(`x',lam * k1) = Eik(k1, `x')/2/lam;
    id DS(`x' + lam * (`y'), lam * k1) = Eik(k1 + `y', `x')/2/lam;
    id DS(`x' + lam * (`y')) = Eik(`y', `x')/2/lam;
    id DS(lam * (`y'), lam * k1) = DS(`y' + k1)/lam^2;
    id DS(lam * (`y')) = DS(`y')/lam^2;
  #enddo
#enddo
id DS(p1 + p2 + lam * p3,lam * k1) = DS(p1 + p2);
id DS(p1 + p2 + p3*lam + k2*lam,k1*lam) = DS(p1 + p2);
id DS( - p1 - p2 - p3*lam - k2*lam,k1*lam) = DS(p1 + p2);
id DS( - p1 - p2 - p3*lam,k1*lam) = DS(p1 + p2);
id DS( - p1 - p2 - k2*lam,k1*lam) = DS(p1 + p2);
id DS( - p1 - p2 - k2*lam) = DS(p1 + p2);
id DS( - p1 - p2 - p3*lam) = DS(p1 + p2);


argument DS;
  if (match(lam));
    exit "Soft limit failed for DS";
  endif;
endargument;

#do x = {k1, k2, p3, p3 + k2, k1 - p3, k1 + p3, k1 - k2, k1 + k2, k1 + p3 + k2, k1 - p3 - k2, k1 - p3 + k2, k1 + p3 - k2}
  id Eik(lam * (`x'), n?) = Eik(`x', n)/lam;
  id Eik(-lam * (`x'), n?) = - Eik(`x', n)/lam;
#enddo

argument Eik;
  if (match(lam));
    print "%t";
    exit "Soft limit failed for Eikonal";
  endif;
endargument;

mul lam^8; * measure
mul lam;

*mul replace_(lam, 0);
.sort

#endprocedure