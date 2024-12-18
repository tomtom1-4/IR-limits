#procedure MomentumShift
.sort

#do k={k1,k2}
if ((match(DV(mu1?, mu2?, `k', 0)) == 0) && (match(DV(mu1?, mu2?, -`k', 0)) == 0));
    splitarg (`k') DV;
    id once DV(v1?, v2?, p?, -`k', 0) = DV(v1, v2, p, -`k', 0) * replace_(`k', -`k');
    id once DV(v1?, v2?, p?, `k', 0) = DV(v1, v2, p, `k', 0) * replace_(`k', `k' - p);
    id DV(v1?, v2?, p?, `k'?, 0) = DV(v1, v2, p + `k', 0);
endif;
#enddo


#endprocedure