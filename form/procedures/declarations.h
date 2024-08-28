****************************************
* Interface to DiaGen

s d;
dimension d;

* because of a bug in form one cannot autodeclare!!!

auto i v;
#if 0
#do i=0,19
  i v`i';
  #do j=0,19
    i v`i'l`j';
    i v`i'l`j'in;
    i v`i'l`j'out;
  #enddo
#enddo
#endif

i mu, nu1, nu2, mu1, mu2, mu3, mu4, mu5;
i ro, ro1, ro2, ro3;
i a1, a2, a3, a4;
v k, k1, k2, k3, k4, k5;
v p, n; * n - light-cone gauge vector (= pbar)
v p1, p2, p3, p4, p5, p12345;

ct e0,...,e4;
ct e0t;
ct t0,...,t4;

s s12, s13, s14, s23, s24, s34, s15, s25, s35, s45;
s si1, si2, si3, sj1, sj2, sj3, sa1, sa2, sa3, sb1, sb2, sb3, sij, sia, sib, sja, sjb, sab;
auto s si;
s m1, m2;
s mi, mj, ma, mb;
s p1n, p2n, p3n, p4n;
s m, n1, n2;
t K(symmetric);
s lam;
cf M;

****************************************
* wave functions

cf VWilson;
cf Eps, EpsStar;
v E0, E3, E4, E5; * Eps(...,p123,...), EpsStar(...,p1,...), etc.
v Et3, Et4;
f UfBar, Vf;
f UBar, V, U, VBar;
nt g;

****************************************
* propagator functions

f SF;       * fermion propagator
cf DV;      * vector propator in light-cone gauge with vector n
cf DS, DG, Eik; * DS(k) = 1/k.k, Eik(k,n) = 1/k.n

****************************************
* vertex functions
*
* VVV(k1, a1, k2, a2, k3, a3) = (k1(a3)-k2(a3))*d_(a1, a2)+...
*
* G(i, a) = g_(i, a)

cf VVV, EikV;
f G;

****************************************
* flavor and color

s nl;                    * number of massless quarks running in a closed loop
s CA, CF, TF, NF, NA, dF;    * S(NF) with fundamental representation trace TF
*dimension NA;            * safeguard for loops of adjoint representation deltas
i c0, c1, c2, c3, c4, c5;        * colors of the external partons
auto i cOli, cOlj, cOlk; * other color indices
cf Tc;
t delta(symmetric), cOlf(antisymmetric), cOlT, cOlTr(cyclic), dsym(symmetric);
f T1, T2;
t TTsym(symmetric);
f T;
*CommuteInSet{T1,T2};
****************************************
* utilities
s marker;
cf summe(symmetric);
v pi, pj, pa, pb;
auto i mui;
auto v pi;
v E1, E2, E3;
v q1, q2, q3;
s i,j, a, b;
auto s i;
s ranktag, numtag, Quad;
cf den, den2, rat, tensor;
s b1,...,b4;

****************************************
* C formating

s p1E3, p2E3, p4E3, p1E4, p2E4, p3E4, E3E4;
s denom1,...,denom5;
cf minkovski;

.global
