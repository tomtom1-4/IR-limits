#-
#: IncDir procedures
#: SmallExtension 200M
#: WorkSpace      4G
#: MaxTermSize 4M
Off Statistics;

format 160;

f J0,J1,J2; * Current expansion coefficients
f Z0,Z1,Z2; * Z operator expansion coefficients
f M0,M1,M2; * Matrix Element Expansion
f F0,F1,F2; * Finite Part Expansion
f MBorn0, MBorn1, MBorn2;
f FBorn0, FBorn1, FBorn2;
s a;
f Gamma0, Gamma1, GammaP0, GammaP1, current;
f anhililate;
s ep;
f T;
s gCusp0, gCusp1;
auto i c;
i b,d;
auto i i,j,k,l,m;
s n;
t cOlf(antisymmetric);
f cOld;
cf log(symmetric), gamma0, gamma1, C, j; * j = -p/p.q
s beta0, beta1, ZUV0, ZUV1;
cf jj(antisymmetric), jjSym(symmetric); *jj(i, j) = p_i/p_i.q - p_j/p_j.q;
cf summeP(symmetric), summe(symmetric);
cf F;
i q;
s g;
s Cg, TF, nl, CF;
cf zeta;
s O; * marker of higher orders in epsilon
auto s marker, Marker;
s sumg0, sumg1;
s sumgCusp0, sumgCusp1;
cf replace(symmetric); *replace helper function
s tripole;

l Finite1 = M1 - Z1*M0;
l Finite2 = M2 - Z1*M1 - (Z2 - Z1*Z1)*M0;
*l test = (1 - a*Z1 + a^2*(Z1^2 - Z2))*(J0(b) + a*J1(b) + a^2*J2(b))*(1 + a*Z1 + a^2*Z2);

id M2 = J0(b)*MBorn2 + J1(b)*MBorn1 + J2(b)*MBorn0; * + (ZUV0*(3/2)*M1 + ZUV1/2*M0);
id M1 = J0(b)*MBorn1 + J1(b)*MBorn0; * + ZUV0/2*M0;
id M0 = J0(b)*MBorn0;

id MBorn2 = FBorn2 + Z1*FBorn1 + Z2*FBorn0;
id MBorn1 = FBorn1 + Z1*FBorn0;
id MBorn0 = FBorn0;

.sort
id once J1(b?) = J1(b) + ZUV0/2*J0(b);
id once J2(b?) = J2(b) + (ZUV1/2 - ZUV0^2/8)*J0(b) + 3/2*ZUV0*J1(b);

if(expression(Finite1) || expression(Finite2));
  if(occurs(FBorn0)==0) discard;
endif;
b FBorn0, FBorn1, FBorn2, a;
print ;
.sort

* last term in Z2 is the interference between UV and IR poles. n - 2 -> n - 0 comes from the UV renormalization of Z1.
id Z2 = GammaP0*GammaP0/32/ep^4 + GammaP0/8/ep^3*(Gamma0 - 3/2*beta0) + Gamma0/8/ep^2*(Gamma0 - 2*beta0) + GammaP1/16/ep^2 + Gamma1/4/ep;
id Z1 = (GammaP0/4/ep^2 + Gamma0/2/ep);

id ZUV0 = -beta0/ep;
id ZUV1 = (beta0^2/ep^2 - beta1/ep/2);
*if(count(GammaP0, 1)==0) discard;

.sort
*l test = + ep^-4 * (
*          + 1/32*J0(b)*GammaP0*GammaP0*FBorn0
*          - 1/16*GammaP0*J0(b)*GammaP0*FBorn0
*          + 1/32*GammaP0*GammaP0*J0(b)*FBorn0
*          )
*
*       + ep^-2 * (
*          + 1/4*J1(b)*GammaP0*FBorn0
*          - 1/4*GammaP0*J1(b)*FBorn0
*          )
*
*       + J2(b)*FBorn0;
.sort
id FBorn0 = current*FBorn0;

#do l=0,1
  #do k=0,1
    id once J0(b?)*current = current*J0(b);
    id once J1(b?)*current = current*J1(b);
    id once J2(b?)*current = current*J2(b);
    id once anhililate*current = current*anhililate;
    #do i={1,3}
      #define j "2"
      #if( `i'== 3 )
        #redefine j "4"
      #endif
      if((occurs(i`i')==0) && (occurs(i`j')==0));
        id once Gamma`l'*current = current*(summeP(i`i', i`j')*T(i`i', c1)*T(i`j', c1)*(-gCusp`l'/2*log(i`i',i`j')) + sumg`l'
          + (summe(i`i')*2*T(i`i', c1)*i_*cOlf(c1, b, d)*(-gCusp`l'/2*log(i`i',q))*replace_(b,d,d,b) + gamma`l'(g))*anhililate);
        sum c1, d;
      endif;
      .sort
    #enddo

    #do i={5,7}
      #define j "6"
      #if( `i' == 7 )
        #redefine j "8"
      #endif

* With color conservation
      if((occurs(i`i')==0) && (occurs(i`j')==0));
        id once GammaP`l'*current = current*(-sumgCusp`l' - Cg*gCusp`l'*anhililate);
* Without color conservation
* id has to be applied from right to left. current serves as a marker for the multiplication
*        id once GammaP`l'*current = current*summeP(i`i', i`j')*T(i`i', c1)*T(i`j', c1)*gCusp`l'
*          + 2*current*summe(i`i')*T(i`i', c1)*i_*cOlf(c1, b, d)*gCusp`l'*anhililate*replace_(b,d,d,b);
*        sum c1, d;
      endif;
      .sort
    #enddo
  #enddo
#enddo
id current = 1;

if(occurs(Gamma0) || occurs(Gamma1));
  print "%t";
  exit "Gamma has not been replaced";
endif;

if(occurs(GammaP0) || occurs(GammaP1));
  print "%t";
  exit "GammaP has not been replaced";
endif;

repeat;
  id once anhililate*T(i1?, b?) = T(i1, b)*anhililate;
endrepeat;

b anhililate, J0, J1;
print+s;
.sort

.sort
#do l=0,2
  id anhililate*FBorn`l' = 0;
  repeat;
    id once anhililate*J`l'(b?) = J`l'(b);
  endrepeat;
#enddo
if(occurs(anhililate));
  print "%t" ;
  exit "Error: anhililate still occurs";
endif;
.sort


*id J0(b?) = summeP(i9, i10)*i_*cOlf(b, c1, c2)*T(i9, c1)*T(i10, c2)*jj(i9, i10)*(-1/Cg);
id J0(b?) = summe(i9)*T(i9, b)*j(i9);
id J1(b?) = summeP(i9, i10)*i_*cOlf(b, c1, c2)*T(i9, c1)*T(i10, c2)*jj(i9, i10)*(1/ep^2 + pi_^2/12 - 7/3*zeta(3)*ep - 13*pi_^4/480*ep^2 + ep^3*O)*(1 + ep*log + ep^2/2*log^2 + ep^3/6*log^3 + ep^4/24*log^4 + ep^5*O);*/12*(pi_^2*log + 2*log^3 -28*zeta(3)) + 0*ep^2/480*(-13*pi_^4 + 20*log*(pi_^2*log + log^3 - 56*zeta(3))) + ep^3*O);
id J2(b?) = summeP(i9, i10)*i_*cOlf(b, c1, c2)*T(i9, c1)*T(i10, c2)*jj(i9, i10)*(1 + ep*log + ep^2/2*log^2 + ep^3/6*log^3 + ep^4/24*log^4 + ep^5*O)^2
  *(-1)*(Cg*(1/2/ep^4 - 11/12/ep^3 + 1/ep^2*(pi_^2/6 - 67/36) - 1/ep*(11/6*zeta(3) + 11/12*pi_^2/6 + 193/54))
        +TF*nl*(1/3/ep^3 + 5/9/ep^2 + 1/ep*(pi_^2/18 + 19/27)))
  + (-1)*summeP(i9, i10, i11)*cOlf(b, c3, c4)*cOlf(c4, c1, c2)*T(i9, c1)*T(i10, c2)*T(i11, c3)*(jj(i9, i11)*F(i9, i10, i11) - jj(i10, i11)*F(i10, i9, i11))*(1 + ep*log + ep^2/2*log^2 + ep^3/6*log^3 + ep^4/24*log^4 + ep^5*O)^2*(1 + ep^2*pi_^2/6);

id F(i1?, i2?, i3?) = -1/ep^2*(log(i2, i3) + log(i1, q) - log(i1, i2) - log(i3, q))*(log(i1, i3) + log(i2, q) - log(i1, i2) - log(i3, q))
                    + 1/3/ep*(-(log(i2, i3) + log(i1, q) - log(i1, i2) - log(i3, q))^2*(log(i1, i3) + log(i2, q) - log(i1, i2) - log(i3, q))
                             -2*(log(i2, i3) + log(i1, q) - log(i1, i2) - log(i3, q))*(log(i1, i3) + log(i2, q) - log(i1, i2) - log(i3, q))^2);
id log = (log(i9, i10) - log(i9, q) - log(i10, q));
id jj(i1?, i2?) = -j(i1) + j(i2);
*id log(i1?, i2?) = 0;
sum c1, c2, c3, c4;
.sort
if(count(ep, 1) > 0) discard;
b summe, summeP;
print+s;
.sort
repeat;
  id once summe(i1?)*summeP(i2?, i3?) = summeP(i1, i2, i3) + summeP(i2, i3)*(replace_(i1, i2) + replace_(i1, i3));
  id once summe(i1?)*summeP(i2?, i3?, i4?) = summeP(i1, i2, i3, i4) + summeP(i2, i3, i4)*(replace_(i1, i2) + replace_(i1, i3) + replace_(i1, i4));
  id once summe(i1?)*summeP(i2?, i3?, i4?, i5?) = summeP(i1, i2, i3, i4, i5) + summeP(i2, i3, i4, i5)*(replace_(i1, i2) + replace_(i1, i3) + replace_(i1, i4) + replace_(i1, i5));
  id once summeP(i1?, i2?)*summeP(i3?, i4?) = summeP(i1, i2, i3, i4) + summeP(i1, i2, i3)*(replace_(i4, i1) + replace_(i4, i2))
                                                                + summeP(i1, i2, i4)*(replace_(i3, i1) + replace_(i3, i2))
                                                                + summeP(i1, i2)*(replace_(i3, i1, i4, i2) + replace_(i3, i2, i4, i1));
  id once summeP(i1?, i2?)*summeP(i3?, i4?, i5?) = summeP(i1, i2, i3, i4, i5)
                                            + summeP(i1, i2, i3, i4)*(replace_(i5, i1) + replace_(i5, i2))
                                            + summeP(i1, i2, i3, i5)*(replace_(i4, i1) + replace_(i4, i2))
                                            + summeP(i1, i2, i4, i5)*(replace_(i3, i1) + replace_(i3, i2))
                                            + summeP(i1, i2, i3)*(replace_(i4, i1, i5, i2) + replace_(i4, i2, i5, i1))
                                            + summeP(i1, i2, i4)*(replace_(i3, i1, i5, i2) + replace_(i3, i2, i5, i1))
                                            + summeP(i1, i2, i5)*(replace_(i3, i1, i4, i2) + replace_(i3, i2, i4, i1));
  id once summe(i1?)*summe(i2?) = summeP(i1, i2) + summe(i1)*replace_(i2, i1);
endrepeat;

if(match(log(i1?, i1?)));
  print "%t";
  exit "Wrong log detected";
endif;

if(count(summe, 1) >= 2);
  print "%t";
  exit "More than one summe detected";
endif;

if(count(summeP, 1) >= 2);
  print "%t";
  exit "More than one summeP detected";
endif;

b summe, summeP;
print+s;

.sort

* rewrite summeP into product form
repeat;
  id once summeP(i10?, i11?, i12?, i13?, i14?) = summeP(i10, i11, i12, i13)*summeP(i14);
  id once summeP(i10?, i11?, i12?, i13?) = summeP(i10, i11, i12)*summeP(i13);
  id once summeP(i10?, i11?, i12?) = summeP(i10, i11)*summeP(i12);
  id once summeP(i10?, i11?) = summeP(i10)*summeP(i11);
endrepeat;
.sort

#do i=1,11
  if(match(summe(i`i')) || match(summeP(i`i')));
    if(match(T(i`i', b?))==0);
      print "%t";
      exit "Error sum without color operator detected";
    endif;
  endif;
#enddo

* rename
#do l=1,3
* find open index slots
#do i=1,10
  #do j=`i'+1,11
    if((match(summe(i`i'))==0) && (match(summeP(i`i'))==0));
      if(match(summe(i`j'))||match(summeP(i`j')));
        mul replace_(i`j', i`i');
        redefine j "12";
        redefine i "1";
      endif;
    endif;
    .sort
  #enddo
#enddo
* switch indices
#do i=1,5
  #do j=`i'+1,6
    if(match(summe(i`i')*summeP(i`j')));
      mul replace_(i`j', i`i', i`i', i`j');
      redefine j "10";
    endif;
    .sort
  #enddo
#enddo
#enddo
.sort


id gCusp0 = 4;
id gCusp1 = (268/9 - 4*pi_^2/3)*Cg - 80/9*TF*nl;
id gamma0(g) = -beta0;
id gamma1(g) = Cg^2*(-692/27 + 11/18*pi_^2 + 2*zeta(3)) + Cg*TF*nl*(256/27 - 2*pi_^2/9) + 4*CF*TF*nl;
id beta0 = 11/3*Cg - 4/3*TF*nl;
id beta1 = 34/3*Cg^2 - 4*CF*TF*nl - 20/3*Cg*TF*nl;
.sort

if(match(log(i1?, i1?)));
  print "%t";
  exit "Wrong log detected";
endif;

.sort


* apply color algebra
#do l=1,4
* rename
* find open index slots
  #do i=1,9
    #do j=`i'+1,10
      if((match(summe(i`i'))==0) && (match(summeP(i`i'))==0));
        if(match(summe(i`j'))||match(summeP(i`j')));
          mul replace_(i`j', i`i');
          redefine j "10";
          redefine i "1";
        endif;
      endif;
      .sort
    #enddo
  #enddo
* switch indices
  #do i=1,5
    #do j=`i'+1,6
      if(match(summe(i`i')*summeP(i`j')));
        mul replace_(i`j', i`i', i`i', i`j');
        redefine j "10";
      endif;
      .sort
    #enddo
  #enddo

  .sort
  #call ColorOrder

* reverse color conservation
  if((occurs(i2)==0) && (occurs(summeP)==0) && (count(summe, 1)==1) && (count(T, 1)==1));
   id once summe(i1)*T(i1, b)*FBorn0 = summeP(i1)*summeP(i2)*2*i_/Cg*cOlf(b, c21, c22)*T(i1, c21)*T(i2, c22)*FBorn0;
   sum c21, c22;
  endif;
*  if(match(jj(i1, i2))==0);
*    id jj(i1, i3?) = jj(i1, i3)*replace_(i3, i2, i2, i3);
*    id jj(i2, i3?) = jj(i2, i3)*replace_(i3, i1, i1, i3);
*  endif;
  .sort
* rename indices
  id once j(i2?!{i1}) = j(i2)*replace_(i2, i1, i1, i2);
  .sort
  #call ColorOrder

  if(count(T, 1)==4);
    id T(i1,c1?)*T(i2,c2?)*T(i3,c3?)*T(i4,c2?) = T(i1,c1)*T(i2,c2)*T(i3,c3)*T(i4,c2)*replace_(i2, i3, i3, i2);
    id T(i1,c1?)*T(i2,c2?)*T(i3,c2?)*T(i4,c3?) = T(i1,c1)*T(i2,c2)*T(i3,c2)*T(i4,c3)*replace_(i2, i4, i4, i2);
*    id T(i1,c1?)*T(i2,c2?)*T(i2,c3?)*T(i2,c4?) = T(i1,c1)*T(i2,c2)*T(i2,c3)*T(i2,c4)*replace_(i1, i2, i2, i1);
  endif;
*  if(count(T, 1)==3) id T(i1, c1?)*T(i2, c2?)*T(i2, c3?) = T(i1, c1)*T(i2, c2)*T(i2, c3)*replace_(i1, i2, i2, i1);
  id j(i2?!{i1}) = j(i2)*replace_(i2, i1, i1, i2);
  #call ColorOrder
  renumber 1;

* jacobi identity
*  id cOlf(b, N2_?,N4_?)*cOlf(N3_?, N1_?, N4_?) = -cOlf(b,N1_?,N4_?)*cOlf(N2_?,N3_?, N4_?) - cOlf(b, N3_?, N4_?)*cOlf(N1_?, N2_?, N4_?);
  id cOlf(b,N1_?,N4_?)*cOlf(N2_?,N3_?, N4_?) = - cOlf(b, N2_?,N4_?)*cOlf(N3_?, N1_?, N4_?) - cOlf(b, N3_?, N4_?)*cOlf(N1_?, N2_?, N4_?);
*  id cOlf(b,N1_?,N3_?)*cOlf(c,N2_?,N3_?) = -cOlf(N1_?, N2_?, N3_?)*cOlf(c, b, N3_?) - cOlf(N2_?, b, N3_?)*cOlf(c, N1_?, N3_?);
  #call ColorOrder
*  id log(i1, i2) = log + log(i1, q) + log(i2, q);
  .sort

* color conservation
  #do j=1,4
*  mul Marker;
* reverse color conservation
  if((occurs(i2)==0) && (occurs(summeP)==0) && (count(summe, 1)==1) && (count(T, 1)==1));
    id once summe(i1)*T(i1, b)*FBorn0 = summeP(i1)*summeP(i2)*2*i_/Cg*cOlf(b, c21, c22)*T(i1, c21)*T(i2, c22)*FBorn0;
    sum c21, c22;
  endif;
  .sort
*  if((count(summe, 1)!=0));
*    print "%t";
*    exit "Illegal summation structure detected";
*  endif;

  if(match(T(i`j', c1?)));
    id T(i`j', c1?) = T(i`j', c1)*marker`j';
    if((match(j(i`j'))==0) && (match(jj(i`j', i12?))==0) && (match(jjSym(i`j', i12?))==0) && (match(log(i`j', i10?))==0) && (match(log)==0)
        && (match(C(i`j'))==0) && (match(gamma0(i`j'))==0) && (match(gamma1(i`j'))==0));
      if(count(marker`j', 1)==1);
        repeat;
          id once summeP(i`j')*summeP(i10?!{i`j'})*T(i`j', c1?)*T(i10?!{i`j'}, c2?) =  summeP(i`j')*summeP(i10)*T(i10, c2)*T(i`j', c1);
        endrepeat;
        if(count(summeP, 1) > 5);
          print "%t";
          exit "To many summations detected.";
        endif;
        #do k=1,5
          #if `j' != `k'
            if(match(summeP(i`k')));
              #message `j' `k'
              if((count(summeP, 1) > 2) && (count(summe, 1) == 0));
                if(occurs(marker99)==0) id once summeP(i`j')*T(i`j', c1?)*FBorn0 = T(i`j', c1)*FBorn0*(summeP(i`j')*marker100 - marker99*replace_(i`j', i`k', marker100, 1));
              else if((count(summeP, 1) == 2) && (count(summe, 1) == 0));
*                if(occurs(marker99)==0) id once summeP(i`j')*T(i`j', c1?)*FBorn0 = T(i`j', c1)*FBorn0*(- marker99*replace_(i`j', i`k', summeP, summe, marker100, 1));
                if(occurs(marker99)==0) id once summeP(i`k')*summeP(i`j')*T(i`j', c1?)*FBorn0 = -summe(i`k')*T(i`k', c1)*FBorn0*marker99;
              endif;
            endif;
          #endif
        #enddo
        id marker100 = 0;
        id marker99 = 1;
        #call ColorOrder
*        if((count(summeP, 1)==1)&&(count(summe, 1)==0)) id summeP(i1?) = summe(i1);
      endif;
    endif;
    id marker`j' = 1;
  endif;
  #enddo
  .sort

  id once T(i1,N1_?)*T(i1,N2_?)*T(i2,N2_?)*T(i2,N3_?) = T(i1,N2_?)*T(i1,N1_?)*T(i2,N2_?)*T(i2,N3_?) + i_*cOlf(N1_?, N2_?, c1)*T(i1, c1)*T(i2,N2_?)*T(i2,N3_?);
  id once T(i1,N1_?)*T(i1,N2_?)*T(i2,N1_?)*T(i3,N3_?) = T(i1,N2_?)*T(i1,N1_?)*T(i2,N1_?)*T(i3,N3_?) + i_*cOlf(N1_?, N2_?, c2)*T(i1, c2)*T(i2,N1_?)*T(i3,N3_?);
  sum c1, c2;
  id once T(i1,N1_?)*T(i1,N2_?)*T(i2,N3_?)*T(i3,N1_?) = T(i1,N1_?)*T(i1,N2_?)*T(i2,N3_?)*T(i3,N1_?)*replace_(i2, i3, i3, i2);
  id once T(i1,N1_?)*T(i2,N2_?)*T(i3,N2_?)*T(i3,N3_?) = T(i1,N1_?)*T(i2,N2_?)*T(i3,N2_?)*T(i3,N3_?)*replace_(i2, i3, i3, i2);
  id once T(i1,N1_?)*T(i2,N2_?)*T(i3,N3_?)*cOlf(b,N3_?,N4_?)*cOlf(N1_?,N2_?,N4_?) = T(i1,N1_?)*T(i2,N2_?)*T(i3,N3_?)*cOlf(b,N3_?,N4_?)*cOlf(N1_?,N2_?,N4_?)*replace_(i2, i3, i3, i2);
  id once cOlf(cOli1?, cOli12?, cOli14?)*cOlf(cOli2?, cOli13?, cOli14?)*cOlf(cOli12?, cOli13?, cOli3?) = Cg/2*cOlf(cOli1, cOli2, cOli3);

  .sort

#enddo
#call ColorOrder
.sort


if(match(log(i1?, i1?)));
  print "%t";
  exit "Wrong log detected";
endif;

*if(count(summeP, 1)==0);
*  if((count(summe, 1)==1) && (count(summeP, 1) == 0))
*  id summe(i1)*T(i1, b)*j(i1) = summeP(i1)*summeP(i2)*i_*cOlf(b, c1, c2)*T(i1, c1)*T(i2, c2)*(-j(i1) - j(i1)*replace_(i1, i2, i2, i1))*(-1/Cg);
*  sum c1, c2;
*endif;

*id j(i2) = (jj(i1, i2) + j(i1));
*id j(i1) = 1/2*j(i1)*(1 + replace_(i1, i2, i2, i1));
#call ColorOrder
*id j(i2) = (jj(i1, i2) + j(i1));

if(count(ep, 1)>0) discard;
id log(q, i1) = log + log(i1, i2) - log(q, i2);
b FBorn0, FBorn1, FBorn2, a, cOlf, cOld, T, i_, summeP, summe, Cg, TF, nl, ep;
*b ep;
print+s;
.sort

.end


