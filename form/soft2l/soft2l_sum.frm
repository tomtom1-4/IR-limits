#-
#: IncDir procedures
#: SmallExtension 100M
*#: MaxTermSize    10M
#: WorkSpace      1G
Off Statistics;
format 160;

#include declarations.h

* get the number of diagrams
#include diagrams/diags.out
.sort
drop;
.sort

#write "`ndiags' diagram(s) for the soft current"
* Load the diagrams computed in soft2l_diags.frm
#do i=1,`ndiags'
    load results/diags/d`i'.sav;
#enddo
* Build sum
g diags2l = d1+...+d`ndiags';
l test = DS(k2)^2*DS( - p3 + k2)*DS(k1 + k2)*Eik(k2,p2)*Eik( - p3 + k2,p1) * (
          - 2*T1(cOli1)*T2(cOli2)*cOlf(c3,cOli1,cOli2)*i_*s12*s13*nl*TF*sumij);
.sort

#call Color

* project onto formfactor
mul replace_(E3, p1);
#call Kinematics

#call NormalizePropagators
.sort
#call NormalizePropagators
#call NormalizePropagators
#call Kinematics
.sort
*id T1(cOli1?) = 1;
*id T2(cOli1?) = 1;
*id cOlf(cOli1?, cOli2?, cOli3?) = 1;
*id s12^m? = 1;
*id s13^m? = 1;
*id s23^m? = 1;
*id sumij = 1;
*id i_ = 1;
*id TF = 1;
*id nl = 1;
*id CA = 1;
*id d = 1;
*id p1?.p2? = 1;
mul replace_(N1_?, cOli1, N2_?, cOli2, N3_?, cOli3, N4_?, cOli4, N5_?, cOli5);
.sort
if(count(Eik, 1) >= 5) discard;
if((match(Eik(p3?, p1))==0) && (match(p1.k1)==0) && (match(p1.k2)==0)) discard;
if((match(Eik(p3?, p2))==0) && (match(p2.k1)==0) && (match(p2.k2)==0)) discard;
*#call Scaleless
.sort
format mathematica;
bracket DS, Eik;
#write <results/integrals.m> "diags2l = (%+E);", diags2l
.sort
if(match(nl)==0) discard;
id DS( - p3 + k2) = DS(-p3 + k2)*replace_(k2, -k2);
id Eik(k2 - p3, p2?) = Eik(k2 - p3, p2)*replace_(k2, -k2);
id DS(k1 - k2) = DS(k1 - k2)*replace_(k1, -k1);
#call NormalizePropagators
#message "Sum of all diagrams"
bracket Eik, DS, k1, k2;
print[] ;
.sort
.store
save results/diags2l.sav diags2l;

.end

.end