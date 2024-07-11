#procedure FullSimplify
.sort
PolyRatFun rat;

#do x = {s12, s13, s14, s15, s23, s24, s25, s34, s35, s45, m1, m2, si1, si2, si3, sj1, sj2, sj3, sa1, sa2, sa3, sb1, sb2, sb3, mi, mj, ma, mb, marker}
    b `x', rat;
    .sort
    keep brackets;
    id 1/`x' = rat(1, `x');
    id `x'^m? = rat(`x'^m, 1);
#enddo
.sort
id den(s12?) = rat(1, s12);
.sort

PolyratFun;
.sort
id rat(s12?, s13?) = s12 * den(s13);
id den(1) = 1;

#endprocedure