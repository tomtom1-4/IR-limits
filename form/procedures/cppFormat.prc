#procedure cppFormat

* rewrite to c++ format
s piq1, piq2, piq3, pjq1, pjq2, pjq3, q1q2, q1q3, q2q3, pipj, mi2, pkq1, pkq2, pkq3, pipk, pjpk;
id pi.q1 = piq1;
id pi.q2 = piq2;
id pi.q3 = piq3;
id pj.q1 = pjq1;
id pj.q2 = pjq2;
id pj.q3 = pjq3;
id pa.q1 = pkq1;
id pa.q2 = pkq2;
id pa.q3 = pkq3;
id q1.q2 = q1q2;
id q1.q3 = q1q3;
id q2.q3 = q2q3;
id pi.pj = pipj;
id pi.pa = pipk;
id pj.pa = pjpk;
id pi.pi = mi2;
argument den;
  id pi.q1 = piq1;
  id pi.q2 = piq2;
  id pi.q3 = piq3;
  id pj.q1 = pjq1;
  id pj.q2 = pjq2;
  id pj.q3 = pjq3;
  id pa.q1 = pkq1;
  id pa.q2 = pkq2;
  id pa.q3 = pkq3;
  id q1.q2 = q1q2;
  id q1.q3 = q1q3;
  id q2.q3 = q2q3;
  id pi.pj = pipj;
  id pi.pa = pipk;
  id pj.pa = pjpk;
  id pi.pi = mi2;
endargument;
*id summe(i?) = 1;
*id T(i?, cOli1?) = 1;
id mi^2 = mi2;

#endprocedure