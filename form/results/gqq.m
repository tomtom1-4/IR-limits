den[x_]:= 1/x;
Sgqqij1:= ( + den[pjq1*pjq3 + pjq1*pjq2]*den[piq1*piq3*q2q3 + piq1*piq2*q2q3]
       * (  - 3/4*pipj^2
       )

       + den[piq3*pjq1 + piq2*pjq1]*den[piq1*piq3*q2q3 + piq1*piq2*q2q3] * ( 1/
         2*pipj*mi2
       )

       + den[piq1*pjq3 + piq1*pjq2]*den[piq1*piq3*q2q3 + piq1*piq2*q2q3] * ( 1/
         2*pipj*mi2
       )

       + den[piq1*pjq3 + piq1*pjq2]*den[piq1*piq3*q2q3^2 + piq1*piq2*q2q3^2]
       * (  - 1/2*piq3*pjq2*mi2
       - 1/2*piq2*pjq3*mi2
       )

       + den[pjq1]*den[pjq3 + pjq2]*den[piq1*piq3*q2q3^2 + piq1*piq2*q2q3^2]
       * ( 3/4*piq3*pjq2*pipj
       + 3/4*piq2*pjq3*pipj
       )

       + den[pjq1]*den[piq3 + piq2]*den[piq1*piq3*q2q3^2 + piq1*piq2*q2q3^2]
       * (  - piq2*piq3*pipj ));
Sgqqij2 = ( + den[piq3 + piq2]*den[piq3*pjq1 + piq2*pjq1]*den[2*piq3*q2q3 + 2*
      piq2*q2q3 + 2*piq1*q2q3] * (  - pipj*mi2
       )

       + den[piq3 + piq2]*den[piq1*pjq3 + piq1*pjq2]*den[2*piq3*q2q3 + 2*piq2*
      q2q3 + 2*piq1*q2q3] * ( pipj*mi2
       )

       + den[piq3*pjq1 + piq2*pjq1]*den[2*q2q3 + 2*q1q3 + 2*q1q2]*den[2*piq3*
      q2q3 + 2*piq2*q2q3 + 2*piq1*q2q3] * (  - 3*pjq3*mi2
       - 3*pjq2*mi2
       + 2*piq1*pipj
       )

       + den[piq3*pjq1 + piq2*pjq1]*den[2*q2q3^2 + 2*q1q3*q2q3 + 2*q1q2*q2q3]*
      den[2*piq3*q2q3 + 2*piq2*q2q3 + 2*piq1*q2q3] * (  - 4*piq3*q1q2*pipj
       - 2*piq3^2*pjq2
       - 4*piq2*q1q3*pipj
       + 6*piq2*piq3*pjq3
       + 6*piq2*piq3*pjq2
       - 2*piq2^2*pjq3
       + 2*piq1*piq3*pjq2
       + 2*piq1*piq2*pjq3
       )

       + den[piq3*pjq1 + piq2*pjq1]*den[2*q1q3*q2q3 + 2*q1q3^2 + 2*q1q2*q1q3]*
      den[2*piq3*q2q3 + 2*piq2*q2q3 + 2*piq1*q2q3] * (  - 2*pjq3*q2q3*mi2
       - pjq3*q1q2*mi2
       - pjq1*q2q3*mi2
       + 4*piq2*piq3*pjq3
       + 2*piq2*piq3*pjq1
       + 2*piq1*piq2*pjq3
       )

       + den[piq3*pjq1 + piq2*pjq1]*den[2*q1q2*q2q3 + 2*q1q2*q1q3 + 2*q1q2^2]*
      den[2*piq3*q2q3 + 2*piq2*q2q3 + 2*piq1*q2q3] * (  - 2*pjq2*q2q3*mi2
       - pjq2*q1q3*mi2
       - pjq1*q2q3*mi2
       + 4*piq2*piq3*pjq2
       + 2*piq2*piq3*pjq1
       + 2*piq1*piq3*pjq2
       )

       + den[piq1*pjq3 + piq1*pjq2]*den[2*q2q3 + 2*q1q3 + 2*q1q2]*den[2*piq3*
      q2q3 + 2*piq2*q2q3 + 2*piq1*q2q3] * ( pjq3*mi2
       + pjq2*mi2
       - 4*pjq1*mi2
       + 2*piq3*pipj
       + 2*piq2*pipj
       + 2*piq1*pipj
       )

       + den[piq1*pjq3 + piq1*pjq2]*den[2*q2q3^2 + 2*q1q3*q2q3 + 2*q1q2*q2q3]*
      den[2*piq3*q2q3 + 2*piq2*q2q3 + 2*piq1*q2q3] * ( 4*pjq3*q1q2*mi2
       + 4*pjq2*q1q3*mi2
       - 2*piq3^2*pjq2
       - 2*piq2*piq3*pjq3
       - 2*piq2*piq3*pjq2
       - 2*piq2^2*pjq3
       - 2*piq1*piq3*pjq2
       - 2*piq1*piq2*pjq3
       )

       + den[piq1*pjq3 + piq1*pjq2]*den[2*q1q3*q2q3 + 2*q1q3^2 + 2*q1q2*q1q3]*
      den[2*piq3*q2q3 + 2*piq2*q2q3 + 2*piq1*q2q3] * ( pjq3*q1q2*mi2
       - pjq1*q2q3*mi2
       + 2*piq3*q2q3*pipj
       - 2*piq3^2*pjq2
       - 2*piq2*piq3*pjq3
       + 2*piq1*q2q3*pipj
       - 2*piq1*piq3*pjq2
       - 2*piq1*piq2*pjq3
       )

       + den[piq1*pjq3 + piq1*pjq2]*den[2*q1q2*q2q3 + 2*q1q2*q1q3 + 2*q1q2^2]*
      den[2*piq3*q2q3 + 2*piq2*q2q3 + 2*piq1*q2q3] * ( pjq2*q1q3*mi2
       - pjq1*q2q3*mi2
       + 2*piq2*q2q3*pipj
       - 2*piq2*piq3*pjq2
       - 2*piq2^2*pjq3
       + 2*piq1*q2q3*pipj
       - 2*piq1*piq3*pjq2
       - 2*piq1*piq2*pjq3
       )

       + den[piq1]*den[piq3*pjq1 + piq2*pjq1]*den[2*piq3*q2q3 + 2*piq2*q2q3 + 
      2*piq1*q2q3] * ( pipj*mi2
       )

       + den[piq1]*den[piq1*pjq3 + piq1*pjq2]*den[2*piq3*q2q3 + 2*piq2*q2q3 + 
      2*piq1*q2q3] * (  - pipj*mi2
       )

       + den[piq1]*den[q2q3]*den[piq3*pjq1 + piq2*pjq1]*den[2*piq3*q2q3 + 2*
      piq2*q2q3 + 2*piq1*q2q3] * (  - 2*piq2*piq3*pipj
       )

       + den[piq1]*den[q2q3]*den[piq1*pjq3 + piq1*pjq2]*den[2*piq3*q2q3 + 2*
      piq2*q2q3 + 2*piq1*q2q3] * ( piq3*pjq2*mi2
       + piq2*pjq3*mi2
       )

       + den[q2q3]*den[piq3 + piq2]*den[piq3*pjq1 + piq2*pjq1]*den[2*piq3*q2q3
       + 2*piq2*q2q3 + 2*piq1*q2q3] * ( 2*piq2*piq3*pipj
       )

       + den[q2q3]*den[piq3 + piq2]*den[piq1*pjq3 + piq1*pjq2]*den[2*piq3*q2q3
       + 2*piq2*q2q3 + 2*piq1*q2q3] * (  - piq3*pjq2*mi2
       - piq2*pjq3*mi2 ));
Sgqqij3AB = ( + den[4*q2q3^2 + 8*q1q3*q2q3 + 4*q1q3^2 + 8*q1q2*q2q3 + 8*q1q2*
      q1q3 + 4*q1q2^2]*den[piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + 
      piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + piq1*pjq2 + piq1*pjq1] * (  - 16*
         pipj
       + 4*d*pipj
       )

       + den[q1q2*q1q3]*den[4*q2q3^2 + 8*q1q3*q2q3 + 4*q1q3^2 + 8*q1q2*q2q3 + 
      8*q1q2*q1q3 + 4*q1q2^2]*den[piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*
      pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + piq1*pjq2 + piq1*pjq1] * ( 8*
         q2q3^2*pipj
       - 8*piq3*pjq2*q2q3
       - 4*piq3*pjq1*q2q3
       - 8*piq2*pjq3*q2q3
       - 4*piq2*pjq1*q2q3
       - 4*piq1*pjq3*q2q3
       - 4*piq1*pjq2*q2q3
       - 16*piq1*pjq1*q2q3
       + 4*d*piq1*pjq1*q2q3
       )

       + den[q1q2]*den[4*q2q3^2 + 8*q1q3*q2q3 + 4*q1q3^2 + 8*q1q2*q2q3 + 8*
      q1q2*q1q3 + 4*q1q2^2]*den[piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3
       + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + piq1*pjq2 + piq1*pjq1] * ( 8*q2q3
         *pipj
       - 4*q1q3*pipj
       - 4*piq3*pjq2
       + 4*piq3*pjq1
       - 4*piq2*pjq3
       + 8*piq2*pjq2
       + 8*piq2*pjq1
       + 4*piq1*pjq3
       + 8*piq1*pjq2
       + 2*d*q1q3*pipj
       - 2*d*piq3*pjq1
       - 2*d*piq2*pjq1
       - 2*d*piq1*pjq3
       - 2*d*piq1*pjq2
       )

       + den[q1q3]*den[4*q2q3^2 + 8*q1q3*q2q3 + 4*q1q3^2 + 8*q1q2*q2q3 + 8*
      q1q2*q1q3 + 4*q1q2^2]*den[piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3
       + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + piq1*pjq2 + piq1*pjq1] * ( 8*q2q3
         *pipj
       - 4*q1q2*pipj
       + 8*piq3*pjq3
       - 4*piq3*pjq2
       + 8*piq3*pjq1
       - 4*piq2*pjq3
       + 4*piq2*pjq1
       + 8*piq1*pjq3
       + 4*piq1*pjq2
       + 2*d*q1q2*pipj
       - 2*d*piq3*pjq1
       - 2*d*piq2*pjq1
       - 2*d*piq1*pjq3
       - 2*d*piq1*pjq2 ));
Sgqqij3NAB = ( + den[4*q2q3]*den[piq3 + piq2]*den[pjq3 + pjq2]*den[piq3*pjq3
       + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3
       + piq1*pjq2 + piq1*pjq1] * (  - 4*pipj^2
       )

       + den[q2q3^2]*den[4*q2q3^2 + 8*q1q3*q2q3 + 4*q1q3^2 + 8*q1q2*q2q3 + 8*
      q1q2*q1q3 + 4*q1q2^2]*den[piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3
       + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + piq1*pjq2 + piq1*pjq1] * ( 32*
         q1q2*q1q3*pipj
       )

       + den[q2q3^2]*den[pjq3 + pjq2]*den[2*q2q3 + 2*q1q3 + 2*q1q2]*den[piq3*
      pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*
      pjq3 + piq1*pjq2 + piq1*pjq1] * (  - 4*pjq3*q1q2*pipj
       - 4*pjq2*q1q3*pipj
       + 4*piq3*pjq2^2
       + 4*piq2*pjq3^2
       + 4*piq1*pjq2*pjq3
       )

       + den[q2q3^2]*den[piq3 + piq2]*den[2*q2q3 + 2*q1q3 + 2*q1q2]*den[piq3*
      pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*
      pjq3 + piq1*pjq2 + piq1*pjq1] * (  - 4*piq3*q1q2*pipj
       + 4*piq3^2*pjq2
       - 4*piq2*q1q3*pipj
       + 4*piq2*piq3*pjq1
       + 4*piq2^2*pjq3
       )

       + den[2*q2q3^2]*den[piq3 + piq2]*den[pjq3 + pjq2]*den[piq3*pjq3 + piq3*
      pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + piq1*
      pjq2 + piq1*pjq1] * ( 2*piq3*pjq2*pipj
       + 2*piq2*pjq3*pipj
       )

       + den[4*q2q3^2 + 8*q1q3*q2q3 + 4*q1q3^2 + 8*q1q2*q2q3 + 8*q1q2*q1q3 + 4
      *q1q2^2]*den[piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2
       + piq2*pjq1 + piq1*pjq3 + piq1*pjq2 + piq1*pjq1] * ( 16*pipj
       - 4*d*pipj
       )

       + den[q1q3*q2q3]*den[pjq3 + pjq2]*den[2*q2q3 + 2*q1q3 + 2*q1q2]*den[
      piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + 
      piq1*pjq3 + piq1*pjq2 + piq1*pjq1] * (  - 2*pjq3*q1q2*pipj
       + piq3*q1q2*mi2
       + 2*piq3*pjq2*pjq3
       + 2*piq2*pjq3^2
       + 2*piq2*pjq1*pjq3
       + 2*piq1*pjq2*pjq3
       )

       + den[q1q3*q2q3]*den[piq3 + piq2]*den[2*q2q3 + 2*q1q3 + 2*q1q2]*den[
      piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + 
      piq1*pjq3 + piq1*pjq2 + piq1*pjq1] * ( pjq3*q1q2*mi2
       - 2*piq3*q1q2*pipj
       + 2*piq3^2*pjq2
       + 2*piq2*piq3*pjq3
       + 2*piq2*piq3*pjq1
       + 2*piq1*piq3*pjq2
       )

       + den[q1q2*q2q3]*den[pjq3 + pjq2]*den[2*q2q3 + 2*q1q3 + 2*q1q2]*den[
      piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + 
      piq1*pjq3 + piq1*pjq2 + piq1*pjq1] * (  - 2*pjq2*q1q3*pipj
       + 2*piq3*pjq2^2
       + 2*piq3*pjq1*pjq2
       + piq2*q1q3*mi2
       + 2*piq2*pjq2*pjq3
       + 2*piq1*pjq2*pjq3
       )

       + den[q1q2*q2q3]*den[piq3 + piq2]*den[2*q2q3 + 2*q1q3 + 2*q1q2]*den[
      piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + 
      piq1*pjq3 + piq1*pjq2 + piq1*pjq1] * ( pjq2*q1q3*mi2
       - 2*piq2*q1q3*pipj
       + 2*piq2*piq3*pjq2
       + 2*piq2*piq3*pjq1
       + 2*piq2^2*pjq3
       + 2*piq1*piq2*pjq3
       )

       + den[q1q2*q1q3]*den[4*q2q3^2 + 8*q1q3*q2q3 + 4*q1q3^2 + 8*q1q2*q2q3 + 
      8*q1q2*q1q3 + 4*q1q2^2]*den[piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*
      pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + piq1*pjq2 + piq1*pjq1] * ( 
          - 8*q2q3^2*pipj
       + 8*piq3*pjq2*q2q3
       + 4*piq3*pjq1*q2q3
       + 8*piq2*pjq3*q2q3
       + 4*piq2*pjq1*q2q3
       + 4*piq1*pjq3*q2q3
       + 4*piq1*pjq2*q2q3
       + 16*piq1*pjq1*q2q3
       - 4*d*piq1*pjq1*q2q3
       )

       + den[piq1]*den[4*q2q3]*den[pjq3 + pjq2]*den[piq3*pjq3 + piq3*pjq2 + 
      piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + piq1*pjq2 + 
      piq1*pjq1] * ( 4*pipj^2
       )

       + den[piq1]*den[q2q3^2]*den[2*q2q3 + 2*q1q3 + 2*q1q2]*den[piq3*pjq3 + 
      piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + 
      piq1*pjq2 + piq1*pjq1] * ( 4*piq3*q1q2*pipj
       - 4*piq3^2*pjq2
       + 4*piq2*q1q3*pipj
       - 4*piq2*piq3*pjq1
       - 4*piq2^2*pjq3
       )

       + den[piq1]*den[2*q2q3^2]*den[pjq3 + pjq2]*den[piq3*pjq3 + piq3*pjq2 + 
      piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + piq1*pjq2 + 
      piq1*pjq1] * (  - 2*piq3*pjq2*pipj
       - 2*piq2*pjq3*pipj
       )

       + den[piq1]*den[q1q3*q2q3]*den[2*q2q3 + 2*q1q3 + 2*q1q2]*den[piq3*pjq3
       + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3
       + piq1*pjq2 + piq1*pjq1] * (  - pjq3*q1q2*mi2
       + 2*piq3*q1q2*pipj
       - 2*piq3^2*pjq2
       - 2*piq2*piq3*pjq3
       - 2*piq2*piq3*pjq1
       - 2*piq1*piq3*pjq2
       )

       + den[piq1]*den[q1q2*q2q3]*den[2*q2q3 + 2*q1q3 + 2*q1q2]*den[piq3*pjq3
       + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3
       + piq1*pjq2 + piq1*pjq1] * (  - pjq2*q1q3*mi2
       + 2*piq2*q1q3*pipj
       - 2*piq2*piq3*pjq2
       - 2*piq2*piq3*pjq1
       - 2*piq2^2*pjq3
       - 2*piq1*piq2*pjq3
       )

       + den[piq1]*den[pjq1]*den[4*q2q3]*den[piq3*pjq3 + piq3*pjq2 + piq3*pjq1
       + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + piq1*pjq2 + piq1*pjq1
      ] * (  - 4*pipj^2
       )

       + den[piq1]*den[pjq1]*den[2*q2q3^2]*den[piq3*pjq3 + piq3*pjq2 + piq3*
      pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + piq1*pjq2 + piq1*
      pjq1] * ( 2*piq3*pjq2*pipj
       + 2*piq2*pjq3*pipj
       )

       + den[piq1]*den[q1q2]*den[2*q2q3 + 2*q1q3 + 2*q1q2]*den[piq3*pjq3 + 
      piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + 
      piq1*pjq2 + piq1*pjq1] * ( pjq1*mi2
       + 2*piq2*pipj
       )

       + den[piq1]*den[q1q3]*den[2*q2q3 + 2*q1q3 + 2*q1q2]*den[piq3*pjq3 + 
      piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + 
      piq1*pjq2 + piq1*pjq1] * ( pjq1*mi2
       + 2*piq3*pipj
       )

       + den[piq1]*den[q2q3]*den[2*q2q3 + 2*q1q3 + 2*q1q2]*den[piq3*pjq3 + 
      piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + 
      piq1*pjq2 + piq1*pjq1] * (  - pjq3*mi2
       - pjq2*mi2
       + 2*pjq1*mi2
       + 4*piq3*pipj
       + 4*piq2*pipj
       - 4*piq1*pipj
       )

       + den[pjq1]*den[4*q2q3]*den[piq3 + piq2]*den[piq3*pjq3 + piq3*pjq2 + 
      piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + piq1*pjq2 + 
      piq1*pjq1] * ( 4*pipj^2
       )

       + den[pjq1]*den[q2q3^2]*den[2*q2q3 + 2*q1q3 + 2*q1q2]*den[piq3*pjq3 + 
      piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + 
      piq1*pjq2 + piq1*pjq1] * ( 4*pjq3*q1q2*pipj
       + 4*pjq2*q1q3*pipj
       - 4*piq3*pjq2^2
       - 4*piq2*pjq3^2
       - 4*piq1*pjq2*pjq3
       )

       + den[pjq1]*den[2*q2q3^2]*den[piq3 + piq2]*den[piq3*pjq3 + piq3*pjq2 + 
      piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + piq1*pjq2 + 
      piq1*pjq1] * (  - 2*piq3*pjq2*pipj
       - 2*piq2*pjq3*pipj
       )

       + den[pjq1]*den[q1q3*q2q3]*den[2*q2q3 + 2*q1q3 + 2*q1q2]*den[piq3*pjq3
       + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3
       + piq1*pjq2 + piq1*pjq1] * ( 2*pjq3*q1q2*pipj
       - piq3*q1q2*mi2
       - 2*piq3*pjq2*pjq3
       - 2*piq2*pjq3^2
       - 2*piq2*pjq1*pjq3
       - 2*piq1*pjq2*pjq3
       )

       + den[pjq1]*den[q1q2*q2q3]*den[2*q2q3 + 2*q1q3 + 2*q1q2]*den[piq3*pjq3
       + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3
       + piq1*pjq2 + piq1*pjq1] * ( 2*pjq2*q1q3*pipj
       - 2*piq3*pjq2^2
       - 2*piq3*pjq1*pjq2
       - piq2*q1q3*mi2
       - 2*piq2*pjq2*pjq3
       - 2*piq1*pjq2*pjq3
       )

       + den[pjq1]*den[q1q2]*den[2*q2q3 + 2*q1q3 + 2*q1q2]*den[piq3*pjq3 + 
      piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + 
      piq1*pjq2 + piq1*pjq1] * ( 2*pjq2*pipj
       + piq1*mi2
       )

       + den[pjq1]*den[q1q3]*den[2*q2q3 + 2*q1q3 + 2*q1q2]*den[piq3*pjq3 + 
      piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + 
      piq1*pjq2 + piq1*pjq1] * ( 2*pjq3*pipj
       + piq1*mi2
       )

       + den[pjq1]*den[q2q3]*den[2*q2q3 + 2*q1q3 + 2*q1q2]*den[piq3*pjq3 + 
      piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + 
      piq1*pjq2 + piq1*pjq1] * ( 4*pjq3*pipj
       + 4*pjq2*pipj
       - 4*pjq1*pipj
       - piq3*mi2
       - piq2*mi2
       + 2*piq1*mi2
       )

       + den[q1q2]*den[2*q2q3 + 2*q1q3 + 2*q1q2]*den[piq3*pjq3 + piq3*pjq2 + 
      piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + piq1*pjq2 + 
      piq1*pjq1] * (  - 8*pipj
       )

       + den[q1q2]*den[4*q2q3^2 + 8*q1q3*q2q3 + 4*q1q3^2 + 8*q1q2*q2q3 + 8*
      q1q2*q1q3 + 4*q1q2^2]*den[piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3
       + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + piq1*pjq2 + piq1*pjq1] * (  - 8*
         q2q3*pipj
       - 4*q1q3*pipj
       + 20*piq3*pjq2
       + 8*piq3*pjq1
       + 20*piq2*pjq3
       - 8*piq2*pjq2
       - 4*piq2*pjq1
       + 8*piq1*pjq3
       - 4*piq1*pjq2
       + 24*piq1*pjq1
       + 2*d*q1q3*pipj
       + 4*d*piq2*pjq1
       + 4*d*piq1*pjq2
       - 4*d*piq1*pjq1
       )

       + den[q1q2]*den[pjq3 + pjq2]*den[2*q2q3 + 2*q1q3 + 2*q1q2]*den[piq3*
      pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*
      pjq3 + piq1*pjq2 + piq1*pjq1] * (  - 2*pjq2*pipj
       - piq1*mi2
       )

       + den[q1q2]*den[piq3 + piq2]*den[2*q2q3 + 2*q1q3 + 2*q1q2]*den[piq3*
      pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*
      pjq3 + piq1*pjq2 + piq1*pjq1] * (  - pjq1*mi2
       - 2*piq2*pipj
       )

       + den[q1q2]*den[q2q3]*den[4*q2q3^2 + 8*q1q3*q2q3 + 4*q1q3^2 + 8*q1q2*
      q2q3 + 8*q1q2*q1q3 + 4*q1q2^2]*den[piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + 
      piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + piq1*pjq2 + piq1*pjq1]
       * ( 12*piq3*pjq2*q1q3
       + 12*piq2*pjq3*q1q3
       - 8*piq2*pjq2*q1q3
       - 12*piq2*pjq1*q1q3
       - 12*piq1*pjq2*q1q3
       - 2*d*piq3*pjq2*q1q3
       - 2*d*piq2*pjq3*q1q3
       - 4*d*piq2*pjq2*q1q3
       + 2*d*piq2*pjq1*q1q3
       + 2*d*piq1*pjq2*q1q3
       )

       + den[q1q3]*den[2*q2q3 + 2*q1q3 + 2*q1q2]*den[piq3*pjq3 + piq3*pjq2 + 
      piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + piq1*pjq2 + 
      piq1*pjq1] * (  - 8*pipj
       )

       + den[q1q3]*den[4*q2q3^2 + 8*q1q3*q2q3 + 4*q1q3^2 + 8*q1q2*q2q3 + 8*
      q1q2*q1q3 + 4*q1q2^2]*den[piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3
       + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + piq1*pjq2 + piq1*pjq1] * (  - 8*
         q2q3*pipj
       - 4*q1q2*pipj
       - 8*piq3*pjq3
       + 20*piq3*pjq2
       - 4*piq3*pjq1
       + 20*piq2*pjq3
       + 8*piq2*pjq1
       - 4*piq1*pjq3
       + 8*piq1*pjq2
       + 24*piq1*pjq1
       + 2*d*q1q2*pipj
       + 4*d*piq3*pjq1
       + 4*d*piq1*pjq3
       - 4*d*piq1*pjq1
       )

       + den[q1q3]*den[pjq3 + pjq2]*den[2*q2q3 + 2*q1q3 + 2*q1q2]*den[piq3*
      pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*
      pjq3 + piq1*pjq2 + piq1*pjq1] * (  - 2*pjq3*pipj
       - piq1*mi2
       )

       + den[q1q3]*den[piq3 + piq2]*den[2*q2q3 + 2*q1q3 + 2*q1q2]*den[piq3*
      pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*
      pjq3 + piq1*pjq2 + piq1*pjq1] * (  - pjq1*mi2
       - 2*piq3*pipj
       )

       + den[q1q3]*den[q2q3]*den[4*q2q3^2 + 8*q1q3*q2q3 + 4*q1q3^2 + 8*q1q2*
      q2q3 + 8*q1q2*q1q3 + 4*q1q2^2]*den[piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + 
      piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + piq1*pjq2 + piq1*pjq1]
       * (  - 8*piq3*pjq3*q1q2
       + 12*piq3*pjq2*q1q2
       - 12*piq3*pjq1*q1q2
       + 12*piq2*pjq3*q1q2
       - 12*piq1*pjq3*q1q2
       - 4*d*piq3*pjq3*q1q2
       - 2*d*piq3*pjq2*q1q2
       + 2*d*piq3*pjq1*q1q2
       - 2*d*piq2*pjq3*q1q2
       + 2*d*piq1*pjq3*q1q2
       )

       + den[q2q3]*den[4*q2q3^2 + 8*q1q3*q2q3 + 4*q1q3^2 + 8*q1q2*q2q3 + 8*
      q1q2*q1q3 + 4*q1q2^2]*den[piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3
       + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + piq1*pjq2 + piq1*pjq1] * ( 16*
         q1q3*pipj
       + 16*q1q2*pipj
       - 16*piq3*pjq3
       + 16*piq3*pjq2
       - 12*piq3*pjq1
       + 16*piq2*pjq3
       - 16*piq2*pjq2
       - 12*piq2*pjq1
       - 12*piq1*pjq3
       - 12*piq1*pjq2
       + 8*piq1*pjq1
       + 2*d*piq3*pjq1
       + 2*d*piq2*pjq1
       + 2*d*piq1*pjq3
       + 2*d*piq1*pjq2
       - 4*d*piq1*pjq1
       )

       + den[q2q3]*den[pjq3 + pjq2]*den[2*q2q3 + 2*q1q3 + 2*q1q2]*den[piq3*
      pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*
      pjq3 + piq1*pjq2 + piq1*pjq1] * (  - 4*pjq3*pipj
       - 4*pjq2*pipj
       + 4*pjq1*pipj
       + piq3*mi2
       + piq2*mi2
       - 2*piq1*mi2
       )

       + den[q2q3]*den[piq3 + piq2]*den[2*q2q3 + 2*q1q3 + 2*q1q2]*den[piq3*
      pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*
      pjq3 + piq1*pjq2 + piq1*pjq1] * ( pjq3*mi2
       + pjq2*mi2
       - 2*pjq1*mi2
       - 4*piq3*pipj
       - 4*piq2*pipj
       + 4*piq1*pipj ));
Sgqqija = ( + den[q1q3*q2q3]*den[2*piq1*pjq3*q2q3*pkq3 + 2*piq1*pjq3*q2q3*pkq2
       + 2*piq1*pjq3*q2q3*pkq1 + 2*piq1*pjq3*q1q3*pkq3 + 2*piq1*pjq3*q1q3*pkq2
       + 2*piq1*pjq3*q1q3*pkq1 + 2*piq1*pjq3*q1q2*pkq3 + 2*piq1*pjq3*q1q2*pkq2
       + 2*piq1*pjq3*q1q2*pkq1 + 2*piq1*pjq2*q2q3*pkq3 + 2*piq1*pjq2*q2q3*pkq2
       + 2*piq1*pjq2*q2q3*pkq1 + 2*piq1*pjq2*q1q3*pkq3 + 2*piq1*pjq2*q1q3*pkq2
       + 2*piq1*pjq2*q1q3*pkq1 + 2*piq1*pjq2*q1q2*pkq3 + 2*piq1*pjq2*q1q2*pkq2
       + 2*piq1*pjq2*q1q2*pkq1] * (  - q1q2*pipj*pkq3
       + pjq3*q1q2*pipk
       + piq3*q1q2*pjpk
       - 2*piq3*pjq3*pkq2
       - 2*piq3*pjq2*pkq3
       - piq3*pjq2*pkq1
       - piq3*pjq1*pkq2
       - piq2*pjq3*pkq1
       + piq2*pjq1*pkq3
       - piq1*pjq3*pkq2
       - piq1*pjq2*pkq3
       )

       + den[q1q2*q2q3]*den[2*piq1*pjq3*q2q3*pkq3 + 2*piq1*pjq3*q2q3*pkq2 + 2*
      piq1*pjq3*q2q3*pkq1 + 2*piq1*pjq3*q1q3*pkq3 + 2*piq1*pjq3*q1q3*pkq2 + 2*
      piq1*pjq3*q1q3*pkq1 + 2*piq1*pjq3*q1q2*pkq3 + 2*piq1*pjq3*q1q2*pkq2 + 2*
      piq1*pjq3*q1q2*pkq1 + 2*piq1*pjq2*q2q3*pkq3 + 2*piq1*pjq2*q2q3*pkq2 + 2*
      piq1*pjq2*q2q3*pkq1 + 2*piq1*pjq2*q1q3*pkq3 + 2*piq1*pjq2*q1q3*pkq2 + 2*
      piq1*pjq2*q1q3*pkq1 + 2*piq1*pjq2*q1q2*pkq3 + 2*piq1*pjq2*q1q2*pkq2 + 2*
      piq1*pjq2*q1q2*pkq1] * ( q1q3*pipj*pkq2
       - pjq2*q1q3*pipk
       + piq3*pjq2*pkq1
       - piq3*pjq1*pkq2
       - piq2*q1q3*pjpk
       + 2*piq2*pjq3*pkq2
       + piq2*pjq3*pkq1
       + 2*piq2*pjq2*pkq3
       + piq2*pjq1*pkq3
       + piq1*pjq3*pkq2
       + piq1*pjq2*pkq3
       )

       + den[q1q2]*den[2*piq1*pjq3*q2q3*pkq3 + 2*piq1*pjq3*q2q3*pkq2 + 2*piq1*
      pjq3*q2q3*pkq1 + 2*piq1*pjq3*q1q3*pkq3 + 2*piq1*pjq3*q1q3*pkq2 + 2*piq1*
      pjq3*q1q3*pkq1 + 2*piq1*pjq3*q1q2*pkq3 + 2*piq1*pjq3*q1q2*pkq2 + 2*piq1*
      pjq3*q1q2*pkq1 + 2*piq1*pjq2*q2q3*pkq3 + 2*piq1*pjq2*q2q3*pkq2 + 2*piq1*
      pjq2*q2q3*pkq1 + 2*piq1*pjq2*q1q3*pkq3 + 2*piq1*pjq2*q1q3*pkq2 + 2*piq1*
      pjq2*q1q3*pkq1 + 2*piq1*pjq2*q1q2*pkq3 + 2*piq1*pjq2*q1q2*pkq2 + 2*piq1*
      pjq2*q1q2*pkq1] * (  - pipj*pkq1
       + pjq1*pipk
       - 2*piq2*pjpk
       - piq1*pjpk
       )

       + den[q1q3]*den[2*piq1*pjq3*q2q3*pkq3 + 2*piq1*pjq3*q2q3*pkq2 + 2*piq1*
      pjq3*q2q3*pkq1 + 2*piq1*pjq3*q1q3*pkq3 + 2*piq1*pjq3*q1q3*pkq2 + 2*piq1*
      pjq3*q1q3*pkq1 + 2*piq1*pjq3*q1q2*pkq3 + 2*piq1*pjq3*q1q2*pkq2 + 2*piq1*
      pjq3*q1q2*pkq1 + 2*piq1*pjq2*q2q3*pkq3 + 2*piq1*pjq2*q2q3*pkq2 + 2*piq1*
      pjq2*q2q3*pkq1 + 2*piq1*pjq2*q1q3*pkq3 + 2*piq1*pjq2*q1q3*pkq2 + 2*piq1*
      pjq2*q1q3*pkq1 + 2*piq1*pjq2*q1q2*pkq3 + 2*piq1*pjq2*q1q2*pkq2 + 2*piq1*
      pjq2*q1q2*pkq1] * ( pipj*pkq1
       - pjq1*pipk
       + 2*piq3*pjpk
       + piq1*pjpk
       )

       + den[q2q3]*den[2*piq1*pjq3*q2q3*pkq3 + 2*piq1*pjq3*q2q3*pkq2 + 2*piq1*
      pjq3*q2q3*pkq1 + 2*piq1*pjq3*q1q3*pkq3 + 2*piq1*pjq3*q1q3*pkq2 + 2*piq1*
      pjq3*q1q3*pkq1 + 2*piq1*pjq3*q1q2*pkq3 + 2*piq1*pjq3*q1q2*pkq2 + 2*piq1*
      pjq3*q1q2*pkq1 + 2*piq1*pjq2*q2q3*pkq3 + 2*piq1*pjq2*q2q3*pkq2 + 2*piq1*
      pjq2*q2q3*pkq1 + 2*piq1*pjq2*q1q3*pkq3 + 2*piq1*pjq2*q1q3*pkq2 + 2*piq1*
      pjq2*q1q3*pkq1 + 2*piq1*pjq2*q1q2*pkq3 + 2*piq1*pjq2*q1q2*pkq2 + 2*piq1*
      pjq2*q1q2*pkq1] * (  - pipj*pkq3
       + pipj*pkq2
       - pjq3*pipk
       + pjq2*pipk
       + piq3*pjpk
       - piq2*pjpk ));
