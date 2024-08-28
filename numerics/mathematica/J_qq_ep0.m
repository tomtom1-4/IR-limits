definitions={l1 -> Log[s1i/(s1i + s2i)], l2 -> Log[s2j/(s1j + s2j)], 
 l3 -> Log[(s12*sij)/(s1i*s2j)]};

Jqq=(2*J*cOlT[a, c2, c]*cOlT[c1, a, b]*Ti[b]*Tj[c])/s12;

J=(2*Jpm*piSlash)/(s1i + s2i) - (2*JpmSwapped*pjSlash)/(s1j + s2j);

Jpm=(18*CF*(2 + ep*(3 + 8*ep)) + 
  (CA*ep*((-66 + 18*l3 + 9*ep*(2*l1^2 + 2*(l1 + l2)*l3 + l3^2) + 
       2*ep*(-76 + 9*l2^2 + 3*Pi^2))*(s1i + s2i)*(s1j + s2j) - 
     (-66 + 18*l3 + ep*(-152 + 9*(2*l1 + l3)^2 + 6*Pi^2))*s12*sij))/
   ((s1i + s2i)*(s1j + s2j) - s12*sij) + 8*ep*(3 + 5*ep)*nl*TF)/(18*CA*ep^2);

JpmSwapped = Jpm /.{s1i -> s2j, s2j -> s1i, s2i -> s1j, s1j -> s2i, l1 -> l2, l2 -> l1};

