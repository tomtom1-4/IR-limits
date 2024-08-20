definitions={l1 -> Log[s1i/(s1i + s2i)], l2 -> Log[s2j/(s1j + s2j)],
 l3 -> Log[(s12*sij)/(s1i*s2j)], L1 -> PolyLog[2, 1 - s2j/(s1j + s2j)],
 L2 -> PolyLog[2, 1 - s1i/(s1i + s2i)],
 L3 -> PolyLog[2, 1 - (s12*sij)/((s1i + s2i)*(s1j + s2j))]};

Jgg=(2*J*cOlf[c1, a, c]*cOlf[c2, b, c]*Ti[a]*Tj[b])/s12;

J=JHat*ghatT + Jtilde*gtilde +
 Jpm*(ghatT + (2*s12*TensorProduct[piT,piT])/(s1i*s2i)) +
 2*Jpp*s12*(TensorProduct[piT,pjT] - TensorProduct[pjT,piT])/(s1i*s2j) +
 JpmSwapped*(ghatT + (2*s12*TensorProduct[pjT, pjT])/(s1j*s2j));

JHat=(-3*(2*l2^2 + 2*l2*l3 + l3^2 + 2*(L2 + L3)) -
  (12*(-(s1j*s2i) + s1i*s2j + s12*sij))/(ep^2*s1i*s2j) -
  (6*l1^2*(-(s1j*s2i) + s1i*s2j + s12*sij))/(s1i*s2j) -
  (6*l3*(-(s1j*s2i) + s1i*s2j + s12*sij))/(ep*s1i*s2j) -
  (6*l1*(l2 + l3)*(-(s1j*s2i) + s1i*s2j + s12*sij))/(s1i*s2j) +
  (-4*CA*s1i*s1j*s2i*s2j + 4*CA*s1i^2*s2j^2 + 6*CA*L1*s1j*s2i*(s1i + s2i)*
     (s1j + s2j) + 6*CA*l2^2*s1j*s2i*(s1i + s2i)*(s1j + s2j) +
    6*CA*L2*s1j*s2i*(s1i + s2i)*(s1j + s2j) + 6*CA*l2*l3*s1j*s2i*(s1i + s2i)*
     (s1j + s2j) + 3*CA*l3^2*s1j*s2i*(s1i + s2i)*(s1j + s2j) +
    6*CA*L3*s1j*s2i*(s1i + s2i)*(s1j + s2j) - 6*CA*L1*s1i*(s1i + s2i)*s2j*
     (s1j + s2j) - 6*CA*L1*s12*(s1i + s2i)*(s1j + s2j)*sij -
    6*CA*l2^2*s12*(s1i + s2i)*(s1j + s2j)*sij - 6*CA*L2*s12*(s1i + s2i)*
     (s1j + s2j)*sij - 6*CA*l2*l3*s12*(s1i + s2i)*(s1j + s2j)*sij -
    3*CA*l3^2*s12*(s1i + s2i)*(s1j + s2j)*sij - 6*CA*L3*s12*(s1i + s2i)*
     (s1j + s2j)*sij + 8*nl*s1i*s2j*(s1j*s2i - s1i*s2j)*TF)/
   (CA*s1i*(s1i + s2i)*s2j*(s1j + s2j)))/6;

Jpp=-1/2*(4 + 2*ep*l3 + ep^2*(2*(l1^2 + l1*l2 + l2^2) + 2*(l1 + l2)*l3 + l3^2 +
     2*(L1 + L2 + L3)))/ep^2;

Jpm=-1/3*((12 + 6*ep*l3 + ep^2*(6*l1^2 + 6*l1*l3 + 3*l3^2 + 2*Pi^2))*s12*s2i*sij*
    (s2i*(s1j + s2j) - s12*sij) + s1i^2*(s1j + s2j)*
    ((12 + 6*ep*l3 + ep^2*(6*l1^2 + 6*(l1 + l2)*l3 + 3*l3^2 +
         2*(3*l2^2 + Pi^2)))*s2i*(s1j - s2j) - 6*ep^2*l1*(l1 + l3)*s12*sij) +
   s1i*((12 + 6*ep*l3 + ep^2*(6*l1^2 + 6*(l1 + l2)*l3 + 3*l3^2 +
         2*(3*l2^2 + Pi^2)))*s2i^2*(s1j - s2j)*(s1j + s2j) +
     2*s12*s2i*(12*s2j + ep*(-6*ep*l1*(l1 + l3)*s1j + 6*l3*s2j +
         ep*(6*l1^2 + 6*l1*l3 + 3*l3^2 + 2*Pi^2)*s2j))*sij +
     6*ep^2*l1*(l1 + l3)*s12^2*sij^2))/(ep^2*(s1i + s2i)*
   (s1j*s2i + s1i*s2j - s12*sij)*((s1i + s2i)*(s1j + s2j) - s12*sij));

JpmSwapped = Jpm /.{s1i -> s2j, s2j -> s1i, s2i -> s1j, s1j -> s2i, l1 -> l2, l2 -> l1};

