#-
#: IncDir procedures
#: SmallExtension 200M
#: WorkSpace      4G
#: MaxTermSize 4M
Off Statistics;

Format 255;


s s12, s23, s13, s14, s24, s34, s123, s124, s1234, s134, s234, z1, z2, z3, z4, z12, z14, z23, z13, z34, z134, z123, z124, z234, z24, D, ep;
s gg;
#do i=1,4
  #do j=1,4
  s k`i'k`j';
  #enddo
#enddo
s CF, CA;
*v E0, E0conj, p1, p2, p3, p4, p12, p123;

#ifndef `split'
  #write "Please define which splitting function you would like to calculate."
  #write "form -Dsplit=... quadruple_collinear.frm"
  #write " (i)   qbp_qp_g_q"
  #write " (ii)  qb_p_g_q"
  #write " (iii) g_g_g_q"
  #write " (iv)  qbp_qp_qb_q"
  .end
#endif

#if (("`split'" == "qbp_qp_g_q") || ("`split'" == "qb_q_g_q") || ("`split'" == "g_g_g_q"))
  #include quark_parent/P_`split'.m
  id D = 4;
  .sort
  Format O1;
  .sort
  Format C;
  #write <quark_parent/P0_`split'.cpp> "#include \"../splitting_functions.hpp\""
  #write <quark_parent/P0_`split'.cpp> "std::vector<std::vector<std::complex<double>>> P0_`split'(LV<double> p1, LV<double> p2, LV<double> p3, LV<double> p4) {"
  #write <quark_parent/P0_`split'.cpp> "  double CA = C_A;"
  #write <quark_parent/P0_`split'.cpp> "  double CF = C_F;"
  #write <quark_parent/P0_`split'.cpp> "  double s12 = 2.*p1*p2;"
  #write <quark_parent/P0_`split'.cpp> "  double s13 = 2.*p1*p3;"
  #write <quark_parent/P0_`split'.cpp> "  double s14 = 2.*p1*p4;"
  #write <quark_parent/P0_`split'.cpp> "  double s23 = 2.*p2*p3;"
  #write <quark_parent/P0_`split'.cpp> "  double s24 = 2.*p2*p4;"
  #write <quark_parent/P0_`split'.cpp> "  double s34 = 2.*p3*p4;"
  #write <quark_parent/P0_`split'.cpp> "  double s123 = s12 + s13 + s23;"
  #write <quark_parent/P0_`split'.cpp> "  double s134 = s13 + s14 + s34;"
  #write <quark_parent/P0_`split'.cpp> "  double s124 = s12 + s14 + s24;"
  #write <quark_parent/P0_`split'.cpp> "  double s234 = s23 + s24 + s34;"
  #write <quark_parent/P0_`split'.cpp> "  double s1234 = s12 + s13 + s14 + s23 + s24 + s34;"
  #write <quark_parent/P0_`split'.cpp> "  LV<double> p = p1 + p2 + p3 + p4;"
  #write <quark_parent/P0_`split'.cpp> "  double z1 = p1.components[0]/p.components[0];"
  #write <quark_parent/P0_`split'.cpp> "  double z2 = p2.components[0]/p.components[0];"
  #write <quark_parent/P0_`split'.cpp> "  double z3 = p3.components[0]/p.components[0];"
  #write <quark_parent/P0_`split'.cpp> "  double z4 = p4.components[0]/p.components[0];"
  #write <quark_parent/P0_`split'.cpp> "  double z12 = z1 + z2;"
  #write <quark_parent/P0_`split'.cpp> "  double z14 = z1 + z4;"
  #write <quark_parent/P0_`split'.cpp> "  double z23 = z2 + z3;"
  #write <quark_parent/P0_`split'.cpp> "  double z13 = z1 + z3;"
  #write <quark_parent/P0_`split'.cpp> "  double z123 = z1 + z2 + z3;"
  #write <quark_parent/P0_`split'.cpp> "  double z134 = z1 + z3 + z4;"
  #write <quark_parent/P0_`split'.cpp> "  double P0`split' = 0.;"

  ExtraSymbols array,Z;
  Format C;
  #optimize P
  #write <quark_parent/P0_`split'.cpp> "  double Z[`optimmaxvar_'+1];"
  #write <quark_parent/P0_`split'.cpp> "%O"
  #write <quark_parent/P0_`split'.cpp> "  P0`split' = %e", P

  #write <quark_parent/P0_`split'.cpp> "  std::vector<std::vector<std::complex<double>>> output = {{P0`split', 0},{0,P0`split'}};"
  #write <quark_parent/P0_`split'.cpp> "  return output;
  }"
  #clearoptimize
  .sort
  drop P;
  .sort
#else if (("`split'"=="qbp_qp_qb_q") || ("`split'"=="qb_q_qb_q")  || ("`split'"=="qb_g_g_q") || ("`split'"=="symm_g_g_g_g"))
  #include gluon_parent/P_`split'.m
  id D = 4;
  .sort

  Format O1;
  .sort
  Format C;
  #write <gluon_parent/P0_`split'.cpp> "#include \"../splitting_functions.hpp\""
  #write <gluon_parent/P0_`split'.cpp> "std::vector<std::vector<std::complex<double>>> P0_`split'(LV<double> p1, LV<double> p2, LV<double> p3, LV<double> p4) {"
  #write <gluon_parent/P0_`split'.cpp> "  double CA = C_A;"
  #write <gluon_parent/P0_`split'.cpp> "  double CF = C_F;"
  #write <gluon_parent/P0_`split'.cpp> "  double s12 = 2.*p1*p2;"
  #write <gluon_parent/P0_`split'.cpp> "  double s13 = 2.*p1*p3;"
  #write <gluon_parent/P0_`split'.cpp> "  double s14 = 2.*p1*p4;"
  #write <gluon_parent/P0_`split'.cpp> "  double s23 = 2.*p2*p3;"
  #write <gluon_parent/P0_`split'.cpp> "  double s24 = 2.*p2*p4;"
  #write <gluon_parent/P0_`split'.cpp> "  double s34 = 2.*p3*p4;"
  #write <gluon_parent/P0_`split'.cpp> "  double s123 = s12 + s13 + s23;"
  #write <gluon_parent/P0_`split'.cpp> "  double s134 = s13 + s14 + s34;"
  #write <gluon_parent/P0_`split'.cpp> "  double s124 = s12 + s14 + s24;"
  #write <gluon_parent/P0_`split'.cpp> "  double s234 = s23 + s24 + s34;"
  #write <gluon_parent/P0_`split'.cpp> "  double s1234 = s12 + s13 + s14 + s23 + s24 + s34;"
  #write <gluon_parent/P0_`split'.cpp> "  LV<double> p = p1 + p2 + p3 + p4;"
  #do i=1,4
    #write <gluon_parent/P0_`split'.cpp> "  double z`i' = p`i'.components[0]/p.components[0];"
    #write <gluon_parent/P0_`split'.cpp> "  LV<double> k`i' = p`i' - z`i'*p;"
  #enddo
  #write <gluon_parent/P0_`split'.cpp> "  double z12 = z1 + z2;"
  #write <gluon_parent/P0_`split'.cpp> "  double z14 = z1 + z4;"
  #write <gluon_parent/P0_`split'.cpp> "  double z23 = z2 + z3;"
  #write <gluon_parent/P0_`split'.cpp> "  double z13 = z1 + z3;"
  #write <gluon_parent/P0_`split'.cpp> "  double z34 = z3 + z4;"
  #write <gluon_parent/P0_`split'.cpp> "  double z24 = z2 + z4;"
  #write <gluon_parent/P0_`split'.cpp> "  double z123 = z1 + z2 + z3;"
  #write <gluon_parent/P0_`split'.cpp> "  double z134 = z1 + z3 + z4;"
  #write <quark_parent/P0_`split'.cpp> "  double z124 = z1 + z2 + z4;"
  #write <quark_parent/P0_`split'.cpp> "  double z234 = z2 + z3 + z4;"

  #write <gluon_parent/P0_`split'.cpp> "  std::complex<double> P0`split' = 0.;"

  #write <gluon_parent/P0_`split'.cpp> "  std::vector<std::vector<std::complex<double>>> output = {{0., 0.}, {0., 0.}};"
  #write <gluon_parent/P0_`split'.cpp> "  for(int i = 0; i <= 1; i++) for(int j = 0; j <= 1; j++) {"
  #write <gluon_parent/P0_`split'.cpp> "    LV<std::complex<double>> E0 = polarization(p, (i?1:-1));"
  #write <gluon_parent/P0_`split'.cpp> "    LV<std::complex<double>> E0conj = polarization(p, (j?1:-1)).conj();"

  #do i=1,4
    #do j=1,4
      #write <gluon_parent/P0_`split'.cpp> "  std::complex<double> k`i'k`j' = (k`i'*E0)*(k`j'*E0conj);"
    #enddo
  #enddo
  #write <gluon_parent/P0_`split'.cpp> "  std::complex<double> gg = E0*E0conj;"

  ExtraSymbols array,Z;
  Format C;
  #optimize P
  #write <gluon_parent/P0_`split'.cpp> "  std::complex<double> Z[`optimmaxvar_'+1];"
  #write <gluon_parent/P0_`split'.cpp> "%O"
  #write <gluon_parent/P0_`split'.cpp> "  P0`split' = %e", P

  #write <gluon_parent/P0_`split'.cpp> "  output[i][j] = P0`split';"
  #write <gluon_parent/P0_`split'.cpp> "  };"
  #write <gluon_parent/P0_`split'.cpp> "  return output;
  }"
  #clearoptimize
  .sort
  drop P;
  .sort
#else
  #write "Unknown splitting function: `split'"
#endif


.end