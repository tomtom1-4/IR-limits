#include "main.hpp"

// Evaluation of input
std::string process_str = "d d~ -> d d~";
std::string unresolved_str = " g";
std::string suffix = "";


const int nBorn = 4;
const int power = 2;
const std::array<unsigned, 3> powerNonQCD = {0,0,0};
const std::string order = "NLO";

void replace(std::string& str, const std::string& from, const std::string& to) {
  if(from.empty())
    return;
  size_t start_pos = 0;
  while((start_pos = str.find(from, start_pos)) != std::string::npos) {
    str.replace(start_pos, from.length(), to);
    start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
  }
}

int main() {
  std::string order_rcl = "LO";
  unsigned loopOrder = 0;
  if(order == "NLO") loopOrder = 1;
  else if(order == "NNLO") loopOrder = 2;
  int delta_power = 0;
  if ((order == "NLO") or (order == "NNLO")) {
    delta_power = 2;
    order_rcl = "NLO";
  }
  std::cout.precision(8);
  // Initialize txt outstream (delete previous file)

  std::ofstream outfile1;
  outfile1.open ("results/" + process_str + " +" + unresolved_str + "_" + order + ".txt");
  outfile1.close();
  outfile1.open ("results/" + process_str + " +" + unresolved_str + "_" + order + ".txt");

  amplitude A(process_str);

  // Stripper Settings
  Stripper::Model::config((powerNonQCD[0]!=0),false,false,false,powerNonQCD[0],powerNonQCD[1],powerNonQCD[2],1.,1.,1.);
  const Stripper::Process process(process_str);
  Stripper::Model::nf = 5;

  // Recola Settings
  Recola::set_reduction_mode_rcl(4);
  Recola::set_print_level_amplitude_rcl(2);
  Recola::set_alphas_rcl(gs * gs/(4 * M_PI), mu, 5);
  Recola::set_mu_ir_rcl(mu);
  Recola::set_delta_ir_rcl(0., M_PI*M_PI/12.);
  Recola::set_momenta_correction_rcl(false);

  Recola::use_alpha0_scheme_rcl(e*e/4./M_PI);
  Recola::set_pole_mass_z_rcl(1.e8, 0.0001);
  Recola::set_pole_mass_w_rcl(1.e15, 0.0001);

  // Define & generate process
  Recola::define_process_rcl(1, process_str, order_rcl);
  Recola::set_otter_mode_rcl(1, "oneloop_qp");

  Recola::generate_processes_rcl();
  // Define initial state
  PSF::PhaseSpace dumn = PSF::Splitting(nBorn - 2, COM);
  PSF::PhaseSpace dumn2 = PSF::Splitting(nBorn - 2, COM);
  PSF::PhaseSpace dumn3 = PSF::Splitting(nBorn - 2, COM);
  //PSF::PhaseSpace dumn4 = PSF::Splitting(nBorn - 2, COM);
  PSF::PhaseSpace pp = PSF::Splitting(nBorn - 2, COM);
  pp.print();
  // Transform to Recola format
  double pp_rcl[nBorn][4], pp_arr[4*nBorn];
  for(int i = 0; i < nBorn; i++) {
    for(int j = 0; j < 4; j++) {
      pp_rcl[i][j] = (i<2?-1.:1.)*pp.momenta[i].components[j];
      pp_arr[i*4+j] = (i<2?-1.:1.)*pp.momenta[i].components[j];
    }
  }
  std::vector<Stripper::Momentum<double>> pp_Stripper;
  for(int i = 0; i < nBorn; i++) {
    Stripper::Momentum<double> p_Stripper = {pp_rcl[i][0], pp_rcl[i][1], pp_rcl[i][2], pp_rcl[i][3],0.,0.};
    pp_Stripper.push_back(p_Stripper);
  }

  XMLDocument xmlDoc;
  std::string process_str_copy = process_str;
  replace(process_str_copy, "-> ", "");
  replace(process_str_copy, "~", "");
  replace(process_str_copy, " " , "");
  replace(process_str_copy, "-" , "");
  replace(process_str_copy, "+" , "");
  XMLNode * pRoot = xmlDoc.NewElement(const_cast<char*>(process_str_copy.c_str()));
  xmlDoc.InsertFirstChild(pRoot);
  XMLElement * ppElement = xmlDoc.NewElement("PhasePhacePoint");
  for(int i = 0; i < A.process.size(); i++) {
    for(int mu = 0; mu < 4; mu++) {
      XMLElement * ppEntry = xmlDoc.NewElement(const_cast<char*>(("p" + std::to_string(i) + "_" + std::to_string(mu)).c_str()));
      ppEntry->SetText(pp.momenta[i].components[mu]);
      ppElement->InsertEndChild(ppEntry);
    }
  }
  pRoot->InsertEndChild(ppElement);

  // Compute process amplitudes
  Recola::compute_process_rcl(1, pp_rcl, order_rcl);

  // Get non-vanishing helicity and color-configurations
  std::vector<std::vector<int> > helicities = non0hel(1, power + delta_power, order_rcl, A);
  std::cout << "Found non zero helicity configurations of Born process" << std::endl;
  std::vector<std::vector<int> > colors = non0col(1, power + delta_power, order_rcl, A, helicities);
  std::cout << "Found non zero color configurations of Born process" << std::endl;

  // Create Hashmap for all non-zero amplitudes that only need to be caluculated once per data set
  std::unordered_map<std::string, std::complex<double>> M0, M1; // hel + col is the key as string
  std::unordered_map<std::string, std::complex<double>> M0_ij, M1_ij, M1_ij_Im, M2_ij, fM0_ijk, fM1_ijk, dM0_ijk, dM_ijk, M0_ijkl, M1_ijkl, M0_ijklab, Q_ijkl, M0_ijkla; // color correlators
  std::vector<std::vector<int>> keys;
  // fill the Hashmaps
  for (int i = 0; i < helicities.size(); i++) {
    // Create hel variables
    std::vector<int> hel_full = helicities[i];
    int hel[A.process.size()];
    std::string hel_string;
    for (int j = 0; j < hel_full.size(); j++) {
      hel[j] = hel_full[j];
      hel_string += std::to_string(hel_full[j]);
    }

    for (int j = 0; j < colors.size(); j++) {
      // Create col variables
      std::vector<int> col_full = colors[j];
      int col[A.process.size()];
      std::string col_string;
      for (int k = 0; k < col_full.size(); k++) {
        col[k] = col_full[k];
        col_string += std::to_string(col_full[k]);
      }
      std::string key = hel_string + col_string;
      // Compute amplitudes
      // Born amplitude
      M0[key] = colorflow2color(hel, col, A, power, "LO", 1);
      if(order_rcl=="NLO") M1[key] = colorflow2color(hel, col, A, power + 2, "NLO", 1);
      if(std::abs(M0[key]) + std::abs(M1[key]) > 1.e-17) {
        std::vector<int> dummy = hel_full;
        dummy.insert(dummy.end(), col_full.begin(), col_full.end());
        keys.push_back(dummy);
        std::cout << key << "\t" << M0[key];
        if(order_rcl=="NLO") std::cout << "\t" << M1[key];
        std::cout << std::endl;
      }
    }
  }
  double average_factor = 1./4.; // Spin of the initial state
  double average_factor_full = 1./4.;
  {
  for(int i = 0; i < 2; i++) {
    average_factor *= 1./double(A.process[i]);
  }
  std::vector<std::string> particles = A.process_particles;
  for(auto s : particles) std::cout << s << ", " ;
  std::cout << std::endl;
  particles.erase(particles.begin(), particles.begin() + 2);
  for(auto s : particles) std::cout << s << ", " ;
  std::cout << std::endl;
  while(particles.size() > 0) {
    std::string particle = particles[0];
    int occurances = 1;
    int i = 1;
    while(i < particles.size()) {
      if(particles[i] == particle) {
        occurances += 1;
        particles.erase(particles.begin() + i);
      }
      else {
        i += 1;
      }
    }
    particles.erase(particles.begin());
    average_factor *= 1./factorial(occurances);
  }
  std::cout << average_factor << std::endl;
  }
  std::cout << "A.particle_type = ";
  for(int i : A.particle_type) {
    std::cout << i << ",";
  }
  std::cout << std::endl;
  std::cout << "A.process = ";
  for(int i : A.process) {
    std::cout << i << ",";
  }
  std::cout << std::endl;

  std::shared_ptr<Stripper::LoopImplementation<double>> matrixelements;
  if(order=="NNLO") {
    matrixelements = std::make_shared<Stripper::Hardcoded<double>>(loopOrder,process,powerNonQCD);
    matrixelements->setKinematics(pp_Stripper);
  }

  // Squared Matrixelements
  double M2_averaged_NNLO;
  double M2_averaged_NLO;
  double M2_averaged_LO;

  Recola::get_squared_amplitude_rcl(1, power, "LO", M2_averaged_LO);
  Recola::get_squared_amplitude_rcl(1, power + delta_power/2, "NLO", M2_averaged_NLO);
  if(order == "NNLO") M2_averaged_NNLO = (matrixelements->eval(mu))*average_factor*std::pow(gs, 2*power + 4);
  std::cout << "M2 (LO) = " << M2_averaged_LO << std::endl;
  std::cout << "M2 (NLO) = " << M2_averaged_NLO << std::endl;
  std::cout << "M2 (NNLO) = " << M2_averaged_NNLO << std::endl;

  XMLElement* M2_Element_LO = xmlDoc.NewElement("M0");
  M2_Element_LO->SetText(M2_averaged_LO);
  pRoot->InsertEndChild(M2_Element_LO);

  XMLElement* M2_Element_NLO = xmlDoc.NewElement("M1");
  M2_Element_NLO->SetText(M2_averaged_NLO);
  pRoot->InsertEndChild(M2_Element_NLO);

  if(order == "NNLO") {
    XMLElement* M2_Element_NNLO = xmlDoc.NewElement("M2");
    M2_Element_NNLO->SetText(M2_averaged_NNLO);
    pRoot->InsertEndChild(M2_Element_NNLO);
  }

  // Get spin correlators
  std::unordered_map<std::string, std::complex<double>> SC0, SC1;
  for(int i = 0; i < A.process.size(); i++) {
    if(A.particle_type[i] == 0) continue;
    for(std::vector<int> hel : helicities) {  // color and helicity of the bra <M|
      std::string hel_string;
      for (int dummy = 0; dummy < hel.size(); dummy++) {
        hel_string += std::to_string(hel[dummy]);
      }
      for(std::vector<int> col : colors) {
        std::string col_string;
        for(int dummy = 0; dummy < col.size(); dummy++) {
          col_string += std::to_string(col[dummy]);
        }
        std::string key = hel_string + col_string;
        std::complex<double> M0_bra = std::conj(M0[key]);
        std::complex<double> M1_bra = std::conj(M1[key]);

        std::vector<int> hel2 = hel;
        for(int hel_i = -1; hel_i <= 1; hel_i += 2) {
          hel2[i] = hel_i;
          std::string hel2_string;
          for(int dummy = 0; dummy < hel.size(); dummy++) {
            hel2_string += std::to_string(hel2[dummy]);
          }
          std::complex<double> M0_ket = M0[hel2_string + col_string];
          std::complex<double> M1_ket = M1[hel2_string + col_string];
          SC0[std::to_string(i) + std::to_string(hel[i]) + std::to_string(hel2[i])] += (M0_bra*M0_ket)*average_factor;
          SC1[std::to_string(i) + std::to_string(hel[i]) + std::to_string(hel2[i])] += (M0_bra*M1_ket + M1_bra*M0_ket)*average_factor;

        }
      }
    }
  }

  // Get Spin corrletors
  XMLElement* SC0_Element = xmlDoc.NewElement("SC0");
  XMLElement* SC1_Element = xmlDoc.NewElement("SC1");
  for (int i = 0; i < A.process.size(); i++) {
    if(A.particle_type[i] == 0) continue;
    for(int s1 = -1; s1 <= 1; s1 += 2) {
      for(int s2 = -1; s2 <= 1; s2 += 2) {
        XMLElement* SC0_iEntry_real = xmlDoc.NewElement(const_cast<char*>(("r" + std::to_string(i) + std::to_string(s1) + std::to_string(s2)).c_str()));
        XMLElement* SC1_iEntry_real = xmlDoc.NewElement(const_cast<char*>(("r" + std::to_string(i) + std::to_string(s1) + std::to_string(s2)).c_str()));
        XMLElement* SC0_iEntry_imag = xmlDoc.NewElement(const_cast<char*>(("i" + std::to_string(i) + std::to_string(s1) + std::to_string(s2)).c_str()));
        XMLElement* SC1_iEntry_imag = xmlDoc.NewElement(const_cast<char*>(("i" + std::to_string(i) + std::to_string(s1) + std::to_string(s2)).c_str()));
        SC0_iEntry_real->SetText(std::real(SC0[std::to_string(i) + std::to_string(s1) + std::to_string(s2)]));
        SC1_iEntry_real->SetText(std::real(SC1[std::to_string(i) + std::to_string(s1) + std::to_string(s2)]));
        SC0_iEntry_imag->SetText(std::imag(SC0[std::to_string(i) + std::to_string(s1) + std::to_string(s2)]));
        SC1_iEntry_imag->SetText(std::imag(SC1[std::to_string(i) + std::to_string(s1) + std::to_string(s2)]));
        std::cout << i << ", " << s1 << ", " << s2 << "\t" << SC0[std::to_string(i) + std::to_string(s1) + std::to_string(s2)] << "\t" <<SC1[std::to_string(i) + std::to_string(s1) + std::to_string(s2)] << std::endl;
        SC0_Element->InsertEndChild(SC0_iEntry_real);
        SC0_Element->InsertEndChild(SC0_iEntry_imag);
        SC1_Element->InsertEndChild(SC1_iEntry_real);
        SC1_Element->InsertEndChild(SC1_iEntry_imag);
      }
    }
  }
  pRoot->InsertEndChild(SC0_Element);
  pRoot->InsertEndChild(SC1_Element);

  // Get color correlators
  // <M|Ti.Tj|M>
  XMLElement* M2_ijElement = xmlDoc.NewElement("M2_ij");
  for (int i = 0; i < A.process.size(); i++) {
    if(A.particle_type[i] == 0) continue;
    for(int j = 0; j < A.process.size(); j++) {
      if(A.particle_type[j] == 0) continue;
      XMLElement* M2_ijEntry = xmlDoc.NewElement(const_cast<char*>(("d" + std::to_string(i) + std::to_string(j)).c_str()));
      double M2ij = 0;
      if(order == "NNLO") {
        M2ij = (matrixelements->evalCC(mu,i,j))*std::pow(gs, 2*power + 4);
        M2_ij[std::to_string(i) + std::to_string(j)] = M2ij*average_factor;
      }

      M2_ijEntry->SetText(std::real(M2_ij[std::to_string(i) + std::to_string(j)]));
      M2_ijElement->InsertEndChild(M2_ijEntry);

      if(order == "NNLO") {
        std::cout << "<M1|T_" << i << ".T_" << j << "|M1> + (<M0|T_" << i << ".T_" << j << "|M2> + c.c.) = " << M2ij*average_factor << "\t" << M2ij*average_factor/M2_averaged_NNLO << std::endl;
      }
    }
  }
  pRoot->InsertEndChild(M2_ijElement);
  std::cout << std::endl;

  XMLElement* M0_ijElement = xmlDoc.NewElement("M0_ij");
  XMLElement* M1_ijElement = xmlDoc.NewElement("M1_ij");
  //XMLElement* M1_ij_ImElement = xmlDoc.NewElement("M1_ij_Im");

  for (int i = 0; i < A.process.size(); i++) {
    if(A.particle_type[i] == 0) continue;
    for(int j = 0; j < A.process.size(); j++) {
      if(A.particle_type[j] == 0) continue;
      XMLElement* M0_ijEntry = xmlDoc.NewElement(const_cast<char*>(("d" + std::to_string(i) + std::to_string(j)).c_str()));
      XMLElement* M1_ijEntry = xmlDoc.NewElement(const_cast<char*>(("d" + std::to_string(i) + std::to_string(j)).c_str()));
      //XMLElement* M1_ij_ImEntry = xmlDoc.NewElement(const_cast<char*>(("d" + std::to_string(i) + std::to_string(j)).c_str()));
      if(std::abs(M0_ij[std::to_string(j) + std::to_string(i)]) + std::abs(M1_ij[std::to_string(j) + std::to_string(i)]) != 0.) {
        M0_ij[std::to_string(i) + std::to_string(j)] = M0_ij[std::to_string(j) + std::to_string(i)];
        M1_ij[std::to_string(i) + std::to_string(j)] = M1_ij[std::to_string(j) + std::to_string(i)];
        //M1_ij_Im[std::to_string(i) + std::to_string(j)] = M1_ij_Im[std::to_string(j) + std::to_string(i)];
        M0_ijEntry->SetText(std::real(M0_ij[std::to_string(i) + std::to_string(j)]));
        M1_ijEntry->SetText(std::real(M1_ij[std::to_string(i) + std::to_string(j)]));
        //M1_ij_ImEntry->SetText(std::real(M1_ij_Im[std::to_string(i) + std::to_string(j)]));
        M0_ijElement->InsertEndChild(M0_ijEntry);
        M1_ijElement->InsertEndChild(M1_ijEntry);
        //M1_ij_ImElement->InsertEndChild(M1_ij_ImEntry);
        continue;
      }

      //Recola::get_squared_amplitude_rcl(1, power + delta_power/2, order, M2_averaged_control);
      double M0ij = 0;
      double M1ij = 0;
      double M1ij_Im = 0;

      for(std::vector<int> helcol : keys) { // color and helicity of the bra <M|
        std::vector<int> hel_full(helcol.begin(), helcol.begin() + A.process.size());
        std::vector<int> col_full(helcol.begin() + A.process.size(), helcol.end());
        std::string hel_string;
        for (int dummy = 0; dummy < hel_full.size(); dummy++) {
          hel_string += std::to_string(hel_full[dummy]);
        }
        std::string col_string;
        for(int dummy = 0; dummy < col_full.size(); dummy++) {
          col_string += std::to_string(col_full[dummy]);
        }
        std::string key = hel_string + col_string;
        std::complex<double> M0_bra = std::conj(M0[key]);
        std::complex<double> M1_bra = std::conj(M1[key]);
        for(int a = 0; a < 8; a++) {
          for(int ci = 0; ci < A.process[i]; ci++) {
            for(int cj = 0; cj < A.process[j]; cj++) {
              int j_color;
              if(i == j) {
                key.replace(key.size() - A.process.size() + j, 1, std::to_string(cj));
                j_color = ci;
              }
              else {
                key.replace(key.size() - A.process.size() + i, 1, std::to_string(ci));
                key.replace(key.size() - A.process.size() + j, 1, std::to_string(cj));
                j_color = col_full[j];
              }
              std::complex<double> col_factor = 1;
              if(A.particle_type[j] == 1)
                col_factor *= lam[a][3*j_color + cj]/2.;
              else if(A.particle_type[j] == -1)
                col_factor *= -lam[a][3*cj + j_color]/2.;
              else if(A.particle_type[j] == 2)
                col_factor *= I*fabc[a][8*cj + j_color];

              if(A.particle_type[i] == 1)
                col_factor *= lam[a][3*col_full[i] + ci]/2.;
              else if(A.particle_type[i] == -1)
                col_factor *= -lam[a][3*ci + col_full[i]]/2.;
              else if(A.particle_type[i] == 2)
                col_factor *= I*fabc[a][8*ci + col_full[i]];
              M0ij += std::real(col_factor*M0_bra*M0[key]);
              M1ij += 2*std::real(col_factor*M1_bra*M0[key]);
              //M1ij_Im += -2*std::imag(col_factor*M1_bra*M0[key]);
            }
          }
        }
      }
      M1_ij[std::to_string(i) + std::to_string(j)] = M1ij*average_factor;
      //M1_ij_Im[std::to_string(i) + std::to_string(j)] = M1ij_Im*average_factor;
      M0_ij[std::to_string(i) + std::to_string(j)] = M0ij*average_factor;
      M0_ijEntry->SetText(std::real(M0_ij[std::to_string(i) + std::to_string(j)]));
      M1_ijEntry->SetText(std::real(M1_ij[std::to_string(i) + std::to_string(j)]));
      //M1_ij_ImEntry->SetText(std::real(M1_ij_Im[std::to_string(i) + std::to_string(j)]));

      M0_ijElement->InsertEndChild(M0_ijEntry);
      M1_ijElement->InsertEndChild(M1_ijEntry);
      //M1_ij_ImElement->InsertEndChild(M1_ij_ImEntry);
      std::cout << "<M0|T_" << i << ".T_" << j << "|M0> = " << M0ij*average_factor << "\t" << M0ij*average_factor/M2_averaged_LO << std::endl;
      if(order_rcl == "NLO") std::cout << "<M0|T_" << i << ".T_" << j << "|M1> + c.c. = " << M1ij*average_factor << "\t" << M1ij*average_factor/M2_averaged_NLO << std::endl;
      //if(order == "NNLO") {
      //  std::cout << "<M0|T_" << i << ".T_" << j << "|M1> - c.c. = " << M1ij_Im*average_factor << "\t" << M1ij_Im*average_factor/M2_averaged_NLO << std::endl;
      //}
      std::cout << std::endl;
    }
  }
  pRoot->InsertEndChild(M0_ijElement);
  pRoot->InsertEndChild(M1_ijElement);
  //pRoot->InsertEndChild(M1_ij_ImElement);

  // <M|Ti.Tj Tk.Tl|M>

  if(unresolved_str==" g g" || unresolved_str==" g g g") {
    XMLElement * M0_ijklElement = xmlDoc.NewElement("M0_ijkl");
    XMLElement * M1_ijklElement = xmlDoc.NewElement("M1_ijkl");
    for (int i = 0; i < A.process.size(); i++) {
      if(A.particle_type[i] == 0) continue;
      for(int j = 0; j < A.process.size(); j++) {
        if(A.particle_type[j] == 0) continue;
        for(int k = 0; k < A.process.size(); k++) {
          if(A.particle_type[k] == 0) continue;
          for(int l = 0; l < A.process.size(); l++) {
            if(A.particle_type[l] == 0) continue;
            std::complex<double> M0ijkl = 0.;
            std::complex<double> M1ijkl = 0.;
            XMLElement * M0_ijklEntry = xmlDoc.NewElement(const_cast<char*>(("d" + std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l)).c_str()));
            XMLElement * M1_ijklEntry = xmlDoc.NewElement(const_cast<char*>(("d" + std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l)).c_str()));
            for(std::vector<int> hel_full : helicities) {  // color and helicity of the bra <M|
              std::string hel_string;
              for (int dummy = 0; dummy < hel_full.size(); dummy++) {
                hel_string += std::to_string(hel_full[dummy]);
              }
              for(std::vector<int> col_full : colors) {
                std::string col_string;
                for(int dummy = 0; dummy < col_full.size(); dummy++) {
                  col_string += std::to_string(col_full[dummy]);
                }
                std::string key = hel_string + col_string;
                std::complex<double> M0_bra = std::conj(M0[key]);
                std::complex<double> M1_bra = std::conj(M1[key]);
                std::vector<int> output_colors = {col_full[i], col_full[j], col_full[k], col_full[l]};
                for(int a = 0; a < 8; a++) for(int b = 0; b < 8; b++) {
                  for(int ci = 0; ci < A.process[i]; ci++) for(int cj = 0; cj < A.process[j]; cj++) for(int ck = 0; ck < A.process[k]; ck++) for(int cl = 0; cl < A.process[l]; cl++) {
                    key.replace(key.size() - A.process.size() + i, 1, std::to_string(ci));
                    key.replace(key.size() - A.process.size() + j, 1, std::to_string(cj));
                    key.replace(key.size() - A.process.size() + k, 1, std::to_string(ck));
                    key.replace(key.size() - A.process.size() + l, 1, std::to_string(cl));

                    if(j == i)
                      output_colors[1] = ci;
                    else if(k == i)
                      output_colors[2] = ci;
                    else if(l == i)
                      output_colors[3] = ci;

                    if(k == j)
                      output_colors[2] = cj;
                    else if(l == j)
                      output_colors[3] = cj;

                    if(l == k)
                      output_colors[3] = ck;

                    std::complex<double> col_factor = 1;
                    if(A.particle_type[l] == 1)
                      col_factor *= lam[a][3*output_colors[3] + cl]/2.;
                    else if(A.particle_type[l] == -1)
                      col_factor *= -lam[a][3*cl + output_colors[3]]/2.;
                    else if(A.particle_type[l] == 2)
                      col_factor *= I*fabc[a][8*cl + output_colors[3]];

                    if(A.particle_type[k] == 1)
                      col_factor *= lam[a][3*output_colors[2] + ck]/2.;
                    else if(A.particle_type[k] == -1)
                      col_factor *= -lam[a][3*ck + output_colors[2]]/2.;
                    else if(A.particle_type[k] == 2)
                      col_factor *= I*fabc[a][8*ck + output_colors[2]];

                    if(A.particle_type[j] == 1)
                      col_factor *= lam[b][3*output_colors[1] + cj]/2.;
                    else if(A.particle_type[j] == -1)
                      col_factor *= -lam[b][3*cj + output_colors[1]]/2.;
                    else if(A.particle_type[j] == 2)
                      col_factor *= I*fabc[b][8*cj + output_colors[1]];

                    if(A.particle_type[i] == 1)
                      col_factor *= lam[b][3*output_colors[0] + ci]/2.;
                    else if(A.particle_type[i] == -1)
                      col_factor *= -lam[b][3*ci + output_colors[0]]/2.;
                    else if(A.particle_type[i] == 2)
                      col_factor *= I*fabc[b][8*ci + output_colors[0]];

                    M0ijkl += col_factor*M0_bra*M0[key];
                    M1ijkl += 2.*std::real(col_factor*M0_bra*M1[key]);
                  }
                }
              }
            }

            M0_ijkl[std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l)] = M0ijkl*average_factor;
            M1_ijkl[std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l)] = M1ijkl*average_factor;
            M0_ijklEntry->SetText(std::real(M0_ijkl[std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l)]));
            M1_ijklEntry->SetText(std::real(M1_ijkl[std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l)]));
            M0_ijklElement->InsertEndChild(M0_ijklEntry);
            M1_ijklElement->InsertEndChild(M1_ijklEntry);
            double M2_averaged;
            Recola::get_squared_amplitude_rcl(1, power, "LO", M2_averaged);
            std::cout << "M2 (LO) = " << M2_averaged << std::endl;
            std::cout << "i = " << i << ", j = " << j << ", k = " << k << ", l = " << l << ": <M|T_i.T_j T_k.T_l|M> = " << M0ijkl*average_factor << "\t" << M0ijkl*average_factor/M2_averaged << std::endl;
          }
        }
      }
    }
    pRoot->InsertEndChild(M0_ijklElement);
    pRoot->InsertEndChild(M1_ijklElement);
  }



  // <M|Ti.Tj Tk.Tl Ta.Tb|M>
  if(unresolved_str==" g g g") {
  for (int i = 0; i < A.process.size(); i++) {
    if(A.particle_type[i] == 0) continue;
    for(int j = 0; j < A.process.size(); j++) {
      if(A.particle_type[j] == 0) continue;
      for(int k = 0; k < A.process.size(); k++) {
        if(A.particle_type[k] == 0) continue;
        for(int l = 0; l < A.process.size(); l++) {
          if(A.particle_type[l] == 0) continue;
          for(int a = 0; a < A.process.size(); a++) {
            if(A.particle_type[a] == 0) continue;
            for(int b = 0; b < A.process.size(); b++) {
              if(A.particle_type[b] == 0) continue;
              std::complex<double> Mijklab = 0.;
              for(std::vector<int> hel_full : helicities) {  // color and helicity of the bra <M|
                std::string hel_string;
                for (int dummy = 0; dummy < hel_full.size(); dummy++) {
                  hel_string += std::to_string(hel_full[dummy]);
                }
                for(std::vector<int> col_full : colors) {
                  std::string col_string;
                  for(int dummy = 0; dummy < col_full.size(); dummy++) {
                    col_string += std::to_string(col_full[dummy]);
                  }
                  std::string key = hel_string + col_string;
                  std::complex<double> M_bra = std::conj(M0[key]);
                  std::vector<int> output_colors = {col_full[i], col_full[j], col_full[k], col_full[l], col_full[a], col_full[b]};
                  for(int m = 0; m < 8; m++) for(int n = 0; n < 8; n++) for(int c = 0; c < 8; c++) {
                    for(int ci = 0; ci < A.process[i]; ci++) for(int cj = 0; cj < A.process[j]; cj++)
                    for(int ck = 0; ck < A.process[k]; ck++) for(int cl = 0; cl < A.process[l]; cl++)
                    for(int ca = 0; ca < A.process[a]; ca++) for(int cb = 0; cb < A.process[b]; cb++) {
                      key.replace(key.size() - A.process.size() + i, 1, std::to_string(ci));
                      key.replace(key.size() - A.process.size() + j, 1, std::to_string(cj));
                      key.replace(key.size() - A.process.size() + k, 1, std::to_string(ck));
                      key.replace(key.size() - A.process.size() + l, 1, std::to_string(cl));
                      key.replace(key.size() - A.process.size() + a, 1, std::to_string(ca));
                      key.replace(key.size() - A.process.size() + b, 1, std::to_string(cb));

                      if(j == i)
                        output_colors[1] = ci;
                      else if(k == i)
                        output_colors[2] = ci;
                      else if(l == i)
                        output_colors[3] = ci;
                      else if(a == i)
                        output_colors[4] = ci;
                      else if(b == i)
                        output_colors[5] = ci;

                      if(k == j)
                        output_colors[2] = cj;
                      else if(l == j)
                        output_colors[3] = cj;
                      else if(a == j)
                        output_colors[4] = cj;
                      else if(b == j)
                        output_colors[5] = cj;

                      if(l == k)
                        output_colors[3] = ck;
                      else if(a == k)
                        output_colors[4] = ck;
                      else if(b == k)
                        output_colors[5] = ck;

                      if(a == l)
                        output_colors[4] = cl;
                      else if(b == l)
                        output_colors[5] = cl;

                      if(b == a)
                        output_colors[5] = ca;

                      std::complex<double> col_factor = 1;
                      if(A.particle_type[b] == 1)
                        col_factor *= lam[c][3*output_colors[5] + cb]/2.;
                      else if(A.particle_type[b] == -1)
                        col_factor *= -lam[c][3*cb + output_colors[5]]/2.;
                      else if(A.particle_type[l] == 2)
                        col_factor *= I*fabc[c][8*cb + output_colors[5]];

                      if(A.particle_type[a] == 1)
                        col_factor *= lam[c][3*output_colors[4] + ca]/2.;
                      else if(A.particle_type[l] == -1)
                        col_factor *= -lam[c][3*ca + output_colors[4]]/2.;
                      else if(A.particle_type[l] == 2)
                        col_factor *= I*fabc[c][8*ca + output_colors[4]];

                      if(A.particle_type[l] == 1)
                        col_factor *= lam[n][3*output_colors[3] + cl]/2.;
                      else if(A.particle_type[l] == -1)
                        col_factor *= -lam[n][3*cl + output_colors[3]]/2.;
                      else if(A.particle_type[l] == 2)
                        col_factor *= I*fabc[n][8*cl + output_colors[3]];

                      if(A.particle_type[k] == 1)
                        col_factor *= lam[n][3*output_colors[2] + ck]/2.;
                      else if(A.particle_type[k] == -1)
                        col_factor *= -lam[n][3*ck + output_colors[2]]/2.;
                      else if(A.particle_type[k] == 2)
                        col_factor *= I*fabc[n][8*ck + output_colors[2]];

                      if(A.particle_type[j] == 1)
                        col_factor *= lam[m][3*output_colors[1] + cj]/2.;
                      else if(A.particle_type[j] == -1)
                        col_factor *= -lam[m][3*cj + output_colors[1]]/2.;
                      else if(A.particle_type[j] == 2)
                        col_factor *= I*fabc[m][8*cj + output_colors[1]];

                      if(A.particle_type[i] == 1)
                        col_factor *= lam[m][3*output_colors[0] + ci]/2.;
                      else if(A.particle_type[i] == -1)
                        col_factor *= -lam[m][3*ci + output_colors[0]]/2.;
                      else if(A.particle_type[i] == 2)
                        col_factor *= I*fabc[m][8*ci + output_colors[0]];

                      Mijklab += col_factor*M_bra*M0[key];
                    }
                  }
                }
              }

              M0_ijklab[std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l) + std::to_string(a) + std::to_string(b)] = Mijklab*average_factor;
              double M2_averaged;
              Recola::get_squared_amplitude_rcl(1, power, "LO", M2_averaged);
              std::cout << "i = " << i << ", j = " << j << ", k = " << k << ", l = " << l << ", a = " << a << ", b = " << b << ": <M|T_i.T_j T_k.T_l Ta.Tb|M> = " << Mijklab*average_factor << "\t" << Mijklab*average_factor/M2_averaged << std::endl;
            }
          }
        }
      }
    }
  }
  }

  // <M|Ti.Tj f^{ck, cl, ca} Tk^ck Tl^cl Ta^ca|M> + c.c.
  if((unresolved_str==" g g g") or (unresolved_str == " g g" and order == "NLO") and suffix=="_EW") {
  XMLElement * M0_ijklaElement = xmlDoc.NewElement("M0_ijkla");
  std::vector<std::vector<int>> configurations;
  for (int i = 0; i < A.process.size(); i++) {
    if(A.particle_type[i] == 0) continue;
    for(int j = 0; j < A.process.size(); j++) {
      if(A.particle_type[j] == 0) continue;
      for(int k = 0; k < A.process.size(); k++) {
        if(A.particle_type[k] == 0) continue;
        for(int l = 0; l < A.process.size(); l++) {
          if(A.particle_type[l] == 0) continue;
          if(l == k) continue;
          for(int a = 0; a < A.process.size(); a++) {
            if(A.particle_type[a] == 0) continue;
            if((a == l) or (a == k)) continue;
            XMLElement * M0_ijklaEntry = xmlDoc.NewElement(const_cast<char*>(("d" + std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l) + std::to_string(a)).c_str()));
            bool permutation = false;
            for(std::vector<int> configuration : configurations) {
              if(configuration == std::vector<int>({i,j,k,l,a})) {
                M0_ijklaEntry->SetText(std::real(M0_ijkla[std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l) + std::to_string(a)]));
                M0_ijklaElement->InsertEndChild(M0_ijklaEntry);
                permutation = true;
                break;
              }
            }
            if(permutation) continue;
            std::complex<double> Mijkla = 0.;
            for(std::vector<int> helcol : keys) { // color and helicity of the bra <M|
              std::vector<int> hel_full(helcol.begin(), helcol.begin() + A.process.size());
              std::vector<int> col_full(helcol.begin() + A.process.size(), helcol.end());
              std::string hel_string;
              for (int dummy = 0; dummy < hel_full.size(); dummy++) {
                hel_string += std::to_string(hel_full[dummy]);
              }
              std::string col_string;
              for(int dummy = 0; dummy < col_full.size(); dummy++) {
                col_string += std::to_string(col_full[dummy]);
              }
              std::string key = hel_string + col_string;
              std::complex<double> M_bra = std::conj(M0[key]);
              std::vector<int> output_colors = {col_full[i], col_full[j], col_full[k], col_full[l], col_full[a]};
              for(int ci = 0; ci < A.process[i]; ci++) for(int cj = 0; cj < A.process[j]; cj++)
              for(int ck = 0; ck < A.process[k]; ck++) for(int cl = 0; cl < A.process[l]; cl++)
              for(int ca = 0; ca < A.process[a]; ca++) {
                key.replace(key.size() - A.process.size() + i, 1, std::to_string(ci));
                key.replace(key.size() - A.process.size() + j, 1, std::to_string(cj));
                key.replace(key.size() - A.process.size() + k, 1, std::to_string(ck));
                key.replace(key.size() - A.process.size() + l, 1, std::to_string(cl));
                key.replace(key.size() - A.process.size() + a, 1, std::to_string(ca));
                if(M0[key]==0.) continue;
                for(int m = 0; m < 8; m++) for(int n = 0; n < 8; n++) for(int b = 0; b < 8; b++) for(int c = 0; c < 8; c++) {
                  if(j == i)
                    output_colors[1] = ci;
                  else if(k == i)
                    output_colors[2] = ci;
                  else if(l == i)
                    output_colors[3] = ci;
                  else if(a == i)
                    output_colors[4] = ci;

                  if(k == j)
                    output_colors[2] = cj;
                  else if(l == j)
                    output_colors[3] = cj;
                  else if(a == j)
                    output_colors[4] = cj;

                  if(l == k)
                    output_colors[3] = ck;
                  else if(a == k)
                    output_colors[4] = ck;

                  if(a == l)
                    output_colors[4] = cl;

                  std::complex<double> col_factor = 1;

                  if(A.particle_type[a] == 1)
                    col_factor *= lam[c][3*output_colors[4] + ca]/2.;
                  else if(A.particle_type[a] == -1)
                    col_factor *= -lam[c][3*ca + output_colors[4]]/2.;
                  else if(A.particle_type[a] == 2)
                    col_factor *= I*fabc[c][8*ca + output_colors[4]];

                  if(A.particle_type[l] == 1)
                    col_factor *= lam[b][3*output_colors[3] + cl]/2.;
                  else if(A.particle_type[l] == -1)
                    col_factor *= -lam[b][3*cl + output_colors[3]]/2.;
                  else if(A.particle_type[l] == 2)
                    col_factor *= I*fabc[b][8*cl + output_colors[3]];

                  if(A.particle_type[k] == 1)
                    col_factor *= lam[n][3*output_colors[2] + ck]/2.;
                  else if(A.particle_type[k] == -1)
                    col_factor *= -lam[n][3*ck + output_colors[2]]/2.;
                  else if(A.particle_type[k] == 2)
                    col_factor *= I*fabc[n][8*ck + output_colors[2]];

                  if(A.particle_type[j] == 1)
                    col_factor *= lam[m][3*output_colors[1] + cj]/2.;
                  else if(A.particle_type[j] == -1)
                    col_factor *= -lam[m][3*cj + output_colors[1]]/2.;
                  else if(A.particle_type[j] == 2)
                    col_factor *= I*fabc[m][8*cj + output_colors[1]];

                  if(A.particle_type[i] == 1)
                    col_factor *= lam[m][3*output_colors[0] + ci]/2.;
                  else if(A.particle_type[i] == -1)
                    col_factor *= -lam[m][3*ci + output_colors[0]]/2.;
                  else if(A.particle_type[i] == 2)
                    col_factor *= I*fabc[m][8*ci + output_colors[0]];

                  Mijkla += 2.*std::real(fabc[n][8*b+c]*col_factor*M_bra*M0[key]);
                }
              }
            }
            M0_ijkla[std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l) + std::to_string(a)] = Mijkla*average_factor;
            M0_ijkla[std::to_string(i) + std::to_string(j) + std::to_string(a) + std::to_string(k) + std::to_string(l)] = Mijkla*average_factor;
            M0_ijkla[std::to_string(i) + std::to_string(j) + std::to_string(l) + std::to_string(a) + std::to_string(k)] = Mijkla*average_factor;
            M0_ijkla[std::to_string(i) + std::to_string(j) + std::to_string(a) + std::to_string(l) + std::to_string(k)] = -Mijkla*average_factor;
            M0_ijkla[std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(a) + std::to_string(l)] = -Mijkla*average_factor;
            M0_ijkla[std::to_string(i) + std::to_string(j) + std::to_string(l) + std::to_string(k) + std::to_string(a)] = -Mijkla*average_factor;
            M0_ijkla[std::to_string(j) + std::to_string(i) + std::to_string(k) + std::to_string(l) + std::to_string(a)] = Mijkla*average_factor;
            M0_ijkla[std::to_string(j) + std::to_string(i) + std::to_string(a) + std::to_string(k) + std::to_string(l)] = Mijkla*average_factor;
            M0_ijkla[std::to_string(j) + std::to_string(i) + std::to_string(l) + std::to_string(a) + std::to_string(k)] = Mijkla*average_factor;
            M0_ijkla[std::to_string(j) + std::to_string(i) + std::to_string(a) + std::to_string(l) + std::to_string(k)] = -Mijkla*average_factor;
            M0_ijkla[std::to_string(j) + std::to_string(i) + std::to_string(k) + std::to_string(a) + std::to_string(l)] = -Mijkla*average_factor;
            M0_ijkla[std::to_string(j) + std::to_string(i) + std::to_string(l) + std::to_string(k) + std::to_string(a)] = -Mijkla*average_factor;
            configurations.push_back(std::vector<int>({i,j,k,l,a}));
            configurations.push_back(std::vector<int>({i,j,a,k,l}));
            configurations.push_back(std::vector<int>({i,j,l,a,k}));
            configurations.push_back(std::vector<int>({i,j,a,l,k}));
            configurations.push_back(std::vector<int>({i,j,k,a,l}));
            configurations.push_back(std::vector<int>({i,j,l,k,a}));
            configurations.push_back(std::vector<int>({j,i,k,l,a}));
            configurations.push_back(std::vector<int>({j,i,a,k,l}));
            configurations.push_back(std::vector<int>({j,i,l,a,k}));
            configurations.push_back(std::vector<int>({j,i,a,l,k}));
            configurations.push_back(std::vector<int>({j,i,k,a,l}));
            configurations.push_back(std::vector<int>({j,i,l,k,a}));

            M0_ijklaEntry->SetText(std::real(M0_ijkla[std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l) + std::to_string(a)]));
            M0_ijklaElement->InsertEndChild(M0_ijklaEntry);
            double M2_averaged;
            Recola::get_squared_amplitude_rcl(1, power, "LO", M2_averaged);
            std::cout << "i = " << i << ", j = " << j << ", k = " << k << ", l = " << l << ", a = " << a <<  ": <M|Ti.Tj f^{ck, cl, ca} Tk^ck Tl^cl Ta^ca|M> = " << Mijkla*average_factor << "\t" << Mijkla*average_factor/M2_averaged << std::endl;
          }
        }
      }
    }
  }
  pRoot->InsertEndChild(M0_ijklaElement);
  }

  // d^{a,b,c} <M|Ti^a Tj^b Tk^c|M>
  // f^{a,b,c} <M|Ti^a Tj^b Tk^c|M>
  if(((unresolved_str==" g g g" or order=="NLO") and (suffix=="_EW")) or (order=="NNLO")) {
  XMLElement * M0_ijkElement = xmlDoc.NewElement("M0_ijk");
  XMLElement * M1_ijkElement = xmlDoc.NewElement("M1_ijk");
  XMLElement * dM0_ijkElement = xmlDoc.NewElement("dM0_ijk");
  std::vector<std::vector<int>> configurations;
  for (int i = 0; i < A.process.size(); i++) {
    if(A.particle_type[i] == 0) continue;
    for(int j = 0; j < A.process.size(); j++) {
      if(A.particle_type[j] == 0) continue;
      //if(j == i) continue;
      for(int k = 0; k < A.process.size(); k++) {
        if(A.particle_type[k] == 0) continue;
        //if((k == i) or (k == j)) continue;
        std::complex<double> fM0ijk = 0.;
        std::complex<double> fM1ijk = 0.;
        std::complex<double> dM0ijk = 0.;
        bool permutation = false;

        XMLElement * M0_ijkEntry = xmlDoc.NewElement(const_cast<char*>(("d" + std::to_string(i) + std::to_string(j) + std::to_string(k)).c_str()));
        XMLElement * M1_ijkEntry = xmlDoc.NewElement(const_cast<char*>(("d" + std::to_string(i) + std::to_string(j) + std::to_string(k)).c_str()));
        XMLElement * dM0_ijkEntry = xmlDoc.NewElement(const_cast<char*>(("d" + std::to_string(i) + std::to_string(j) + std::to_string(k)).c_str()));
        for(std::vector<int> configuration : configurations) {
          if(configuration==std::vector<int>({i,j,k})) {
            if((i != j) and (i != k) and (j != k)) {
              M0_ijkEntry->SetText(std::real(fM0_ijk[std::to_string(i) + std::to_string(j) + std::to_string(k)]));
              M1_ijkEntry->SetText(std::real(fM1_ijk[std::to_string(i) + std::to_string(j) + std::to_string(k)]));
            }
            else {
              M0_ijkEntry->SetText(0);
              M1_ijkEntry->SetText(0);
            }
            M0_ijkElement->InsertEndChild(M0_ijkEntry);
            M1_ijkElement->InsertEndChild(M1_ijkEntry);
            dM0_ijkEntry->SetText(std::real(dM0_ijk[std::to_string(i) + std::to_string(j) + std::to_string(k)]));
            dM0_ijkElement->InsertEndChild(dM0_ijkEntry);
            permutation = true;
            break;
          }
        }
        if(permutation) continue;
        for(std::vector<int> helcol : keys) { // color and helicity of the bra <M|
          std::vector<int> hel_full(helcol.begin(), helcol.begin() + A.process.size());
          std::vector<int> col_full(helcol.begin() + A.process.size(), helcol.end());
          std::string hel_string;
          for (int dummy = 0; dummy < hel_full.size(); dummy++) {
            hel_string += std::to_string(hel_full[dummy]);
          }
          std::string col_string;
          for(int dummy = 0; dummy < col_full.size(); dummy++) {
            col_string += std::to_string(col_full[dummy]);
          }
          std::string key = hel_string + col_string;
          std::complex<double> M0_bra = std::conj(M0[key]);
          std::vector<int> output_colors = {col_full[i], col_full[j], col_full[k]};
          for(int ci = 0; ci < A.process[i]; ci++) for(int cj = 0; cj < A.process[j]; cj++) for(int ck = 0; ck < A.process[k]; ck++) {
            key.replace(key.size() - A.process.size() + i, 1, std::to_string(ci));
            key.replace(key.size() - A.process.size() + j, 1, std::to_string(cj));
            key.replace(key.size() - A.process.size() + k, 1, std::to_string(ck));

            if(j == i)
              output_colors[1] = ci;
            else if(k == i)
              output_colors[2] = ci;

            if(k == j)
              output_colors[2] = cj;

            if(M0[key] == 0. and M1[key] == 0.) continue;
            for(int a = 0; a < 8; a++) for(int b = 0; b < 8; b++) for(int c = 0; c < 8; c++) {
              std::complex<double> col_factor = 1;

              if(A.particle_type[k] == 1)
                col_factor *= lam[c][3*output_colors[2] + ck]/2.;
              else if(A.particle_type[k] == -1)
                col_factor *= -lam[c][3*ck + output_colors[2]]/2.;
              else if(A.particle_type[k] == 2)
                col_factor *= I*fabc[c][8*ck + output_colors[2]];
              else
                std::cout << "unknown particle type" << std::endl;

              if(A.particle_type[j] == 1)
                col_factor *= lam[b][3*output_colors[1] + cj]/2.;
              else if(A.particle_type[j] == -1)
                col_factor *= -lam[b][3*cj + output_colors[1]]/2.;
              else if(A.particle_type[j] == 2)
                col_factor *= I*fabc[b][8*cj + output_colors[1]];
              else
                std::cout << "unknown particle type" << std::endl;

              if(A.particle_type[i] == 1)
                col_factor *= lam[a][3*output_colors[0] + ci]/2.;
              else if(A.particle_type[i] == -1)
                col_factor *= -lam[a][3*ci + output_colors[0]]/2.;
              else if(A.particle_type[i] == 2)
                col_factor *= I*fabc[a][8*ci + output_colors[0]];
              else
                std::cout << "unknown particle type" << std::endl;
              fM0ijk += col_factor*M0_bra*M0[key]*fabc[a][8*b+c];
              fM1ijk += 2.*std::real(col_factor*M0_bra*M1[key])*fabc[a][8*b+c];
              dM0ijk += col_factor*M0_bra*M0[key]*dabc[a][8*b+c];
            }
          }
        }
        fM0_ijk[std::to_string(i) + std::to_string(j) + std::to_string(k)] = fM0ijk*average_factor;
        fM1_ijk[std::to_string(i) + std::to_string(j) + std::to_string(k)] = fM1ijk*average_factor;
        dM0_ijk[std::to_string(i) + std::to_string(j) + std::to_string(k)] = dM0ijk*average_factor;
        if((i != j) and (i != k) and (j != k)) {
          M0_ijkEntry->SetText(std::real(fM0_ijk[std::to_string(i) + std::to_string(j) + std::to_string(k)]));
          M1_ijkEntry->SetText(std::real(fM1_ijk[std::to_string(i) + std::to_string(j) + std::to_string(k)]));
        }
        else {
          M0_ijkEntry->SetText(0);
          M1_ijkEntry->SetText(0);
        }
        M0_ijkElement->InsertEndChild(M0_ijkEntry);
        M1_ijkElement->InsertEndChild(M1_ijkEntry);
        dM0_ijkEntry->SetText(std::real(dM0_ijk[std::to_string(i) + std::to_string(j) + std::to_string(k)]));
        dM0_ijkElement->InsertEndChild(dM0_ijkEntry);
        // permutations
        fM0_ijk[std::to_string(i) + std::to_string(k) + std::to_string(j)] = -fM0ijk*average_factor;
        fM0_ijk[std::to_string(j) + std::to_string(i) + std::to_string(k)] = -fM0ijk*average_factor;
        fM0_ijk[std::to_string(j) + std::to_string(k) + std::to_string(i)] = fM0ijk*average_factor;
        fM0_ijk[std::to_string(k) + std::to_string(j) + std::to_string(i)] = -fM0ijk*average_factor;
        fM0_ijk[std::to_string(k) + std::to_string(i) + std::to_string(j)] = fM0ijk*average_factor;

        fM1_ijk[std::to_string(i) + std::to_string(k) + std::to_string(j)] = -fM1ijk*average_factor;
        fM1_ijk[std::to_string(j) + std::to_string(i) + std::to_string(k)] = -fM1ijk*average_factor;
        fM1_ijk[std::to_string(j) + std::to_string(k) + std::to_string(i)] = fM1ijk*average_factor;
        fM1_ijk[std::to_string(k) + std::to_string(j) + std::to_string(i)] = -fM1ijk*average_factor;
        fM1_ijk[std::to_string(k) + std::to_string(i) + std::to_string(j)] = fM1ijk*average_factor;

        dM0_ijk[std::to_string(i) + std::to_string(k) + std::to_string(j)] = dM0ijk*average_factor;
        dM0_ijk[std::to_string(j) + std::to_string(i) + std::to_string(k)] = dM0ijk*average_factor;
        dM0_ijk[std::to_string(j) + std::to_string(k) + std::to_string(i)] = dM0ijk*average_factor;
        dM0_ijk[std::to_string(k) + std::to_string(j) + std::to_string(i)] = dM0ijk*average_factor;
        dM0_ijk[std::to_string(k) + std::to_string(i) + std::to_string(j)] = dM0ijk*average_factor;
        configurations.push_back(std::vector<int>({i,j,k}));
        configurations.push_back(std::vector<int>({i,k,j}));
        configurations.push_back(std::vector<int>({j,i,k}));
        configurations.push_back(std::vector<int>({j,k,i}));
        configurations.push_back(std::vector<int>({k,j,i}));
        configurations.push_back(std::vector<int>({k,i,j}));

        dM0_ijk[std::to_string(i) + std::to_string(j) + std::to_string(k)] = dM0ijk*average_factor;

        double M2_averaged;
        Recola::get_squared_amplitude_rcl(1, power, "LO", M2_averaged);
        std::string lock = std::to_string(i) + std::to_string(j) + std::to_string(k);
        std::cout << lock << ": f^{abc} <M0|Ti Tj Tk|M0> = " << fM0_ijk[lock] << "\t" << fM0ijk*average_factor << "\t" << fM0_ijk[lock]/M2_averaged << std::endl;
        std::cout << lock << ": f^{abc} 2 Re <M0|Ti Tj Tk|M1> = " << fM1_ijk[lock] << "\t" << fM1ijk*average_factor << "\t" << fM1_ijk[lock]/M2_averaged << std::endl;
        std::cout << "i = " << i << ", j = " << j << ", k = " << k << ": d^{abc} <M|Ti Tj Tk|M> = " << dM0ijk*average_factor << "\t" << dM0ijk*average_factor/M2_averaged << std::endl;

      }
    }
  }
  pRoot->InsertEndChild(M0_ijkElement);
  pRoot->InsertEndChild(M1_ijkElement);
  pRoot->InsertEndChild(dM0_ijkElement);
  }

  // f^{a,d;b,c} <M|Ti^a {Tj^b, Tk^c} Tl^d|M> + h.c.
  if((order=="NLO" && unresolved_str==" g g") || (unresolved_str==" g g g") || (order=="NNLO")) {
  XMLElement * Q_ijklElement = xmlDoc.NewElement("Q_ijkl");
  std::vector<std::vector<int>> configurations;
  for (int i = 0; i < A.process.size(); i++) {
    if(A.particle_type[i] == 0) continue;
    for(int j = 0; j < A.process.size(); j++) {
      if(A.particle_type[j] == 0) continue;
      for(int k = 0; k < A.process.size(); k++) {
        if(A.particle_type[k] == 0) continue;
        if(k == j) continue;
        for(int l = 0; l < A.process.size(); l++) {
          if(A.particle_type[l] == 0) continue;
          if(l == i) continue;
          XMLElement * Q_ijklEntry = xmlDoc.NewElement(("d" + std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l)).c_str());
          bool permutation = false;
          for(std::vector<int> configuration : configurations) {
            if(configuration==std::vector<int>({i,j,k,l})) {
              Q_ijklEntry->SetText(std::real(Q_ijkl[std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l)]));
              Q_ijklElement->InsertEndChild(Q_ijklEntry);
              permutation = true;
              break;
            }
          }
          if(permutation) continue;
          std::complex<double> ffMijkl = 0.;

          for(std::vector<int> helcol : keys) { // color and helicity of the bra <M|
            std::vector<int> hel_full(helcol.begin(), helcol.begin() + A.process.size());
            std::vector<int> col_full(helcol.begin() + A.process.size(), helcol.end());
            std::string hel_string;
            for (int dummy = 0; dummy < hel_full.size(); dummy++) {
              hel_string += std::to_string(hel_full[dummy]);
            }
            std::string col_string;
            for(int dummy = 0; dummy < col_full.size(); dummy++) {
              col_string += std::to_string(col_full[dummy]);
            }
            std::string key = hel_string + col_string;
            std::complex<double> M_bra = std::conj(M0[key]);
            std::vector<int> output_colors = {col_full[i], col_full[j], col_full[k],  col_full[l]};
            std::vector<int> output_colors2 = {col_full[i], col_full[k], col_full[j],  col_full[l]};
            std::string key2 = key;
            for(int ci = 0; ci < A.process[i]; ci++) for(int cj = 0; cj < A.process[j]; cj++) for(int ck = 0; ck < A.process[k]; ck++) for(int cl = 0; cl < A.process[l]; cl++) {
              key.replace(key.size() - A.process.size() + i, 1, std::to_string(ci));
              key.replace(key.size() - A.process.size() + j, 1, std::to_string(cj));
              key.replace(key.size() - A.process.size() + k, 1, std::to_string(ck));
              key.replace(key.size() - A.process.size() + l, 1, std::to_string(cl));

              key2.replace(key2.size() - A.process.size() + i, 1, std::to_string(ci));
              key2.replace(key2.size() - A.process.size() + k, 1, std::to_string(ck));
              key2.replace(key2.size() - A.process.size() + j, 1, std::to_string(cj));
              key2.replace(key2.size() - A.process.size() + l, 1, std::to_string(cl));
              if((M0[key] == 0.) and (M0[key2] == 0.)) continue;
              for(int a = 0; a < 8; a++) for(int b = 0; b < 8; b++) for(int c = 0; c < 8; c++) for(int d = 0; d < 8; d++) {
                if(j == i)
                  output_colors[1] = ci;
                else if(k == i)
                  output_colors[2] = ci;
                else if(l == i)
                  output_colors[3] = ci;

                if(k == j)
                  output_colors[2] = cj;
                else if(l == j)
                  output_colors[3] = cj;

                if(l == k)
                  output_colors[3] = ck;

                if(k == i)
                  output_colors2[1] = ci;
                else if(j == i)
                  output_colors2[2] = ci;
                else if(l == i)
                  output_colors2[3] = ci;

                if(j == k)
                  output_colors2[2] = ck;
                else if(l == k)
                  output_colors2[3] = ck;

                if(l == j)
                  output_colors2[3] = cj;

                std::complex<double> col_factor1 = 1;
                if(A.particle_type[l] == 1)
                  col_factor1 *= lam[b][3*output_colors[3] + cl]/2.;
                else if(A.particle_type[l] == -1)
                  col_factor1 *= -lam[b][3*cl + output_colors[3]]/2.;
                else if(A.particle_type[l] == 2)
                  col_factor1 *= I*fabc[b][8*cl + output_colors[3]];

                if(A.particle_type[k] == 1)
                  col_factor1 *= lam[d][3*output_colors[2] + ck]/2.;
                else if(A.particle_type[k] == -1)
                  col_factor1 *= -lam[d][3*ck + output_colors[2]]/2.;
                else if(A.particle_type[k] == 2)
                  col_factor1 *= I*fabc[d][8*ck + output_colors[2]];

                if(A.particle_type[j] == 1)
                  col_factor1 *= lam[c][3*output_colors[1] + cj]/2.;
                else if(A.particle_type[j] == -1)
                  col_factor1 *= -lam[c][3*cj + output_colors[1]]/2.;
                else if(A.particle_type[j] == 2)
                  col_factor1 *= I*fabc[c][8*cj + output_colors[1]];

                if(A.particle_type[i] == 1)
                  col_factor1 *= lam[a][3*output_colors[0] + ci]/2.;
                else if(A.particle_type[i] == -1)
                  col_factor1 *= -lam[a][3*ci + output_colors[0]]/2.;
                else if(A.particle_type[i] == 2)
                  col_factor1 *= I*fabc[a][8*ci + output_colors[0]];

                std::complex<double> col_factor2 = 1;
                if(A.particle_type[l] == 1)
                  col_factor2 *= lam[b][3*output_colors2[3] + cl]/2.;
                else if(A.particle_type[l] == -1)
                  col_factor2 *= -lam[b][3*cl + output_colors2[3]]/2.;
                else if(A.particle_type[l] == 2)
                  col_factor2 *= I*fabc[b][8*cl + output_colors2[3]];

                if(A.particle_type[j] == 1)
                  col_factor2 *= lam[c][3*output_colors2[2] + cj]/2.;
                else if(A.particle_type[j] == -1)
                  col_factor2 *= -lam[c][3*cj + output_colors2[2]]/2.;
                else if(A.particle_type[j] == 2)
                  col_factor2 *= I*fabc[c][8*cj + output_colors2[2]];

                if(A.particle_type[k] == 1)
                  col_factor2 *= lam[d][3*output_colors2[1] + ck]/2.;
                else if(A.particle_type[k] == -1)
                  col_factor2 *= -lam[d][3*ck + output_colors2[1]]/2.;
                else if(A.particle_type[k] == 2)
                  col_factor2 *= I*fabc[d][8*ck + output_colors2[1]];

                if(A.particle_type[i] == 1)
                  col_factor2 *= lam[a][3*output_colors2[0] + ci]/2.;
                else if(A.particle_type[i] == -1)
                  col_factor2 *= -lam[a][3*ci + output_colors2[0]]/2.;
                else if(A.particle_type[i] == 2)
                  col_factor2 *= I*fabc[a][8*ci + output_colors2[0]];

                for(int s = 0; s < 8; s++)
                  ffMijkl += 2*std::real(M_bra*(col_factor1*M0[key] + col_factor2*M0[key2]))*fabc[a][8*b + s]*fabc[c][8*d + s];
              }
            }
          }
          Q_ijkl[std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l)] = ffMijkl*average_factor;
          Q_ijkl[std::to_string(i) + std::to_string(k) + std::to_string(j) + std::to_string(l)] = -ffMijkl*average_factor;
          Q_ijkl[std::to_string(l) + std::to_string(j) + std::to_string(k) + std::to_string(i)] = -ffMijkl*average_factor;
          Q_ijkl[std::to_string(l) + std::to_string(k) + std::to_string(j) + std::to_string(i)] = ffMijkl*average_factor;

          configurations.push_back(std::vector<int>({i,j,k,l}));
          configurations.push_back(std::vector<int>({l,j,k,i}));
          configurations.push_back(std::vector<int>({i,k,j,l}));
          configurations.push_back(std::vector<int>({l,k,j,i}));

          Q_ijklEntry->SetText(std::real(Q_ijkl[std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l)]));
          Q_ijklElement->InsertEndChild(Q_ijklEntry);

          double M2_averaged;
          Recola::get_squared_amplitude_rcl(1, power, "LO", M2_averaged);
          std::cout << "i = " << i << ", j = " << j << ", k = " << k << ", l = " << l << ": f^{abc} f^{cde} <M|Ti^a {Tj^b, Tk^d} Tl^e|M> = " << ffMijkl*average_factor << "\t" << ffMijkl*average_factor/M2_averaged << std::endl;
        }
      }
    }
  }
  pRoot->InsertEndChild(Q_ijklElement);
  }

  std::cout << "Filled the Hashmaps" << std::endl;
  std::cout << "Save results to " << "results/ColorCorrelators/" + process_str + suffix + ".xml" << std::endl;
  XMLError eResult = xmlDoc.SaveFile(const_cast<char*>(("results/ColorCorrelators/" + process_str + suffix + ".xml").c_str()));
  return 0;
}
