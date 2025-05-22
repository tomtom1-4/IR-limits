#include "main.hpp"

//using namespace std;
//using namespace Stripper;

// Evaluation of input
std::string process_str = "u u~ -> A A";
std::string unresolved_str = " g";
std::string limit = "soft";
//std::string process_full_str = process_str + unresolved_str;
std::string process_full_str = "u u~ -> A A g";

const int nBorn = 4;
const int power = 0;
const std::array<unsigned, 3> powerNonQCD = {2,0,0};
const std::vector<int> flavor = {0,0,1,1};
const std::string order = "NNLO";
const std::string suffix = "_QED"; // "" or "_EW" or "_QED"
std::vector<double> M2_custom;

int main() {
  // switch floating point type here
  using Real = double;
  int loopOrder = 0;
  if(order == "LO") loopOrder = 0;
  if(order == "NLO") loopOrder = 1;
  if(order == "NNLO") loopOrder = 2;

  // user defined matrix element configuration
  Stripper::Model::config((powerNonQCD[0]!=0),false,false,false,powerNonQCD[0],powerNonQCD[1],powerNonQCD[2],1.,0.,0.);
  const Stripper::Process process_full(process_full_str);
  const Stripper::Process process(process_str);
  Stripper::Model::nf = n_f;
  Stripper::Model::print(std::cout);

  // phase space point
  const unsigned n = 3*(nBorn - 2) - 4;

  amplitude A(process_str);
  amplitude A_full(process_full_str);
  const int nUnresolved = A_full.process.size() - A.process.size();

  bool custom = false;
  if((process_full_str == "u u~ -> u u~ e- e+ g g") and (order == "NLO")) {
    custom = true;
    M2_custom.push_back(-4.58613727787852e-26);
    M2_custom.push_back(-3.05566873598841e-22);
    M2_custom.push_back(-6.78415558513149e-20);
    M2_custom.push_back(-1.2213501061979e-17);
    M2_custom.push_back(-1.93961648681395e-15);
    M2_custom.push_back(-2.82993140226372e-13);
    M2_custom.push_back(-3.88880026114949e-11);
    M2_custom.push_back(-5.11393342703531e-09);
    M2_custom.push_back(-6.5040578409581e-07);
    M2_custom.push_back(-8.05858355886322e-05);
    M2_custom.push_back( -0.0097845849231643);
    M2_custom.push_back( -1.16517247862666);
    M2_custom.push_back( -136.406428204994);
  }
  else if((process_full_str == "e- e+ -> u u~ u u~ g g") and (order == "NLO")) {
    custom = true;
    M2_custom.push_back(8.42714370716502e-25);
    M2_custom.push_back(-1.07948493039382e-22);
    M2_custom.push_back(-3.1589420159308e-20);
    M2_custom.push_back(-1.25043991989418e-17);
    M2_custom.push_back(-1.01797400996888e-15);
    M2_custom.push_back(-1.50917187473323e-13);
    M2_custom.push_back(-2.09150310626981e-11);
    M2_custom.push_back(-2.76442260780936e-09);
    M2_custom.push_back(-3.52770343464854e-07);
    M2_custom.push_back(-4.38128658522702e-05);
    M2_custom.push_back( -0.00532339243608503);
    M2_custom.push_back( -0.632583283928193);
    M2_custom.push_back( -72.7129589447317);
  }
  else if((process_full_str == "u u~ -> A A g") and (order == "NNLO")) {
    custom = true;
    M2_custom.push_back(std::pow(gs, 2*(power + nUnresolved) + 4)*4.83538613801978e-06);
    M2_custom.push_back(std::pow(gs, 2*(power + nUnresolved) + 4)*1.54476748572392e-05);
    M2_custom.push_back(std::pow(gs, 2*(power + nUnresolved) + 4)*5.49200500812953e-05);
    M2_custom.push_back(std::pow(gs, 2*(power + nUnresolved) + 4)*0.000191162117465791);
    M2_custom.push_back(std::pow(gs, 2*(power + nUnresolved) + 4)*0.000630408722842142);
    M2_custom.push_back(std::pow(gs, 2*(power + nUnresolved) + 4)*0.00197007173283094);
    M2_custom.push_back(std::pow(gs, 2*(power + nUnresolved) + 4)*0.0058781023986536);
    M2_custom.push_back(std::pow(gs, 2*(power + nUnresolved) + 4)*0.016875328259249);
    M2_custom.push_back(std::pow(gs, 2*(power + nUnresolved) + 4)*0.0469193849876823);
    M2_custom.push_back(std::pow(gs, 2*(power + nUnresolved) + 4)*0.126995124044529);
    M2_custom.push_back(std::pow(gs, 2*(power + nUnresolved) + 4)*0.335987583035483);
    M2_custom.push_back( std::pow(gs, 2*(power + nUnresolved) + 4)*0.871716694258973);
    M2_custom.push_back( std::pow(gs, 2*(power + nUnresolved) + 4)*2.22409119793593);
    M2_custom.push_back( std::pow(gs, 2*(power + nUnresolved) + 4)*5.59076343769503);
    M2_custom.push_back( std::pow(gs, 2*(power + nUnresolved) + 4)*13.8716256713624);
    M2_custom.push_back( std::pow(gs, 2*(power + nUnresolved) + 4)*34.0212692292854);
    M2_custom.push_back( std::pow(gs, 2*(power + nUnresolved) + 4)*82.5789232014353);
    M2_custom.push_back( std::pow(gs, 2*(power + nUnresolved) + 4)*198.578361773121);
    M2_custom.push_back( std::pow(gs, 2*(power + nUnresolved) + 4)*473.502601645154);
    M2_custom.push_back( std::pow(gs, 2*(power + nUnresolved) + 4)*1120.39590157353);
    M2_custom.push_back( std::pow(gs, 2*(power + nUnresolved) + 4)*2632.49086045922);
  }
  else if((process_full_str == "d d~ -> g g g g") and (order == "NLO")) {
    custom = true;
    M2_custom.push_back( -1.07529890579258e-06);
    M2_custom.push_back( -2.60045451018945e-06);
    M2_custom.push_back( -6.03844054006154e-06);
    M2_custom.push_back( -1.36442201108769e-05);
    M2_custom.push_back( -3.02934353742203e-05);
    M2_custom.push_back( -6.65012200836531e-05);
    M2_custom.push_back( -0.00014489308977931);
    M2_custom.push_back( -0.000314054902497589);
    M2_custom.push_back( -0.000678129355886073);
    M2_custom.push_back( -0.00145996160406465);
    M2_custom.push_back( -0.00313561407655795);
    M2_custom.push_back( -0.00672047760786235);
    M2_custom.push_back( -0.0143769182369063);
    M2_custom.push_back( -0.0307029256774041);
    M2_custom.push_back( -0.0654610931772924);
    M2_custom.push_back( -0.139349808857436);
    M2_custom.push_back( -0.296191670082831);
    M2_custom.push_back( -0.628641187086145);
    M2_custom.push_back( -1.33233894317267);
    M2_custom.push_back( -2.81985020103629);
    M2_custom.push_back( -5.96010598030109);
    M2_custom.push_back( -12.5810152286265);
    M2_custom.push_back( -26.5233468681044);
    M2_custom.push_back( -55.8481778022263);
    M2_custom.push_back( -117.455931890383);
    M2_custom.push_back( -246.742072703789);
    M2_custom.push_back( -517.762403376488);
    M2_custom.push_back( -1085.3088410625);
    M2_custom.push_back( -2272.629051835855);
    M2_custom.push_back( -4754.1273305696);
    M2_custom.push_back( -9935.61738850097);
    M2_custom.push_back( -20745.3830428635);
    M2_custom.push_back( -43275.7168515828);
    M2_custom.push_back( -90202.9844854199);
    M2_custom.push_back( -187857.783735964);
    M2_custom.push_back( -390930.416613994);
    M2_custom.push_back( -812924.988423128);
    M2_custom.push_back( -1688056.64421381);
    M2_custom.push_back( -3500582.94057637);
    M2_custom.push_back( -7254873.16785345);
    M2_custom.push_back( -15029752.4463456);
    M2_custom.push_back( -34014957.0833781);
    M2_custom.push_back( -32301871.466692);
    M2_custom.push_back( +201602230.12205);
    M2_custom.push_back( 241911353.271096);
    M2_custom.push_back( -9861511028.22394);
    M2_custom.push_back( -229797749693.066);
  }

  int delta_power = 0;
  if (order == "NLO") delta_power = 2;
  std::cout.precision(8);
  // Initialize txt outstream (delete previous file)

  std::ofstream outfile1;
  outfile1.open ("results/" + process_str + " +" + unresolved_str + "_" + order + suffix + ".txt");
  outfile1.close();
  outfile1.open ("results/" + process_str + " +" + unresolved_str + "_" + order + suffix + ".txt");

  if(!custom and (order=="NLO" or order=="LO")) {
    std::cout << "Using Recola to compute tree-level and one-loop amplitudes" << std::endl;
    // Recola Settings
    Recola::set_reduction_mode_rcl(4);
    Recola::set_print_level_amplitude_rcl(2);
    Recola::set_alphas_rcl(gs * gs/(4 * M_PI), mu, n_f);
    Recola::set_mu_ir_rcl(mu);
    Recola::set_mu_uv_rcl(mu);
    Recola::set_delta_ir_rcl(0., M_PI*M_PI/12.);
    Recola::set_momenta_correction_rcl(true);

    Recola::use_alpha0_scheme_rcl(e*e/4./M_PI);
    //Recola::set_pole_mass_z_rcl(1.e8, 0.0001);
    //Recola::set_pole_mass_w_rcl(1.e15, 0.0001);

    // Define & generate process
    Recola::define_process_rcl(1, process_str, order);
    Recola::define_process_rcl(2, process_full_str, order);
    Recola::set_otter_mode_rcl(1, "oneloop_qp");
    Recola::set_otter_mode_rcl(2, "oneloop_qp");

    Recola::generate_processes_rcl();
  }

  XMLDocument xmlDoc;
  XMLError eResult = xmlDoc.LoadFile(const_cast<char*>(("results/ColorCorrelators/" + process_str + suffix + ".xml").c_str()));
  XMLNode * pRoot = xmlDoc.FirstChild();
  if(pRoot == nullptr) {
    std::cout << "Could not find file: " << "results/ColorCorrelators/" + process_str + suffix + ".xml" << std::endl;
    return -1;
  }
  // Define initial state
  PSF::PhaseSpace dumn = PSF::Splitting(nBorn - 2, COM);
  PSF::PhaseSpace pp;// = Splitting(nBorn - 2, COM);
  XMLElement * ppElement = pRoot->FirstChildElement("PhasePhacePoint");
  for(int i = 0; i < nBorn; i++) {
    PSF::Momentum pi;
    for(int mu = 0; mu < 4; mu++) {
      XMLElement * ppEntry = ppElement->FirstChildElement(const_cast<char*>(("p" + std::to_string(i) + "_" + std::to_string(mu)).c_str()));
      if(ppEntry != nullptr) eResult = ppEntry->QueryDoubleText(&pi.components[mu]);
    }
    pp.momenta.push_back(pi);
  }
  pp.print();
  double pp_rcl[nBorn][4], ppFull_rcl[nBorn + nUnresolved][4], pp_arr[4*nBorn];
  for(int i = 0 ; i < nBorn; i++) {
    for(int j = 0; j < 4; j++)
      pp_rcl[i][j] = (i<2?-1.:1.)*pp.momenta[i].components[j];
  }

  // Get non-vanishing helicity and color-configurations
  std::unordered_map<std::string, double> M0_ij, M1_ij, M2_ij, M0_ijk, M1_ijk, dM0_ijk, M0_ijkl, M1_ijkl, M0_ijklab, Q_ijkl, M0_ijkla; // color correlators

  double average_factor = 1./4.; // Spin of the initial state
  double average_factor_full = 1./4.;
  {
  for(int i = 0; i < 2; i++) {
    average_factor *= 1./double(A.process[i]);
    average_factor_full *= 1./double(A_full.process[i]);
  }
  std::vector<std::string> particles = A.process_particles;
  std::vector<std::string> particles_full = A_full.process_particles;
  for(auto s : particles) std::cout << s << ", " ;
  std::cout << std::endl;
  particles.erase(particles.begin(), particles.begin() + 2);
  for(auto s : particles) std::cout << s << ", " ;
  std::cout << std::endl;
  particles_full.erase(particles_full.begin(), particles_full.begin() + 2);
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
  while(particles_full.size() > 0) {
    std::string particle = particles_full[0];
    int occurances = 1;
    int i = 1;
    while(i < particles_full.size()) {
      if(particles_full[i] == particle) {
        occurances += 1;
        particles_full.erase(particles_full.begin() + i);
      }
      else {
        i += 1;
      }
    }
    particles_full.erase(particles_full.begin());
    average_factor_full *= 1./factorial(occurances);
  }
  }
  std::cout << "average_factor = " << average_factor << std::endl;
  std::cout << "average_factor_full = " << average_factor_full << std::endl;
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
  // Get matrix elements
  double M0, M1, M2;
  XMLElement * M0_Element = pRoot->FirstChildElement("M0");
  eResult = M0_Element->QueryDoubleText(&M0);
  XMLElement * M1_Element = pRoot->FirstChildElement("M1");
  eResult = M1_Element->QueryDoubleText(&M1);
  if(order=="NNLO") {
    XMLElement * M2_Element = pRoot->FirstChildElement("M2");
    eResult = M2_Element->QueryDoubleText(&M2);
  }
  // Get spin correlator
  std::unordered_map<std::string, std::complex<double>> SC0, SC1, SC2;
  if(limit == "collinear") {
    XMLElement * SC0_iElement = pRoot->FirstChildElement("SC0");
    XMLElement * SC1_iElement = pRoot->FirstChildElement("SC1");
    XMLElement * SC2_iElement = pRoot->FirstChildElement("SC2");
    for(int i = 0; i < A.process.size(); i++) {
      if(A.particle_type[i] == 0) continue;
      for(int s1 = -1; s1 <= 1; s1 += 2) {
        for(int s2 = -1; s2 <= 1; s2 += 2) {
          XMLElement * SC0_iEntry_real = SC0_iElement->FirstChildElement(const_cast<char*>(("r" + std::to_string(i) + std::to_string(s1) + std::to_string(s2)).c_str()));
          XMLElement * SC0_iEntry_imag = SC0_iElement->FirstChildElement(const_cast<char*>(("i" + std::to_string(i) + std::to_string(s1) + std::to_string(s2)).c_str()));
          XMLElement * SC1_iEntry_real = SC1_iElement->FirstChildElement(const_cast<char*>(("r" + std::to_string(i) + std::to_string(s1) + std::to_string(s2)).c_str()));
          XMLElement * SC1_iEntry_imag = SC1_iElement->FirstChildElement(const_cast<char*>(("i" + std::to_string(i) + std::to_string(s1) + std::to_string(s2)).c_str()));

          XMLElement * SC2_iEntry_real = SC2_iElement->FirstChildElement(const_cast<char*>(("r" + std::to_string(i) + std::to_string(s1) + std::to_string(s2)).c_str()));
          XMLElement * SC2_iEntry_imag = SC2_iElement->FirstChildElement(const_cast<char*>(("i" + std::to_string(i) + std::to_string(s1) + std::to_string(s2)).c_str()));
          double real0, imag0, real1, imag1, real2, imag2;
          if(SC0_iEntry_real != nullptr) eResult = SC0_iEntry_real->QueryDoubleText(&real0);
          if(SC0_iEntry_imag != nullptr) eResult = SC0_iEntry_imag->QueryDoubleText(&imag0);
          if(SC1_iEntry_real != nullptr) eResult = SC1_iEntry_real->QueryDoubleText(&real1);
          if(SC1_iEntry_imag != nullptr) eResult = SC1_iEntry_imag->QueryDoubleText(&imag1);
          if(SC2_iEntry_real != nullptr) eResult = SC2_iEntry_real->QueryDoubleText(&real2);
          if(SC2_iEntry_imag != nullptr) eResult = SC2_iEntry_imag->QueryDoubleText(&imag2);
          SC0[std::to_string(i) + std::to_string(s1) + std::to_string(s2)] = real0 + I*imag0;
          SC1[std::to_string(i) + std::to_string(s1) + std::to_string(s2)] = real1 + I*imag1;
          SC2[std::to_string(i) + std::to_string(s1) + std::to_string(s2)] = real2 + I*imag2;

          std::cout << std::to_string(i) + std::to_string(s1) + std::to_string(s2) << "\t" << SC0[std::to_string(i) + std::to_string(s1) + std::to_string(s2)]
                                                                                           << SC1[std::to_string(i) + std::to_string(s1) + std::to_string(s2)]
                                                                                           << SC2[std::to_string(i) + std::to_string(s1) + std::to_string(s2)] << std::endl;
        }
      }
    }
  }
  else if(limit == "soft") {
    // Get color correlators
    // <M|Ti.Tj|M>
    XMLElement * M0_ijElement = pRoot->FirstChildElement("M0_ij");
    XMLElement * M1_ijElement = pRoot->FirstChildElement("M1_ij");
    XMLElement * M2_ijElement = pRoot->FirstChildElement("M2_ij");
    for (int i = 0; i < A.process.size(); i++) {
      if(A.particle_type[i] == 0) continue;
      for(int j = 0; j < A.process.size(); j++) {
        if(A.particle_type[j] == 0) continue;
        if(M0_ijElement != nullptr) {
          XMLElement * M0_ijEntry = M0_ijElement->FirstChildElement(const_cast<char*>(("d" + std::to_string(i) + std::to_string(j)).c_str()));
          if(M0_ijEntry != nullptr) eResult = M0_ijEntry->QueryDoubleText(&M0_ij[std::to_string(i) + std::to_string(j)]);
        }
        if(M1_ijElement != nullptr) {
          XMLElement * M1_ijEntry = M1_ijElement->FirstChildElement(const_cast<char*>(("d" + std::to_string(i) + std::to_string(j)).c_str()));
          if(M1_ijEntry != nullptr) eResult = M1_ijEntry->QueryDoubleText(&M1_ij[std::to_string(i) + std::to_string(j)]);
        }
        if(M2_ijElement != nullptr) {
          XMLElement * M2_ijEntry = M2_ijElement->FirstChildElement(const_cast<char*>(("d" + std::to_string(i) + std::to_string(j)).c_str()));
          if(M2_ijEntry != nullptr) eResult = M2_ijEntry->QueryDoubleText(&M2_ij[std::to_string(i) + std::to_string(j)]);
        }
      }
    }
    // <M|Ti.Tj Tk.Tl|M>
    if(unresolved_str==" g g" || unresolved_str==" g g g") {
      XMLElement * M0_ijklElement = pRoot->FirstChildElement("M0_ijkl");
      XMLElement * M1_ijklElement = pRoot->FirstChildElement("M1_ijkl");
      for (int i = 0; i < A.process.size(); i++) {
        if(A.particle_type[i] == 0) continue;
        for(int j = 0; j < A.process.size(); j++) {
          if(A.particle_type[j] == 0) continue;
          for(int k = 0; k < A.process.size(); k++) {
            if(A.particle_type[k] == 0) continue;
            for(int l = 0; l < A.process.size(); l++) {
              if(A.particle_type[l] == 0) continue;
              if(M0_ijklElement != nullptr) {
                XMLElement * M0_ijklEntry = M0_ijklElement->FirstChildElement(const_cast<char*>(("d" + std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l)).c_str()));
                if(M0_ijklEntry != nullptr) eResult = M0_ijklEntry->QueryDoubleText(&M0_ijkl[std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l)]);
              }
              if(M1_ijklElement != nullptr) {
                XMLElement * M1_ijklEntry = M1_ijklElement->FirstChildElement(const_cast<char*>(("d" + std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l)).c_str()));
                if(M1_ijklEntry != nullptr) eResult = M1_ijklEntry->QueryDoubleText(&M1_ijkl[std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l)]);
              }
            }
          }
        }
      }
    }

    // f^{a,d;b,c} <M|Ti^a {Tj^b, Tk^c} Tl^d|M> + h.c.
    if((order=="NLO" && unresolved_str==" g g") || unresolved_str==" g g g" or order=="NNLO") {
      XMLElement * Q_ijklElement = pRoot->FirstChildElement("Q_ijkl");
      if(Q_ijklElement != nullptr) {
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
                XMLElement * Q_ijklEntry = Q_ijklElement->FirstChildElement(const_cast<char*>(("d" + std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l)).c_str()));
                if(Q_ijklEntry != nullptr) Q_ijklEntry->QueryDoubleText(&Q_ijkl[std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l)]);
              }
            }
          }
        }
      }
    }

    // f^{a,b,c} <M|Ti^a Tj^b Tk^c|M>
    if(unresolved_str==" g g g" or order=="NLO" or order=="NNLO") {
      XMLElement * M0_ijkElement = pRoot->FirstChildElement("M0_ijk");
      XMLElement * M1_ijkElement = pRoot->FirstChildElement("M1_ijk");
      if(M0_ijkElement != nullptr) {
        for(int i = 0; i < A.process.size(); i++) {
          if(A.particle_type[i] == 0) continue;
          for(int j = 0; j < A.process.size(); j++) {
            if(A.particle_type[j] == 0) continue;
            if(j == i) continue;
            for(int k = 0; k < A.process.size(); k++) {
              if(A.particle_type[k] == 0) continue;
              if(k == i or k == j) continue;
              XMLElement * M0_ijkEntry = M0_ijkElement->FirstChildElement(const_cast<char*>(("d" + std::to_string(i) + std::to_string(j) + std::to_string(k)).c_str()));
              if(M0_ijkEntry != nullptr) M0_ijkEntry->QueryDoubleText(&M0_ijk[std::to_string(i) + std::to_string(j) + std::to_string(k)]);
              XMLElement * M1_ijkEntry = M1_ijkElement->FirstChildElement(const_cast<char*>(("d" + std::to_string(i) + std::to_string(j) + std::to_string(k)).c_str()));
              if(M1_ijkEntry != nullptr) M1_ijkEntry->QueryDoubleText(&M1_ijk[std::to_string(i) + std::to_string(j) + std::to_string(k)]);
            }
          }
        }
      }
    }

    // d^{a,b,c} <M|Ti^a Tj^b Tk^c|M>
    if(unresolved_str==" g g g" or order=="NLO") {
      XMLElement * dM0_ijkElement = pRoot->FirstChildElement("dM0_ijk");
      if(dM0_ijkElement != nullptr) {
        for(int i = 0; i < A.process.size(); i++) {
          if(A.particle_type[i] == 0) continue;
          for(int j = 0; j < A.process.size(); j++) {
            if(A.particle_type[j] == 0) continue;
            for(int k = 0; k < A.process.size(); k++) {
              if(A.particle_type[k] == 0) continue;
              XMLElement * dM0_ijkEntry = dM0_ijkElement->FirstChildElement(const_cast<char*>(("d" + std::to_string(i) + std::to_string(j) + std::to_string(k)).c_str()));
              if(dM0_ijkEntry != nullptr) dM0_ijkEntry->QueryDoubleText(&dM0_ijk[std::to_string(i) + std::to_string(j) + std::to_string(k)]);
            }
          }
        }
      }
    }

    // <M|Ti.Tj f^{ck, cl, ca} Tk^ck Tl^cl Ta^ca|M> + c.c.
    if((order=="NLO") and (unresolved_str== " g g") and (suffix == "_EW")) {
      XMLElement * M0_ijklaElement = pRoot->FirstChildElement("M0_ijkla");
      if(M0_ijklaElement != nullptr) {
        for(int i = 0; i < A.process.size(); i++) {
          if(A.particle_type[i] == 0) continue;
          for(int j = 0; j < A.process.size(); j++) {
            if(A.particle_type[j] == 0) continue;
            for(int k = 0; k < A.process.size(); k++) {
              if(A.particle_type[k] == 0) continue;
              for(int l = 0; l < A.process.size(); l++) {
                if(A.particle_type[l] == 0) continue;
                for(int a = 0; a < A.process.size(); a++) {
                  if(A.particle_type[a] == 0) continue;
                  XMLElement * M0_ijklaEntry = M0_ijklaElement->FirstChildElement(const_cast<char*>(("d" + std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l) + std::to_string(a)).c_str()));
                  if(M0_ijklaEntry != nullptr) M0_ijklaEntry->QueryDoubleText(&M0_ijkla[std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l) + std::to_string(a)]);
                }
              }
            }
          }
        }
      }
    }
  }

  std::cout << "Filled the Hashmaps" << std::endl;

  // Define ClusterTree
  std::vector<PSF::Tree<PSF::Cluster>> trees = PSF::GenTrees(nUnresolved);
  for(int tree_counter = 0; tree_counter < trees.size(); tree_counter++) {
    PSF::Tree<PSF::Cluster> clusterTree = trees[tree_counter];
    std::cout << "tree_counter = " << tree_counter << std:: endl;
    clusterTree.print();
    std::cout << "#########################################################################" << std::endl;
  }
  PSF::Tree<PSF::Cluster> tree;
  if(nUnresolved == 1) tree = trees[0];
  else if(nUnresolved == 2) tree = trees[1];
  else if(nUnresolved == 3) tree = trees[2];
  else std::cout << "Define Tree for " << nUnresolved << " unresolved partons" << std::endl;
  tree.print();
  std::vector<PSF::Tree<PSF::Cluster>> sectors = PSF::GenSectors(flavor, tree, nBorn);

  for(int sec_counter = 0; sec_counter < 1; sec_counter++) {
  PSF::Tree<PSF::Cluster> clusterTree = sectors[sec_counter];
  std::cout << "reference = " << clusterTree.getRoot()->children[0]->data.reference << std::endl;

  // Generate phase-space points
  double scale = std::sqrt(0.1);
  std::vector<double> etas(nUnresolved);
  std::vector<double> phis(nUnresolved);
  double dummy = PSF::rnd(0., 1.);
  for(int i = 0; i < nUnresolved; i++) {
    etas[i] = PSF::rnd(0.1, 0.9);
    phis[i] = PSF::rnd(0., 1.);
  }
  //etas[0] = 1./3.;
  //etas[1] = 0.5;
  //phis[0] = 0.;
  //phis[1] = 1./6.;
  double increment = std::sqrt(0.5);
  int counter = 0;
  while (scale > 3.e-4) {
    scale *= increment;
    std::vector<std::vector<std::vector<double>>> xParFull;
    int level_int = 1;
    std::vector<PSF::TreeNode<PSF::Cluster>*> level = clusterTree.getLevel(level_int);
    int unresolved_counter = 0;
    while(level.size() > 0) {
      int unresolved_level = 0;
      for(PSF::TreeNode<PSF::Cluster>* node : level){
        unresolved_level += node->data.unresolved;

        std::vector<std::vector<double>> xPar;
        for(int c = 0; c < node->data.unresolved; c++) {
          double xi, eta;
          if (limit == "soft") {
            xi = std::pow(scale, 1);
            eta = etas[unresolved_counter];
          }
          else if(limit == "collinear") {
            xi = etas[unresolved_counter];
            eta = std::pow(scale, 1);
          }
          double phi = phis[unresolved_counter];
          unresolved_counter++;
          std::vector<double> xPar_c = {eta, xi, phi};
          xPar.push_back(xPar_c);
        }
        if(node->children.size() > 1)
          xParFull.push_back(xPar);
      }
      level_int++;
      level = clusterTree.getLevel(level_int);
    }
    PSF::PhaseSpace ppFull = PSF::GenMomenta2(pp, clusterTree, xParFull);

    ppFull.print();

    double pp_full[4*(nBorn+nUnresolved)];
    for(int i = 0; i < ppFull.momenta.size(); i++) {
      for(int j = 0; j < 4; j++) {
        pp_full[4*i+j] = (i<2?-1.:1.)*ppFull.momenta[i].components[j];
        ppFull_rcl[i][j] = (i<2?-1.:1.)*ppFull.momenta[i].components[j];
      }
    }
    std::vector<Stripper::Momentum<Real>> ppFull_Stripper, pp_Stripper;
    for(int i = 0; i < nBorn; i++) {
      Stripper::Momentum<Real> p_Stripper = {pp_rcl[i][0], pp_rcl[i][1], pp_rcl[i][2], pp_rcl[i][3],0.,0.};
      pp_Stripper.push_back(p_Stripper);
    }
    for(int i = 0; i < nBorn+nUnresolved; i++) {
      Stripper::Momentum<Real> p_Stripper = {ppFull_rcl[i][0], ppFull_rcl[i][1], ppFull_rcl[i][2], ppFull_rcl[i][3],0.,0.};
      //Stripper::Momentum<double> p_Stripper = {ppFull.momenta[i].components[0], ppFull.momenta[i].components[1], ppFull.momenta[i].components[2], ppFull.momenta[i].components[3],0.,0.};
      ppFull_Stripper.push_back(p_Stripper);
    }

    double M2_Full = 0.;
    if(!custom) {
      if(order == "LO" or order == "NLO") {
        double M2_test = 0;
        // Recola one-loop
        Recola::compute_process_rcl(1, pp_rcl, order);
        Recola::get_squared_amplitude_rcl(1, power + delta_power/2, order, M2_test);
        Recola::compute_process_rcl(2, ppFull_rcl, order);
        Recola::get_squared_amplitude_rcl(2, power + nUnresolved + delta_power/2, order, M2_Full);
        //Stripper one-loop
        //double M2_Full_Stripper = 0.;
        //Stripper::Born<Real> me(process_full, ppFull_Stripper);
        //M2_Full_Stripper = Stripper::toDouble(me())*average_factor_full*std::pow(gs, 2*(power));
        //M2_Full = M2_Full_Stripper;
        //Stripper::OneLoop<double> me(process_full,ppFull_Stripper);
        //for(int i = 0; i <= 2; i++) M2_Full_Stripper += ((me()[i])*std::pow(std::log(mu*mu/COM/COM), i))*average_factor_full*std::pow(gs, 2*(power + nUnresolved) + 4);
        //std::cout << "gs = " << gs << std::endl;
        //std::cout << M2_Full/average_factor << "\t" << M2_Full_Stripper/average_factor << "\t" << M2_Full/M2_Full_Stripper << std::endl;
        //M2_Full = M2_Full_Stripper*2.;
        //std::cout << "M2_Stripper_full = " << M2_Full_Stripper << "\t" << M2_Full << "\t" << M2_Full/M2_Full_Stripper << std::endl;
      }
      else if(order == "NNLO") {
        Stripper::Born<Real> me0(process_full, ppFull_Stripper);
        double M2_Born_Full = me0();
        Stripper::TwoLoop<Real> me(process_full, ppFull_Stripper);
        //Stripper::OneLoop<double> me(process_full,ppFull_Stripper);
        for(int i = 0; i <= 4; i++) M2_Full += ((Stripper::toDouble(me()[i]))*std::pow(std::log(mu*mu/COM/COM), i))*average_factor_full;
        //for(int i = 4; i <= 4; i++) M2_Full += (Stripper::toDouble(me()[i]))*average_factor_full*24.;
      }
    }
    else {
      M2_Full = M2_custom[counter];
      counter++;
    }

    double M2_approx, M2_test;
    if(limit == "soft") {
      if(unresolved_str == " g") {
        if(order=="LO") M2_approx = soft_g_squared(pp_full, M0_ij, A)*average_factor_full/average_factor;
        else if(order=="NLO") M2_approx = soft_g_squared_1l(pp_full, M0_ij, M0_ijk, M1_ij, A)*average_factor_full/average_factor;
        else if(order=="NNLO") {
          //M2_test = soft_g_squared_1l(pp_full, M1_ij, M1_ijk, M2_ij, A)*average_factor_full/average_factor;
          M2_approx = soft_g_squared_2l(pp_full, M0_ij, M0_ijk, Q_ijkl, M1_ij, M1_ijk, M2_ij, A, n_f)*average_factor_full/average_factor; // *(-0.5)?
          M2_test = soft_g_squared_2l_reducible(pp_full, M0_ij, M0_ijk, Q_ijkl, M1_ij, M1_ijk, M2_ij, A, n_f)*average_factor_full/average_factor;
        }
      }
      else if(unresolved_str == " d d~") {
        if(order == "LO") M2_approx = soft_qq_squared(pp_full, M0_ij, A)*average_factor_full/average_factor;
        else if(order == "NLO") M2_approx = soft_qq_squared_1l(pp_full, M0_ij, M1_ij, M0_ijk, dM0_ijk, A, n_f);
      }
      else if(unresolved_str == " g g") {
        if(order=="LO") M2_approx = soft_gg_squared(pp_full, M0_ij, M0_ijkl, A)*average_factor_full/average_factor;
        else if(order=="NLO") M2_approx = soft_gg_squared_1l(pp_full, M0_ij, M1_ij, M0_ijkl, M1_ijkl, M0_ijk, M0_ijkla, Q_ijkl, A, n_f)*average_factor_full/average_factor;
      }
      else if(unresolved_str == " g d d~")
        M2_approx = soft_gqq_squared(pp_full, M0_ij, M0_ijkl, dM0_ijk, A)*average_factor_full/average_factor;
      //else if(unresolved_str == " g g g")
      //  M2_approx = soft_ggg_squared(pp_full, M_ij, M0_ijkl, M_ijklab, dM0_ijk, fM_ijkl)*average_factor_full/average_factor;

    }
    else if(limit == "collinear") {
      int i_reference = clusterTree.getRoot()->children[0]->data.reference;
      switch (nUnresolved)
      {
      case 1:
        if(order=="LO") M2_approx = collinear_squared(pp_full, SC0, A, A_full, i_reference)*average_factor_full/average_factor;
        else if(order=="NLO") {
          M2_test = collinear_squared(pp_full, SC1, A, A_full, i_reference)*average_factor_full/average_factor;
          M2_approx = collinear_squared_1l(pp_full, SC0, SC1, A, A_full, i_reference)*average_factor_full/average_factor;
        }
        else if(order=="NNLO") {
          M2_test = collinear_squared(pp_full, SC2, A, A_full, i_reference)*average_factor_full/average_factor;
          M2_test += collinear_squared_1l(pp_full, SC1, SC2, A, A_full, i_reference)*average_factor_full/average_factor;
          M2_approx = collinear_squared_2l(pp_full, SC0, SC1, SC2, A, A_full, i_reference)*average_factor_full/average_factor;
        }
        break;
      case 2:
        if(order=="LO") M2_approx = triple_collinear_squared(pp_full, SC0, A, A_full, i_reference)*average_factor_full/average_factor;
        else if(order=="NLO") {
          M2_test = triple_collinear_squared(pp_full, SC1, A, A_full, i_reference)*average_factor_full/average_factor;
          M2_approx = triple_collinear_squared_1l(pp_full, SC0, SC1, A, A_full, i_reference)*average_factor_full/average_factor;
        }
        break;
      case 3:
        if(order=="LO") M2_approx = quadruple_collinear_squared(pp_full, SC0, A, A_full, i_reference)*average_factor_full/average_factor;
        break;
      default:
        break;
      }
    }
    std::cout << scale << "\t" << M2_Full << "\t" << M2_approx << "\t" << M2_test << "\t" << M2_approx + M2_test << "\t" << (M2_Full - M2_test)/M2_approx << "\t" << std::abs(1. - M2_Full/(M2_approx + M2_test)) << std::endl;
    outfile1 << scale << ", " << std::abs(1. - M2_Full/(M2_approx + M2_test)) << std::endl;
  }
  }
  outfile1.close();
}