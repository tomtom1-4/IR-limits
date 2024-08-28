#include "main.hpp"

// Evaluation of input
std::string process_str = "u u~ -> u u~ e- e+";
std::string unresolved_str = " d d~";
std::string process_full_str = process_str + unresolved_str;


const int nBorn = 6;
const int power = 2;
const std::vector<int> flavor = {0,0,1,1,0,0};
const std::string order = "NLO";
const std::string suffix = "_EW"; // "" or "_EW"
double M2_custom[13];

int main() {
  bool custom = false;
  if((process_full_str == "u u~ -> u u~ e- e+ g g") and (order == "NLO")) {
    custom = true;
    M2_custom[0] = -4.58613727787852e-26;
    M2_custom[1] = -3.05566873598841e-22;
    M2_custom[2] = -6.78415558513149e-20;
    M2_custom[3] = -1.2213501061979e-17;
    M2_custom[4] = -1.93961648681395e-15;
    M2_custom[5] = -2.82993140226372e-13;
    M2_custom[6] = -3.88880026114949e-11;
    M2_custom[7] = -5.11393342703531e-09;
    M2_custom[8] = -6.5040578409581e-07;
    M2_custom[9] = -8.05858355886322e-05;
    M2_custom[10] = -0.0097845849231643;
    M2_custom[11] = -1.16517247862666;
    M2_custom[12] = -136.406428204994;
  }
  else if((process_full_str == "e- e+ -> u u~ u u~ g g") and (order == "NLO")) {
    custom = true;
    M2_custom[0] = 8.42714370716502e-25;
    M2_custom[1] = -1.07948493039382e-22;
    M2_custom[2] = -3.1589420159308e-20;
    M2_custom[3] = -1.25043991989418e-17;
    M2_custom[4] = -1.01797400996888e-15;
    M2_custom[5] = -1.50917187473323e-13;
    M2_custom[6] = -2.09150310626981e-11;
    M2_custom[7] = -2.76442260780936e-09;
    M2_custom[8] = -3.52770343464854e-07;
    M2_custom[9] = -4.38128658522702e-05;
    M2_custom[10] = -0.00532339243608503;
    M2_custom[11] = -0.632583283928193;
    M2_custom[12] = -72.7129589447317;
  }

  int delta_power = 0;
  if (order == "NLO") delta_power = 2;
  std::cout.precision(8);
  // Initialize txt outstream (delete previous file)

  std::ofstream outfile1;
  outfile1.open ("results/" + process_str + " +" + unresolved_str + "_" + order + suffix + ".txt");
  outfile1.close();
  outfile1.open ("results/" + process_str + " +" + unresolved_str + "_" + order + suffix + ".txt");

  amplitude A(process_str);
  amplitude A_full(process_full_str);
  const int nUnresolved = A_full.process.size() - A.process.size();

  if(!custom) {
    // Recola Settings
    Recola::set_reduction_mode_rcl(4);
    Recola::set_print_level_amplitude_rcl(2);
    Recola::set_alphas_rcl(gs * gs/(4 * M_PI), mu, 5);
    Recola::set_mu_ir_rcl(mu);
    Recola::set_delta_ir_rcl(0, M_PI * M_PI/12.);
    Recola::set_momenta_correction_rcl(false);

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
  PhaseSpace dumn = Splitting(nBorn - 2, COM);
  PhaseSpace pp = Splitting(nBorn - 2, COM);
  PhaseSpace ppFull_dummy = Splitting(nBorn + nUnresolved - 2, COM);
  pp.print();
  // Transform to Recola format
  double pp_rcl[nBorn][4], ppFull_rcl[nBorn + nUnresolved][4], pp_arr[4*nBorn];
  for(int i = 0; i < nBorn; i++) {
    for(int j = 0; j < 4; j++) {
      pp_rcl[i][j] = (i<2?-1.:1.)*pp.momenta[i].components[j];
      pp_arr[i*4+j] = (i<2?-1.:1.)*pp.momenta[i].components[j];
    }
  }
  for(int i = 0; i < nBorn + nUnresolved; i++) {
    for(int j = 0; j < 4; j++) {
      ppFull_rcl[i][j] = (i<2?-1.:1.)*ppFull_dummy.momenta[i].components[j];
    }
  }


  if(!custom) {
    // Compute process amplitudes
    Recola::compute_process_rcl(1, pp_rcl, order);
    //Recola::compute_all_colour_correlations_rcl(1, pp_rcl);
    Recola::compute_process_rcl(2, ppFull_rcl, order);
  }

  // Get non-vanishing helicity and color-configurations
  std::unordered_map<std::string, double> M0_ij, M1_ij, M0_ijk, dM0_ijk, M0_ijkl, M1_ijkl, M0_ijklab, Q_ijkl, M0_ijkla; // color correlators

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
  std::cout << average_factor << std::endl;
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
  std::cout << "average_factor_full/average_factor = " << average_factor_full/average_factor << std::endl;
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
  // Get color correlators
  // <M|Ti.Tj|M>
  XMLElement * M0_ijElement = pRoot->FirstChildElement("M0_ij");
  XMLElement * M1_ijElement = pRoot->FirstChildElement("M1_ij");
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
  if((order=="NLO" && unresolved_str==" g g") || unresolved_str==" g g g") {
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
  if(unresolved_str==" g g g" or order=="NLO") {
    XMLElement * M0_ijkElement = pRoot->FirstChildElement("M0_ijk");
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

  std::cout << "Filled the Hashmaps" << std::endl;

  // Define ClusterTree
  std::vector<Tree<Cluster>> trees = GenTrees(nUnresolved);
  for(int tree_counter = 0; tree_counter < trees.size(); tree_counter++) {
    Tree<Cluster> clusterTree = trees[tree_counter];
    std::cout << "tree_counter = " << tree_counter << std:: endl;
    clusterTree.print();
    std::cout << "#########################################################################" << std::endl;
  }
  Tree<Cluster> tree;
  if(nUnresolved == 1) tree = trees[0];
  else if(nUnresolved == 2) tree = trees[1];
  tree.print();
  std::vector<Tree<Cluster>> sectors = GenSectors(flavor, tree, nBorn);

  for(int sec_counter = 0; sec_counter < 1; sec_counter++) {
  Tree<Cluster> clusterTree = sectors[sec_counter];
  std::cout << "reference = " << clusterTree.getRoot()->children[0]->data.reference << std::endl;

  // Generate phase-space points
  double scale = 1;
  std::vector<double> etas(nUnresolved);
  std::vector<double> phis(nUnresolved);
  for(int i = 0; i < nUnresolved; i++) {
    etas[i] = rnd(0.1, 0.9);
    phis[i] = rnd(0., 1.);
  }
  double increment = std::sqrt(0.1);
  int counter = 0;
  while (scale > 1.e-8) {
    scale *= increment;
    std::vector<std::vector<std::vector<double>>> xParFull;
    int level_int = 1;
    std::vector<TreeNode<Cluster>*> level = clusterTree.getLevel(level_int);
    int unresolved_counter = 0;
    while(level.size() > 0) {
      int unresolved_level = 0;
      for(TreeNode<Cluster>* node : level){
        unresolved_level += node->data.unresolved;

        std::vector<std::vector<double>> xPar;
        for(int c = 0; c < node->data.unresolved; c++) {
          double xi = std::pow(scale, 1);
          double eta = etas[unresolved_counter];
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
    PhaseSpace ppFull = GenMomenta2(pp, clusterTree, xParFull);
    ppFull.print();
    double pp_full[4*(nBorn+nUnresolved)];
    for(int i = 0; i < ppFull.momenta.size(); i++) {
      for(int j = 0; j < 4; j++) {
        pp_full[4*i+j] = (i<2?-1.:1.)*ppFull.momenta[i].components[j];
        ppFull_rcl[i][j] = (i<2?-1.:1.)*ppFull.momenta[i].components[j];
      }
    }

    double M2_Full;
    if(!custom) {
      Recola::compute_process_rcl(2, ppFull_rcl, order);
      Recola::get_squared_amplitude_rcl(2, power + nUnresolved + delta_power/2, order, M2_Full);
    }
    else {
      M2_Full = M2_custom[counter];
      counter++;
    }

    double M2_approx;
    if(unresolved_str == " g") {
      if(order=="LO") M2_approx = soft_g_squared(pp_full, M0_ij, A)*average_factor_full/average_factor;
      else if(order=="NLO") M2_approx = soft_g_squared_1l(pp_full, M0_ij, M0_ijk, M1_ij, A)*average_factor_full/average_factor;
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
    double M2_test = soft_qq_squared_1l_oneLoop(pp_full, M0_ij, M1_ij, M0_ijk, dM0_ijk, A, n_f);
    std::cout << scale << "\t" << M2_Full << "\t" << M2_approx  << "\t" << (M2_Full - M2_approx)/M2_test << "\t" << std::abs(1. - M2_Full/M2_approx) << std::endl;
    outfile1 << scale << ", " << std::abs(1. - M2_Full/M2_approx) << std::endl;
  }
  }
  outfile1.close();
}