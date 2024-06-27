#include "main.hpp"

// Evaluation of input
std::string process_str = "d d~ -> g g";
std::string unresolved_str = " g";
std::string process_full_str = process_str + unresolved_str;


const int nBorn = 4;
const int nUnresolved = 1;
const int power = 2;
const std::vector<int> flavor = {1,1,1,1};

int main() {
  std::cout.precision(8);
  // Initialize txt outstream (delete previous file)

  std::ofstream outfile1;
  outfile1.open ("results/" + process_full_str + ".txt");
  outfile1.close();
  outfile1.open ("results/" + process_full_str + ".txt");

  amplitude A(process_str);
  amplitude A_full(process_full_str);

  // Recola Settings
  Recola::set_reduction_mode_rcl(4);
  Recola::set_print_level_amplitude_rcl(2);
  Recola::set_alphas_rcl(gs * gs/(4 * M_PI), mu, 5);
  Recola::set_mu_ir_rcl(mu);
  Recola::set_delta_ir_rcl(0, M_PI * M_PI/12.);
  Recola::set_momenta_correction_rcl(false);

  Recola::use_alpha0_scheme_rcl(e*e/4./M_PI);
  Recola::set_pole_mass_z_rcl(1.e8, 0.0001);
  Recola::set_pole_mass_w_rcl(1.e15, 0.0001);

  // Define & generate process
  Recola::define_process_rcl(1, process_str, "LO");
  Recola::define_process_rcl(2, process_full_str, "LO");

  Recola::generate_processes_rcl();
  // Define initial state
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



  // Compute process amplitudes
  Recola::compute_process_rcl(1, pp_rcl, std::string("LO"));
  Recola::compute_all_colour_correlations_rcl(1, pp_rcl);
  Recola::compute_process_rcl(2, ppFull_rcl, std::string("LO"));

  // Get non-vanishing helicity and color-configurations
  std::vector<std::vector<int> > helicities = non0hel(1, power, "LO", A);
  std::cout << "Found non zero helicity configurations of Born process" << std::endl;
  std::vector<std::vector<int> > colors = non0col(1, power, "LO", A, helicities);
  std::cout << "Found non zero color configurations of Born process" << std::endl;
  std::vector<std::vector<int> > helicities_full = non0hel(2, power + nUnresolved, "LO", A_full);
  std::cout << "Found non zero helicity configurations of Full process" << std::endl;
  std::vector<std::vector<int> > colors_full = non0col(2, power + nUnresolved, "LO", A_full, helicities_full);
  std::cout << "Found non zero color configurations of Full process" << std::endl;
  // Create Hashmap for all non-zero amplitudes that only need to be caluculated once per data set
  std::unordered_map<std::string, std::complex<double>> M; // hel + col is the key as string
  std::unordered_map<std::string, std::complex<double>> M_ij, fM_ijk, dM_ijk; // color correlators

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
      M[key] = colorflow2color(hel, col, A, power, "LO", 1);
      std::cout << key << "\t" << M[key] << std::endl;
    }
  }
  double average_factor = 1./4.; // Spin of the initial state
  double average_factor_full = 1./4.;
  for(int i = 0; i < 2; i++) {
    average_factor *= 1./double(A.process[i]);
    average_factor_full *= 1./double(A_full.process[i]);
  }
  std::vector<std::string> particles = A.process_particles;
  std::vector<std::string> particles_full = A_full.process_particles;
  particles.erase(particles.begin(), particles.begin() + 2);
  particles_full.erase(particles_full.begin(), particles_full.begin() + 2);

  while(particles.size() > 0) {
    std::string particle = particles[0];
    int occurances = 1;
    int i = 1;
    while(i < particles.size()) {
      std::cout << particles[i] << std::endl;
      if(particles[i] == particle) {
        occurances += 1;
        particles.erase(particles.begin(), particles.begin() + i);
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
      std::cout << particles_full[i] << std::endl;
      if(particles_full[i] == particle) {
        occurances += 1;
        particles_full.erase(particles_full.begin(), particles_full.begin() + i);
      }
      else {
        i += 1;
      }
    }
    particles_full.erase(particles_full.begin());
    average_factor_full *= 1./factorial(occurances);
  }
  for (int i = 0; i < A.process.size(); i++) {
    if(A.particle_type[i] == 0) continue;
    for(int j = 0; j < A.process.size(); j++) {
      if(A.particle_type[j] == 0) continue;
      double control;
      double M2_averaged_control;
      double M2_averaged = 0;
      Recola::get_colour_correlation_rcl(1, power, i + 1, j + 1, control);
      if(A.process[i] == 3) control*= 4./3.;
      else if(A.particle_type[i] == 8) control *= 3.;
      Recola::get_squared_amplitude_rcl(1, power, "LO", M2_averaged_control);
      std::complex<double> Mij = 0;
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
          std::complex<double> M_bra = std::conj(M[key]);
          M2_averaged += std::pow(std::abs(M[key]), 2)*average_factor;
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
                Mij += col_factor*M_bra*M[key];
              }
            }
          }
        }
      }
      M_ij[std::to_string(i) + std::to_string(j)] = Mij*average_factor;

      std::cout << "i = " << i << ", j = " << j << ": <M|T_i.T_j|M> = " << Mij*average_factor << "\t" << control << "\t" << Mij*average_factor/control << std::endl;
      std::cout << "M2 = " << M2_averaged << "\tcontrol = " << M2_averaged_control << "\tratio = " << M2_averaged/M2_averaged_control << std::endl;
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
  for(int tree_counter = 0; tree_counter < 1; tree_counter++) {
  Tree<Cluster>& tree = trees[tree_counter];
  tree.print();
  if(tree.getRoot()->children.size() > 1) continue; // Remove multiple reference, as they are trivial
  std::vector<Tree<Cluster>> sectors = GenSectors(flavor, tree, nBorn);

  for(int sec_counter = 2; sec_counter < 3; sec_counter++) {
  Tree<Cluster> clusterTree = sectors[sec_counter];
  std::cout << "reference = " << clusterTree.getRoot()->children[0]->data.reference << std::endl;

  // Generate phase-space points
  double scale = 1;
  double increment = 0.1;
  while (scale > 1.e-7) {
    scale *= increment;
    std::vector<std::vector<std::vector<double>>> xParFull;
    int level_int = 1;
    std::vector<TreeNode<Cluster>*> level = clusterTree.getLevel(level_int);
    while(level.size() > 0) {
      int unresolved_level = 0;
      for(TreeNode<Cluster>* node : level){
        unresolved_level += node->data.unresolved;

        std::vector<std::vector<double>> xPar;
        for(int c = 0; c < node->data.unresolved; c++) {
          double xi = std::pow(scale, 1);
          double eta = 0.5;
          double phi = rnd(0., 1.);
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
    std::vector<std::unordered_map<std::string, std::complex<double>>> J_g, J_gg, J_qq, J_qqg;

    if((unresolved_str == " g") or (unresolved_str == " g g")) {
      for(int i = 0; i < nUnresolved; i++) J_g.push_back(J_g_eikonal(pp_arr, ppFull_rcl[A.process.size() + i], A));
      for(int i = 0; i < nUnresolved - 1; i++)
        for(int j = i + 1; j < nUnresolved; j++)
          J_gg.push_back(J_gg_eikonal(pp_arr, ppFull_rcl[A.process.size() + i], ppFull_rcl[A.process.size() + j], A));
    }
    else if ((unresolved_str == " d d~")) {
      for(int i = 0; i < nUnresolved - 1; i++)
        for(int j = i + 1; j < nUnresolved; j++)
          J_qq.push_back(J_qq_eikonal(pp_arr, ppFull_rcl[A.process.size() + i], ppFull_rcl[A.process.size() + j], A));
    }
    else if ((unresolved_str == " d d~ g")) {
      J_g.push_back(J_g_eikonal(pp_arr, ppFull_rcl[A.process.size() + 2], A));
      J_qq.push_back(J_qq_eikonal(pp_arr, ppFull_rcl[A.process.size()], ppFull_rcl[A.process.size() + 1], A));
      J_qqg.push_back(J_qqg_eikonal(pp_arr, ppFull_rcl[A.process.size()], ppFull_rcl[A.process.size() + 1], ppFull_rcl[A.process.size() + 2], A));
    }

    Recola::compute_process_rcl(2, ppFull_rcl, std::string("LO"));
    double M2_Full;
    Recola::get_squared_amplitude_rcl(2, power + nUnresolved, "LO", M2_Full);
    double M2_approx = soft_g_squared(pp_full, M_ij, A)*average_factor_full/average_factor;
    std::cout << M2_Full << "\t" << M2_approx << "\t" << std::abs(1. - M2_Full/M2_approx) << "\t" << M2_Full/M2_approx << "\t" << scale << std::endl;
    double diff = 0.;
    int counter = 0;
    for (int i = 0; i < helicities_full.size(); i++) {
      // Create hel variables
      int *hel_full = &helicities_full[i][0];
      for (int j = 0; j < colors_full.size(); j++) {
        // Create col variables
        int *col_full = &colors_full[j][0];
        // Compute amplitudes
        // Born amplitude
        //std::complex<double> M_eikonal = tree_lp(pp_full, col_full, hel_full, A, M, 0);
        std::complex<double> M_eikonal;
        if(unresolved_str == " g")
          M_eikonal = soft_g(pp_full, J_g[0],  M, helicities_full[i], colors_full[j], A);
        else if(unresolved_str == " g g") {
          M_eikonal = soft_gg(pp_full, J_g[0], J_g[1], J_gg[0],  M, helicities_full[i], colors_full[j], A);
          //M_eikonal = double_soft_tree(pp_full, J_g[0], J_g[1], M, helicities_full[i],  colors_full[j], A);
        }
        else if(unresolved_str == " d d~") {
          //M_eikonal = double_soft_tree_qq(pp_full, M, helicities_full[i], colors_full[j], A);
          M_eikonal = soft_qq(pp_full, J_qq[0], M, helicities_full[i], colors_full[j], A);
        }
        else if (unresolved_str == " d d~ g") {
          M_eikonal = soft_qqg(pp_full, J_qq[0], J_g[0], J_qqg[0], M, helicities_full[i], colors_full[j], A);
        }
        std::complex<double> M_full = colorflow2color(hel_full, col_full, A_full, power + nUnresolved, "LO", 2);
        std::string hel_string;
        std::string col_string;
        for (int i = 0; i < A_full.process.size(); i++) {
            hel_string += std::to_string(hel_full[i]);
            col_string += std::to_string(col_full[i]);
        }
        std::string key = hel_string + col_string;
        if(std::abs(M_full)/std::sqrt(M2_Full) > scale) {
          diff += std::abs((M_full - M_eikonal)/M_full);
        }

        std::cout << key << ":\t" << std::abs(M_full) << "\t" << std::abs(M_eikonal) << "\t" << std::abs(M_eikonal/M_full) << "\t" << std::abs((M_eikonal - M_full)/M_full) << "\t" << scale << std::endl;

        counter++;
      }
    }
    outfile1 << scale << ", " << diff/double(counter) << std::endl;
  }
  }}
  outfile1.close();
}