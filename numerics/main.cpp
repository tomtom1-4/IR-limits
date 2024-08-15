#include "main.hpp"

// Evaluation of input
std::string process_str = "e- e+ -> u u~";
std::string unresolved_str = " g";
std::string process_full_str = process_str + unresolved_str;


const int nBorn = 4;
const int nUnresolved = 1;
const int power = 0;
const std::vector<int> flavor = {0,0,1,1};
const std::string order = "NLO";

int main() {
  int delta_power = 0;
  if (order == "NLO") delta_power = 2;
  std::cout.precision(8);
  // Initialize txt outstream (delete previous file)

  std::ofstream outfile1;
  outfile1.open ("results/" + process_full_str + "_" + order + ".txt");
  outfile1.close();
  outfile1.open ("results/" + process_full_str + "_" + order + ".txt");

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
  //Recola::set_pole_mass_z_rcl(1.e8, 0.0001);
  //Recola::set_pole_mass_w_rcl(1.e15, 0.0001);

  // Define & generate process
  Recola::define_process_rcl(1, process_str, order);
  Recola::define_process_rcl(2, process_full_str, order);
  Recola::set_otter_mode_rcl(1, "oneloop_qp");
  Recola::set_otter_mode_rcl(2, "oneloop_qp");

  Recola::generate_processes_rcl();
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



  // Compute process amplitudes
  Recola::compute_process_rcl(1, pp_rcl, order);
  //Recola::compute_all_colour_correlations_rcl(1, pp_rcl);
  Recola::compute_process_rcl(2, ppFull_rcl, order);

  // Get non-vanishing helicity and color-configurations
  std::vector<std::vector<int> > helicities = non0hel(1, power + delta_power, order, A);
  std::cout << "Found non zero helicity configurations of Born process" << std::endl;
  std::vector<std::vector<int> > colors = non0col(1, power + delta_power, order, A, helicities);
  std::cout << "Found non zero color configurations of Born process" << std::endl;
  /*std::vector<std::vector<int> > helicities_full = non0hel(2, power + nUnresolved, "LO", A_full);
  std::cout << "Found non zero helicity configurations of Full process" << std::endl;
  std::vector<std::vector<int> > colors_full = non0col(2, power + nUnresolved, "LO", A_full, helicities_full);*/
  std::cout << "Found non zero color configurations of Full process" << std::endl;
  // Create Hashmap for all non-zero amplitudes that only need to be caluculated once per data set
  std::unordered_map<std::string, std::complex<double>> M0, M1; // hel + col is the key as string
  std::unordered_map<std::string, std::complex<double>> M0_ij, M1_ij, fM_ijk, dM_ijk, M_ijkl, M_ijklab, Q_ijkl; // color correlators
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
      if(order=="NLO") M1[key] = colorflow2color(hel, col, A, power + 2, "NLO", 1);
      if(std::abs(M0[key]) + std::abs(M1[key]) > 1.e-17) {
        std::vector<int> dummy = hel_full;
        dummy.insert(dummy.end(), col_full.begin(), col_full.end());
        keys.push_back(dummy);
        std::cout << key << "\t" << M0[key];
        if(order=="NLO") std::cout << "\t" << M1[key];
        std::cout << std::endl;
      }
    }
  }
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
  // Get color correlators
  // <M|Ti.Tj|M>
  for (int i = 0; i < A.process.size(); i++) {
    if(A.particle_type[i] == 0) continue;
    for(int j = 0; j < A.process.size(); j++) {
      if(A.particle_type[j] == 0) continue;
      if(std::abs(M0_ij[std::to_string(j) + std::to_string(i)]) + std::abs(M1_ij[std::to_string(j) + std::to_string(i)]) != 0.) {
        M0_ij[std::to_string(i) + std::to_string(j)] = M0_ij[std::to_string(j) + std::to_string(i)];
        M1_ij[std::to_string(i) + std::to_string(j)] = M1_ij[std::to_string(j) + std::to_string(i)];
        continue;
      }
      double control;
      double M2_averaged_control;
      double M2_averaged = 0;
      //Recola::get_colour_correlation_rcl(1, power, i + 1, j + 1, control);

      Recola::get_squared_amplitude_rcl(1, power, "LO", M2_averaged_control);
      //Recola::get_squared_amplitude_rcl(1, power + delta_power/2, order, M2_averaged_control);
      std::cout << "M2_averaged_control = " << M2_averaged_control << std::endl;
      std::complex<double> M0ij = 0;
      std::complex<double> M1ij = 0;
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
        if(order=="LO") M2_averaged += std::pow(std::abs(M0[key]), 2)*average_factor;
        if(order=="NLO") M2_averaged += 2.*std::real(M0[key]*std::conj(M1[key]))*average_factor;
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
              M0ij += col_factor*M0_bra*M0[key];
              M1ij += 2*std::real(col_factor*M1_bra*M0[key]);
            }
          }
        }
      }
      M1_ij[std::to_string(i) + std::to_string(j)] = M1ij*average_factor;
      M0_ij[std::to_string(i) + std::to_string(j)] = M0ij*average_factor;
      std::cout << "M2Control = " << M2_averaged_control << "\t" << M2_averaged << std::endl;
      std::cout << "i = " << i << ", j = " << j << ": <M|T_i.T_j|M> = " << M0ij*average_factor << "\t" << control;
      if(order == "NLO") std::cout << "\t" << M1ij*average_factor;
      std::cout << std::endl;
    }
  }
  // <M|Ti.Tj Tk.Tl|M>
  if(unresolved_str==" g g" || unresolved_str==" g g g") {
  for (int i = 0; i < A.process.size(); i++) {
    if(A.particle_type[i] == 0) continue;
    for(int j = 0; j < A.process.size(); j++) {
      if(A.particle_type[j] == 0) continue;
      for(int k = 0; k < A.process.size(); k++) {
        if(A.particle_type[k] == 0) continue;
        for(int l = 0; l < A.process.size(); l++) {
          if(A.particle_type[l] == 0) continue;
          std::complex<double> Mijkl = 0.;
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

                  Mijkl += col_factor*M_bra*M0[key];
                }
              }
            }
          }

          M_ijkl[std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l)] = Mijkl*average_factor;
          double M2_averaged;
          Recola::get_squared_amplitude_rcl(1, power, "LO", M2_averaged);
          std::cout << "i = " << i << ", j = " << j << ", k = " << k << ", l = " << l << ": <M|T_i.T_j T_k.T_l|M> = " << Mijkl*average_factor << "\t" << Mijkl*average_factor/M2_averaged << std::endl;
        }
      }
    }
  }
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

              M_ijkl[std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l) + std::to_string(a) + std::to_string(b)] = Mijklab*average_factor;
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

  // d^{a,b,c} <M|Ti^a Tj^b Tk^c|M>
  // f^{a,b,c} <M|Ti^a Tj^b Tk^c|M>
  if(unresolved_str==" g g g" or order=="NLO") {
  std::vector<std::vector<int>> configurations;
  for (int i = 0; i < A.process.size(); i++) {
    if(A.particle_type[i] == 0) continue;
    for(int j = 0; j < A.process.size(); j++) {
      if(A.particle_type[j] == 0) continue;
      if(j == i) continue;
      for(int k = 0; k < A.process.size(); k++) {
        if(A.particle_type[k] == 0) continue;
        if((k == i) or (k == j)) continue;
        std::complex<double> fMijk = 0.;
        std::complex<double> dMijk = 0.;
        bool permutation = false;
        for(std::vector<int> configuration : configurations) {
          if(configuration==std::vector<int>({i,j,k})) {
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
          std::complex<double> M_bra = std::conj(M0[key]);
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

            if(std::abs(M0[key]) < 1.e-17) continue;
            for(int a = 0; a < 8; a++) for(int b = 0; b < 8; b++) for(int c = 0; c < 8; c++) {
              std::complex<double> col_factor = 1;

              if(A.particle_type[k] == 1)
                col_factor *= lam[a][3*output_colors[2] + ck]/2.;
              else if(A.particle_type[k] == -1)
                col_factor *= -lam[a][3*ck + output_colors[2]]/2.;
              else if(A.particle_type[k] == 2)
                col_factor *= I*fabc[a][8*ck + output_colors[2]];
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
                col_factor *= lam[c][3*output_colors[0] + ci]/2.;
              else if(A.particle_type[i] == -1)
                col_factor *= -lam[c][3*ci + output_colors[0]]/2.;
              else if(A.particle_type[i] == 2)
                col_factor *= I*fabc[c][8*ci + output_colors[0]];
              else
                std::cout << "unknown particle type" << std::endl;
              fMijk += col_factor*M_bra*M0[key]*fabc[a][8*b+c];
              dMijk += col_factor*M_bra*M0[key]*dabc[a][8*b+c];
            }
          }
        }
        fM_ijk[std::to_string(i) + std::to_string(j) + std::to_string(k)] = fMijk*average_factor;
        // permutations
        fM_ijk[std::to_string(i) + std::to_string(k) + std::to_string(j)] = -fMijk*average_factor;
        fM_ijk[std::to_string(j) + std::to_string(i) + std::to_string(k)] = -fMijk*average_factor;
        fM_ijk[std::to_string(j) + std::to_string(k) + std::to_string(i)] = fMijk*average_factor;
        fM_ijk[std::to_string(k) + std::to_string(j) + std::to_string(i)] = -fMijk*average_factor;
        fM_ijk[std::to_string(k) + std::to_string(i) + std::to_string(j)] = fMijk*average_factor;
        configurations.push_back(std::vector<int>({i,j,k}));
        configurations.push_back(std::vector<int>({i,k,j}));
        configurations.push_back(std::vector<int>({j,i,k}));
        configurations.push_back(std::vector<int>({j,k,i}));
        configurations.push_back(std::vector<int>({k,j,i}));
        configurations.push_back(std::vector<int>({k,i,j}));

        dM_ijk[std::to_string(i) + std::to_string(j) + std::to_string(k)] = dMijk*average_factor;

        double M2_averaged;
        Recola::get_squared_amplitude_rcl(1, power, "LO", M2_averaged);
        std::string lock = std::to_string(i) + std::to_string(j) + std::to_string(k);
        std::cout << lock << ": f^{abc} <M|Ti Tj Tk|M> = " << fM_ijk[lock] << "\t" << fMijk*average_factor << "\t" << fM_ijk[lock]/M2_averaged << std::endl;
        //std::cout << "i = " << i << ", j = " << j << ", k = " << k << ": d^{abc} <M|Ti Tj Tk|M> = " << dMijk*average_factor << "\t" << dMijk*average_factor/M2_averaged << std::endl;

      }
    }
  }
  }

  // f^{a,b;c,d} <M|Ti^a {Tj^b, Tk^c} Tl^d|M> + h.c.
  if((order=="NLO" && unresolved_str==" g g") || unresolved_str==" g g g") {
  for (int i = 0; i < A.process.size(); i++) {
    if(A.particle_type[i] == 0) continue;
    for(int j = 0; j < A.process.size(); j++) {
      if(A.particle_type[j] == 0) continue;
      for(int k = 0; k < A.process.size(); k++) {
        if(A.particle_type[k] == 0) continue;
        for(int l = 0; l < A.process.size(); l++) {
          if(A.particle_type[l] == 0) continue;
          std::complex<double> ffMijkl = 0.;

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
              std::vector<int> output_colors = {col_full[i], col_full[j], col_full[k],  col_full[l]};
              for(int a = 0; a < 8; a++) for(int b = 0; b < 8; b++) for(int c = 0; c < 8; c++) for(int d = 0; d < 8; d++) {
                for(int ci = 0; ci < A.process[i]; ci++) for(int cj = 0; cj < A.process[j]; cj++) for(int ck = 0; ck < A.process[k]; ck++) for(int cl = 0; cl < A.process[l]; cl++) {
                  std::string key2 = key;

                  key.replace(key.size() - A.process.size() + i, 1, std::to_string(ci));
                  key.replace(key.size() - A.process.size() + j, 1, std::to_string(cj));
                  key.replace(key.size() - A.process.size() + k, 1, std::to_string(ck));
                  key.replace(key.size() - A.process.size() + l, 1, std::to_string(cl));

                  key2.replace(key2.size() - A.process.size() + i, 1, std::to_string(ci));
                  key2.replace(key2.size() - A.process.size() + k, 1, std::to_string(ck));
                  key2.replace(key2.size() - A.process.size() + j, 1, std::to_string(cj));
                  key2.replace(key2.size() - A.process.size() + l, 1, std::to_string(cl));

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
                    col_factor2 *= lam[b][3*output_colors[3] + cl]/2.;
                  else if(A.particle_type[l] == -1)
                    col_factor2 *= -lam[b][3*cl + output_colors[3]]/2.;
                  else if(A.particle_type[l] == 2)
                    col_factor2 *= I*fabc[b][8*cl + output_colors[3]];

                  if(A.particle_type[k] == 1)
                    col_factor2 *= lam[d][3*output_colors[2] + ck]/2.;
                  else if(A.particle_type[k] == -1)
                    col_factor2 *= -lam[d][3*ck + output_colors[2]]/2.;
                  else if(A.particle_type[k] == 2)
                    col_factor2 *= I*fabc[d][8*ck + output_colors[2]];

                  if(A.particle_type[j] == 1)
                    col_factor2 *= lam[c][3*output_colors[1] + cj]/2.;
                  else if(A.particle_type[j] == -1)
                    col_factor2 *= -lam[c][3*cj + output_colors[1]]/2.;
                  else if(A.particle_type[j] == 2)
                    col_factor2 *= I*fabc[c][8*cj + output_colors[1]];

                  if(A.particle_type[i] == 1)
                    col_factor2 *= lam[a][3*output_colors[0] + ci]/2.;
                  else if(A.particle_type[i] == -1)
                    col_factor2 *= -lam[a][3*ci + output_colors[0]]/2.;
                  else if(A.particle_type[i] == 2)
                    col_factor2 *= I*fabc[a][8*ci + output_colors[0]];

                  for(int s = 0; s < 8; s++)
                    ffMijkl += 2*std::real(M_bra*(col_factor1*M0[key] + col_factor2*M0[key2]))*fabc[a][8*b + s]*fabc[s][8*c + d];
                }
              }
            }
            Q_ijkl[std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l)] = ffMijkl*average_factor;

            double M2_averaged;
            Recola::get_squared_amplitude_rcl(1, power, "LO", M2_averaged);
            std::cout << "i = " << i << ", j = " << j << ", k = " << k << ", l = " << l << ": f^{abc} f^{cde} <M|Ti^a {Tj^b, Tk^d} Tl^e|M> = " << ffMijkl*average_factor << "\t" << ffMijkl*average_factor/M2_averaged << std::endl;

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
  Tree<Cluster>& tree = trees[0];
  tree.print();
  std::vector<Tree<Cluster>> sectors = GenSectors(flavor, tree, nBorn);

  for(int sec_counter = 0; sec_counter < 1; sec_counter++) {
  Tree<Cluster> clusterTree = sectors[sec_counter];
  std::cout << "reference = " << clusterTree.getRoot()->children[0]->data.reference << std::endl;

  // Generate phase-space points
  double scale = 1;
  double increment = 0.1;//std::sqrt(0.1);
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

    Recola::compute_process_rcl(2, ppFull_rcl, order);
    double M2_Full;
    Recola::get_squared_amplitude_rcl(2, power + nUnresolved + delta_power/2, order, M2_Full);

    double M2_approx;
    if(unresolved_str == " g") {
      if(order=="LO") M2_approx = soft_g_squared(pp_full, M0_ij, A)*average_factor_full/average_factor;
      else if(order=="NLO") M2_approx = soft_g_squared_1l(pp_full, M0_ij, fM_ijk, M1_ij, A)*average_factor_full/average_factor;
    }
    else if(unresolved_str == " d d~")
      M2_approx = soft_qq_squared(pp_full, M0_ij, A)*average_factor_full/average_factor;
    else if(unresolved_str == " g g")
      M2_approx = soft_gg_squared(pp_full, M0_ij, M_ijkl, A)*average_factor_full/average_factor;
    else if(unresolved_str == " g d d~")
      M2_approx = soft_gqq_squared(pp_full, M0_ij, M_ijkl, dM_ijk, A)*average_factor_full/average_factor;
    //else if(unresolved_str == " g g g")
    //  M2_approx = soft_ggg_squared(pp_full, M_ij, M_ijkl, M_ijklab, dM_ijk, fM_ijkl)*average_factor_full/average_factor;

    std::cout << scale << "\t" << M2_Full << "\t" << M2_approx << "\t" << M2_Full/M2_approx << "\t" << std::abs(1. - M2_Full/M2_approx) << std::endl;
    outfile1 << scale << ", " << std::abs(1. - M2_Full/M2_approx) << std::endl;
  }
  }
  outfile1.close();
}