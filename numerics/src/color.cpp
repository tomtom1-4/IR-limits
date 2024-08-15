#include "color.hpp"


std::complex<double> f_abc(int a, int b, int c) {
    std::complex<double> TbTc[9];
    complex_mult(lam[b], lam[c], TbTc, 3);
    std::complex<double> TcTb[9];
    complex_mult(lam[c], lam[b], TcTb, 3);
    std::complex<double> TaTbTc[9];
    complex_mult(lam[a], TbTc, TaTbTc, 3);
    std::complex<double> TaTcTb[9];
    complex_mult(lam[a], TcTb, TaTcTb, 3);
    std::complex<double> fabc = -I/T_F * (TaTbTc[0] + TaTbTc[4] + TaTbTc[8] - TaTcTb[0] - TaTcTb[4] - TaTcTb[8])/8. ;
    return fabc;
}

std::complex<double> cOlT(int a, int b, int c) {
    return lam[c][3 * a + b]/2.;
}

std::complex<double> kronecker(int c1, int c2) {
    if (c1 == c2) {
        return 1.;
    }
    else {
        return 0.;
    }
}

std::string get_key(int* hel, int* col, int length) {
    std::string hel_string;
    std::string col_string;
    for (int i = 0; i < length; i++) {
        hel_string += std::to_string(hel[i]);
        col_string += std::to_string(col[i]);
    }
    std::string key = hel_string + col_string;
    return key;
}

int coeff2baseN(int *coeff, int base, int length) {
    int loc = 0;
    for (int i = 0; i < length; i++) {
        loc = loc + coeff[i] * pow(base, i);
    }
    return loc;
}

std::vector<int> baseN2coeff(int value, int base, int length) {
    int loc = 0;
    std::vector<int> coeff;
    coeff.assign(length, 0);
    for (int i = 0; i < length; i++) {
        coeff[length - i - 1] = (value - loc)/pow(base, length - i - 1);
        loc = loc + coeff[length - i - 1] * pow(base, length - i - 1);
    }
    return coeff;
}

std::vector<int> remove0(int* in, int length) {
    std::vector<int> out;
    for (int i = 0; i < length; i++) {
        if (in[i] != 0) {
            out.push_back(in[i]);
        }
    }
    return out;
}

std::vector<int> remove0(std::vector<int> in) {
    std::vector<int> out;
    for (int i = 0; i < in.size(); i++) {
        if (in[i] != 0) {
            out.push_back(in[i]);
        }
    }
    return out;
}

bool contain(std::vector<int> list, int element) {
    // Checks if element is contained in list
    for (int i = 0; i < list.size(); i++) {
        if (list[i] == element) {
            return true;
            break;
        }
    }
    return false;
}

std::vector<int> missing_and_gluon(std::vector<int> colb, std::vector<int> process) {
    //creates list with all missing indices (initial anti-quarks & final quarks) and gluons
    std::vector<int> out;
    for (int i = 1; i <= process.size(); i++) {
        if (contain(colb, i) and process[i - 1] == 8) { // position i is a gluon so it is added to the list
            out.push_back(i);
        }
        else if (!contain(colb, i)) { // position i is not contained in colb so it is added to the list
            out.push_back(i);
        }
    }
    return out;
}

std::vector<int> gluon_order(std::vector<int> process, std::vector<int> colb) {
    // Returns a list of ordered gluons, e.g. if colb = {4,1,3} with 1 = quark and 4, 3 = gluon then the output is {1, 0}
    std::vector<int> out;
    // remove quarks
    for (int i = 0; i < colb.size(); i++) {
        if (process[colb[i] - 1] == 8) {
            out.push_back(colb[i]);
        }
    }
    // out = {4, 3}
    // lower indices to 0,1,...
    int counter = 0;
    for (int i = 1; i <= process.size(); i++) {
        for (int j = 0; j < out.size(); j++) {
            if (out[j] == i) {
                out[j] = counter;
                counter++;
            }
        }
    }
    return out;
}

std::vector<int> partial_color(std::vector<int> colb, int* colors, std::vector<int> gluon_colors, std::vector<int> process) {
    // Return list of colors of the respective colorbasis.
    std::vector<int> out;
    int gluon_counter = 0;
    std::vector<int> gluon_list = gluon_order(process, colb);
    for (int i = 0; i < colb.size(); i++) {
        if (process[colb[i] - 1] == 3) {
            out.push_back(colors[colb[i] - 1]);
        }
        else if (process[colb[i] - 1] == 8) {
            out.push_back(gluon_colors[gluon_list[gluon_counter]]);
            gluon_counter++;
        }
    }
    return out;
}

std::vector<std::vector<int> > split_gluon_color(std::vector<int> gluon_col) {
    // Return a list of list with spitted gluon colors, e.g. [5, 5, 6, 3, 1, 4] -> [[5, 6, 1], [5, 3, 4]]
    std::vector<std::vector<int> > out;
    std::vector<int> out1, out2;
    for (int i = 0; 2 * i + 1 < gluon_col.size(); i++) {
        out1.push_back(gluon_col[2 * i]);
        out2.push_back(gluon_col[2 * i + 1]);
    }
    out.push_back(out1);
    out.push_back(out2);
    return out;
}

void hel_plus(int *hel, int *hel_out, int length) {
    int i = 0;
    while (i < length) {
        if (hel[i] == -1) {
            hel_out[i] = 1;
            break;
        }
        else {
            i = i + 1;
            for (int j = 0; j < i; j++) {
                hel_out[j] = -1;
            }
        }
    }
}

void col_plus(int *col, int *col_out, amplitude& A) {
    for (int i = 0; i < A.process.size(); i++) {
        if (col[i] < A.process[i] - 1) {
             col_out[i] = col[i] + 1;
             break;
        }
        else {
            col_out[i] = 0;
        }
    }
}

std::vector<std::vector<int> > non0hel(int id , int power, std::string loop_order, amplitude& A) {
    std::vector<std::vector<int>> helnon0;
    int ncolb, colb[A.process.size()];
    Recola::get_n_colour_configurations_rcl(id, ncolb);
    // Sum over all helicities
    // Create start array
    int hel[A.process.size()], hel_stop[A.process.size()];
    for (int i = 0; i < A.process.size(); i++) {
        hel[i] = -1;  // Let's hope we ain't never gonna do gravitons o7
        hel_stop[i] = 1;
    }
    bool flag = true;
    while (flag) {
        std::complex<double> M_tot;
		double M2;
		Recola::get_squared_amplitude_rcl(id, power, loop_order, M2);
        for (int i = 1; i <= ncolb; i++) {
            Recola::get_colour_configuration_rcl(id, i, colb);  // Get the basis configuration

            // Get Coefficient infront of the basis vectors
            std::complex<double> M_component;
            Recola::get_amplitude_rcl(id, power, loop_order, colb, hel, M_component);
            M_tot += std::abs(M_component);
            std::cout << A.process_str << ", " << power << ", " << loop_order << ": ";
            for (int j = 0; j < A.process.size(); j++) std::cout << colb[j] << ", ";
            std::cout << "  ";
            for (int j = 0; j < A.process.size(); j++) std::cout << hel[j] << ", ";
            std::cout << "\t" << M_component << std::endl;
        }
        if (M_tot != 0.) {
            std::vector<int> hel_vec;
            for (int j = 0; j < A.process.size(); j++) {
                hel_vec.push_back(hel[j]);
            }
            helnon0.push_back(hel_vec);
        }
        if (arrayeq(hel, hel_stop, A.process.size())) {
            flag = false;
        }
        hel_plus(hel, hel, A.process.size());
    }
    return helnon0;
}

std::complex<double> colorflow2color(int *hel, int *col, amplitude& A, int power, std::string loop_order, int id) {
    // Computes the amplitude of some process with helicity hel and color col
    std::complex<double> M = 0;
    // Basis transformation to usual basis
    int ncolb, colb[A.process.size()];
    Recola::get_n_colour_configurations_rcl(id, ncolb);
    // Extract gluon color
    std::vector<int> gluon_color;
    for (int j = 0; j < A.process.size(); j++) {
        if (A.process[j] == 8) {
            gluon_color.push_back(col[j]);
        }
    }
    // Sum over gluon colors as two quark indices (c1,d1, c2, d2, ...)
    for (int i = 0; i < pow(3, 2 * A.gluon_num); i++) {
        std::vector<int> gluon_col = baseN2coeff(i, 3,  2 * A.gluon_num);  // Get from base 10 to the color coefficients
        //Split Gluon color in upper and lower indices
        std::vector<int> gluon_col_up = split_gluon_color(gluon_col)[0];
        std::vector<int> gluon_col_down = split_gluon_color(gluon_col)[1];
        std::complex<double> M_colorflow = 0; // Result in colorflow basis
        // compute colorfactor, e.g. multiply by the gellmann matrices
        std::complex<double> col_fac = 1;
        for (int k = 0; k < A.gluon_num; k++) {
            col_fac *= lam[gluon_color[k]][gluon_col_up[k] * 3 + gluon_col_down[k]]/std::sqrt(2);
        }
        if (col_fac != 0.) {
            // Sum over all basis vectors
            for (int j = 1; j <= ncolb; j++) {
                Recola::get_colour_configuration_rcl(id, j, colb);  // Get the basis configuration
                std::vector<int> colb_vec = remove0(colb, A.process.size()); // Remove zeros
                std::vector<int> dual_colb = missing_and_gluon(colb_vec, A.process); // Get dual basis

                // Get Coefficient infront of the basis vectors
                std::complex<double> M_component;
                Recola::get_amplitude_rcl(id, power, loop_order, colb, hel, M_component);
                //std::cout << id << ", " << loop_order << "\t" << hel[0] << ", " << hel[1] << ", " << hel[2] << ", " << hel[3] << "\t"
                // << colb[0] << "\t" << colb[1] << "\t" << colb[2] << "\t" << colb[3] << "\t" << M_component << std::endl;

                // Get colors of quarks and gluons, where gluons have color fixed by summation
                std::vector<int> up_col = partial_color(colb_vec, col, gluon_col_up, A.process);
                std::vector<int> down_col = partial_color(dual_colb, col, gluon_col_down, A.process);

                // Kronecker delta and summation
                if (up_col == down_col) {
                    M_colorflow += M_component;
                }
            }
        }
        M += M_colorflow * col_fac;
    }
    return M;
}

std::vector<std::vector<int> > non0col(int id, int power, std::string loop_order, amplitude& A, std::vector<std::vector<int> > helnon0) {
    // Returns a list with all non-vanishing color configurations
    std::vector<std::vector<int> > colnon0;
    // Failsafe for empty helnon0
    if (helnon0.size() == 0) {
        std::cout << "process impossible" << std::endl;
        std::vector<int> empty = {};
        colnon0.push_back(empty);
        return colnon0;
    }
    int col[A.process.size()];
    int col_stop[A.process.size()];
    for (int i = 0; i < A.process.size(); i++) {
        col[i] = 0;
        col_stop[i] = A.process[i] - 1;
    }
    // Loop over non-vanishing helicity configurations
    int i = 0;
    while (i < helnon0.size()) {
        bool col_flag = true;
        std::vector<int> hel_vec = helnon0[i];
        int hel[A.process.size()];
        vector2arr(hel_vec, hel);
        std::cout << i << ": " << hel[0] << ", " << hel[1] << ", " << hel[2] << ", " << hel[3] << std::endl;
        while (col_flag) {
            std::complex<double> m_full = colorflow2color(hel, col, A, power, loop_order, id);
            if (A.process_str == "u u~ -> g g") {
                std::cout << helnon0.size() << "\t";
                std::cout << hel[0] << ", " << hel[1] << ", " << hel[2] << ", " << hel[3] << "\t";
                std::cout << col[0] << ", " << col[1] << ", " << col[2] << ", " << col[3] << "\t";
                std::cout << m_full << std::endl;
            }
            if ( std::abs(m_full) > 1.e-17 ) {
                std::vector<int> col_vec;
                for (int k = 0; k < A.process.size(); k++) {
                    col_vec.push_back(col[k]);
                }
                colnon0.push_back(col_vec);
            }
            if (arrayeq(col, col_stop, A.process.size())) {
                col_flag = false;
            }
            col_plus(col, col, A);
        }
        if (colnon0.size() == 0) {
            i += 1;
        }
        else {
            i = helnon0.size();
        }
    }
    return colnon0;
}