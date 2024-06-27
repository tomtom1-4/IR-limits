#ifndef AMPLITUDE_HEADER
#define AMPLITUDE_HEADER

#include <string>
#include <vector>

class amplitude {
  public:
    std::string process_str;
    std::vector<int> process; // Contains a list with 3 for every anti-quark, 3 for every quark and 8 for every gluon; 1 for any uncharged particle
    std::vector<int> particle_type; // Contains a list with 1 for every final state quark/initial state anti-quark
                                    //                     -1 for every initial state quark/final state anti-quark
                                    //                      2 for every gluon
                                    //                      0 for any uncharged particle
    int quark_num = 0;
    int gluon_num = 0;
    int uncol_num = 0;
    int col_num = 0;
    std::vector<int> gluon_positions;
    std::vector<int> up_index_pos; // Contains a list with the index of all up index particles, e.g. (1, 1000, 1000, 2, 3) for u u~ -> d d~ g

    std::vector<std::string> process_particles;
    std::vector<std::vector<int>> quark_pairs; // Contains list of quark positions
    std::vector<std::string> off_diagonals;
    std::vector<int> off_diagonals_ids;
    amplitude(std::string);
    ~amplitude(){}

};

#endif