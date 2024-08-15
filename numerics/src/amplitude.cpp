#include "amplitude.hpp"

int str2int(std::string character) {
    // input: input string is a particle, like u, u~, d, d~ or g
    // ouput: returns 8 for glouns, 3 for quarks and 0 for anything else
    if (character.compare("g") == 0) {
        return 8;
    }
    else if (character.compare("u") == 0 or character.compare("d") == 0 or character.compare("c") == 0 or character.compare("s") == 0
          or character.compare("u~") == 0 or character.compare("d~") == 0 or character.compare("c~") == 0 or character.compare("s~") == 0){
        return 3;
    }
    else {
        return 1;
    }
}

int str2int_type(std::string character) {
    // input: input string is a particle, like u, u~, d, d~ or g
    // ouput: returns 8 for glouns, 3 for quarks and 0 for anything else
    if (character.compare("g") == 0) {
        return 2;
    }
    else if (character.compare("u") == 0 or character.compare("d") == 0 or character.compare("c") == 0 or character.compare("s") == 0){
        return 1;
    }
    else if (character.compare("u~") == 0 or character.compare("d~") == 0 or character.compare("c~") == 0 or character.compare("s~") == 0){
        return -1;
    }
    else {
        return 0;
    }
}

int scan_next_number(std::string name_str, int index, int *next_index) {
    // input: input string, and an index where to look for
    // return the next number after index as an integer. Also updates next_index to be the length of that number
    std::string number;
    int output;
    std::string space = " ";
    for (int i = 0; i < 3; i++){
        std::string character {name_str[index + i]};
        if (character.compare(space) != 0 and index + i != name_str.length()) {
            number.insert(i, character);
        }
        else {
            *next_index = i;
            output = str2int(number);
            return output;
        }
    }
    return 0;
}

int scan_next_number_type(std::string name_str, int index, int *next_index) {
    // input: input string, and an index where to look for
    // return the next number after index as an integer. Also updates next_index to be the length of that number
    std::string number;
    int output;
    std::string space = " ";
    for (int i = 0; i < 3; i++){
        std::string character {name_str[index + i]};
        if (character.compare(space) != 0 and index + i != name_str.length()) {
            number.insert(i, character);
        }
        else {
            *next_index = i;
            output = str2int_type(number);
            return output;
        }
    }
    return 0;
}

std::vector<int> scan_process(std::string name_str){
    // Delete "->" from neme_str
    int pos = name_str.find("->");
    name_str.erase(pos, 3);

    // Declare output vector
    std::vector<int> out;
    // Run through the input string
    int i = 0;
    while (i < name_str.length()) {
        int length_of_name;
        int particle_name = scan_next_number(name_str, i, &length_of_name);
        out.push_back(particle_name);
        if (length_of_name == 0) {
            i = i + 1;
        }
        i = i + length_of_name + 1;
    }
    return out;
}

std::vector<int> scan_process_type(std::string name_str){
    // Delete "->" from neme_str
    int pos = name_str.find("->");
    name_str.erase(pos, 3);

    // Declare output vector
    std::vector<int> out;
    // Run through the input string
    int i = 0;
    while (i < name_str.length()) {
        int length_of_name;
        int particle_name = scan_next_number_type(name_str, i, &length_of_name);
        if (out.size() < 2) {
            if(particle_name==2) particle_name = -2;
            out.push_back(-particle_name);
        }
        else {
            out.push_back(particle_name);
        }
        if (length_of_name == 0) {
            i = i + 1;
        }
        i = i + length_of_name + 1;
    }
    return out;
}


amplitude::amplitude(std::string process_str):process_str(process_str) {
    process = scan_process(process_str);
    particle_type = scan_process_type(process_str);
    int counter;
    for (int i = 0; i < process.size(); i++) {
        if (process[i] == 3) {
            quark_num++;
            if (particle_type[i] == -1) {
                up_index_pos.push_back(counter);
                counter++;
            }
            else if (particle_type[i] == 1) {
                up_index_pos.push_back(1000); // should never be be seen
            }
        }
        else if (process[i] == 8) {
            gluon_num++;
            gluon_positions.push_back(i);
            up_index_pos.push_back(counter);
            counter++;
        }
    }
    col_num = quark_num + gluon_num;
    uncol_num = process.size() - quark_num - gluon_num;

    // Get process_particles
    std::string space = " ";
    int i = 1;
    int index = 0;
    while (index + i < process_str.size()) {
        std::string character = process_str.substr(index + i, 1);
        if (character.compare(space) != 0) {
            i++;
        }
        else {
            std::string particle_string = process_str.substr(index, i);
            process_particles.push_back(particle_string);
            index = index + i + 1;
            i = 1;
        }
    }
    process_particles.push_back(process_str.substr(index, i));
    std::string arrow = "->";
    for (int i = 0; i < process_particles.size(); i++) {
        if (process_particles[i].compare(arrow) == 0) {
            process_particles.erase(process_particles.begin() + i);
        }
    }

    // Get quark pairs
    for (int i = 0; i < process_particles.size() - 1; i++) {
        std::string jTarget = "";
        if (process[i] == 3) {
            jTarget = process_particles[i];
        }
        std::string jTarget_saved = jTarget;
        for (int j = i + 1; j < process_particles.size(); j++) {
            jTarget = jTarget_saved;
            if ((i < 2 and j < 2) or (i >= 2 and j >= 2)){
                jTarget += "~";
            }
            // Replace "~~" -> ""
            if (jTarget.size() >= 3){
                jTarget = jTarget.substr(0, 1);
            }

            if (process_particles[j].compare(jTarget) == 0) {
                std::vector<int> quark_pair = {i, j};
                if (particle_type[i] == -1) {
                    quark_pair = {j, i};
                }
                quark_pairs.push_back(quark_pair);
            }
        }
    }

    // Get off-diagonals
    for (int i = 0; i < quark_pairs.size(); i++) {
        std::string off_diagonal = "";
        for (int j = 0; j < process_particles.size(); j++) {
            if (j == 2) {
                off_diagonal += "-> ";
            }
            if ((j == quark_pairs[i][0]) or (j == quark_pairs[i][1])) {
                off_diagonal += "g";
            }
            else {
                off_diagonal += process_particles[j];
            }
            if (j != process_particles.size() - 1) {
                off_diagonal += " ";
            }
        }
        off_diagonals.push_back(off_diagonal);
    }

    for (int i = 0; i < quark_pairs.size(); i++) {
        off_diagonals_ids.push_back(i+3);
    }
}


