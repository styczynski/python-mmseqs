#ifndef API_UTILS_H
#define API_UTILS_H

#include <string>
#include <vector>
#include <fstream>

#include <mmseqs/api/rand.h>

std::string write_temp_fasta(std::vector<std::string> sequences, std::string temp_path) {
    std::string file_path = temp_path + get_uuid() + "query_inputs.fasta";
    std::ofstream out(file_path, std::ios_base::binary);
    for (int i=0;i<sequences.size();++i) {
        out << ">query" + SSTR(i) + "\n" + sequences[i] + "\n";
    }
    return file_path;
}

#endif