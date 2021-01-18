#ifndef OUTPUT_H
#define OUTPUT_H

#include <string>
#include <map>

struct mmseqs_output {
    mmseqs_output() {
        vars_str = std::map<std::string, std::string>();
    }
    std::map<std::string, std::string> vars_str;

    void output_string(std::string name, std::string value) {
        vars_str[name] = value;
    }
};

#endif