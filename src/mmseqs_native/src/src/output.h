#ifndef OUTPUT_H
#define OUTPUT_H

#include <iostream>
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

    void print() {
        std::cout << "Output{\n";
        for(auto it = vars_str.cbegin(); it != vars_str.cend(); ++it) {
            std::cout << "  " << it->first << ": [" << it->second << "]\n";
        }
        std::cout << "}\n";
    }
};

#endif