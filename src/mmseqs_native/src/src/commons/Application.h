struct mmseqs_call_args {
    mmseqs_call_args() {
        cli_args = std::vector<std::string>();
    }
    std::vector<std::string> cli_args;
};

int call_mmseqs(mmseqs_call_args args);