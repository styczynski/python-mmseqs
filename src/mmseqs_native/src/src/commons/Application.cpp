#include "Debug.h"
#include "Parameters.h"
#include "Util.h"
#include "Command.h"
//include "DistanceCalculator.h"
#include "Timer.h"
#include "Application.h"
#include "output.h"

#include <iomanip>

extern const char *binary_name;
extern const char *tool_name;
extern const char *tool_introduction;
extern const char *main_author;
extern const char *version;
extern const char *show_extended_help;
extern const char *show_bash_info;
extern bool hide_base_commands;

extern std::vector<Command> commands;
extern std::vector<Command> baseCommands;
extern std::vector<Categories> categories;

Command *getCommandByName(const char *s) {
    for (size_t i = 0; i < commands.size(); i++) {
        Command &p = commands[i];
        if (!strcmp(s, p.cmd))
            return &p;
    }
    for (size_t i = 0; i < baseCommands.size(); i++) {
        Command &p = baseCommands[i];
        if (!strcmp(s, p.cmd))
            return &p;
    }
    return NULL;
}

int runCommand(mmseqs_output* out, Command* c, Parameters &par) {
    Timer timer;
    int status = c->commandFunction(out, par);
    Debug(Debug::INFO) << "Time for processing: " << timer.lap() << "\n";
    return status;
}

void printUsage(bool showExtended) {
    std::stringstream usage;

    usage << tool_introduction << "\n\n";
    usage << tool_name << " Version: " << version << "\n";
    usage << "© " << main_author << "\n\n";
    usage << "usage: " << binary_name << " <command> [<args>]" << "\n";

    std::vector<int> showCategoryHeader(categories.size(), 0);
    for (size_t i = 0; i < categories.size(); ++i) {
        for (size_t j = 0; j < commands.size(); j++) {
            Command &p = commands[j];
            if (p.mode & categories[i].mode) {
                showCategoryHeader[i] = 1;
                break;
            }
        }
        if (hide_base_commands) {
            continue;
        }
        for (size_t j = 0; j < baseCommands.size(); j++) {
            Command &p = baseCommands[j];
            if (p.mode & categories[i].mode) {
                showCategoryHeader[i] = 1;
                break;
            }
        }
    }


    for (size_t i = 0; i < categories.size(); ++i) {
        if (showExtended == false
            && (categories[i].mode & COMMAND_MAIN) == 0
            && (categories[i].mode & COMMAND_EASY) == 0
            && (categories[i].mode & COMMAND_DATABASE_CREATION) == 0
            && (categories[i].mode & COMMAND_FORMAT_CONVERSION) == 0
            ) {
            continue;
        }

        if (showCategoryHeader[i] == 0) {
            continue;
        }

        usage << "\n" << std::setw(20) << categories[i].title << "\n";
        for (size_t j = 0; j < commands.size(); j++) {
            struct Command &p = commands[j];
            if (showExtended == false && (p.mode & COMMAND_EXPERT) != 0) {
                continue;
            }
            if (p.mode & categories[i].mode) {
                usage << std::left << std::setw(20) << "  " + std::string(p.cmd) << "\t" << p.description << "\n";
            }
        }
        if (hide_base_commands) {
            continue;
        }
        for (size_t j = 0; j < baseCommands.size(); j++) {
            struct Command &p = baseCommands[j];
            if (showExtended == false && (p.mode & COMMAND_EXPERT) != 0) {
                continue;
            }
            if (p.mode & categories[i].mode) {
                usage << std::left << std::setw(20) << "  " + std::string(p.cmd) << "\t" << p.description << "\n";
            }
        }
    }

    if (show_extended_help != NULL) {
        if (showExtended == false) {
            usage << "\nAn extended list of all modules can be obtained by calling '" << binary_name << " -h'.\n";
        }
    }
    if (show_bash_info != NULL) {
        usage  << "\nBash completion for modules and parameters can be installed by adding \"source MMSEQS_HOME/util/bash-completion.sh\" to your \"$HOME/.bash_profile\".\nInclude the location of the " << tool_name << " binary in your \"$PATH\" environment variable.";
    }
    Debug(Debug::INFO) << usage.str() << "\n";
}

int shellcompletion(int argc, const char **argv) {
    // mmseqs programs
    if (argc == 0) {
        for (size_t i = 0; i < commands.size(); i++) {
            struct Command &p = commands[i];
            if (p.mode & COMMAND_HIDDEN)
                continue;
            Debug(Debug::INFO) << p.cmd << " ";
        }
        if (hide_base_commands == false) {
            for (size_t i = 0; i < baseCommands.size(); i++) {
                struct Command &p = baseCommands[i];
                if (p.mode & COMMAND_HIDDEN)
                    continue;
                Debug(Debug::INFO) << p.cmd << " ";
            }
        }
        Debug(Debug::INFO) << "\n";
    }

    // mmseqs parameters for given program
    if (argc == 1) {
        for (size_t i = 0; i < commands.size(); i++) {
            struct Command &p = commands[i];
            if (strcmp(p.cmd, argv[0]) != 0) {
                continue;
            }
            if (p.params == NULL) {
                continue;
            }
            for (std::vector<MMseqsParameter *>::const_iterator it = p.params->begin(); it != p.params->end(); ++it) {
                Debug(Debug::INFO) << (*it)->name << " ";
            }
            Debug(Debug::INFO) << "\n";
            break;
        }
        if (hide_base_commands == false) {
            for (size_t i = 0; i < baseCommands.size(); i++) {
                struct Command &p = baseCommands[i];
                if (strcmp(p.cmd, argv[0]) != 0) {
                    continue;
                }
                if (p.params == NULL) {
                    continue;
                }
                for (std::vector<MMseqsParameter *>::const_iterator it = p.params->begin();
                     it != p.params->end(); ++it) {
                    Debug(Debug::INFO) << (*it)->name << " ";
                }
                Debug(Debug::INFO) << "\n";
                break;
            }
        }
        Debug(Debug::INFO) << "\n";
    }
    return EXIT_SUCCESS;
}


void subcall_mmseqs(mmseqs_output* out, std::string command_name, Parameters args) {

    Parameters par(args);

//    const int argc = args.cli_args.size();
//    char** argvp = (char**) malloc(argc * sizeof(char*));
//
//    for (int i=0; i<argc; ++i) {
//        const int len = args.cli_args[i].size();
//        argvp[i] = (char*) malloc((len+1) * sizeof(char));
//        for (int j=0;j<len;++j) {
//            argvp[i][j] = args.cli_args[i][j];
//        }
//        argvp[i][len] = '\0';
//    }
//    const char** argv = (const char**) argvp;
//
//    std::cout << "We will execute:\n";
//    for (int i=0;i<argc;++i) {
//        std::cout << "[" << std::string(argv[i]) << "] ";
//    }
//    std::cout << "\n-- END --\n";
//
//    if (argc < 2) {
//        printUsage(false);
//        return out;
//    }
//
//    if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0){
//        printUsage(true);
//        return out;
//    }
//
//    setenv("MMSEQS", argv[0], true);
    Command *c = getCommandByName(command_name.c_str());
    runCommand(out, c, par);

//
//    if (strncmp(argv[1], "shellcompletion", strlen("shellcompletion")) == 0) {
//        shellcompletion(argc - 2, argv + 2);
//        return out;
//    } else if ((c = getCommandByName(argv[1])) != NULL) {
//        runCommand(&out, c, argc - 2, argv + 2);
//        return out;
//    } else {
//        printUsage(true);
//        Debug(Debug::INFO) << "\nInvalid Command: " << argv[1] << "\n";
//        return out;
//    }
//
//    return out;
}

mmseqs_output call_mmseqs(std::string command_name, Parameters args) {
    mmseqs_output out;
    subcall_mmseqs(&out, command_name, args);
    return out;
}