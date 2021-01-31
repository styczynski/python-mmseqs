#include <mmseqs/commons/commandCaller.h>
#include <mmseqs/output.h>
#include <mmseqs/commons/util.h>

#include <strings.h>
#include <unistd.h>
#include <cstdlib>
#include <sstream>

#ifdef OPENMP
#include <omp.h>
#endif

CommandCaller::CommandCaller(mmseqs_output* output): out(output) {
#ifdef OPENMP
#if _OPENMP >= 201307
  if (omp_get_proc_bind() != omp_proc_bind_false) {
#else
  char* procBind = getenv("OMP_PROC_BIND");
  if (procBind != NULL && strcasecmp(procBind, "false") != 0 &&
      strcasecmp(procBind, "0") != 0) {
#endif
    out->failure("Calling program has OMP_PROC_BIND set in its environment. Please unset OMP_PROC_BIND.");
  }
#endif

  std::string depth = SSTR(getCallDepth());
  addVariable("MMSEQS_CALL_DEPTH", depth.c_str());
}

unsigned int CommandCaller::getCallDepth() {
  char* currentCallDepth = getenv("MMSEQS_CALL_DEPTH");
  if (currentCallDepth == NULL) {
    return 0;
  }

  char* rest;
  int depth = strtol(currentCallDepth, &rest, 10);
  if (rest == currentCallDepth || errno == ERANGE) {
    out->failure("Invalid non-numeric value for environment variable MMSEQS_CALL_DEPTH");
  }
  return depth + 1;
}

void CommandCaller::addVar(std::string key, std::string value) {
  addVariable(key.c_str(), value.c_str());
}

void CommandCaller::addVariable(const char* key, const char* value) {
  if (value == NULL) {
    unsetenv(key);
  } else {
    setenv(key, value, true);
  }
}

int CommandCaller::callProgram(const char* program, size_t argc,
                               const char** argv) {
  std::stringstream argStream(program);
  for (size_t i = 0; i < argc; i++) {
    argStream << " " << argv[i];
  }

  std::string argString = argStream.str();
  if (std::system(argString.c_str()) != EXIT_SUCCESS) {
    EXIT(EXIT_FAILURE);
  }

  return 0;
}

void CommandCaller::execProgram(const char* program,
                                const std::vector<std::string>& argv) {
  std::cerr.flush();
  std::cout.flush();
  // hack: our argv string does not contain a program name anymore, readd it
  const char** pArgv = new const char*[argv.size() + 2];
  pArgv[0] = program;
  for (size_t i = 0; i < argv.size(); ++i) {
    pArgv[i + 1] = argv[i].c_str();
  }
  pArgv[argv.size() + 1] = NULL;

  int res = execvp(program, (char* const*)pArgv);

  if (res == -1) {
    delete[] pArgv;
    out->failure("Failed to execute {} with error {}", program, errno);
  }

  // should not be reached in the normal case
  delete[] pArgv;
  EXIT(EXIT_FAILURE);
}

int CommandCaller::callProgram(const char* program,
                               const std::vector<std::string>& argv) {
  std::stringstream argStream;
  argStream << "bash " << std::string(program);
  for (size_t i = 0; i < argv.size(); i++) {
    argStream << " " << argv[i];
  }

  std::string argString = argStream.str();
  if (std::system(argString.c_str()) != EXIT_SUCCESS) {
    EXIT(EXIT_FAILURE);
  }

  return 0;
}
