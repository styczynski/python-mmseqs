#include "Command.h"
#include "Parameters.h"
#include "Debug.h"
#include "Util.h"
#include "output.h"

extern const char* version;

int versionstring(mmseqs_output* out, Parameters &par) {
    Debug(Debug::INFO) << version << "\n";
    EXIT(EXIT_SUCCESS);
}
