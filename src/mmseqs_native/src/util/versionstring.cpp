#include "Command.h"
#include "Parameters.h"
#include "Debug.h"
#include "Util.h"
#include "output.h"

int versionstring(mmseqs_output* out, Parameters &par) {
    Debug(Debug::INFO) << version << "\n";
    EXIT(EXIT_SUCCESS);
}
