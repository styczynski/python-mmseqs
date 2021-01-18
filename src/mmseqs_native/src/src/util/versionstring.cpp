#include "Command.h"
#include "Debug.h"
#include "Util.h"
#include "output.h"

extern const char* version;

int versionstring(mmseqs_output* out, int, const char**, const Command&) {
    Debug(Debug::INFO) << version << "\n";
    EXIT(EXIT_SUCCESS);
}
