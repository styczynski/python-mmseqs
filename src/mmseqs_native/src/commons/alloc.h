#ifndef ALLOC_H
#define ALLOC_H

#include <stdlib.h>

//define NO_VLA

#ifdef NO_VLA
#define VLA(TYPE, NAME, SIZE) TYPE* NAME = (TYPE*) malloc((SIZE) * sizeof(TYPE))
#else
#if _MSC_VER && !__INTEL_COMPILER
#define VLA(TYPE, NAME, SIZE) TYPE* NAME = (TYPE*) _malloca((SIZE) * sizeof(TYPE))
#else
#define VLA(TYPE, NAME, SIZE) TYPE NAME[(SIZE)]
#endif
#endif // NO_VLA

#endif