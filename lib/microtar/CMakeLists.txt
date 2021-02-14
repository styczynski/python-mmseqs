add_library(microtar microtar.h microtar.c)
set_target_properties(microtar PROPERTIES COMPILE_FLAGS "${BIOSNAKE_C_FLAGS}" LINK_FLAGS "${BIOSNAKE_C_FLAGS}")
