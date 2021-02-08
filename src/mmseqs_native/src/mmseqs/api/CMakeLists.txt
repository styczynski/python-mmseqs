set(api_header_files
        mmseqs/api/client.h
        mmseqs/api/database.h
        mmseqs/api/rand.h
        mmseqs/api/utils.h
        mmseqs/api/parameters.h
        mmseqs/api/search_results.h
        PARENT_SCOPE
        )

set(api_source_files
        mmseqs/api/client.cpp
        mmseqs/api/database.cpp
        mmseqs/api/search_results.cpp
        PARENT_SCOPE
        )
