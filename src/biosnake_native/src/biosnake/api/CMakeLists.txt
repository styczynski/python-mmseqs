set(api_header_files
        biosnake/api/client.h
        biosnake/api/database.h
        biosnake/api/rand.h
        biosnake/api/utils.h
        biosnake/api/parameters.h
        biosnake/api/search_results.h
        PARENT_SCOPE
        )

set(api_source_files
        biosnake/api/client.cpp
        biosnake/api/database.cpp
        biosnake/api/search_results.cpp
        PARENT_SCOPE
        )
