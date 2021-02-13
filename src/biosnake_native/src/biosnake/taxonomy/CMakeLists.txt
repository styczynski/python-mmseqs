set(taxonomy_header_files
        biosnake/taxonomy/ncbiTaxonomy.h
        PARENT_SCOPE
        )


set(taxonomy_source_files
        biosnake/taxonomy/lca.cpp
        biosnake/taxonomy/addtaxonomy.cpp
        biosnake/taxonomy/ncbiTaxonomy.cpp
        biosnake/taxonomy/filtertaxdb.cpp
        biosnake/taxonomy/filtertaxseqdb.cpp
        biosnake/taxonomy/aggregatetax.cpp
        biosnake/taxonomy/createtaxdb.cpp
        biosnake/taxonomy/createbintaxonomy.cpp
        biosnake/taxonomy/taxonomyreport.cpp
        biosnake/taxonomy/taxonomyExpression.h
        PARENT_SCOPE
        )
