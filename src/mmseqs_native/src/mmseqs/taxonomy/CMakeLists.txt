set(taxonomy_header_files
        mmseqs/taxonomy/ncbiTaxonomy.h
        PARENT_SCOPE
        )


set(taxonomy_source_files
        mmseqs/taxonomy/lca.cpp
        mmseqs/taxonomy/addtaxonomy.cpp
        mmseqs/taxonomy/ncbiTaxonomy.cpp
        mmseqs/taxonomy/filtertaxdb.cpp
        mmseqs/taxonomy/filtertaxseqdb.cpp
        mmseqs/taxonomy/aggregatetax.cpp
        mmseqs/taxonomy/createtaxdb.cpp
        mmseqs/taxonomy/createbintaxonomy.cpp
        mmseqs/taxonomy/taxonomyreport.cpp
        mmseqs/taxonomy/taxonomyExpression.h
        PARENT_SCOPE
        )
