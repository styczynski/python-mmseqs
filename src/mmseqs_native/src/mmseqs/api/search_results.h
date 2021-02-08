#ifndef API_SEARCH_RESULTS_H
#define API_SEARCH_RESULTS_H

#include <mmseqs/output.h>

class SearchResults {
public:

    std::vector<std::vector<mmseqs_blast_tab_record>> _records;
    std::vector<std::string> _headers;

    SearchResults();

    explicit SearchResults(
        const std::vector<std::vector<mmseqs_blast_tab_record>>& records,
        const std::vector<std::string>& headers
    );

    const std::vector<std::vector<mmseqs_blast_tab_record>>& getRecords() const;
};

#endif