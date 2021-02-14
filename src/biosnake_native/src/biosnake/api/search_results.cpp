#include <biosnake/api/search_results.h>

SearchResults::SearchResults() {}

SearchResults::SearchResults(
    const std::vector<std::vector<biosnake_blast_tab_record>>& records,
    const std::vector<std::string>& headers
): _records(records), _headers(headers) {}

const std::vector<std::vector<biosnake_blast_tab_record>>& SearchResults::getRecords() const {
    return _records;
}
