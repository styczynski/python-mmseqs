#include "state.h"

int main() {
    State state;
    state.version = "1.0.0";

    StateDatabase db;
    db.name = "db1";
    db.description = "My description";
    db.database_type = "nucleotide";
    state.databases["db1"] = db;

    state.dump(std::ofstream("config.toml", std::ios_base::binary));
    state = State::load(std::ifstream("config.toml", std::ios_base::binary), "config.toml");
    return 0;
}