#ifndef __STATE_H
#define __STATE_H

#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <toml.hpp>

class StateDatabase {
public:
   std::string name;
   std::string description;
   toml::local_datetime created_on;
   std::string database_type;

   StateDatabase(): created_on(toml::local_datetime(std::chrono::system_clock::now())) {}

   StateDatabase(
     std::string arg_name,
     std::string arg_description,
     std::string arg_database_type,
     toml::local_datetime arg_created_on = toml::local_datetime(std::chrono::system_clock::now())
   ): name(arg_name),
      description(arg_description),
      database_type(arg_database_type),
      created_on(arg_created_on) {}

   toml::value dump_value() {
       return toml::value({
            {"name", name},
            {"description", description},
            {"created_on", created_on},
            {"database_type", database_type},
       });
   }

   static StateDatabase load_value(toml::value val) {
       StateDatabase state;
       state.name = toml::get<std::string>(toml::find(val, "name"));
       state.description = toml::get<std::string>(toml::find(val, "description"));
       state.created_on = toml::get<toml::local_datetime>(toml::find(val, "created_on"));
       state.database_type = toml::get<std::string>(toml::find(val, "database_type"));

       return state;
   }
};

class State {
public:
    std::map<std::string, StateDatabase> databases;
    std::vector<std::string> databases_ids;
    std::string version;

    explicit State(std::string arg_version): version(arg_version) {}

    StateDatabase set_database_description(std::string name, std::string description) {
        databases[name].description = description;
        return databases[name];
    }

    StateDatabase rename_database(std::string name, std::string new_name) {
        StateDatabase db_state = databases[name];
        remove_database(name);
        db_state.name = new_name;
        databases[new_name] = db_state;
        databases_ids.push_back(new_name);
        return db_state;
    }

    void remove_database(std::string name) {
        databases.erase(name);
        databases_ids.erase(std::remove(databases_ids.begin(), databases_ids.end(), name), databases_ids.end());
    }

    StateDatabase create_database(
        std::string arg_name,
        std::string arg_description,
        std::string arg_database_type,
        toml::local_datetime arg_created_on = toml::local_datetime(std::chrono::system_clock::now())
    ) {
        StateDatabase new_state(arg_name, arg_description, arg_database_type, arg_created_on);
        databases[arg_name] = new_state;
        databases_ids.push_back(arg_name);
        return new_state;
    }

    toml::value dump_value() {
        toml::table databases_val;
        for(auto& [k, v]: databases) {
            databases_val[k] = v.dump_value();
        }
        databases_val["count"] = databases.size();
        return toml::value({
            {"version", version},
            {"databases", databases_val},
        });
    }

    static State load_value(toml::value data) {
        State state("0.0.0");

        if (data.contains("databases")) {
            const auto& table_databases = toml::find(data, "databases");
            for(auto& [k, v] : table_databases.as_table()) {
                if (k != "count") {
                    auto dbState = StateDatabase::load_value(v);
                    state.databases[dbState.name] = dbState;
                }
            }
        }
        state.version = toml::get<std::string>(toml::find(data, "version"));

        std::vector<std::string> all_databases_ids;
        for(auto& [k, v] : state.databases) {
            all_databases_ids.push_back(k);
        }
        std::stable_sort(all_databases_ids.begin(), all_databases_ids.end());
        state.databases_ids = all_databases_ids;

        return state;
    }

    void dump(std::ofstream output) {
        output << dump_value();
    }

    static State load(std::ifstream input, std::string file_path) {
        return State::load_value(toml::parse(input, file_path));
    }
};

#endif