#include <snt/dip/cursor.h>
#include <snt/dip/dip.h>
#include <snt/dip/environment.h>
#include <snt/puq/quantity.h>
#include <snt/val/values_number.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numbers>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace fs = std::filesystem;
using snt::dip::Environment;

namespace {

struct Nuclide {
    std::string id;
    std::string label;
    std::string daughter;
    std::string mode;
    double atomic_mass_g_mol{};
    double half_life_s{};
    double branching{};
    double energy_mev{};
    bool stable{};
    double lambda() const { return stable ? 0.0 : std::numbers::ln2 / half_life_s; }
};

using BatemanTerms = std::map<double, double>;

// The bundled schema is a serial chain (one explicit daughter per nuclide).
// Representing each inventory as a sum of exp(-lambda*t) terms is the exact
// Bateman solution and avoids the stiffness/conditioning trap of timestepping
// microsecond and gigayear half-lives on one reporting grid.
std::vector<BatemanTerms> bateman_terms(const std::vector<Nuclide>& chain, double initial_atoms) {
    std::vector<BatemanTerms> terms(chain.size());
    terms[0].emplace(chain[0].lambda(), initial_atoms);
    for (std::size_t i = 1; i < chain.size(); ++i) {
        const double daughter_lambda = chain[i].lambda();
        const double feed = chain[i - 1].branching * chain[i - 1].lambda();
        for (const auto& [parent_lambda, coefficient] : terms[i - 1]) {
            const double denominator = daughter_lambda - parent_lambda;
            if (std::abs(denominator) < 1e-30)
                throw std::runtime_error("equal decay constants require a repeated-root Bateman term");
            terms[i][parent_lambda] += feed * coefficient / denominator;
            terms[i][daughter_lambda] -= feed * coefficient / denominator;
        }
    }
    return terms;
}

std::vector<double> inventory_at(const std::vector<BatemanTerms>& terms, double time_s) {
    if (time_s == 0.0) {
        std::vector<double> initial(terms.size());
        for (const auto& [_, coefficient] : terms.front()) initial[0] += coefficient;
        return initial;
    }
    std::vector<double> inventory;
    inventory.reserve(terms.size());
    for (const auto& isotope_terms : terms) {
        double atoms = 0.0;
        for (const auto& [lambda, coefficient] : isotope_terms) atoms += coefficient * std::exp(-lambda * time_s);
        inventory.push_back(std::max(0.0, atoms));
    }
    return inventory;
}

double quantity_as(const Environment& env, const std::string& path, const std::string& unit) {
    const auto cursor = env[path];
    const auto source_unit = cursor.get_units();
    if (!source_unit) throw std::runtime_error("parameter '" + path + "' has no physical unit");
    const snt::puq::Quantity quantity(cursor.as<double>(), source_unit->to_string());
    const auto converted = quantity.convert(unit);
    const auto* values = dynamic_cast<const snt::val::ArrayValue<double>*>(converted.measurement.result.estimate.get());
    if (!values || values->get_size() != 1) throw std::runtime_error("parameter '" + path + "' is not a scalar float");
    return values->get_value(0);
}

Nuclide read_nuclide(const Environment& env, const std::string& id) {
    const std::string root = "nuclear.isotopes." + id;
    Nuclide n;
    n.id = env[root + ".id"].as<std::string>();
    n.label = env[root + ".label"].as<std::string>();
    n.daughter = env[root + ".decay.daughter"].as<std::string>();
    n.mode = env[root + ".decay.mode"].as<std::string>();
    n.atomic_mass_g_mol = quantity_as(env, root + ".atomic_mass", "g/mol");
    n.half_life_s = quantity_as(env, root + ".half_life", "s");
    n.branching = env[root + ".decay.branching"].as<double>();
    n.energy_mev = quantity_as(env, root + ".decay.energy", "MeV");
    n.stable = env[root + ".stable"].as<bool>();
    if (!n.stable && n.half_life_s <= 0.0) throw std::runtime_error(id + " has a non-positive half-life");
    if (n.branching < 0.0 || n.branching > 1.0) throw std::runtime_error(id + " has an invalid branching ratio");
    return n;
}

void append_chain(const Environment& env, const std::string& id, std::vector<Nuclide>& chain) {
    if (std::any_of(chain.begin(), chain.end(), [&](const Nuclide& n) { return n.id == id; }))
        throw std::runtime_error("cycle detected at " + id);
    chain.push_back(read_nuclide(env, id));
    if (!chain.back().stable && !chain.back().daughter.empty()) append_chain(env, chain.back().daughter, chain);
}

void print_nuclide(const Nuclide& n) {
    std::cout << n.id << " — " << n.label << '\n'
              << "  atomic mass: " << n.atomic_mass_g_mol << " g/mol\n"
              << "  half-life:   " << (n.stable ? std::string("stable") : std::to_string(n.half_life_s) + " s") << '\n'
              << "  decay:       " << n.mode << (n.daughter.empty() ? "" : " -> " + n.daughter) << '\n'
              << "  Q value:     " << n.energy_mev << " MeV\n";
}

void usage() {
    std::cout << "Usage: nuclide-atlas [--scenario FILE] [--output FILE] [--data DIR] [--isotope ID] [--list]\n";
}

} // namespace

int main(int argc, char** argv) {
    try {
        fs::path data_dir = SNT_NUCLEAR_DATA_DIR;
        fs::path scenario = SNT_NUCLEAR_DEFAULT_SCENARIO;
        fs::path output;
        std::string inspect;
        bool list = false;
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            auto value = [&](const char* name) -> std::string {
                if (++i >= argc) throw std::runtime_error(std::string("missing value for ") + name);
                return argv[i];
            };
            if (arg == "--scenario") scenario = value("--scenario");
            else if (arg == "--output") output = value("--output");
            else if (arg == "--data") data_dir = value("--data");
            else if (arg == "--isotope") inspect = value("--isotope");
            else if (arg == "--list") list = true;
            else if (arg == "--help" || arg == "-h") { usage(); return 0; }
            else throw std::runtime_error("unknown option: " + arg);
        }

        data_dir = fs::absolute(data_dir);
        scenario = fs::absolute(scenario);
        const fs::path original_working_directory = fs::current_path();
        // DIPL source paths are evaluated by SNT relative to the working directory.
        // Make the database self-contained even when the executable is launched from build/.
        fs::current_path(data_dir);
        snt::dip::DIP dip;
        dip.add_file(data_dir / "nuclear.dip", {}, true);
        if (inspect.empty() && !list) dip.add_file(scenario, {}, true);
        const Environment env = dip.parse();
        fs::current_path(original_working_directory);
        if (list) {
            for (const auto& id : {"U238", "Th234", "Pa234m", "U234", "Th230", "Ra226", "Rn222", "Po218", "Pb214", "Bi214", "Po214", "Pb210", "Bi210", "Po210", "Pb206", "Pu239", "U235", "Pb207"})
                std::cout << id << '\n';
            return 0;
        }
        if (!inspect.empty()) { print_nuclide(read_nuclide(env, inspect)); return 0; }

        const std::string parent = env["simulation.sample.isotope"].as<std::string>();
        const double sample_mass_g = quantity_as(env, "simulation.sample.mass", "g");
        const double duration_s = quantity_as(env, "simulation.duration", "s");
        const std::size_t points = static_cast<std::size_t>(env["simulation.points"].as<int64_t>());
        const std::string time_grid = env["simulation.time_grid"].as<std::string>();
        const double minimum_time_s = quantity_as(env, "simulation.minimum_time", "s");
        if (minimum_time_s >= duration_s) throw std::runtime_error("simulation.minimum_time must be shorter than simulation.duration");
        const double avogadro = quantity_as(env, "nuclear.constants.avogadro", "1/mol");
        std::vector<Nuclide> chain;
        append_chain(env, parent, chain);
        std::vector<double> inventory(chain.size());
        inventory[0] = sample_mass_g / chain[0].atomic_mass_g_mol * avogadro;
        const auto terms = bateman_terms(chain, inventory[0]);
        if (output.empty()) output = env["simulation.output.csv"].as<std::string>();
        std::ofstream csv(output);
        if (!csv) throw std::runtime_error("cannot open output: " + output.string());
        csv << "time_s,time_yr,total_activity_bq";
        for (const auto& n : chain) csv << ',' << n.id << "_atoms," << n.id << "_activity_bq";
        csv << '\n' << std::setprecision(12);
        const double seconds_per_year = 31557600.0;
        for (std::size_t p = 0; p < points; ++p) {
            double t = duration_s * static_cast<double>(p) / static_cast<double>(points - 1);
            if (time_grid == "logarithmic" && p > 0) {
                if (points == 2) t = duration_s;
                else {
                    const double fraction = static_cast<double>(p - 1) / static_cast<double>(points - 2);
                    t = minimum_time_s * std::exp(fraction * std::log(duration_s / minimum_time_s));
                }
            }
            const auto values = inventory_at(terms, t);
            double activity = 0.0;
            for (std::size_t i = 0; i < values.size(); ++i) activity += chain[i].lambda() * std::max(0.0, values[i]);
            csv << t << ',' << t / seconds_per_year << ',' << activity;
            for (std::size_t i = 0; i < values.size(); ++i) csv << ',' << std::max(0.0, values[i]) << ',' << chain[i].lambda() * std::max(0.0, values[i]);
            csv << '\n';
        }
        std::cout << "Nuclide Atlas\n"
                  << "  scenario: " << env["simulation.title"].as<std::string>() << '\n'
                  << "  chain:    ";
        for (std::size_t i = 0; i < chain.size(); ++i) std::cout << (i ? " -> " : "") << chain[i].id;
        std::cout << "\n  atoms at t=0: " << std::setprecision(6) << inventory[0]
                  << "\n  output: " << output << '\n';
    } catch (const std::exception& error) {
        std::cerr << "nuclide-atlas: " << error.what() << '\n';
        return 1;
    }
}
