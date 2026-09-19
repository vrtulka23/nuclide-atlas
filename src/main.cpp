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
#include <limits>
#include <map>
#include <numbers>
#include <numeric>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
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
        double absolute_sum = 0.0;
        for (const auto& [lambda, coefficient] : isotope_terms) {
            const double contribution = coefficient * std::exp(-lambda * time_s);
            atoms += contribution;
            absolute_sum += std::abs(contribution);
        }
        // Near t=0, a distant descendant is the difference of large nearly
        // equal Bateman terms.  Values below this floating-point cancellation
        // bound have no physical significance; report an honest zero instead
        // of a spurious trace inventory.
        const double cancellation_bound = 64.0 * std::numeric_limits<double>::epsilon() * absolute_sum;
        inventory.push_back(std::abs(atoms) <= cancellation_bound ? 0.0 : std::max(0.0, atoms));
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
    const std::string root = "nuclear.isotopes[" + id + "]";
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

#if SNT_NUCLEAR_ENABLE_ACTIVATION
using Matrix = std::vector<std::vector<double>>;

struct ActivationNetwork {
    std::vector<Nuclide> nuclides;
    // This absorbing state receives omitted radioactive branches.  It closes
    // the rate matrix, so numerical propagation can enforce atom balance
    // without pretending that omitted branches are still reported isotopes.
    std::size_t loss_index{};
};

Matrix identity_matrix(std::size_t size) {
    Matrix result(size, std::vector<double>(size));
    for (std::size_t i = 0; i < size; ++i) result[i][i] = 1.0;
    return result;
}

double matrix_norm(const Matrix& matrix) {
    double result = 0.0;
    for (const auto& row : matrix) {
        double sum = 0.0;
        for (const double value : row) sum += std::abs(value);
        result = std::max(result, sum);
    }
    return result;
}

Matrix multiply(const Matrix& left, const Matrix& right) {
    const std::size_t size = left.size();
    Matrix result(size, std::vector<double>(size));
    for (std::size_t i = 0; i < size; ++i)
        for (std::size_t k = 0; k < size; ++k) {
            const double coefficient = left[i][k];
            if (coefficient == 0.0) continue;
            for (std::size_t j = 0; j < size; ++j) result[i][j] += coefficient * right[k][j];
        }
    return result;
}

void normalize_stochastic_columns(Matrix& matrix) {
    for (std::size_t column = 0; column < matrix.size(); ++column) {
        double sum = 0.0;
        for (std::size_t row = 0; row < matrix.size(); ++row) {
            double& value = matrix[row][column];
            if (value < 0.0) {
                if (value < -1e-13) throw std::runtime_error("activation propagator lost positivity");
                value = 0.0;
            }
            sum += value;
        }
        if (!std::isfinite(sum) || sum <= 0.0) throw std::runtime_error("activation propagator lost atom balance");
        for (std::size_t row = 0; row < matrix.size(); ++row) matrix[row][column] /= sum;
    }
}

std::vector<double> multiply(const Matrix& matrix, const std::vector<double>& values) {
    std::vector<double> result(matrix.size());
    for (std::size_t i = 0; i < matrix.size(); ++i)
        for (std::size_t j = 0; j < values.size(); ++j) result[i] += matrix[i][j] * values[j];
    return result;
}

// Scaling and squaring makes the Taylor series well conditioned even when the
// network combines microsecond daughters with years of cooling.  The augmented
// generator has zero column sums, so its propagator has unit column sums;
// normalizing after each squaring keeps
// the finite-precision propagator positive and atom-conserving.
Matrix matrix_exponential(Matrix matrix, const double time_s) {
    const std::size_t size = matrix.size();
    for (auto& row : matrix)
        for (double& value : row) value *= time_s;
    const double norm = matrix_norm(matrix);
    const int squarings = norm > 0.5 ? static_cast<int>(std::ceil(std::log2(norm / 0.5))) : 0;
    const double scale = std::ldexp(1.0, -squarings);
    for (auto& row : matrix)
        for (double& value : row) value *= scale;

    Matrix result = identity_matrix(size);
    Matrix term = result;
    for (int order = 1; order <= 128; ++order) {
        term = multiply(term, matrix);
        for (auto& row : term)
            for (double& value : row) value /= static_cast<double>(order);
        for (std::size_t i = 0; i < size; ++i)
            for (std::size_t j = 0; j < size; ++j) result[i][j] += term[i][j];
        if (matrix_norm(term) <= 1e-15 * std::max(1.0, matrix_norm(result))) break;
    }
    normalize_stochastic_columns(result);
    for (int step = 0; step < squarings; ++step) {
        result = multiply(result, result);
        normalize_stochastic_columns(result);
    }
    return result;
}

ActivationNetwork activation_network(
    const Environment& env, const std::string& sample, const std::string& capture_product
) {
    ActivationNetwork network;
    std::unordered_map<std::string, std::size_t> index;
    const auto visit = [&](const auto& self, const std::string& id) -> void {
        if (index.contains(id)) return;
        index.emplace(id, network.nuclides.size());
        network.nuclides.push_back(read_nuclide(env, id));
        const auto& nuclide = network.nuclides.back();
        if (!nuclide.stable && !nuclide.daughter.empty() && nuclide.daughter != "none") self(self, nuclide.daughter);
    };
    visit(visit, sample);
    visit(visit, capture_product);
    network.loss_index = network.nuclides.size();
    return network;
}

Matrix transition_matrix(
    const ActivationNetwork& network,
    const std::string& capture_target,
    const std::string& capture_product,
    const double capture_rate_s
) {
    const std::size_t size = network.nuclides.size() + 1;
    Matrix result(size, std::vector<double>(size));
    std::unordered_map<std::string, std::size_t> index;
    for (std::size_t i = 0; i < network.nuclides.size(); ++i) index.emplace(network.nuclides[i].id, i);
    for (std::size_t parent = 0; parent < network.nuclides.size(); ++parent) {
        const auto& nuclide = network.nuclides[parent];
        const double decay_rate = nuclide.lambda();
        result[parent][parent] -= decay_rate;
        const bool has_declared_daughter = !nuclide.stable && !nuclide.daughter.empty() && nuclide.daughter != "none";
        if (has_declared_daughter) {
            result[index.at(nuclide.daughter)][parent] += nuclide.branching * nuclide.lambda();
        }
        if (!nuclide.stable) result[network.loss_index][parent] += (has_declared_daughter ? 1.0 - nuclide.branching : 1.0) * decay_rate;
    }
    if (capture_rate_s > 0.0) {
        const auto source = index.at(capture_target);
        result[source][source] -= capture_rate_s;
        result[index.at(capture_product)][source] += capture_rate_s;
    }
    return result;
}

std::vector<double> propagate_activation(const Matrix& transition, const double time_s, const std::vector<double>& initial) {
    std::vector<double> values = multiply(matrix_exponential(transition, time_s), initial);
    const double expected_atoms = std::accumulate(initial.begin(), initial.end(), 0.0);
    double actual_atoms = 0.0;
    for (double& value : values) {
        if (value < 0.0) {
            if (value < -1e-12 * expected_atoms) throw std::runtime_error("activation propagation produced a negative inventory");
            value = 0.0;
        }
        actual_atoms += value;
    }
    if (std::abs(actual_atoms - expected_atoms) > 2e-12 * expected_atoms)
        throw std::runtime_error("activation propagation failed its atom-balance check");
    return values;
}

void run_activation(const Environment& env, fs::path output) {
    const std::string sample = env["activation.sample.isotope"].as<std::string>();
    const double mass_g = quantity_as(env, "activation.sample.mass", "g");
    const double irradiation_s = quantity_as(env, "activation.irradiation.duration", "s");
    const std::size_t irradiation_points = static_cast<std::size_t>(env["activation.irradiation.points"].as<int64_t>());
    const double flux = quantity_as(env, "activation.irradiation.neutron_flux", "1/(cm2*s)");
    const std::string capture_target = env["activation.irradiation.capture.target"].as<std::string>();
    const std::string capture_product = env["activation.irradiation.capture.product"].as<std::string>();
    const double cross_section = quantity_as(env, "activation.irradiation.capture.cross_section", "cm2");
    const double cooldown_s = quantity_as(env, "activation.cooldown.duration", "s");
    const std::size_t cooldown_points = static_cast<std::size_t>(env["activation.cooldown.points"].as<int64_t>());
    const double avogadro = quantity_as(env, "nuclear.constants.avogadro", "1/mol");
    const double seconds_per_year = quantity_as(env, "nuclear.constants.seconds_per_year", "s");
    const double joules_per_mev = quantity_as(env, "nuclear.constants.joules_per_mev", "J/MeV");
    const bool include_activity = env["activation.output.include_activity"].as<bool>();
    const bool include_q_power = env["activation.output.include_q_power"].as<bool>();
    const auto network = activation_network(env, sample, capture_product);
    const auto sample_it = std::find_if(network.nuclides.begin(), network.nuclides.end(), [&](const Nuclide& n) { return n.id == sample; });
    if (sample_it == network.nuclides.end()) throw std::runtime_error("activation sample is absent from its network");
    if (std::none_of(network.nuclides.begin(), network.nuclides.end(), [&](const Nuclide& n) { return n.id == capture_target; }))
        throw std::runtime_error("activation capture target is absent from its network");
    const double capture_rate_s = flux * cross_section;
    const Matrix during_irradiation = transition_matrix(network, capture_target, capture_product, capture_rate_s);
    const Matrix during_cooldown = transition_matrix(network, capture_target, capture_product, 0.0);
    std::vector<double> initial(network.nuclides.size() + 1);
    initial[static_cast<std::size_t>(std::distance(network.nuclides.begin(), sample_it))] = mass_g / sample_it->atomic_mass_g_mol * avogadro;
    const auto end_of_irradiation = propagate_activation(during_irradiation, irradiation_s, initial);
    if (output.empty()) output = env["activation.output.csv"].as<std::string>();
    std::ofstream csv(output);
    if (!csv) throw std::runtime_error("cannot open output: " + output.string());
    csv << "time_s,time_yr,phase,total_atoms";
    if (include_activity) csv << ",total_activity_bq";
    if (include_q_power) csv << ",total_q_power_w";
    csv << ",capture_reactions_per_s,untracked_loss_atoms";
    for (const auto& nuclide : network.nuclides) {
        csv << ',' << nuclide.id << "_atoms";
        if (include_activity) csv << ',' << nuclide.id << "_activity_bq";
    }
    csv << '\n' << std::setprecision(12);
    const auto write_snapshot = [&](const double time_s, const std::string_view phase, std::vector<double> values) {
        double total_atoms = 0.0;
        double total_activity = 0.0;
        double q_power = 0.0;
        for (std::size_t i = 0; i < network.nuclides.size(); ++i) {
            const double activity = network.nuclides[i].lambda() * values[i];
            total_atoms += values[i];
            total_activity += activity;
            q_power += activity * network.nuclides[i].energy_mev * joules_per_mev;
        }
        double reactions_per_second = 0.0;
        if (phase == "irradiation") {
            const auto target = std::find_if(network.nuclides.begin(), network.nuclides.end(), [&](const Nuclide& n) { return n.id == capture_target; });
            reactions_per_second = capture_rate_s * values[static_cast<std::size_t>(std::distance(network.nuclides.begin(), target))];
        }
        csv << time_s << ',' << time_s / seconds_per_year << ',' << phase << ',' << total_atoms;
        if (include_activity) csv << ',' << total_activity;
        if (include_q_power) csv << ',' << q_power;
        csv << ',' << reactions_per_second << ',' << values[network.loss_index];
        for (std::size_t i = 0; i < network.nuclides.size(); ++i) {
            csv << ',' << values[i];
            if (include_activity) csv << ',' << network.nuclides[i].lambda() * values[i];
        }
        csv << '\n';
    };
    for (std::size_t point = 0; point < irradiation_points; ++point) {
        const double time_s = irradiation_s * static_cast<double>(point) / static_cast<double>(irradiation_points - 1);
        write_snapshot(time_s, "irradiation", propagate_activation(during_irradiation, time_s, initial));
    }
    for (std::size_t point = 1; point < cooldown_points; ++point) {
        const double elapsed_s = cooldown_s * static_cast<double>(point) / static_cast<double>(cooldown_points - 1);
        write_snapshot(irradiation_s + elapsed_s, "cooldown", propagate_activation(during_cooldown, elapsed_s, end_of_irradiation));
    }
    std::cout << "Nuclide Atlas activation\n"
              << "  scenario: " << env["activation.title"].as<std::string>() << '\n'
              << "  reaction: " << capture_target << " (n,gamma) " << capture_product << '\n'
              << "  flux:     " << flux << " 1/(cm2*s)\n"
              << "  output:   " << output << '\n';
}
#endif

void print_nuclide(const Nuclide& n) {
    std::cout << n.id << " — " << n.label << '\n'
              << "  atomic mass: " << n.atomic_mass_g_mol << " g/mol\n"
              << "  half-life:   " << (n.stable ? std::string("stable") : std::to_string(n.half_life_s) + " s") << '\n'
              << "  decay:       " << n.mode << (n.daughter.empty() ? "" : " -> " + n.daughter) << '\n'
              << "  Q value:     " << n.energy_mev << " MeV\n";
}

void usage() {
    std::cout << "Usage: nuclide-atlas [--scenario FILE] [--output FILE] [--data DIR] [--isotope ID] [--list]";
#if SNT_NUCLEAR_ENABLE_ACTIVATION
    std::cout << " [--activation]";
#endif
    std::cout << '\n';
}

} // namespace

int main(int argc, char** argv) {
    try {
        fs::path data_dir = SNT_NUCLEAR_DATA_DIR;
        fs::path scenario = SNT_NUCLEAR_DEFAULT_SCENARIO;
        fs::path output;
        std::string inspect;
        bool list = false;
        bool activation = false;
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
#if SNT_NUCLEAR_ENABLE_ACTIVATION
            else if (arg == "--activation") activation = true;
#else
            else if (arg == "--activation") throw std::runtime_error("activation support was disabled in dip/build.dip");
#endif
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
            std::vector<std::string> ids;
            for (const auto& [id, _] : env["nuclear.isotopes"].items()) ids.push_back(id);
            std::sort(ids.begin(), ids.end());
            for (const auto& id : ids)
                std::cout << id << '\n';
            return 0;
        }
        if (!inspect.empty()) { print_nuclide(read_nuclide(env, inspect)); return 0; }

#if SNT_NUCLEAR_ENABLE_ACTIVATION
        if (activation) {
            run_activation(env, output);
            return 0;
        }
#endif

        const std::string parent = env["simulation.sample.isotope"].as<std::string>();
        const double sample_mass_g = quantity_as(env, "simulation.sample.mass", "g");
        const double duration_s = quantity_as(env, "simulation.duration", "s");
        const std::size_t points = static_cast<std::size_t>(env["simulation.points"].as<int64_t>());
        const std::string time_grid = env["simulation.time_grid"].as<std::string>();
        const double minimum_time_s = quantity_as(env, "simulation.minimum_time", "s");
        if (minimum_time_s >= duration_s) throw std::runtime_error("simulation.minimum_time must be shorter than simulation.duration");
        const double avogadro = quantity_as(env, "nuclear.constants.avogadro", "1/mol");
        const double seconds_per_year = quantity_as(env, "nuclear.constants.seconds_per_year", "s");
        const double joules_per_mev = quantity_as(env, "nuclear.constants.joules_per_mev", "J/MeV");
        const bool include_activity = env["simulation.output.include_activity"].as<bool>();
        const bool include_q_power = env["simulation.output.include_q_power"].as<bool>();
        std::vector<Nuclide> chain;
        append_chain(env, parent, chain);
        std::vector<double> inventory(chain.size());
        inventory[0] = sample_mass_g / chain[0].atomic_mass_g_mol * avogadro;
        const auto terms = bateman_terms(chain, inventory[0]);
        if (output.empty()) output = env["simulation.output.csv"].as<std::string>();
        std::ofstream csv(output);
        if (!csv) throw std::runtime_error("cannot open output: " + output.string());
        csv << "time_s,time_yr";
        if (include_activity) csv << ",total_activity_bq";
        if (include_q_power) csv << ",total_q_power_w";
        for (const auto& n : chain) {
            csv << ',' << n.id << "_atoms";
            if (include_activity) csv << ',' << n.id << "_activity_bq";
        }
        csv << '\n' << std::setprecision(12);
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
            double q_power = 0.0;
            for (std::size_t i = 0; i < values.size(); ++i) {
                const double atoms = std::max(0.0, values[i]);
                const double isotope_activity = chain[i].lambda() * atoms;
                activity += isotope_activity;
                q_power += isotope_activity * chain[i].energy_mev * joules_per_mev;
            }
            csv << t << ',' << t / seconds_per_year;
            if (include_activity) csv << ',' << activity;
            if (include_q_power) csv << ',' << q_power;
            for (std::size_t i = 0; i < values.size(); ++i) {
                const double atoms = std::max(0.0, values[i]);
                csv << ',' << atoms;
                if (include_activity) csv << ',' << chain[i].lambda() * atoms;
            }
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
