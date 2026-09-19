#!/usr/bin/env bash
# Small local workflow helper for Nuclide Atlas.
set -euo pipefail

root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
build_dir="${root_dir}/build"
scenario="${root_dir}/scenarios/u238_age.dip"
prefix_path="${CMAKE_PREFIX_PATH:-${HOME}/.local}"
python_bin="${PYTHON:-python3}"

usage() {
    cat <<'EOF'
Usage: ./dev.sh [options]

  -b, --build       Configure the CMake build directory.
  -c, --compile     Compile Nuclide Atlas.
  -t, --test        Run CTest.
  -r, --run         Run the selected scenario and write build/inventory.csv.
  -p, --plot        Generate build/inventory.png (requires matplotlib).
  -s, --scenario    Scenario DIPL file for --run (default: scenarios/u238_age.dip).
      --clean       Remove this project's build directory.
  -h, --help        Show this help.

Examples:
  ./dev.sh -bc              # configure and compile
  ./dev.sh -bctr            # configure, compile, test, then run
  ./dev.sh --clean -bc -r   # clean rebuild and run

Environment:
  CMAKE_PREFIX_PATH  SciNumTools3 install prefix (default: $HOME/.local)
  PYTHON             Python interpreter used by the optional plot target
EOF
}

configure() {
    cmake -S "${root_dir}" -B "${build_dir}" -G Ninja \
        -DCMAKE_PREFIX_PATH="${prefix_path}" \
        -DPython3_EXECUTABLE="${python_bin}"
}

compile() { cmake --build "${build_dir}" --parallel; }
test_project() { ctest --test-dir "${build_dir}" --output-on-failure; }
run_project() { "${build_dir}/nuclide-atlas" --scenario "${scenario}" --output "${build_dir}/inventory.csv"; }
plot_project() { cmake --build "${build_dir}" --target plot; }

do_build=false do_compile=false do_test=false do_run=false do_plot=false do_clean=false
while (($#)); do
    case "$1" in
        -[bctrp]*)
            flags="${1#-}"
            for ((index = 0; index < ${#flags}; ++index)); do
                case "${flags:index:1}" in
                    b) do_build=true ;;
                    c) do_compile=true ;;
                    t) do_test=true ;;
                    r) do_run=true ;;
                    p) do_plot=true ;;
                    *) echo "Unknown short option: -${flags:index:1}" >&2; exit 2 ;;
                esac
            done
            ;;
        -b|--build) do_build=true ;;
        -c|--compile) do_compile=true ;;
        -t|--test) do_test=true ;;
        -r|--run) do_run=true ;;
        -p|--plot) do_plot=true ;;
        --clean) do_clean=true ;;
        -s|--scenario)
            shift
            (($#)) || { echo "Missing scenario after --scenario" >&2; exit 2; }
            scenario="$1"
            ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
    shift
done

if ! $do_build && ! $do_compile && ! $do_test && ! $do_run && ! $do_plot && ! $do_clean; then
    usage
    exit 2
fi

if $do_clean; then
    [[ -d "${build_dir}" ]] && rm -rf "${build_dir}"
fi

# Later stages imply their prerequisites.
if $do_compile || $do_test || $do_run || $do_plot; then do_build=true; fi
if $do_test || $do_run || $do_plot; then do_compile=true; fi

$do_build && configure
$do_compile && compile
$do_test && test_project
$do_run && run_project
$do_plot && plot_project
