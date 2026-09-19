#!/usr/bin/env bash
# Small local workflow helper for Nuclide Atlas.
set -euo pipefail

root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
build_dir="${root_dir}/build"
scenario="${root_dir}/dip/scenarios/u238_age.dip"
venv_dir="${root_dir}/.venv"
requirements_file="${root_dir}/requirements.txt"
bootstrap_python="${PYTHON:-python3}"
python_bin="${venv_dir}/bin/python"

usage() {
    cat <<'EOF'
Usage: ./setup.sh [options]

  -b, --build       Configure the CMake build directory.
  -c, --compile     Compile Nuclide Atlas.
  -t, --test        Run CTest.
  -r, --run         Run the selected scenario and write build/inventory.csv.
  -p, --plot        Generate build/inventory.png.
  -a, --activation  Treat the selected scenario as a neutron-activation run.
  -d, --deuteron-activation
                    Treat the selected scenario as a one-group deuteron-reaction run.
  -s, --scenario    Scenario DIPL file for --run (default: dip/scenarios/u238_age.dip).
      --clean       Remove this project's build directory.
  -h, --help        Show this help.

Examples:
  ./setup.sh -bc              # configure and compile
  ./setup.sh -bctr            # configure, compile, test, then run
  ./setup.sh -ar -s dip/scenarios/u238_activation.dip
  ./setup.sh -dr -s dip/scenarios/h2_deuteron_activation.dip
  ./setup.sh --clean -bc -r   # clean rebuild and run

Environment:
  CMAKE_PREFIX_PATH  Optional SciNumTools3 prefix when it is not installed in
                     a standard system location
  PYTHON             Python interpreter used to create .venv (default: python3)
EOF
}

prepare_python() {
    if [[ ! -x "${python_bin}" ]]; then
        "${bootstrap_python}" -m venv "${venv_dir}"
    fi

    "${python_bin}" -m pip install --requirement "${requirements_file}"
}

configure() {
    local -a cmake_args=(
        -S "${root_dir}" -B "${build_dir}" -G Ninja
        -DPython3_EXECUTABLE="${python_bin}"
    )

    # Use CMake's normal system-prefix search unless an override is requested.
    if [[ -n "${CMAKE_PREFIX_PATH:-}" ]]; then
        cmake_args+=("-DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}")
    fi

    cmake "${cmake_args[@]}"
}

compile() { cmake --build "${build_dir}" --parallel; }
test_project() { ctest --test-dir "${build_dir}" --output-on-failure; }
run_project() {
    local -a run_args=(--scenario "${scenario}" --output "${build_dir}/inventory.csv")
    $do_activation && run_args+=(--activation)
    $do_deuteron_activation && run_args+=(--deuteron-activation)
    "${build_dir}/nuclide-atlas" "${run_args[@]}"
}
plot_project() {
    if $do_activation || $do_deuteron_activation; then
        $do_run || run_project
        cmake -E make_directory "${build_dir}/.plot-cache"
        cmake -E env "MPLCONFIGDIR=${build_dir}/.plot-cache" "XDG_CACHE_HOME=${build_dir}/.plot-cache" \
            "${python_bin}" "${root_dir}/python/plot_inventory.py" "${build_dir}/inventory.csv" "${build_dir}/inventory.png"
    else
        cmake --build "${build_dir}" --target plot
    fi
}

do_build=false do_compile=false do_test=false do_run=false do_plot=false do_clean=false do_activation=false do_deuteron_activation=false
while (($#)); do
    case "$1" in
        -[bctrdpa]*)
            flags="${1#-}"
            for ((index = 0; index < ${#flags}; ++index)); do
                case "${flags:index:1}" in
                    b) do_build=true ;;
                    c) do_compile=true ;;
                    t) do_test=true ;;
                    r) do_run=true ;;
                    p) do_plot=true ;;
                    a) do_activation=true ;;
                    d) do_deuteron_activation=true ;;
                    *) echo "Unknown short option: -${flags:index:1}" >&2; exit 2 ;;
                esac
            done
            ;;
        -b|--build) do_build=true ;;
        -c|--compile) do_compile=true ;;
        -t|--test) do_test=true ;;
        -r|--run) do_run=true ;;
        -p|--plot) do_plot=true ;;
        -a|--activation) do_activation=true ;;
        -d|--deuteron-activation) do_deuteron_activation=true ;;
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

if $do_activation && $do_deuteron_activation; then
    echo "Choose either --activation or --deuteron-activation" >&2
    exit 2
fi

if $do_clean; then
    [[ -d "${build_dir}" ]] && rm -rf "${build_dir}"
fi

# Later stages imply their prerequisites.
if $do_compile || $do_test || $do_run || $do_plot; then do_build=true; fi
if $do_test || $do_run || $do_plot; then do_compile=true; fi

# Keep the optional plotting dependency isolated from the system Python. CMake
# receives this interpreter, so its plot target uses the same environment.
$do_build && prepare_python
$do_build && configure
$do_compile && compile
$do_test && test_project
$do_run && run_project
$do_plot && plot_project
