cmake_minimum_required(VERSION 3.22)

project(nuclide_atlas VERSION 0.1.0 LANGUAGES CXX)

set(SNT_NUCLEAR_DEFAULT_SCENARIO "${CMAKE_CURRENT_SOURCE_DIR}/scenarios/u238_age.dip" CACHE FILEPATH "Default DIPL scenario")

# SNT is a first-class package dependency. Its package config exposes both the
# C++ targets and a host `snt` CLI; the latter evaluates build.dip at configure
# time. This deliberately keeps the CMake feature visible and reproducible.
find_package(snt CONFIG REQUIRED)
snt_dip_get(
    FILE "${CMAKE_CURRENT_SOURCE_DIR}/cmake/build.dip"
    PATH build.cxx_standard
    TYPE INTEGER
    OUTPUT SNT_NUCLEAR_CXX_STANDARD
)
snt_dip_get(
    FILE "${CMAKE_CURRENT_SOURCE_DIR}/cmake/build.dip"
    PATH build.enable_plot_target
    TYPE BOOL
    OUTPUT SNT_NUCLEAR_ENABLE_PLOT_TARGET
)

add_executable(nuclide-atlas src/main.cpp)
target_compile_features(nuclide-atlas PRIVATE cxx_std_${SNT_NUCLEAR_CXX_STANDARD})
target_link_libraries(nuclide-atlas PRIVATE snt-dip snt-puq)
target_compile_definitions(nuclide-atlas PRIVATE
    SNT_NUCLEAR_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data"
    SNT_NUCLEAR_DEFAULT_SCENARIO="${SNT_NUCLEAR_DEFAULT_SCENARIO}"
)

install(TARGETS nuclide-atlas RUNTIME DESTINATION bin)
install(DIRECTORY data scenarios DESTINATION share/nuclide-atlas)
install(DIRECTORY python DESTINATION share/nuclide-atlas)

include(CTest)
if(BUILD_TESTING)
    add_test(
        NAME u238_inventory
        COMMAND nuclide-atlas --scenario "${CMAKE_CURRENT_SOURCE_DIR}/scenarios/u238_age.dip"
                --output "${CMAKE_CURRENT_BINARY_DIR}/u238_inventory.csv"
    )
endif()

find_package(Python3 COMPONENTS Interpreter QUIET)
if(SNT_NUCLEAR_ENABLE_PLOT_TARGET AND Python3_Interpreter_FOUND)
    add_custom_target(plot
        COMMAND $<TARGET_FILE:nuclide-atlas> --output "${CMAKE_CURRENT_BINARY_DIR}/inventory.csv"
        COMMAND "${CMAKE_COMMAND}" -E make_directory "${CMAKE_CURRENT_BINARY_DIR}/.plot-cache"
        COMMAND "${CMAKE_COMMAND}" -E env
                "MPLCONFIGDIR=${CMAKE_CURRENT_BINARY_DIR}/.plot-cache"
                "XDG_CACHE_HOME=${CMAKE_CURRENT_BINARY_DIR}/.plot-cache"
                "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/python/plot_inventory.py"
                "${CMAKE_CURRENT_BINARY_DIR}/inventory.csv"
                "${CMAKE_CURRENT_BINARY_DIR}/inventory.png"
        DEPENDS nuclide-atlas
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        COMMENT "Plotting decay-chain inventory"
        VERBATIM
    )
endif()
