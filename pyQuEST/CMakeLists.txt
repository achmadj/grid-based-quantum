cmake_minimum_required(VERSION 3.5)

project(
    pyQuEST_Project
    LANGUAGES C CXX)

if(SKBUILD)
    message(STATUS "Building project " ${CMAKE_PROJECT_NAME} " using scikit-build")
endif()
message(STATUS "CMake version: ${CMAKE_VERSION}")
message(STATUS "CMake build type: ${CMAKE_BUILD_TYPE}")

# ==== QuEST ====
# We want the QuEST library to be shared between all modules of pyQuEST.
# Because building with static libraries would introduce bugs, we force
# building a shared lib regardless of the cache variable.
set(BUILD_SHARED_LIBS TRUE)
# Build the QuEST library.
add_subdirectory(QuEST/QuEST)
# To get proper exception/error handling, we must override `invalidQuESTInputError()`
# in QuEST, which we can do by adding `pyquest/quest_exception.c` to its sources.
target_include_directories(QuEST PUBLIC ${CMAKE_SOURCE_DIR}/pyquest)
target_sources(QuEST PRIVATE ${CMAKE_SOURCE_DIR}/pyquest/quest_exception.cpp)
# For Win32 `dll`s every exported symbol must usually be marked explicitly
# as such. Fortunately, CMake supports automatic export of all symbols (yay).
# On non-Windows systems, this is a NOP.
set_property(TARGET QuEST PROPERTY WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
# The shared library goes directly into the 'pyquest' folder, so it can
# easily be found by the modules.
install(TARGETS QuEST DESTINATION pyquest)

# ==== pyQuEST ====
if(APPLE)
    # The @rpath is embedded in the Python modules and tells the loader where to
    # look for dynamic libraries. In our case, it will be in the same directory as
    # the modules (because QuEST is installed into the 'pyquest' directory, see
    # above). On macOS, the module path is available as the variable '@loader_path'.
    set(CMAKE_INSTALL_RPATH "@loader_path")
elseif(UNIX)
    # The same as for macOS goes for other Unix systems, but the path of the
    # loader module is stored in '$ORIGIN'.
    set(CMAKE_INSTALL_RPATH "$ORIGIN")
endif()
# pyquest installs all its targets automatically in the 'pyquest' folder,
# so this call is all that is needed to add it.
add_subdirectory(pyquest)
