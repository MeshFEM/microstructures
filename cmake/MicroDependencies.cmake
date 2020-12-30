# Prepare dependencies
#
# For each third-party library, if the appropriate target doesn't exist yet,
# download it via external project, and add_subdirectory to build it alongside
# this project.

### Configuration
set(MICRO_ROOT     "${CMAKE_CURRENT_LIST_DIR}/..")
set(MICRO_EXTERNAL "${MICRO_ROOT}/3rdparty")

# Download and update 3rdparty libraries
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
list(REMOVE_DUPLICATES CMAKE_MODULE_PATH)
include(MicroDownloadExternal)

################################################################################
# Required libraries
################################################################################

# TBB library; must be brought in before MeshFEM to override! We need tbbmalloc,
# which MeshFEM chooses not to build.
# There are also some segfaults on shutdown with TBB 2017 (the version in wjakob's
# repository), so we need to use a more recent version of TBB.

# While wjakob has been updated to TBB 2019 recently, it seems to hang Travis
# at the linking stage for some reason, so we'll just use the upstream version
# for now.
if(MICRO_WITH_UBUNTU AND NOT TARGET tbb::tbb)
    micro_download_tbb()
    list(APPEND CMAKE_MODULE_PATH ${MICRO_EXTERNAL}/tbb/cmake)
    include(TBBBuild)
    tbb_build(TBB_ROOT ${MICRO_EXTERNAL}/tbb CONFIG_DIR TBB_DIR)
    find_package(TBB REQUIRED tbb tbbmalloc)
    add_library(tbb_tbb INTERFACE)
    add_library(tbb::tbb ALIAS tbb_tbb)
    target_link_libraries(tbb_tbb INTERFACE TBB::tbb TBB::tbbmalloc)
endif()

if(NOT TARGET tbb::tbb)
    set(TBB_BUILD_STATIC ON CACHE BOOL " " FORCE)
    set(TBB_BUILD_SHARED OFF CACHE BOOL " " FORCE)
    set(TBB_BUILD_TBBMALLOC ON CACHE BOOL " " FORCE)
    set(TBB_BUILD_TBBMALLOC_PROXY OFF CACHE BOOL " " FORCE)
    set(TBB_BUILD_TESTS OFF CACHE BOOL " " FORCE)
    set(TBB_NO_DATE ON CACHE BOOL " " FORCE)

    micro_download_tbb()
    add_subdirectory(${MICRO_EXTERNAL}/tbb tbb EXCLUDE_FROM_ALL)
    set_property(TARGET tbb_static tbb_def_files PROPERTY FOLDER "dependencies")
    if(NOT MSVC)
        set_target_properties(tbb_static PROPERTIES COMPILE_FLAGS "-Wno-implicit-fallthrough -Wno-missing-field-initializers -Wno-unused-parameter -Wno-keyword-macro")
    endif()

    add_library(tbb_tbb INTERFACE)
    target_include_directories(tbb_tbb SYSTEM INTERFACE ${MICRO_EXTERNAL}/tbb/include)
    target_link_libraries(tbb_tbb INTERFACE tbb_static tbbmalloc_static)
    add_library(tbb::tbb ALIAS tbb_tbb)

    micro_target_hide_warnings(tbb_tbb tbb_static tbbmalloc_static)
endif()

if(NOT TARGET micro::tbb)
    add_library(micro_tbb INTERFACE)
    if(MICRO_WITH_TBB)
        target_link_libraries(micro_tbb INTERFACE tbb::tbb)
        target_compile_definitions(micro_tbb INTERFACE -DMICRO_WITH_TBB)
    endif()
    add_library(micro::tbb ALIAS micro_tbb)
endif()

# MeshFEM library
if(NOT TARGET MeshFEM)
    micro_download_meshfem()
    option(MESHFEM_WITH_CERES "Compile MeshFEM with Ceres" ${MICRO_WITH_CERES})
    add_subdirectory(${MICRO_EXTERNAL}/MeshFEM MeshFEM)
endif()

# CLI11 library
if(NOT TARGET CLI11::CLI11)
    micro_download_cli11()
    add_library(CLI11 INTERFACE)
    target_include_directories(CLI11 SYSTEM INTERFACE ${MICRO_EXTERNAL}/CLI11/include)
    add_library(CLI11::CLI11 ALIAS CLI11)
endif()

# CGAL library
if(MICRO_WITH_CGAL AND NOT TARGET CGAL::CGAL)
    micro_download_cgal()
    set(CGAL_DIR ${MICRO_EXTERNAL}/cgal)
    find_package(CGAL CONFIG REQUIRED COMPONENTS PATHS ${CGAL_DIR} NO_DEFAULT_PATH)
endif()

# Nanoflann
if(NOT TARGET nanoflann::nanoflann)
    micro_download_nanoflann()
    add_library(nanoflann INTERFACE)
    add_library(nanoflann::nanoflann ALIAS nanoflann)
    target_include_directories(nanoflann INTERFACE ${MICRO_EXTERNAL}/nanoflann/include)
endif()

# Quickhull
if(NOT TARGET quickhull::quickhull)
    micro_download_quickhull()
    file(GLOB QUICKHULL_SOURCES "${MICRO_EXTERNAL}/quickhull/*.cpp")
    add_library(quickhull ${QUICKHULL_SOURCES})
    add_library(quickhull::quickhull ALIAS quickhull)
    target_include_directories(quickhull SYSTEM PUBLIC
        ${MICRO_EXTERNAL}/quickhull/
        ${MICRO_EXTERNAL}/quickhull/Structs
    )
endif()

# Catch2
if(NOT TARGET Catch2::Catch2 AND (CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR))
    micro_download_catch()
    add_subdirectory(${MICRO_EXTERNAL}/Catch2)
    list(APPEND CMAKE_MODULE_PATH ${MICRO_EXTERNAL}/Catch2/contrib)
endif()

# Cotire
if(MICRO_WITH_COTIRE)
    micro_download_cotire()
endif()

# Unit test data
if(MICRO_TOPLEVEL_PROJECT AND NOT TARGET micro::data)
    micro_download_microstructures_data()
    add_library(micro_data INTERFACE)
    add_library(micro::data ALIAS micro_data)
    set(DATA_DIR "${MICRO_EXTERNAL}/microstructures_data/")
    target_compile_definitions(micro_data INTERFACE -DDATA_DIR=\"${DATA_DIR}\")
endif()

################################################################################
# Optional libraries
################################################################################

# Ceres library
if(NOT TARGET micro::ceres)
    # Target ceres::ceres should have been defined above if the option MICRO_WITH_CERES was given
    if(MICRO_WITH_CERES AND TARGET ceres::ceres)
        add_library(micro_ceres INTERFACE)
        add_library(micro::ceres ALIAS micro_ceres)
        target_compile_definitions(micro_ceres INTERFACE -DHAS_CERES)
        target_link_libraries(micro_ceres INTERFACE ceres::ceres)
    else()
        if(MICRO_WITH_CERES)
            message(STATUS "Google's ceres-solver not found; levenberg-marquardt disabled")
        endif()
        add_library(micro::ceres INTERFACE IMPORTED)
    endif()
endif()

# Dlib Library
if(NOT TARGET micro::dlib)
    find_package(Dlib QUIET)
    if(DLIB_FOUND)
        # Complains when Dlib is compiled in Release and the main app is compiled in Debug
        add_library(micro_dlib INTERFACE)
        target_include_directories(micro_dlib SYSTEM INTERFACE ${DLIB_INCLUDE_DIR})
        target_link_libraries(micro_dlib ${DLIB_LIBRARIES})
        target_compile_definitions(micro_dlib INTERFACE -DHAS_DLIB)
        add_library(micro::dlib ALIAS micro_dlib)
    else()
        message(STATUS "DLib not found; bfgs optimizer is disabled")
        add_library(micro::dlib INTERFACE IMPORTED)
    endif()
endif()

# Knitro library
if(NOT TARGET micro::knitro)
    find_package(Knitro QUIET)
    if(KNITRO_FOUND)
        message(STATUS "Found Knitro; active_set enabled")
        add_library(micro_knitro INTERFACE)
        add_library(micro::knitro ALIAS micro_knitro)
        target_link_libraries(micro_knitro INTERFACE knitro::knitro)
        target_compile_definitions(micro_knitro INTERFACE -DHAS_KNITRO)
    else()
        message(STATUS "Knitro not found; active_set disabled")
        add_library(micro::knitro INTERFACE IMPORTED)
    endif()
endif()

# libigl library
if(NOT TARGET micro::libigl)
    if(NOT TARGET igl::core)
        micro_download_libigl()
        find_package(LIBIGL QUIET)
    endif()
    if(LIBIGL_FOUND)
        add_library(micro_libigl INTERFACE)
        target_link_libraries(micro_libigl INTERFACE igl::core)
        target_compile_definitions(micro_libigl INTERFACE -DHAS_LIBIGL)
        add_library(micro::libigl ALIAS micro_libigl)
    else()
        add_library(micro::libigl INTERFACE IMPORTED)
    endif()
endif()

# NLopt library
if(NOT TARGET micro::nlopt)
    micro_download_nlopt()
    find_package(NLopt QUIET)
    if(NLOPT_FOUND)
        add_library(micro_nlopt INTERFACE)
        target_link_libraries(micro_nlopt INTERFACE ${NLOPT_LIBRARIES})
        target_compile_definitions(micro_nlopt INTERFACE -DHAS_NLOPT)
        add_library(micro::nlopt ALIAS micro_nlopt)
    else()
        message(STATUS "NLopt not found; slsqp disabled")
        add_library(micro::nlopt INTERFACE IMPORTED)
    endif()
endif()

# PyMesh Wires library
if(NOT TARGET micro::pymesh)
    # find_package(PyMesh QUIET)
    if(PYMESH_FOUND)
        add_library(micro_pymesh INTERFACE)
        target_link_libraries(micro_pymesh INTERFACE PyMesh::core PyMesh::wires)
        target_compile_definitions(micro_pymesh INTERFACE -DHAS_PYMESH)
        add_library(micro::pymesh ALIAS micro_pymesh)
    else()
        message(STATUS "PyMesh not found; disabling James' Inflator")
        add_library(micro::pymesh INTERFACE IMPORTED)
    endif()
endif()

# Sanitizers
if(MICRO_WITH_SANITIZERS)
    if(NOT COMMAND add_sanitizers)
        micro_download_sanitizers()
    endif()
    find_package(Sanitizers)
endif()
