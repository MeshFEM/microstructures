################################################################################
include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
    set(MICRO_EXTRA_OPTIONS "GIT_CONFIG advice.detachedHead=false")
else()
    set(MICRO_EXTRA_OPTIONS "")
endif()

# Shortcut function
function(micro_download_project name)
    download_project(
        PROJ         ${name}
        SOURCE_DIR   ${MICRO_EXTERNAL}/${name}
        DOWNLOAD_DIR ${MICRO_EXTERNAL}/.cache/${name}
        QUIET
        ${MICRO_EXTRA_OPTIONS}
        ${ARGN}
    )
endfunction()

################################################################################

## MeshFEM
function(micro_download_meshfem)
    micro_download_project(MeshFEM
        GIT_REPOSITORY https://github.com/MeshFEM/MeshFEM.git
        GIT_TAG        c9fa2effae174e2783cf0de32454c6267bd1dee5
    )
endfunction()

## TBB
function(micro_download_tbb)
    if(MICRO_WITH_UBUNTU)
        micro_download_project(tbb
            GIT_REPOSITORY https://github.com/01org/tbb.git
            GIT_TAG        2019_U1
        )
    else()
        micro_download_project(tbb
            GIT_REPOSITORY https://github.com/wjakob/tbb.git
            GIT_TAG        344fa84f34089681732a54f5def93a30a3056ab9
        )
    endif()
endfunction()

## CGAL
function(micro_download_cgal)
    micro_download_project(cgal
        URL     https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-4.12/CGAL-4.12.tar.xz
        URL_MD5 b12fd24dedfa889a04abfaea565a88bd
    )
endfunction()

## Catch2
function(micro_download_catch)
    micro_download_project(Catch2
        URL     https://github.com/catchorg/Catch2/archive/v2.3.0.tar.gz
        URL_MD5 1fc90ff3b7b407b83057537f4136489e
    )
endfunction()

## CLI11
function(micro_download_cli11)
    micro_download_project(CLI11
        URL     https://github.com/CLIUtils/CLI11/archive/v1.6.0.tar.gz
        URL_MD5 c8e3dc70e3b7ebf6b01f618f7cdcc85f
    )
endfunction()

## nlopt
function(micro_download_nlopt)
    micro_download_project(nlopt
        GIT_REPOSITORY https://github.com/stevengj/nlopt.git
        GIT_TAG        37b74a8c2037eea5dc72fea7eeb9b850fa978913
    )
endfunction()

## libigl
function(micro_download_libigl)
    micro_download_project(libigl
        GIT_REPOSITORY https://github.com/libigl/libigl.git
        GIT_TAG        75d60e40a8edc6868571fbdca2e74f97d5dddab8
    )
endfunction()

## nanoflann
function(micro_download_nanoflann)
    micro_download_project(nanoflann
        GIT_REPOSITORY https://github.com/jlblancoc/nanoflann
        GIT_TAG v1.3.0
    )
endfunction()

## quickhull
function(micro_download_quickhull)
    micro_download_project(quickhull
        GIT_REPOSITORY https://github.com/akuukka/quickhull.git
        GIT_TAG 4f65e0801b8f60c9a97da2dadbe63c2b46397694
    )
endfunction()

## Sanitizers
function(micro_download_sanitizers)
    micro_download_project(sanitizers-cmake
        GIT_REPOSITORY https://github.com/arsenm/sanitizers-cmake.git
        GIT_TAG        99e159ec9bc8dd362b08d18436bd40ff0648417b
    )
endfunction()

## Cotire
function(micro_download_cotire)
    micro_download_project(cotire
        GIT_REPOSITORY https://github.com/sakra/cotire.git
        GIT_TAG        391bf6b7609e14f5976bd5247b68d63cbf8d4d12
    )
endfunction()

## unit test data
function(micro_download_microstructures_data)
    micro_download_project(microstructures_data
        GIT_REPOSITORY https://github.com/MeshFEM/microstructures_data.git
        GIT_TAG 30438a36558b88f79b1335f9f1003e0122a552e0
    )
endfunction()
