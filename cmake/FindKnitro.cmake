# - Find Knitro
#
#  Searches for includes/libraries using environment variable $KNITRO_PATH, $KNITRO_DIR, $KNITRO_HOME or #KNITRO_ROOT
#
#  KNITRO_FOUND        - True if Knitro found.
#  knitro::knitro      - Imported target for Knitro.
#

if(KNITRO_INCLUDE_DIRS)
    # Already in cache, be silent
    set(knitro_FIND_QUIETLY TRUE)
endif(KNITRO_INCLUDE_DIRS)

find_path(KNITRO_INCLUDE_DIR knitro.h
	HINTS
        ENV KNITRO_PATH
        ENV KNITRO_DIR
        ENV KNITRO_HOME
        ENV KNITRO_ROOT
        ENV KNITRO_INC
        "/home/jdumas/.knitro"
    PATH_SUFFIXES
        include
)

find_library(KNITRO_LIBRARY NAMES knitro
	HINTS
        ENV KNITRO_PATH
        ENV KNITRO_DIR
        ENV KNITRO_HOME
        ENV KNITRO_ROOT
        ENV KNITRO_LIB
        "/home/jdumas/.knitro"
    PATH_SUFFIXES
        lib
)

# handle the QUIETLY and REQUIRED arguments and set KNITRO_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(KNITRO DEFAULT_MSG KNITRO_LIBRARY KNITRO_INCLUDE_DIR)

if(KNITRO_FOUND)
    add_library(knitro::knitro UNKNOWN IMPORTED)
    message("include dirs: ${KNITRO_INCLUDE_DIR}")
    set_target_properties(knitro::knitro PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
        IMPORTED_LOCATION "${KNITRO_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${KNITRO_INCLUDE_DIR};${KNITRO_INCLUDE_DIR}/../examples/C++/include"
    )
endif()

mark_as_advanced(KNITRO_LIBRARY KNITRO_INCLUDE_DIR)
