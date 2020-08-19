# - Find Dlib
# Find Dlib includes
#
#  DLIB_INCLUDE_DIR - where to find dlib/optimization.h, etc.
#  DLIB_FOUND       - True if Dlib found.


IF (DLIB_INCLUDE_DIR)
  # Already in cache, be silent
  SET (dlib_FIND_QUIETLY TRUE)
ENDIF (DLIB_INCLUDE_DIR)

FIND_PATH(DLIB_INCLUDE_DIR dlib/optimization.h
	HINTS
		${EXTERNAL_LIBS_PREFIX}
        ${DLIB_PREFIX_PATH}
	PATHS
        $ENV{DLIB_INC}
)

# handle the QUIETLY and REQUIRED arguments and set DLIB_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE (FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS (DLIB DEFAULT_MSG
    DLIB_INCLUDE_DIR)

MARK_AS_ADVANCED (DLIB_INCLUDE_DIR)
