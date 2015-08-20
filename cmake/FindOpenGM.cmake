FIND_PATH(OPENGM_INCLUDE_DIR
          opengm/opengm.hxx
          HINTS "~/usr/local/include/opengm"
          "~/usr/include/opengm"
          "$ENV{OPENGM_ROOT}/include"
		  "${CMAKE_CURRENT_LIST_DIR}/../dependencies/opengm/include"
)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(OpenGM DEFAULT_MSG OPENGM_INCLUDE_DIR)

MARK_AS_ADVANCED( OPENGM_INCLUDE_DIR )
