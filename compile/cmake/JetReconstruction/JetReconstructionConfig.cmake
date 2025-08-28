# Config file for the JetReconstruction.jl package
# Manualy adjusted from standard cmake generated config file 
get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

# - Create relocatable paths to headers.
# NOTE: Do not strictly need paths as all usage requirements are encoded in
# the imported targets created later.
set_and_check(JetReconstruction_INCLUDE_DIR "${PACKAGE_PREFIX_DIR}/include")

# - Create path to installed read-only data files (e.g. yaml files)
set_and_check(JetReconstruction_DATA_DIR "${PACKAGE_PREFIX_DIR}/share/julia")

# - Include the targets file to create the imported targets that a client can
# link to (libraries) or execute (programs)
include("${CMAKE_CURRENT_LIST_DIR}/JetReconstructionTargets.cmake")

# print the default "Found:" message and check library location
include(FindPackageHandleStandardArgs)
get_property(TEST_JETRECONSTRUCTION_LIBRARY TARGET JetReconstruction::JetReconstruction PROPERTY LOCATION)
find_package_handle_standard_args(JetReconstruction DEFAULT_MSG CMAKE_CURRENT_LIST_FILE TEST_JETRECONSTRUCTION_LIBRARY)
