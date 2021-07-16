if(ENABLE_PYTHON)
  set(PYTHON_SOURCES __init__.py ecl_type.py ecl_util.py)
  add_python_package("python.ecl" ${PYTHON_INSTALL_PREFIX}/ecl
                     "${PYTHON_SOURCES}" True)

  add_subdirectory(eclfile)
  add_subdirectory(grid)
  add_subdirectory(rft)
  add_subdirectory(gravimetry)
  add_subdirectory(summary)
  add_subdirectory(util)
  add_subdirectory(well)

  configure_file(
    ecl_lib_info_build.py.in
    ${PROJECT_BINARY_DIR}/${PYTHON_INSTALL_PREFIX}/ecl/__ecl_lib_info.py)
  configure_file(
    ecl_lib_info_install.py.in
    ${PROJECT_BINARY_DIR}/${PYTHON_INSTALL_PREFIX}/ecl_lib_info_install.py)
  install(
    FILES ${PROJECT_BINARY_DIR}/${PYTHON_INSTALL_PREFIX}/ecl_lib_info_install.py
    DESTINATION ${PYTHON_INSTALL_PREFIX}/ecl
    RENAME __ecl_lib_info.py)
endif()

if(NOT SKBUILD)
  return()
endif()

execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c
          "import pybind11; print(pybind11.get_cmake_dir())"
  OUTPUT_VARIABLE _tmp_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ECHO STDOUT)
list(APPEND CMAKE_PREFIX_PATH "${_tmp_dir}")

find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(_native MODULE _native/__init__.cpp _native/summary.cpp)
target_link_libraries(_native PRIVATE ecl-static)
target_include_directories(
  _native
  PRIVATE "${CMAKE_SOURCE_DIR}/lib/include"
          "${CMAKE_SOURCE_DIR}/lib/private-include"
          "${CMAKE_BINARY_DIR}/lib/include")
target_compile_definitions(_native PRIVATE VERSION_INFO=${PROJECT_VERSION})
