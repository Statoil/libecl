if(ENABLE_PYTHON)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")
  include(init_python)
  init_python()

  configure_file(test_env.py.in
                 ${PROJECT_BINARY_DIR}/${PYTHON_INSTALL_PREFIX}/test_env.py)

  if(RST_DOC)
    add_subdirectory(docs)
  endif()
endif()

add_subdirectory(ecl)
