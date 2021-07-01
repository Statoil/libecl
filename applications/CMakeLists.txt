project(libecl-applications)

function(target_link_libecl target)
    target_link_libraries(${target} ecl)
    if(SKBUILD)
        if(APPLE)
            set_target_properties(${target} PROPERTIES INSTALL_RPATH "$executable_path/../.libs")
        else()
            set_target_properties(${target} PROPERTIES INSTALL_RPATH "$ORIGIN/../.libs")
        endif()
    endif()
endfunction()

if (BUILD_APPLICATIONS)
    add_executable(sum_write   ecl/sum_write.cpp)
    add_executable(make_grid   ecl/make_grid.c)
    add_executable(grdecl_grid ecl/grdecl_grid.c)
    add_executable(summary     ecl/view_summary.cpp)
    add_executable(kw_extract  ecl/kw_extract.cpp)
    target_link_libecl(sum_write)
    target_link_libecl(make_grid)
    target_link_libecl(grdecl_grid)
    target_link_libecl(summary)
    target_link_libecl(kw_extract)

    list(APPEND apps make_grid grdecl_grid summary kw_extract)

    foreach (app ecl_pack
                 ecl_unpack
                 grid_info
                 grid_dump
                 grid_dump_ascii
                 select_test
                 load_test
            )
        add_executable(${app} ecl/${app}.c)
        target_link_libecl(${app})
        list(APPEND apps ${app})
        if (ERT_LINUX)
            # The .x extension creates problems on windows
            set_target_properties(${app} PROPERTIES SUFFIX ".x")
        endif ()
    endforeach ()


    if (ERT_LINUX)
        foreach (app convert grdecl_test kw_list)
            add_executable(${app} ecl/${app}.c)
            target_link_libecl(${app})
            # The .x extension creates problems on windows
            set_target_properties(${app} PROPERTIES SUFFIX ".x")
            list(APPEND apps ${app})
        endforeach ()

        set_target_properties(summary PROPERTIES SUFFIX ".x")
    endif ()


    # This minor eclipse application depends on the config parser library
    # and should be moved to the libres repository - it can not be
    # built here.
    if (BUILD_ERT)
       add_executable( ecl_quantile ecl/ecl_quantile.c )
       target_link_libraries( ecl_quantile config )
       list(APPEND apps ecl_quantile)
    endif()
    if (PTHREAD_LIBRARY)
        add_executable(bls block_fs/bls.c)
        target_link_libecl(bls)
        list(APPEND apps bls)
    endif ()

    add_executable(segment_info well/segment_info.c)
    add_executable(CF_dump      well/well_CF_dump.c)
    add_executable(ri_well_test well/ri_well_test.c)

    target_link_libecl(segment_info)
    target_link_libecl(CF_dump)
    target_link_libecl(ri_well_test)

    list(APPEND apps segment_info CF_dump ri_well_test)
endif()

if (BUILD_ECL_SUMMARY)
   add_executable(ecl_summary ecl/view_summary.cpp)
   target_link_libecl(ecl_summary)
   list(APPEND apps ecl_summary)

   install(FILES ${PROJECT_SOURCE_DIR}/man/man1/ecl_summary.1
           DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)
endif()

if(apps)
     if (ERT_USE_OPENMP)
        foreach(app ${apps})
            target_compile_options(${app} PUBLIC ${OpenMP_C_FLAGS})
            set_property(TARGET ${app} APPEND PROPERTY LINK_FLAGS ${OpenMP_C_FLAGS})
        endforeach()
    endif()
    install(TARGETS ${apps}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
