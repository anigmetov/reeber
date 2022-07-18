# src files
foreach(real   double)

    # include files
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include
            ${CMAKE_CURRENT_SOURCE_DIR}/../amr-merge-tree/include
            ${AMREX_INSTALL_DIR}/include
#            /home/narn/code/LowFive/include
            ${DIY_INCLUDE_DIR})

    # AMR
    find_package(AMReX REQUIRED HINTS ${AMREX_INSTALL_DIR})

    set (mylibraries ${libraries}
            ${CMAKE_THREAD_LIBS_INIT}
            ${AMREX_INSTALL_DIR}/lib/libamrex.so
            HighFive
#            /home/narn/code/LowFive/build_debug/src/liblowfive.so
#            /home/narn/code/LowFive/build_debug/src/liblowfive-dist.a
            )

    add_executable(amr_connected_components_${real}
            ${CMAKE_CURRENT_SOURCE_DIR}/src/amr-connected-components.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/../amr-merge-tree/src/amr-plot-reader.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/../amr-merge-tree/src/hdf5-reader.cpp)

    if($real STREQUAL "float")
        target_compile_definitions(amr_connected_components_${real} PRIVATE REEBER_REAL=${real} AMREX_SPACEDIM=3 BL_USE_MPI=ON BL_USE_FLOAT=YES)
    else()
        target_compile_definitions(amr_connected_components_${real} PRIVATE REEBER_REAL=${real} AMREX_SPACEDIM=3 BL_USE_MPI=ON BL_USE_DOUBLE=YES)
    endif()

    target_link_libraries(amr_connected_components_${real} PUBLIC ${mylibraries})

    # for Henson

    add_library(henson_amr_connected_components_${real} SHARED
            ${CMAKE_CURRENT_SOURCE_DIR}/src/amr-connected-components.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/../amr-merge-tree/src/amr-plot-reader.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/../amr-merge-tree/src/hdf5-reader.cpp)

    if($real STREQUAL "float")
        target_compile_definitions(henson_amr_connected_components_${real} PRIVATE REEBER_REAL=${real} AMREX_SPACEDIM=3 BL_USE_MPI=ON BL_USE_FLOAT=YES)
    else()
        target_compile_definitions(henson_amr_connected_components_${real} PRIVATE REEBER_REAL=${real} AMREX_SPACEDIM=3 BL_USE_MPI=ON BL_USE_DOUBLE=YES)
    endif()


    set_target_properties(henson_amr_connected_components_${real} PROPERTIES PREFIX "")
    set_target_properties(henson_amr_connected_components_${real} PROPERTIES SUFFIX ".hx")

    target_link_libraries(henson_amr_connected_components_${real} ${HENSON_BUILD_DIR}/libhenson-pmpi.so)
    target_link_libraries(henson_amr_connected_components_${real} ${HENSON_BUILD_DIR}/libhenson.a)
    target_link_libraries(henson_amr_connected_components_${real} ${mylibraries})

    # Henson-specific
    set(linker_flags "-Wl,--export-dynamic")
    set(linker_flags "${linker_flags} -Wl,-u,henson_set_contexts,-u,henson_set_namemap")
    set_target_properties(henson_amr_connected_components_${real} PROPERTIES LINK_FLAGS ${linker_flags})

endforeach()
