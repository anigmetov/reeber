# src files
foreach(real   double)

    # include files
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include
            ${CMAKE_CURRENT_SOURCE_DIR}/../amr-merge-tree/include
            ${LOW_FIVE_DIR}/include
            ${LOW_FIVE_DIR}/ext/fmt/include
            ${DIY_INCLUDE_DIR})

    set (mylibraries ${libraries}
            ${CMAKE_THREAD_LIBS_INIT}
            HighFive
            ${LOW_FIVE_BUILD_DIR}/src/liblowfive.so
            ${LOW_FIVE_BUILD_DIR}/src/liblowfive-dist.a
            )

    add_library(amr_connected_components_${real}_for_lowfive SHARED
            ${CMAKE_CURRENT_SOURCE_DIR}/src/amr-connected-components-lib.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/../amr-merge-tree/src/hdf5-reader.cpp)

    target_include_directories(amr_connected_components_${real}_for_lowfive PUBLIC
            ${LOW_FIVE_DIR}/ext/fmt/include
            ${LOW_FIVE_DIR}/ext/spdlog/include
            ${CMAKE_CURRENT_SOURCE_DIR}/include
            ${CMAKE_CURRENT_SOURCE_DIR}/../amr-merge-tree/include
            ${LOW_FIVE_DIR}/include
            ${DIY_INCLUDE_DIR})

    target_compile_definitions(amr_connected_components_${real}_for_lowfive PRIVATE REEBER_REAL=${real} AMREX_SPACEDIM=3 BL_USE_MPI=ON ENABLE_DP=OFF BL_USE_FLOAT=0 SPDLOG_FMT_EXTERNAL=ON)

    # target_link_libraries(amr_merge_tree PUBLIC ${libraries})
    target_link_libraries(amr_connected_components_${real}_for_lowfive PUBLIC ${mylibraries})

    # AMR
    find_package(AMReX REQUIRED HINTS ${AMREX_INSTALL_DIR})

    set (mylibraries ${libraries}
            ${CMAKE_THREAD_LIBS_INIT}
            AMReX::amrex
            HighFive
            )



add_executable(amr_connected_components_${real}
        ${CMAKE_CURRENT_SOURCE_DIR}/src/amr-connected-components.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/../amr-merge-tree/src/amr-plot-reader.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/../amr-merge-tree/src/hdf5-reader.cpp)

target_compile_definitions(amr_connected_components_${real} PRIVATE REEBER_REAL=${real} AMREX_SPACEDIM=3 BL_USE_MPI=ON ENABLE_DP=OFF BL_USE_FLOAT=0)

# target_link_libraries(amr_merge_tree PUBLIC ${libraries})
target_link_libraries(amr_connected_components_${real} PUBLIC ${mylibraries})

endforeach()
