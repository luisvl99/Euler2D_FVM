cmake_minimum_required(VERSION 3.19)
project(Euler2D_FVM LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra)
endif()

find_package(Qt6 6.5 REQUIRED COMPONENTS Core Widgets)

qt_standard_project_setup()

qt_add_executable(Euler2D_FVM
    WIN32 MACOSX_BUNDLE
    main.cpp
    mainwindow.cpp
    mainwindow.h
    mainwindow.ui
    mesh_structured.h mesh_structured.cpp
    physics_primitivestate.h
    physics_eulerphysics.h physics_eulerphysics.cpp
    gui_meshviewer.h gui_meshviewer.cpp
    mesh_node.h
    mesh_face.h
    mesh_cell.h
    mesh_boundarytype.h
    gui_fieldtype.h
    gui_colormap.h
    gui_fieldextractor.h
    physics_boundarycondition.h
    physics_flowparameters.h
    physics_solver.h physics_solver.cpp
    physics_eulerstate.h
    gui_exportmanager.h gui_exportmanager.cpp
    gui_residualplotter.h gui_residualplotter.cpp
    sod_exact.h sod_exact.cpp
    physics_testcases.h physics_testcases.cpp
    gui_sodviewer.h gui_sodviewer.cpp
)

target_link_libraries(Euler2D_FVM
    PRIVATE
        Qt::Core
        Qt::Widgets
)

# Headless regression tests (no Qt): exact Riemann solver + test cases
enable_testing()
add_executable(solver_tests tests.cpp
    mesh_structured.cpp physics_eulerphysics.cpp physics_solver.cpp sod_exact.cpp
    physics_testcases.cpp)
add_test(NAME solver_tests COMMAND solver_tests)

include(GNUInstallDirs)

install(TARGETS Euler2D_FVM
    BUNDLE  DESTINATION .
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

qt_generate_deploy_app_script(
    TARGET Euler2D_FVM
    OUTPUT_SCRIPT deploy_script
    NO_UNSUPPORTED_PLATFORM_ERROR
)
install(SCRIPT ${deploy_script})
