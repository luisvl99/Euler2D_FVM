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

option(EULER_BUILD_GUI "Build the Qt 6 GUI (needs Qt >= 6.5)" ON)

# ---------------------------------------------------------------------------
# Solver library: mesh, physics, test cases, exact Riemann solver (no Qt)
# ---------------------------------------------------------------------------
add_library(euler_core STATIC
    src/core/mesh_boundarytype.h
    src/core/mesh_cell.h
    src/core/mesh_face.h
    src/core/mesh_node.h
    src/core/mesh_structured.h          src/core/mesh_structured.cpp
    src/core/physics_boundarycondition.h
    src/core/physics_eulerphysics.h     src/core/physics_eulerphysics.cpp
    src/core/physics_eulerstate.h
    src/core/physics_flowparameters.h
    src/core/physics_primitivestate.h
    src/core/physics_solver.h           src/core/physics_solver.cpp
    src/core/physics_testcases.h        src/core/physics_testcases.cpp
    src/core/sod_exact.h                src/core/sod_exact.cpp
)
target_include_directories(euler_core PUBLIC src/core)

# ---------------------------------------------------------------------------
# Headless regression tests
# ---------------------------------------------------------------------------
enable_testing()
add_executable(solver_tests tests/tests.cpp)
target_link_libraries(solver_tests PRIVATE euler_core)
add_test(NAME solver_tests COMMAND solver_tests)

# ---------------------------------------------------------------------------
# Qt GUI
# ---------------------------------------------------------------------------
if(EULER_BUILD_GUI)
    find_package(Qt6 6.5 REQUIRED COMPONENTS Core Widgets)
    qt_standard_project_setup()

    qt_add_executable(Euler2D_FVM
        WIN32 MACOSX_BUNDLE
        src/gui/main.cpp
        src/gui/mainwindow.h            src/gui/mainwindow.cpp
        src/gui/mainwindow.ui
        src/gui/gui_colormap.h
        src/gui/gui_exportmanager.h     src/gui/gui_exportmanager.cpp
        src/gui/gui_fieldextractor.h
        src/gui/gui_fieldtype.h
        src/gui/gui_meshviewer.h        src/gui/gui_meshviewer.cpp
        src/gui/gui_residualplotter.h   src/gui/gui_residualplotter.cpp
        src/gui/gui_sodviewer.h         src/gui/gui_sodviewer.cpp
    )
    target_link_libraries(Euler2D_FVM PRIVATE euler_core Qt::Core Qt::Widgets)

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
endif()
