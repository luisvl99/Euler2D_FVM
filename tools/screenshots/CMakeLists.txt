# Standalone helper that regenerates the GUI screenshots in docs/screenshots.
# It is not part of the main build.  From the repository root:
#
#   cmake -S tools/screenshots -B build-screenshots -DCMAKE_BUILD_TYPE=Release
#         -DCMAKE_PREFIX_PATH=/path/to/Qt/6.x/<kit>
#   cmake --build build-screenshots
#   build-screenshots/screenshots docs/screenshots
#
# To render without opening windows, set QT_QPA_PLATFORM=offscreen (and on
# Windows QT_QPA_FONTDIR=C:/Windows/Fonts so text renders).

cmake_minimum_required(VERSION 3.19)
project(Euler2D_FVM_screenshots LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Qt6 6.5 REQUIRED COMPONENTS Core Widgets)
qt_standard_project_setup()

set(SRC "${CMAKE_CURRENT_SOURCE_DIR}/../../src")
file(GLOB CORE_SOURCES "${SRC}/core/*.cpp")

qt_add_executable(screenshots
    screenshots.cpp
    ${CORE_SOURCES}
    ${SRC}/gui/mainwindow.h            ${SRC}/gui/mainwindow.cpp
    ${SRC}/gui/mainwindow.ui
    ${SRC}/gui/gui_exportmanager.h     ${SRC}/gui/gui_exportmanager.cpp
    ${SRC}/gui/gui_meshviewer.h        ${SRC}/gui/gui_meshviewer.cpp
    ${SRC}/gui/gui_residualplotter.h   ${SRC}/gui/gui_residualplotter.cpp
    ${SRC}/gui/gui_sodviewer.h         ${SRC}/gui/gui_sodviewer.cpp
)
target_include_directories(screenshots PRIVATE ${SRC}/core ${SRC}/gui)
target_link_libraries(screenshots PRIVATE Qt::Core Qt::Widgets)
