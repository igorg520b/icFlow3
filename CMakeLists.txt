cmake_minimum_required(VERSION 3.5)

project(fl3 LANGUAGES CXX C)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

if(UNIX)
    find_package(spdlog REQUIRED)
    find_package(HDF5 REQUIRED COMPONENTS CXX)
endif()

find_package(Boost HINTS "C:\\Users\\s\\Projects\\libraries\\boost_1_79_0\\stage\\lib\\cmake\\Boost-1.79.0" REQUIRED)


find_package(VTK 9.0 HINTS "C:\\Users\\s\\Projects\\libraries\\vtk-build\\lib\\cmake\\vtk-9.1" REQUIRED COMPONENTS
  CommonColor
  CommonCore
  FiltersSources
  FiltersGeometry
  InteractionStyle
  RenderingContextOpenGL2
  RenderingCore
  RenderingFreeType
  RenderingGL2PSOpenGL2
  RenderingOpenGL2
  ViewsQt
  IOGeometry
  RenderingAnnotation
  RenderingLabel
  IOImage
)

find_package(Qt5 COMPONENTS Widgets REQUIRED)
find_package(Qt5Charts)

add_compile_definitions(QT_NO_KEYWORDS)


add_executable(fl3
    gui/main.cpp
    gui/windows/mainwindow.cpp
    gui/windows/mainwindow.h
    gui/windows/mainwindow.ui
    gui/icons.qrc
    gui/preferences_gui.h
    simulation/parameters_sim.h
    simulation/parameters_beam.h


    simulation/linearsystem.cpp
    simulation/linearsystem.h
    simulation/serializer.cpp
    simulation/serializer.h
    simulation/SimpleObjectPool.h

    geometry/baseelement.h

    geometry/node.cpp
    geometry/node.h
    geometry/edge.cpp
    geometry/edge.h

    geometry/boundaryedge.h

    geometry/element.cpp
    geometry/element.h
    geometry/mesh.cpp
    geometry/mesh.h
    simulation/modelcontrollerinterface.h
    simulation/model.cpp
    simulation/model.h
    geometry/floevisualization.h
    geometry/floevisualization.cpp

    simulation/frameinfo.h
    gui/backgroundworker.cpp
    gui/backgroundworker.h
)

target_link_libraries(fl3 PRIVATE Qt5::Widgets Qt5::Charts Qt5::Core)


target_include_directories(fl3 PRIVATE geometry )
target_include_directories(fl3 PRIVATE gui gui/windows)
target_include_directories(fl3 PRIVATE simulation simulation/material)


target_include_directories(fl3 PRIVATE ${Boost_INCLUDE_DIRS})


# OpenMP
if(WIN32)
    target_compile_options(fl3 PRIVATE /openmp)
else()
    target_compile_options(fl3 PRIVATE -fopenmp -Wno-stringop-overflow -Wno-invalid-offsetof)
    target_link_libraries(fl3 PRIVATE gomp)
endif()

# spdlog
if (WIN32)
    target_link_libraries(fl3 PRIVATE spdlog)
    target_include_directories(fl3 PRIVATE "C:\\Users\\s\\Projects\\libraries\\spdlog\\include")
    target_link_directories(fl3 PRIVATE "C:\\Users\\s\\Projects\\libraries\\spdlog-build\\Release")
else()
    target_link_libraries(fl3 PRIVATE spdlog::spdlog)
endif()

# VTK
target_link_libraries(fl3 PRIVATE ${VTK_LIBRARIES})

# TBB and MKL
if (WIN32)

    target_include_directories(fl3 PRIVATE "C:\\Program Files (x86)\\Intel\\oneAPI\\tbb\\2021.5.2\\include")
    target_link_directories(fl3 PRIVATE "C:\\Program Files (x86)\\Intel\\oneAPI\\tbb\\2021.5.2\\lib\\intel64\\vc_mt")
    target_link_libraries(fl3 PRIVATE tbb)

    target_include_directories(fl3 PRIVATE "C:\\Program Files (x86)\\Intel\\oneAPI\\mkl\\2022.0.3\\include")
    target_link_directories(fl3 PRIVATE "C:\\Program Files (x86)\\Intel\\oneAPI\\mkl\\2022.0.3\\lib\\intel64")
    target_link_libraries(fl3 PRIVATE mkl_intel_lp64_dll mkl_tbb_thread_dll mkl_core_dll)

else()
    target_include_directories(fl3 PRIVATE /usr/include/mkl)
    target_link_libraries(fl3 PRIVATE tbb mkl_intel_lp64 mkl_gnu_thread mkl_core gomp m dl)
    target_link_libraries(fl3 PRIVATE ${TBB_IMPORTED_TARGETS})
endif ()


# Gmsh
target_include_directories(fl3 PRIVATE "C:\\Users\\s\\Projects\\libraries\\gmsh-4.9.5-source\\api")
target_link_directories(fl3 PRIVATE "C:\\Users\\s\\Projects\\libraries\\gmsh-build\\Release")
target_link_libraries(fl3 PRIVATE gmsh)

# Eigen
find_package(Eigen3 REQUIRED NO_MODULE HINTS "C:\\Users\\s\\Projects\\libraries\\eigen-build")
target_link_libraries(fl3 PRIVATE Eigen3::Eigen)

# HDF5
if(WIN32)
    add_compile_definitions(H5_BUILT_AS_DYNAMIC_LIB)
    target_include_directories(fl3 PRIVATE "C:\\Users\\s\\Projects\\libraries\\hdf5-build\\src")
    target_include_directories(fl3 PRIVATE "C:\\Users\\s\\Projects\\libraries\\hdf5-develop\\src")
    target_include_directories(fl3 PRIVATE "C:\\Users\\s\\Projects\\libraries\\hdf5-develop\\c++\\src")
    target_link_directories(fl3 PRIVATE "C:\\Users\\s\\Projects\\libraries\\hdf5-build\\bin\\Release")
    target_link_libraries(fl3 PRIVATE hdf5_cpp hdf5 )
else()
    target_link_libraries(fl3 PRIVATE ${HDF5_LIBRARIES})
    target_include_directories(fl3 PRIVATE ${HDF5_INCLUDE_DIRS})
endif()

#PBrowser
if(WIN32)
    target_link_directories(fl3 PRIVATE "C:\\Users\\s\\Projects\\libraries\\pbrowser-build")
    target_include_directories(fl3 PRIVATE "C:\\Users\\s\\Projects\\libraries\\PBrowser")
    target_include_directories(fl3 PRIVATE "C:\\Users\\s\\Projects\\libraries\\PBrowser\\propertybrowser")
    target_include_directories(fl3 PRIVATE "C:\\Users\\s\\Projects\\libraries\\PBrowser\\objectpropertybrowser")
else()
    target_link_directories(fl3 PRIVATE "")
    target_include_directories(fl3 PRIVATE /home/s2/Projects-Qt/PBrowser/objectpropertybrowser)
    target_include_directories(fl3 PRIVATE /home/s2/Projects-Qt/PBrowser/propertybrowser)
    target_include_directories(fl3 PRIVATE /home/s2/Projects-Qt/PBrowser)
endif()
target_link_libraries(fl3 PRIVATE PBrowser)
