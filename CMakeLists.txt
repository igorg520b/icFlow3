cmake_minimum_required(VERSION 3.5)

project(fl3 LANGUAGES CXX C)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)
find_package(Qt5Charts)

find_package(HDF5 REQUIRED COMPONENTS C)
find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(spdlog REQUIRED)

find_package(VTK 9.0 REQUIRED COMPONENTS
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


add_executable(fl3
    #    gui/main_test.cpp
    gui/main.cpp
    gui/windows/mainwindow.cpp
    gui/windows/mainwindow.h
    gui/windows/mainwindow.ui
    gui/icons.qrc
    gui/preferences_gui.h
    simulation/parameters_sim.h

    #property browser
    gui/propertybrowser/qtvariantproperty.h
    gui/propertybrowser/qtbuttonpropertybrowser.cpp
    gui/propertybrowser/qtgroupboxpropertybrowser.cpp
    gui/propertybrowser/qtgroupboxpropertybrowser.h
    gui/propertybrowser/qteditorfactory.cpp
    gui/propertybrowser/qteditorfactory.h
    gui/propertybrowser/qtpropertybrowserutils.cpp
    gui/propertybrowser/qtvariantproperty.cpp
    gui/propertybrowser/qtpropertymanager.cpp
    gui/propertybrowser/qtpropertymanager.h
    gui/propertybrowser/qtpropertybrowserutils_p.h
    gui/propertybrowser/qtpropertybrowser.qrc
    gui/propertybrowser/qttreepropertybrowser.h
    gui/propertybrowser/qttreepropertybrowser.cpp
    gui/propertybrowser/qtpropertybrowser.h
    gui/propertybrowser/qtbuttonpropertybrowser.h
    gui/propertybrowser/qtpropertybrowser.cpp
    gui/objectpropertybrowser/objectpropertybrowser.cpp
    gui/objectpropertybrowser/objectpropertybrowser.h

    simulation/linearsystem.cpp
    simulation/linearsystem.h
    simulation/serializer.cpp
    simulation/serializer.h
    simulation/SimpleObjectPool.h

    geometry/solid3d.cpp
    geometry/solid3d.h

    geometry/baseelement.h

    geometry/node.cpp
    geometry/node.h
    geometry/edge.cpp
    geometry/edge.h

    geometry/boundaryedge.h

    geometry/element.cpp
    geometry/element.h
    geometry/geometry.cpp
    geometry/geometry_stable.cpp
    geometry/geometry.h
    geometry/model.cpp
    geometry/model.h
    geometry/floevisualization.h
    geometry/floevisualization.cpp

    simulation/frameinfo.h
    simulation/modelcontroller.cpp
    simulation/modelcontroller.h
    gui/backgroundworker.cpp
    gui/backgroundworker.h
)


target_include_directories(fl3 PRIVATE geometry gui/propertybrowser)
target_include_directories(fl3 PRIVATE gui gui/objectpropertybrowser gui/windows)
target_include_directories(fl3 PRIVATE simulation simulation/material)

# MKL path is set manually
target_include_directories(fl3 PRIVATE /home/s2/intel/mkl/include)
target_link_directories(fl3 PRIVATE /home/s2/intel/mkl/lib/intel64_lin)

target_include_directories(fl3 PRIVATE ${Boost_INCLUDE_DIRS})
target_include_directories(fl3 PRIVATE ${HDF5_INCLUDE_DIRS})

target_link_libraries(fl3 PRIVATE Qt5::Widgets Qt5::Charts Qt5::Core)
target_link_libraries(fl3 PRIVATE ${VTK_LIBRARIES})
target_link_libraries(fl3 PRIVATE gmsh)
target_link_libraries(fl3 PRIVATE tbb)
target_link_libraries(fl3 PRIVATE mkl_intel_lp64 mkl_gnu_thread mkl_core gomp m dl)
#target_link_libraries(fl3 PRIVATE pthread)
target_link_libraries(fl3 PRIVATE ${HDF5_C_LIBRARIES})
target_link_libraries(fl3 PRIVATE z boost_serialization)
target_link_libraries(fl3 PRIVATE ${Boost_LIBRARIES})
target_link_libraries(fl3 PRIVATE Eigen3::Eigen)
target_link_libraries(fl3 PRIVATE spdlog::spdlog)

target_compile_options(fl3 PRIVATE -fopenmp)
