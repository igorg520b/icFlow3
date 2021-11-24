#include "mainwindow.h"

#include <QApplication>
#include <QSurfaceFormat>
#include <QCommandLineParser>
#include <gmsh.h>
#include <omp.h>
#include <iostream>
#include "spdlog/spdlog.h"

#include "linearsystem.h"
#include "element.h"

int main(int argc, char *argv[])
{
    std::cout << "num_threads " << omp_get_max_threads() << std::endl;
    std::cout << "testing threads" << std::endl;
    int nthreads, tid;
#pragma omp parallel
    { std::cout << omp_get_thread_num(); }
    std::cout << std::endl;

    spdlog::set_pattern("%v");
    gmsh::initialize();

    QApplication a(argc, argv);
    QApplication::setApplicationName("icFlow");
    QApplication::setApplicationVersion("3.0");


    QSurfaceFormat fmt = QVTKOpenGLNativeWidget::defaultFormat();
    QSurfaceFormat::setDefaultFormat(fmt);

    MainWindow w;
    w.showMaximized();
    return a.exec();
}
