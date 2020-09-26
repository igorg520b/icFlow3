#include "mainwindow.h"

#include <QApplication>
#include <QSurfaceFormat>
#include <QCommandLineParser>
#include <gmsh.h>
#include <omp.h>
#include <iostream>
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

//    icy::Element elem; elem.Test(); return 0;
//    icy::LinearSystem ls; ls.TestSolve(); ls.TestSolve(); return 0;
/*
    // h5cpp test
     double *buf;
     unsigned long N1 = 30;
     unsigned long N2 = 3;
     buf = new double[N1*N2];
     for(unsigned long i=0;i<N1;i++)
         for(unsigned long j=0;j<N2;j++) {
             buf[i*N2+j] = i*10+j;
         }

     h5::fd_t fd = h5::create("test.h5", H5F_ACC_TRUNC);
     h5::ds_t ds = h5::create<double>(fd, "test",
                                      h5::current_dims{10,N2},
                                      h5::max_dims{H5S_UNLIMITED,N2},
                                      h5::chunk{5,N2});

     h5::write<double>(ds,  buf, h5::count{10,N2});
     h5::set_extent(ds, h5::current_dims{20,N2});
*/
    gmsh::initialize();
//    double value;
//    gmsh::option::getNumber("Mesh.MaxNumThreads2D", value);
//    qDebug() << "Mesh.MaxNumThreads2D " << value;

    QApplication a(argc, argv);
    QApplication::setApplicationName("icFlow");
    QApplication::setApplicationVersion("3.0");

    // parse command line options
    QCommandLineParser parser;
    parser.addHelpOption();
    parser.addVersionOption();
    // An option with a value
    QCommandLineOption initialConfigOption(QStringList() << "f" << "file",
            QCoreApplication::translate("main", "Load initial configuration from a JSON file <fileName>"),
            QCoreApplication::translate("main", "fileName"));

    parser.addOption(initialConfigOption);
    parser.process(a);

    QSurfaceFormat fmt = QVTKOpenGLNativeWidget::defaultFormat();
//    fmt.setAlphaBufferSize(0);
    QSurfaceFormat::setDefaultFormat(fmt);

    MainWindow w;
    w.showMaximized();
    return a.exec();
}
