#include <QApplication>
#include <omp.h>
#include <iostream>
//#include "linearsystem.h"

int main()
{
    std::cout << "num_threads " << omp_get_max_threads() << std::endl;
    std::cout << "testing threads" << std::endl;
    int nthreads, tid;
#pragma omp parallel
    { std::cout << omp_get_thread_num(); }
    std::cout << std::endl;

//    icy::Element elem; elem.Test(); return 0;
//    icy::LinearSystem ls; ls.TestSolve(); ls.TestSolve(); return 0;
    return 0;

}
