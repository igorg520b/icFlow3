#if !defined(Q_MOC_RUN) // MOC has a glitch when parsing concurrent_set.h
#ifndef LINEARSYSTEM_H
#define LINEARSYSTEM_H
#include <tbb/concurrent_vector.h>
#include <Eigen/Core>
#include <initializer_list>
#include <vector>
#include <chrono>

namespace icy { class LinearSystem; }

class icy::LinearSystem
{
public:
    constexpr static int DOFS = 5;
    constexpr static int DOFS_SQ = DOFS*DOFS;

    // initializing and creating structure
    long ClearAndResize(unsigned N_);     // size N must be set; return execution time
    long CreateStructure(); // return execution_time
    void AddEntriesToStructure(const int* idx_begin, const int* idx_end); // insert nxn matrix of indices of non-zero entries


    // add values to non-zero elements
    void AddToEquation(const double *linearEntries, const double *quadraticEntries, const std::initializer_list<int> ids);

    // creating the values array
//    void SubtractRHS(const int idx, const Eigen::Matrix<double,DOFS,1> &vec); // we subtract, because RHS in N-R has opposite sign
//    void AddLHS(const int row, const int column, const Eigen::Matrix<double,DOFS,DOFS> &mat);
    long Solve(int verbosity = 0);  // return execution time
    void AdjustCurrentGuess(int idx, Eigen::Matrix<double,DOFS,1> &vec);  // solution => convenient vector form

    double SqNormOfDx();

    // testing and benchmarking
    void Assert();
    void TestSolve(); // test solver with sample data (below)

private:
    std::vector<double> vals, rhs, sln;
    std::vector<int> csr_rows, csr_cols;
    unsigned N, nnz;

    // concurrent set allows combining the computation of forces with creation of structure
    std::vector<std::unique_ptr<tbb::concurrent_vector<unsigned>>> rows_neighbors;

    void AddNNZEntry(int row, int column);    // reserve non-zero positions one-by-one (thread-safe)
    void AddToH(const int row, const int column, const double *v);
    void AddToC(const int idx, const double *v);
    unsigned get_offset(const int row, const int column) const;


};

#endif // LINEARSYSTEM_H
#endif
