#if !defined(Q_MOC_RUN) // MOC has a glitch when parsing concurrent_set.h
#ifndef LINEARSYSTEM_H
#define LINEARSYSTEM_H
#include <concurrent_vector.h>
#include <Eigen/Core>
#include <map>
#include <unordered_map>
#include <vector>
#include <chrono>
#define DOFS 5

namespace icy { class LinearSystem; }

class icy::LinearSystem
{
public:
    double *vals=nullptr, *rhs=nullptr, *dx=nullptr;   // non-zero values, right-hand side and solution
    int *csr_rows=nullptr, *csr_cols=nullptr;  // structure arrays of the sparse matrix
    int N;      // number of variables
    int nnz;    // number of non-zero entries

    LinearSystem();

    // initializing and creating structure
    long ClearAndResize(std::size_t N);     // size N must be set; return execution time
    void AddElementToStructure(int row, int column);    // reserve non-zero positions one-by-one (thread-safe)
    long CreateStructure(); // return execution_time

    // creating the values array
    void SubtractRHS(const int idx, const Eigen::Matrix<double,DOFS,1> &vec); // we subtract, because RHS in N-R has opposite sign
    void AddLHS(const int row, const int column, const Eigen::Matrix<double,DOFS,DOFS> &mat);
    long Solve(int verbosity = 0);  // return execution time
    void AdjustCurrentGuess(int idx, Eigen::Matrix<double,DOFS,1> &vec);  // solution => convenient vector form

    double SqNormOfDx();

    // testing and benchmarking
    void Assert();
    void TestSolve(); // test solver with sample data (below)

private:
    // concurrent set allows combining the computation of forces with creation of structure
    std::vector<tbb::concurrent_vector<int>*> rows_Neighbors;
    int dx_length = 0;      // number of currently allocated elements for dx and rhs
    int vals_length = 0;    // currently allocated vals
//    std::vector<std::vector<std::pair<int,int>>*> rows_pcsr;   // per row mappings between columns and offset in "values"
    std::vector<std::vector<int>*> rows_pcsr;   // per row mappings between columns and offset in "values"
    int csr_rows_size = 0;
    int csr_cols_size = 0;
    int dvalsSize() { return nnz*DOFS*DOFS; }
    int dxSize() { return N*DOFS; }

    void ResizeRows();
};

#endif // LINEARSYSTEM_H
#endif
