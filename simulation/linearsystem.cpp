#include "linearsystem.h"
#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl.h"
#include <stdexcept>

#include <algorithm>
#include <spdlog/spdlog.h>

long icy::LinearSystem::ClearAndResize(unsigned N_)
{
    auto t1 = std::chrono::high_resolution_clock::now();

    // clear the list of non-zero element indices
    this->N = N_;
    if(rhs.size() < N*DOFS)
    {
        rhs.resize(N*DOFS*2);
        sln.resize(N*DOFS*2);
    }

    if(csr_rows.size() < N+1) csr_rows.resize(N*2+1);

    std::fill(rhs.begin(), rhs.begin()+N*DOFS, 0);

    if(rows_neighbors.size()<N)
    {
        rows_neighbors.reserve(N*2);
        while(rows_neighbors.size()<N*2)
            rows_neighbors.push_back(std::make_unique<tbb::concurrent_vector<unsigned>>(20));
    }

    for(auto &row_vec : rows_neighbors) row_vec->clear();

    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
}


void icy::LinearSystem::AddNNZEntry(int row, int column)
{
    if(row < 0 || column < 0) return; // the element does not belong in the matrix
    if(row > column) std::swap(row,column);    // enforce upper-triangular matrix
    if((unsigned)column >= N) throw std::runtime_error("AddNNZEntry: trying to insert an element beyond the matrix size");
    rows_neighbors[row]->push_back(column);
}

void icy::LinearSystem::AddEntriesToStructure(const int* idx_begin, const int* idx_end)
{
    for(auto iter=(idx_begin+1); iter!=idx_end; ++iter)
        for(auto j=idx_begin; j!=iter; ++j)
            AddNNZEntry(*iter,*j);
}





long icy::LinearSystem::CreateStructure()
{
    auto t1 = std::chrono::high_resolution_clock::now();
    // CREATE STRUCTURE ARRAYS

    nnz = 0;

    // sort the neighbor list of each row
#pragma omp parallel for reduction(+:nnz)
    for(unsigned i=0;i<N;i++)
    {
        auto &rn = *rows_neighbors[i];
        rn.push_back(i);    // add diagonal entry
        std::sort(rn.begin(),rn.end());
        rn.resize(std::distance(rn.begin(),std::unique(rn.begin(), rn.end())));     // remove duplicates
        nnz+=(rn.size());
    }

    csr_rows[N] = nnz;
    if(csr_cols.size() < nnz) csr_cols.resize(nnz*2);
    nnz *= DOFS_SQ;

    if(vals.size() < nnz) vals.resize(nnz*2);
    std::fill(vals.begin(), vals.begin()+nnz, 0);

    // enumerate entries
    unsigned count=0;
    for(unsigned row=0;row<N;row++)
    {
        csr_rows[row] = count;

        auto &rn = *rows_neighbors[row];

        for(unsigned const &local_column : rn)
        {
            if(row > local_column) throw std::runtime_error("CreateStructure: matrix is not upper-triangular");

            csr_cols[count] = local_column;
            count++;
        }
    }

    if(nnz != DOFS_SQ*count) throw std::runtime_error("CreateStructure: nnz != count");

    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
}




unsigned icy::LinearSystem::get_offset(const int row, const int column) const
{
    int col_offset_begin = csr_rows[row];
    int col_offset_end = csr_rows[row+1];

    const int *start_pt = &csr_cols[col_offset_begin];
    const int *end_pt = &csr_cols[col_offset_end];

    auto it = std::lower_bound(start_pt,end_pt,column);
    if(it == end_pt || *it != column)
    {
        throw std::runtime_error("get_offset(): (i,j) index not found");
    }

    unsigned offset = std::distance(start_pt,it)+col_offset_begin;
    offset *= DOFS_SQ;

    return offset;
}


void icy::LinearSystem::AddToH(const int row, const int column, const double *v)
{
    if (row < 0 || column < 0 || row > column) return;
    else if((unsigned)row >= N || (unsigned)column >= N) throw std::runtime_error("AddToH: out of range");
    int offset = get_offset(row,column);

    for(unsigned i=0;i<DOFS_SQ;i++)
    {
#pragma omp atomic
        vals[offset + i] += v[i];
    }
}

void icy::LinearSystem::AddToC(const int idx, const double *v)
{
    if(idx < 0) return;
    if((unsigned)idx >= N) throw std::runtime_error("AddToC: index out of range");

    for(unsigned i=0;i<DOFS;i++)
    {
#pragma omp atomic
        rhs[idx*DOFS + i] += v[i];
    }
}


void icy::LinearSystem::AddToEquation(const double *lE, const double *qE, const std::initializer_list<int> ids)
{
    unsigned n = ids.size();

    Eigen::Map<const Eigen::MatrixXd> levec(lE, 1, n*DOFS);
    Eigen::Map<const Eigen::MatrixXd> qevec(qE, n*DOFS, n*DOFS);

    unsigned idx_i=0;
    for(auto iter_row=ids.begin(); iter_row!=ids.end(); ++iter_row,++idx_i)
    {
        int row = *iter_row;
        if(row < 0) continue;

        Eigen::Matrix<double, 1, DOFS> vec = -levec.block(idx_i*DOFS,0, DOFS,1);
        AddToC(row, vec.data());

        unsigned idx_j=0;
        for(auto iter_col=ids.begin(); iter_col!=ids.end(); ++iter_col,++idx_j)
        {
            int col = *iter_col;
            if(col < 0) continue;
            Eigen::Matrix<double,DOFS,DOFS> mat = qevec.block(idx_i*DOFS, idx_j*DOFS, DOFS, DOFS);
            mat.transposeInPlace();
            AddToH(row,col,mat.data());
        }
    }


}





long icy::LinearSystem::Solve(int verbosity)
{
    auto t1 = std::chrono::high_resolution_clock::now();

    int n = N;
    MKL_INT mtype = -2;       // Real symmetric matrix: -2;  real unsymmetric: 11
    MKL_INT nrhs = 1;     // Number of right hand sides.
    void *pt[64] = {};
    MKL_INT iparm[64] = {};
    MKL_INT maxfct, mnum, phase, error, msglvl;
    MKL_INT idum;
    iparm[0] = 1;       // No solver default
    iparm[1] = 3;       // Fill-in reordering from METIS (was 2)
    iparm[3] = 0;       // No iterative-direct algorithm
    iparm[4] = 0;       // No user fill-in reducing permutation
    iparm[5] = 0;       // Write solution into x
    iparm[6] = 0;       // Not in use
    iparm[7] = 0;       // Max numbers of iterative refinement steps
    iparm[8] = 0;
    iparm[9] = 8;       // Perturb the pivot elements with 1E-iparm[9];
    iparm[10] = 0;      // Use nonsymmetric permutation and scaling MPS
    iparm[11] = 0;
    iparm[12] = 0;      // Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy
    iparm[13] = 0;      // Output: Number of perturbed pivots
    iparm[14] = 0;
    iparm[15] = 0;
    iparm[16] = 0;
    iparm[17] = 1;      // 1 - disable report; Output: Number of nonzeros in the factor LU
    iparm[18] = 1;		// 1- disable report; output number of operations
    iparm[19] = 0;
#ifdef QT_DEBUG
    iparm[26] = 1;      // check matrix structure for errors
#else
    iparm[26] = 0;      // check matrix structure for errors
#endif
    iparm[27] = 0;      // 0 double; 1 single
    iparm[34] = 1;      // zero-base index
    iparm[36] = DOFS;    // BSR with block size DOFS
    maxfct = 1;
    mnum = 1;
    msglvl = verbosity; // use 1 for verbose output
    error = 0;
    phase = 13;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, vals.data(), csr_rows.data(), csr_cols.data(),
            &idum, &nrhs, iparm, &msglvl, rhs.data(), sln.data(), &error);

    if(error != 0) throw std::runtime_error("MKL solver error");

    phase = -1; //clean up
    double ddum;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, csr_rows.data(), csr_cols.data(),
            &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if(error != 0) throw std::runtime_error("MKL solver error");

    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
}

void icy::LinearSystem::AdjustCurrentGuess(int idx, Eigen::Matrix<double,DOFS,1> &vec)
{
    int b_idx = idx*DOFS;
    for(int i=0;i<DOFS;i++) vec(i)+=sln[b_idx+i];
}

double icy::LinearSystem::SqNormOfDx()
{
    double result = 0;
    int max_count = N*DOFS;
#pragma omp parallel for reduction(+:result)
    for (int i = 0; i < max_count; i++) result += sln[i]*sln[i];
    return result;
}

