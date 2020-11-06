#include "linearsystem.h"
#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl.h"
#include <stdexcept>
#include <cmath>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <random>
#include <QDebug>

icy::LinearSystem::LinearSystem()
{
    // typical values for medium-sized FEM geometry
    int reserve = 1300000;
    rows_Neighbors.reserve(reserve);
    rows_pcsr.reserve(reserve);
}

long icy::LinearSystem::ClearAndResize(std::size_t N_)
{
    auto t1 = std::chrono::high_resolution_clock::now();

    // clear the list of non-zero element indices
    this->N = (int)N_;
    if((int)rows_Neighbors.size() < N)
    {
        while((int)rows_Neighbors.size()<N*1.2)
        rows_Neighbors.push_back(new tbb::concurrent_vector<int>(10));
    }

    if((int)rows_pcsr.size() < N)
    {
        while((int)rows_pcsr.size()<N*1.2)
            rows_pcsr.push_back(new std::vector<int>(10));
    }

#pragma omp parallel for
    for(int i=0;i<N;i++)
    {
        rows_Neighbors[i]->clear();
        rows_Neighbors[i]->push_back(i);    // diagonal elements must be non-zero
        rows_pcsr[i]->clear(); // clear the mapping of (i,j)->offset
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
}


void icy::LinearSystem::AddElementToStructure(int row, int column)
{
    if(row < 0 || column < 0) return; // the element does not belong in the matrix
    // since we are creatigng upper-triangular matrix, row<=column always
    if(row<=column) {
        if (row >= N) throw std::runtime_error("trying to insert element beyond matrix size");
        rows_Neighbors[row]->push_back(column);
    } else
    {
        if (column >= N) throw std::runtime_error("trying to insert element beyond matrix size");
        rows_Neighbors[column]->push_back(row);
    }
}

long icy::LinearSystem::CreateStructure()
{
    auto t1 = std::chrono::high_resolution_clock::now();
    // CREATE STRUCTURE ARRAYS

    // sort the neighbor list of each row
#pragma omp parallel for
    for(int i=0;i<N;i++)
    {
        tbb::concurrent_vector<int> &rn = *rows_Neighbors[i];
        std::sort(rn.begin(),rn.end());
        auto unique_res = std::unique(rn.begin(), rn.end());
        unsigned newSize = std::distance(rn.begin(),unique_res);
        rn.resize(newSize);
    }

    // count non-zero entries
    nnz = 0;
#pragma omp parallel for reduction(+:nnz)
    for(int i=0;i<N;i++) nnz+=rows_Neighbors[i]->size();

    // allocate structure arrays
    if(csr_rows_size < N+1)
    {
        csr_rows_size = N*1.3;
        delete csr_rows;
        csr_rows = new int[csr_rows_size];
        if(csr_rows == nullptr) throw std::runtime_error("csr_rows allocation error");
    }

    if(csr_cols_size < nnz)
    {
        csr_cols_size = nnz*1.5;
        delete csr_cols;
        csr_cols = new int[csr_cols_size];
        if(csr_cols == nullptr) throw std::runtime_error("csr_cols allocation error");
    }
    csr_rows[N] = nnz;

    // enumerate entries
    int count=0;
    for(int i=0;i<N;i++)
    {
        csr_rows[i] = count;
        if(rows_Neighbors[i]->size() == 0) throw std::runtime_error("matrix row contains no entries");
        tbb::concurrent_vector<int> &sorted_vec = *rows_Neighbors[i];

        int column_for_assertion = -1;
        for(int const &local_column : sorted_vec)
        {
            rows_pcsr[i]->push_back(count);
            csr_cols[count] = local_column;
            count++;

            if(local_column <= column_for_assertion) throw std::runtime_error("rn not sorted");
            column_for_assertion = local_column;
        }
    }

    if(nnz!=count) {
        qDebug()<<"csr_rows["<<N<<"]="<<nnz;
        qDebug()<<"nnz="<<nnz<<"; count="<<count;
        throw std::runtime_error("nnz != count");
    }

    // ALLOCATE AND CLEAR VALUE ARRAYS
    // allocate value arrays
    if (vals_length < dvalsSize())
    {
        delete vals;
        vals_length = dvalsSize()*1.5;
        vals = new double[vals_length];
    }

    if(dx_length < dxSize())
    {
        delete rhs;
        delete dx;
        dx_length = dxSize()*1.3;
        rhs = new double[dx_length];
        dx = new double[dx_length];
    }
    memset(rhs, 0, sizeof(double)*dxSize());
    memset(dx, 0, sizeof(double)*dxSize());
    memset(vals, 0, sizeof(double)*dvalsSize());

    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
}

void icy::LinearSystem::SubtractRHS(const int idx, const Eigen::Matrix<double,DOFS,1> &vec)
{
    if (idx < 0) return;
#ifdef QT_DEBUG
    if(idx >= N) throw std::runtime_error("ADDRHS: index out of range");
#endif
    int i3 = idx*DOFS;
    for(int i=0;i<DOFS;i++) {
#pragma omp atomic
        rhs[i3+i] -= vec.coeff(i);
    }
}

void icy::LinearSystem::AddLHS(const int row, const int column, const Eigen::Matrix<double,DOFS,DOFS> &mat)
{
    if (row < 0 || column < 0 || row > column) return;
#ifdef QT_DEBUG
    if(row >= N) throw std::runtime_error("LHS: row out of range");
    if(column >= N) throw std::runtime_error("LHS: column out of range");
#endif

    int offset=-1;
//    std::vector<std::pair<int,int>> &entry = *rows_pcsr[row];
//    for(std::pair<int,int> &p : entry) if(p.first == column) {offset=p.second; break;}
    tbb::concurrent_vector<int>*vec = rows_Neighbors[row];
    for(unsigned count = 0;count<vec->size();count++)
    {
        if(vec->at(count) == column) {offset = rows_pcsr[row]->at(count); break;}
    }
    if(offset<0) throw std::runtime_error("offset < 0");

#ifdef QT_DEBUG
    if(offset >= nnz) throw std::runtime_error("offset >= nnz");
#endif

    offset *= (DOFS*DOFS);
    for(int i=0;i<DOFS;i++)
        for(int j=0;j<DOFS;j++) {
#pragma omp atomic
            vals[offset + (j+i*DOFS)] += mat.coeff(i,j);
        }
}

void icy::LinearSystem::Assert()
{
    if(csr_rows[0] != 0) throw std::runtime_error("rows[0] != 0");
    if(csr_rows[N] != nnz) throw std::runtime_error("rows[N] != nnz");
    for (int i = 1; i < N + 1; i++)
        if (csr_rows[i] <= csr_rows[i - 1]) {
            std::cout << csr_rows[0] << "; " << csr_rows[1]<< "; " << csr_rows[2] << std::endl;
            throw std::runtime_error("rows[i] is not increasing");
        }

    // verify columns array, upper triangular
    for (int i = 0; i < N; i++)
    {
        if (csr_cols[csr_rows[i]] != i) throw std::runtime_error("structure not UT");
        for (int j = csr_rows[i]; j < csr_rows[i + 1] - 1; j++) {
            if (csr_cols[j + 1] < csr_cols[j])
            {
                qDebug() << "i="<<i<<"; j="<<j;
                qDebug()<<"csr_rows[i]="<<csr_rows[i]<<";csr_rows[i+1]="<<csr_rows[i+1];
                qDebug()<<"csr_cols[j]="<<csr_cols[j]<<";csr_cols[j+1]="<<csr_cols[j+1];
                throw std::runtime_error("cols in same row decreasing");
            }
            if (csr_cols[j + 1] == csr_cols[j]) throw std::runtime_error("cols in same row equal");
        }
    }

    for(int i=0;i<dxSize();i++)
        if(std::isnan(rhs[i])) throw std::runtime_error("rhs constains NaN");
    for(int i=0;i<dvalsSize();i++)
        if(std::isnan(vals[i])) throw std::runtime_error("matrix contains NaN");
}

long icy::LinearSystem::Solve(int verbosity)
{
    auto t1 = std::chrono::high_resolution_clock::now();
#ifdef QT_DEBUG
    Assert();
#endif
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
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, vals, csr_rows, csr_cols,
            &idum, &nrhs, iparm, &msglvl, rhs, dx, &error);

    if(error != 0) throw std::runtime_error("MKL solver error");

    phase = -1; //clean up
    double ddum;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, csr_rows, csr_cols,
            &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if(error != 0) throw std::runtime_error("MKL solver error");

    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
}

void icy::LinearSystem::AdjustCurrentGuess(int idx, Eigen::Matrix<double,DOFS,1> &vec)
{
    int b_idx = idx*DOFS;
    for(int i=0;i<DOFS;i++) vec(i)+=dx[b_idx+i];
}

double icy::LinearSystem::SqNormOfDx()
{
    double result = 0;
    int max_count = dxSize();
#pragma omp parallel for reduction(+:result)
    for (int i = 0; i < max_count; i++) result += dx[i]*dx[i];
    return result;
}

void icy::LinearSystem::TestSolve()
{
    std::vector<std::pair<int, int>> nnzs;

    int testing_size = 5000;
    int testing_nnz = 30000;
    nnzs.reserve(testing_size+testing_nnz);
    for(int i=0;i<testing_size;i++) nnzs.push_back(std::make_pair(i,i));
    for(int i=0;i<testing_nnz;i++) {
        int row = std::rand()%testing_size;
        int col = std::rand()%testing_size;
        if(row <= col) nnzs.push_back(std::make_pair(row, col));
        else nnzs.push_back(std::make_pair(col, row));
    }
    this->ClearAndResize((std::size_t)testing_size);

#pragma omp parallel for
    for(std::size_t i=0;i<nnzs.size();i++) {
        const std::pair<int,int> &t = nnzs[i];
        this->AddElementToStructure(t.first, t.second);
    }

    CreateStructure();

    double lower_bound = 0.1;
    double upper_bound = 10;
    std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
    std::default_random_engine re;

    for(std::size_t i=0;i<nnzs.size();i++) {
        const std::pair<int,int> &t = nnzs[i];
        Eigen::Matrix<double,DOFS,DOFS> mat = Eigen::Matrix<double,DOFS,DOFS>::Random();
        this->AddLHS(t.first, t.second, mat);
    }

    for(std::size_t i=0;(int)i<testing_size;i++)
    {
        Eigen::Matrix<double,DOFS,1> vec = Eigen::Matrix<double,DOFS,1>::Random();
        this->SubtractRHS(i, vec);
    }
    std::cout << "starting solve with N=" << N << "; nnz=" << nnz << std::endl;
    this->Assert();
    this->Solve(1);
    std::cout << "Sq dx = " << SqNormOfDx() << std::endl;
}
