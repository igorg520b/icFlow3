#include "element.h"
#include "edge.h"
#include "model.h"
#include "node.h"

#include <cstdlib>
#include <algorithm>

// barycentric coords of Gauss points
double icy::Element::N[3][3] = {
       {2.0/3.0, 1.0/6.0, 1.0/6.0},
       {1.0/6.0, 2.0/3.0, 1.0/6.0},
       {1.0/6.0, 1.0/6.0, 2.0/3.0}};

void icy::Element::ComputeInitialNormal()
{
    // translate the element
    Eigen::Vector3d p1, p2;
    p1 = nds[1]->x_initial - nds[0]->x_initial;
    p2 = nds[2]->x_initial - nds[0]->x_initial;

    normal_initial = p1.cross(p2);
    area_initial=normal_initial.norm()/2;
    if(area_initial<1e-8) {
        qDebug() << "element " << nds[0]->locId << ", " << nds[1]->locId << ", " << nds[2]->locId;
        qDebug() << "area" << area_initial;
        qDebug() << "nd0 " << nds[0]->x_initial.x() << ", " << nds[0]->x_initial.y();
        qDebug() << "nd1 " << nds[1]->x_initial.x() << ", " << nds[1]->x_initial.y();
        qDebug() << "nd2 " << nds[2]->x_initial.x() << ", " << nds[2]->x_initial.y();
        throw std::runtime_error("degenerate element created");
    }
    normal_initial.normalize();
}


void icy::Element::UpdateSparseSystemEntries(LinearSystem &ls)
{
    if(nds[0]->lsId < 0 && nds[1]->lsId < 0 && nds[2]->lsId < 0) return;

    // reserve non-zero entries in the sparse matrix structure
    ls.AddElementToStructure(nds[0]->lsId, nds[1]->lsId);
    ls.AddElementToStructure(nds[0]->lsId, nds[2]->lsId);
    ls.AddElementToStructure(nds[1]->lsId, nds[2]->lsId);
}



void icy::Element::DistributeStresses()
{
    // distribute the values from elements to nodes
    for(int i=0;i<3;i++) {
        Node *nd = nds[i];
        double coeff1 = 1.0/nd->adjacent_elems.size();//area_initial/nd->area;
        double coeff2 = coeff1/3;//area_initial/(3*nd->area);

        for(int j=0;j<3;j++)
        {
#pragma omp atomic
            nd->str_b[j] += a_str_b[j]*coeff1;
#pragma omp atomic
            nd->str_m[j] += a_str_m[j]*coeff1;
#pragma omp atomic
            nd->str_b_top[j] += a_str_b_top[j]*coeff1;
#pragma omp atomic
            nd->str_b_bottom[j] -= a_str_b_top[j]*coeff1;
        }

        float str_s_combined[2] = {};
        for(int j=0;j<3;j++) {
            str_s_combined[0] += a_str_s[j][0]*coeff2*N[i][j];
            str_s_combined[1] += a_str_s[j][1]*coeff2*N[i][j];
        }
#pragma omp atomic
            nd->str_s[0] += str_s_combined[0];
#pragma omp atomic
            nd->str_s[1] += str_s_combined[1];

//        nd->normal_n += normal_n;
    }

    if(principal_stress_exceeds_threshold)
        for(int k=0;k<3;k++) nds[k]->potentially_can_fracture=true;
}
/*
void icy::Element::rotationMatrix(Eigen::Vector3d &p1, Eigen::Vector3d &p2, Eigen::Matrix3d &result,
                                   double &area, Eigen::Vector3d &normal)
{
    Eigen::Vector3d r1 = p1.normalized();
    normal = p1.cross(p2); area=normal.norm()/2; normal.normalize();
    Eigen::Vector3d r2 = normal.cross(r1).normalized();
    result << r1, r2, normal;
}

void icy::Element::rotationMatrix_alt(Eigen::Vector3d &p12, Eigen::Vector3d &p13,
                        Eigen::Matrix3d &result,
                        double &area, Eigen::Vector3d &normal)
{
    normal = p12.cross(p13); area=normal.norm()/2; normal.normalize();
    double nx = normal.coeff(0);
    double ny = normal.coeff(0);
    double nz = normal.coeff(2);

    Eigen::Vector3d r1;
    if(nx == 0 && ny == 0)
    {
        r1 << 1, 0, 0;
    }
    else if(nx != 0 && nz != 0)
    {
        double c1 = nx/nz;
        double c1sq = c1*c1;
        r1 << 1.0/sqrt(1.0+c1sq), 0, -1.0/sqrt(1.0+1.0/c1sq);
        r1.normalize();
    }
    else {
        r1 = p12.normalized();
    }

    Eigen::Vector3d r2 = normal.cross(r1).normalized();
    result << r1, r2, normal;
}
*/

void icy::Element::ComputeNormal()
{
    Eigen::Vector3d u = (nds[1]->xn - nds[0]->xn).block(0,0,3,1);
    Eigen::Vector3d v = (nds[2]->xn - nds[0]->xn).block(0,0,3,1);
    normal_n = u.cross(v);
}



//=============== fracture helpers
icy::Node* icy::Element::getOppositeNode(icy::Edge edge)
{
    icy::Node *nd0 = edge.nds[0];
    icy::Node *nd1 = edge.nds[1];

    for(int i=0;i<3;i++)
    {
        int idx_next = (i+1)%3;
        if((nds[i] == nd0 && nds[idx_next] == nd1)||
                (nds[i] == nd1 && nds[idx_next] == nd0))
            return nds[(i+2)%3];
    }
    throw std::runtime_error("opposite node not found");
}

icy::Node* icy::Element::getOppositeNode(Node *nd0, Node* nd1)
{
    for(int i=0;i<3;i++)
    {
        int idx_next = (i+1)%3;
        if((nds[i] == nd0 && nds[idx_next] == nd1)||
                (nds[i] == nd1 && nds[idx_next] == nd0))
            return nds[(i+2)%3];
    }
    throw std::runtime_error("opposite node not found 2");
}


void icy::Element::ReplaceNode(icy::Node *replaceWhat, icy::Node *replaceWith)
{
    if(nds[0] == replaceWhat) nds[0] = replaceWith;
    else if(nds[1] == replaceWhat) nds[1] = replaceWith;
    else if(nds[2] == replaceWhat) nds[2] = replaceWith;
    else throw std::runtime_error("replaceWhat is not in nds[]");
    ComputeInitialNormal();
}

Eigen::Vector3d icy::Element::getCenter()
{
    return (nds[0]->x_initial + nds[1]->x_initial + nds[2]->x_initial)/3.0;
}


void icy::Element::getIdxs(const icy::Node*nd, short &thisIdx, short &CWIdx, short &CCWIdx) const
{
    if(nd==nds[0]) thisIdx=0;
    else if(nd==nds[1]) thisIdx=1;
    else if(nd==nds[2]) thisIdx=2;
    else {
        std::cout << "getIdxs; node does not belong to the element" << std::endl;
        throw std::runtime_error("getIdxs; node does not belong to the element");
    }

    CWIdx = (thisIdx+1)%3;
    CCWIdx = (thisIdx+2)%3;
}

icy::Edge icy::Element::CWEdge(const Node* nd) const
{
    short thisIdx, CWIdx, CCWIdx;
    getIdxs(nd, thisIdx, CWIdx, CCWIdx);
    return edges[CCWIdx];
}

icy::Edge icy::Element::CCWEdge(const Node* nd) const
{
    short thisIdx, CWIdx, CCWIdx;
    getIdxs(nd, thisIdx, CWIdx, CCWIdx);
    return edges[CWIdx];

}

icy::Edge icy::Element::OppositeEdge(const Node* nd) const
{
    short thisIdx, CWIdx, CCWIdx;
    getIdxs(nd, thisIdx, CWIdx, CCWIdx);
    return edges[thisIdx];
}


icy::Edge icy::Element::getEdgeOppositeToNode(icy::Node *nd)
{
    int thisIdx;
    if(nd==nds[0]) thisIdx=0;
    else if(nd==nds[1]) thisIdx=1;
    else if(nd==nds[2]) thisIdx=2;
    else throw std::runtime_error("getIdxs; node does not belong to the element");
    return edges[thisIdx];
}

short icy::Element::getNodeIdx(Node *nd)
{
    if(nds[0]==nd) return 0;
    else if(nds[1]==nd) return 1;
    else if(nds[2]==nd) return 2;
    else
    {
        nd->PrintoutFan();
        std::cout << "trying to find the index of the node " << nd->locId << "\n";
        std::cout << "in the element: " << nds[0]->locId << ", " << nds[1]->locId << ", " << nds[2]->locId << std::endl;
        throw std::runtime_error("getNodeIdx");
    }
}

icy::Element* icy::Element::getAdjacentElementOppositeToNode(Node *nd)
{
    if(!ContainsNode(nd)) throw std::runtime_error("getAdjacentElementOppositeToNode");
    return adj_elems[getNodeIdx(nd)];
}

void icy::Element::AssertEdges()
{
    // verify that edes are correctly initialized

    for(int i=0;i<3;i++)
    {
        if(!ContainsNode(edges[i].nds[0])) throw std::runtime_error("AssertEdges");
        if(!ContainsNode(edges[i].nds[1])) throw std::runtime_error("AssertEdges");
    }
}


//=================

void icy::Element::ComputeMatrices(SimParams &prms,
                     Eigen::Matrix3d &elasticityMatrix,
                     Eigen::Matrix2d &D_mats,
                     Eigen::Matrix<double,3,DOFS*3> &bmat_b,
                     Eigen::Matrix<double,2,DOFS*3> (&bmat_s)[3],
                     Eigen::Matrix<double,3,DOFS*3> &bmat_m,
                     Eigen::Matrix<double,DOFS*3,DOFS*3> &K)
{
    // translate the element
    Eigen::Vector3d p1, p2;
    p1 = nds[1]->x_initial - nds[0]->x_initial;
    p2 = nds[2]->x_initial - nds[0]->x_initial;

    normal_initial = p1.cross(p2);
    area_initial=normal_initial.norm()/2;
    if(area_initial<1e-8) {
        qDebug() << "element " << nds[0]->locId << ", " << nds[1]->locId << ", " << nds[2]->locId;
        qDebug() << "area" << area_initial;
        qDebug() << "nd0 " << nds[0]->x_initial.x() << ", " << nds[0]->x_initial.y();
        qDebug() << "nd1 " << nds[1]->x_initial.x() << ", " << nds[1]->x_initial.y();
        qDebug() << "nd2 " << nds[2]->x_initial.x() << ", " << nds[2]->x_initial.y();
        throw std::runtime_error("degenerate element created");
    }
    normal_initial.normalize();

    // (use 1-based indices, yij = yi-yj)
    double x1 = 0;
    double x2 = p1.x();
    double x3 = p2.x();
    double y1 = 0;
    double y2 = p1.y();
    double y3 = p2.y();

    double A2 = 2*area_initial;
    double y23, y31, y12, x32, x13, x21;
    y23 = y2-y3;
    y31 = y3-y1;
    y12 = y1-y2;
    x32 = x3-x2;
    x13 = x1-x3;
    x21 = x2-x1;

    // derivatives of shape functions
    double dN1dx = y23/A2;
    double dN2dx = y31/A2;
    double dN3dx = y12/A2;
    double dN1dy = x32/A2;
    double dN2dy = x13/A2;
    double dN3dy = x21/A2;

    // bending strain-displacement matrix
    bmat_b <<
         0,0,0,dN1dx,0,       0,0,0,dN2dx,0,        0,0,0,dN3dx,0,
         0,0,0,0,dN1dy,       0,0,0,0,dN2dy,        0,0,0,0,dN3dy,
         0,0,0,dN1dy,dN1dx,  0,0,0,dN2dy,dN2dx,   0,0,0,dN3dy,dN3dx;

    // membrane strain-displacement matrix
    bmat_m <<
         dN1dx,0,0,0,0,       dN2dx,0,0,0,0,      dN3dx,0,0,0,0,
         0,dN1dy,0,0,0,       0,dN2dy,0,0,0,      0,dN3dy,0,0,0,
         dN1dy,dN1dx,0,0,0,   dN2dy,dN2dx,0,0,0,  dN3dy,dN3dx,0,0,0;

    // shear
    for(int i=0;i<3;i++)
        bmat_s[i] <<
           0,0,dN1dx,-N[i][0],0,   0,0,dN2dx,-N[i][1],0,     0,0,dN3dx,-N[i][2], 0,
           0,0,dN1dy,0,-N[i][0],   0,0,dN2dy,0,-N[i][1],     0,0,dN3dy,0,-N[i][2];

    double thickness = prms.Thickness;
    // K and M depend on rho, Young's modulus and Poisson's ratio,
    // therefore they are computed after these parameters are set

//    Eigen::Matrix<double,DOFS*3,DOFS*3> K_b, K_m, K_s;
//    K = Eigen::Matrix<double,DOFS*3,DOFS*3>::Zero();
    K = bmat_m.transpose()*elasticityMatrix*bmat_m*(area_initial*thickness);

    double coeff = area_initial*thickness*thickness*thickness/12.0;
    K += bmat_b.transpose()*elasticityMatrix*bmat_b*coeff;

    for(int i=0;i<3;i++) K += bmat_s[i].transpose()*D_mats*bmat_s[i]*(thickness*area_initial/3.0);
}

void icy::Element::ComputeElasticForce(icy::LinearSystem &ls, icy::SimParams &prms, double,
                                       Eigen::Matrix3d &elasticityMatrix,
                                       Eigen::Matrix2d &D_mats)
{
    if(nds[0]->lsId < 0 && nds[1]->lsId < 0 && nds[2]->lsId < 0) return;

    Eigen::Matrix<double,3,DOFS*3> bmat_b;
    Eigen::Matrix<double,2,DOFS*3> bmat_s[3];   // 3 gauss points
    Eigen::Matrix<double,3,DOFS*3> bmat_m;
    Eigen::Matrix<double,DOFS*3,DOFS*3> K;    // element stiffness matrix (3 gauss points)
    ComputeMatrices(prms, elasticityMatrix, D_mats, bmat_b, bmat_s, bmat_m, K);

    Eigen::Matrix<double,DOFS*3,1> un;
    Eigen::Matrix<double,DOFS*3,1> ut;
    un << nds[0]->un, nds[1]->un, nds[2]->un;
    ut << nds[0]->ut, nds[1]->ut, nds[2]->ut;

    // calculate elastic forces and Hessian at step n+1

    // absolute position is not important, the origin is set at node 0
    Eigen::Matrix<double,DOFS*3,1> Fn, Fnp1;     // internal force at steps n and n+1
    Eigen::Matrix<double,DOFS*3,DOFS*3> &dFnp1=K;

//    FdF(un, Fn, nullptr);
//    FdF(ut, Fnp1, &dFnp1);
    Fn = K*un;
    Fnp1 = K*ut;

    // combine the results into a linearized equation of motion with HHT-alpha integration scheme
    double alpha = prms.HHTalpha;
    Eigen::Matrix<double,DOFS*3,1> F;      // right-hand side of the equation is equal to -F
    Eigen::Matrix<double,DOFS*3,DOFS*3> dF;

    F = Fn*alpha + Fnp1*(1-alpha);
    dF= dFnp1*(1-alpha);

#ifdef QT_DEBUG
    // assert
    for(int i=0;i<DOFS*3;i++)
        for(int j=0;j<DOFS*3;j++)
            if(std::isnan(dF(i,j)))
                throw std::runtime_error("elem.ComputeElasticForce: dF contains NaN");
#endif

    // assemble
    for(int i=0;i<3;i++) {
        int row = nds[i]->lsId;
        Eigen::Matrix<double,DOFS,1> locF = F.block(i*DOFS,0,DOFS,1);
        ls.SubtractRHS(row, locF);
        for(int j=0;j<3;j++) {
            int col = nds[j]->lsId;
            Eigen::Matrix<double,DOFS,DOFS> loc_dF = dF.block(i*DOFS,j*DOFS,DOFS,DOFS);
            ls.AddLHS(row, col, loc_dF);
        }
    }
}


void icy::Element::EvaluateStresses(icy::SimParams &prms,
                                    Eigen::Matrix3d &elasticityMatrix,
                                    Eigen::Matrix2d &D_mats)
{
    Eigen::Matrix<double,DOFS*3,1> u;
    u << nds[0]->ut, nds[1]->ut, nds[2]->ut;

    Eigen::Matrix<double,3,DOFS*3> bmat_b;
    Eigen::Matrix<double,2,DOFS*3> bmat_s[3];   // 3 gauss points
    Eigen::Matrix<double,3,DOFS*3> bmat_m;
    Eigen::Matrix<double,DOFS*3,DOFS*3> K;    // element stiffness matrix (3 gauss points)
    ComputeMatrices(prms, elasticityMatrix, D_mats, bmat_b, bmat_s, bmat_m, K);

    Eigen::Vector3d str_b, str_m, str_b_top;
    Eigen::Vector2d str_s[3];

    double thickness = prms.Thickness;
    str_b_top = str_b = elasticityMatrix*bmat_b*u;
    str_b*= (thickness*thickness*thickness/12.0); // !!! this is bening moment, not stress
    str_b_top *= (-thickness/2);
    str_m = elasticityMatrix*bmat_m*u*thickness;
    for(int i=0;i<3;i++) str_s[i] = D_mats*bmat_s[i]*u*(thickness/3.0);

    // transfer to arrays
    for(int i=0;i<3;i++)
    {
        a_str_b[i] = (float)str_b.coeff(i);
        a_str_m[i] = (float)str_m.coeff(i);
        a_str_b_top[i] = (float)str_b_top.coeff(i);
        a_str_s[i][0] = (float)str_s[i].coeff(0);
        a_str_s[i][1] = (float)str_s[i].coeff(1);
    }

    // stress on top/bottom of the plate and corresponding principal stresses

    float sx = str_b_top.coeff(0) + str_m.coeff(0);
    float sy = str_b_top.coeff(1) + str_m.coeff(1);
    float txy = str_b_top.coeff(2) + str_m.coeff(2);
    str_top << sx, txy, txy, sy;

    float coeff1 = sqrt((sx-sy)*(sx-sy)+txy*txy*4.0);
    float s1 = 0.5*(sx+sy+coeff1);
    float s2 = 0.5*(sx+sy-coeff1);
    float max_principal_top = std::max(s1,s2);

    sx = -str_b_top.coeff(0) + str_m.coeff(0);
    sy = -str_b_top.coeff(1) + str_m.coeff(1);
    txy = -str_b_top.coeff(2) + str_m.coeff(2);
    str_bottom << sx, txy, txy, sy;

    coeff1 = sqrt((sx-sy)*(sx-sy)+txy*txy*4.0);
    s1 = 0.5*(sx+sy+coeff1);
    s2 = 0.5*(sx+sy-coeff1);
    float max_principal_bottom = std::max(s1,s2);

    principal_stress_exceeds_threshold =
            std::max(max_principal_top, max_principal_bottom) > prms.normal_traction_threshold*prms.cutoff_coefficient;

    ComputeNormal();
}
