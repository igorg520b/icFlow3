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

void icy::Element::InitializePersistentVariables()
{
    x_initial << nds[0]->x_initial, nds[1]->x_initial, nds[2]->x_initial;
    // translate the element
    Eigen::Vector3d p1, p2;
    p1 = nds[1]->x_initial.block(0,0,3,1) - nds[0]->x_initial.block(0,0,3,1);
    p2 = nds[2]->x_initial.block(0,0,3,1) - nds[0]->x_initial.block(0,0,3,1);

    rotationMatrix(p1, p2, R0, area_initial, normal_initial);
    initial_normal_up = normal_initial.z() > 0;
    for(int j=0;j<3;j++) nds[j]->area += area_initial/3; // distribute area to adjacent nodes

    R0t = R0.transpose();
    pr1_initial = p1;
    pr2_initial = p2;

    // (use 1-based indices, yij = yi-yj)
    double x1 = 0;
    double x2 = pr1_initial.x();
    double x3 = pr2_initial.x();
    double y1 = 0;
    double y2 = pr1_initial.y();
    double y3 = pr2_initial.y();

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
    for(int i=0;i<3;i++) {
        bmat_s[i] <<
           0,0,dN1dx,-N[i][0],0,   0,0,dN2dx,-N[i][1],0,     0,0,dN3dx,-N[i][2], 0,
           0,0,dN1dy,0,-N[i][0],   0,0,dN2dy,0,-N[i][1],     0,0,dN3dy,0,-N[i][2];
    }
}

void icy::Element::PrecomputeStiffnessMatrix(icy::SimParams &prms,
                                             Eigen::Matrix3d &elasticityMatrix,
                                             Eigen::Matrix2d &D_mats)
{
    double thickness = prms.Thickness;
    // K and M depend on rho, Young's modulus and Poisson's ratio,
    // therefore they are computed after these parameters are set

    Eigen::Matrix<double,DOFS*3,DOFS*3> K_b, K_m, K_s;
    K = Eigen::Matrix<double,DOFS*3,DOFS*3>::Zero();
    K += bmat_m.transpose()*elasticityMatrix*bmat_m*(area_initial*thickness);

    double coeff = area_initial*thickness*thickness*thickness/12.0;
    K += bmat_b.transpose()*elasticityMatrix*bmat_b*coeff;

    for(int i=0;i<3;i++) K += bmat_s[i].transpose()*D_mats*bmat_s[i]*(thickness*area_initial/3.0);
}

void icy::Element::FdF(
        Eigen::Matrix<double,DOFS*3,1> &u,
        Eigen::Matrix<double,DOFS*3,1> &Fo,
        Eigen::Matrix<double,DOFS*3,DOFS*3> *dFo)
{
    Fo = K*u;
    if(dFo != nullptr) *dFo = K;
}

void icy::Element::ComputeElasticForce(icy::LinearSystem &ls, icy::SimParams &prms, double)
{
    if(nds[0]->lsId < 0 && nds[1]->lsId < 0 && nds[2]->lsId < 0) return;

    // reserve non-zero entries in the sparse matrix structure
    ls.AddElementToStructure(nds[0]->lsId, nds[1]->lsId);
    ls.AddElementToStructure(nds[0]->lsId, nds[2]->lsId);
    ls.AddElementToStructure(nds[1]->lsId, nds[2]->lsId);

    Eigen::Matrix<double,DOFS*3,1> un;
    Eigen::Matrix<double,DOFS*3,1> ut;
    un << nds[0]->un, nds[1]->un, nds[2]->un;
    ut << nds[0]->ut, nds[1]->ut, nds[2]->ut;

    // calculate elastic forces and Hessian at step n+1

    // absolute position is not important, the origin is set at node 0
    Eigen::Matrix<double,DOFS*3,1> Fn, Fnp1;     // internal force at steps n and n+1
    Eigen::Matrix<double,DOFS*3,DOFS*3> dFnp1;   // Hessian of the force function

    FdF(un, Fn, nullptr);
    FdF(ut, Fnp1, &dFnp1);

    // combine the results into a linearized equation of motion with HHT-alpha integration scheme
    double alpha = prms.HHTalpha;
    F = Fn*alpha + Fnp1*(1-alpha);
    dF= dFnp1*(1-alpha);

#ifdef QT_DEBUG
    // assert
    for(int i=0;i<DOFS*3;i++)
        for(int j=0;j<DOFS*3;j++)
            if(std::isnan(dF(i,j)))
                throw std::runtime_error("elem.ComputeElasticForce: dF contains NaN");
#endif
}

void icy::Element::Assemble(icy::LinearSystem &ls) const
{
    if(nds[0]->lsId < 0 && nds[1]->lsId < 0 && nds[2]->lsId < 0) return;
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

    double thickness = prms.Thickness;
    double coeff = thickness*thickness*thickness/12.0;
    str_b = elasticityMatrix*bmat_b*u*coeff;
    str_m = -elasticityMatrix*bmat_m*u*thickness;
    for(int i=0;i<3;i++) str_s[i] = D_mats*bmat_s[i]*u*(thickness/3.0);

    // stress on top/bottom of the plate and corresponding principal stresses
    str_b_top = elasticityMatrix*bmat_b*u*area_initial*(-thickness/2);

    float sx = (float)str_b_top.coeff(0);
    float sy = str_b_top.coeff(1);
    float txy = str_b_top.coeff(2);
    str_top << sx, txy, txy, sy;

    str_b_bottom = elasticityMatrix*bmat_b*u*area_initial*(thickness/2);

    sx = str_b_bottom.coeff(0) + str_m(0);
    sy = str_b_bottom.coeff(1) + str_m(1);
    txy = str_b_bottom.coeff(2) + str_m(2);
    str_bottom << sx, txy, txy, sy;

    ComputeNormal();
}

void icy::Element::DistributeStresses()
{
    // distribute the values from elements to nodes
    for(int i=0;i<3;i++) {
        Node *nd = nds[i];
        double coeff1 = area_initial/nd->area;
        double coeff2 = area_initial/(3*nd->area);

        nd->str_b += str_b*coeff1;
        nd->str_m += str_m*coeff1;

        for(int j=0;j<3;j++) nd->str_s += str_s[j]*coeff2*N[i][j];

        nd->str_b_top += str_b_top*coeff1;
        nd->str_b_bottom += str_b_bottom*coeff1;

        nd->normal_n += normal_n;
    }
}

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


void icy::Element::ComputeNormal()
{
    Eigen::Vector3d u = (nds[1]->xn - nds[0]->xn).block(0,0,3,1);
    Eigen::Vector3d v = (nds[2]->xn - nds[0]->xn).block(0,0,3,1);
    normal_n = u.cross(v);
}



//=============== fracture helpers
icy::Node* icy::Element::getOppositeNode(icy::Edge *edge)
{
    icy::Node *nd0 = edge->nds[0];
    icy::Node *nd1 = edge->nds[1];

    for(int i=0;i<3;i++)
    {
        int idx_next = (i+1)%3;
        if((nds[i] == nd0 && nds[idx_next] == nd1)||
                (nds[i] == nd1 && nds[idx_next] == nd0))
            return nds[(i+2)%3];
    }
    throw std::runtime_error("opposite node not found");
}

std::pair<int,int> icy::Element::getOppositeEdge(icy::Node *nd)
{
    int nd_idx;
    if(nd == nds[0]) nd_idx = 0;
    else if(nd == nds[1]) nd_idx = 1;
    else if(nd == nds[2]) nd_idx = 2;
    else throw std::runtime_error("node not found");

    int edge_idx0 = nds[(nd_idx+1)%3]->locId;
    int edge_idx1 = nds[(nd_idx+2)%3]->locId;
    if(edge_idx0<edge_idx1) return std::make_pair(edge_idx0, edge_idx1);
    else return std::make_pair(edge_idx1, edge_idx0);
}

void icy::Element::ReplaceNode(icy::Node *replaceWhat, icy::Node *replaceWith)
{
    if(nds[0] == replaceWhat) { nds[0] = replaceWith; return; }
    if(nds[1] == replaceWhat) { nds[1] = replaceWith; return; }
    if(nds[2] == replaceWhat) { nds[2] = replaceWith; return; }
    throw std::runtime_error("replaceWhat is not in nds[]");
}

void icy::Element::Initialize(icy::Node* nd0, icy::Node* nd1, icy::Node* nd2, bool orientation)
{
    nds[0] = nd0;
    if(orientation) {
        nds[1] = nd1;
        nds[2] = nd2;
    } else {
        nds[2] = nd1;
        nds[1] = nd2;
    }
    InitializePersistentVariables();
}

bool icy::Element::Orientation(icy::Node* nd0, icy::Node* nd1)
{
    if(nds[0] == nd0 && nds[1] == nd1) return true;
    if(nds[1] == nd0 && nds[0] == nd1) return false;

    if(nds[1] == nd0 && nds[2] == nd1) return true;
    if(nds[2] == nd0 && nds[1] == nd1) return false;

    if(nds[2] == nd0 && nds[0] == nd1) return true;
    if(nds[0] == nd0 && nds[2] == nd1) return false;

    throw std::runtime_error("nodes do not belong to this element");
}

Eigen::Vector3d icy::Element::getCenter()
{
    return (nds[0]->x_initial + nds[1]->x_initial + nds[2]->x_initial).block(0,0,3,1)/3.0;
}

icy::Node* icy::Element::getCCWNode(icy::Node* nd)
{
    if(normal_initial.z()<0) {
        if(nd==nds[0]) return nds[1];
        if(nd==nds[1]) return nds[2];
        if(nd==nds[2]) return nds[0];
    } else
    {
        if(nd==nds[0]) return nds[2];
        if(nd==nds[1]) return nds[0];
        if(nd==nds[2]) return nds[1];
    }
    throw std::runtime_error("nd not found");
}

icy::Node* icy::Element::getCWNode(icy::Node* nd)
{
    if(normal_initial.z()<0) {

        if(nd==nds[0]) return nds[2];
        if(nd==nds[1]) return nds[0];
        if(nd==nds[2]) return nds[1];
    } else {
        if(nd==nds[0]) return nds[1];
        if(nd==nds[1]) return nds[2];
        if(nd==nds[2]) return nds[0];
    }
    throw std::runtime_error("nd not found");
}

short icy::Element::getCCWIdx(icy::Node* nd)
{
    if(!initial_normal_up) {
        if(nd==nds[0]) return 1;
        if(nd==nds[1]) return 2;
        if(nd==nds[2]) return 0;
    } else
    {
        if(nd==nds[0]) return 2;
        if(nd==nds[1]) return 0;
        if(nd==nds[2]) return 1;
    }
    throw std::runtime_error("nd not found");
}

short icy::Element::getCWIdx(icy::Node* nd)
{
    if(!initial_normal_up) {
        if(nd==nds[0]) return 2;
        if(nd==nds[1]) return 0;
        if(nd==nds[2]) return 1;
    } else {
        if(nd==nds[0]) return 1;
        if(nd==nds[1]) return 2;
        if(nd==nds[2]) return 0;
    }
    throw std::runtime_error("nd not found");
}
