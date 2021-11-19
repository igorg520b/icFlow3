#include "element.h"
#include "edge.h"
#include "model.h"
#include "node.h"
#include "boundaryedge.h"

#include <cstdlib>
#include <algorithm>

#include <spdlog/spdlog.h>

// barycentric coords of Gauss points
double icy::Element::N[3][3] = {
       {2.0/3.0, 1.0/6.0, 1.0/6.0},
       {1.0/6.0, 2.0/3.0, 1.0/6.0},
       {1.0/6.0, 1.0/6.0, 2.0/3.0}};


void icy::Element::Reset()
{
    nds[0] = nds[1] = nds[2] = nullptr;
    incident_elems[0] = incident_elems[1] = incident_elems[2] = nullptr;
    area_initial = 0;
    normal_initial.setZero();
    normal_n.setZero();
}


void icy::Element::Initialize(Node *nd0, Node *nd1, Node *nd2)
{
    nds[0] = nd0;
    nds[1] = nd1;
    nds[2] = nd2;
}


void icy::Element::PrecomputeInitialArea()
{
    // translate the element
    Eigen::Vector3d p1, p2;
    p1 = nds[1]->x_initial - nds[0]->x_initial;
    p2 = nds[2]->x_initial - nds[0]->x_initial;

    normal_initial = p1.cross(p2);
    area_initial = normal_initial.norm()/2;
    if(area_initial<degenerate_area_threshold)
    {
        spdlog::critical("element {} - {} - {}",nds[0]->locId,nds[1]->locId,nds[2]->locId);
        spdlog::critical("area {}", area_initial);
        spdlog::critical("nd0 ({},{}); nd1 ({},{}); nd2 ({},{})",nds[0]->x_initial.x(),nds[0]->x_initial.y(),
                nds[1]->x_initial.x(),nds[1]->x_initial.y(),
                nds[2]->x_initial.x(),nds[2]->x_initial.y());
        throw std::runtime_error("degenerate element created");
    }
    normal_initial.normalize();
}

void icy::Element::UpdateSparseSystemEntries(LinearSystem &ls)
{
    int idxs[3] {nds[0]->lsId, nds[1]->lsId, nds[2]->lsId};
    ls.AddEntriesToStructure(std::begin(idxs),std::end(idxs));
}

void icy::Element::ComputeElasticForce(icy::LinearSystem &ls, icy::SimParams &prms, double,
                                       Eigen::Matrix3d &elasticityMatrix,
                                       Eigen::Matrix2d &D_mats)
{
    if(nds[0]->lsId < 0 && nds[1]->lsId < 0 && nds[2]->lsId < 0) return;

    Eigen::Matrix<double,3,LinearSystem::DOFS*3> bmat_b;
    Eigen::Matrix<double,2,LinearSystem::DOFS*3> bmat_s[3];   // 3 gauss points
    Eigen::Matrix<double,3,LinearSystem::DOFS*3> bmat_m;
    Eigen::Matrix<double,LinearSystem::DOFS*3,LinearSystem::DOFS*3> K;    // element stiffness matrix (3 gauss points)
    ComputeMatrices(prms, elasticityMatrix, D_mats, bmat_b, bmat_s, bmat_m, K);

    Eigen::Matrix<double,LinearSystem::DOFS*3,1> un;
    Eigen::Matrix<double,LinearSystem::DOFS*3,1> ut;
    un << nds[0]->un, nds[1]->un, nds[2]->un;
    ut << nds[0]->ut, nds[1]->ut, nds[2]->ut;

    // calculate elastic forces and Hessian at step n+1

    // absolute position is not important, the origin is set at node 0
    Eigen::Matrix<double,LinearSystem::DOFS*3,1> Fn, Fnp1;     // internal force at steps n and n+1
    Eigen::Matrix<double,LinearSystem::DOFS*3,LinearSystem::DOFS*3> &dFnp1 = K;

    Fn = K*un;
    Fnp1 = K*ut;

    // combine the results into a linearized equation of motion with HHT-alpha integration scheme
    double alpha = prms.HHTalpha;
    Eigen::Matrix<double,LinearSystem::DOFS*3,1> F;      // right-hand side of the equation is equal to -F
    Eigen::Matrix<double,LinearSystem::DOFS*3,LinearSystem::DOFS*3> dF;

    F = Fn*alpha + Fnp1*(1-alpha);
    dF= dFnp1*(1-alpha);

    // assert
    for(int i=0;i<LinearSystem::DOFS_SQ*3*3;i++)
            if(std::isnan(dF.data()[i]))
                throw std::runtime_error("icy::Element::ComputeElasticForce: dF contains NaN");

    // assemble
    ls.AddToEquation(F.data(),dF.data(),{nds[0]->lsId,nds[1]->lsId,nds[2]->lsId});
}


void icy::Element::ComputeMatrices(SimParams &prms,
                     Eigen::Matrix3d &elasticityMatrix,
                     Eigen::Matrix2d &D_mats,
                     Eigen::Matrix<double,3,LinearSystem::DOFS*3> &bmat_b,
                     Eigen::Matrix<double,2,LinearSystem::DOFS*3> (&bmat_s)[3],
                     Eigen::Matrix<double,3,LinearSystem::DOFS*3> &bmat_m,
                     Eigen::Matrix<double,LinearSystem::DOFS*3,LinearSystem::DOFS*3> &K)
{
    // translate the element
    Eigen::Vector3d p1, p2;
    p1 = nds[1]->x_initial - nds[0]->x_initial;
    p2 = nds[2]->x_initial - nds[0]->x_initial;

    normal_initial = p1.cross(p2);
    area_initial = normal_initial.norm()/2;
    normal_initial.normalize();
    if(area_initial<degenerate_area_threshold)
        throw std::runtime_error("icy::Element::ComputeMatrices: degenerate element detected");

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

    K = bmat_m.transpose()*elasticityMatrix*bmat_m*(area_initial*thickness);

    double coeff = area_initial*thickness*thickness*thickness/12.0;
    K += bmat_b.transpose()*elasticityMatrix*bmat_b*coeff;

    for(int i=0;i<3;i++) K += bmat_s[i].transpose()*D_mats*bmat_s[i]*(thickness*area_initial/3.0);
}


void icy::Element::EvaluateStresses(icy::SimParams &prms,
                                    Eigen::Matrix3d &elasticityMatrix,
                                    Eigen::Matrix2d &D_mats)
{
    Eigen::Matrix<double,LinearSystem::DOFS*3,1> u;
    u << nds[0]->ut, nds[1]->ut, nds[2]->ut;

    Eigen::Matrix<double,3,LinearSystem::DOFS*3> bmat_b;
    Eigen::Matrix<double,2,LinearSystem::DOFS*3> bmat_s[3];   // 3 gauss points
    Eigen::Matrix<double,3,LinearSystem::DOFS*3> bmat_m;
    Eigen::Matrix<double,LinearSystem::DOFS*3,LinearSystem::DOFS*3> K;    // element stiffness matrix (3 gauss points)
    ComputeMatrices(prms, elasticityMatrix, D_mats, bmat_b, bmat_s, bmat_m, K);

    Eigen::Vector3d str_b, str_m, str_b_top;
    Eigen::Vector2d str_s[3];

    double thickness = prms.Thickness;
    str_b_top = str_b = elasticityMatrix*bmat_b*u;
    str_b *= (thickness*thickness*thickness/12.0); // !!! this is bening moment, not stress
    str_b_top *= (-thickness/2);
    str_m = elasticityMatrix*bmat_m*u*thickness;
    for(int i=0;i<3;i++) str_s[i] = D_mats*bmat_s[i]*u*(thickness/3.0);

    // transfer to arrays
    for(int i=0;i<3;i++)
    {
        a_str_b[i] = str_b[i];
        a_str_m[i] = str_m[i];
        a_str_b_top[i] = str_b_top[i];
        a_str_s[i][0] = str_s[i].coeff(0);
        a_str_s[i][1] = str_s[i].coeff(1);
    }

    // stress on top/bottom of the plate and corresponding principal stresses

    double sx = str_b_top.coeff(0) + str_m.coeff(0);
    double sy = str_b_top.coeff(1) + str_m.coeff(1);
    double txy = str_b_top.coeff(2) + str_m.coeff(2);
    str_top << sx, txy, txy, sy;

    double coeff1 = sqrt((sx-sy)*(sx-sy)+txy*txy*4.0);
    double s1 = 0.5*(sx+sy+coeff1);
    double s2 = 0.5*(sx+sy-coeff1);
    double max_principal_top = std::max(s1,s2);

    sx = -str_b_top.coeff(0) + str_m.coeff(0);
    sy = -str_b_top.coeff(1) + str_m.coeff(1);
    txy = -str_b_top.coeff(2) + str_m.coeff(2);
    str_bottom << sx, txy, txy, sy;

    coeff1 = sqrt((sx-sy)*(sx-sy)+txy*txy*4.0);
    s1 = 0.5*(sx+sy+coeff1);
    s2 = 0.5*(sx+sy-coeff1);
    double max_principal_bottom = std::max(s1,s2);

    principal_stress_exceeds_threshold =
            std::max(max_principal_top, max_principal_bottom) > prms.FractureTractionThreshold*prms.cutoff_coefficient;

    ComputeNormal();
}

void icy::Element::ComputeNormal()
{
    Eigen::Vector3d u = nds[1]->xn3() - nds[0]->xn3();
    Eigen::Vector3d v = nds[2]->xn3() - nds[0]->xn3();
    normal_n = u.cross(v);
}


void icy::Element::DistributeStresses()
{
    // distribute the values from elements to nodes
    for(int i=0;i<3;i++)
    {
        Node *nd = nds[i];
        double coeff1 = 1.0/nd->adj_elems.size();//area_initial/nd->area;
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

        double str_s_combined[2] = {};
        for(int j=0;j<3;j++)
        {
            str_s_combined[0] += a_str_s[j][0]*coeff2*N[i][j];
            str_s_combined[1] += a_str_s[j][1]*coeff2*N[i][j];
        }
#pragma omp atomic
            nd->str_s[0] += str_s_combined[0];
#pragma omp atomic
            nd->str_s[1] += str_s_combined[1];
    }

    if(principal_stress_exceeds_threshold)
        for(int k=0;k<3;k++) nds[k]->potentially_can_fracture = true;
}





// GEOMETRICAL HELPER FUNCTIONS

uint8_t icy::Element::getNodeIdx(const Node *nd) const
{
    if(nds[0]==nd) return 0;
    else if(nds[1]==nd) return 1;
    else if(nds[2]==nd) return 2;
    else
    {
        spdlog::info("getNodeIdx: node not found; *nd = {}", (void*)nd);
        spdlog::info("getNodeIdx: nd->locId {}", nd->locId);
        spdlog::info("getNodeIdx: elem nds {} - {} - {}", (void*)nds[0],(void*)nds[1],(void*)nds[2]);
        spdlog::info("getNodeIdx: elem nds locid {} - {} - {}", nds[0]->locId,nds[1]->locId,nds[2]->locId);
        throw std::runtime_error("getNodeIdx");
    }
}

uint8_t icy::Element::getEdgeIdx(const Node *nd1, const Node *nd2) const
{
    for(uint8_t i=0;i<3;i++)
        if((nds[(i+1)%3]==nd1 && nds[(i+2)%3]==nd2) || (nds[(i+2)%3]==nd1 && nds[(i+1)%3]==nd2))
            return i;
    throw std::runtime_error("icy::Element::getEdgeIdx: edge not found");
}

std::pair<icy::Node*,icy::Node*> icy::Element::CW_CCW_Node(const Node* nd) const
{
    uint8_t idx = getNodeIdx(nd);
    return {nds[(idx+1)%3],nds[(idx+2)%3]};
}

bool icy::Element::isOnBoundary(const Node* nd) const
{
    uint8_t idx = getNodeIdx(nd);
    uint8_t cw_idx = (idx+1)%3;
    uint8_t ccw_idx = (idx+2)%3;
    return isBoundaryEdge(cw_idx) || isBoundaryEdge(ccw_idx);
}

bool icy::Element::isCWBoundary(const Node* nd) const
{
    uint8_t idx = getNodeIdx(nd);
    uint8_t cw_idx = (idx+2)%3;
    bool result = isBoundaryEdge(cw_idx);
    return result;
}

bool icy::Element::isCCWBoundary(const Node* nd) const
{
    uint8_t idx = getNodeIdx(nd);
    uint8_t ccw_idx = (idx+1)%3;
    return isBoundaryEdge(ccw_idx);
}

bool icy::Element::isEdgeCW(const Node *nd1, const Node *nd2) const
{
    for(int i=0;i<3;i++)
    {
        if(nds[i%3]==nd1 && nds[(i+2)%3]==nd2) return true;
        if(nds[i%3]==nd1 && nds[(i+1)%3]==nd2) return false;
    }
    throw std::runtime_error("icy::Element::isEdgeCW: edge not found");
}

bool icy::Element::containsEdge(const Node *nd1, const Node *nd2) const
{
    for(uint8_t i=0;i<3;i++)
        if((nds[(i+1)%3]==nd1 && nds[(i+2)%3]==nd2) || (nds[(i+2)%3]==nd1 && nds[(i+1)%3]==nd2))
            return true;
    return false;
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

    spdlog::critical("getOppositeNode: trying to find edge {}-{} in element {}-{}-{}",
                     nd0->locId,nd1->locId,nds[0]->locId,nds[1]->locId,nds[2]->locId);
    throw std::runtime_error("getOppositeNode: opposite node not found");
}

icy::BaseElement* icy::Element::getIncidentElementOppositeToNode(Node *nd)
{
    return incident_elems[getNodeIdx(nd)];
}


// FRACTURE ALGORITHM


void icy::Element::ReplaceNode(Node *replaceWhat, Node *replaceWith)
{
    uint8_t nd_idx = getNodeIdx(replaceWhat);
    nds[nd_idx] = replaceWith;

    PrecomputeInitialArea();

    uint8_t cw_idx = (nd_idx+1)%3;
    uint8_t ccw_idx = (nd_idx+2)%3;

    // update incident elements after replacing the node
    for(uint8_t idx : {cw_idx,ccw_idx}) incident_elems[idx]->UpdateNodes();
}

void icy::Element::ReplaceIncidentElem(const BaseElement* which, BaseElement* withWhat)
{
    if(incident_elems[0] == which) incident_elems[0] = withWhat;
    else if(incident_elems[1] == which) incident_elems[1] = withWhat;
    else if(incident_elems[2] == which) incident_elems[2] = withWhat;
    else throw std::runtime_error("ReplaceIncidentElem: incident elem not found");
}

void icy::Element::ReplaceAdjacentElem(const Element* originalElem, Element* insertedElem, uint8_t idx)
{
    ReplaceIncidentElem(originalElem,insertedElem);
}




std::pair<icy::Node*,icy::Node*> icy::Element::SplitElem(Node *nd, Node *nd0, Node *nd1, double where)
{
    BaseElement* incident_elem = getIncidentElementOppositeToNode(nd);

    if(incident_elem->type == ElementType::BEdge)
        return {SplitBoundaryElem(nd, nd0, nd1, where),nullptr};
    else if(incident_elem->type == ElementType::TElem)
        return {SplitNonBoundaryElem(nd, nd0, nd1, where),nullptr};
    else throw std::runtime_error("SplitElem: unknown incident element type");
}


icy::Node* icy::Element::SplitBoundaryElem(Node *nd, Node *nd0, Node *nd1, double where)
{
    MeshFragment *fragment = nd->fragment;

    uint8_t ndIdx = getNodeIdx(nd);
    uint8_t nd0Idx = getNodeIdx(nd0);
    uint8_t nd1Idx = getNodeIdx(nd1);


    // insert element
    Element *insertedElem = fragment->AddElement();
    nd->adj_elems.push_back(insertedElem);

    // insert the node between nd0 and nd1; initialize its coordinates; connect to adjacent elements
    Node *split = fragment->AddNode();
    split->InitializeLERP(nd0, nd1, where);
    split->isBoundary = true;
    split->adj_elems.push_back(this);
    split->adj_elems.push_back(insertedElem);

    // modify the original element
    nds[nd1Idx] = split;

    // initialize the inserted element's nodes
    insertedElem->nds[ndIdx] = nd;
    insertedElem->nds[nd1Idx] = nd1;
    insertedElem->nds[nd0Idx] = split;

    // if the original element had a boundary at nd0Idx, the inserted element now takes that boundary
    incident_elems[nd0Idx]->ReplaceAdjacentElem(this, insertedElem, nd0Idx);


    // add the boundary that has just appeared
    // automatically initialize the inserted element's adjacency data (insertedElem->incident_elems[ndIdx])
    fragment->AddBoundary(insertedElem,ndIdx,3);
    incident_elems[ndIdx]->UpdateNodes();

    insertedElem->incident_elems[nd0Idx] = incident_elems[nd0Idx];
    incident_elems[nd0Idx] = insertedElem;
    insertedElem->incident_elems[nd1Idx] = this;

    // from node "nd1", disconnect the original element and replace it with the inserted element
    nd1->ReplaceAdjacentElement(this, insertedElem);

    // compute the new area and reference shape matrix
    this->PrecomputeInitialArea();
    insertedElem->PrecomputeInitialArea();

    // re-evaluate PiMultiplier on both elements to maintain consistent plasticity
    this->RecalculatePiMultiplierFromDeformationGradient(F_orig);
    insertedElem->RecalculatePiMultiplierFromDeformationGradient(F_orig);

    // "fix" the fan for the node, whose element was just replaced
    nd1->PrepareFan();
    return split;
}

icy::Node* icy::Element::SplitNonBoundaryElem(Node *nd, Node *nd0, Node *nd1, double where)
{
    MeshFragment *fragment = nd->fragment;

    icy::Element *adjElem = dynamic_cast<icy::Element*>(getIncidentElementOppositeToNode(nd));
    if(adjElem == nullptr) throw std::runtime_error("icy::Element::SplitNonBoundaryElem dynamic cast issue");

    uint8_t ndIdx_orig = getNodeIdx(nd);
    uint8_t nd0Idx_orig = getNodeIdx(nd0);
    uint8_t nd1Idx_orig = getNodeIdx(nd1);

    // preserve deformation gradient
    Eigen::Matrix2d F_orig = getF_at_n();
    Eigen::Matrix2d F_adj = adjElem->getF_at_n();

    Node *oppositeNode = adjElem->getOppositeNode(nd0, nd1);
    uint8_t nd0Idx_adj = adjElem->getNodeIdx(nd0);
    uint8_t nd1Idx_adj = adjElem->getNodeIdx(nd1);
    uint8_t oppIdx_adj = adjElem->getNodeIdx(oppositeNode);

    // insert "main" element
    Element *insertedElem = fragment->AddElement();
    nd->adj_elems.push_back(insertedElem);

    // insert "adjacent" element
    Element *insertedElem_adj = fragment->AddElement();

    // insert the "split" node between nd0 and nd1
    Node *split = fragment->AddNode();
    split->InitializeLERP(nd0, nd1, where);
    split->adj_elems.insert(split->adj_elems.end(),{this,insertedElem,adjElem,insertedElem_adj});
    split->isBoundary = false;

    // modify the original element
    nds[nd1Idx_orig] = split;

    nd1->ReplaceAdjacentElement(this,insertedElem);

    // initialize the inserted "main" element
    insertedElem->nds[ndIdx_orig] = nd;
    insertedElem->nds[nd1Idx_orig] = nd1;
    insertedElem->nds[nd0Idx_orig] = split;

    // if the original element had a boundary at nd0Idx, the inserted element now takes that boundary
    incident_elems[nd0Idx_orig]->ReplaceAdjacentElem(this, insertedElem, nd0Idx_orig);


    insertedElem->incident_elems[ndIdx_orig] = insertedElem_adj;
    insertedElem->incident_elems[nd0Idx_orig] = incident_elems[nd0Idx_orig];
    insertedElem->incident_elems[nd1Idx_orig] = this;
    incident_elems[nd0Idx_orig] = insertedElem;

    // similarly, modify the existing adjacent element
    adjElem->nds[nd1Idx_adj] = split;
    insertedElem_adj->nds[oppIdx_adj] = oppositeNode;
    insertedElem_adj->nds[nd1Idx_adj] = nd1;
    insertedElem_adj->nds[nd0Idx_adj] = split;

    // whichever element was attached at nd0Idx is now attached to insertedElem_adj
    adjElem->incident_elems[nd0Idx_adj]->ReplaceAdjacentElem(adjElem, insertedElem_adj, nd0Idx_adj);


    insertedElem_adj->incident_elems[oppIdx_adj] = insertedElem;
    insertedElem_adj->incident_elems[nd1Idx_adj] = adjElem;
    insertedElem_adj->incident_elems[nd0Idx_adj] = adjElem->incident_elems[nd0Idx_adj];
    adjElem->incident_elems[nd0Idx_adj] = insertedElem_adj;

    oppositeNode->adj_elems.push_back(insertedElem_adj);
    nd1->ReplaceAdjacentElement(adjElem,insertedElem_adj);

    this->PrecomputeInitialArea();
    insertedElem->PrecomputeInitialArea();
    adjElem->PrecomputeInitialArea();
    insertedElem_adj->PrecomputeInitialArea();

    // "fix" palsticity on all four elements
    this->RecalculatePiMultiplierFromDeformationGradient(F_orig);
    insertedElem->RecalculatePiMultiplierFromDeformationGradient(F_orig);
    adjElem->RecalculatePiMultiplierFromDeformationGradient(F_adj);
    insertedElem_adj->RecalculatePiMultiplierFromDeformationGradient(F_adj);

    oppositeNode->PrepareFan();
    nd1->PrepareFan();
    return split;
}





