#ifndef ELEMENT123_H
#define ELEMENT123_H

#include <iostream>
#include <vector>
#include <utility>

#include <QtDebug>
#include <QtGlobal>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "parameters_sim.h"
#include "linearsystem.h"
#include "node.h"

namespace icy { class Element; class Node; class SimParams; class Edge; }

class icy::Element
{
public:
    icy::Node* nds[3];          // initialized when the geometry is loaded or remeshed
    icy::Edge edges[3];        // element's edges 0-1; 1-2; 2-0;
    icy::Element* adj_elems[3]; // nullptr of no adjacent element
    unsigned short region;
    unsigned short traversal;             // for traversal when identifying region connectivity

    // at initial state
    double area_initial;
    Eigen::Vector3d normal_initial, normal_n;
    bool initial_normal_up; // normal_initial.z() > 0
//    Eigen::Vector3d pr1_initial, pr2_initial;
//    Eigen::Matrix<double,DOFS*3,1> x_initial;

    // strain-displacement for bending, shear and membrane
    Eigen::Matrix<double,3,DOFS*3> bmat_b;
    Eigen::Matrix<double,2,DOFS*3> bmat_s[3];   // 3 gauss points
    Eigen::Matrix<double,3,DOFS*3> bmat_m;
    Eigen::Matrix<double,DOFS*3,DOFS*3> K;    // element stiffness matrix (3 gauss points)

    // sress values / visualization
    Eigen::Vector3d str_b, str_m, str_b_top;
    Eigen::Vector2d str_s[3];
    Eigen::Matrix2f str_top, str_bottom;    // stress on top and bottom surface of the plate, as 3x3 matrix
    bool principal_stress_exceeds_threshold;

    void InitializePersistentVariables(); // everything that depends on initial position (not K and M)
    // K and M are computed after rho, Y and nu are set/changed
    void PrecomputeStiffnessMatrix(SimParams &prms,
                                   Eigen::Matrix3d &elasticityMatrix,
                                   Eigen::Matrix2d &D_mats);
    // compute forces and insert nz-entries to sparse structure
    void UpdateSparseSystemEntries(LinearSystem &ls);
    void ComputeElasticForce(LinearSystem &ls, SimParams &prms, double timeStep);

    // this is done after the new displacement values are accepted
    void EvaluateStresses(SimParams &prms,
                          Eigen::Matrix3d &elasticityMatrix,
                          Eigen::Matrix2d &D_mats);
    void DistributeStresses();

    // helper functions for fracture
    icy::Node* getOppositeNode(Edge edge);    // return the node across from a given edge
    icy::Node* getOppositeNode(Node *nd0, Node* nd1);
//    std::pair<int,int> getOppositeEdge(Node *nd);
    Eigen::Vector3d getCenter();

    void getIdxs(Node*nd, short &thisIdx, short &CWIdx, short &CCWIdx);
    Edge getEdgeOppositeToNode(Node *nd);
    Element* getAdjacentElementOppositeToNode(Node *nd);
    short getNodeIdx(Node *nd);


    bool ContainsNode(Node *nd){return (nds[0]==nd || nds[1]==nd || nds[2]==nd);}
    void ReplaceNode(Node *replaceWhat, Node *replaceWith);
//    void Initialize(Node* nd0, Node* nd1, Node* nd2, bool orientation);
//    bool Orientation(const Node* nd0, const Node* nd1);
    void ComputeNormal();
    void AssertEdges();

private:
    static double N[3][3];
    void rotationMatrix(Eigen::Vector3d &p1, Eigen::Vector3d &p2, Eigen::Matrix3d &result,
                        double &area, Eigen::Vector3d &normal);
    void rotationMatrix_alt(Eigen::Vector3d &p1, Eigen::Vector3d &p2,
                            Eigen::Matrix3d &result,
                            double &area, Eigen::Vector3d &normal);

    // compute the elastic force and Hessian based on given location of vertices
    // if dFo is null, don't compute it
    void FdF(Eigen::Matrix<double,DOFS*3,1> &p,
            Eigen::Matrix<double,DOFS*3,1> &Fo,
            Eigen::Matrix<double,DOFS*3,DOFS*3> *dFo);

};

#endif // ELEMENT123_H
