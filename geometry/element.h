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
    unsigned region;
    bool traversed;             // for traversal when identifying region connectivity

    // at initial state
    double area_initial;
    Eigen::Vector3d normal_initial, normal_n;
    bool initial_normal_up; // normal_initial.z() > 0
    Eigen::Vector3d pr1_initial, pr2_initial;
    Eigen::Matrix3d R0, R0t; // for testing
    Eigen::Matrix<double,DOFS*3,1> x_initial;

    // strain-displacement for bending, shear and membrane
    Eigen::Matrix<double,3,DOFS*3> bmat_b;
    Eigen::Matrix<double,2,DOFS*3> bmat_s[3];   // 3 gauss points
    Eigen::Matrix<double,3,DOFS*3> bmat_m;
    Eigen::Matrix<double,DOFS*3,DOFS*3> K;    // element stiffness matrix (3 gauss points)

    // linear system
    Eigen::Matrix<double,DOFS*3,1> F;      // right-hand side of the equation is equal to -F
    Eigen::Matrix<double,DOFS*3,DOFS*3> dF;

    // sress values / visualization
    Eigen::Vector3d str_b, str_m, str_b_top, str_b_bottom;
    Eigen::Vector2d str_s[3];
    Eigen::Matrix2f str_top, str_bottom;    // stress on top and bottom surface of the plate, as 3x3 matrix
    bool principal_stress_exceeds_threshold;

    void InitializePersistentVariables(); // everything that depends on initial position (not K and M)
    // K and M are computed after rho, Y and nu are set/changed
    void PrecomputeStiffnessMatrix(icy::SimParams &prms,
                                   Eigen::Matrix3d &elasticityMatrix,
                                   Eigen::Matrix2d &D_mats);
    // compute forces and insert nz-entries to sparse structure
    void ComputeElasticForce(icy::LinearSystem &ls, icy::SimParams &prms, double timeStep);
    void Assemble(icy::LinearSystem &ls) const;   // distributed computed F and dF into linear system

    // this is done after the new displacement values are accepted
    void EvaluateStresses(icy::SimParams &prms,
                          Eigen::Matrix3d &elasticityMatrix,
                          Eigen::Matrix2d &D_mats);
    void DistributeStresses();

    // helper functions for fracture
    icy::Node* getOppositeNode(icy::Edge *edge);    // return the node across from a given edge
    std::pair<int,int> getOppositeEdge(icy::Node *nd);
    Eigen::Vector3d getCenter();
//    icy::Node* getCWNode(icy::Node* nd);
//    icy::Node* getCCWNode(icy::Node* nd);
//    short getCWIdx(icy::Node* nd);
//    short getCCWIdx(icy::Node* nd);
    void getIdxs(icy::Node*nd, short &thisIdx, short &CWIdx, short &CCWIdx);
    icy::Edge getEdgeOppositeToNode(icy::Node *nd);

    bool ContainsNode(icy::Node *nd){return (nds[0]==nd || nds[1]==nd || nds[2]==nd);}
    void ReplaceNode(icy::Node *replaceWhat, icy::Node *replaceWith);
    void Initialize(icy::Node* nd0, icy::Node* nd1, icy::Node* nd2, bool orientation);
    bool Orientation(icy::Node* nd0, icy::Node* nd1);
    void ComputeNormal();

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
