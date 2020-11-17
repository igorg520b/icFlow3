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
    icy::Edge edges[3];        // element's edges opposite to nd 0,1,2
    icy::Element* adj_elems[3]; // nullptr of no adjacent element
    unsigned short region;
    unsigned short traversal;             // for traversal when identifying region connectivity

    // at initial state
    double area_initial;
    Eigen::Vector3d normal_initial, normal_n;

    // sress values / visualization
    float a_str_b[3], a_str_m[3], a_str_b_top[3];// str_b_top;
    float a_str_s[3][2];

    Eigen::Matrix2f str_top, str_bottom;    // stress on top and bottom surface of the plate, as 3x3 matrix
    bool principal_stress_exceeds_threshold;

    void ComputeInitialNormal();

    void UpdateSparseSystemEntries(LinearSystem &ls);
    void ComputeElasticForce(LinearSystem &ls, SimParams &prms, double timeStep,
                             Eigen::Matrix3d &elasticityMatrix,
                             Eigen::Matrix2d &D_mats);

    // this is done after the new displacement values are accepted
    void EvaluateStresses(SimParams &prms,
                          Eigen::Matrix3d &elasticityMatrix,
                          Eigen::Matrix2d &D_mats);
    void DistributeStresses();

    void ComputeMatrices(SimParams &prms,
                         Eigen::Matrix3d &elasticityMatrix,
                         Eigen::Matrix2d &D_mats,
                         Eigen::Matrix<double,3,DOFS*3> &bmat_b,
                         Eigen::Matrix<double,2,DOFS*3> (&bmat_s)[3],
                         Eigen::Matrix<double,3,DOFS*3> &bmat_m,
                         Eigen::Matrix<double,DOFS*3,DOFS*3> &K);

    // helper functions for fracture
    icy::Node* getOppositeNode(Edge edge);    // return the node across from a given edge
    icy::Node* getOppositeNode(Node *nd0, Node* nd1);
    Eigen::Vector3d getCenter();

    void getIdxs(const Node*nd, short &thisIdx, short &CWIdx, short &CCWIdx) const;
    Edge getEdgeOppositeToNode(Node *nd);
    Element* getAdjacentElementOppositeToNode(Node *nd);
    short getNodeIdx(Node *nd);

    Edge CWEdge(const Node* nd) const;
    Edge CCWEdge(const Node* nd) const;
    Edge OppositeEdge(const Node* nd) const;

    bool ContainsNode(Node *nd){return (nds[0]==nd || nds[1]==nd || nds[2]==nd);}
    void ReplaceNode(Node *replaceWhat, Node *replaceWith);
    void ComputeNormal();
    void AssertEdges();

private:
    static double N[3][3];
    /*
    void rotationMatrix(Eigen::Vector3d &p1, Eigen::Vector3d &p2, Eigen::Matrix3d &result,
                        double &area, Eigen::Vector3d &normal);
    void rotationMatrix_alt(Eigen::Vector3d &p1, Eigen::Vector3d &p2,
                            Eigen::Matrix3d &result,
                            double &area, Eigen::Vector3d &normal);
                            */
};

#endif // ELEMENT123_H
