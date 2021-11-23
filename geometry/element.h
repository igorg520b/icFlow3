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
#include "baseelement.h"


namespace icy { class Element; class Node; class SimParams; class Edge; }

class icy::Element : public icy::BaseElement
{
public:
    icy::Node* nds[3];          // initialized when the geometry is loaded or remeshed
    icy::BaseElement* incident_elems[3]; // nullptr or an adjacent element
    unsigned short region;
    unsigned short traversal;             // for traversal when identifying region connectivity

    // at initial state
    double area_initial;
    Eigen::Vector3d normal_initial, normal_n;

    // sress values / visualization
    double a_str_b[3], a_str_m[3], a_str_b_top[3];
    double a_str_s[3][2];

    Eigen::Matrix2d str_top, str_bottom;    // stress on top and bottom surface of the plate, as 3x3 matrix
    bool principal_stress_exceeds_threshold;

    Element() { type = ElementType::TElem; }
    ~Element() = default;
    Element& operator=(const Element&) = delete;

    void Reset();
    void Initialize(Node *nd0, Node *nd1, Node *nd2);
    void PrecomputeInitialArea();


    // ELASTIC FORCES AND STRESSES
public:

    void UpdateSparseSystemEntries(LinearSystem &ls);
    void ComputeElasticForce(LinearSystem &ls, SimParams &prms, double timeStep,
                             Eigen::Matrix3d &elasticityMatrix,
                             Eigen::Matrix2d &D_mats);


    // this is done after the new displacement values are accepted
    void EvaluateStresses(SimParams &prms,
                          Eigen::Matrix3d &elasticityMatrix,
                          Eigen::Matrix2d &D_mats);
    void DistributeStresses();

private:
    static double N[3][3]; // barycentric coords of Gauss points
    constexpr static double degenerate_area_threshold = 1e-8;
    void ComputeNormal();
    void ComputeMatrices(SimParams &prms,
                         Eigen::Matrix3d &elasticityMatrix,
                         Eigen::Matrix2d &D_mats,
                         Eigen::Matrix<double,3,LinearSystem::DOFS*3> &bmat_b,
                         Eigen::Matrix<double,2,LinearSystem::DOFS*3> (&bmat_s)[3],
                         Eigen::Matrix<double,3,LinearSystem::DOFS*3> &bmat_m,
                         Eigen::Matrix<double,LinearSystem::DOFS*3,LinearSystem::DOFS*3> &K);


    // FRACTURE ALGORITHM
    public:
        std::pair<Node*,Node*> SplitElem(Node *nd, Node *nd0, Node *nd1, double where); // split the element by inserting a node between nd0 and nd1

        bool isBoundary() {return std::any_of(std::begin(nds),std::end(nds),[](Node *nd){return nd->isBoundary;});}

        uint8_t getNodeIdx(const Node* nd) const;
        uint8_t getEdgeIdx(const Node *nd1, const Node *nd2) const;
        std::pair<Node*,Node*> CW_CCW_Node(const Node* nd) const;

        bool isBoundaryEdge(const uint8_t idx) const {return incident_elems[idx]->type != ElementType::TElem;}
        bool isOnBoundary(const Node* nd) const;
        bool isCWBoundary(const Node* nd) const;
        bool isCCWBoundary(const Node* nd) const;
        bool isEdgeCW(const Node *nd1, const Node *nd2) const; // true if nd1-nd2 is oriented clockwise
        bool containsEdge(const Node *nd1, const Node *nd2) const;

        bool containsNode(const Node* nd) const {return (nds[0]==nd || nds[1]==nd || nds[2]==nd);}
        Eigen::Vector3d getCenter() const {return (nds[0]->x_initial + nds[1]->x_initial + nds[2]->x_initial)/3.0;};
        icy::Node* getOppositeNode(Node* nd0, Node* nd1);

        void ReplaceNode(Node* replaceWhat, Node* replaceWith);
        void ReplaceIncidentElem(const BaseElement* which, BaseElement* withWhat);
        void ReplaceAdjacentElem(const Element* originalElem, Element* insertedElem, uint8_t idx) override;

    private:
        constexpr static double threshold_area = 1e-7;
        Node* SplitBoundaryElem(Node *nd, Node *nd0, Node *nd1, double where);
        Node* SplitNonBoundaryElem(Node *nd, Node *nd0, Node *nd1, double where);

        BaseElement* getIncidentElementOppositeToNode(Node* nd);

};

#endif // ELEMENT123_H
