#ifndef EDGE_H
#define EDGE_H

#include <Eigen/Core>

namespace icy {class Edge; class Node; class Element;}

class icy::Edge
{
public:
    void Initialize(icy::Node* nd0, icy::Node* nd1, icy::Element* elem0, icy::Element* elem1);

    icy::Node* nds[2];
    bool isBoundary;    // belongs to only one element
    icy::Element* elems[2];
    double angle0_initial; // [-pi, +pi] from nds[0] to nds[1]
    double angle1_initial; // [-pi, +pi] from nds[1] to nds[0]

    bool elem0_isCCW, elem1_isCCW;
    Eigen::Vector3d res1, res2;

    double getAngle(icy::Node *center_node) const;
    icy::Element* get_CCW_Element(icy::Node *center_node) const; // adjacent element from ccw side; deprecated
    Eigen::Vector2f getVec(icy::Node *center_node) const;     // as vector at step n
    icy::Node* getOtherNode(icy::Node *center_node) const;
    // icy::Node* getOppositeNode(int idx);
    icy::Element* getTheOnlyElement();

    icy::Element* getElementWithNode(icy::Node *nd);
    icy::Node* getFarNode(icy::Node *nd);

private:
    bool ElementBoundaryFollowsEdge(icy::Element* elem, icy::Node* &opposite_node);
};

#endif // EDGE_H
