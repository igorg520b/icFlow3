#ifndef EDGE_H
#define EDGE_H

#include <Eigen/Core>

namespace icy {class Edge; class Node; class Element;}

class icy::Edge
{
public:
    Edge() {}
    Edge(icy::Node* nd0, icy::Node* nd1);
    void AddElement(icy::Element* elem, short idx);

    icy::Node* nds[2];
    bool isBoundary;    // belongs to only one element
    icy::Element* elems[2];
    short edge_in_elem_idx[2];  // index [0,3] where the edge is located in the corresponding element 0 and 1
//    double angle0_initial; // [-pi, +pi] from nds[0] to nds[1]
//    double angle1_initial; // [-pi, +pi] from nds[1] to nds[0]
    bool toSplit = false; // for mesh splitting

//    double getAngle(icy::Node *center_node) const;
    icy::Element* get_CCW_Element(icy::Node *center_node) const; // adjacent element from ccw side; deprecated
    Eigen::Vector2f getVec(icy::Node *center_node) const;     // as vector at step n
    icy::Node* getOtherNode(icy::Node *center_node) const;
    icy::Element* getTheOnlyElement();
    icy::Element* getOtherElement(const icy::Element* elem) const;

    icy::Element* getElementWithNode(icy::Node *nd);
    icy::Node* getFarNode(icy::Node *nd);
    bool sameAs(Edge other);
    bool containsNode(Node *nd);

private:
    bool ElementBoundaryFollowsEdge(icy::Element* elem, icy::Node* &opposite_node);
};

#endif // EDGE_H
