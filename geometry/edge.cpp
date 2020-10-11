#include "edge.h"
#include "node.h"
#include "element.h"

void icy::Edge::Initialize(icy::Node* nd0, icy::Node* nd1, icy::Element* elem0, icy::Element* elem1)
{
    elems[0] = elems[1] = nullptr;
    nds[0] = nd0;
    nds[1] = nd1;

    // calculate edge angle
    Eigen::Matrix<double,DOFS,1> e = nds[1]->x_initial - nds[0]->x_initial;
    angle0_initial = atan2(e.y(), e.x());
    angle1_initial = (angle0_initial > 0) ? angle0_initial - M_PI : angle0_initial + M_PI;

    Eigen::Vector3d u = e.block(0,0,3,1);
    // get oppisite node for elem0
    icy::Node *opposite_node0;
    ElementBoundaryFollowsEdge(elem0, opposite_node0);
    Eigen::Vector3d v0 = (opposite_node0->x_initial - nds[0]->x_initial).block(0,0,3,1);
    res1 = u.cross(v0);
    elem0_isCCW = res1.z() > 0;
    elems[elem0_isCCW ? 0 : 1] = elem0;

    isBoundary = (elem1==nullptr);
    if(!isBoundary)
    {
        icy::Node *opposite_node1;
        ElementBoundaryFollowsEdge(elem1, opposite_node1);
        Eigen::Vector3d v1 = (opposite_node1->x_initial - nds[0]->x_initial).block(0,0,3,1);
        res2 = u.cross(v1);
        elem1_isCCW = res2.z() > 0;
        if(elem1_isCCW == elem0_isCCW) throw std::runtime_error("two elems on the same side of edge");
        elems[elem1_isCCW ? 0 : 1] = elem1;
    } else
    {
        nd0->isBoundary = nd1->isBoundary = true;
    }
}

// determine if nds[0],nds[1] are contained in elem.nds in forward order
bool icy::Edge::ElementBoundaryFollowsEdge(icy::Element* elem, icy::Node* &opposite_node)
{
    icy::Node* n0 = nds[0];
    icy::Node* n1 = nds[1];
    if(elem->nds[0] == n0 && elem->nds[1] == n1) {opposite_node = elem->nds[2]; return true;}
    if(elem->nds[1] == n0 && elem->nds[0] == n1) {opposite_node = elem->nds[2]; return false;}

    if(elem->nds[1] == n0 && elem->nds[2] == n1) {opposite_node = elem->nds[0]; return true;}
    if(elem->nds[2] == n0 && elem->nds[1] == n1) {opposite_node = elem->nds[0]; return false;}

    if(elem->nds[2] == n0 && elem->nds[0] == n1) {opposite_node = elem->nds[1]; return true;}
    if(elem->nds[0] == n0 && elem->nds[2] == n1) {opposite_node = elem->nds[1]; return false;}

    throw std::runtime_error("element does not contain the edge");
}


double icy::Edge::getAngle(icy::Node *center_node) const
{
    if(center_node == nds[0]) return angle0_initial;
    else if(center_node == nds[1]) return angle1_initial;
    else throw std::runtime_error("node does not belong to the edge");
}

icy::Element* icy::Edge::get_CCW_Element(icy::Node *center_node) const
{
    if(center_node == nds[0]) return elems[0];
    else if(center_node == nds[1]) return elems[1];
    else throw std::runtime_error("node does not belong to the edge");
}

Eigen::Vector2f icy::Edge::getVec(icy::Node *center_node) const
{
    Node* other = getOtherNode(center_node);
    float x = (float)(other->xt.x()-center_node->xt.x());
    float y = (float)(other->xt.y()-center_node->xt.y());
    return Eigen::Vector2f(x,y);
}

icy::Node* icy::Edge::getOtherNode(icy::Node *center_node) const
{
    if(center_node == nds[0]) return nds[1];
    else if(center_node == nds[1]) return nds[0];
    else throw std::runtime_error("center node does not belong to the edge");
}

icy::Element* icy::Edge::getTheOnlyElement()
{
    return elems[0]==nullptr ? elems[1] : elems[0];
}

icy::Element* icy::Edge::getElementWithNode(icy::Node *nd)
{
    if(elems[0] != nullptr && elems[0]->ContainsNode(nd)) return elems[0];
    else if(elems[1] != nullptr && elems[1]->ContainsNode(nd)) return elems[1];
    else throw std::runtime_error("cannot find element with a given node");
}

icy::Node* icy::Edge::getFarNode(icy::Node *nd)
{
    if(elems[0]->ContainsNode(nd)) return elems[1]->getOppositeNode(this);
    else if(elems[1]->ContainsNode(nd)) return elems[0]->getOppositeNode(this);
    else throw std::runtime_error("cannot find element with a given node");
}
