#include "edge.h"
#include "node.h"
#include "element.h"


icy::Edge::Edge(icy::Node* nd0, icy::Node* nd1)
{
    if(nd0->locId > nd1->locId) std::swap(nd0,nd1);
    nds[0] = nd0;
    nds[1] = nd1;
    elems[0] = elems[1] = nullptr;
}

void icy::Edge::AddElement(icy::Element* elem, short idx)
{
    Eigen::Vector3d u = nds[1]->x_initial - nds[0]->x_initial;
    icy::Node *opposite_node;
    ElementBoundaryFollowsEdge(elem, opposite_node);
    Eigen::Vector3d v1 = opposite_node->x_initial - nds[0]->x_initial;
    Eigen::Vector3d res2 = u.cross(v1);
    bool elem1_isCCW = res2.z() > 0;

    if(elem1_isCCW && elems[0] == nullptr) { elems[0] = elem; edge_in_elem_idx[0] = idx; }
    else if(!elem1_isCCW && elems[1] == nullptr) { elems[1] = elem; edge_in_elem_idx[1] = idx; }
    else throw std::runtime_error("mesh topology error");
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

/*
double icy::Edge::getAngle(icy::Node *center_node) const
{
    if(center_node == nds[0]) return angle0_initial;
    else if(center_node == nds[1]) return angle1_initial;
    else throw std::runtime_error("node does not belong to the edge");
}
*/

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
    else
    {
        qDebug() << "center node " << center_node->locId;
        qDebug() << "edge " << nds[0]->locId << " - " << nds[1]->locId;
        throw std::runtime_error("center node does not belong to the edge");
    }
}

icy::Element* icy::Edge::getTheOnlyElement()
{
    return elems[0]==nullptr ? elems[1] : elems[0];
}

icy::Element* icy::Edge::getOtherElement(const icy::Element* elem) const
{
    if (elems[0]==elem) return elems[1];
    else if(elems[1]==elem) return elems[0];
    else throw std::runtime_error("getOtherElement can't find other element");
}


icy::Element* icy::Edge::getElementWithNode(icy::Node *nd)
{
    if(elems[0] != nullptr && elems[0]->ContainsNode(nd)) return elems[0];
    else if(elems[1] != nullptr && elems[1]->ContainsNode(nd)) return elems[1];
    else {
        qDebug() << "edge " << nds[0]->locId << "-" << nds[1]->locId;
        qDebug() << "boundary " << isBoundary;
        qDebug() << "searching for element with node " << nd->locId;
        if(elems[0] == nullptr) qDebug() << "elem0 is null";
        else qDebug() << "elem0: " << elems[0]->nds[0]->locId << ", " << elems[0]->nds[1]->locId << ", " << elems[0]->nds[2]->locId;
        if(elems[1] == nullptr) qDebug() << "elem1 is null";
        else qDebug() << "elem1: " << elems[1]->nds[0]->locId << ", " << elems[1]->nds[1]->locId << ", " << elems[1]->nds[2]->locId;
        nd->PrintoutFan();
        throw std::runtime_error("getElementWithNode: cannot find element with a given node");
    }
}

bool icy::Edge::sameAs(Edge other)
{
    return (nds[0]==other.nds[0] && nds[1]==other.nds[1]);
}

bool icy::Edge::containsNode(Node *nd)
{
    return (nd==nds[0] || nd==nds[1]);
}
