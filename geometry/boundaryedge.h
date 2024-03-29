#ifndef BOUNDARYEDGE_H
#define BOUNDARYEDGE_H


#include "baseelement.h"
#include "element.h"

namespace icy {class BoundaryEdge; class Node;}


class icy::BoundaryEdge : public BaseElement
{
public:
    Node* nds[2];
    Element *elem;                      // element to which this boundary belongs
    uint8_t edge_idx;                   // side in the element to which this boundary is attached
    uint8_t status;                     // color/type for visualization

    BoundaryEdge() {type = ElementType::BEdge;};

    BoundaryEdge* Initialize(Element *elem_, uint8_t edge_idx_, uint8_t status_ = 0)
    {
        status = status_;
        elem = elem_;
        edge_idx = edge_idx_;
        UpdateNodes();
        return this;
    }

    BoundaryEdge* Initialize(Node *nd0, Node *nd1)
    {
        status = 1;
        elem = nullptr;
        nds[0] = nd0;
        nds[1] = nd1;
        return this;
    }

    void ReplaceAdjacentElem(const Element* originalElem, Element* insertedElem, uint8_t idx) override
    {
        Initialize(insertedElem, idx);
    }


    // infer nds[2] from elem and edge_idx
    void UpdateNodes()
    {
        if(elem!=nullptr)
        {
            uint8_t nd1Idx = (edge_idx+2)%3;
            uint8_t nd2Idx = (edge_idx+1)%3;
//            MeshFragment *fr = elem->nds[0]->fragment;
            nds[0] = elem->nds[nd1Idx];
            nds[1] = elem->nds[nd2Idx];
        }
        else throw std::runtime_error("BoundaryEdge::UpdateNodes called with nullptr element");
    }
};

#endif // BOUNDARYEDGE_H
