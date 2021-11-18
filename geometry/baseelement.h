#ifndef BASEELEMENT_H
#define BASEELEMENT_H

#include <cstdint>


namespace icy { class BaseElement; class Element; }

class icy::BaseElement
{
public:
    enum ElementType {BEdge, TElem};
    ElementType type;

    virtual ~BaseElement() = default;

    virtual void ReplaceAdjacentElem(const Element* originalElem, Element* insertedElem, uint8_t idx) = 0;

    // if BaseElement is BoudnaryEdge, then UpdateNodes re-creates the list of nodes who which this element is conected
    virtual void UpdateNodes() {};
};

#endif // BASEELEMENT_H
