#if !defined(Q_MOC_RUN) // MOC has a glitch when parsing TBB headers

#ifndef FL333_H
#define FL333_H

#include <gmsh.h>

#include <vector>
#include <memory>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <set>
#include <utility>
#include <algorithm>
#include <chrono>
#include <tuple>
#include <queue>
#include <set>

#include <QObject>
#include <QString>
#include <QDebug>

#include "linearsystem.h"
#include "element.h"
#include "edge.h"
#include "SimpleObjectPool.h"

#include <hdf5.h>

#include <concurrent_unordered_map.h>

namespace icy { class Geometry; class Node; class Element; class Edge;}

class icy::Geometry
{
public:
    Geometry();

    // nodes, elements and edges
    std::unique_ptr<std::vector<icy::Node*>> nodes = std::make_unique<std::vector<icy::Node*>>();
    std::unique_ptr<std::vector<icy::Element*>> elems = std::make_unique<std::vector<icy::Element*>>();

    double length, width, area;

    // at the "setting up the scene" stage - remesh 2d floe if needed
    void Reset();
    void ImportFloePatch(QString fileName, double CharacteristicLengthMax);
    void Remesh(double CharacteristicLengthMax);

    void RecomputeElasticityMatrix(SimParams &prms);
    void AssignLsIds();
    void CreateEdges2();     // from the list of elements, infer inner edges and boundary
    long IdentifyDisconnectedRegions(); // used to deal with small fragments
//    long RemoveDegenerateFragments();
    std::vector<std::tuple<unsigned, double, unsigned>> regions; // region#, area, element count

    unsigned getElemCount() {return elems->size();}
    unsigned getNodeCount() {return nodes->size();}
    unsigned getFreeNodeCount() { return nodes->size(); }
//        return std::count_if(nodes->begin(), nodes->end(), [](Node* &nd){return !nd->prescribed;}); }

    void EvaluateStresses(SimParams &prms, std::vector<Element*> &elems_range);     // needed for ComputeFractureDirections
    void DistributeStresses();                  // needed for visualization
    long ComputeFractureDirections(SimParams &prms, double timeStep = 0, bool startingFracture = false); // sets maxNode to breakable node
    long InferLocalSupport(SimParams &prms);
    void EvaluateAllNormalTractions(SimParams &prms); // for visualization when reloading from file

    long SplitNodeAlt(SimParams &prms);
    std::vector<Node*> new_crack_tips; // populated by SplitNodeAlt

    // save/load
    void WriteToHD5(unsigned offset_nodes, unsigned offset_elems,
                    hid_t ds_nodes_handle, hid_t ds_elems_handle);
    void RestoreFromHD5(unsigned offset_nodes, unsigned offset_elems,
                        unsigned nNodes, unsigned nElems,
                        hid_t ds_nodes_handle, hid_t ds_elems_handle);

    // fracture
    tbb::concurrent_vector<Node*> breakable_range_concurrent;
    std::vector<Node*> breakable_range;
    std::vector<Node*> neighbors_of_crack_tip, local_support;
    std::vector<Element*> local_elems; // elems corresponding to breakable_range;
    //std::unordered_set<Element*> local_elems_set; // for computing local_elems
    icy::Node *maxNode = nullptr;

    std::vector<Edge> boundaryEdges;

    // matrices for elements, recomputed when rho/Y/nu change
    Eigen::Matrix3d elasticityMatrix;        // this has to be pre-computed whenever Y and nu change
    Eigen::Matrix2d D_mats;

private:
    void ResizeNodes(std::size_t newSize);
    void ResizeElems(std::size_t newSize);

    icy::Node* AddNode(icy::Node *otherNd=nullptr);
    icy::Element* AddElement(); // makes a new element

    std::unordered_set<Element*> affected_elements_during_split; // a list of elements that were affected by SplitNode
    void UpdateEdges();

    void EstablishSplittingEdge(Edge &splitEdge, Node* nd,
                                const float phi, const float theta, const float fracture_epsilon,
                                const Edge e0, const Edge e1, Element *elem);
    void Fix_X_Topology(Node *nd);
    // preserve boundaries and orientation
    void CarefulSplitBoundaryElem(Element *originalElem, Node *nd, Node *nd0, Node *nd1, float where, Edge &insertedEdge);
    void CarefulSplitNonBoundaryElem(Element *originalElem, Element *adjElem, Node *nd,
                                     Node *nd0, Node *nd1, float where, Edge &insertedEdge);

    void MeshingStepTwo(double CharacteristicLengthMax);

    icy::SimpleObjectPool<Node> s_pool_nodes;
    icy::SimpleObjectPool<Element> s_pool_elems;

    void CreateSupportRange(int neighborLevel, std::vector<Element*> &initial_set);
    std::vector<Element*>local_elems2;

//    void RemoveRegion(unsigned idx);

};
#endif
#endif // Q_MOC_RUN
