#if !defined(Q_MOC_RUN) // MOC has a glitch when parsing TBB headers

#ifndef FL333_H
#define FL333_H

#include <gmsh.h>

#include <vector>
#include <utility>
#include <algorithm>
#include <chrono>
#include <tuple>
#include <queue>
#include <set>

#include "element.h"
#include "linearsystem.h"
#include "edge.h"
#include "SimpleObjectPool.h"

#include <hdf5.h>


namespace icy { class Mesh; class Node; class Element; class Edge; class BoundaryEdge;}

class icy::Mesh
{
public:
    Mesh();
    ~Mesh();
    Mesh& operator=(Mesh&) = delete;

    // nodes, elements and edges
    std::vector<icy::Node*> nodes;
    std::vector<icy::Element*> elems;
    std::vector<icy::BoundaryEdge*> edges;

    double length, width, area;

    // SCENE
    void Reset(unsigned typeOfSetup_ = 0);
    void ImportFloePatch(std::string fileName, double ElemSizeMax);
//    void Remesh(double CharacteristicLengthMax);

private:
    void SetupLShapedBeam();
    constexpr static unsigned reserve_param = 10000;

    icy::Node* AddNode();
    icy::Element* AddElement();

    void ResizeNodes(std::size_t newSize);
    void ResizeElems(std::size_t newSize);


    // SERIALIZATION
    void WriteToHD5(unsigned offset_nodes, unsigned offset_elems,
                    hid_t ds_nodes_handle, hid_t ds_elems_handle);
    void RestoreFromHD5(unsigned offset_nodes, unsigned offset_elems,
                        unsigned nNodes, unsigned nElems,
                        hid_t ds_nodes_handle, hid_t ds_elems_handle);


    // FRACTURE
private:
    std::vector<Node*> breakable_range;     // populated in ComputeFractureDirections() when startingFracture==true
    std::vector<Node*> new_crack_tips;      // populated in SplitNode(), then used when startingFracture==false
    icy::Node *maxNode;
    constexpr static double fracture_epsilon = 0.1;   // if an edge splits too close to its vertex, then just go through the vertex
    void ComputeFractureDirections(const SimParams &prms, double timeStep, bool startingFracture);
    void PropagateCrack(const SimParams &prms);
    void EstablishSplittingEdge(Node* nd, const double phi, const double theta, Element *elem, Node* &adjacentNode);
    Node* Fix_X_Topology(Node *nd_to_split, Node *alignment_node);

    void InferLocalSupport(SimParams &prms);
    void ResetFractureTimer(SimParams &prms);
    void CreateSupportRange(const int neighborLevel);     // result is in local_elems, local_czs, local_support


    std::vector<Element*> local_elems; // elems corresponding to breakable_range;
    std::vector<Node*> local_support;

    long IdentifyDisconnectedRegions(); // used to deal with small fragments
    std::vector<std::tuple<unsigned, double, unsigned>> regions; // region#, area, element count

    /*

    void AssignLsIds();
    void CreateEdges2();     // from the list of elements, infer inner edges and boundary


    void EvaluateStresses(SimParams &prms, std::vector<Element*> &elems_range);     // needed for ComputeFractureDirections
    void DistributeStresses();                  // needed for visualization
    long ComputeFractureDirections(SimParams &prms, double timeStep = 0, bool startingFracture = false); // sets maxNode to breakable node
    long InferLocalSupport(SimParams &prms);
    void EvaluateAllNormalTractions(SimParams &prms); // for visualization when reloading from file


    // FRACTURE
    icy::Node *maxNode = nullptr;

*/











    void CreateSupportRange(int neighborLevel, std::vector<Element*> &initial_set);


    static icy::SimpleObjectPool<icy::Node> s_pool_nodes;
    static icy::SimpleObjectPool<icy::Element> s_pool_elems;
    static icy::SimpleObjectPool<icy::BoundaryEdge> s_pool_edges;


    friend class icy::Element;
    friend class icy::Model;
};

#endif
#endif // Q_MOC_RUN
