// code from icy::Geometry class that changes less often

#include "geometry.h"
#include <bits/stdc++.h>
#include <algorithm>

icy::Geometry::Geometry()
{
    tmp_range0->reserve(1000);
    tmp_range1->reserve(1000);
    breakable_range.reserve(100);
    neighbors_of_crack_tip.reserve(100);
    local_support.reserve(100);

    const std::size_t expected_size = 16384;
    nodes->reserve(expected_size);
    nodes2->reserve(expected_size);
    elems->reserve(expected_size);
    elems2->reserve(expected_size);
    //edges->reserve(expected_size);

    regions.reserve(100);
    length = width = area = 0;
}

void icy::Geometry::Reset()
{
    qDebug() << "icy::Geometry::Reset()";
    ResizeNodes(0);
    ResizeElems(0);
    //ResizeEdges(0);
    SwapCurrentAndTmp();
    ResizeNodes(0);
    ResizeElems(0);

    length = width = area = 0;
}

void icy::Geometry::SwapCurrentAndTmp()
{
    std::swap(nodes, nodes2);
    std::swap(elems, elems2);
}

void icy::Geometry::ResizeNodes(std::size_t newSize)
{
    std::size_t nNodes = nodes->size();
    if(newSize > nNodes)
    {
        do {
            icy::Node* newNode = pool_nodes.construct();
            newNode->locId = nodes->size();
            nodes->push_back(newNode);
        } while(nodes->size() < newSize);
    }
    else if(newSize < nNodes)
    {
        do {
//            pool_nodes.destroy(nodes->back());
            pool_nodes.free(nodes->back());
            nodes->pop_back();
        } while(nodes->size() > newSize);
    }
}

void icy::Geometry::ResizeElems(std::size_t newSize)
{
    std::size_t nElems = elems->size();
    if(newSize > nElems)
    {
        do {
            elems->push_back(pool_elems.malloc());
        } while(elems->size() < newSize);
    }
    else if(newSize < nElems)
    {
        do {
            pool_elems.free(elems->back());
            elems->pop_back();
        } while(elems->size() > newSize);
    }
}

/*
void icy::Geometry::ResizeEdges(std::size_t newSize)
{
    std::size_t nEdges = edges->size();
    if(newSize > nEdges)
    {
        do {
            edges->push_back(pool_edges.malloc());
        } while(edges->size() < newSize);
    }
    else if(newSize < nEdges)
    {
        do {
            pool_edges.free(edges->back());
            edges->pop_back();
        } while(edges->size() > newSize);
    }
}
*/

icy::Node* icy::Geometry::AddNode(icy::Node *otherNd)
{
    icy::Node* result = pool_nodes.construct();
    result->locId = nodes->size();
    nodes->push_back(result);
    if(otherNd!=nullptr) result->InitializeFromAnother(otherNd);
    return result;
}

icy::Element* icy::Geometry::AddElement()
{
    icy::Element* result = pool_elems.malloc();
    elems->push_back(result);
    return result;
}

void icy::Geometry::ImportFloePatch(QString fileName, double CharacteristicLengthMax)
{
    qDebug() << "importing floe patch from STL file " << fileName;
    if(fileName.isEmpty()) return;

    Reset();
    gmsh::clear();
    gmsh::option::setNumber("General.Terminal", 0);
    gmsh::open(fileName.toStdString());

    MeshingStepTwo(CharacteristicLengthMax);

    qDebug() << "patch import successful";
    // for testing
//    for(std::size_t i=0;i<boundary.size();i++) nodes[boundary[i]].prescribed = true;
}

void icy::Geometry::Remesh(double CharacteristicLengthMax)
{
    gmsh::clear();
    gmsh::option::setNumber("General.Terminal", 0);
    gmsh::model::add("floe1");

    const int dim = 2;
    int entityTag = gmsh::model::addDiscreteEntity(dim);
    //std::cout << "addDiscreteEntity return entity tag " << entityTag << std::endl;

    std::vector<std::size_t> ndTags(nodes->size());
    std::iota(ndTags.begin(), ndTags.end(), 1);
    std::vector<double> ndCoord(nodes->size()*3);
    for(std::size_t i=0;i<nodes->size();i++)
    {
        ndCoord[i*3+0] = (*nodes)[i]->x_initial.x();
        ndCoord[i*3+1] = (*nodes)[i]->x_initial.y();
        ndCoord[i*3+2] = 0;
    }
    gmsh::model::mesh::addNodes(dim, entityTag, ndTags, ndCoord);

    std::vector<int> elementTypes;
    elementTypes.push_back(2); // 2 == triangle

    std::vector<std::vector<std::size_t>> elementTags(1);
    std::vector<std::size_t> &elementTagsType2 = elementTags[0];
    elementTagsType2.resize(elems->size());
    std::iota(elementTagsType2.begin(), elementTagsType2.end(), 1);

    std::vector<std::vector<std::size_t>> nodeTags(1);
    std::vector<std::size_t> &nodeTags0 = nodeTags[0];
    nodeTags0.resize(elems->size()*3);

    for(std::size_t i=0;i<elems->size();i++)
        for(int j=0;j<3;j++) nodeTags0[i*3+j] = (*elems)[i]->nds[j]->locId+1;
    gmsh::model::mesh::addElements(dim,entityTag, elementTypes, elementTags, nodeTags);

    MeshingStepTwo(CharacteristicLengthMax);
}

void icy::Geometry::MeshingStepTwo(double CharacteristicLengthMax)
{
    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", CharacteristicLengthMax);
    gmsh::option::setNumber("Mesh.CharacteristicLengthExtendFromBoundary",0);

    double angle_threshold = 0.01*M_PI/180.0;
    gmsh::model::mesh::classifySurfaces(10.000*M_PI/180.0, false, false, angle_threshold);
    gmsh::model::mesh::createGeometry();

    gmsh::model::mesh::generate(2);

    // get nodes
    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, parametricCoords;
    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, parametricCoords);

    // get elems
    std::vector<std::size_t> trisTags, nodeTagsInTris;
    gmsh::model::mesh::getElementsByType(2, trisTags, nodeTagsInTris);
    gmsh::clear();

    // nodeTags, nodeCoords, nodeTagsInTris => nodes, elems
    std::map<std::size_t, int> nodeTagsMap1; // nodeTag -> its sequential position in nodeTag
    for(std::size_t i=0;i<nodeTags.size();i++) nodeTagsMap1[nodeTags[i]] = i;

    std::unordered_set<std::size_t> tagSet; // only keep nodes from the tris
    for(std::size_t &tag : nodeTagsInTris) tagSet.insert(tag);

    int count = 0;
    std::map<std::size_t, int> mtags; // nodeTag -> sequential position
    ResizeNodes(tagSet.size());

    double xmax, xmin, ymax, ymin;
    xmax = ymax = -DBL_MAX;
    xmin = ymin = DBL_MAX;
    for(const std::size_t &tag : tagSet)
    {
        int idx1 = nodeTagsMap1[tag];
        double x = nodeCoords[idx1*3+0];
        double y = nodeCoords[idx1*3+1];
        if(xmax < x) xmax = x;
        if(ymax < y) ymax = y;
        if(xmin > x) xmin = x;
        if(ymin > y) ymin = y;

        icy::Node* nd = (*nodes)[count];
        nd->Reset();
        nd->x_initial << x, y, 0, 0, 0;
        nd->xt = nd->xn = nd->x_initial;
        nd->locId = count;
        mtags[tag] = count;
        count++;
    }

    length = xmax-xmin;
    width = ymax-ymin;

    area = 0;
    ResizeElems(nodeTagsInTris.size()/3);

    for(std::size_t i=0;i<nodeTagsInTris.size()/3;i++)
    {
        icy::Element *elem = (*elems)[i];
        for(int j=0;j<3;j++) elem->nds[j] = (*nodes)[mtags[nodeTagsInTris[i*3+j]]];
        elem->InitializePersistentVariables();
        area += elem->area_initial;
        for(int j=0;j<3;j++) elem->nds[j]->normal_n+=elem->normal_initial;
    }

    for(unsigned i=0;i<nodes->size();i++) (*nodes)[i]->normal_n.normalize();

    CreateEdges2();
}


long icy::Geometry::CreateEdges2()
{
    auto t1 = std::chrono::high_resolution_clock::now();

    std::size_t nNodes = nodes->size();

#pragma omp parallel for
    for(std::size_t i=0;i<nNodes;i++)
    {
        icy::Node* nd = (*nodes)[i];
        nd->adjacent_nodes.clear();
        nd->fan.clear();
        nd->isBoundary = false;
        nd->adjacent_edges_map.clear();
        nd->area = 0;
    }

    // edges_map will hold all edges and their connected elements
    edges_map2.clear();
    area = 0;

    // associate edges with one or two adjacent elements

    std::size_t nElems = elems->size();
#pragma omp parallel for
    for(std::size_t k=0;k<nElems;k++)
    {
        icy::Element *elem = (*elems)[k];
        for(int i=0;i<3;i++)
        {
            // process edges
            int nd0 = elem->nds[i]->locId;
            int nd1 = elem->nds[(i+1)%3]->locId;

            if(nd0 > nd1) std::swap(nd0, nd1);
            uint64_t key = ((uint64_t)nd0 << 32) | nd1;

            Edge edge;
            edge.elems[0] = elem;
            edge.elems[1] = nullptr;
            edge.nds[0] = (*nodes)[nd0];
            edge.nds[1] = (*nodes)[nd1];

            auto insertion_result = edges_map2.insert({key,edge});
            if(insertion_result.second == false)
            {
                icy::Edge &existing_edge = edges_map2.at(key);
                if(existing_edge.elems[1]!=nullptr) throw std::runtime_error("edge is attached to 3+ elems");
                existing_edge.elems[1] = elem;
            }
        }
    }

    // distribute edges to nodes
    for(auto &kvpair : edges_map2)
    {
        icy::Edge &edge = kvpair.second;
        edge.RepairElementOrder();

        uint64_t key = kvpair.first;
        int nd0idx = key >> 32;
        int nd1idx = (uint32_t)key;// & 0xffffffff;

        icy::Node* nd0 = (*nodes)[nd0idx];
        icy::Node* nd1 = (*nodes)[nd1idx];

        // nodes keep the list of adjacent edges
        nd0->adjacent_edges_map.insert({nd1idx, edge});
        nd1->adjacent_edges_map.insert({nd0idx, edge});

        // populate the list of adjacent nodes
        nd0->adjacent_nodes.push_back(nd1);
        nd1->adjacent_nodes.push_back(nd0);
    }

    for(icy::Element* elem : *elems)
    {
        area += elem->area_initial;
        for(int i=0;i<3;i++)
        {
            // process adjacent elements
            icy::Node *nd = elem->nds[i];
            nd->area += elem->area_initial/3;
            icy::Node::Sector f(elem, nd);
            nd->fan.push_back(f);
        }
    }

#pragma omp parallel for
    for(std::size_t i=0;i<nodes->size();i++) {
        icy::Node* nd = (*nodes)[i];
        nd->PrepareFan();
        if(nd->fan.size()==0) throw std::runtime_error("disconnected node");
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
}


void icy::Geometry::AssignLsIds()
{
    std::size_t nNodes = nodes->size();
    int count = 0;
    for(std::size_t i=0;i<nNodes;i++)
    {
        icy::Node *nd = (*nodes)[i];
        if(nd->prescribed) nd->lsId = -1;
        else nd->lsId = count++;
    }
}


void icy::Geometry::PrecomputePersistentVariables(SimParams &prms)
{
    // compute elasticityMatrix and D_mats
    RecomputeElasticityMatrix(prms);
    std::size_t nElems = elems->size();
#pragma omp parallel for
    for(std::size_t i=0;i<nElems;i++)
        (*elems)[i]->PrecomputeStiffnessMatrix(prms, elasticityMatrix, D_mats);
}

void icy::Geometry::RecomputeElasticityMatrix(SimParams &prms)
{
    elasticityMatrix = Eigen::Matrix3d::Zero();
    double k = prms.YoungsModulus / (1-prms.PoissonsRatio*prms.PoissonsRatio);
    elasticityMatrix(0,0) = elasticityMatrix(1,1) = k;
    elasticityMatrix(0,1) = elasticityMatrix(1,0) = k*prms.PoissonsRatio;
    elasticityMatrix(2,2) = prms.YoungsModulus/((1+prms.PoissonsRatio)*2.0);
    D_mats = Eigen::Matrix2d::Identity();
    D_mats *= ((5.0/6.0)*prms.YoungsModulus/((1+prms.PoissonsRatio)*2.0));
}

void icy::Geometry::WriteToSerializationBuffers()
{
    std::size_t nNodes = nodes->size();
    std::size_t nElems = elems->size();

    node_buffer.resize(nNodes*icy::Node::NumberOfSerializedFields);
    elems_buffer.resize(nElems*3);

#pragma omp parallel for
    for(std::size_t i=0;i<nNodes;i++)
    {
        icy::Node *nd = (*nodes)[i];
        std::size_t idx = i*icy::Node::NumberOfSerializedFields;
        node_buffer[idx+0] = nd->prescribed ? 1.0 : 0.0;
        node_buffer[idx+1] = nd->x_initial.x();
        node_buffer[idx+2] = nd->x_initial.y();
        node_buffer[idx+3] = nd->timeLoadedAboveThreshold;

        for(int j=0;j<5;j++)
        {
            node_buffer[idx+4+j] = nd->un[j];
            node_buffer[idx+4+5+j] = nd->vn[j];
            node_buffer[idx+4+10+j] = nd->an[j];
        }
    }

#pragma omp parallel for
    for(std::size_t i=0;i<nElems;i++)
        for(int j=0;j<3;j++) elems_buffer[i*3+j] = (*elems)[i]->nds[j]->locId;
}

void icy::Geometry::RestoreFromSerializationBuffers()
{
    std::size_t nNodes = node_buffer.size()/icy::Node::NumberOfSerializedFields;
    std::size_t nElems = elems_buffer.size()/3;

    ResizeNodes(nNodes);
    ResizeElems(nElems);

#pragma omp parallel for
    for(std::size_t i=0;i<nNodes;i++)
    {
        icy::Node *nd = (*nodes)[i];
        nd->Reset();
        std::size_t idx = i*icy::Node::NumberOfSerializedFields;

        nd->prescribed = node_buffer[idx+0]==0 ? false : true;
        nd->x_initial << node_buffer[idx+1], node_buffer[idx+2], 0, 0, 0;
        nd->timeLoadedAboveThreshold = node_buffer[idx+3];

        for(int j=0;j<5;j++)
        {
            nd->un(j)=node_buffer[idx+4+j];
            nd->vn(j)=node_buffer[idx+4+5+j];
            nd->an(j)=node_buffer[idx+4+10+j];
        }
        nd->xt = nd->xn = nd->x_initial + nd->un;
        nd->ut = nd->un;
        nd->vt = nd->vn;
        nd->at = nd->an;
    }

#pragma omp parallel for
    for(std::size_t i=0;i<nElems;i++)
    {
        icy::Element *elem = (*elems)[i];
        for(int j=0;j<3;j++) elem->nds[j] = (*nodes)[elems_buffer[i*3+j]];
        elem->InitializePersistentVariables();
        elem->ComputeNormal();
    }

    area=0;
    for(std::size_t i=0;i<nElems;i++)
    {
        icy::Element *elem = (*elems)[i];
        for(int j=0;j<3;j++) elem->nds[j]->normal_n+=elem->normal_n;
        area+= elem->area_initial;
    }

#pragma omp parallel for
    for(unsigned i=0;i<nodes->size();i++)
        (*nodes)[i]->normal_n.normalize();


    double xmax, xmin, ymax, ymin;
    xmax = ymax = -DBL_MAX;
    xmin = ymin = DBL_MAX;
    for(unsigned i=0;i<nodes->size();i++) {
        icy::Node *nd = (*nodes)[i];
        double x = nd->x_initial.x();
        double y = nd->x_initial.y();
        if(xmax < x) xmax = x;
        if(ymax < y) ymax = y;
        if(xmin > x) xmin = x;
        if(ymin > y) ymin = y;
    }
    length = xmax-xmin;
    width = ymax-ymin;

    CreateEdges2();
}

icy::Edge icy::Geometry::getEdgeByNodalIdx(int idx1, int idx2)
{
    if(idx1 > idx2) std::swap(idx1, idx2);
    uint64_t edgeIdx = ((uint64_t)idx1 << 32) | idx2;
    //return edges_map.at(edgeIdx).edge;
    return edges_map2.at(edgeIdx);
}
