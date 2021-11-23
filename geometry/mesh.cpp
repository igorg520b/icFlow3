#include "mesh.h"
#include "model.h"
#include <numeric>
#include <algorithm>
#include <iterator>
#include <spdlog/spdlog.h>

icy::SimpleObjectPool<icy::Node> icy::Mesh::s_pool_nodes(icy::Mesh::reserve_param);
icy::SimpleObjectPool<icy::Element> icy::Mesh::s_pool_elems(icy::Mesh::reserve_param*2);
icy::SimpleObjectPool<icy::BoundaryEdge> icy::Mesh::s_pool_edges(icy::Mesh::reserve_param/10);


icy::Mesh::Mesh()
{
    nodes.reserve(reserve_param);
    elems.reserve(reserve_param*2);
    edges.reserve(reserve_param/10);
}

icy::Mesh::~Mesh()
{
    Reset();
}

void icy::Mesh::Reset(unsigned typeOfSetup_)
{
    s_pool_nodes.release(nodes);
    s_pool_elems.release(elems);
    s_pool_edges.release(edges);

    length = width = area = 0;
}

/*
long icy::Geometry::ComputeFractureDirections(SimParams &prms, double timeStep, bool startingFracture)
{

    auto t1 = std::chrono::high_resolution_clock::now();

    maxNode=nullptr;
    double temporal_attenuation = prms.temporal_attenuation;

    float threashold = prms.normal_traction_threshold;

    std::size_t nNodes = nodes->size();


    if(startingFracture)
    {
        // evaluate all nodes to compute breakable range
        breakable_range_concurrent.clear();

        EvaluateStresses(prms, (*elems));
        DistributeStresses();

#pragma omp parallel for
        for(std::size_t i=0;i<nNodes;i++)
        {
            icy::Node *nd = nodes->at(i);
            if(nd->potentially_can_fracture)
            {
                if(nd->time_loaded_above_threshold >= temporal_attenuation)
                {
                    nd->ComputeFanVariables(prms);
                    if(nd->max_normal_traction > threashold) breakable_range_concurrent.push_back(nd);
                    else nd->time_loaded_above_threshold = 0;
                }
                else
                {
                    nd->time_loaded_above_threshold += timeStep;
                    nd->max_normal_traction = 0;
                    nd->dir.setZero();
                }
            }
            else
            {
                nd->time_loaded_above_threshold = 0;
                nd->max_normal_traction = 0;
                nd->dir.setZero();
            }
        }
        breakable_range.clear();
        std::copy(breakable_range_concurrent.begin(), breakable_range_concurrent.end(), std::back_inserter(breakable_range));
        std::sort(breakable_range.begin(), breakable_range.end(), [](Node *nd1, Node *nd2)
        {return nd1->max_normal_traction > nd2->max_normal_traction;});

        const unsigned max_breakable_range = 100;
        if(breakable_range.size() > max_breakable_range) breakable_range.resize(max_breakable_range);
    }
    else
    {

        // insert the recently created crack tips into the breakable range
        for(Node *nct : new_crack_tips)
        {
//            nct->PrepareFan2();
            nct->ComputeFanVariablesAlt(prms);
            nct->timeLoadedAboveThreshold = temporal_attenuation;
            auto find_result = std::find(breakable_range.begin(), breakable_range.end(),nct);
            bool already_contains = find_result!=breakable_range.end();

            if(!already_contains)
                breakable_range.push_back(nct);
        }
        new_crack_tips.clear();

        // remove the nodes that were affected by the crack on the previous step
        breakable_range.erase(std::remove_if(breakable_range.begin(), breakable_range.end(),
                                          [temporal_attenuation](Node *nd)
                {return nd->max_normal_traction==0 || (nd->timeLoadedAboveThreshold < temporal_attenuation && !nd->crack_tip);}),
                breakable_range.end());

        // update Sector in case if topology changed around this node
        for(Node *nd : breakable_range)
        {
            //nd->PrepareFan2();
            nd->ComputeFanVariablesAlt(prms);
        }

    }

    if(breakable_range.size() > 0)
    {
        // take out maximal node from breakable_range
        auto it_nd = std::max_element(breakable_range.begin(), breakable_range.end(),
                                      [](Node *nd1, Node *nd2) {
                if(nd2->crack_tip && nd2->max_normal_traction>0 && !nd1->crack_tip) return true;
                return nd1->max_normal_traction < nd2->max_normal_traction; });

        if((*it_nd)->max_normal_traction > 0)
        {
            maxNode = *it_nd;

            // make sure that the Sector information is updated

            //maxNode->PrepareFan2();
            //maxNode->ComputeFanVariablesAlt(prms);
            maxNode->timeLoadedAboveThreshold = 0;

#ifdef QT_DEBUG
            std::cout << "\n\nselected node " << maxNode->locId << std::endl;
            std::cout << "breakable range " << breakable_range.size() << "\n";
            for(Node *nd : breakable_range)
                std::cout << nd->locId << "; " << nd->max_normal_traction << (nd->crack_tip ? " *" : "") << std::endl;
#endif
            breakable_range.erase(it_nd);

        }
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
}





//========================== ALTERNATIVE ALGORITHM FOR SPLITTING


// code from icy::Geometry class that changes less often

#include "geometry.h"
#include <bits/stdc++.h>
#include <algorithm>




void icy::Geometry::ResizeNodes(std::size_t newSize)
{
    std::size_t nNodes = nodes->size();
    if(newSize == 0)
    {
        nodes->clear();
        s_pool_nodes.releaseAll();
    }
    else if(newSize > nNodes)
    {
        do {
            icy::Node* newNode = s_pool_nodes.take();
            newNode->Reset();
            newNode->locId = nodes->size();
            nodes->push_back(newNode);
        } while(nodes->size() < newSize);
    }
    else if(newSize < nNodes)
    {
        do {
            s_pool_nodes.release(nodes->back());
            nodes->pop_back();
        } while(nodes->size() > newSize);
    }
}

void icy::Geometry::ResizeElems(std::size_t newSize)
{
    std::size_t nElems = elems->size();
    if(newSize == 0)
    {
        elems->clear();
        s_pool_elems.releaseAll();
    }
    else if(newSize > nElems)
    {
        do {
            elems->push_back(s_pool_elems.take());
        } while(elems->size() < newSize);
    }
    else if(newSize < nElems)
    {
        do {
            s_pool_elems.release(elems->back());
            elems->pop_back();
        } while(elems->size() > newSize);
    }
}

icy::Node* icy::Geometry::AddNode()
{
    icy::Node* result = s_pool_nodes.take();
    result->Reset();
    result->locId = nodes->size();
    nodes->push_back(result);
    return result;
}

icy::Element* icy::Geometry::AddElement()
{
    icy::Element* elem = s_pool_elems.take();
    elem->Reset();
    elems->push_back(elem);
    return elem;
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
        nd->x_initial << x, y, 0;
        nd->xn << x, y, 0, 0, 0;
        nd->xt = nd->xn;
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
        elem->ComputeInitialNormal();
        if(elem->normal_initial.z()<0)
        {
            for(int j=0;j<3;j++) elem->nds[j] = (*nodes)[mtags[nodeTagsInTris[i*3+(2-j)]]];
            elem->ComputeInitialNormal();
        }
        else if(elem->normal_initial.z() < 0) throw std::runtime_error("normals inconsistent");

        area += elem->area_initial;
        for(int j=0;j<3;j++) elem->nds[j]->normal_n+=elem->normal_initial;
    }

    for(unsigned i=0;i<nodes->size();i++) (*nodes)[i]->normal_n.normalize();

    CreateEdges2();
    IdentifyDisconnectedRegions();
}


void icy::Geometry::CreateEdges2()
{
    std::size_t nNodes = nodes->size();

#pragma omp parallel for
    for(std::size_t i=0;i<nNodes;i++)
    {
        icy::Node* nd = (*nodes)[i];
        //nd->adjacent_nodes.clear();
        nd->adjacent_elems.clear();
        nd->isBoundary = false;
    }

    // edges_map will hold all edges and their connected elements
    tbb::concurrent_unordered_map<uint64_t, Edge> edges_map2;

    area = 0;

    // associate edges with one or two adjacent elements

    std::size_t nElems = elems->size();
#pragma omp parallel for
    for(std::size_t k=0;k<nElems;k++)
    {
        icy::Element *elem = (*elems)[k];
        for(int i=0;i<3;i++)
        {
            elem->adj_elems[i] = nullptr;
            elem->nds[i]->adjacent_elems.push_back(elem);
            // process edges
            int nd0idx = elem->nds[i]->locId;
            int nd1idx = elem->nds[(i+1)%3]->locId;

            if(nd0idx > nd1idx) std::swap(nd0idx, nd1idx);
            uint64_t key = ((uint64_t)nd0idx << 32) | nd1idx;

            icy::Node *nd0 = (*nodes)[nd0idx];
            icy::Node *nd1 = (*nodes)[nd1idx];

            Edge edge(nd0, nd1);

            edges_map2.insert({key,edge});
        }
    }

#pragma omp parallel for
    for(std::size_t k=0;k<nElems;k++)
    {
        icy::Element *elem = (*elems)[k];
        for(int i=0;i<3;i++)
        {
            int nd0idx = elem->nds[i]->locId;
            int nd1idx = elem->nds[(i+1)%3]->locId;

            if(nd0idx > nd1idx) std::swap(nd0idx, nd1idx);
            uint64_t key = ((uint64_t)nd0idx << 32) | nd1idx;

            icy::Edge &existing_edge = edges_map2.at(key);
            existing_edge.AddElement(elem, (i+2)%3);
        }
    }

    std::vector<Edge> allEdges;
    allEdges.resize(edges_map2.size());
    boundaryEdges.clear();

    std::size_t count = 0;
    for(auto &kvpair : edges_map2)
    {
        icy::Edge &e = kvpair.second;
        e.isBoundary = (e.elems[0] == nullptr || e.elems[1] == nullptr);
        if(e.isBoundary) boundaryEdges.push_back(e);
        allEdges[count++] = e;
    }

#pragma omp parallel for
    for(std::size_t i=0;i<allEdges.size();i++)
    {
        icy::Edge &existing_edge = allEdges[i];
        icy::Element *elem_of_edge0 = existing_edge.elems[0];
        icy::Element *elem_of_edge1 = existing_edge.elems[1];
        short idx0 = existing_edge.edge_in_elem_idx[0];
        short idx1 = existing_edge.edge_in_elem_idx[1];

        if(elem_of_edge0 == nullptr && elem_of_edge1 == nullptr) throw std::runtime_error("disconnected edge?");

        if(elem_of_edge0 != nullptr) elem_of_edge0->edges[idx0] = existing_edge;
        if(elem_of_edge1 != nullptr) elem_of_edge1->edges[idx1] = existing_edge;

        if(!existing_edge.isBoundary)
        {
            elem_of_edge0->adj_elems[idx0] = elem_of_edge1;
            elem_of_edge1->adj_elems[idx1] = elem_of_edge0;
        }
    }

#pragma omp parallel for
    for(std::size_t i=0;i<nodes->size();i++) (*nodes)[i]->PrepareFan2();
}

void icy::Geometry::AssignLsIds()
{
    std::size_t nNodes = nodes->size();
    int count = 0;
    for(std::size_t i=0;i<nNodes;i++)
    {
        icy::Node *nd = (*nodes)[i];
        //if(nd->prescribed) nd->lsId = -1;
        nd->lsId = count++;
    }
}


void icy::Geometry::WriteToHD5(unsigned offset_nodes, unsigned offset_elems,
                               hid_t ds_nodes_handle, hid_t ds_elems_handle)
{

    std::size_t nNodes = nodes->size();

    unsigned long nodes_extent = offset_nodes+nNodes;

    // set extent
    hsize_t nodes_dims[2] = {nodes_extent,icy::Node::NumberOfSerializedFields};
    H5Dset_extent(ds_nodes_handle, nodes_dims);

    // write

    hid_t file_space_id = H5Dget_space(ds_nodes_handle);
    constexpr unsigned long buffer_rows = 100;

    // write in blocks sized buffer_rows
    unsigned long written_node_count = 0;
    do
    {
        double node_buffer[buffer_rows][icy::Node::NumberOfSerializedFields];
        unsigned long writing_now = std::min(buffer_rows, nNodes-written_node_count);
        hsize_t mem_dims[2] = {writing_now,icy::Node::NumberOfSerializedFields};
        hid_t mem_space_id = H5Screate_simple(2, mem_dims, mem_dims);


        for(unsigned k=0;k<writing_now;k++)
        {
            unsigned idx = written_node_count+k;
            icy::Node *nd = (*nodes)[idx];
            //node_buffer[k][0] = nd->prescribed ? 1.0 : 0.0;
            node_buffer[k][1] = nd->x_initial.x();
            node_buffer[k][2] = nd->x_initial.y();
            node_buffer[k][3] = nd->timeLoadedAboveThreshold;

            for(int j=0;j<5;j++)
            {
                node_buffer[k][4+j] = nd->un[j];
                node_buffer[k][4+5+j] = nd->vn[j];
                node_buffer[k][4+10+j] = nd->an[j];
            }
        }
        hsize_t offset_nds[2] = {offset_nodes+written_node_count,0};
        hsize_t count[2] = {writing_now, icy::Node::NumberOfSerializedFields};
        H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, offset_nds, NULL, count, NULL);
        H5Dwrite(ds_nodes_handle, H5T_NATIVE_DOUBLE, mem_space_id, file_space_id, H5P_DEFAULT,(void*)node_buffer);
        H5Sclose(mem_space_id);

        written_node_count+=writing_now;
    }while(written_node_count < nNodes);

    H5Sclose(file_space_id);

    std::size_t nElems = elems->size();
    unsigned long elems_extent = offset_elems+nElems;

    // set extent
    hsize_t elems_dims[2] = {elems_extent,3};
    H5Dset_extent(ds_elems_handle, elems_dims);

    // write
    hsize_t mem_dims2[2] = {1,3};
    hid_t mem_space_id2 = H5Screate_simple(2, mem_dims2, mem_dims2);
    hid_t file_space_id2 = H5Dget_space(ds_elems_handle);
    int elems_buffer[3];
    for(unsigned i=0;i<nElems;i++)
    {
        for(int j=0;j<3;j++) elems_buffer[j] = (*elems)[i]->nds[j]->locId;
        hsize_t offset_e[2] = {offset_elems+i,0};
        H5Sselect_hyperslab(file_space_id2, H5S_SELECT_SET, offset_e, NULL, mem_dims2, NULL);
        H5Dwrite(ds_elems_handle, H5T_NATIVE_INT, mem_space_id2, file_space_id2, H5P_DEFAULT,(void*)elems_buffer);
    }

    H5Sclose(mem_space_id2);
    H5Sclose(file_space_id2);
}

void icy::Geometry::RestoreFromHD5(unsigned offset_nodes, unsigned offset_elems,
                                   unsigned nNodes, unsigned nElems,
                                   hid_t ds_nodes_handle, hid_t ds_elems_handle)
{
    ResizeNodes(nNodes);
    ResizeElems(nElems);

    hid_t file_space_id = H5Dget_space(ds_nodes_handle);
    hsize_t count[2] = {1, icy::Node::NumberOfSerializedFields};
    hid_t mem_space_id = H5Screate_simple(2, count, NULL);

    double node_buffer[icy::Node::NumberOfSerializedFields];
    for(unsigned i=0;i<nNodes;i++)
    {
        hsize_t offset[2] = {offset_nodes+i,0};
        H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        H5Dread(ds_nodes_handle, H5T_NATIVE_DOUBLE, mem_space_id, file_space_id, H5P_DEFAULT, (void*)node_buffer);

        icy::Node *nd = (*nodes)[i];
        nd->Reset();

        //nd->prescribed = node_buffer[0]==0 ? false : true;
        nd->x_initial << node_buffer[1], node_buffer[2], 0;
        nd->timeLoadedAboveThreshold = node_buffer[3];

        for(int j=0;j<5;j++)
        {
            nd->un(j)=node_buffer[4+j];
            nd->vn(j)=node_buffer[4+5+j];
            nd->an(j)=node_buffer[4+10+j];
        }
        nd->xn = nd->un;
        nd->xn.x()+=nd->x_initial.x();
        nd->xn.y()+=nd->x_initial.y();

        nd->xt = nd->xn;
        nd->ut = nd->un;
        nd->vt = nd->vn;
        nd->at = nd->an;
    }
    H5Sclose(file_space_id);
    H5Sclose(mem_space_id);


    int elems_buffer[3];

    hid_t file_space_id2 = H5Dget_space(ds_elems_handle);
    hsize_t count2[2] = {1,3};
    hid_t mem_space_id2 = H5Screate_simple(2, count2, NULL);

    for(unsigned i=0;i<nElems;i++)
    {
        hsize_t offset2[2] = {offset_elems+i,0};
        H5Sselect_hyperslab(file_space_id2, H5S_SELECT_SET, offset2, NULL, count2, NULL);
        H5Dread(ds_elems_handle, H5T_NATIVE_INT, mem_space_id2,
                file_space_id2, H5P_DEFAULT, (void*)elems_buffer);

        icy::Element *elem = (*elems)[i];
        for(int j=0;j<3;j++) elem->nds[j] = (*nodes)[elems_buffer[j]];

        elem->ComputeInitialNormal();
        if(elem->normal_initial.z() <0 ) {
            qDebug() << "negative elem normal";
            throw std::runtime_error("RestoreFromHD5: negative normals ");
        }
        elem->ComputeNormal();
    }

    H5Sclose(file_space_id2);
    H5Sclose(mem_space_id2);

    // -
    for(std::size_t i=0;i<nElems;i++)
    {
        icy::Element *elem = (*elems)[i];
        for(int j=0;j<3;j++) elem->nds[j]->normal_n+=elem->normal_n;
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
    IdentifyDisconnectedRegions();
}

//===========================================================

void icy::Geometry::EvaluateStresses(SimParams &prms, std::vector<Element*> &elems_range)
{
#pragma omp parallel for
    for(std::size_t i=0;i<elems_range.size();i++)
        elems_range[i]->EvaluateStresses(prms, elasticityMatrix, D_mats);
}

void icy::Geometry::DistributeStresses()
{
    std::size_t nElems = elems->size();
    std::size_t nNodes = nodes->size();
#pragma omp parallel for
    for(std::size_t i=0;i<nNodes;i++)
    {
        icy::Node *nd = (*nodes)[i];
        for(int k=0;k<3;k++) nd->str_b[k] = nd->str_m[k] = nd->str_b_top[k] = nd->str_b_bottom[k] = 0;
        nd->str_s[0] = nd->str_s[1] = 0;
        nd->potentially_can_fracture = false;
    }

#pragma omp parallel for
    for(std::size_t i=0;i<nElems;i++) (*elems)[i]->DistributeStresses();
}


void icy::Geometry::EvaluateAllNormalTractions(SimParams &prms)
{
#pragma omp parallel for
    for(std::size_t i=0;i<nodes->size();i++) (*nodes)[i]->ComputeFanVariablesAlt(prms);
}

long icy::Geometry::IdentifyDisconnectedRegions()
{
    auto t1 = std::chrono::high_resolution_clock::now();
    regions.clear();
    for(icy::Element *e : *elems) e->traversal = 0;  // set to not-traversed

    unsigned short current_region = 0;
    std::vector<Element*> wave;
    wave.reserve(elems->size());
    area = 0;
    for(icy::Element *e : *elems)
    {
        if(e->traversal != 0) continue;

        wave.push_back(e);
        unsigned count_elems = 0;
        double region_area = 0;
        while(wave.size() > 0)
        {
            icy::Element *elem = wave.back();
            wave.pop_back();
            count_elems++;
            region_area += elem->area_initial;
            elem->traversal = 1;
            elem->region = current_region;
            for(int i=0;i<3;i++)
            {
                icy::Element *adj_e = elem->adj_elems[i];
                if(adj_e!= nullptr && adj_e->traversal==0) wave.push_back(adj_e);
            }
        }
        regions.push_back(std::make_tuple(current_region, region_area, count_elems));
        current_region++;

        area+=region_area;
    }

    // for testing
    //    std::cout << "printing regions:\n";
    //    for(std::tuple<unsigned, double, unsigned> &r : regions)
    //        std::cout << std::get<0>(r) << ": " << std::get<1>(r) << "; " << std::get<2>(r) <<  std::endl;
    //    std::cout << "============= \n";
    // std::cout << "Regions " << regions.size() << std::endl;

    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
}


    long icy::Geometry::InferLocalSupport(SimParams &prms)
{
    auto t1 = std::chrono::high_resolution_clock::now();
    if(maxNode==nullptr) throw std::runtime_error("CreateSupportRange nullptr");
    local_elems.clear();
    std::copy(maxNode->adjacent_elems.begin(),maxNode->adjacent_elems.end(),std::back_inserter(local_elems));
    CreateSupportRange(prms.substep_radius, local_elems);

    std::unordered_set<Node*> local_support_set;
    for(Element *elem : local_elems) for(int k=0;k<3;k++) local_support_set.insert(elem->nds[k]);
    local_support.clear();
    std::copy(local_support_set.begin(), local_support_set.end(),std::back_inserter(local_support));

    // reset the loading timer in the vicinity of the crack
    local_elems2.clear();
    std::copy(maxNode->adjacent_elems.begin(),maxNode->adjacent_elems.end(),std::back_inserter(local_elems2));
    CreateSupportRange(prms.substep_radius2, local_elems2);
    for(icy::Node *nd : *nodes) nd->reset_timing = nd->support_node = false;
    for(Element *e : local_elems2) for(int k=0;k<3;k++)
        {
            e->nds[k]->timeLoadedAboveThreshold=0;
            e->nds[k]->reset_timing=true;
        }

    // for visualization - mark support range (stored in breakable_range)
    for(icy::Node *nd : local_support) nd->support_node = true; // for visualization


    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
}


void icy::Geometry::CreateSupportRange(int neighborLevel, std::vector<Element*> &initial_set)
{
#pragma omp parallel for
    for(unsigned i=0;i<elems->size();i++) (*elems)[i]->traversal=0;

    std::queue<Element*> q_wave;
    for(Element *e : initial_set)
    {
        e->traversal=1;
        q_wave.push(e);
    }
    initial_set.clear();

    while(q_wave.size() > 0)
    {
        icy::Element *elem = q_wave.front();
        q_wave.pop();
        initial_set.push_back(elem);

        unsigned short level = elem->traversal;
        if(level < neighborLevel)
        {
            for(int i=0;i<3;i++)
            {
                icy::Element *adj_e = elem->adj_elems[i];
                if(adj_e!= nullptr && adj_e->traversal==0)
                {
                    adj_e->traversal=level+1;
                    q_wave.push(adj_e);
                }
            }
        }
    }
}
*/



/*
void icy::GeneratorTool::GenerateBeam(BeamParams *beamParams, Mesh *outMesh)
{
    outMesh->isDeformable = true;
    double CharacteristicLengthMax = beamParams->CharacteristicLengthMax;
    double rm = beamParams->RefinementMultiplier;
    double a = beamParams->beamA; // beamA
    double b = beamParams->beamB; // beamB
    double l1 = beamParams->beamL1;
    double l2 = beamParams->beamL2;
    double c = beamParams->beamGap; // beam gap
    double d = beamParams->beamMargin; // beam margin
    double h = beamParams->beamThickness; // thickness

gmsh::clear();
gmsh::option::setNumber("General.Terminal", 1);
model::add("beam1");

double dy = l1 + c + d;
double dx = c+d;

int point1 = factory::addPoint(dx + 0, dy + 0, 0, rm);
int point2 = factory::addPoint(dx + l2, dy + 0, 0, rm);
int point22 = factory::addPoint(dx + l2 / 2,dy + 0, 0, rm);
int point3 = factory::addPoint(dx + l2, dy - a, 0, rm);
int point4 = factory::addPoint(dx + 0,dy - l1 + c / 2, 0, rm);
int point5 = factory::addPoint(dx - c,dy - l1 + c / 2, 0, rm);
int point6 = factory::addPoint(dx - c / 2,dy - l1 + c / 2, 0, rm);
int point7 = factory::addPoint(dx - c, dy + c, 0, 1.0);
int point8 = factory::addPoint(dx + l2 + c, dy + c, 0, 1.0);
int point9 = factory::addPoint(dx + l2 + c, dy - a - c, 0, 1.0);
int point10 = factory::addPoint(dx + b + 2 * c,dy - a, 0, rm);
int point11 = factory::addPoint(dx + b,dy - a - 2 * c, 0, rm);
int point12 = factory::addPoint(dx + b + 2 * c,dy - a - 2 * c,0, rm);
int point13 = factory::addPoint(dx + b + 2 * c,dy - a - c, 0, rm);
int point14 = factory::addPoint(dx + b + c,dy - a - 2 * c, 0, rm);
int point15 = factory::addPoint(dx + b,dy - l1 + c / 2, 0, rm);
int point16 = factory::addPoint(dx + b + c,dy - l1 + c / 2, 0, rm);
int point17 = factory::addPoint(dx + b + c / 2,dy - l1 + c / 2, 0, rm);
int point18 = factory::addPoint(-d, dy + c + d, 0, 1.0);
int point19 = factory::addPoint(dx + l2 + c + d, dy + c + d, 0, 1.0);
int point20 = factory::addPoint(dx + l2 + c + d, dy - l1 - c - 2*d, 0, 1.0);
int point21 = factory::addPoint(-d, -d, 0, rm);
int point23 = factory::addPoint(dx/4 + l2/4 + c/4 + d/4, -d , 0, rm);
int point24 = factory::addPoint(-d, dy/4 + c/4 + d/4, 0, rm);
int point25 = factory::addPoint(-d , dy/2 + c/2 + d/2, 0, 1.0);
int point26 = factory::addPoint(dx/2 + l2/2 + c/2 + d/2, -d, 0, 1.0);

int circle1 = factory::addCircleArc(point4, point6, point5);
int circle8 = factory::addCircleArc(point10, point12, point11);
int circle11 = factory::addCircleArc(point13, point12, point14);
int circle14 = factory::addCircleArc(point16, point17, point15);

int line2 = factory::addLine(point1, point22);
int line19 = factory::addLine(point22, point2);
int line3 = factory::addLine(point2, point3);
int line4 = factory::addLine(point4, point1);
int line5 = factory::addLine(point5, point7);
int line6 = factory::addLine(point7, point8);
int line7 = factory::addLine(point8, point9);
int line9 = factory::addLine(point3, point10);
int line10 = factory::addLine(point9, point13);
int line12 = factory::addLine(point11, point15);
int line13 = factory::addLine(point14, point16);
int line15 = factory::addLine(point18, point19);
int line16 = factory::addLine(point19, point20);
int line17 = factory::addLine(point20, point26);
int line18 = factory::addLine(point26, point23);
int line20 = factory::addLine(point23, point21);
int line21 = factory::addLine(point21, point24);
int line22 = factory::addLine(point24, point25);
int line23 = factory::addLine(point25, point18);

std::vector<int> curveTags;
curveTags.push_back(line15);
curveTags.push_back(line16);
curveTags.push_back(line17);
curveTags.push_back(line18);
curveTags.push_back(line20);
curveTags.push_back(line21);
curveTags.push_back(line22);
curveTags.push_back(line23);
int loop2 = factory::addCurveLoop(curveTags);

curveTags.clear();
curveTags.push_back(line6);
curveTags.push_back(line7);
curveTags.push_back(line10);
curveTags.push_back(circle11);
curveTags.push_back(line13);
curveTags.push_back(circle14);
curveTags.push_back(-line12);
curveTags.push_back(-circle8);
curveTags.push_back(-line9);
curveTags.push_back(-line3);
curveTags.push_back(-line19);
curveTags.push_back(-line2);
curveTags.push_back(-line4);
curveTags.push_back(circle1);
curveTags.push_back(line5);
int loop3 = factory::addCurveLoop(curveTags);

std::vector<int> loops;
loops.push_back(loop2);
loops.push_back(loop3);
factory::addPlaneSurface(loops);

factory::synchronize();

gmsh::vectorpair vp;
gmsh::vectorpair vpOut;

model::getEntities(vp, 2);
factory::extrude(vp, 0, 0, h, vpOut);

factory::synchronize();
gmsh::option::setNumber("Mesh.CharacteristicLengthMax", CharacteristicLengthMax);
model::mesh::generate(3);

// process the result

std::vector<std::size_t> nodeTags3;
std::vector<double> nodeCoords3, parametricCoords3;
model::mesh::getNodes(nodeTags3, nodeCoords3, parametricCoords3, -1, -1, false, false);

// retrieve elements
std::vector<std::size_t> elementTags, nodeTagsInElems;
model::mesh::getElementsByType(4, elementTags, nodeTagsInElems);
outMesh->elems.resize(elementTags.size());

// get triangles
std::vector<std::size_t> faceTags, nodeTagsInFaces;
model::mesh::getElementsByType(2, faceTags, nodeTagsInFaces);
outMesh->faces.resize(faceTags.size());

// compile a set of node tags that are present in elements
std::unordered_set<int> usedNodes;
for(int i=0;i<(int)nodeTagsInElems.size();i++) usedNodes.insert(nodeTagsInElems[i]);
for(int i=0;i<(int)nodeTagsInFaces.size();i++) usedNodes.insert(nodeTagsInFaces[i]);

std::map<int,int> nodeTagMap; // "gmsh tag" -> "sequential tag"
int count = 0;
outMesh->nodes.resize(usedNodes.size());
for(int i=0;i<(int)nodeTags3.size();i++) {
    int tag = nodeTags3[i];
    if(usedNodes.find(tag)!=usedNodes.end() && nodeTagMap.find(tag)==nodeTagMap.end())
    {
        if(count >= (int)usedNodes.size()) throw std::runtime_error("generator tool error");
        outMesh->nodes[count].Initialize(nodeCoords3[i*3+0], nodeCoords3[i*3+1], nodeCoords3[i*3+2], count);
        nodeTagMap[tag]=count;
        count++;
    }
}

// elements & faces
for(int i=0;i<(int)elementTags.size();i++)
{
    Element *elem = &(outMesh->elems[i]);
    for(int j=0;j<4;j++) {
        int ndidx = nodeTagsInElems[i*4+j];
        int newIdx = nodeTagMap[ndidx];
        elem->vrts[j] = &(outMesh->nodes[newIdx]);
    }
}

for(int i=0;i<(int)faceTags.size();i++)
{
    Face &fc = outMesh->faces[i];
    for(int j=0;j<3;j++) {
        int ndidx = nodeTagsInFaces[i*3+j];
        int newIdx = nodeTagMap[ndidx];
        fc.vrts[j]=&(outMesh->nodes[newIdx]);
    }
}
outMesh->ComputeBoundingBox();

// region of cz insertion
std::vector<p2d> region;
region.push_back(std::make_pair(dx + l2*0.6, dy+c/2));
region.push_back(std::make_pair(dx + b + 2*c, dy -a-c/2));
region.push_back(std::make_pair(dx + b + c/2, dy - a - 2*c));
region.push_back(std::make_pair(dx + b + c/2, dy - l1 - c - d/2));
region.push_back(std::make_pair(dx - c/2, dy - l1 - c - d/2));
region.push_back(std::make_pair(dx - c/2, dy - l1 / 2));
region.push_back(std::make_pair(dx - c / 2, dy + c / 2));


// exterior point for the region
double maxX = -DBL_MAX;
double maxY = -DBL_MAX;
for(auto &pt : region) {
    if(pt.first > maxX) maxX=pt.first;
    if(pt.second > maxY) maxY=pt.second;
}
p2d exteriorPt = std::make_pair(maxX+1, maxY+1);

count = 2;
for(auto &elem : outMesh->elems) {
    double centerX = 0, centerY = 0;
    for(int i=0;i<4;i++) {
        centerX+=elem.vrts[i]->x0/4;
        centerY+=elem.vrts[i]->y0/4;
    }
    p2d elemCenter = std::make_pair(centerX, centerY);
    if(PointInsideLoop(elemCenter, region, exteriorPt)) elem.tag = count++;
    else elem.tag = 1;
}
}
*/
