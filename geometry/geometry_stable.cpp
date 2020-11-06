// code from icy::Geometry class that changes less often

#include "geometry.h"
#include <bits/stdc++.h>
#include <algorithm>

icy::Geometry::Geometry()
{
    breakable_range.reserve(100);
    neighbors_of_crack_tip.reserve(100);
    local_support.reserve(100);

    const std::size_t expected_size = 16384;
    nodes->reserve(expected_size);
    elems->reserve(expected_size);

    regions.reserve(100);
    length = width = area = 0;

    s_pool_nodes.reserve(5000);
    s_pool_elems.reserve(5000);
}

void icy::Geometry::Reset()
{
    qDebug() << "icy::Geometry::Reset()";
    ResizeNodes(0);
    ResizeElems(0);
    length = width = area = 0;
}


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

icy::Node* icy::Geometry::AddNode(icy::Node *otherNd)
{
    icy::Node* result = s_pool_nodes.take();
    result->Reset();
    result->locId = nodes->size();
    nodes->push_back(result);
    if(otherNd!=nullptr) result->InitializeFromAnother(otherNd);
    return result;
}

icy::Element* icy::Geometry::AddElement()
{
    icy::Element* result = s_pool_elems.take();
    elems->push_back(result);
    for(int i=0;i<3;i++) result->adj_elems[i]=nullptr;
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
/*
long icy::Geometry::RemoveDegenerateFragments()
{
    auto t1 = std::chrono::high_resolution_clock::now();
    double avg_elem_area = area/elems->size();
    bool region_was_removed;
    do
    {
        region_was_removed = false;
        auto iter = std::min_element(regions.begin(), regions.end(),
                                      [](std::tuple<unsigned, double, unsigned> r1,
                                      std::tuple<unsigned, double, unsigned> r2)
        {return std::get<2>(r1) < std::get<2>(r2);});
        if(std::get<2>(*iter) <= 2)
        {
            unsigned idx = std::get<0>(*iter);
            RemoveRegion(idx);
            regions.erase(iter);
            region_was_removed = true;
        }

        iter = std::min_element(regions.begin(), regions.end(),
                                [](std::tuple<unsigned, double, unsigned> r1,
                                std::tuple<unsigned, double, unsigned> r2)
        {return std::get<1>(r1) < std::get<1>(r2);});

        if(std::get<1>(*iter) < avg_elem_area*1.5)
        {
            unsigned idx = std::get<0>(*iter);
            RemoveRegion(idx);
            regions.erase(iter);
            region_was_removed = true;
        }

    } while(region_was_removed && regions.size() > 0);

    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
}

void icy::Geometry::RemoveRegion(unsigned idx)
{
    std::unordered_set<icy::Node*>nds;
    for(icy::Element *elem : *elems)
    {
        if(elem->region == idx) {
            for(int k=0;k<3;k++) nds.insert(elem->nds[k]);
            pool_elems.free(elem);
        }
    }

    elems->erase(std::remove_if(elems->begin(), elems->end(),
                                [idx](icy::Element *elem){return elem->region==idx;}), elems->end());

    for(icy::Node *nd : nds) pool_nodes.destroy(nd);

    nodes->erase(std::remove_if(nodes->begin(), nodes->end(),
                                [nds](icy::Node *nd){return nds.find(nd)!=nds.end();}), nodes->end());

    for(std::size_t i=0;i<nodes->size();i++) (*nodes)[i]->locId=i;

    CreateEdges2();
}
*/
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
