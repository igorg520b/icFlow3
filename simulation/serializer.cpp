#include "serializer.h"
#include <QDebug>
//#include <iostream>
#include "frameinfo.h"


void icy::Serializer::CreateFile(std::string fileName, unsigned params_size)
{
    if(fileIsOpen) CloseFile();
    qDebug() << "icy::Serializer::CreateFile";

    file_handle = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    const int deflate_val = 5;
    // nodes
    hsize_t current_dims_nodes[2] = {0,NodeDataFields}; // can this be replaced by NULL ?
    hsize_t max_dims_nodes[2] = {H5S_UNLIMITED,NodeDataFields};
    hid_t space_nodes = H5Screate_simple(2, current_dims_nodes, max_dims_nodes);
    hid_t dcpl_nodes = H5Pcreate (H5P_DATASET_CREATE);
    H5Pset_deflate (dcpl_nodes, deflate_val);
    hsize_t chunk_nodes[2] = {2000,NodeDataFields};
    H5Pset_chunk(dcpl_nodes, 2, chunk_nodes);
    ds_nodes_handle = H5Dcreate(file_handle, "nodes", H5T_NATIVE_DOUBLE, space_nodes,
                                H5P_DEFAULT, dcpl_nodes, H5P_DEFAULT);
    H5Pclose(dcpl_nodes);
    H5Sclose(space_nodes);


    // elems
    hsize_t current_dims_elems[2] = {0,3}; // can this be replaced by NULL ?
    hsize_t max_dims_elems[2] = {H5S_UNLIMITED,3};
    hid_t space_elems = H5Screate_simple(2, current_dims_elems, max_dims_elems);
    hid_t dcpl_elems = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_deflate(dcpl_elems, deflate_val);
    hsize_t chunk_elems[2] = {2000,3};
    H5Pset_chunk(dcpl_elems, 2, chunk_elems);
    ds_elems_handle = H5Dcreate(file_handle, "elems", H5T_NATIVE_INT, space_elems,
                                H5P_DEFAULT, dcpl_elems, H5P_DEFAULT);
    H5Pclose(dcpl_elems);
    H5Sclose(space_elems);


    // steps
    hsize_t current_dims_steps[2] = {0,sizeof(icy::FrameInfo)}; // can this be replaced by NULL ?
    hsize_t max_dims_steps[2] = {H5S_UNLIMITED,sizeof(icy::FrameInfo)};
    hid_t space_steps = H5Screate_simple(2, current_dims_steps, max_dims_steps);
    hid_t dcpl_steps = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_deflate (dcpl_steps, deflate_val);
    hsize_t chunk_steps[2] = {200,sizeof(icy::FrameInfo)};
    H5Pset_chunk(dcpl_steps, 2, chunk_steps);
    ds_steps_handle = H5Dcreate(file_handle, "steps", H5T_NATIVE_CHAR, space_steps,
                                H5P_DEFAULT, dcpl_steps, H5P_DEFAULT);
    H5Pclose(dcpl_steps);
    H5Sclose(space_steps);

    // params
    hsize_t current_dims_params[1] = {params_size}; // can this be replaced by NULL ?
    hsize_t max_dims_params[1] = {params_size};
    hid_t space_params = H5Screate_simple(1, current_dims_params, max_dims_params);
    hid_t dcpl_params = H5Pcreate(H5P_DATASET_CREATE);
    ds_params_handle = H5Dcreate(file_handle, "params", H5T_NATIVE_CHAR, space_params,
                                H5P_DEFAULT, dcpl_params, H5P_DEFAULT);
    H5Pclose(dcpl_params);
    H5Sclose(space_params);

    fileIsOpen = true;
    qDebug() << "finished icy::Serializer::CreateFile";
}

void icy::Serializer::OpenFile(std::string fileName)
{
    if(fileIsOpen) CloseFile();

    file_handle = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    ds_nodes_handle = H5Dopen(file_handle, "nodes", H5P_DEFAULT);
    ds_elems_handle = H5Dopen(file_handle, "elems", H5P_DEFAULT);
    ds_steps_handle = H5Dopen(file_handle, "steps", H5P_DEFAULT);
    ds_params_handle = H5Dopen(file_handle, "params", H5P_DEFAULT);

    fileIsOpen = true;
}

void icy::Serializer::CloseFile()
{
    if(!fileIsOpen) return;
    qDebug() << "icy::Serializer::CloseFile()";
    H5Dclose(ds_nodes_handle);
    H5Dclose(ds_elems_handle);
    H5Dclose(ds_steps_handle);
    H5Dclose(ds_params_handle);
    H5Fclose(file_handle);
    fileIsOpen = false;
}

void icy::Serializer::LoadParams(void *data, unsigned length)
{
    if(!fileIsOpen) return;

    hid_t file_space_id = H5Dget_space(ds_params_handle);
    hsize_t offset[1] = {0};
    hsize_t count[1] = {length};
    H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    hid_t mem_space_id = H5Screate_simple(1, count, NULL);
    H5Dread(ds_params_handle, H5T_NATIVE_CHAR, mem_space_id, file_space_id,
            H5P_DEFAULT, (void*)data);
    H5Sclose(file_space_id);
    H5Sclose(mem_space_id);

}

void icy::Serializer::SaveParams(void *data, unsigned length)
{
    if(!fileIsOpen) return;

    hsize_t current_dims[1] = {length};
    hid_t mem_space_id = H5Screate_simple(1, current_dims, current_dims);
    H5Dwrite(ds_params_handle, H5T_NATIVE_CHAR, mem_space_id, H5S_ALL, H5P_DEFAULT,(void*)data);
    H5Sclose(mem_space_id);
}
/*
void icy::Serializer::Write(std::vector<double> &node_buffer,
                            std::vector<int> &elems_buffer,
                            unsigned offset_nodes, unsigned offset_elems)
{
    if(!fileIsOpen) return;
    std::size_t nNodes = node_buffer.size()/NodeDataFields;
    std::size_t nElems = elems_buffer.size()/3;

    unsigned long nodes_extent = offset_nodes+nNodes;

    // set extent
    hsize_t nodes_dims[2] = {nodes_extent,NodeDataFields};
    H5Dset_extent(ds_nodes_handle, nodes_dims);

    // write
    hsize_t mem_dims[2] = {nNodes,NodeDataFields};
    hid_t mem_space_id = H5Screate_simple(2, mem_dims, mem_dims);

    hid_t file_space_id = H5Dget_space(ds_nodes_handle);
    hsize_t offset[2] = {offset_nodes,0};
    hsize_t count[2] = {nNodes,NodeDataFields};
    H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, offset, NULL, count, NULL);

    H5Dwrite(ds_nodes_handle, H5T_NATIVE_DOUBLE, mem_space_id, file_space_id,
             H5P_DEFAULT,(void*)node_buffer.data());
    H5Sclose(mem_space_id);
    H5Sclose(file_space_id);

    unsigned long elems_extent = offset_elems+nElems;

    // set extent
    hsize_t elems_dims[2] = {elems_extent,3};
    H5Dset_extent(ds_elems_handle, elems_dims);

    // write
    hsize_t mem_dims2[2] = {nElems,3};
    hid_t mem_space_id2 = H5Screate_simple(2, mem_dims2, mem_dims2);

    hid_t file_space_id2 = H5Dget_space(ds_elems_handle);
    hsize_t offset2[2] = {offset_elems,0};
    H5Sselect_hyperslab(file_space_id2, H5S_SELECT_SET, offset2, NULL, mem_dims2, NULL);

    H5Dwrite(ds_elems_handle, H5T_NATIVE_INT, mem_space_id2, file_space_id2,
             H5P_DEFAULT,(void*)elems_buffer.data());
    H5Sclose(mem_space_id2);
    H5Sclose(file_space_id2);

}

void icy::Serializer::Read(std::vector<double> &node_buffer,
                            std::vector<int> &elems_buffer,
                           unsigned offset_nodes, unsigned offset_elems,
                           unsigned nNodes, unsigned nElems)
{
    node_buffer.resize(nNodes*NodeDataFields);
    elems_buffer.resize(nElems*3);

    hid_t file_space_id = H5Dget_space(ds_nodes_handle);
    hsize_t offset[2] = {offset_nodes,0};
    hsize_t count[2] = {nNodes,NodeDataFields};
    H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    hid_t mem_space_id = H5Screate_simple(2, count, NULL);
    H5Dread(ds_nodes_handle, H5T_NATIVE_DOUBLE, mem_space_id,
            file_space_id, H5P_DEFAULT, (void*)node_buffer.data());
    H5Sclose(file_space_id);
    H5Sclose(mem_space_id);

    hid_t file_space_id2 = H5Dget_space(ds_elems_handle);
    hsize_t offset2[2] = {offset_elems,0};
    hsize_t count2[2] = {nElems,3};
    H5Sselect_hyperslab(file_space_id2, H5S_SELECT_SET, offset2, NULL, count2, NULL);
    hid_t mem_space_id2 = H5Screate_simple(2, count2, NULL);
    H5Dread(ds_elems_handle, H5T_NATIVE_INT, mem_space_id2,
            file_space_id2, H5P_DEFAULT, (void*)elems_buffer.data());
    H5Sclose(file_space_id2);
    H5Sclose(mem_space_id2);
}
*/

void icy::Serializer::Trim(unsigned steps_extent, unsigned nodes_extent, unsigned elems_extent)
{
    if(!fileIsOpen) return;

    // trim file
    hsize_t steps_dims[2] = {steps_extent,sizeof(icy::FrameInfo)};
    H5Dset_extent(ds_steps_handle, steps_dims);

    hsize_t nodes_dims[2] = {nodes_extent,NodeDataFields};
    H5Dset_extent(ds_nodes_handle, nodes_dims);

    hsize_t elems_dims[2] = {elems_extent,3};
    H5Dset_extent(ds_elems_handle, elems_dims);
}

void icy::Serializer::WriteSteps(unsigned steps_extent, FrameInfo *f)
{
    if(!fileIsOpen) return;
    hsize_t current_dims[2] = {steps_extent,sizeof(icy::FrameInfo)};

    H5Dset_extent(ds_steps_handle, current_dims);

    hsize_t frameinfo_dims[2] = {1,sizeof(icy::FrameInfo)};
    hid_t mem_space_id = H5Screate_simple(2, frameinfo_dims, frameinfo_dims);

    //hid_t file_space_id = H5Screate_simple(2, current_dims, NULL);
    hid_t file_space_id = H5Dget_space(ds_steps_handle);

    hsize_t offset[2] = {steps_extent-1,0};
    hsize_t count[2] = {1, sizeof(icy::FrameInfo)};
    H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, offset, NULL, count, NULL);

    H5Dwrite(ds_steps_handle, H5T_NATIVE_CHAR, mem_space_id, file_space_id, H5P_DEFAULT,(void*)f);
    H5Sclose(mem_space_id);
    H5Sclose(file_space_id);

}

void icy::Serializer::ReadSteps(std::vector<icy::FrameInfo> &stepStats)
{
    if(!fileIsOpen) return;

    hid_t filespace = H5Dget_space(ds_steps_handle);
    hsize_t dims[2]; // dataset dimensions
    H5Sget_simple_extent_dims(filespace, dims, NULL);
    H5Sclose(filespace);
    unsigned size = dims[0];
    stepStats.resize(size);
    void *ptr = (void*)stepStats.data();
    H5Dread(ds_steps_handle, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, ptr);
}

/*
void icy::Serializer::WriteAll(std::vector<double> &node_buffer,
                               std::vector<int> &elems_buffer,
              unsigned offset_nodes, unsigned offset_elems,
              unsigned steps_extent, FrameInfo *f)
{
    if(!fileIsOpen) return;
    Write(node_buffer, elems_buffer, offset_nodes, offset_elems);
    WriteSteps(steps_extent, f);
}
*/
