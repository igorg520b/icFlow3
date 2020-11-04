#ifndef SERIALIZER_H
#define SERIALIZER_H

#include <string>
#include <vector>
#include "hdf5.h"


// the reason for having a separate Serializer class is to
// separate the "stable" chunk of code and avoid recompilation of h5cpp

namespace icy { class Serializer; class FrameInfo; }

class icy::Serializer
{
public:
    bool fileIsOpen=false;

    hid_t file_handle, ds_nodes_handle, ds_elems_handle, ds_steps_handle, ds_params_handle;

    const unsigned NodeDataFields = 19;    // from node.h

    // save/load mesh
    void CreateFile(std::string fileName, unsigned params_size);
    void OpenFile(std::string fileName);
    void CloseFile();
    void SaveParams(void *data, unsigned length);
    void LoadParams(void *data, unsigned max_length);
/*
    void Write(std::vector<double> &node_buffer,
               std::vector<int> &elems_buffer,
               unsigned offset_nodes, unsigned offset_elems);

    void Read(std::vector<double> &node_buffer,
              std::vector<int> &elems_buffer,
             unsigned offset_nodes, unsigned offset_elems,
             unsigned nNodes, unsigned nElems);

    void WriteAll(std::vector<double> &node_buffer,
                  std::vector<int> &elems_buffer,
                  unsigned offset_nodes, unsigned offset_elems,
                  unsigned steps_extent, FrameInfo *f);

*/
    void WriteSteps(unsigned steps_extent, FrameInfo *f);
    void ReadSteps(std::vector<icy::FrameInfo> &stepStats);

    void Trim(unsigned steps_extent, unsigned nodes_extent, unsigned elems_extent);
};

#endif // SERIALIZER_H
