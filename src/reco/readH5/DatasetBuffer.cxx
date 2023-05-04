
#include "DatasetBuffer.h"

namespace cafmaker
{
  DatasetBufferBase::DatasetBufferBase(const H5::H5File &f, const std::string &dsName)
    : ds(f.openDataSet(dsName)), dsp(ds.getSpace())
  {
    hsize_t dimsMax = dsp.getSimpleExtentNdims();
    auto dims = new hsize_t[dimsMax];
    dsp.getSimpleExtentDims(dims);
    nEntries = dims[dimsMax - 1];
    delete[] dims;
  }


}
