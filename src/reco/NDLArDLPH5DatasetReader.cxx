#include "NDLArDLPH5DatasetReader.h"

namespace cafmaker
{
  // -----------------------------------------------------------
  // -----------------------------------------------------------

  NDLArDLPH5DatasetReader::NDLArDLPH5DatasetReader(const std::string &h5filename,
                                                   const std::unordered_map<std::type_index, std::string> &datasetNames)
    : fInputFile(h5filename, H5F_ACC_RDONLY), fDatasetNames(datasetNames)
  {}

  // -----------------------------------------------------------

  std::string NDLArDLPH5DatasetReader::InputFileName() const
  {
    return fInputFile.getFileName();
  }

  // -----------------------------------------------------------





} // namespace cafmaker
