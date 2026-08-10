/// \file DatasetBuffer.h
///
/// Utility for storing an HDF5 Dataset object
/// and the associated buffer we'll use to read it
///
/// \author  J. Wolcott <jwolcott@fnal.gov>
/// \date    May 2023

#ifndef ND_CAFMAKER_DATASETBUFFER_H
#define ND_CAFMAKER_DATASETBUFFER_H

#include <functional>
#include <vector>

#include "H5Cpp.h"

namespace cafmaker
{
  /// Base for the dataset buffer storage containing non-templated shared stuff
  struct DatasetBufferBase
  {
    DatasetBufferBase(const H5::H5File &f, const std::string &dsName);
    virtual ~DatasetBufferBase() = default;

    H5::DataSet ds;
    H5::DataSpace dsp;

    std::size_t nEntries;   //< loaded from dataset
  };

  /// Storage class for the buffer used for an HDF structured datatype,
  /// which is mapped to the C++ class specified by the template argument.
  ///
  /// One DatasetBuffer owns exactly one read's worth of data.
  /// It is meant to be filled once (by the reader) and thereafter treated as immutable, 
  /// then access to it handed out via a shared_ptr inside an H5DataView.  
  /// Because the mapped type T may contain variable-length members (BufferView<>s / char* strings)
  /// whose backing memory is allocated by HDF5 during read(), 
  /// the destructor reclaims that memory -- so the whole read result is cleaned up 
  /// as soon as the last view referencing it goes away.
  template<typename T>
  class DatasetBuffer : public DatasetBufferBase
  {
    public:
      //todo: make a static_assert<> checking that T is an appropriate type

      DatasetBuffer(const H5::H5File &f,
                    const std::string &dsName,
                    const std::function<H5::CompType()> &compTypeBuilder)
        : DatasetBufferBase(f, dsName), fCompType(compTypeBuilder())
      {}

      ~DatasetBuffer() override
      {
        // HDF5 allocates the memory behind any variable-length members
        // (vlen buffers *and* variable-length strings) during read().
        // We own it, so we must give it back.  vlenReclaim() walks the buffer
        // according to the compound type and frees those allocations; 
        // it is a no-op for types with no variable-length members, so it's always safe.
        if (fBuffer.empty())
          return;
        hsize_t dim = fBuffer.size();
        H5::DataSpace space(1, &dim);
        H5::DataSet::vlenReclaim(fBuffer.data(), fCompType, space);
      }

      // owns its data -- no copies (copying would double-free on reclaim)
      DatasetBuffer(const DatasetBuffer &) = delete;
      DatasetBuffer &operator=(const DatasetBuffer &) = delete;

      /// Get a H5 Compound Type instance corresponding to this buffer.
      const H5::CompType &compType() const { return fCompType; }

      /// Read-only access to the underlying std::vector buffer (for H5DataView)
      const std::vector <T> & buffer() const { return fBuffer; }

      /// ensure any vectors within type T are synchronized with the HDF5 'handles'
      /// call this after loading data into the buffer...
      void syncVectors()
      {
        for (auto & e : fBuffer)
          e.SyncVectors();
      }

      // ape the std::vector interface where relevant
      T *  data()                    { return fBuffer.data(); }
      std::size_t size() const       { return fBuffer.size(); }
      void resize(std::size_t count) { fBuffer.resize(count); }

    private:
      H5::CompType fCompType;
      std::vector <T> fBuffer;
  };
}
#endif //ND_CAFMAKER_DATASETBUFFER_H
