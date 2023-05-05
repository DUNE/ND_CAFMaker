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

  /// Class that tests if a given type has a method SyncVectors().
  /// Used in DatasetBuffer<T>::syncVectors().
  /// (Adapted from https://retroscience.net/cpp-detect-functions-template.html)
  template<typename T>
  class component_has_syncvector
  {
      template<typename> static std::false_type test();
      template<typename U> static auto test() -> decltype(std::declval<U>().SyncVectors(), std::true_type());
    public:
      static constexpr bool value = std::is_same<decltype(test<T>(0)), std::true_type>::value;
  };

  /// Storage class for the buffer used for an HDF structured datatype,
  /// which is mapped to the C++ class specified by the template argument.
  template<typename T>
  class DatasetBuffer : public DatasetBufferBase
  {
    public:
      DatasetBuffer(const H5::H5File &f,
                    const std::string &dsName,
                    const std::function<H5::CompType()> &compTypeBuilder)
        : DatasetBufferBase(f, dsName), fCompType(compTypeBuilder())
      {}

      /// Get a H5 Compound Type instance corresponding to this buffer.
      const H5::CompType &compType() const { return fCompType; }

      /// Get the address of the underlying std::vector buffer
      const std::vector <T> * bufferaddr() const { return &fBuffer; }

      /// ensure any vectors within type T are synchroniz with the HDF5 'handles'
      /// call this after loading data into the buffer...
      void syncVectors()
      {
        if constexpr (component_has_syncvector<T>::value)
        {
          for (const auto & e : fBuffer)
            e.SyncVectors();
        }
      }

      // ape the std::vector interface where relevant
      T *  data()                    { return fBuffer.data(); }
      void resize(std::size_t count) { fBuffer.resize(count); }

    private:
      H5::CompType fCompType;
      std::vector <T> fBuffer;
  };
}
#endif //ND_CAFMAKER_DATASETBUFFER_H
