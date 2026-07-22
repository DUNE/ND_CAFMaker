/// \file H5DataView.h
///
/// Utility for handing out access to the contents of an HDF5 Dataset read.
///
/// \author  J. Wolcott <jwolcott@fnal.gov>
/// \date    May 2023

#ifndef ND_CAFMAKER_H5DATAVIEW_H
#define ND_CAFMAKER_H5DATAVIEW_H

#include <memory>
#include <vector>

#include "DatasetBuffer.h"

namespace cafmaker
{

  /// Wrapper class for viewing the contents of an HDF5 dataset
  /// that consists of a sequence of a structured type,
  /// which is mapped to the C++ class passed as the template argument.
  //
  // A view shares ownership of the immutable DatasetBuffer it was made from
  // (via shared_ptr), so the data it points at is guaranteed to stay alive and unchanged 
  // for as long as the view -- or any copy of it -- exists.
  // There is consequently no notion of a view becoming "invalid": if you hold one, it's good.  
  // Independent reads produce independent buffers,
  // so it's fine to hold several views of the same type at once.
  template <typename T>
  class H5DataView
  {
    public:
      H5DataView() = default;

      explicit H5DataView(std::shared_ptr<const DatasetBuffer<T>> buffer)
        : fBuffer(std::move(buffer))
      {}

      const std::vector<T> & operator*() const { return fBuffer->buffer(); }

      // enable use with range-based for
      auto begin() const { return fBuffer->buffer().begin(); }
      auto end()   const { return fBuffer->buffer().end(); }

      // other vector-like operations
      std::size_t size() const { return fBuffer->buffer().size(); }

      const T & operator[](std::size_t idx) const { return fBuffer->buffer()[idx]; }

    private:
      std::shared_ptr<const DatasetBuffer<T>> fBuffer;
  };

}

#endif //ND_CAFMAKER_H5DATAVIEW_H
