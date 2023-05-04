/// \file H5DataView.h
///
/// Utility for creating a stateful "view" into an HDF5 Dataset
///
/// \author  J. Wolcott <jwolcott@fnal.gov>
/// \date    May 2023

#ifndef ND_CAFMAKER_H5DATAVIEW_H
#define ND_CAFMAKER_H5DATAVIEW_H

#include <stdexcept>
#include <vector>

namespace cafmaker
{

  class IH5Viewer;

  /// Base type, non-templated common stuff for H5DataViews
  class H5DataViewBase
  {
    public:
      H5DataViewBase(const IH5Viewer * viewer);

      // the IH5Viewer needs to find out about the new object when copying
      H5DataViewBase(const H5DataViewBase&);

      virtual ~H5DataViewBase();

      bool valid() const  { return fValid; };    ///< Check if this view is valid.  If not, it should be discarded
      void invalidate()   { fValid = false; }    ///< Set this view status to invalid

    private:
      bool fValid = true;
      const IH5Viewer * fViewer;
  };

  // -----------------------------------------------------------

  /// Wrapper class for viewing the contents of an HDF5 dataset
  /// that consists of a sequence of a structured type,
  /// which is mapped to the C++ class passed as the template argument.
  /// This viewer maintains an internal state corresponding to whether the view is valid.
  template <typename T>
  class H5DataView : public H5DataViewBase
  {
    // only IH5Viewers should be making H5DataViews
    friend class IH5Viewer;

    public:
      const std::vector<T> & operator*() const
      {
        if(valid())
          return *fBuffer;
        throw std::runtime_error("H5DataView is invalid");
      }

      // enable use with range-based for
      const auto begin() const { return fBuffer->begin(); }
      const auto end()   const { return fBuffer->end(); }

      const T & operator[](std::size_t idx) const
      {
        if(valid())
          return (*fBuffer)[idx];
        throw std::runtime_error("H5DataView is invalid");
      }

    private:
      H5DataView(const IH5Viewer * reader, const std::vector<T> * vec)
      : H5DataViewBase(reader), fBuffer(vec)
      {}

      const std::vector<T> * fBuffer;
  };

}

#endif //ND_CAFMAKER_H5DATAVIEW_H
