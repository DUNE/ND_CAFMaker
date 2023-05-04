/// \file IH5Viewer.h
///
/// Base class for HDF5 readers that want to use the H5DataView tool
/// to track views into their buffers
///
/// \author  J. Wolcott <jwolcott@fnal.gov>
/// \date    May 2023

#ifndef ND_CAFMAKER_IH5VIEWER_H
#define ND_CAFMAKER_IH5VIEWER_H

#include <unordered_set>

namespace cafmaker
{
  class H5DataViewBase;

  class IH5Viewer
  {
    public:
      virtual ~IH5Viewer();

      /// Remove a view from the list.  Public so that H5DataViewBase instances can call it.
      void RemoveView(H5DataViewBase* view) const;

    protected:
      /// Add a view from the list.  Only IH5Viewer derived types should be able to make these
      void AddView(H5DataViewBase* view) const;

    private:
      mutable std::unordered_set<H5DataViewBase*> fCurrentViews;
  };


}

#endif //ND_CAFMAKER_IH5VIEWER_H
