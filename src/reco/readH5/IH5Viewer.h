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
  template <typename T>
  class H5DataView;

  class H5DataViewBase;

  class IH5Viewer
  {
    public:
      virtual ~IH5Viewer();

      /// Remove a view from the list.  Public so that H5DataViewBase instances can call it.
      void RemoveView(H5DataViewBase* view) const;

    protected:
      /// Add a view and store it in the cache.
      /// Only IH5Viewer derived types are supposed to be able to make views,
      /// so it's protected.
      template <typename T, typename ...Args>
      H5DataView<T> NewView(Args&&...args) const
      {
        H5DataView<T> view(this, std::forward<Args>(args)...);
        fCurrentViews.insert(&view);
        return view;
      }

      /// Add an already-created view into the cache.
      /// Should only be used by the H5DataView copy constructor.
      void AddView(H5DataViewBase* view) const;

    private:
      mutable std::unordered_set<H5DataViewBase*> fCurrentViews;
  };


}

#endif //ND_CAFMAKER_IH5VIEWER_H
