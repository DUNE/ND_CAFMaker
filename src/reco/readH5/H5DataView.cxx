#include "H5DataView.h"
#include "IH5Viewer.h"

namespace cafmaker
{
  // -----------------------------------------------------------

  H5DataViewBase::H5DataViewBase(const IH5Viewer * viewer)
  : fViewer(viewer)
  {}

  // -----------------------------------------------------------

  H5DataViewBase::~H5DataViewBase()
  {
    // notify our parent reader object that we're not around any longer
    if (fValid)
      fViewer->RemoveView(this);
  }


}
