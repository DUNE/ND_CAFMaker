#include "IH5Viewer.h"
#include "H5DataView.h"

namespace cafmaker
{
  IH5Viewer::~IH5Viewer()
  {
    // make sure any child views know they can't be used any more
    for (H5DataViewBase *view: fCurrentViews)
      view->invalidate();
  }

  // -----------------------------------------------------------


  void IH5Viewer::RemoveView(H5DataViewBase* view) const
  {
    fCurrentViews.erase(view);
  }


}
