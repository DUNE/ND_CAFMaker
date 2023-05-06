/// \file BufferView.h
///
/// A viewer type for variable-length HDF5 buffers
///
/// \author  J. Wolcott <jwolcott@fnal.gov>
/// \date    May 2023

#ifndef ND_CAFMAKER_BUFFERVIEW_H
#define ND_CAFMAKER_BUFFERVIEW_H

#include "H5Cpp.h"

namespace cafmaker
{

  /// \brief Viewer type to work with
  template <typename T>
  class BufferView
  {
    public:
      BufferView()
        : fHandle(nullptr)
      {}

      explicit BufferView(const hvl_t* handle)
        : fHandle(handle)
      {}

      class forward_iterator: public std::iterator<std::forward_iterator_tag, T>
      {
        public:
          forward_iterator(const BufferView<T> * buffer, long idx = 0) : fBuffer(buffer), fIdx(idx) {}
          forward_iterator& operator++() {fIdx++; return *this;}
          forward_iterator operator++(int) {forward_iterator retval = *this; ++(*this); return retval;}
          bool operator==(forward_iterator other) const {return fBuffer == other.fBuffer && fIdx == other.fIdx;}
          bool operator!=(forward_iterator other) const {return !(*this == other);}
          const T& operator*() const {return (*fBuffer)[fIdx];}

        private:
          const BufferView<T> * fBuffer;
          std::size_t fIdx;
      };
      forward_iterator begin() const { return forward_iterator(this, 0); }
      forward_iterator end() const { return forward_iterator(this, fHandle->len); }

      const T& operator[](std::size_t i) const  { return *(static_cast<T*>(fHandle->p) + i); }

      void reset(const hvl_t* handle)  { fHandle = handle; }
      std::size_t size() const { return fHandle->len; }


    private:
      const hvl_t * fHandle;
  };

} // cafmaker

#endif //ND_CAFMAKER_BUFFERVIEW_H
