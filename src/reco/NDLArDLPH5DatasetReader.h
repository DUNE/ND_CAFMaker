/// \file NDLArDLPH5DatasetReader.h
///
/// Wrapper for retrieving products from a Deep Learn Physics ML Reco HDF5 file
///
/// \author  J. Wolcott <jwolcott@fnal.gov>
/// \date    May 2023


#ifndef ND_CAFMAKER_NDLARDLPH5DATASETREADER_H
#define ND_CAFMAKER_NDLARDLPH5DATASETREADER_H

#include <memory>
#include <string>
#include <typeinfo>
#include <typeindex>
#include <unordered_map>

#include "H5Cpp.h"

#include "DLP_h5_classes.h"
#include "readH5/DatasetBuffer.h"
#include "readH5/H5DataView.h"
#include "readH5/IH5Viewer.h"

namespace cafmaker
{
  class NDLArDLPH5DatasetReader;

  // -----------------------------------------------------------

//  std::type_info Get

  /// Reader for the Deep Learn Physics ML Reco chain's HDF5 files.
  /// Not quite a generic HDF5 reader because
  /// it knows that the 'events' dataset is special
  /// (it holds region references to the other datasets)
  /// and because it knows the structured types it's expecting
  class NDLArDLPH5DatasetReader : public IH5Viewer
  {
    public:
      struct DatasetNames
      {
        std::string event;
        std::string particle;
        std::string interaction;
        std::string trueParticle;
        std::string trueInteraction;
      };

      NDLArDLPH5DatasetReader(const std::string & h5filename,
                              const std::unordered_map<std::type_index, std::string> & datasetNames);

      template <typename T>
      const std::string & GetDatasetName() const
      {
        auto it = this->fDatasetNames.find(std::type_index(typeid(T)));
        if (it == this->fDatasetNames.end())
          throw std::invalid_argument(std::string("Could not find dataset for type: ") + typeid(T).name());

        return it->second;
      }

      /// Retrieve all of the products for a given event index (or all events if given -1)
      template <typename T>
      H5DataView<T> GetProducts(long int evtIdx=-1) const
      {
        // todo: implement a caching mechanism so repeated requests for the same evtIdx don't cause re-reads from the file

        if (fDatasetBuffers.find(typeid(T)) == fDatasetBuffers.end())
          fDatasetBuffers.emplace(typeid(T), std::make_unique<DatasetBuffer<T>>(fInputFile,
                                                                                GetDatasetName<T>(),
                                                                                cafmaker::types::dlp::BuildCompType<T>));

        auto dsBuffer = dynamic_cast<DatasetBuffer<T>*>(fDatasetBuffers.at(typeid(T)).get());

        // the easy case is if the user wants all entries.  no filtering then...
        if (evtIdx < 0)
        {
          dsBuffer->resize(dsBuffer->nEntries);
          dsBuffer->ds.read(dsBuffer->data(), dsBuffer->compType(), H5::DataSpace::ALL, H5::DataSpace::ALL);
        }
        else
        {
          // when it's just the Event object they want, it's a tad simpler
          if constexpr (std::is_same_v<T, cafmaker::types::dlp::Event>)
          {
            dsBuffer->resize(1);
            hsize_t start[1] = {0};
            hsize_t count[1] = {1};
            dsBuffer->dsp.selectHyperslab(H5S_SELECT_SET, count, start);
            // I think the 3rd argument is likely suspect.  it should specify how to map events into memory,
            // but I'm not sure what happens when you give ALL for that and a subselection for the subsequent argument
            dsBuffer->ds.read(dsBuffer->data(), dsBuffer->compType(), H5::DataSpace::ALL, dsBuffer->dsp);

          } // if (T == Event)
          else
          {
            H5DataView<cafmaker::types::dlp::Event> evts = GetProducts<cafmaker::types::dlp::Event>(evtIdx);

            H5::DataSet ds_ref;
            ds_ref.dereference(fInputFile, evts[0].GetRef<T>(),
                               H5R_DATASET_REGION);  // event 0 because we selected using the evtIdx...
            H5::DataSpace ref_region = fInputFile.getRegion(evts[0].GetRef<T>());

            std::size_t newSize = ref_region.getSelectNpoints();
            dsBuffer->resize(newSize);
            H5::DataSpace memspace(ref_region);
            hsize_t start[1] = {0};
            hsize_t count[1] = {newSize};
            memspace.selectHyperslab(H5S_SELECT_SET, count, start);
            // we need two DataSpaces here because the first one specifies how to map the elements into memory
            // (you could in principle want to rearrange them from the input, though we don't want to here)
            // and the second one specifies which elements to pull from the file.
            // in our case we want to use the reference region selection to pull from the file,
            // and to just stuff them all into the vector from the beginning
            ds_ref.read(dsBuffer->data(), dsBuffer->compType(), memspace, ref_region);
          } // else if (T != Event)
        } // else if (evtIdx >= 0)

        H5DataView<T> view = AddView(dsBuffer->bufferaddr());;

        return view;
      } // H5DataView<T> NDLArDLPH5DatasetReader::GetProducts()


      std::string InputFileName() const;

    private:
      H5::H5File  fInputFile;

      std::unordered_map<std::type_index, std::string> fDatasetNames;

      mutable std::unordered_map<std::type_index, std::unique_ptr<DatasetBufferBase>> fDatasetBuffers;
  };
}

#endif //ND_CAFMAKER_NDLARDLPH5DATASETREADER_H
