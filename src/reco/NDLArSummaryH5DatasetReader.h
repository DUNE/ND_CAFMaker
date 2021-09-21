/// \file NDLArSummaryH5.h
///
/// Wrapper for retrieving products from a ND-LAr summary H5 file
///
/// \author  J. Wolcott <jwolcott@fnal.gov>
/// \date    Sept. 2021


#ifndef ND_CAFMAKER_NDLARSUMMARYH5DATASETREADER_H
#define ND_CAFMAKER_NDLARSUMMARYH5DATASETREADER_H

#include <array>
#include <cassert>
#include <string>
#include <set>
#include <vector>

#include "H5Cpp.h"
#include "TVector3.h"

#include "duneanaobj/StandardRecord/SRShower.h"
#include "duneanaobj/StandardRecord/SRTrack.h"

namespace cafmaker
{
  class NDLArSummaryH5DatasetReader
  {
    public:
      NDLArSummaryH5DatasetReader(const std::string & h5filename,
                                  const std::string & h5dataset,
                                  const std::string & column_name_attr = "column_names",
                                  const std::string & evt_col_name = "Event");

      const std::vector<std::string> & ColumnNames() const { return fColumnNames; }

      /// What are the known events in this file?
      std::set<std::size_t> Events() const;

      /// What row numbers in the dataset are the upper & lower bounds for the given event?
      std::pair<std::size_t, std::size_t> EventRowEdges(std::size_t event) const;

      /// Read a rectangular grid of numbers from the h5 file.
      /// ColumnValues() is a shortcut if you want all the values in a column
      const float * GridValues(std::size_t startRow, std::size_t endRow,
                               std::size_t startCol, std::size_t endCol) const;

      std::string InputDatasetName() const  { return fInputDataset.getObjName(); }
      std::string InputFileName() const     { return fInputDataset.getFileName(); }

      /// Return the indices of the product information columns
      std::size_t ProductFirstColumn() const;

    private:
      /// Get the values contained in an entire column
      std::vector<float> ColumnValues(std::size_t colIdx, std::size_t firstRow = 0, int lastRow = -1) const;

      /// Get the event number corresponding to each row
      const std::vector<std::size_t> & EventRowMap() const;

      void ReadColumnNames();

      H5::H5File  fInputFile;
      H5::DataSet fInputDataset;
      H5::DataSpace fInputDataspace;

      std::string fColumnNameAttr;     ///<  name of the attribute inside the dataset that gives the dataset's column names
      std::string fEventColumnName;    ///<  what's the string stored in the column names that tells us which one is the event number?

      std::vector<std::string> fColumnNames;
      mutable std::vector<std::size_t> fRowEvents;    ///<  which event each row corresponds to

      mutable int fEvtColumnIdx = -1;
      mutable std::vector<float> fReadBuffer;
  };
}

#endif //ND_CAFMAKER_NDLARSUMMARYH5DATASETREADER_H
