/// \file NDLArSummaryH5.h
///
/// Wrapper for a ND-LAr summary H5 file
///

#ifndef ND_CAFMAKER_NDLARSUMMARYH5_H
#define ND_CAFMAKER_NDLARSUMMARYH5_H

#include <array>
#include <cassert>
#include <string>
#include <set>
#include <vector>

#include "H5Cpp.h"
#include "TVector3.h"

#include "duneanaobj/StandardRecord/SRTrack.h"

namespace cafmaker
{
  class NDLArSummaryH5
  {
    public:
      NDLArSummaryH5(const std::string & h5filename,
                     const std::string & h5dataset,
                     const std::string & column_name_attr = "column_names");

      /// What are the known events in this file?
      std::set<std::size_t> Events() const;

      /// Get the tracks from a particular event within the file.
      std::vector<caf::SRTrack> EventTracks(std::size_t event) const;

      const std::string & GetEventColumnName() const   { return fEventColumnName; };
      void SetEventColumnName(const std::string & ecn) { fEventColumnName = ecn; };

    private:
      std::vector<float> ColumnValues(std::size_t colIdx, std::size_t firstRow = 0, int lastRow = -1) const;

      const std::vector<std::size_t> & EventRows() const;

      /// Read a rectangular grid of numbers from the h5 file.
      /// ColumnValues() is a shortcut if you want all the values in a column
      void GridValues(float *buffer,
                      std::size_t startRow, std::size_t endRow,
                      std::size_t startCol, std::size_t endCol) const;

      void ReadColumnNames();

      H5::H5File  fInputFile;
      H5::DataSet fInputDataset;
      H5::DataSpace fInputDataspace;

      std::string fColumnNameAttr;     ///<  name of the attribute inside the dataset that gives the dataset's column names
      std::string fEventColumnName;    ///<  what's the string stored in the column names that tells us which one is the event number?

      std::vector<std::string> fColumnNames;
      mutable std::vector<std::size_t> fRowEvents;    ///<  which event each row corresponds to
  };
}

#endif //ND_CAFMAKER_NDLARSUMMARYH5_H
