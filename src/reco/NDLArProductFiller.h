/// \file NDLArProductFiller.h
///
/// Fill CAF reco branches using DeepLearnPhysics machine learning based reconstruction.
///
/// \author  J. Wolcott <jwolcott@fnal.gov>
/// \date    Sept. 2021

#ifndef ND_CAFMAKER_NDLARPRODUCTFILLER_H
#define ND_CAFMAKER_NDLARPRODUCTFILLER_H

#include <iostream>
#include <tuple>

#include "NDLArSummaryH5DatasetReader.h"

namespace caf
{
  class StandardRecord;
}

namespace cafmaker
{
  // note: these specialized in the .cxx
  template <typename T>
  const std::vector<std::string> EXPECTED_COLUMN_NAMES;

  // -----------------------------------------

  /// Class to copy info from array read out of hdf5 into CAF StanardRecord branches
  template<typename T>
  class NDLArProductFiller
  {
    public:
      NDLArProductFiller(const std::string &h5filename,
                         const std::string &product_name,
                         const std::string &evt_col_name = "Event",
                         const std::string &column_name_attr = "column_names")
        : fDSReader(h5filename, product_name , column_name_attr, evt_col_name)
      {
        this->ValidateColumns(EXPECTED_COLUMN_NAMES<T>);
      }

      /// This will be specialized for each product type that can be handled.
      /// See the implementations in the .cxx.
      void FillSR(caf::StandardRecord &sr, std::size_t evtIdx) const;

    private:
      std::vector<T> EventProducts(std::size_t evt) const
      {
        std::size_t firstRow, lastRow;
        std::tie(firstRow, lastRow) = fDSReader.EventRowEdges(evt);
        if (firstRow < 0)
          return {};

        const float * buffer = fDSReader.GridValues(firstRow,
                                                    lastRow,
                                                    fDSReader.ProductFirstColumn(),
                                                    fDSReader.ColumnNames().size() - 1);

        std::size_t nRows = lastRow - firstRow + 1;
        std::size_t nCols = fDSReader.ColumnNames().size() - fDSReader.ProductFirstColumn() + 1;

        return ParseVals(buffer, nRows, nCols);
      }

      /// This will be specialized for each product type that can be handled.
      /// See the implementations in the .cxx.
      static std::vector<T> ParseVals(const float * buffer, std::size_t nRows, std::size_t nCols);

      /// Compare the columns read back from the dataset to those in an expectation list
      void ValidateColumns(const std::vector<std::string> & expectedColumns) const
      {
        if (this->fDSReader.ColumnNames() != expectedColumns)
        {
          std::cerr << "Column names read from dataset '" << this->fDSReader.InputDatasetName()
                    << "' in file '" << this->fDSReader.InputFileName() << "'"
                    << " don't match expected names!" << std::endl;
          std::cerr << "Expected columns:";
          for (const auto & c : expectedColumns)
            std::cerr << "  " << c;
          std::cerr << std::endl;
          std::cerr << "Columns read from dataset:" << std::endl;
          for (const auto & c : this->fDSReader.ColumnNames())
            std::cerr << "  " << c;
          std::cerr << std::endl;

          abort();
        } // if(columns don't match)
      } // ValidateColumns()

      NDLArSummaryH5DatasetReader fDSReader;

  };
}

#endif //ND_CAFMAKER_NDLARPRODUCTFILLER_H
