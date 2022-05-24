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

  template <>
  const std::vector<std::string> EXPECTED_COLUMN_NAMES<caf::SRTrack>
  {
      "trk_start_x",
      "trk_start_y",
      "trk_start_z",
      "trk_end_x",
      "trk_end_y",
      "trk_end_z",
      "trk_end_dir_x",
      "trk_end_dir_y",
      "trk_end_dir_z",
      "trk_visE"
  };

  // -------------------------------------------------------------
  template <>
  const std::vector<std::string> EXPECTED_COLUMN_NAMES<caf::SRShower>
  {
      "shw_start_x",
      "shw_start_y",
      "shw_start_z",
      "shw_dir_x",
      "shw_dir_y",
      "shw_dir_z",
      "shw_visE"
  };

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

      /// Get the internal dataset reader.  For diagnostics only; prefer the FillSR() interface where possible
      const NDLArSummaryH5DatasetReader & DatasetReader() const  { return fDSReader; }

      /// This will be specialized for each product type that can be handled.
      /// See the implementations in the .cxx.
      void FillSR(caf::StandardRecord &sr, std::size_t evtIdx) const;

    private:
      std::vector<T> EventProducts(std::size_t evt) const
      {
        int firstRow, lastRow;
        std::tie(firstRow, lastRow) = fDSReader.EventRowEdges(evt);
        std::cout << "  event " << evt << " spans rows " << firstRow << ", " << lastRow << std::endl;
        if (firstRow < 0)
          return {};

        const float * buffer = fDSReader.GridValues(firstRow,
                                                    lastRow,
                                                    fDSReader.ProductFirstColumn(),
                                                    fDSReader.ColumnNames().size() - 1);

        std::size_t nRows = lastRow - firstRow + 1;
        std::size_t nCols = fDSReader.ColumnNames().size() - fDSReader.ProductFirstColumn();

        return ParseVals(buffer, nRows, nCols);
      }

      /// This will be specialized for each product type that can be handled.
      /// See the implementations in the .cxx.
      static std::vector<T> ParseVals(const float * buffer, std::size_t nRows, std::size_t nCols);

      /// Compare the columns read back from the dataset to those in an expectation list
      void ValidateColumns(const std::vector<std::string> & expectedColumns) const
      {
        bool match = true;
        std::size_t startCol = this->fDSReader.ProductFirstColumn();
        if (expectedColumns.size() != this->fDSReader.ColumnNames().size() - this->fDSReader.ProductFirstColumn())
          match = false;
        else
        {
          for (std::size_t col = startCol; col < this->fDSReader.ColumnNames().size(); col++)
          {
            std::cout << "fDSReader.ColumnNames()[" << col << "] = " << this->fDSReader.ColumnNames()[col]
                      << "; expectedColumns[" << col - startCol << "] = " << expectedColumns[col - startCol]
                      << std::endl;
            if (this->fDSReader.ColumnNames()[col] != expectedColumns[col - startCol])
            {
              match = false;
              break;
            }
          } // for (col)
        } // else (column sizes do match)

        if (!match)
        {
          std::cerr << "Column names read from dataset '" << this->fDSReader.InputDatasetName()
                    << "' in file '" << this->fDSReader.InputFileName() << "'"
                    << " don't match expected names!" << std::endl;
          std::cerr << "Expected " << expectedColumns.size() << " columns:";
          for (const auto & c : expectedColumns)
            std::cerr << "  " << c;
          std::cerr << std::endl;
          std::cerr << this->fDSReader.ColumnNames().size() - this->fDSReader.ProductFirstColumn()
                    << " columns read from dataset:" << std::endl;
          for (std::size_t col = startCol; col < this->fDSReader.ColumnNames().size(); col++)
            std::cerr << "  " << this->fDSReader.ColumnNames()[col];
          std::cerr << std::endl;

          abort();
        } // if(columns don't match)
      } // ValidateColumns()

      NDLArSummaryH5DatasetReader fDSReader;

  };
}

#endif //ND_CAFMAKER_NDLARPRODUCTFILLER_H
