#include "NDLArSummaryH5.h"

#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <set>

#include <iostream>


namespace cafmaker
{
  NDLArSummaryH5::NDLArSummaryH5(const std::string & h5filename,
                                 const std::string & h5dataset,
                                 const std::string & column_name_attr)
      : fInputFile(h5filename.c_str(), H5F_ACC_RDONLY),
        fInputDataset(fInputFile.openDataSet(h5dataset)),
        fInputDataspace(fInputDataset.getSpace()),
        fColumnNameAttr(column_name_attr)
  {
    ReadColumnNames();
  }

  // -----------------------------------------------------------

  std::vector<float> NDLArSummaryH5::ColumnValues(std::size_t colIdx, std::size_t firstRow, int lastRow) const
  {
    std::cout << "Reading column " << colIdx << std::endl;

    // in principle just taking the number of rows verbatim from the file
    // could cause us to run out of memory if it's a REALLY large number of rows.
    // however that seems pretty unlikely when handed input files corresponding to single subruns,
    // which is what we assume here...
    if (lastRow < 0)
    {
      hsize_t dims[2];
      fInputDataspace.getSimpleExtentDims(dims);
      lastRow = dims[0] - 1;  // because dims[0] is the *count* of rows
    }
    std::size_t nRows = lastRow - firstRow + 1;

    std::cout << "Reading rows " << firstRow << " through " << lastRow << std::endl;

    // if we somehow wound up with an empty array?
    if (static_cast<std::size_t>(lastRow) <= firstRow)
      return {};

    // if the file layout is chunked we should only read one chunk at a time
    std::size_t chunkSize = nRows;
    H5::DSetCreatPropList cparms = fInputDataset.getCreatePlist();
    if (cparms.getLayout() == H5D_CHUNKED)
    {
      hsize_t chunk_dims[2];
      cparms.getChunk( 2, chunk_dims);

      if (chunk_dims[0] < chunkSize)
        chunkSize = chunk_dims[0];

      assert(chunk_dims[1] >= colIdx);
    }

    std::vector<float> column(nRows);
    hsize_t col_dims[1];
    col_dims[0] = chunkSize;
    hsize_t coll_offset[1];
    coll_offset[0] = 0;
    H5::DataSpace memspace(1, col_dims);  // rank = 1 (1-D) == 1 column

    std::cout << " memspace has " << memspace.getSimpleExtentNpoints() << " points" << std::endl;
    for (std::size_t chunkIdx = 0; chunkIdx < nRows / chunkSize + 1; chunkIdx++)
    {
      std::size_t startRow = firstRow + chunkIdx * chunkSize;
      std::size_t chunkRows = std::min(chunkSize, nRows - startRow);
      std::cout << "  reading chunk: rows " << startRow << " through " << startRow + chunkRows - 1 << std::endl;
      hsize_t offset[2]{startRow, colIdx};  // i.e., start from the first row, and the column of interest
      hsize_t count[2]{chunkRows, 1};       // i.e., read all the rows in this chunk, and a single column

      if (chunkRows != chunkSize)
      {
        col_dims[0] = chunkRows;
        memspace.selectHyperslab(H5S_SELECT_SET, col_dims, coll_offset);
      }

      fInputDataspace.selectHyperslab(H5S_SELECT_SET, count, offset);
      std::cout << "   input dataset has " << fInputDataspace.getSelectNpoints() << " points selected" << std::endl;

      fInputDataset.read(&column[startRow],
                         H5Tget_native_type(fInputDataset.getDataType().getId(), H5T_DIR_DEFAULT),
                         memspace,
                         fInputDataspace);
    }

    return column;
  } // NDLArSummaryH5::ColumnValues

  // -----------------------------------------------------------

  void NDLArSummaryH5::GridValues(float * buffer,
                                  std::size_t startRow, std::size_t endRow,
                                  std::size_t startCol, std::size_t endCol) const
  {
    std::cout << "Reading 2D grid from (row, col) = (" << startRow << "," << startCol << ") "
              << " to (" << endRow << "," << endCol << ") "
              << std::endl;

    std::size_t nRows = endRow - startRow + 1;
    std::size_t nCols = endCol - startCol + 1;

    // if the file layout is chunked we should only read one chunk at a time.
    // n.b. we only work with chunking along the row axis...
    std::size_t chunkSize = nRows;
    H5::DSetCreatPropList cparms = fInputDataset.getCreatePlist();
    if (cparms.getLayout() == H5D_CHUNKED)
    {
      hsize_t chunk_dims[2];
      cparms.getChunk( 2, chunk_dims);

//      std::cout << "Chunk dimensions: " << chunk_dims[0] << " x " << chunk_dims[1] << std::endl;
      if (chunk_dims[0] < chunkSize)
        chunkSize = chunk_dims[0];

      if (chunk_dims[1] < endCol)
      {
        std::cerr << "Warning: chunk's horizontal size is " << chunk_dims[1] << ", which is less than the last column to read: "
                  << endCol << std::endl;
        std::cout << "This may impact data read speed!" << std::endl;
      }
    }

    hsize_t * grid_dims;
    int rank = 2;
    if (chunkSize == 1)
    {
      rank = 1;
      grid_dims = new hsize_t[1];
      grid_dims[0] = endCol - startCol + 1;
      std::cout << "   memspace is 1-dimensional with size " << grid_dims[0] << std::endl;
    }
    else if (endCol - startCol == 0)
    {
      rank = 1;
      grid_dims = new hsize_t[1];
      grid_dims[0] = chunkSize;
      std::cout << "   memspace is 1-dimensional with size " << grid_dims[0] << std::endl;
    }
    else
    {
      grid_dims = new hsize_t[2];
      grid_dims[0] = chunkSize;
      grid_dims[1] = endCol - startCol + 1;
      std::cout << "   memspace has dims. " << grid_dims[0] << " x " << grid_dims[1] << std::endl;
    }
    H5::DataSpace memspace(rank, grid_dims);

    std::cout << "    memspace is rank " << rank << std::endl;
    std::cout << "    memspace has " << memspace.getSimpleExtentNpoints() << " points" << std::endl;
    std::cout << "    total n rows = " << nRows << "; initial row chunk size = " << chunkSize << std::endl;
    for (std::size_t chunkIdx = 0; ; chunkIdx++)
    {
      // i.e.: either startRow, or the aligned start of a chunk block.
      // note this is integer division (i.e. floor())!
      std::size_t firstRow = std::max(startRow, (startRow/chunkSize + chunkIdx) * chunkSize);
      if (firstRow > endRow)
        break;

      // i.e.: from the firstRow until the end of the chunk block or the last row, whichever is first.
      // again note integer division (== floor())
      std::size_t chunkRows = std::min( (firstRow/chunkSize + 1) * chunkSize - startRow, endRow - firstRow + 1);

      std::cout << "     this chunk nrows = " << chunkRows << std::endl;
      std::cout << "     reading chunk: rows " << firstRow << " through " << firstRow + chunkRows - 1 << std::endl;
      hsize_t offset[2]{firstRow, startCol};
      hsize_t count[2]{chunkRows, nCols};

      // ensure the memory space map is the same size as the chunk we're going to read.
      // since we originally built it to be the full size, if this block is smaller
      // (perhaps the first block starts or the last one ends in the middle of a chunk),
      // we might need to resize
      hsize_t grid_start[2];
      hsize_t grid_end[2];
      memspace.getSelectBounds(grid_start, grid_end);
      if (rank > 1 && chunkRows != (grid_end[0] - grid_start[0] + 1))
      {
        grid_dims[0] = chunkRows;
        grid_dims[1] = count[1];
        hsize_t grid_offset[2];
        grid_offset[0] = grid_offset[1] = 0;
//        std::cout << " resizing memspace to " << grid_dims[0] << " x " << grid_dims[1] << std::endl;
        memspace.selectHyperslab(H5S_SELECT_SET, grid_dims, grid_offset);
        std::cout << " memspace now has " << memspace.getSimpleExtentNpoints() << " points" << std::endl;
      }

      fInputDataspace.selectHyperslab(H5S_SELECT_SET, count, offset);
      std::cout << "   input dataset has " << fInputDataspace.getSelectNpoints() << " points selected" << std::endl;

      fInputDataset.read(&buffer[(firstRow - startRow) * nCols],
                         H5Tget_native_type(fInputDataset.getDataType().getId(), H5T_DIR_DEFAULT),
                         memspace,
                         fInputDataspace);

    } // for (chunkIdx)

    delete [] grid_dims;

  } // NDLArSummaryH5::GridValues

  // -----------------------------------------------------------

  std::set<std::size_t> NDLArSummaryH5::Events() const
  {
    return { EventRows().begin(), EventRows().end() };
  }

  // -----------------------------------------------------------

  const std::vector<size_t> & NDLArSummaryH5::EventRows() const
  {
    // only do it once per file
    if (fRowEvents.empty())
    {
      // which column tells me what event the other items below came from
      std::size_t evtColumnIdx = std::distance(fColumnNames.begin(),
                                               std::find(fColumnNames.begin(), fColumnNames.end(), "Event"));

      // read them all
      auto colVals = ColumnValues(evtColumnIdx);
      fRowEvents.insert(fRowEvents.end(), colVals.begin(), colVals.end()); // convert to size_t here

      std::cout << "fRowEvents thinks it is " << fRowEvents.size() << "elements long now." << std::endl;
      std::cout << "  its first element is " << fRowEvents[0] << std::endl;
      std::cout << "  its last element (idx " << fRowEvents.size() - 1 << ") is " << fRowEvents[fRowEvents.size() - 1] << std::endl;
      std::cout << "     ... if you use .back() to obtain it: " << fRowEvents.back() << std::endl;
    }

    return fRowEvents;
  }


  // -----------------------------------------------------------

  std::vector<Track> NDLArSummaryH5::EventTracks(std::size_t event) const
  {
    // determine the range of rows corresponding to this event.
    // we could do something smarter and keep track of it across calls
    // to this method, but I seriously doubt this will ever
    // be the limiting step.
    int firstRow = -1;
    int lastRow = -1;

    // search for the range of rows corresponding to this event.
    // look from both ends.  stop when we find the first instance
    for (std::size_t counter = 0; counter < EventRows().size() / 2; counter++)
    {
      if (EventRows()[counter] == event)
      {
        std::cout << " index " << counter << " starts block for event " << event << std::endl;
        firstRow = counter;
        break;
      }

      if (EventRows()[EventRows().size() - counter - 1] == event)
      {
        std::cout << " index " << EventRows().size() - counter - 1 << " ends block for event " << event << std::endl;
        std::cout << "    EventRows()[" << EventRows().size() << " - " << counter << " - 1] = " << EventRows()[EventRows().size() - counter - 1] << std::endl;
        lastRow = EventRows().size() - counter - 1;
        break;
      }
    }

    // we didn't find this event in the collection
    if (firstRow == -1 && lastRow == -1)
      return {};

    // now walk through until you find the other end of the range
    int step = firstRow >= 0 ? 1 : -1;
    int start = firstRow >= 0 ? firstRow : lastRow;
    int & otherEnd = firstRow >= 0 ? lastRow : firstRow;
    otherEnd = start;
    while (EventRows()[otherEnd] == event && otherEnd >= 0 && otherEnd < static_cast<int>(fRowEvents.size()))
      otherEnd += step;
    otherEnd -= step;  // undo the last step that failed

    if (firstRow < 0 || lastRow < 0 || lastRow < firstRow)
    {
      std::cerr << "firstRow and lastRow were not determined correctly." << std::endl;
      std::cerr << "  firstRow = " << firstRow << std::endl;
      std::cerr << "  lastRow = " << lastRow << std::endl;
      abort();
    }

    // assume that the track info starts after the "Event" column and reaches to the end
    std::size_t evtColumnIdx = std::distance(fColumnNames.begin(),
                                             std::find(fColumnNames.begin(), fColumnNames.end(), "Event"));

    // just double-check things are as we expect
    if (fColumnNames[evtColumnIdx + 1] != "trk_start_x")
    {
      std::cerr << "fColumnNames[" <<  evtColumnIdx+1 <<  + "] = " << fColumnNames[evtColumnIdx + 1] << std::endl;
      assert(false);
    }
    assert(fColumnNames[evtColumnIdx + 4] == "trk_end_x");
    assert(fColumnNames[evtColumnIdx + 7] == "trk_end_dir_x");


    // read the whole group of numbers together for efficiency
    std::size_t nRows = lastRow - firstRow + 1;
    std::size_t nCols = fColumnNames.size() - evtColumnIdx - 1;
    std::vector<Track> tracksOut;
    tracksOut.reserve(nRows);

    float vals[nRows * nCols];
    GridValues(vals, firstRow, lastRow, evtColumnIdx + 1, evtColumnIdx + nCols);

    std::cout << "  buffer now contains:" << std::endl;
    for (std::size_t idx = 0; idx < nRows * nCols; idx++)
      std::cout << "  " << vals[idx];
    std::cout << std::endl;


    //  fill the Tracks object
    std::cout << nRows << " tracks for event " << event << ":" << std::endl;
   for (std::size_t rowIdx = 0; rowIdx < nRows; rowIdx++)
   {
      Track tr;
      std::size_t rowOffset = rowIdx * nCols;
      std::cout << "row = " << rowIdx << "; rowoffset = " << rowOffset << std::endl;
      tr.start = {vals[rowOffset + 0], vals[rowOffset + 1], vals[rowOffset + 2]};
      tr.end = {vals[rowOffset + 3], vals[rowOffset + 4], vals[rowOffset + 5]};
      tr.end_dir = {vals[rowOffset + 6], vals[rowOffset + 7], vals[rowOffset + 8]};
      std::cout << "  " << tr << std::endl;

      tracksOut.emplace_back(std::move(tr));
    }

    //for (std::size_t i = 0; i < nRows; i++)
    //  delete [] vals[i];
    //delete [] vals;

    return tracksOut;
  }

  // -----------------------------------------------------------

  void NDLArSummaryH5::ReadColumnNames()
  {
    H5::Attribute attr = fInputDataset.openAttribute(fColumnNameAttr);
    auto dtype = attr.getDataType();
    auto dSpace = attr.getSpace();
    if (dtype.getClass() != H5T_STRING)
      throw std::runtime_error("Unexpected type for 'column_names' attribute in HDF5 file: " + std::to_string(dtype.getClass()));


    // the following modified from https://stackoverflow.com/questions/43722194/reading-a-string-array-hdf5-attribute-in-c
    hsize_t dim = 0;
    dSpace.getSimpleExtentDims(&dim);   // the number of strings in the array...

    fColumnNames.resize(dim);
    char **rdata = new char*[dim];
    attr.read(dtype, (void*)rdata);
    for(std::size_t iStr=0; iStr < dim; ++iStr)
    {
      fColumnNames[iStr].assign(rdata[iStr]);
      delete[] rdata[iStr];
    }
    delete[] rdata;

  }

} // namespace cafmaker


std::ostream &operator<<(std::ostream &stream, const TVector3 &vec)
{
  stream << "(" << vec.x() << "," << vec.y() << "," << vec.z() << ")";
  return stream;
}

std::ostream &operator<<(std::ostream &stream, const cafmaker::Track &tr)
{
  stream << "cafmaker::Track with start=" << tr.start << "  end=" << tr.end << "  end_dir=" << tr.end_dir << std::endl;
  return stream;
}
