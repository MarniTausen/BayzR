//
//  simpleMatrix.h
//  Class to hold a matrix of doubles, quite basic functionality but designed for performance.
//  - all memory allocated as one block
//  - organized column-major so computations go fastest working by column
//  - accessible with double index ->data[col][row] - NOTE col first!
//    or alternatively, ->data[col] is a pointer to the data of column `col`
//  - constructor specifies number of rows first!
// Note: the actual elements need to be accessed as object->data[] (column) or ->data[][] (element),
//    at least the first one can be made more fancy by adding a function for operator[]?
//
//  Created by Luc Janss on 23/10/2020.
//

#ifndef simpleMatrix_h
#define simpleMatrix_h

#include <algorithm>

class simpleMatrix {

public:
   simpleMatrix(size_t nr, size_t nc) {
      data0 = new double[nr*nc];    // a block of memory nrows * ncolumns
      data  = new double*[nc];      // a list of column-pointers to index it by column
      for(size_t i=0; i<nc; i++)
         data[i] = data0 + i*nr;
      nrow = nr; ncol = nc;
      std::fill_n(data0, nr*nc, 0.0l);
   }

   ~simpleMatrix() {
      delete[] data;
      delete[] data0;
   }

   double *data0;
   double **data;
   size_t nrow,ncol;
   
};

#endif /* simpleMatrix_h */
