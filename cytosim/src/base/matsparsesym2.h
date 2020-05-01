//RCS: $Id: matsparsesym2.h,v 2.8 2005/04/24 22:48:02 nedelec Exp $

#ifndef MATSPARSESYM2_H
#define MATSPARSESYM2_H

#include <cstdio>
#include "matrix.h"


///real symmetric sparse Matrix, faster equivalent to MatrixSparseSymmetric
class MatrixSparseSymmetric2 : public Matrix
{

private:
  
  ///Element describes an element in a sparse matrix
  struct Element {
    real val;    ///< The value of the element
    int  indx;   ///< The index of the line (the elements are stored per columns, that is how the column index is known)
  };
  
  static const int SPCHUNCK = 8;

  ///accessory function for qsort:
  static int mat_comp( const void *, const void * );
    
protected:
 
  /// array [ size ][ indx ] holding the Element for the 'indx' column
  Element ** Vcol;
  
  ///number of elements in each column: Vsize[ indx ] is the number of elements in column 'indx'
  int  * Vsize;
  
  ///numer of element allocated in each column
  int  * Vallo;
  
  ///nextCol[ii] is the index of the first non-empty column of index strictly greater than ii
  int  * nextCol;
  
  ///index of the first non-empty column
  int    firstCol;
  
  /// allocate column to hold nb values
  void allocateColumn( int column_index, int nb_values );

public:
    
  /// base for constructor
  void allocate();
  
  /// base for destructor
  void deallocate();
  
  /// default constructor
  MatrixSparseSymmetric2()   { allocate(); } 
  
  /// default destructor
  ~MatrixSparseSymmetric2()  { deallocate(); }
  
  /// set all the element to zero
  void setToZero();
  
  /// allocate the matrix to hold ( sz * sz )
  void allocate( const int sz );
  
  /// allocate ( sz * sz ), specifying the number of full diagonals
  void allocate(const int sz, const int diagonals ) {
    allocate( sz ); 
  }
  
  /// returns the address of element at (x, y), no allocation is done
  real * addr( int x, int y ) const;
  
  /// returns the address of element at (x, y), allocating if necessary
  real & operator()( int x, int y );
   
  /// scale the matrix by a scalar factor
  void scale( const real a );
    
  /// add the diagonal block ( x, x, x+sx, x+sx ) from this matrix to M
  void addDiagonalBlock( real* M, const int x, const int sx) const;
  
  /// add the upper triagular block ( x, x, x+sx, x+sx ) from this matrix to M
  void addTriangularBlock( real* M, const int x, const int sx) const;
    
  ///optional optimization that may accelerate multiplications by a vector
  void optimizeForMultiply();
 
  /// multiplication of a vector: Y = Y + M * X, dim(X) = dim(M)
  void multiplyVectorAdd( const real * X, real * Y ) const;
  
  /// isotropic multiplication of a vector: Y = Y + M * X, dim(X) = DIM * dim(M)
  void multiplyVectorAddIsotropic( const real * X, real * Y ) const;
  
  /// true if matrix is non-zero
  bool nonZero() const;
    
  /// number of element which are non-zero
  int  nbNonZeroElements() const;
  
  /// returns a string which a description of the type of matrix
  char * shortDescription() const;
  
  /// printf debug function in sparse mode: i, j : value
  void printSparse(FILE * f=stdout) const;
    
  /// debug function
  void checkDataConsistency() const;
};


#endif
