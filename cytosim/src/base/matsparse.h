//RCS: $Id: matsparse.h,v 2.1 2004/12/06 18:36:27 nedelec Exp $
//======================================================================
//===                          matsparse.h                           ===
//===                real symmetric sparse matrix,                   ===
//======================================================================

#ifndef MATSPARSE_H
#define MATSPARSE_H

#include <cstdio>
#include "matrix.h"

/// a real (non-symmetric) sparse Matrix
class MatrixSparse : public Matrix
{

private:
  
  static const int SPCHUNCK = 8;

protected:
 
  // array [ size ][ ? ] holding the values for each column
  real ** Vcol;
  
  // array [ size ][ ? ] holding the line index for each column
  int  ** Vrow;

  // allocate column to hold nb values
  void allocateColumn( int column_index, int nb_values );

public:
    
  /// base for constructor
  void allocate();
  
  /// base for destructor
  void deallocate();
  
  /// default constructor
  MatrixSparse()   { allocate(); } 
  
  /// default destructor
  ~MatrixSparse()  { deallocate(); }
  
  /// set all the element to zero
  void setToZero();
  
  /// allocate the matrix to hold ( sz * sz )
  void allocate( const int sz );
  
  /// allocate ( sz * sz ), specifying the number of full diagonals
  void allocate( const int sz, const int diagonals ) {
    allocate( sz ); 
  }
  
  /// returns the address of element at (x, y), no allocation is done
  real * addr( int x, int y ) const;
  
  /// returns the address of element at (x, y), allocating if necessary
  real & operator()( int x, int y );
   
  /// scale the matrix by a scalar factor
  void scale( const real a );
    
  /// add the diagonal block ( x, x, x+sx, x+sx ) from this matrix to M
  void addDiagonalBlock(real* M, const int x, const int sx) const;
  
  /// add the upper triagular block ( x, x, x+sx, x+sx ) from this matrix to M
  void addTriangularBlock(real* M, const int x, const int sx) const;
    
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
