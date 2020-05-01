//RCS: $Id: matrix.h,v 2.0 2004/08/16 15:55:08 nedelec Exp $

#ifndef MATRIX_H
#define MATRIX_H

#include <cstdio>
#include "types.h"

/// class Matrix is the common parent for all the big matrices
class Matrix
{

protected:
  
  /// size of the matrix
  int size;
  
  /// size of memory which has been allocated
  int allocated;

public:
    
  /// base for constructor
  virtual void allocate() { }
  
  /// base for destructor
  virtual void deallocate() { }
  
  /// default constructor
  Matrix()           { allocate(); } 
  
  /// default destructor
  virtual ~Matrix()  { deallocate(); }
  
  //----------------------------------------------------------------------

  
  /// allocate the matrix to hold ( sz * sz ), all values are lost
  virtual void allocate(const int sz ) = 0;
  
  /// allocate ( sz * sz ), specifying the number of full diagonals
  void allocate(const int sz, const int diagonals ) {
    allocate( sz ); 
  }

  /// returns the address of element at (x, y), no allocation is done
  virtual real * addr( int x, int y ) const = 0;
  
  /// returns the address of element at (x, y), allocating if necessary
  virtual real & operator()( int x, int y) = 0;
  
  /// returns the value of element at (x, y) or zero if not allocated
  real value( int x, int y) const;
  
  //----------------------------------------------------------------------

  /// returns the size of the matrix
  virtual int getSize() const { return size; }
  
  /// set all the elements to zero
  virtual void setToZero() = 0;

  /// scale the matrix by a scalar factor
  virtual void scale( const real a ) = 0;
  
  /// copy the block ( x, y, x+sx, y+sy ) from this matrix into M
  virtual void copyBlock(real* M, const int x, const int sx, const int y, const int sy) const;
  
  /// add the diagonal block ( x, x, x+sx, x+sx ) from this matrix to M
  virtual void addDiagonalBlock(real* M, const int x, const int sx) const;
  
  /// add the upper triagular block ( x, x, x+sx, x+sx ) from this matrix to M
  virtual void addTriangularBlock(real* M, const int x, const int sx) const;
  
  //----------------------------------------------------------------------
  
  ///optional optimization to accelerate multiplications below
  virtual void optimizeForMultiply() {}
  
  /// multiplication of a vector: Y = Y + M * X, dim(X) = dim(M)
  virtual void multiplyVectorAdd( const real * X, real * Y ) const = 0;
  
  /// isotropic multiplication of a vector: Y = Y + M * X, dim(X) = DIM * dim(M)
  virtual void multiplyVectorAddIsotropic( const real * X, real * Y ) const = 0;
  
  /// multiplication of a vector: Y = M * X, dim(X) = dim(M)
  virtual void multiplyVector( const real * X, real * Y ) const;
  
  /// isotropic multiplication of a vector, Y = M * X, dim(X) = DIM * dim(M)
  virtual void multiplyVectorIsotropic( const real * X, real * Y ) const;

  //----------------------------------------------------------------------
  
  /// maximum of absolute value on all the elements
  virtual real maxNorm() const;
  
  /// true if the matrix is non-zero
  virtual bool nonZero() const;
  
  /// number of element which are non-zero
  virtual int  nbNonZeroElements() const;
  
  /// returns a string which a description of the type of matrix
  virtual char * shortDescription() const { return "Matrix"; }
  
  /// printf debug function in sparse mode: i, j : value
  virtual void printSparse(FILE * f=stdout) const;
  
  /// printf debug function in full lines, all columns
  virtual void printFull(FILE * f=stdout) const;
  
  /// debug function
  virtual void checkDataConsistency() const {}
};


#endif
