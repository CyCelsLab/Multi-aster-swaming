//RCS: $Id: matsym.h,v 2.1 2004/12/06 18:36:50 nedelec Exp $
//======================================================================
//===                          matsym.h                              ===
//======================================================================
#ifndef MATSYM_H
#define MATSYM_H

#include <cstdio>
#include "matrix.h"

///a real symmetric Matrix, the full upper triangular is stored
class MatrixSymmetric : public Matrix
{
  
protected:

  // full upper triangle:
  real * val;

  // will not call delete[]
  int do_not_delete_array;
  
public:

  /// base for constructor
  void allocate();
  
  /// base for destructor
  void deallocate();
  
  /// default constructor
  MatrixSymmetric()   { allocate(); } 
  
  /// constructor from an existing array
  MatrixSymmetric( int sz, real * array ) {
    size = sz;
    val = array;
    do_not_delete_array = 1;
  }
  
  /// default destructor
  ~MatrixSymmetric()  { deallocate(); }
  
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
};  

#endif //MAT_H
