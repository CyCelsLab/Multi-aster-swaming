//RCS: $Id: matsparsebandsym.h,v 2.1 2004/12/06 18:36:32 nedelec Exp $


#ifndef MATRIXSPARSEBANDSYM_H
#define MATRIXSPARSEBANDSYM_H

#include <cstdio>
#include "matrix.h"

/// a sparse real symmetric Matrix, with a variable number of full diagonals
class MatrixSparseBandSymmetric : public Matrix
{
  
  static const int SPCHUNCK=8;
  
private:
  
  int allocated_diag;
  
  // number of off-diagonals stored in val
  int diag;
  
  // the diagonal's values
  real * val; 
  
  // the off-diagonal's values: normal sparse
  
  // array [ size ][ ? ] holding the values for each column
  real ** Vcol;
  
  // array [ size ][ ? ] holding the line index for each column
  int  ** Vrow;
  
  void allocateColumn( int, int );
 
public:
   
 /// base for constructor
 void allocate();
 
 /// base for destructor
 void deallocate();
 
 /// default constructor
 MatrixSparseBandSymmetric()   { allocate(); } 
 
 /// default destructor
 ~MatrixSparseBandSymmetric()  { deallocate(); }
 
 /// set all the element to zero
 void setToZero();
 
 /// allocate the matrix, specifying the number of off diagonal
 void allocate( const int sz, const int diagonals );

 /// allocate the matrix to hold ( sz * sz )
 void allocate( const int sz ) {
   allocate( sz, 1 );
 }
  
 /// returns the address of element at (x, y), no allocation is done
 real * addr( int x, int y ) const;
 
 /// returns the address of element at (x, y), allocating if necessary
 real & operator()( int x, int y );
  
 /// scale the matrix by a scalar factor
 void scale( const real a );
  
 /// add the diagonal block ( x, x, x+sx, x+sx ) from this matrix to M
 void addDiagonalBlock(real* M, const int x, const int sx ) const;
 
 /// add the upper triagular block ( x, x, x+sx, x+sx ) from this matrix to M
 void addTriangularBlock(real* M, const int x, const int sx ) const;
  
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
